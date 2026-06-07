"""Skill 主入口与管线编排.

VCF → Prefilter → VEP Annotate → Sensory Filter → Impact Assess →
Specialized Logic → API Enrichment → Report Context → Markdown/JSON Generator → Output
"""

import argparse
import asyncio
import json
import logging
import os
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from src.assessment.engine import ImpactEngine
from src.config_loader import load_config
from src.enrichment.cache import CacheManager
from src.enrichment.clinvar import ClinVarClient
from src.enrichment.gnomad import GnomADClient
from src.enrichment.gtex import GTExClient
from src.enrichment.uniprot import UniProtClient
from src.gene_sets.bed_mapper import BedMapper
from src.gene_sets.filter import SensoryFilter
from src.gene_sets.loader import GeneSetLoader
from src.logger import get_logger, setup_logger
from src.models import (
    AnalysisConfig,
    DataAvailability,
    ExecutiveSummary,
    GeneCard,
    SensoryReport,
    Variant,
)
from src.report.json_generator import JsonReportGenerator
from src.report.markdown_generator import MarkdownReportGenerator
from src.specialized.mitochondrial import MitochondrialAnnotator
from src.specialized.or_tiers import ORTierClassifier
from src.specialized.tas2r38 import TAS2R38Analyzer
from src.vcf.parser import VcfParser
from src.vcf.prefilter import Prefilter
from src.vep.client import VepClient
from src.vep.hybrid_client import HybridVepClient
from src.vep.local_client import LocalVepClient

logger = get_logger(__name__)


class SensoryPipeline:
    """感官基因组学分析管线."""

    def __init__(self, config: AnalysisConfig) -> None:
        """初始化管线.

        Args:
            config: 分析配置。
        """
        self.config = config
        self.skill_config = load_config()
        self.vcf_parser = VcfParser(config.vcf_path, config.sex)
        self.prefilter = Prefilter(
            min_qual=self.skill_config.filter.min_qual,
            min_dp=self.skill_config.filter.min_dp,
            pass_only=self.skill_config.filter.pass_only,
        )
        self.gene_sets = GeneSetLoader()
        self.sensory_filter = SensoryFilter(self.gene_sets)
        self.impact_engine = ImpactEngine()
        self.bed_mapper: Optional[BedMapper] = None

        # 专用逻辑模块
        self.tas2r38_analyzer = TAS2R38Analyzer()
        self.mt_annotator = MitochondrialAnnotator()
        self.or_classifier = ORTierClassifier()

        # 缓存与 API 富集客户端（可选）
        self.cache = CacheManager(
            db_path=self.skill_config.cache.db_path,
            default_ttl_days=self.skill_config.cache.default_ttl_days,
        )

        # 预计算数据：gene -> canonical_transcript_id 映射
        self._gene_transcript_map: Dict[str, str] = {}
        self._load_precomputed_data()

        # 根据配置选择 VEP 客户端
        vep_source = self.skill_config.vep.source
        if vep_source == "local_docker":
            logger.info("Using LocalVepClient (Docker offline)")
            self.vep_client = LocalVepClient(
                cache_dir=self.skill_config.vep.cache_dir,
                species="homo_sapiens",
                assembly="GRCh38",
                cache_version=115,
            )
        elif vep_source == "rest_api":
            logger.info("Using VepClient (REST API)")
            self.vep_client = VepClient(
                base_url=self.skill_config.vep.base_url,
                batch_size=self.skill_config.vep.batch_size,
                rate_limit=self.skill_config.vep.rate_limit,
                max_retries=self.skill_config.vep.max_retries,
                timeout=self.skill_config.vep.timeout,
                cache=self.cache,
                canonical_only=True,
                transcript_ids=list(self._gene_transcript_map.values()) if self._gene_transcript_map else None,
                precompute_db=self.skill_config.precompute_db,
            )
        else:
            # auto（默认）：HybridVepClient，自动检测本地环境，fallback REST
            logger.info("Using HybridVepClient (auto-detect local → fallback REST)")
            self.vep_client = HybridVepClient(
                local_config={
                    "cache_dir": self.skill_config.vep.cache_dir,
                    "species": "homo_sapiens",
                    "assembly": "GRCh38",
                    "cache_version": 115,
                },
                rest_config={
                    "base_url": self.skill_config.vep.base_url,
                    "batch_size": self.skill_config.vep.batch_size,
                    "rate_limit": self.skill_config.vep.rate_limit,
                    "max_retries": self.skill_config.vep.max_retries,
                    "timeout": self.skill_config.vep.timeout,
                    "cache": self.cache,
                    "canonical_only": True,
                    "transcript_ids": list(self._gene_transcript_map.values()) if self._gene_transcript_map else None,
                    "precompute_db": self.skill_config.precompute_db,
                },
            )
        self.uniprot_client = UniProtClient(
            cache=self.cache,
            rate_limit=self.skill_config.rate_limits.uniprot,
        )
        self.gnomad_client = GnomADClient(
            cache=self.cache,
            rate_limit=self.skill_config.rate_limits.ensembl,
        )
        self.clinvar_client = ClinVarClient(
            cache=self.cache,
            rate_limit=self.skill_config.rate_limits.ncbi,
        )
        self.gtex_client = GTExClient(
            cache=self.cache,
            rate_limit=self.skill_config.rate_limits.uniprot,
        )

    def _load_precomputed_data(self) -> None:
        """加载预计算的基因坐标和转录本映射数据."""
        data_dir = Path(__file__).resolve().parent.parent / "data"

        # 1. gene -> canonical transcript 映射
        transcript_map_path = data_dir / "gene_transcript_map.json"
        if transcript_map_path.exists():
            try:
                with open(transcript_map_path, "r") as f:
                    self._gene_transcript_map = json.load(f)
                logger.info("Loaded precomputed transcript map: %d genes", len(self._gene_transcript_map))
            except Exception as exc:
                logger.warning("Failed to load transcript map: %s", exc)

        # 2. 精确外显子 BED（用于 strict 模式预过滤）
        exon_bed_path = data_dir / "sensory_genes_exons.bed"
        if exon_bed_path.exists():
            self._exon_bed_mapper = BedMapper(str(exon_bed_path))
            logger.info("Loaded exon BED for strict filtering")
        else:
            self._exon_bed_mapper = None

    async def run(self) -> SensoryReport:
        """执行完整分析管线.

        Returns:
            SensoryReport 报告对象。
        """
        logger.info("Pipeline started for %s", self.config.vcf_path)

        # Stage 1: VCF 解析与预过滤
        raw_variants = self._parse_and_prefilter()
        logger.info("Stage 1 complete: %d variants after prefilter", len(raw_variants))

        # Stage 2: VEP 批量注释
        annotated_variants = await self._vep_annotate(raw_variants)
        logger.info("Stage 2 complete: %d variants annotated", len(annotated_variants))

        # Stage 2.5: BED 映射补充基因名（当 VEP 不可用时）
        if self.bed_mapper:
            mapped = 0
            for variant in annotated_variants:
                if not variant.gene_symbol:
                    gene = self.bed_mapper.lookup(variant.chrom, variant.pos)
                    if gene:
                        variant.gene_symbol = gene
                        mapped += 1
            logger.info("Stage 2.5 complete: BED mapped %d variants to genes", mapped)

        # Stage 3: 感官基因精确筛选
        sensory_variants = self.sensory_filter.filter_variants(annotated_variants)
        logger.info("Stage 3 complete: %d sensory variants", len(sensory_variants))

        # Stage 4: 按基因分组并评估
        gene_groups = self.sensory_filter.group_by_gene(sensory_variants)
        gene_cards = self._assess_genes(gene_groups)
        logger.info("Stage 4 complete: %d gene cards", len(gene_cards))

        # Stage 5: 专用逻辑分析
        tas2r38_result = self._analyze_tas2r38(gene_cards)
        mt_results = self._analyze_mitochondrial(gene_cards)
        or_tier_results = self._classify_or_genes(gene_cards)
        logger.info(
            "Stage 5 complete: TAS2R38=%s, MT=%d, OR=%d",
            "yes" if tas2r38_result else "no",
            len(mt_results) if mt_results else 0,
            len(or_tier_results) if or_tier_results else 0,
        )

        # Stage 6: API 富集（异步并发）
        enrichment_data = await self._enrich_genes(gene_cards)
        logger.info("Stage 6 complete: enriched %d genes", len(enrichment_data))

        # Stage 7: 构建报告上下文
        report = self._build_report(
            gene_cards=gene_cards,
            tas2r38=tas2r38_result,
            mt_results=mt_results,
            or_tier_results=or_tier_results,
            enrichment_data=enrichment_data,
        )
        logger.info("Stage 7 complete: report built")

        logger.info("Pipeline completed successfully")
        return report

    def _parse_and_prefilter(self) -> List[Variant]:
        """解析 VCF 并应用预过滤."""
        variants = []
        for variant in self.vcf_parser.iter_variants():
            if self.prefilter.apply(variant):
                variants.append(variant)
        return variants

    async def _vep_annotate(self, variants: List[Variant]) -> List[Variant]:
        """VEP 注释并合并结果."""
        if not variants:
            return []

        vep_results = await self.vep_client.annotate(variants)
        annotated = []
        for variant, vep_data in zip(variants, vep_results):
            merged = self._merge_vep(variant, vep_data)
            annotated.append(merged)
        return annotated

    # 后果严重度分级（高 → 低），用于跨基因选最优转录本
    _CONSEQUENCE_RANK: Dict[str, int] = {
        "transcript_ablation": 40, "splice_acceptor_variant": 39,
        "splice_donor_variant": 38, "stop_gained": 37,
        "frameshift_variant": 36, "stop_lost": 35,
        "start_lost": 34, "transcript_amplification": 33,
        "inframe_insertion": 32, "inframe_deletion": 31,
        "missense_variant": 30, "protein_altering_variant": 29,
        "splice_region_variant": 28, "splice_donor_5th_base_variant": 27,
        "splice_donor_region_variant": 26, "splice_polypyrimidine_tract_variant": 25,
        "incomplete_terminal_codon_variant": 24, "start_retained_variant": 23,
        "stop_retained_variant": 22, "synonymous_variant": 21,
        "coding_sequence_variant": 20, "mature_miRNA_variant": 19,
        "5_prime_UTR_variant": 18, "3_prime_UTR_variant": 17,
        "non_coding_transcript_exon_variant": 16, "intron_variant": 15,
        "NMD_transcript_variant": 14, "non_coding_transcript_variant": 13,
        "upstream_gene_variant": 12, "downstream_gene_variant": 11,
        "TFBS_ablation": 10, "TFBS_amplification": 9,
        "TF_binding_site_variant": 8, "regulatory_region_ablation": 7,
        "regulatory_region_amplification": 6, "feature_elongation": 5,
        "regulatory_region_variant": 4, "feature_truncation": 3,
        "intergenic_variant": 2,
    }

    @classmethod
    def _score_transcript(cls, tc: Dict[str, Any], gene: str = "") -> int:
        """对单个转录本打分，用于跨基因选择最佳转录本.

        评分规则:
            - 基因匹配: +100
            - 典型转录本 (canonical): +10
            - protein_coding 生物型: +5
            - 后果严重度 (_CONSEQUENCE_RANK): +0-40
        """
        score = 0
        if gene and tc.get("gene_symbol", "") == gene:
            score += 100
        flags = tc.get("flags") or []
        if tc.get("canonical") or "CANONICAL" in flags:
            score += 10
        biotype = tc.get("biotype", "")
        if "protein_coding" in biotype:
            score += 5
        max_rank = 0
        for term in tc.get("consequence_terms", []):
            max_rank = max(max_rank, cls._CONSEQUENCE_RANK.get(term, 0))
        score += max_rank
        return score

    @classmethod
    def _select_best_transcript(
        cls, consequences: List[Dict[str, Any]], gene: str = ""
    ) -> Optional[Dict[str, Any]]:
        """从转录本列表中选出最优的一个.

        Args:
            consequences: VEP transcript_consequences 列表
            gene: 目标基因符号（空字符串 = 跨基因选择，不偏袒任何基因）

        Returns:
            最优转录本 dict，若列表为空则返回 None
        """
        if not consequences:
            return None
        return max(consequences, key=lambda tc: cls._score_transcript(tc, gene))

    @classmethod
    def _merge_vep(cls, variant: Variant, vep_data: Dict[str, Any]) -> Variant:
        """将 VEP 结果合并到 Variant 模型.

        对基因重叠区域的变异（如 NLRP3/OR2B11 共区），选择后果最严重的
        代表性转录本，而非盲目取 transcript_consequences[0]。
        """
        if not vep_data:
            return variant

        consequences = vep_data.get("transcript_consequences", [])
        if consequences:
            # 选出整体最优转录本（不偏袒任何基因，按严重度排名）
            best_tc = cls._select_best_transcript(consequences, gene="")

            if best_tc:
                variant.gene_symbol = best_tc.get("gene_symbol", variant.gene_symbol)
                cons_terms = best_tc.get("consequence_terms", [])
                variant.consequence = ",".join(cons_terms) if cons_terms else ""
                variant.hgvsc = best_tc.get("hgvsc")
                variant.hgvsp = best_tc.get("hgvsp")
                variant.protein_domain = (
                    best_tc.get("protein_domains", [None])[0]
                    if best_tc.get("protein_domains")
                    else None
                )

        # 提取 colocated variants（gnomAD AF, ClinVar 等）
        colocated = vep_data.get("colocated_variants", [])
        for cv in colocated:
            if "gnomad" in str(cv.get("allele_frequency", "")).lower():
                variant.af_gnomad = cv.get("allele_frequency")

        variant.raw_vep = vep_data
        return variant

    def _assess_genes(self, gene_groups: Dict[str, List[Variant]]) -> List[GeneCard]:
        """对每个基因进行功能影响评估."""
        gene_cards = []
        for gene_symbol, variants in gene_groups.items():
            if not variants:
                continue

            # 取影响程度最高的变异作为代表
            assessments = [
                self.impact_engine.assess(v, self.config.sex, self.gene_sets)
                for v in variants
            ]
            best_assessment = max(
                assessments,
                key=lambda a: self._level_rank(a.level),
            )

            gene_card = GeneCard(
                gene_symbol=gene_symbol,
                subsystem=self.gene_sets.get_subsystem(gene_symbol),
                sensory_function_zh=self.gene_sets.get_gene_function(gene_symbol),
                variants=variants,
                assessment=best_assessment,
            )
            gene_cards.append(gene_card)

        return gene_cards

    def _analyze_tas2r38(self, gene_cards: List[GeneCard]) -> Optional[Any]:
        """分析 TAS2R38 Haplotype."""
        tas2r38_variants = []
        for card in gene_cards:
            if card.gene_symbol == "TAS2R38":
                tas2r38_variants.extend(card.variants)

        if not tas2r38_variants:
            return None

        try:
            return self.tas2r38_analyzer.analyze(tas2r38_variants)
        except Exception as exc:
            logger.error("TAS2R38 analysis failed: %s", exc)
            return None

    def _analyze_mitochondrial(
        self, gene_cards: List[GeneCard]
    ) -> Optional[List[Any]]:
        """分析线粒体耳聋变异."""
        mt_variants = []
        for card in gene_cards:
            if card.gene_symbol.startswith("MT-"):
                mt_variants.extend(card.variants)

        if not mt_variants:
            return None

        try:
            return self.mt_annotator.annotate(mt_variants)
        except Exception as exc:
            logger.error("Mitochondrial annotation failed: %s", exc)
            return None

    def _classify_or_genes(self, gene_cards: List[GeneCard]) -> Optional[List[Any]]:
        """对 OR 基因进行分级分类."""
        or_variants = []
        for card in gene_cards:
            if card.gene_symbol.startswith("OR"):
                or_variants.extend(card.variants)

        if not or_variants:
            return None

        try:
            return self.or_classifier.classify(
                or_variants, self.config.sex, self.gene_sets
            )
        except Exception as exc:
            logger.error("OR tier classification failed: %s", exc)
            return None

    # Impact levels that trigger API enrichment (部分影响及以上)
    _ENRICH_LEVELS = frozenset({"完全丧失", "显著影响", "部分影响"})

    # Per-API concurrency limits to avoid rate limiting
    _API_SEMAPHORES: Dict[str, asyncio.Semaphore] = {}

    async def _enrich_genes(
        self, gene_cards: List[GeneCard]
    ) -> Dict[str, Dict[str, Any]]:
        """异步并发查询外部 API 富集基因信息（仅对部分影响及以上的基因）."""
        if not self.config.show_reference_info:
            return {}

        # Filter: only enrich genes at or above 部分影响
        priority_genes = [
            card.gene_symbol for card in gene_cards
            if card.gene_symbol and card.assessment.level in self._ENRICH_LEVELS
        ]
        if not priority_genes:
            logger.info("Stage 6: no genes meet impact threshold, skipping API enrichment")
            return {}

        logger.info(
            "Stage 6: enriching %d/%d genes (impact >= 部分影响)",
            len(priority_genes), len(gene_cards),
        )

        # Per-API concurrency: gnomAD=3, GTEx=5, ClinVar=5, UniProt=10
        api_limits = {"gnomad": 3, "gtex": 5, "clinvar": 5, "uniprot": 10}

        async def rate_limited_query(client, api_name: str, gene: str):
            sem = self._API_SEMAPHORES.setdefault(
                api_name, asyncio.Semaphore(api_limits.get(api_name, 5))
            )
            async with sem:
                return await self._safe_query(client.query, gene)

        all_tasks = []
        task_meta = []

        for gene in priority_genes:
            all_tasks.append(rate_limited_query(self.uniprot_client, "uniprot", gene))
            all_tasks.append(rate_limited_query(self.gnomad_client, "gnomad", gene))
            all_tasks.append(rate_limited_query(self.clinvar_client, "clinvar", gene))
            all_tasks.append(rate_limited_query(self.gtex_client, "gtex", gene))
            task_meta.extend([
                (gene, "uniprot"),
                (gene, "gnomad"),
                (gene, "clinvar"),
                (gene, "gtex"),
            ])

        results = await asyncio.gather(*all_tasks, return_exceptions=True)

        enrichment: Dict[str, Dict[str, Any]] = {gene: {} for gene in priority_genes}
        for (gene, api_name), result in zip(task_meta, results):
            if isinstance(result, Exception):
                enrichment[gene][api_name] = {"error": str(result)}
            else:
                enrichment[gene][api_name] = result

        return enrichment

    @staticmethod
    async def _safe_query(query_func, key: str) -> Any:
        """安全查询，异常时返回降级结果."""
        try:
            return await query_func(key)
        except Exception as exc:
            return {"error": str(exc), "api": "unknown", "key": key}

    def _build_report(
        self,
        gene_cards: List[GeneCard],
        tas2r38: Optional[Any],
        mt_results: Optional[List[Any]],
        or_tier_results: Optional[List[Any]],
        enrichment_data: Dict[str, Dict[str, Any]],
    ) -> SensoryReport:
        """构建完整报告."""
        # 注入富集数据
        for card in gene_cards:
            gene = card.gene_symbol
            if gene in enrichment_data:
                card.enrichment_data = enrichment_data[gene]

        # 构建执行摘要
        executive_summary = self._build_executive_summary(gene_cards)

        # 构建数据可用性
        data_availability: Dict[str, DataAvailability] = {}
        for card in gene_cards:
            data = enrichment_data.get(card.gene_symbol, {})
            data_availability[card.gene_symbol] = DataAvailability(
                gnomad_af="可用" if data.get("gnomad", {}).get("found") else "N/A",
                clinvar="可用" if data.get("clinvar", {}).get("found") else "N/A",
                spliceai="N/A",
                cadd="N/A",
                topology="可用" if data.get("uniprot", {}).get("found") else "N/A",
            )

        report = SensoryReport(
            sample_id=self.vcf_parser.sample_name,
            sex=self.config.sex,
            ref_genome="GRCh38",
            subsystems=self.config.subsystems,
            gene_cards=gene_cards,
            tas2r38=tas2r38,
            mitochondrial=mt_results,
            or_tiers=or_tier_results,
            executive_summary=executive_summary,
            data_availability=data_availability,
            disclaimer_zh=self._default_disclaimer(),
        )
        return report

    @staticmethod
    def _build_executive_summary(gene_cards: List[GeneCard]) -> ExecutiveSummary:
        """构建执行摘要统计."""
        subsystem_counts: Dict[str, Dict[str, int]] = {}
        key_findings = []

        for card in gene_cards:
            level = card.assessment.level
            if level in ("完全丧失", "显著影响"):
                key_findings.append(
                    f"{card.gene_symbol}: {level}（{card.assessment.rationale_zh[:50]}...）"
                )

        return ExecutiveSummary(
            subsystem_counts=subsystem_counts,
            key_findings=key_findings,
        )

    # 预计算的影响程度排序字典（避免 _level_rank 每次重建）
    _LEVEL_RANK: Dict[str, int] = {
        "无影响": 0,
        "可能轻微影响": 1,
        "部分影响": 2,
        "显著影响": 3,
        "完全丧失": 4,
    }

    @classmethod
    def _level_rank(cls, level: str) -> int:
        """影响程度排序（越大越严重）."""
        return cls._LEVEL_RANK.get(level, 0)

    @staticmethod
    def _default_disclaimer() -> str:
        """默认免责声明文本."""
        return (
            "【重要声明】本报告仅基于基因组变异数据进行描述性分析，不构成医学诊断，"
            "不可替代专业医疗建议或临床检查。报告中的\"影响\"描述仅表示基因层面的功能预测，"
            "不等同于临床表现或疾病状态。感官能力受基因、环境、生活方式等多因素影响，"
            "基因型与表型之间不存在简单的一一对应关系。如报告中有高风险发现，"
            "建议咨询遗传咨询师或相关专科医生。未成年人应在监护人陪同下解读报告。"
        )

    async def close(self) -> None:
        """释放资源（带超时保护，防止 aiohttp session 挂死）."""
        _CLOSE_TIMEOUT = 5.0
        clients = [
            ("vep", self.vep_client),
            ("uniprot", self.uniprot_client),
            ("gnomad", self.gnomad_client),
            ("clinvar", self.clinvar_client),
            ("gtex", self.gtex_client),
        ]
        for name, client in clients:
            try:
                await asyncio.wait_for(client.close(), timeout=_CLOSE_TIMEOUT)
            except asyncio.TimeoutError:
                logger.warning("Close timeout for %s client (%.1fs), skipping", name, _CLOSE_TIMEOUT)
            except Exception as exc:
                logger.warning("Close error for %s client: %s", name, exc)


async def run_analysis(config: AnalysisConfig, bed_path: Optional[str] = None) -> SensoryReport:
    """运行分析管线并返回报告.

    Args:
        config: 分析配置。
        bed_path: BED 文件路径，用于坐标到基因名的映射。

    Returns:
        SensoryReport 报告对象。
    """
    pipeline = SensoryPipeline(config)
    if bed_path:
        pipeline.bed_mapper = BedMapper(bed_path)
    try:
        report = await pipeline.run()
    finally:
        await pipeline.close()
    return report


def _filter_vcf_by_bed(vcf_path: str, strict_filter: bool = False) -> tuple[str, Optional[str]]:
    """用 BED 文件筛选感官基因区域的变异，返回临时 VCF 路径和使用的 BED 路径.

    Args:
        vcf_path: 输入 VCF 路径。
        strict_filter: 若 True，优先使用精确外显子/CDS BED 做硬过滤。

    Returns:
        (filtered_vcf_path, bed_path_used)
    """
    script_dir = Path(__file__).resolve().parent
    work_dir = Path(__file__).resolve().parent.parent

    # 精确坐标 BED（外显子/CDS）—— strict 模式优先
    precise_beds = []
    if strict_filter:
        precise_beds = [
            script_dir / "data" / "sensory_genes_exons.bed",
            script_dir.parent / "assets" / "data" / "sensory_genes_exons.bed",
            work_dir / "sensory_genes_exons.bed",
            work_dir / "assets" / "data" / "sensory_genes_exons.bed",
        ]

    # 基因全区域 BED —— fallback
    broad_beds = [
        script_dir / "data" / "sensory_gene_regions.bed",
        script_dir.parent / "assets" / "data" / "sensory_gene_regions.bed",
        work_dir / "sensory_genes.bed",
        work_dir / "assets" / "data" / "sensory_gene_regions.bed",
    ]

    # 外部 BED（通用路径，无用户特定硬编码）
    external_beds = [
        "/tmp/core_genes.bed",
    ]

    bed_path = None
    bed_source = "broad"

    # 先找精确 BED（strict 模式）
    if strict_filter:
        for bed in precise_beds:
            if bed.exists():
                bed_path = str(bed)
                bed_source = "precise (exon/CDS)"
                logger.info("Using strict filter BED: %s", bed_path)
                break
        if bed_path is None:
            logger.warning("Strict filter requested but no exon BED found, falling back to broad BED")

    # 再找基因全区域 BED
    if bed_path is None:
        for bed in broad_beds + [Path(p) for p in external_beds]:
            if bed.exists():
                bed_path = str(bed)
                bed_source = "broad (gene region)"
                break

    if bed_path is None:
        logger.warning("No BED file found, using full VCF")
        return vcf_path, None

    # 确保 VCF 有 .tbi 索引（bcftools view -R 需要）
    if not os.path.exists(vcf_path + ".tbi") and not os.path.exists(vcf_path + ".csi"):
        logger.info("VCF index missing, creating .tbi for %s", vcf_path)
        try:
            subprocess.run(
                ["bcftools", "index", "-t", vcf_path],
                stderr=subprocess.PIPE,
                check=True,
                timeout=120,
            )
            logger.info("Created .tbi index for %s", vcf_path)
        except subprocess.CalledProcessError as exc:
            logger.error("Failed to create .tbi index: %s", exc.stderr.decode() if exc.stderr else str(exc))
            return vcf_path, None
        except subprocess.TimeoutExpired:
            logger.error("Index creation timed out for %s", vcf_path)
            return vcf_path, None

    # 检查 VCF 染色体命名
    has_chr = False
    try:
        result = subprocess.run(
            ["bcftools", "view", vcf_path],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, timeout=10
        )
        for line in result.stdout.split("\n"):
            if line.startswith("chr"):
                has_chr = True
                break
            elif line and not line.startswith("#"):
                break
    except Exception:
        pass

    # 创建带 chr 前缀的 BED 文件（如果需要）
    with open(bed_path, "r") as f:
        first_line = f.readline()
    bed_has_chr = first_line.startswith("chr")

    if has_chr and not bed_has_chr:
        # 为 BED 添加 chr 前缀
        tmp_bed = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
        with open(bed_path, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if parts:
                    parts[0] = "chr" + parts[0]
                    tmp_bed.write("\t".join(parts) + "\n")
        tmp_bed.close()
        bed_path = tmp_bed.name

    # 筛选 VCF
    tmp_vcf = tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False)
    tmp_vcf.close()

    try:
        subprocess.run(
            ["bcftools", "view", "-Oz", "-o", tmp_vcf.name, "-R", bed_path, vcf_path],
            stderr=subprocess.PIPE,
            check=True,
            timeout=60,
        )
        # 创建索引
        subprocess.run(
            ["bcftools", "index", "-t", tmp_vcf.name],
            stderr=subprocess.PIPE,
            check=True,
            timeout=30,
        )
        count = subprocess.run(
            ["bcftools", "view", tmp_vcf.name],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, timeout=30
        )
        variant_count = sum(1 for line in count.stdout.split("\n") if line and not line.startswith("#"))
        logger.info("BED-filtered VCF: %s -> %d variants", tmp_vcf.name, variant_count)
        if variant_count == 0:
            logger.warning("No variants in sensory gene regions, using full VCF")
            os.unlink(tmp_vcf.name)
            return vcf_path, None
        return tmp_vcf.name, bed_path
    except subprocess.CalledProcessError as exc:
        logger.error("BED filter failed: %s", exc)
        if os.path.exists(tmp_vcf.name):
            os.unlink(tmp_vcf.name)
        return vcf_path, None
    except subprocess.TimeoutExpired:
        logger.error("BED filter timeout")
        if os.path.exists(tmp_vcf.name):
            os.unlink(tmp_vcf.name)
        return vcf_path, None


async def main(config: Optional[AnalysisConfig] = None) -> Dict[str, str]:
    """Skill 主入口.

    Args:
        config: 分析配置，若为 None 则解析命令行参数。

    Returns:
        输出文件路径字典 {"markdown": ..., "json": ...}。
    """
    if config is None:
        config = _parse_args()

    # 初始化日志
    setup_logger(level=logging.INFO)

    # BED 筛选感官基因区域
    filtered_vcf, bed_path_used = _filter_vcf_by_bed(config.vcf_path, strict_filter=config.strict_filter)
    if filtered_vcf != config.vcf_path:
        config.vcf_path = filtered_vcf

    # 运行分析管线
    report = await run_analysis(config, bed_path=bed_path_used)

    # 生成报告
    md_generator = MarkdownReportGenerator()
    json_generator = JsonReportGenerator()

    markdown_content = md_generator.generate(report)
    json_content = json_generator.generate(report)

    # 保存报告
    output_dir = config.output_dir or os.path.expanduser(
        "~/.workbuddy/skills/sensory-genomics/output/"
    )
    os.makedirs(output_dir, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    md_path = os.path.join(output_dir, f"sensory_report_{timestamp}.md")
    json_path = os.path.join(output_dir, f"sensory_report_{timestamp}.json")

    with open(md_path, "w", encoding="utf-8") as f:
        f.write(markdown_content)

    with open(json_path, "w", encoding="utf-8") as f:
        f.write(json_content)

    logger.info("Reports saved: markdown=%s, json=%s", md_path, json_path)
    return {"markdown": md_path, "json": json_path}


def _parse_args() -> AnalysisConfig:
    """解析命令行参数."""
    parser = argparse.ArgumentParser(description="Sensory Genomics Analysis Skill")
    parser.add_argument("--vcf", required=True, help="VCF file path")
    parser.add_argument("--sex", required=True, choices=["M", "F"], help="Sample sex")
    parser.add_argument(
        "--subsystems",
        default="vision,hearing,olfaction,taste,somatosensation",
        help="Comma-separated subsystems",
    )
    parser.add_argument("--output-dir", default=None, help="Output directory")
    parser.add_argument("--known-phenotype", default=None, help="Known phenotype text")
    parser.add_argument(
        "--no-reference-info",
        action="store_true",
        help="Skip API reference info enrichment",
    )
    parser.add_argument(
        "--strict-filter",
        action="store_true",
        help="Use precise exon/CDS coordinates for pre-filtering (excludes intronic variants)",
    )

    args = parser.parse_args()
    return AnalysisConfig(
        vcf_path=args.vcf,
        sex=args.sex,
        subsystems=args.subsystems.split(","),
        output_dir=args.output_dir,
        known_phenotype=args.known_phenotype,
        show_reference_info=not args.no_reference_info,
        strict_filter=args.strict_filter,
    )


if __name__ == "__main__":
    asyncio.run(main())
