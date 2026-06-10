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
import sys
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
from src.enrichment.local_gnomad import LocalGnomADClient
from src.enrichment.uniprot import UniProtClient
from src.gene_sets.bed_mapper import BedMapper
from src.gene_sets.filter import SensoryFilter
from src.gene_sets.loader import GeneSetLoader
from src.key_snps import KeySNPInferrer
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
from src.vcf.sex_detector import detect_sex
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
        self.bed_mapper = BedMapper()

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
        self.local_gnomad = LocalGnomADClient()

    def _load_precomputed_data(self) -> None:
        """加载预计算的基因坐标和转录本映射数据."""
        data_dir = Path(__file__).resolve().parent.parent / "assets" / "data"

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

        # Stage 5.5: 关键性状 SNP 基因型推断（使用原始 VCF，非 BED 过滤版）
        key_snp_results = None
        snp_vcf_path = self.config.original_vcf_path or self.config.vcf_path
        try:
            inferrer = KeySNPInferrer(snp_vcf_path)
            key_snp_results = inferrer.infer_all()
            logger.info(
                "Stage 5.5 complete: %d key SNPs inferred (%d found in VCF, %d REF/REF) [VCF: %s]",
                len(key_snp_results),
                sum(1 for r in key_snp_results if r.found_in_vcf),
                sum(1 for r in key_snp_results if not r.found_in_vcf),
                "original" if self.config.original_vcf_path else "filtered",
            )
        except Exception as exc:
            logger.warning("Key SNP inference failed: %s", exc)

        # Stage 6: API 富集（异步并发）
        enrichment_data = await self._enrich_genes(gene_cards)
        logger.info("Stage 6 complete: enriched %d genes", len(enrichment_data))

        # Stage 6.3: gnomAD AF 降级评估（对所有级别生效）
        downgraded_count = self._apply_gnomad_af_downgrade(gene_cards, enrichment_data)
        logger.info("Stage 6.3 complete: %d genes downgraded by gnomAD AF", downgraded_count)

        # Stage 6.5: OR 基因分级重分类（使用 enrichment 数据）
        or_tier_results = self._reclassify_or_genes_with_enrichment(
            gene_cards, enrichment_data, or_tier_results
        )
        logger.info("Stage 6.5 complete: OR tiers reclassified with enrichment data")

        # Stage 7: 构建报告上下文
        report = self._build_report(
            gene_cards=gene_cards,
            tas2r38=tas2r38_result,
            mt_results=mt_results,
            or_tier_results=or_tier_results,
            enrichment_data=enrichment_data,
            key_snp_results=key_snp_results,
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

    def _merge_vep(self, variant: Variant, vep_data: Dict[str, Any]) -> Variant:
        """将 VEP 结果合并到 Variant 模型.

        对基因重叠区域的变异（如 NLRP3/OR2B11 共区、HERC2/OCA2 重叠），
        收集所有涉及的感官基因符号，确保变异能正确归属到每个相关基因。
        """
        if not vep_data:
            return variant

        consequences = vep_data.get("transcript_consequences", [])
        if consequences:
            # 收集该变异涉及的所有感官基因符号
            all_sensory_genes: set = set()
            for tc in consequences:
                gs = tc.get("gene_symbol", "")
                if gs and self.gene_sets.is_sensory_gene(gs):
                    all_sensory_genes.add(gs)

            # BED 映射补充：当 VEP 未返回某个感官基因（如 HERC2/OCA2 重叠区）
            bed_gene = self.bed_mapper.lookup(variant.chrom, variant.pos)
            if bed_gene and self.gene_sets.is_sensory_gene(bed_gene):
                if bed_gene not in all_sensory_genes:
                    all_sensory_genes.add(bed_gene)
                    logger.debug(
                        "BED-mapped gene added: %s:%d -> %s (VEP had %s)",
                        variant.chrom, variant.pos, bed_gene,
                        sorted(all_sensory_genes - {bed_gene}) or "none",
                    )

            # 如果涉及多个感官基因，存储到 gene_symbols
            if len(all_sensory_genes) > 1:
                variant.gene_symbols = sorted(all_sensory_genes)
                logger.debug(
                    "Multi-gene variant: %s:%d %s>%s -> %s",
                    variant.chrom, variant.pos, variant.ref, variant.alt,
                    variant.gene_symbols,
                )

            # 选出整体最优转录本（不偏袒任何基因，按严重度排名）
            best_tc = self._select_best_transcript(consequences, gene="")

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
                # 提取蛋白位置
                if best_tc.get("protein_start"):
                    variant.protein_position = best_tc.get("protein_start")
                # 提取 SIFT / PolyPhen
                if best_tc.get("sift"):
                    variant.sift = str(best_tc.get("sift"))
                if best_tc.get("polyphen"):
                    variant.polyphen = str(best_tc.get("polyphen"))
                # 提取氨基酸替换（兼容异常格式如 GGGDTRA*527X、GGG*/527X 等）
                if best_tc.get("amino_acids") and best_tc.get("protein_start"):
                    aa = str(best_tc.get("amino_acids"))
                    pos = str(best_tc.get("protein_start", ""))
                    if "/" in aa and aa != "/":
                        parts = aa.split("/", 1)
                        ref_aa = parts[0].strip()
                        alt_aa = parts[1].strip()
                        # 处理涉及终止密码子的异常格式
                        if "*" in alt_aa and len(alt_aa) > 3:
                            variant.amino_acid_change = f"{ref_aa}{pos}{alt_aa}"[:30]
                        else:
                            variant.amino_acid_change = f"{ref_aa}{pos}{alt_aa}"

        # 提取 colocated variants（gnomAD AF, ClinVar, rsID 等）
        colocated = vep_data.get("colocated_variants", [])
        for cv in colocated:
            # rsID — VEP 使用 "id" 字段
            cid = cv.get("id", "")
            if cid and str(cid).startswith("rs"):
                variant.rsid = str(cid)
            # LoF flags（直接从 cv 级提取）
            if cv.get("lof"):
                variant.lof_flags = str(cv.get("lof"))
            if cv.get("lof_flags"):
                variant.lof_flags = str(cv.get("lof_flags"))
            # gnomAD 频率 — 在 frequencies.{alt_allele} 下
            frequencies = cv.get("frequencies", {})
            if isinstance(frequencies, dict) and variant.alt:
                alt_allele_data = frequencies.get(variant.alt)
                if alt_allele_data and isinstance(alt_allele_data, dict):
                    # gnomAD genomes overall AF
                    gnomadg = alt_allele_data.get("gnomadg")
                    if gnomadg is not None:
                        variant.gnomad_af_genome = float(gnomadg)
                    # gnomAD exomes overall AF
                    gnomade = alt_allele_data.get("gnomade")
                    if gnomade is not None:
                        variant.gnomad_af_exome = float(gnomade)
                    # 综合 AF 取两者中的较大值
                    g_vals = [v for v in [gnomadg, gnomade] if v is not None]
                    if g_vals:
                        variant.af_gnomad = max(float(v) for v in g_vals)

        variant.raw_vep = vep_data
        return variant

    def _assess_genes(self, gene_groups: Dict[str, List[Variant]]) -> List[GeneCard]:
        """对每个基因进行功能影响评估，并筛选关键功能变异."""
        gene_cards = []
        for gene_symbol, variants in gene_groups.items():
            if not variants:
                continue

            # 评估所有变异
            assessments = [
                self.impact_engine.assess(v, self.config.sex, self.gene_sets)
                for v in variants
            ]
            best_assessment = max(
                assessments,
                key=lambda a: self._level_rank(a.level),
            )

            # 筛选关键功能变异（基因维度分析：只保留影响蛋白功能的）
            key_variants = [v for v in variants if self._is_key_variant(v)]
            # 标记关键变异
            for v in key_variants:
                v.is_key_variant = True

            # 如果没有关键变异但基因有影响评估，保留最高影响的一个用于展示
            if not key_variants and best_assessment.level in (
                "完全丧失", "显著影响", "部分影响", "可能轻微影响"
            ):
                top_var = max(variants, key=lambda v: self._variant_impact_score(v))
                top_var.is_key_variant = True
                key_variants = [top_var]

            # 构建蛋白影响综合摘要
            protein_summary = self._build_protein_impact_summary(gene_symbol, key_variants)

            gene_card = GeneCard(
                gene_symbol=gene_symbol,
                subsystem=self.gene_sets.get_subsystem(gene_symbol),
                sensory_function_zh=self.gene_sets.get_gene_function(gene_symbol),
                variants=variants,
                key_variants=key_variants,
                assessment=best_assessment,
                protein_impact_summary=protein_summary,
            )
            gene_cards.append(gene_card)

        return gene_cards

    @staticmethod
    def _is_key_variant(variant: Variant) -> bool:
        """判断变异是否为关键功能变异.

        关键变异标准：
        - 必须是高后果类型（frameshift, stop_gained, splice, missense 等）
        - 排除同义、内含子、UTR、调控区、非编码等低影响类型
        - missense 额外要求：有害预测或功能域注释或纯合
        - 个人特征基因（pigmentation/metabolism/muscle/hair）的纯合调控区变异也保留，
          因为这类变异可能是功能调控位点（如 HERC2 rs12913832 增强子变异）
        """
        cons = (variant.consequence or "").lower()
        gene = variant.gene_symbol or ""

        # 个人特征基因列表：调控区变异可能具有功能意义
        trait_genes = {
            "HERC2", "SLC45A2", "SLC24A4", "SLC24A5", "IRF4",
            "CYP1A2", "LCT", "ALDH2", "ADH1B",
            "ACTN3", "EDAR", "OR6A2",
        }
        is_trait_gene = gene in trait_genes

        # 对个人特征基因，纯合的调控区/UTR 变异也视为关键变异
        if is_trait_gene and variant.is_homozygous:
            if any(t in cons for t in [
                "upstream_gene_variant", "downstream_gene_variant",
                "regulatory_region_variant", "tf_binding_site_variant",
                "5_prime_utr_variant", "3_prime_utr_variant",
            ]):
                return True

        # 排除的低影响类型
        low_impact_keywords = [
            "synonymous", "intron", "utr", "upstream", "downstream",
            "intergenic", "regulatory_region", "non_coding", "nmd_transcript",
            "feature_elongation", "feature_truncation", "mature_mirna",
        ]
        for low in low_impact_keywords:
            if low in cons:
                return False

        # 必须是功能相关的后果类型
        high_impact_keywords = [
            "frameshift", "stop_gained", "stop_lost", "start_lost",
            "splice", "missense", "inframe", "protein_altering",
            "transcript_ablation", "transcript_amplification",
        ]
        has_high = any(hi in cons for hi in high_impact_keywords)
        if not has_high:
            return False

        # missense 额外过滤：需要有害预测、功能域注释、或纯合
        if "missense" in cons:
            has_damaging = (
                (variant.sift and "deleterious" in variant.sift.lower())
                or (variant.polyphen and "damaging" in variant.polyphen.lower())
                or bool(variant.protein_domain)
                or variant.is_homozygous
                or (variant.af_gnomad is not None and variant.af_gnomad < 0.01)
            )
            if not has_damaging:
                return False

        return True

    @staticmethod
    def _variant_impact_score(variant: Variant) -> int:
        """计算变异的蛋白影响分数（用于排序）."""
        score = 0
        cons = (variant.consequence or "").lower()

        consequence_scores = {
            "frameshift": 100, "stop_gained": 95, "stop_lost": 90,
            "start_lost": 85, "splice": 80, "missense": 50,
            "inframe": 40, "protein_altering": 35,
        }
        for key, val in consequence_scores.items():
            if key in cons:
                score = max(score, val)

        if variant.is_homozygous:
            score += 30
        if variant.sift and "deleterious" in variant.sift.lower():
            score += 20
        if variant.polyphen and "probably_damaging" in variant.polyphen.lower():
            score += 15
        if variant.polyphen and "possibly_damaging" in variant.polyphen.lower():
            score += 10
        if variant.protein_domain:
            score += 5
        score += min(int(variant.qual / 10), 10)
        return score

    @staticmethod
    def _build_protein_impact_summary(gene_symbol: str, key_variants: List[Variant]) -> str:
        """基于关键变异构建综合蛋白影响摘要."""
        if not key_variants:
            return "未发现影响蛋白功能的关键变异。"

        lof_types = ["frameshift", "stop_gained", "stop_lost", "start_lost", "splice"]
        lof_variants = [v for v in key_variants if any(t in (v.consequence or "").lower() for t in lof_types)]
        missense_variants = [v for v in key_variants if "missense" in (v.consequence or "").lower()]
        damaging_missense = [v for v in missense_variants
            if (v.sift and "deleterious" in v.sift.lower())
            or (v.polyphen and "damaging" in v.polyphen.lower())]
        hom_variants = [v for v in key_variants if v.is_homozygous]

        parts = []
        if lof_variants:
            lof_info = f"发现 {len(lof_variants)} 个功能丧失型变异"
            if hom_variants:
                hom_lof = [v for v in lof_variants if v.is_homozygous]
                if hom_lof:
                    lof_info += f"（含 {len(hom_lof)} 个纯合）"
            parts.append(lof_info)

        if missense_variants:
            mis_info = f"发现 {len(missense_variants)} 个错义变异"
            if damaging_missense:
                mis_info += f"，其中 {len(damaging_missense)} 个被预测为可能有害"
            if hom_variants:
                hom_mis = [v for v in missense_variants if v.is_homozygous]
                if hom_mis:
                    mis_info += f"（含 {len(hom_mis)} 个纯合）"
            parts.append(mis_info)

        if not lof_variants and not missense_variants:
            parts.append(f"发现 {len(key_variants)} 个蛋白结构改变变异")

        summary = "；".join(parts) + "。"

        # 添加遗传模式推断
        if hom_variants:
            summary += "存在纯合变异，若为隐性遗传模式则可能显著影响蛋白功能。"
        elif len(key_variants) >= 2:
            summary += "多个杂合变异共存，需考虑复合杂合的可能性。"
        else:
            summary += "杂合状态下，野生型等位基因通常可维持基本功能。"

        return summary

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

    def _apply_gnomad_af_downgrade(
        self,
        gene_cards: List[GeneCard],
        enrichment_data: Dict[str, Dict[str, Any]],
    ) -> int:
        """用 gnomAD AF 对所有基因进行可能性评估降级.

        对 enrichment 数据中 gnomAD AF > 30% 的纯合蛋白影响变异，
        无论基因当前评估级别如何，都降级为 "无影响"。
        仅影响具有纯合蛋白影响变异的基因。

        Returns:
            被降级的基因数量。
        """
        protein_csq = {"frameshift_variant", "stop_gained",
                       "splice_acceptor_variant", "splice_donor_variant",
                       "missense_variant", "inframe_deletion", "inframe_insertion"}
        downgraded = 0

        for card in gene_cards:
            ed = enrichment_data.get(card.gene_symbol, {})
            gnmd_vars = ed.get("gnomad_variants", [])
            if not gnmd_vars:
                continue

            # Find homozygous protein-affecting variants with AF > 30%
            has_high_af = False
            af_values = []
            for v in card.variants:
                if not v.is_homozygous:
                    continue
                if v.consequence not in protein_csq:
                    continue
                # Match to gnomAD data
                for gv in gnmd_vars:
                    vi = gv.get("variant", {})
                    if vi.get("pos") == v.pos and vi.get("ref") == v.ref and vi.get("alt") == v.alt:
                        af = gv.get("result", {}).get("gnomad_af")
                        if af is not None:
                            af_values.append(af)
                            if af > 0.30:
                                has_high_af = True
                        break

            if not has_high_af:
                continue

            # All homozygous protein-affecting variants have AF > 30% → common polymorphism
            # Downgrade assessment to "无影响"
            old_level = card.assessment.level
            old_level_rank = self._level_rank(old_level)
            if old_level_rank <= self._level_rank("无影响"):
                continue  # Already minimal

            max_af = max(af_values) * 100 if af_values else 0
            card.assessment.level = "无影响"
            card.assessment.rationale_zh = (
                f"gnomAD 频率校验：检出纯合蛋白影响变异但人群频率 {max_af:.0f}% (常见多态)，"
                f"参考基因组可能携带罕见等位基因。原评估: {old_level}"
            )
            downgraded += 1
            logger.debug("AF downgrade: %s %s → 无影响 (max AF=%.1f%%)",
                        card.gene_symbol, old_level, max_af)

        return downgraded

    def _reclassify_or_genes_with_enrichment(
        self,
        gene_cards: List[GeneCard],
        enrichment_data: Dict[str, Dict[str, Any]],
        fallback_tiers: Optional[List[Any]] = None,
    ) -> Optional[List[Any]]:
        """使用 enrichment 数据重新分类 OR 基因（gnomAD AF + ClinVar）."""
        or_variants = []
        for card in gene_cards:
            if card.gene_symbol.startswith("OR"):
                or_variants.extend(card.variants)

        if not or_variants:
            return None

        try:
            return self.or_classifier.classify_with_enrichment(
                or_variants, self.config.sex, self.gene_sets, enrichment_data
            )
        except Exception as exc:
            logger.warning("OR enrichment-aware reclassification failed: %s, using fallback", exc)
            return fallback_tiers

    # 预计算的影响程度排序字典（用于按需 API 富集阈值判断）
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

    async def _enrich_genes(
        self, gene_cards: List[GeneCard]
    ) -> Dict[str, Dict[str, Any]]:
        """异步并发查询外部 API 富集基因信息.

        双路径策略:
            Path A (primary): 所有子系统评估级别 >= "部分影响" 的基因
                → 基因级: UniProt + gnomAD(gene) + ClinVar + GTEx
            Path B (secondary): OR 基因中检出纯合蛋白影响变异的基因
                → 基因级: UniProt + ClinVar(gene) + GTEx
                → 变异级: gnomAD(per variant) → enrichment[gene]["gnomad_variants"]

        两条路径去重后统一并发查询。
        """
        if not self.config.show_reference_info:
            return {}

        # ── Path A: 影响评估阈值筛选 ──
        min_rank = self._LEVEL_RANK.get("部分影响", 2)
        impactful_cards = [
            card for card in gene_cards
            if card.gene_symbol and self._level_rank(card.assessment.level) >= min_rank
        ]
        path_a_genes = {card.gene_symbol for card in impactful_cards}

        # ── Path B: OR 基因纯合蛋白影响变异 ──
        or_protein_csq = {"frameshift_variant", "stop_gained",
                          "splice_acceptor_variant", "splice_donor_variant",
                          "missense_variant", "inframe_deletion", "inframe_insertion"}
        path_b_genes: Dict[str, list] = {}  # gene → list of homozygous protein-affecting variants
        for card in gene_cards:
            gene = card.gene_symbol
            if not gene or not gene.startswith("OR"):
                continue
            if gene in path_a_genes:
                continue
            homo_variants = [
                v for v in card.variants
                if v.gt == "1/1" and v.consequence in or_protein_csq
            ]
            if homo_variants:
                path_b_genes[gene] = homo_variants

        all_genes = sorted(set(path_a_genes) | set(path_b_genes.keys()))
        if not all_genes:
            logger.info("Stage 6: no genes match enrichment criteria, skipping")
            return {}

        logger.info(
            "Stage 6: enriching %d genes (Path A: %d impactful, Path B: %d OR-homozygous)",
            len(all_genes), len(path_a_genes), len(path_b_genes),
        )

        # ── 逐个基因顺序执行，带详细日志 ──
        enrichment: Dict[str, Dict[str, Any]] = {}
        for idx, gene in enumerate(all_genes):
            logger.info("Stage 6: [%d/%d] enriching %s", idx + 1, len(all_genes), gene)
            if gene not in enrichment:
                enrichment[gene] = {}

            # UniProt
            try:
                logger.debug("  -> UniProt query for %s", gene)
                enrichment[gene]["uniprot"] = await self._safe_query(self.uniprot_client.query, gene)
                logger.debug("  <- UniProt done for %s", gene)
            except Exception as exc:
                logger.error("UniProt error for %s: %s", gene, exc)
                enrichment[gene]["uniprot"] = {"error": str(exc)}

            # ClinVar
            try:
                logger.debug("  -> ClinVar query for %s", gene)
                enrichment[gene]["clinvar"] = await self._safe_query(self.clinvar_client.query, gene)
                logger.debug("  <- ClinVar done for %s", gene)
            except Exception as exc:
                logger.error("ClinVar error for %s: %s", gene, exc)
                enrichment[gene]["clinvar"] = {"error": str(exc)}

            # GTEx
            try:
                logger.debug("  -> GTEx query for %s", gene)
                enrichment[gene]["gtex"] = await self._safe_query(self.gtex_client.query, gene)
                logger.debug("  <- GTEx done for %s", gene)
            except Exception as exc:
                logger.error("GTEx error for %s: %s", gene, exc)
                enrichment[gene]["gtex"] = {"error": str(exc)}

            # gnomAD
            if gene in path_a_genes:
                try:
                    logger.debug("  -> gnomAD gene query for %s", gene)
                    enrichment[gene]["gnomad"] = await self._safe_query(self.gnomad_client.query, gene)
                    logger.debug("  <- gnomAD done for %s", gene)
                except Exception as exc:
                    logger.error("gnomAD error for %s: %s", gene, exc)
                    enrichment[gene]["gnomad"] = {"error": str(exc)}
                enrichment[gene]["gnomad_variants"] = []
            elif gene in path_b_genes:
                variants_list = enrichment[gene].setdefault("gnomad_variants", [])
                for v in path_b_genes[gene]:
                    try:
                        logger.debug("  -> gnomAD variant query for %s %s:%d", gene, v.chrom, v.pos)
                        result = await self._safe_query_variant(
                            self.local_gnomad.query_variant_async,
                            str(v.chrom), v.pos, v.ref, v.alt,
                        )
                        variants_list.append({"variant": {"pos": v.pos, "ref": v.ref, "alt": v.alt, "consequence": v.consequence}, "result": result})
                        logger.debug("  <- gnomAD variant done for %s", gene)
                    except Exception as exc:
                        logger.error("gnomAD variant error for %s: %s", gene, exc)
                        variants_list.append({"variant": {"pos": v.pos, "ref": v.ref, "alt": v.alt, "consequence": v.consequence}, "error": str(exc)})

        logger.info("Stage 6: enrichment loop completed for %d genes", len(enrichment))
        return enrichment

        # ── 按基因归类 ──
        enrichment: Dict[str, Dict[str, Any]] = {}
        for (gene, api_name, variant_info), result in zip(task_meta, results):
            if gene not in enrichment:
                enrichment[gene] = {}

            if api_name == "gnomad_variant":
                # 变异性 gnomAD → 聚合到 gnomad_variants 列表
                variants_list = enrichment[gene].setdefault("gnomad_variants", [])
                entry = {
                    "variant": variant_info,
                }
                if isinstance(result, Exception):
                    entry["error"] = str(result)
                else:
                    entry["result"] = result
                variants_list.append(entry)
            else:
                if isinstance(result, Exception):
                    enrichment[gene][api_name] = {"error": str(result)}
                else:
                    enrichment[gene][api_name] = result

        # 同样为 Path A 基因补充一个空的 gnomad_variants（如果没有的话）
        # 维持 report 模板的兼容性
        for gene in path_a_genes:
            if gene in enrichment and "gnomad_variants" not in enrichment[gene]:
                enrichment[gene]["gnomad_variants"] = []

        return enrichment

    @staticmethod
    async def _safe_query(query_func, key: str) -> Any:
        """安全查询，异常时返回降级结果."""
        try:
            return await query_func(key)
        except Exception as exc:
            return {"error": str(exc), "api": "unknown", "key": key}

    @staticmethod
    async def _safe_query_variant(query_func, *args, **kwargs) -> Any:
        """安全查询变异级别的 API，异常时返回降级结果."""
        try:
            return await query_func(*args, **kwargs)
        except Exception as exc:
            return {"error": str(exc)}

    def _build_report(
        self,
        gene_cards: List[GeneCard],
        tas2r38: Optional[Any],
        mt_results: Optional[List[Any]],
        or_tier_results: Optional[List[Any]],
        enrichment_data: Dict[str, Dict[str, Any]],
        key_snp_results: Optional[List[Any]] = None,
    ) -> SensoryReport:
        """构建完整报告."""
        # 注入富集数据
        for card in gene_cards:
            gene = card.gene_symbol
            if gene in enrichment_data:
                card.enrichment_data = enrichment_data[gene]

        # 构建执行摘要
        executive_summary = self._build_executive_summary(gene_cards, key_snp_results)

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
            key_snps=key_snp_results,
            executive_summary=executive_summary,
            data_availability=data_availability,
            disclaimer_zh=self._default_disclaimer(),
        )
        return report

    @staticmethod
    def _build_executive_summary(
        gene_cards: List[GeneCard],
        key_snp_results: Optional[List[Any]] = None,
    ) -> ExecutiveSummary:
        """构建执行摘要统计（含个人特征定性预测）."""
        subsystem_counts: Dict[str, Dict[str, int]] = {}
        key_findings = []

        for card in gene_cards:
            level = card.assessment.level
            if level in ("完全丧失", "显著影响"):
                key_findings.append(
                    f"{card.gene_symbol}: {level}（{card.assessment.rationale_zh[:50]}...）"
                )

        # 生成个人特征定性预测
        from src.trait_predictor import PersonalTraitPredictor
        predictor = PersonalTraitPredictor(
            key_snp_results=key_snp_results,
            gene_cards=gene_cards,
        )
        personal_traits = predictor.predict_all()

        return ExecutiveSummary(
            subsystem_counts=subsystem_counts,
            key_findings=key_findings,
            personal_traits=personal_traits,
        )

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
        """释放资源."""
        await self.vep_client.close()
        await self.uniprot_client.close()
        await self.gnomad_client.close()
        await self.clinvar_client.close()
        await self.gtex_client.close()


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


def _resolve_vcf_with_index(vcf_path: str) -> str:
    """确保 VCF 有可用的 tabix 索引，必要时通过多重策略解决路径问题.

    策略优先级:
        1. 索引已存在 → 直接返回原路径
        2. 在原始路径创建索引
        3. 创建 /tmp symlink（绕过路径编码问题）→ 在 symlink 上索引
        4. 复制到 /tmp（绕过只读/iCloud 问题）→ 在副本上索引

    Returns:
        带可用索引的 VCF 路径（可能指向 tmp 副本）。
    Raises:
        RuntimeError: 所有策略均失败。
    """
    idx_path = vcf_path + ".tbi"
    idx_path_csi = vcf_path + ".csi"

    # 策略 1: 索引已存在
    if os.path.exists(idx_path) or os.path.exists(idx_path_csi):
        return vcf_path

    # 策略 2: 在原始路径直接索引
    logger.info("VCF index not found, attempting direct indexing...")
    try:
        subprocess.run(
            ["bcftools", "index", "-t", vcf_path],
            stderr=subprocess.PIPE, check=True, timeout=120,
        )
        logger.info("VCF index created at original path: %s", idx_path)
        return vcf_path
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        logger.warning("Direct indexing failed (%s), trying symlink strategy...", exc)

    # 策略 3: 创建 /tmp symlink（绕过路径中的特殊字符）
    tmp_symlink = os.path.join(tempfile.gettempdir(), f"sensory_vcf_{os.getpid()}.vcf.gz")
    try:
        if os.path.exists(tmp_symlink) or os.path.islink(tmp_symlink):
            os.unlink(tmp_symlink)
        os.symlink(os.path.abspath(vcf_path), tmp_symlink)
        logger.info("Created symlink: %s -> %s", tmp_symlink, vcf_path)
        subprocess.run(
            ["bcftools", "index", "-t", tmp_symlink],
            stderr=subprocess.PIPE, check=True, timeout=120,
        )
        logger.info("VCF index created via symlink: %s.tbi", tmp_symlink)
        return tmp_symlink
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        logger.warning("Symlink strategy failed (%s), trying copy strategy...", exc)
        if os.path.exists(tmp_symlink) or os.path.islink(tmp_symlink):
            os.unlink(tmp_symlink)

    # 策略 4: 复制到 /tmp（最终 fallback）
    tmp_copy = os.path.join(tempfile.gettempdir(), f"sensory_vcf_{os.getpid()}_copy.vcf.gz")
    try:
        logger.info("Copying VCF to %s ...", tmp_copy)
        subprocess.run(
            ["cp", vcf_path, tmp_copy],
            stderr=subprocess.PIPE, check=True, timeout=300,
        )
        subprocess.run(
            ["bcftools", "index", "-t", tmp_copy],
            stderr=subprocess.PIPE, check=True, timeout=120,
        )
        logger.info("VCF copied and indexed at: %s", tmp_copy)
        return tmp_copy
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        logger.error("Copy strategy also failed: %s", exc)
        if os.path.exists(tmp_copy):
            os.unlink(tmp_copy)
        raise RuntimeError(
            f"Failed to create VCF index via all strategies. "
            f"Original path: {vcf_path}. "
            f"Please check file permissions and bcftools installation."
        )


def _filter_vcf_by_bed(vcf_path: str, strict_filter: bool = False) -> tuple[str, Optional[str]]:
    """用 BED 文件筛选感官基因区域的变异，返回临时 VCF 路径和使用的 BED 路径.

    索引策略：多重 fallback（直接索引 → symlink → 复制到 /tmp），
    只有所有策略都失败才放弃 BED 过滤。不再因路径编码或权限问题静默跳过。

    Args:
        vcf_path: 输入 VCF 路径。
        strict_filter: 若 True，优先使用精确外显子/CDS BED 做硬过滤。

    Returns:
        (filtered_vcf_path, bed_path_used)
    """
    script_dir = Path(__file__).resolve().parent
    work_dir = Path(__file__).resolve().parent.parent
    skill_root = Path(__file__).resolve().parent.parent.parent

    # 精确坐标 BED（外显子/CDS）—— strict 模式优先
    precise_beds = []
    if strict_filter:
        precise_beds = [
            script_dir.parent / "assets" / "data" / "sensory_genes_exons.bed",
            work_dir / "sensory_genes_exons.bed",
            work_dir / "assets" / "data" / "sensory_genes_exons.bed",
            skill_root / "assets" / "data" / "sensory_genes_exons.bed",
            Path("/tmp/sensory_precompute/sensory_genes_exons.bed"),
        ]

    # 基因全区域 BED —— fallback
    broad_beds = [
        script_dir.parent / "assets" / "data" / "sensory_gene_regions.bed",
        work_dir / "sensory_genes.bed",
        work_dir / "assets" / "data" / "sensory_gene_regions.bed",
        skill_root / "assets" / "data" / "sensory_gene_regions.bed",
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
        for bed in broad_beds:
            if bed.exists():
                bed_path = str(bed)
                bed_source = "broad (gene region)"
                break

    if bed_path is None:
        # 彻底没有 BED — 这是 hard error，不能静默跳过
        raise RuntimeError(
            "No sensory gene BED file found. "
            "Please run precompute first: python scripts/precompute.py"
        )

    # 多重策略确保 VCF 索引可用
    try:
        indexed_vcf = _resolve_vcf_with_index(vcf_path)
    except RuntimeError as exc:
        logger.critical("BED filter failed — all indexing strategies exhausted: %s", exc)
        raise  # 不再静默跳过，直接报错让用户知道

    # 检查 VCF 染色体命名
    has_chr = False
    try:
        result = subprocess.run(
            ["bcftools", "view", indexed_vcf],
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
            ["bcftools", "view", "-Oz", "-o", tmp_vcf.name, "-R", bed_path, indexed_vcf],
            stderr=subprocess.PIPE,
            check=True,
            timeout=120,
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
        logger.info("BED-filtered VCF: %s -> %d variants (%s)", tmp_vcf.name, variant_count, bed_source)
        if variant_count == 0:
            logger.warning("No variants in sensory gene regions (0 variants after BED filter)")
        return tmp_vcf.name, bed_path
    except subprocess.CalledProcessError as exc:
        logger.error("BED filter failed: %s", exc)
        if os.path.exists(tmp_vcf.name):
            os.unlink(tmp_vcf.name)
        raise RuntimeError(f"bcftools view -R failed: {exc}") from exc
    except subprocess.TimeoutExpired as exc:
        logger.error("BED filter timeout")
        if os.path.exists(tmp_vcf.name):
            os.unlink(tmp_vcf.name)
        raise RuntimeError("bcftools view -R timed out") from exc


def _auto_infer_sex(vcf_path: str) -> Optional[str]:
    """自动推断样本性别（从 chrX/chrY 基因型模式）."""
    logger.info("Auto-detecting sample sex from VCF...")
    sex = detect_sex(vcf_path)
    if sex:
        logger.info("Inferred sex: %s (%s)", sex, "男性" if sex == "M" else "女性")
    else:
        logger.warning("Could not determine sex from VCF (chrX/chrY may be absent)")
    return sex


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

    # 保存原始 VCF 路径（BED 过滤前），供 KeySNPInferrer 使用
    original_vcf_path = config.vcf_path
    config.original_vcf_path = original_vcf_path

    # 自动推断性别（当用户使用 --auto-sex 未提供 --sex 时）
    if hasattr(config, 'auto_sex') and config.auto_sex:
        inferred = _auto_infer_sex(original_vcf_path)
        if inferred:
            config.sex = inferred
        else:
            logger.error("Cannot auto-detect sex. Please provide --sex M|F explicitly.")
            sys.exit(1)

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
    parser.add_argument("--sex", default=None, choices=["M", "F"],
                        help="Sample sex (M/F). Use --auto-sex to auto-detect from VCF.")
    parser.add_argument("--auto-sex", action="store_true",
                        help="Auto-detect sex from chrX/chrY genotype patterns")
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

    # Sex resolution: --auto-sex or --sex required
    sex = args.sex
    auto_sex = args.auto_sex
    if not sex and not auto_sex:
        parser.error("Either --sex M|F or --auto-sex is required")

    return AnalysisConfig(
        vcf_path=args.vcf,
        sex=sex or "M",  # placeholder, overwritten by auto_sex if enabled
        subsystems=args.subsystems.split(","),
        output_dir=args.output_dir,
        known_phenotype=args.known_phenotype,
        show_reference_info=not args.no_reference_info,
        strict_filter=args.strict_filter,
        auto_sex=auto_sex,
    )


if __name__ == "__main__":
    asyncio.run(main())
