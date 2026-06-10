"""报告上下文构建器.

汇总核心管线、专用模块和 API 富集数据，构建 SensoryReport 对象。
"""

from typing import Any, Dict, List, Optional

from src.gene_sets.loader import GeneSetLoader
from src.logger import get_logger
from src.models import (
    DataAvailability,
    ExecutiveSummary,
    GeneCard,
    KeySNPResult,
    MitochondrialResult,
    ORTierResult,
    PersonalTraitPrediction,
    SensoryReport,
    TAS2R38Result,
)
from src.specialized.mitochondrial import MitochondrialAnnotator
from src.specialized.or_tiers import ORTierClassifier
from src.specialized.tas2r38 import TAS2R38Analyzer
from src.trait_predictor import PersonalTraitPredictor

logger = get_logger(__name__)


class ReportContextBuilder:
    """报告上下文构建器."""

    def __init__(self, gene_sets: GeneSetLoader) -> None:
        """初始化构建器.

        Args:
            gene_sets: 感官基因集加载器。
        """
        self.gene_sets = gene_sets
        self.tas2r38_analyzer = TAS2R38Analyzer()
        self.mt_annotator = MitochondrialAnnotator()
        self.or_classifier = ORTierClassifier()

    def build(
        self,
        gene_cards: List[GeneCard],
        enrichment_data: Optional[Dict[str, Dict[str, Any]]] = None,
        key_snp_results: Optional[List[KeySNPResult]] = None,
    ) -> SensoryReport:
        """构建完整的报告上下文.

        Args:
            gene_cards: 全部基因卡片。
            enrichment_data: API 富集数据（可选）。

        Returns:
            SensoryReport 报告对象。
        """
        if enrichment_data is None:
            enrichment_data = {}

        # 1. 执行专用逻辑分析
        tas2r38_result = self._analyze_tas2r38(gene_cards)
        mt_results = self._analyze_mitochondrial(gene_cards)
        or_tier_results = self._classify_or_genes(gene_cards)

        # 2. 将富集数据注入基因卡片
        for card in gene_cards:
            gene = card.gene_symbol
            if gene in enrichment_data:
                card.enrichment_data = enrichment_data[gene]

        # 3. 构建执行摘要（包含个人特征预测）
        executive_summary = self._build_executive_summary(
            gene_cards, key_snp_results=key_snp_results
        )

        # 4. 构建数据可用性总览
        data_availability = {}
        for card in gene_cards:
            data_availability[card.gene_symbol] = self._build_data_availability(
                card.gene_symbol, enrichment_data
            )

        report = SensoryReport(
            gene_cards=gene_cards,
            tas2r38=tas2r38_result,
            mitochondrial=mt_results if mt_results else None,
            or_tiers=or_tier_results if or_tier_results else None,
            key_snps=key_snp_results if key_snp_results else None,
            executive_summary=executive_summary,
            data_availability=data_availability,
            disclaimer_zh=self._default_disclaimer(),
        )

        return report

    def _analyze_tas2r38(self, gene_cards: List[GeneCard]) -> Optional[TAS2R38Result]:
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
    ) -> Optional[List[MitochondrialResult]]:
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

    def _classify_or_genes(self, gene_cards: List[GeneCard]) -> Optional[List[ORTierResult]]:
        """对 OR 基因进行分级分类."""
        or_variants = []
        for card in gene_cards:
            if card.gene_symbol.startswith("OR"):
                or_variants.extend(card.variants)

        if not or_variants:
            return None

        try:
            # 使用默认 sex，实际应由外部传入
            return self.or_classifier.classify(or_variants, "M", self.gene_sets)
        except Exception as exc:
            logger.error("OR tier classification failed: %s", exc)
            return None

    @staticmethod
    def _build_executive_summary(
        gene_cards: List[GeneCard],
        key_snp_results: Optional[List[KeySNPResult]] = None,
    ) -> ExecutiveSummary:
        """构建执行摘要（含个人特征定性预测）."""
        subsystem_counts: Dict[str, Dict[str, int]] = {}
        key_findings = []

        for card in gene_cards:
            subsystem = card.subsystem or "unknown"
            level = card.assessment.level
            if level in ("完全丧失", "显著影响"):
                key_findings.append(
                    f"{card.gene_symbol}: {level}（{card.assessment.rationale_zh[:50]}...）"
                )
            # 统计各子系统 × 影响级别计数
            if subsystem not in subsystem_counts:
                subsystem_counts[subsystem] = {}
            subsystem_counts[subsystem][level] = subsystem_counts[subsystem].get(level, 0) + 1

        # 生成个人特征定性预测
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
    def _build_data_availability(
        gene: str, enrichment_data: Dict[str, Dict[str, Any]]
    ) -> DataAvailability:
        """构建数据可用性状态."""
        data = enrichment_data.get(gene, {})
        return DataAvailability(
            gnomad_af="可用" if data.get("gnomad", {}).get("found") else "N/A",
            clinvar="可用" if data.get("clinvar", {}).get("found") else "N/A",
            spliceai="N/A",
            cadd="N/A",
            topology="可用" if data.get("uniprot", {}).get("found") else "N/A",
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
