"""端到端冒烟测试.

验证管线骨架能正确组装、主要模块可导入、核心流程可执行（使用 mock 数据）。
"""

import json
import sys
from pathlib import Path

import pytest

# 确保测试覆盖生产入口 scripts/main.py
_SCRIPTS_DIR = Path(__file__).parent.parent / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from src.assessment.engine import ImpactEngine
from src.config_loader import load_config
from src.gene_sets.filter import SensoryFilter
from src.gene_sets.loader import GeneSetLoader
from src.models import AnalysisConfig, GeneCard, ImpactAssessment, SensoryReport, Variant
from src.report.json_generator import JsonReportGenerator
from src.report.markdown_generator import MarkdownReportGenerator
from src.report.report_context import ReportContextBuilder
from src.specialized.mitochondrial import MitochondrialAnnotator
from src.specialized.or_tiers import ORTierClassifier
from src.specialized.tas2r38 import TAS2R38Analyzer
from src.vcf.prefilter import Prefilter


class TestImports:
    """模块导入测试."""

    def test_all_main_modules_importable(self) -> None:
        """验证所有核心模块可导入."""
        import src.models
        import src.config_loader
        import src.gene_sets.loader
        import src.gene_sets.filter
        import src.assessment.engine
        import src.assessment.rules
        import src.assessment.inheritance
        import src.specialized.tas2r38
        import src.specialized.mitochondrial
        import src.specialized.or_tiers
        import src.report.markdown_generator
        import src.report.json_generator
        import src.report.report_context
        import src.exceptions
        import src.logger

    def test_main_module_importable(self) -> None:
        """生产入口 scripts/main.py 应可导入."""
        from main import SensoryPipeline, run_analysis
        assert SensoryPipeline is not None
        assert run_analysis is not None


class TestPipelineAssembly:
    """管线组装冒烟测试."""

    def test_gene_set_loading(self) -> None:
        loader = GeneSetLoader()
        assert len(loader.get_all_genes()) > 0

    def test_config_loading(self) -> None:
        cfg = load_config()
        assert cfg.vep.source in ("rest_api", "local", "local_docker")

    def test_prefilter_creation(self) -> None:
        p = Prefilter(min_qual=30, min_dp=10, pass_only=True)
        assert p.min_qual == 30

    def test_impact_engine_creation(self) -> None:
        engine = ImpactEngine()
        assert len(engine.rules) == 2
        assert engine.inheritance_matcher is not None

    def test_specialized_modules_creation(self) -> None:
        tas2r38 = TAS2R38Analyzer()
        mt = MitochondrialAnnotator()
        or_cls = ORTierClassifier()
        assert tas2r38 is not None
        assert mt is not None
        assert or_cls is not None

    def test_report_generators_creation(self) -> None:
        md_gen = MarkdownReportGenerator()
        json_gen = JsonReportGenerator()
        assert md_gen is not None
        assert json_gen is not None


class TestMockPipeline:
    """使用 mock 数据模拟完整分析流程."""

    @pytest.fixture
    def gene_sets(self) -> GeneSetLoader:
        return GeneSetLoader()

    @pytest.fixture
    def mock_variants(self) -> list[Variant]:
        return [
            Variant(chrom="7", pos=141972755, ref="G", alt="C", gt="1/1", gene_symbol="TAS2R38", consequence="missense_variant"),
            Variant(chrom="13", pos=100, ref="A", alt="T", gt="1/1", gene_symbol="GJB2", consequence="frameshift_variant"),
            Variant(chrom="1", pos=100, ref="A", alt="T", gt="1/1", gene_symbol="OR2T11", consequence="frameshift_variant"),
            Variant(chrom="MT", pos=1555, ref="A", alt="G", gt="1", gene_symbol="MT-RNR1", consequence="synonymous_variant"),
        ]

    def test_sensory_filter(self, mock_variants: list[Variant], gene_sets: GeneSetLoader) -> None:
        f = SensoryFilter(gene_sets)
        filtered = f.filter_variants(mock_variants)
        assert len(filtered) == 4
        genes = [v.gene_symbol for v in filtered]
        assert "TAS2R38" in genes
        assert "GJB2" in genes

    def test_group_by_gene(self, mock_variants: list[Variant], gene_sets: GeneSetLoader) -> None:
        f = SensoryFilter(gene_sets)
        groups = f.group_by_gene(mock_variants)
        assert "TAS2R38" in groups
        assert "GJB2" in groups
        assert "OR2T11" in groups

    def test_impact_assessment(self, mock_variants: list[Variant], gene_sets: GeneSetLoader) -> None:
        engine = ImpactEngine()
        for v in mock_variants:
            result = engine.assess(v, "M", gene_sets)
            assert result.level in ("完全丧失", "显著影响", "部分影响", "可能轻微影响", "无影响")
            assert result.rationale_zh != ""

    def test_gene_card_building(self, mock_variants: list[Variant], gene_sets: GeneSetLoader) -> None:
        engine = ImpactEngine()
        f = SensoryFilter(gene_sets)
        groups = f.group_by_gene(mock_variants)

        gene_cards = []
        for gene_symbol, variants in groups.items():
            assessments = [engine.assess(v, "M", gene_sets) for v in variants]
            best = max(assessments, key=lambda a: {"无影响": 0, "可能轻微影响": 1, "部分影响": 2, "显著影响": 3, "完全丧失": 4}.get(a.level, 0))
            card = GeneCard(
                gene_symbol=gene_symbol,
                sensory_function_zh=gene_sets.get_gene_function(gene_symbol),
                variants=variants,
                assessment=best,
            )
            gene_cards.append(card)

        assert len(gene_cards) == 4

    def test_tas2r38_analysis(self, mock_variants: list[Variant]) -> None:
        analyzer = TAS2R38Analyzer()
        tas2r38_vars = [v for v in mock_variants if v.gene_symbol == "TAS2R38"]
        result = analyzer.analyze(tas2r38_vars)
        assert result is not None
        assert result.diplotype != ""

    def test_mitochondrial_analysis(self, mock_variants: list[Variant]) -> None:
        annotator = MitochondrialAnnotator()
        mt_vars = [v for v in mock_variants if v.gene_symbol.startswith("MT-")]
        results = annotator.annotate(mt_vars)
        assert len(results) >= 1
        assert results[0].variant_name == "m.1555A>G"

    def test_or_tier_classification(self, mock_variants: list[Variant], gene_sets: GeneSetLoader) -> None:
        classifier = ORTierClassifier()
        or_vars = [v for v in mock_variants if v.gene_symbol.startswith("OR")]
        results = classifier.classify(or_vars, "M", gene_sets)
        assert len(results) >= 1
        assert results[0].tier in ("A", "B", "C")

    def test_report_context_building(self, mock_variants: list[Variant], gene_sets: GeneSetLoader) -> None:
        engine = ImpactEngine()
        f = SensoryFilter(gene_sets)
        groups = f.group_by_gene(mock_variants)

        gene_cards = []
        for gene_symbol, variants in groups.items():
            assessments = [engine.assess(v, "M", gene_sets) for v in variants]
            best = max(assessments, key=lambda a: {"无影响": 0, "可能轻微影响": 1, "部分影响": 2, "显著影响": 3, "完全丧失": 4}.get(a.level, 0))
            gene_cards.append(GeneCard(
                gene_symbol=gene_symbol,
                sensory_function_zh=gene_sets.get_gene_function(gene_symbol),
                variants=variants,
                assessment=best,
            ))

        builder = ReportContextBuilder(gene_sets)
        report = builder.build(gene_cards)

        assert isinstance(report, SensoryReport)
        assert report.tas2r38 is not None
        assert report.mitochondrial is not None
        assert report.or_tiers is not None
        assert report.disclaimer_zh != ""
        assert "不构成医学诊断" in report.disclaimer_zh

    def test_markdown_generation(self, mock_variants: list[Variant], gene_sets: GeneSetLoader) -> None:
        engine = ImpactEngine()
        f = SensoryFilter(gene_sets)
        groups = f.group_by_gene(mock_variants)

        gene_cards = []
        for gene_symbol, variants in groups.items():
            assessments = [engine.assess(v, "M", gene_sets) for v in variants]
            best = max(assessments, key=lambda a: {"无影响": 0, "可能轻微影响": 1, "部分影响": 2, "显著影响": 3, "完全丧失": 4}.get(a.level, 0))
            gene_cards.append(GeneCard(
                gene_symbol=gene_symbol,
                sensory_function_zh=gene_sets.get_gene_function(gene_symbol),
                variants=variants,
                assessment=best,
            ))

        builder = ReportContextBuilder(gene_sets)
        report = builder.build(gene_cards)

        md_gen = MarkdownReportGenerator()
        md = md_gen.generate(report)
        assert len(md) > 100
        assert "TEST" not in md  # 未设置 sample_id

    def test_json_generation(self, mock_variants: list[Variant], gene_sets: GeneSetLoader) -> None:
        engine = ImpactEngine()
        f = SensoryFilter(gene_sets)
        groups = f.group_by_gene(mock_variants)

        gene_cards = []
        for gene_symbol, variants in groups.items():
            assessments = [engine.assess(v, "M", gene_sets) for v in variants]
            best = max(assessments, key=lambda a: {"无影响": 0, "可能轻微影响": 1, "部分影响": 2, "显著影响": 3, "完全丧失": 4}.get(a.level, 0))
            gene_cards.append(GeneCard(
                gene_symbol=gene_symbol,
                sensory_function_zh=gene_sets.get_gene_function(gene_symbol),
                variants=variants,
                assessment=best,
            ))

        builder = ReportContextBuilder(gene_sets)
        report = builder.build(gene_cards)

        json_gen = JsonReportGenerator()
        json_str = json_gen.generate(report)
        data = json.loads(json_str)
        assert "gene_cards" in data
        assert "tas2r38" in data
