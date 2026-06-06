"""测试报告生成器与模板渲染.

覆盖 MarkdownReportGenerator、JsonReportGenerator 及 Jinja2 模板渲染。
"""

import json
from datetime import datetime

import pytest

from src.models import (
    DataAvailability,
    ExecutiveSummary,
    GeneCard,
    ImpactAssessment,
    MitochondrialResult,
    ORTierResult,
    SensoryReport,
    TAS2R38Result,
    Variant,
)
from src.report.json_generator import JsonReportGenerator
from src.report.markdown_generator import MarkdownReportGenerator


@pytest.fixture
def mock_report() -> SensoryReport:
    """构建一个包含完整数据的 mock 报告."""
    v1 = Variant(
        chrom="7", pos=141972755, ref="G", alt="C", gt="C/C",
        gene_symbol="TAS2R38", consequence="missense_variant",
    )
    card1 = GeneCard(
        gene_symbol="TAS2R38",
        sensory_function_zh="苦味受体",
        variants=[v1],
        assessment=ImpactAssessment(
            level="无影响",
            protein_impact="蛋白结构改变（错义/非框内缺失）",
            gene_certainty="高",
            zygosity_match=False,
            inheritance_pattern="隐性纯合",
            rationale_zh="蛋白影响: ...",
        ),
    )
    return SensoryReport(
        sample_id="TEST001",
        sex="M",
        ref_genome="GRCh38",
        analysis_date=datetime(2024, 6, 1, 12, 0, 0),
        subsystems=["vision", "hearing", "taste"],
        gene_cards=[card1],
        tas2r38=TAS2R38Result(
            rs713598_gt="C/C",
            rs1726866_gt="T/T",
            rs10246939_gt="G/G",
            diplotype="PAV/PAV",
            phenotype_zh="对苦味化合物高度敏感",
            phenotype_level="苦味高敏感型",
        ),
        mitochondrial=[
            MitochondrialResult(
                variant_name="m.1555A>G",
                gene="MT-RNR1",
                heteroplasmy=1.0,
                drug_warning_zh="氨基糖苷类抗生素高风险",
                risk_level="高风险",
            )
        ],
        or_tiers=[
            ORTierResult(
                tier="A",
                gene_symbol="OR2T11",
                known_ligand_zh="甲硫醇",
                variant=Variant(chrom="1", pos=100, ref="A", alt="T", gt="1/1", gene_symbol="OR2T11", consequence="frameshift_variant"),
                assessment=ImpactAssessment(level="完全丧失"),
            )
        ],
        executive_summary=ExecutiveSummary(
            subsystem_counts={},
            key_findings=["OR2T11: 完全丧失"],
        ),
        data_availability={
            "TAS2R38": DataAvailability(gnomad_af="可用", clinvar="N/A", topology="可用"),
        },
        disclaimer_zh="【重要声明】本报告不构成医学诊断...",
    )


class TestMarkdownReportGenerator:
    """Markdown 报告生成器测试."""

    def test_generate_basic(self, mock_report: SensoryReport) -> None:
        gen = MarkdownReportGenerator()
        md = gen.generate(mock_report)
        assert "# 个体感官能力基因分析报告" in md or "个体感官" in md or "感官能力" in md
        assert "TEST001" in md
        assert "TAS2R38" in md

    def test_level_badge_filter(self) -> None:
        gen = MarkdownReportGenerator()
        assert "完全丧失" in gen._level_badge_filter("完全丧失")
        assert "显著影响" in gen._level_badge_filter("显著影响")
        assert gen._level_badge_filter("未知") == "未知"

    def test_risk_badge_filter(self) -> None:
        gen = MarkdownReportGenerator()
        assert "高风险" in gen._risk_badge_filter("高风险")
        assert "⚪ 未知" == gen._risk_badge_filter("未知")

    def test_template_rendering_with_tas2r38(self, mock_report: SensoryReport) -> None:
        gen = MarkdownReportGenerator()
        md = gen.generate(mock_report)
        assert "PAV/PAV" in md or "苦味" in md

    def test_template_rendering_with_mitochondrial(self, mock_report: SensoryReport) -> None:
        gen = MarkdownReportGenerator()
        md = gen.generate(mock_report)
        assert "m.1555A>G" in md or "MT-RNR1" in md or "氨基糖苷" in md

    def test_disclaimer_in_report(self, mock_report: SensoryReport) -> None:
        gen = MarkdownReportGenerator()
        md = gen.generate(mock_report)
        assert "不构成医学诊断" in md

    def test_empty_report(self) -> None:
        gen = MarkdownReportGenerator()
        report = SensoryReport(sample_id="EMPTY", sex="F")
        md = gen.generate(report)
        assert "EMPTY" in md


class TestJsonReportGenerator:
    """JSON 报告生成器测试."""

    def test_generate_json(self, mock_report: SensoryReport) -> None:
        gen = JsonReportGenerator()
        json_str = gen.generate(mock_report)
        data = json.loads(json_str)
        assert data["sample_id"] == "TEST001"
        assert data["sex"] == "M"
        assert data["tas2r38"]["diplotype"] == "PAV/PAV"

    def test_generate_dict(self, mock_report: SensoryReport) -> None:
        gen = JsonReportGenerator()
        data = gen.generate_dict(mock_report)
        assert data["sample_id"] == "TEST001"
        assert isinstance(data["analysis_date"], str)

    def test_json_contains_all_fields(self, mock_report: SensoryReport) -> None:
        gen = JsonReportGenerator()
        data = gen.generate_dict(mock_report)
        assert "gene_cards" in data
        assert "tas2r38" in data
        assert "mitochondrial" in data
        assert "or_tiers" in data
        assert "executive_summary" in data
        assert "data_availability" in data
        assert "disclaimer_zh" in data

    def test_json_unicode(self, mock_report: SensoryReport) -> None:
        gen = JsonReportGenerator()
        json_str = gen.generate(mock_report)
        # 确保中文未被转义
        assert "苦味" in json_str
        assert "\\u" not in json_str  # 不应出现 unicode escape
