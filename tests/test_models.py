"""测试 Pydantic 数据模型.

覆盖 Variant, ImpactAssessment, GeneCard, TAS2R38Result, MitochondrialResult,
ORTierResult, ExecutiveSummary, SensoryReport, AnalysisConfig 等核心模型的
构造、序列化、属性计算和边界条件。
"""

import json
from datetime import datetime

import pytest
from pydantic import ValidationError

from src.models import (
    AnalysisConfig,
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


class TestVariant:
    """Variant 模型测试."""

    def test_basic_construction(self) -> None:
        v = Variant(chrom="7", pos=141972755, ref="G", alt="C", gt="0/1")
        assert v.chrom == "7"
        assert v.pos == 141972755
        assert v.ref == "G"
        assert v.alt == "C"
        assert v.gt == "0/1"
        assert v.dp == 0
        assert v.qual == 0.0
        assert v.filter_status == "PASS"

    def test_required_fields(self) -> None:
        with pytest.raises(ValidationError):
            Variant()  # type: ignore[call-arg]

    def test_vcf_id_property(self) -> None:
        v = Variant(chrom="7", pos=100, ref="A", alt="T", gt="0/1")
        assert v.vcf_id == "7:100:A:T"

    def test_is_homozygous(self) -> None:
        assert Variant(chrom="1", pos=1, ref="A", alt="T", gt="1/1").is_homozygous is True
        assert Variant(chrom="1", pos=1, ref="A", alt="T", gt="1|1").is_homozygous is True
        assert Variant(chrom="1", pos=1, ref="A", alt="T", gt="0/1").is_homozygous is False

    def test_is_heterozygous(self) -> None:
        assert Variant(chrom="1", pos=1, ref="A", alt="T", gt="0/1").is_heterozygous is True
        assert Variant(chrom="1", pos=1, ref="A", alt="T", gt="0|1").is_heterozygous is True
        assert Variant(chrom="1", pos=1, ref="A", alt="T", gt="1/0").is_heterozygous is True
        assert Variant(chrom="1", pos=1, ref="A", alt="T", gt="1|0").is_heterozygous is True
        assert Variant(chrom="1", pos=1, ref="A", alt="T", gt="1/1").is_heterozygous is False

    def test_is_hemizygous(self) -> None:
        assert Variant(chrom="X", pos=1, ref="A", alt="T", gt="1").is_hemizygous is True
        assert Variant(chrom="X", pos=1, ref="A", alt="T", gt="0/1").is_hemizygous is False
        assert Variant(chrom="X", pos=1, ref="A", alt="T", gt="1/1").is_hemizygous is False

    def test_hash_and_eq(self) -> None:
        v1 = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1")
        v2 = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1")
        v3 = Variant(chrom="1", pos=100, ref="A", alt="T", gt="1/1")
        assert hash(v1) == hash(v2)
        assert v1 == v2
        assert v1 != v3
        assert v1 != "not a variant"

    def test_optional_vep_fields(self) -> None:
        v = Variant(
            chrom="1",
            pos=100,
            ref="A",
            alt="T",
            gt="0/1",
            gene_symbol="GJB2",
            consequence="missense_variant",
            hgvsc="c.100A>T",
            hgvsp="p.Lys34Met",
            af_gnomad=0.001,
        )
        assert v.gene_symbol == "GJB2"
        assert v.consequence == "missense_variant"
        assert v.hgvsc == "c.100A>T"
        assert v.hgvsp == "p.Lys34Met"
        assert v.af_gnomad == 0.001


class TestImpactAssessment:
    """ImpactAssessment 模型测试."""

    def test_defaults(self) -> None:
        a = ImpactAssessment()
        assert a.level == "无影响"
        assert a.protein_impact == ""
        assert a.gene_certainty == ""
        assert a.zygosity_match is False
        assert a.inheritance_pattern == "未知"

    def test_all_levels(self) -> None:
        for level in ("完全丧失", "显著影响", "部分影响", "可能轻微影响", "无影响"):
            a = ImpactAssessment(level=level)
            assert a.level == level

    def test_invalid_level(self) -> None:
        with pytest.raises(ValidationError):
            ImpactAssessment(level="invalid")  # type: ignore[arg-type]


class TestGeneCard:
    """GeneCard 模型测试."""

    def test_construction(self) -> None:
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1")
        card = GeneCard(
            gene_symbol="GJB2",
            sensory_function_zh="缝隙连接蛋白",
            variants=[v],
            assessment=ImpactAssessment(level="显著影响"),
        )
        assert card.gene_symbol == "GJB2"
        assert card.sensory_function_zh == "缝隙连接蛋白"
        assert len(card.variants) == 1
        assert card.assessment.level == "显著影响"


class TestTAS2R38Result:
    """TAS2R38Result 模型测试."""

    def test_construction(self) -> None:
        r = TAS2R38Result(
            rs713598_gt="C/C",
            rs1726866_gt="T/T",
            rs10246939_gt="G/G",
            diplotype="PAV/PAV",
            phenotype_zh="对苦味化合物高度敏感",
            phenotype_level="苦味高敏感型",
        )
        assert r.diplotype == "PAV/PAV"
        assert r.phenotype_level == "苦味高敏感型"


class TestMitochondrialResult:
    """MitochondrialResult 模型测试."""

    def test_construction(self) -> None:
        r = MitochondrialResult(
            variant_name="m.1555A>G",
            gene="MT-RNR1",
            heteroplasmy=0.85,
            drug_warning_zh="氨基糖苷类抗生素高风险",
            risk_level="高风险",
        )
        assert r.variant_name == "m.1555A>G"
        assert r.heteroplasmy == 0.85

    def test_defaults(self) -> None:
        r = MitochondrialResult(variant_name="m.1555A>G", gene="MT-RNR1")
        assert r.heteroplasmy == 0.0
        assert r.drug_warning_zh == ""
        assert r.risk_level == ""


class TestORTierResult:
    """ORTierResult 模型测试."""

    def test_valid_tiers(self) -> None:
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="1/1")
        for tier in ("A", "B", "C"):
            r = ORTierResult(tier=tier, gene_symbol="OR2T11", variant=v)
            assert r.tier == tier

    def test_invalid_tier(self) -> None:
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="1/1")
        with pytest.raises(ValidationError):
            ORTierResult(tier="D", gene_symbol="OR2T11", variant=v)  # type: ignore[arg-type]


class TestSensoryReport:
    """SensoryReport 模型测试."""

    def test_default_construction(self) -> None:
        report = SensoryReport()
        assert report.sample_id == ""
        assert report.sex == "M"
        assert report.ref_genome == "GRCh38"
        assert isinstance(report.analysis_date, datetime)
        assert report.gene_cards == []

    def test_serialization(self) -> None:
        report = SensoryReport(sample_id="S001", sex="F")
        data = report.model_dump(mode="json")
        assert data["sample_id"] == "S001"
        assert data["sex"] == "F"
        assert "analysis_date" in data

    def test_json_roundtrip(self) -> None:
        report = SensoryReport(sample_id="S001", sex="F")
        json_str = report.model_dump_json()
        parsed = json.loads(json_str)
        assert parsed["sample_id"] == "S001"


class TestAnalysisConfig:
    """AnalysisConfig 模型测试."""

    def test_required_fields(self) -> None:
        with pytest.raises(ValidationError):
            AnalysisConfig()  # type: ignore[call-arg]

    def test_default_subsystems(self) -> None:
        cfg = AnalysisConfig(vcf_path="test.vcf", sex="M")
        assert cfg.subsystems == ["vision", "hearing", "olfaction", "taste", "somatosensation"]

    def test_invalid_sex(self) -> None:
        with pytest.raises(ValidationError):
            AnalysisConfig(vcf_path="test.vcf", sex="X")  # type: ignore[arg-type]

    def test_invalid_subsystem(self) -> None:
        with pytest.raises(ValidationError):
            AnalysisConfig(vcf_path="test.vcf", sex="M", subsystems=["invalid"])  # type: ignore[list-item]
