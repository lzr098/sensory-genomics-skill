"""测试功能影响评估体系.

覆盖 ImpactRule（ProteinImpactRule, GeneCertaintyRule）、
InheritanceMatcher（遗传模式匹配）、ImpactEngine（五级评估引擎）。
"""

import pytest

from src.assessment.engine import ImpactEngine
from src.assessment.inheritance import InheritanceMatcher
from src.assessment.rules import GeneCertaintyRule, ProteinImpactRule
from src.gene_sets.loader import GeneSetLoader
from src.models import Variant


@pytest.fixture
def gene_sets() -> GeneSetLoader:
    return GeneSetLoader()


@pytest.fixture
def matcher() -> InheritanceMatcher:
    return InheritanceMatcher()


@pytest.fixture
def engine() -> ImpactEngine:
    return ImpactEngine()


class TestProteinImpactRule:
    """蛋白影响规则测试."""

    def test_high_impact(self, gene_sets: GeneSetLoader) -> None:
        rule = ProteinImpactRule()
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1", consequence="stop_gained")
        result = rule.evaluate(v, "M", gene_sets)
        assert "严重受损" in result

    def test_moderate_impact(self, gene_sets: GeneSetLoader) -> None:
        rule = ProteinImpactRule()
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1", consequence="missense_variant")
        result = rule.evaluate(v, "M", gene_sets)
        assert "结构改变" in result

    def test_low_impact(self, gene_sets: GeneSetLoader) -> None:
        rule = ProteinImpactRule()
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1", consequence="synonymous_variant")
        result = rule.evaluate(v, "M", gene_sets)
        assert "保守" in result

    def test_modifier_impact(self, gene_sets: GeneSetLoader) -> None:
        rule = ProteinImpactRule()
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1", consequence="intron_variant")
        result = rule.evaluate(v, "M", gene_sets)
        assert "非编码" in result or "调控区域" in result

    def test_multiple_consequences(self, gene_sets: GeneSetLoader) -> None:
        """多个 consequence 取最高影响等级."""
        rule = ProteinImpactRule()
        v = Variant(
            chrom="1", pos=100, ref="A", alt="T", gt="0/1",
            consequence="missense_variant,synonymous_variant"
        )
        result = rule.evaluate(v, "M", gene_sets)
        assert "结构改变" in result

    def test_unknown_consequence(self, gene_sets: GeneSetLoader) -> None:
        rule = ProteinImpactRule()
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1", consequence="totally_unknown")
        result = rule.evaluate(v, "M", gene_sets)
        assert result == "非编码或调控区域变异"


class TestGeneCertaintyRule:
    """基因确定性规则测试."""

    def test_high_certainty(self, gene_sets: GeneSetLoader) -> None:
        rule = GeneCertaintyRule()
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1", gene_symbol="GJB2")
        result = rule.evaluate(v, "M", gene_sets)
        assert result == "高"

    def test_or_gene_medium(self, gene_sets: GeneSetLoader) -> None:
        rule = GeneCertaintyRule()
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1", gene_symbol="OR2T11")
        result = rule.evaluate(v, "M", gene_sets)
        assert result == "中"

    def test_unknown_gene_medium(self, gene_sets: GeneSetLoader) -> None:
        rule = GeneCertaintyRule()
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1", gene_symbol="UNKNOWN")
        result = rule.evaluate(v, "M", gene_sets)
        assert result == "中"


class TestInheritanceMatcher:
    """遗传模式匹配器测试."""

    def test_detect_mitochondrial(self, matcher: InheritanceMatcher) -> None:
        assert matcher.detect_pattern("MT-RNR1", "M") == "线粒体"
        assert matcher.detect_pattern("MT-TS1", "F") == "线粒体"

    def test_detect_x_linked(self, matcher: InheritanceMatcher) -> None:
        assert matcher.detect_pattern("OPN1LW", "M") == "X-连锁"
        assert matcher.detect_pattern("OPN1MW", "F") == "X-连锁"

    def test_detect_dominant(self, matcher: InheritanceMatcher) -> None:
        assert matcher.detect_pattern("SCN9A", "M") == "显性"
        assert matcher.detect_pattern("SCN10A", "F") == "显性"

    def test_detect_recessive(self, matcher: InheritanceMatcher) -> None:
        assert matcher.detect_pattern("GJB2", "M") == "隐性纯合"
        assert matcher.detect_pattern("MYO7A", "F") == "隐性纯合"

    def test_detect_or_recessive(self, matcher: InheritanceMatcher) -> None:
        assert matcher.detect_pattern("OR2T11", "M") == "隐性纯合"
        assert matcher.detect_pattern("OR7D4", "F") == "隐性纯合"

    def test_detect_unknown(self, matcher: InheritanceMatcher) -> None:
        assert matcher.detect_pattern("UNKNOWN", "M") == "未知"

    def test_zygosity_match_mitochondrial(self, matcher: InheritanceMatcher) -> None:
        v = Variant(chrom="MT", pos=1555, ref="A", alt="G", gt="1")
        assert matcher.zygosity_match(v, "M", "线粒体") is True
        v2 = Variant(chrom="7", pos=100, ref="A", alt="T", gt="0/1")
        assert matcher.zygosity_match(v2, "M", "线粒体") is False

    def test_zygosity_match_x_linked_male(self, matcher: InheritanceMatcher) -> None:
        # 男性半合子
        v = Variant(chrom="X", pos=100, ref="A", alt="T", gt="1")
        assert matcher.zygosity_match(v, "M", "X-连锁") is True
        # 男性纯合
        v2 = Variant(chrom="X", pos=100, ref="A", alt="T", gt="1/1")
        assert matcher.zygosity_match(v2, "M", "X-连锁") is True
        # 男性杂合
        v3 = Variant(chrom="X", pos=100, ref="A", alt="T", gt="0/1")
        assert matcher.zygosity_match(v3, "M", "X-连锁") is False

    def test_zygosity_match_x_linked_female(self, matcher: InheritanceMatcher) -> None:
        # 女性纯合
        v = Variant(chrom="X", pos=100, ref="A", alt="T", gt="1/1")
        assert matcher.zygosity_match(v, "F", "X-连锁") is True
        # 女性杂合
        v2 = Variant(chrom="X", pos=100, ref="A", alt="T", gt="0/1")
        assert matcher.zygosity_match(v2, "F", "X-连锁") is True

    def test_zygosity_match_dominant(self, matcher: InheritanceMatcher) -> None:
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1")
        assert matcher.zygosity_match(v, "M", "显性") is True
        v2 = Variant(chrom="1", pos=100, ref="A", alt="T", gt="1/1")
        assert matcher.zygosity_match(v2, "M", "显性") is True

    def test_zygosity_match_recessive(self, matcher: InheritanceMatcher) -> None:
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="1/1")
        assert matcher.zygosity_match(v, "M", "隐性纯合") is True
        v2 = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1")
        assert matcher.zygosity_match(v2, "M", "隐性纯合") is False

    def test_zygosity_match_unknown(self, matcher: InheritanceMatcher) -> None:
        v = Variant(chrom="1", pos=100, ref="A", alt="T", gt="0/1")
        assert matcher.zygosity_match(v, "M", "未知") is True

    def test_explain_pattern(self, matcher: InheritanceMatcher) -> None:
        assert "母系遗传" in matcher.explain_pattern("MT-RNR1", "M")
        assert "男性半合子" in matcher.explain_pattern("OPN1LW", "M")
        assert "女性纯合或杂合均可表型" in matcher.explain_pattern("OPN1LW", "F")
        assert "常染色体显性" in matcher.explain_pattern("SCN9A", "M")
        assert "常染色体隐性" in matcher.explain_pattern("GJB2", "M")


class TestImpactEngine:
    """五级评估引擎测试."""

    def test_or_homozygous_lof_complete_loss(self, engine: ImpactEngine, gene_sets: GeneSetLoader) -> None:
        """OR 基因纯合 LoF → 完全丧失."""
        v = Variant(
            chrom="1", pos=100, ref="A", alt="T", gt="1/1",
            gene_symbol="OR2T11", consequence="frameshift_variant"
        )
        result = engine.assess(v, "M", gene_sets)
        assert result.level == "完全丧失"
        assert result.zygosity_match is True
        assert result.inheritance_pattern == "隐性纯合"

    def test_or_homozygous_missense_complete_loss(self, engine: ImpactEngine, gene_sets: GeneSetLoader) -> None:
        """OR 基因纯合错义 → 完全丧失."""
        v = Variant(
            chrom="1", pos=100, ref="A", alt="T", gt="1/1",
            gene_symbol="OR2T11", consequence="missense_variant"
        )
        result = engine.assess(v, "M", gene_sets)
        assert result.level == "完全丧失"

    def test_or_heterozygous_no_effect(self, engine: ImpactEngine, gene_sets: GeneSetLoader) -> None:
        """OR 基因杂合错义 → 无影响."""
        v = Variant(
            chrom="1", pos=100, ref="A", alt="T", gt="0/1",
            gene_symbol="OR2T11", consequence="missense_variant"
        )
        result = engine.assess(v, "M", gene_sets)
        assert result.level == "无影响"

    def test_gjb2_homozygous_frameshift_significant(self, engine: ImpactEngine, gene_sets: GeneSetLoader) -> None:
        """GJB2 纯合 frameshift → 显著影响."""
        v = Variant(
            chrom="13", pos=100, ref="A", alt="T", gt="1/1",
            gene_symbol="GJB2", consequence="frameshift_variant"
        )
        result = engine.assess(v, "M", gene_sets)
        assert result.level == "显著影响"

    def test_gjb2_homozygous_missense_partial(self, engine: ImpactEngine, gene_sets: GeneSetLoader) -> None:
        """GJB2 纯合错义 → 部分影响."""
        v = Variant(
            chrom="13", pos=100, ref="A", alt="T", gt="1/1",
            gene_symbol="GJB2", consequence="missense_variant"
        )
        result = engine.assess(v, "M", gene_sets)
        assert result.level == "部分影响"

    def test_scn9a_homozygous_frameshift_complete_loss(self, engine: ImpactEngine, gene_sets: GeneSetLoader) -> None:
        """SCN9A 纯合 LoF → 完全丧失（痛觉）."""
        v = Variant(
            chrom="2", pos=100, ref="A", alt="T", gt="1/1",
            gene_symbol="SCN9A", consequence="stop_gained"
        )
        result = engine.assess(v, "M", gene_sets)
        assert result.level == "完全丧失"

    def test_dominant_heterozygous_lof(self, engine: ImpactEngine, gene_sets: GeneSetLoader) -> None:
        """显性基因杂合 LoF → 显著影响."""
        v = Variant(
            chrom="1", pos=100, ref="A", alt="T", gt="0/1",
            gene_symbol="SCN9A", consequence="frameshift_variant"
        )
        result = engine.assess(v, "M", gene_sets)
        assert result.level == "显著影响"

    def test_opn1lw_heterozygous_missense(self, engine: ImpactEngine, gene_sets: GeneSetLoader) -> None:
        """OPN1LW 杂合错义 → 可能轻微影响."""
        v = Variant(
            chrom="X", pos=100, ref="A", alt="T", gt="0/1",
            gene_symbol="OPN1LW", consequence="missense_variant"
        )
        result = engine.assess(v, "F", gene_sets)
        assert result.level == "可能轻微影响"
        assert result.limitation_note is not None
        assert "CNV" in result.limitation_note

    def test_synonymous_no_effect(self, engine: ImpactEngine, gene_sets: GeneSetLoader) -> None:
        """同义变异 → 无影响."""
        v = Variant(
            chrom="1", pos=100, ref="A", alt="T", gt="0/1",
            gene_symbol="GJB2", consequence="synonymous_variant"
        )
        result = engine.assess(v, "M", gene_sets)
        assert result.level == "无影响"

    def test_utr_no_effect(self, engine: ImpactEngine, gene_sets: GeneSetLoader) -> None:
        """UTR 变异 → 无影响."""
        v = Variant(
            chrom="1", pos=100, ref="A", alt="T", gt="0/1",
            gene_symbol="GJB2", consequence="3_prime_UTR_variant"
        )
        result = engine.assess(v, "M", gene_sets)
        assert result.level == "无影响"

    def test_rationale_building(self, engine: ImpactEngine, gene_sets: GeneSetLoader) -> None:
        v = Variant(
            chrom="1", pos=100, ref="A", alt="T", gt="1/1",
            gene_symbol="GJB2", consequence="missense_variant"
        )
        result = engine.assess(v, "M", gene_sets)
        assert "蛋白影响" in result.rationale_zh
        assert "基因功能确定性" in result.rationale_zh
        assert "遗传模式" in result.rationale_zh

    def test_trpv1_heterozygous_missense(self, engine: ImpactEngine, gene_sets: GeneSetLoader) -> None:
        """TRPV1 杂合错义 → 可能轻微影响."""
        v = Variant(
            chrom="17", pos=100, ref="A", alt="T", gt="0/1",
            gene_symbol="TRPV1", consequence="missense_variant"
        )
        result = engine.assess(v, "M", gene_sets)
        assert result.level == "可能轻微影响"

    def test_trpv1_homozygous_missense(self, engine: ImpactEngine, gene_sets: GeneSetLoader) -> None:
        """TRPV1 纯合错义 → 可能轻微影响."""
        v = Variant(
            chrom="17", pos=100, ref="A", alt="T", gt="1/1",
            gene_symbol="TRPV1", consequence="missense_variant"
        )
        result = engine.assess(v, "M", gene_sets)
        assert result.level == "可能轻微影响"
