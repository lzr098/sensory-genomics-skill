"""功能影响评估引擎.

编排多条 ImpactRule，综合蛋白影响 × 基因确定性 × 遗传模式，输出五级评估。
"""

from typing import Any, List

from src.assessment.inheritance import InheritanceMatcher
from src.assessment.rules import GeneCertaintyRule, ImpactRule, ProteinImpactRule
from src.gene_sets.loader import GeneSetLoader
from src.logger import get_logger
from src.models import ImpactAssessment, ImpactLevel, Sex, Variant

logger = get_logger(__name__)

# 后果类型到蛋白影响的映射（简化版）
_CONSEQUENCE_IMPACT = {
    "transcript_ablation": "high",
    "splice_acceptor_variant": "high",
    "splice_donor_variant": "high",
    "stop_gained": "high",
    "frameshift_variant": "high",
    "stop_lost": "high",
    "start_lost": "high",
    "transcript_amplification": "high",
    "inframe_insertion": "moderate",
    "inframe_deletion": "moderate",
    "missense_variant": "moderate",
    "protein_altering_variant": "moderate",
    "splice_region_variant": "low",
    "splice_donor_5th_base_variant": "low",
    "splice_donor_region_variant": "low",
    "splice_polypyrimidine_tract_variant": "low",
    "incomplete_terminal_codon_variant": "low",
    "start_retained_variant": "low",
    "stop_retained_variant": "low",
    "synonymous_variant": "low",
    "coding_sequence_variant": "modifier",
    "mature_miRNA_variant": "modifier",
    "5_prime_UTR_variant": "modifier",
    "3_prime_UTR_variant": "modifier",
    "non_coding_transcript_exon_variant": "modifier",
    "intron_variant": "modifier",
    "NMD_transcript_variant": "modifier",
    "non_coding_transcript_variant": "modifier",
    "upstream_gene_variant": "modifier",
    "downstream_gene_variant": "modifier",
    "TFBS_ablation": "modifier",
    "TFBS_amplification": "modifier",
    "TF_binding_site_variant": "modifier",
    "regulatory_region_ablation": "modifier",
    "regulatory_region_amplification": "modifier",
    "feature_elongation": "modifier",
    "regulatory_region_variant": "modifier",
    "feature_truncation": "modifier",
    "intergenic_variant": "modifier",
}

# 基因确定性等级
_GENE_CERTAINTY = {
    "GJB2": "高", "GJB6": "高", "MYO7A": "高", "CDH23": "高", "OTOF": "高",
    "SLC26A4": "高", "KCNQ4": "高",
    "OPN1LW": "高", "OPN1MW": "高", "OPN1SW": "高", "RHO": "高",
    "ABCA4": "高", "EYS": "高", "GJA8": "高", "CRYAA": "高",
    "OCA2": "高", "TYRP1": "高",
    "TAS2R38": "高", "TAS2R16": "高", "TAS2R19": "高",
    "TAS1R2": "高", "TAS1R3": "高", "TAS1R1": "高",
    "SCNN1A": "高", "SCNN1B": "高", "SCNN1G": "高",
    "OTOP1": "高", "OTOP2": "高", "OTOP3": "高",
    "SCN9A": "高", "SCN10A": "高", "SCN11A": "高",
    "TRPV1": "高", "TRPM8": "高",
    "PIEZO1": "高", "PIEZO2": "高",
    "NGF": "高", "NTRK1": "高", "GRPR": "高", "IL31RA": "高",
    "CNGA2": "高", "ADCY3": "高", "CABP2": "高",
    "HERC2": "高", "SLC45A2": "高", "SLC24A4": "高", "SLC24A5": "高", "IRF4": "高",
    "CYP1A2": "高", "LCT": "高", "ALDH2": "高", "ADH1B": "高",
    "ACTN3": "高", "EDAR": "高", "OPRM1": "高", "COMT": "高",
}

# 特殊基因处理
_SPECIAL_GENE_RULES = {
    "TAS2R38": "taste_special",
    "TAS2R16": "taste_special",
    "TAS2R19": "taste_special",
    "MT-RNR1": "mitochondrial",
    "MT-TS1": "mitochondrial",
    "HERC2": "pigmentation_special",
    "LCT": "metabolism_special",
    "CYP1A2": "metabolism_special",
    "ALDH2": "metabolism_special",
    "ADH1B": "metabolism_special",
    "ACTN3": "muscle_special",
    "EDAR": "hair_special",
    "COMT": "pain_special",
    "OPRM1": "pain_special",
}


class ImpactEngine:
    """功能影响评估引擎.

    综合蛋白影响、基因确定性、遗传模式进行五级评估。
    对继承模式和基因确定性结果按基因缓存，避免同基因多变异重复计算。
    """

    def __init__(self) -> None:
        """初始化评估引擎，注册默认规则."""
        self.rules: List[ImpactRule] = [
            ProteinImpactRule(),
            GeneCertaintyRule(),
        ]
        self.inheritance_matcher = InheritanceMatcher()
        # 按 (gene, sex) 缓存继承模式结果
        self._inheritance_cache: Dict[Any, Any] = {}
        # 按基因缓存确定性等级
        self._certainty_cache: Dict[str, str] = {}

    def _get_inheritance(self, gene: str, sex: Sex) -> Any:
        """获取遗传模式（带缓存）."""
        key = (gene, sex)
        if key not in self._inheritance_cache:
            self._inheritance_cache[key] = self.inheritance_matcher.detect_pattern(gene, sex)
        return self._inheritance_cache[key]

    def _get_gene_certainty(self, gene: str) -> str:
        """获取基因确定性等级（带缓存）."""
        if gene not in self._certainty_cache:
            self._certainty_cache[gene] = _GENE_CERTAINTY.get(gene, "中")
        return self._certainty_cache[gene]

    def assess(
        self, variant: Variant, sex: Sex, gene_sets: GeneSetLoader
    ) -> ImpactAssessment:
        """对单个变异进行功能影响评估.

        Args:
            variant: 已注释的变异。
            sex: 样本性别。
            gene_sets: 感官基因集加载器。

        Returns:
            ImpactAssessment 评估结果。
        """
        gene = variant.gene_symbol

        # 1. 蛋白影响评估
        protein_impact = self._assess_protein_impact(variant)

        # 2. 基因确定性（缓存）
        gene_certainty = self._get_gene_certainty(gene)

        # 3. 遗传模式匹配（缓存）
        inheritance = self._get_inheritance(gene, sex)
        zygosity_match = self.inheritance_matcher.zygosity_match(variant, sex, inheritance)

        # 4. 综合五级评估
        level = self._determine_level(
            variant, protein_impact, gene_certainty, zygosity_match, inheritance, gene_sets
        )

        # 5. 生成依据文本
        rationale = self._build_rationale(
            variant, protein_impact, gene_certainty, zygosity_match, inheritance
        )

        limitation = None
        if gene in ("OPN1LW", "OPN1MW"):
            limitation = "本分析仅覆盖 SNV/indel，未检测拷贝数变异（CNV）和基因重组，色觉分析不完整。"

        return ImpactAssessment(
            level=level,
            protein_impact=protein_impact,
            gene_certainty=gene_certainty,
            zygosity_match=zygosity_match,
            inheritance_pattern=inheritance,
            rationale_zh=rationale,
            limitation_note=limitation,
        )

    @staticmethod
    def _assess_protein_impact(variant: Variant) -> str:
        """评估蛋白影响."""
        consequence = variant.consequence or ""
        # VEP consequence 可能包含多个，取最高影响
        impacts = []
        for term in consequence.split(","):
            term = term.strip()
            impact = _CONSEQUENCE_IMPACT.get(term, "modifier")
            impacts.append(impact)

        if not impacts:
            return "未知"

        priority = ["high", "moderate", "low", "modifier"]
        highest = "modifier"
        for p in priority:
            if p in impacts:
                highest = p
                break

        mapping = {
            "high": "蛋白功能严重受损（移码/无义/剪接位点）",
            "moderate": "蛋白结构改变（错义/非框内缺失）",
            "low": "蛋白序列保守（同义/UTR）",
            "modifier": "非编码或调控区域变异",
        }
        return mapping.get(highest, "未知")

    @staticmethod
    def _determine_level(
        variant: Variant,
        protein_impact: str,
        gene_certainty: str,
        zygosity_match: bool,
        inheritance: str,
        gene_sets: GeneSetLoader,
    ) -> ImpactLevel:
        """确定五级影响程度.

        评分原则（保守优先）：
        - 完全丧失：仅用于明确的纯合功能丧失（frameshift/stop_gained/splice_donor/splice_acceptor）
          且该基因功能确定性高、无冗余替代通路。
        - 显著影响：纯合 LoF 但基因冗余度较高，或强致病基因杂合 LoF（显性）。
        - 部分影响：纯合错义/剪接区变异，或重要基因的功能改变。
        - 可能轻微影响：调控区/UTR纯合变异，或冗余基因家族的错义变异。
        - 无影响：同义、杂合错义（隐性遗传）、非编码背景变异。
        """
        gene = variant.gene_symbol
        consequence = variant.consequence or ""

        # ── 1. 明确的功能丧失（LoF）─────────────────────────────────────
        is_lof = any(t in consequence for t in [
            "frameshift_variant", "stop_gained", "stop_lost", "start_lost",
            "splice_acceptor_variant", "splice_donor_variant",
            "transcript_ablation",
        ])

        if is_lof:
            if variant.is_homozygous:
                # 纯合 LoF：只有高确定性关键基因才判"完全丧失"
                if gene in ("SCN9A", "NTRK1"):
                    return "完全丧失"
                # 重要感官基因 → "显著影响"
                if gene in ("GJB2", "GJB6", "OTOF", "SLC26A4", "MYO7A", "CDH23",
                            "OPN1LW", "OPN1MW", "OPN1SW", "RHO", "ABCA4", "EYS"):
                    return "显著影响"
                # OR 基因虽有纯合 LoF，但嗅觉受体高度冗余（~400个），不应判"完全丧失"
                if gene.startswith("OR"):
                    return "可能轻微影响"
                # 其他基因保守处理
                return "部分影响"

            if variant.is_heterozygous:
                # 杂合 LoF + 显性遗传 → 显著影响
                if inheritance == "显性":
                    return "显著影响"
                # 杂合 LoF + 隐性遗传 → 通常携带者状态，无影响
                return "无影响"

        # ── 2. 错义变异（missense）───────────────────────────────────────
        is_missense = "missense_variant" in consequence or "protein_altering_variant" in consequence

        if is_missense:
            if variant.is_homozygous:
                # 纯合错义：不等于功能丧失，最多"部分影响"
                if gene in ("GJB2", "GJB6", "OTOF", "SLC26A4", "MYO7A", "CDH23"):
                    return "部分影响"
                if gene in ("SCN9A", "NTRK1"):
                    return "部分影响"
                # OR 基因纯合错义非常常见，大量是多态性 → "可能轻微影响"
                if gene.startswith("OR"):
                    return "可能轻微影响"
                # 其他基因保守处理
                return "可能轻微影响"

            if variant.is_heterozygous:
                # 杂合错义：隐性遗传下是携带者，不表现表型
                if gene in ("OPN1LW", "OPN1MW", "OPN1SW"):
                    return "可能轻微影响"
                if gene in ("TRPV1", "TRPM8"):
                    return "可能轻微影响"
                # 默认：杂合错义无表型影响
                return "无影响"

        # ── 3. 同义变异 ────────────────────────────────────────────────
        if "synonymous_variant" in consequence:
            return "无影响"

        # ── 4. 剪接区域变异（非经典剪接位点）─────────────────────────────
        is_splice_region = any(t in consequence for t in [
            "splice_region_variant", "splice_donor_5th_base_variant",
            "splice_donor_region_variant", "splice_polypyrimidine_tract_variant",
        ])
        if is_splice_region:
            if variant.is_homozygous:
                return "可能轻微影响"
            return "无影响"

        # ── 5. UTR 变异 ────────────────────────────────────────────────
        if "UTR" in consequence:
            if variant.is_homozygous:
                return "可能轻微影响"
            return "无影响"

        # ── 6. 调控区变异 ──────────────────────────────────────────────
        if any(t in consequence for t in ["upstream_gene", "downstream_gene", "regulatory_region", "TFBS", "TF_binding_site"]):
            if variant.is_homozygous:
                return "可能轻微影响"
            return "无影响"

        # ── 7. 内含子 / 基因间区 ───────────────────────────────────────
        if "intron" in consequence or "intergenic" in consequence:
            return "无影响"

        # 默认保守处理
        return "无影响"

    @staticmethod
    def _build_rationale(
        variant: Variant,
        protein_impact: str,
        gene_certainty: str,
        zygosity_match: bool,
        inheritance: str,
    ) -> str:
        """生成影响依据文本."""
        parts = [
            f"蛋白影响: {protein_impact}",
            f"基因功能确定性: {gene_certainty}",
            f"遗传模式: {inheritance}",
        ]
        if zygosity_match:
            parts.append("基因型与遗传模式匹配")
        else:
            parts.append("基因型与遗传模式不匹配，但变异本身可能影响蛋白功能")
        return "；".join(parts)
