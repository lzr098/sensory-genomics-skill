"""ImpactRule 接口及具体规则实现.

- ProteinImpactRule: 基于 VEP consequence 评估蛋白影响
- GeneCertaintyRule: 基于基因-表型关联确定性评估
"""

from abc import ABC, abstractmethod

from src.gene_sets.loader import GeneSetLoader
from src.models import Sex, Variant


class ImpactRule(ABC):
    """功能影响规则接口."""

    name: str = ""

    @abstractmethod
    def evaluate(self, variant: Variant, sex: Sex, gene_sets: GeneSetLoader) -> str:
        """评估变异并返回结果文本.

        Args:
            variant: 已注释的变异。
            sex: 样本性别。
            gene_sets: 感官基因集加载器。

        Returns:
            评估结果描述文本。
        """
        raise NotImplementedError


class ProteinImpactRule(ImpactRule):
    """蛋白影响规则.

    基于 VEP consequence 术语评估变异对蛋白功能的影响程度。
    """

    name = "protein_impact"

    # VEP consequence 影响等级映射
    IMPACT_LEVELS = {
        "high": [
            "transcript_ablation",
            "splice_acceptor_variant",
            "splice_donor_variant",
            "stop_gained",
            "frameshift_variant",
            "stop_lost",
            "start_lost",
            "transcript_amplification",
        ],
        "moderate": [
            "inframe_insertion",
            "inframe_deletion",
            "missense_variant",
            "protein_altering_variant",
        ],
        "low": [
            "splice_region_variant",
            "splice_donor_5th_base_variant",
            "splice_donor_region_variant",
            "splice_polypyrimidine_tract_variant",
            "incomplete_terminal_codon_variant",
            "start_retained_variant",
            "stop_retained_variant",
            "synonymous_variant",
        ],
        "modifier": [
            "coding_sequence_variant",
            "mature_miRNA_variant",
            "5_prime_UTR_variant",
            "3_prime_UTR_variant",
            "non_coding_transcript_exon_variant",
            "intron_variant",
            "NMD_transcript_variant",
            "non_coding_transcript_variant",
            "upstream_gene_variant",
            "downstream_gene_variant",
            "TFBS_ablation",
            "TFBS_amplification",
            "TF_binding_site_variant",
            "regulatory_region_ablation",
            "regulatory_region_amplification",
            "feature_elongation",
            "regulatory_region_variant",
            "feature_truncation",
            "intergenic_variant",
        ],
    }

    def evaluate(self, variant: Variant, sex: Sex, gene_sets: GeneSetLoader) -> str:
        """评估蛋白影响.

        Args:
            variant: 已注释的变异。
            sex: 样本性别。
            gene_sets: 感官基因集加载器。

        Returns:
            蛋白影响评估文本。
        """
        consequence = variant.consequence or ""
        terms = [t.strip() for t in consequence.split(",")]

        highest = "modifier"
        for term in terms:
            for level, term_list in self.IMPACT_LEVELS.items():
                if term in term_list:
                    if self._level_priority(level) < self._level_priority(highest):
                        highest = level
                    break

        mapping = {
            "high": "蛋白功能严重受损（移码/无义/剪接位点）",
            "moderate": "蛋白结构改变（错义/非框内缺失）",
            "low": "蛋白序列保守（同义/UTR）",
            "modifier": "非编码或调控区域变异",
        }
        return mapping.get(highest, "未知")

    @staticmethod
    def _level_priority(level: str) -> int:
        """返回影响等级优先级数值（越小越高）."""
        order = {"high": 0, "moderate": 1, "low": 2, "modifier": 3}
        return order.get(level, 99)


class GeneCertaintyRule(ImpactRule):
    """基因功能确定性规则.

    基于已知的基因-表型关联强度评估该基因在感官系统中的功能确定性。
    """

    name = "gene_certainty"

    # 高确定性基因：OMIM/ClinVar 中已明确与感官表型关联
    HIGH_CERTAINTY_GENES = {
        "GJB2", "GJB6", "MYO7A", "CDH23", "OTOF", "CABP2", "SLC26A4", "KCNQ4",
        "OPN1LW", "OPN1MW", "OPN1SW", "RHO", "ABCA4", "EYS", "GJA8", "CRYAA",
        "OCA2", "TYRP1",
        "TAS2R38", "TAS1R2", "TAS1R3", "TAS1R1", "SCNN1A", "SCNN1B", "SCNN1G",
        "OTOP1", "OTOP2", "OTOP3",
        "SCN9A", "SCN10A", "SCN11A", "TRPV1", "TRPM8", "PIEZO1", "PIEZO2",
        "NGF", "NTRK1", "GRPR", "IL31RA",
        "CNGA2", "ADCY3",
    }

    def evaluate(self, variant: Variant, sex: Sex, gene_sets: GeneSetLoader) -> str:
        """评估基因功能确定性.

        Args:
            variant: 已注释的变异。
            sex: 样本性别。
            gene_sets: 感官基因集加载器。

        Returns:
            基因确定性评估文本（高/中/低）。
        """
        gene = variant.gene_symbol
        if gene in self.HIGH_CERTAINTY_GENES:
            return "高"
        # OR 基因大多数有功能研究但表型关联证据较少
        if gene.startswith("OR"):
            return "中"
        return "中"
