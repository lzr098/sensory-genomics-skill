"""遗传模式匹配器.

支持隐性（纯合/复合杂合）、X-连锁（男性半合子/女性纯合/杂合）、显性、线粒体四种模式。
"""

from typing import Optional

from src.logger import get_logger
from src.models import InheritancePattern, Sex, Variant

logger = get_logger(__name__)

# X-连锁基因（OPN1SW 位于 7 号染色体，非 X-连锁）
_X_LINKED_GENES = {"OPN1LW", "OPN1MW"}

# 线粒体基因
_MITOCHONDRIAL_GENES = {"MT-RNR1", "MT-TS1", "MT-CO1", "MT-ATP6", "MT-CYB"}

# 已知隐性遗传感官基因
_RECESSIVE_GENES = {
    "GJB2", "GJB6", "MYO7A", "CDH23", "OTOF", "CABP2", "SLC26A4", "KCNQ4",
    "RHO", "PDE6A", "ABCA4", "EYS", "GJA8", "CRYAA", "OCA2", "TYRP1",
    "TAS1R2", "TAS1R3", "TAS1R1", "SCNN1A", "SCNN1B", "SCNN1G",
    "OTOP1", "OTOP2", "OTOP3",
    "TRPV1", "TRPM8", "PIEZO1", "PIEZO2", "NGF", "NTRK1", "GRPR", "IL31RA",
    "CNGA2", "ADCY3",
}

# 已知显性遗传感官基因
_DOMINANT_GENES = {
    "KCNQ4",
    "GJA8",
    "CRYAA",
    "ABCA4",
    "SCN9A", "SCN10A", "SCN11A",
}

# OR 基因视为隐性（纯合 LoF 才表型）
_OR_RECESSIVE = True


class InheritanceMatcher:
    """遗传模式匹配器."""

    def __init__(self) -> None:
        """初始化匹配器."""
        pass

    def detect_pattern(self, gene_symbol: str, sex: Sex) -> InheritancePattern:
        """检测基因的预期遗传模式.

        Args:
            gene_symbol: 基因符号。
            sex: 样本性别。

        Returns:
            遗传模式。
        """
        if gene_symbol in _MITOCHONDRIAL_GENES:
            return "线粒体"
        if gene_symbol in _X_LINKED_GENES:
            return "X-连锁"
        if gene_symbol in _DOMINANT_GENES:
            return "显性"
        if gene_symbol.startswith("OR"):
            return "隐性纯合"
        if gene_symbol in _RECESSIVE_GENES:
            return "隐性纯合"
        return "未知"

    def zygosity_match(
        self, variant: Variant, sex: Sex, expected_pattern: str
    ) -> bool:
        """判断变异基因型是否与预期遗传模式匹配.

        Args:
            variant: 变异记录。
            sex: 样本性别。
            expected_pattern: 预期遗传模式。

        Returns:
            True 表示基因型与遗传模式匹配。
        """
        if expected_pattern == "线粒体":
            return variant.chrom == "MT"

        if expected_pattern == "X-连锁":
            if variant.chrom not in ("X", "chrX"):
                return False
            if sex == "M":
                return variant.is_hemizygous or variant.is_homozygous
            else:
                return variant.is_homozygous or variant.is_heterozygous

        if expected_pattern == "显性":
            return variant.is_heterozygous or variant.is_homozygous

        if expected_pattern in ("隐性纯合", "隐性复合杂合"):
            return variant.is_homozygous

        # 未知模式：只要有变异即视为可能匹配
        return True

    def explain_pattern(self, gene_symbol: str, sex: Sex) -> str:
        """生成遗传模式的说明文本.

        Args:
            gene_symbol: 基因符号。
            sex: 样本性别。

        Returns:
            遗传模式中文说明。
        """
        pattern = self.detect_pattern(gene_symbol, sex)
        explanations = {
            "线粒体": "线粒体母系遗传",
            "X-连锁": f"X-连锁遗传（{'男性半合子或女性纯合才表型' if sex == 'M' else '女性纯合或杂合均可表型'}）",
            "显性": "常染色体显性遗传，杂合即可表型",
            "隐性纯合": "常染色体隐性遗传，通常需纯合才表型",
            "隐性复合杂合": "常染色体隐性遗传，复合杂合可表型",
            "未知": "遗传模式未知",
        }
        return explanations.get(pattern, "遗传模式未知")
