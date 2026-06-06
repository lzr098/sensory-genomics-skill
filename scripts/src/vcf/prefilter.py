"""VCF 预过滤规则.

基于 QUAL、DP、FILTER 字段进行质量过滤，减少下游分析负荷。
"""

from src.logger import get_logger
from src.models import Variant

logger = get_logger(__name__)


class Prefilter:
    """预过滤规则.

    默认规则:
    - QUAL >= 30
    - DP >= 10
    - FILTER == PASS（若 pass_only=True）
    """

    def __init__(self, min_qual: int = 30, min_dp: int = 10, pass_only: bool = True) -> None:
        """初始化预过滤器.

        Args:
            min_qual: 最小质量值阈值。
            min_dp: 最小测序深度阈值。
            pass_only: 是否仅保留 FILTER=PASS 的记录。
        """
        self.min_qual = min_qual
        self.min_dp = min_dp
        self.pass_only = pass_only

    def apply(self, variant: Variant) -> bool:
        """应用预过滤规则.

        Args:
            variant: 待过滤的变异记录。

        Returns:
            True 表示通过过滤，False 表示被过滤掉。
        """
        if variant.qual < self.min_qual:
            return False
        if variant.dp < self.min_dp:
            return False
        if self.pass_only and variant.filter_status not in ("PASS", ".", ""):
            return False
        return True

    def apply_with_logging(self, variant: Variant) -> bool:
        """应用预过滤规则并记录被过滤的原因.

        Args:
            variant: 待过滤的变异记录。

        Returns:
            True 表示通过过滤，False 表示被过滤掉。
        """
        if variant.qual < self.min_qual:
            logger.debug(
                "Filtered %s:%d (QUAL=%.1f < %d)",
                variant.chrom,
                variant.pos,
                variant.qual,
                self.min_qual,
            )
            return False
        if variant.dp < self.min_dp:
            logger.debug(
                "Filtered %s:%d (DP=%d < %d)",
                variant.chrom,
                variant.pos,
                variant.dp,
                self.min_dp,
            )
            return False
        if self.pass_only and variant.filter_status not in ("PASS", ".", ""):
            logger.debug(
                "Filtered %s:%d (FILTER=%s != PASS)",
                variant.chrom,
                variant.pos,
                variant.filter_status,
            )
            return False
        return True
