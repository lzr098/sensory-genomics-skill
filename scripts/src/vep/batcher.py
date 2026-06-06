"""变异批量切分器.

将变异列表切分为适配 VEP REST API limit 的批次。
"""

from typing import Iterator, List

from src.models import Variant


class VariantBatcher:
    """变异批量切分器."""

    def __init__(self, batch_size: int = 200) -> None:
        """初始化切分器.

        Args:
            batch_size: 每批最大变异数量。
        """
        self.batch_size = batch_size

    def split(self, variants: List[Variant]) -> Iterator[List[Variant]]:
        """将变异列表切分为多个批次.

        Args:
            variants: 待切分的变异列表。

        Yields:
            每批变异列表。
        """
        for i in range(0, len(variants), self.batch_size):
            yield variants[i : i + self.batch_size]

    @staticmethod
    def chunk_list(items: List[Variant], n: int) -> Iterator[List[Variant]]:
        """通用列表切分工具.

        Args:
            items: 待切分列表。
            n: 每批大小。

        Yields:
            每批子列表。
        """
        for i in range(0, len(items), n):
            yield items[i : i + n]
