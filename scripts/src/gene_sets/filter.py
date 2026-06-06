"""感官基因精确筛选.

基于 VEP 返回的 gene_symbol 筛选落在感官基因集上的变异。
"""

from typing import Dict, List

from src.gene_sets.loader import GeneSetLoader
from src.logger import get_logger
from src.models import Variant

logger = get_logger(__name__)


class SensoryFilter:
    """感官基因精确筛选器."""

    def __init__(self, gene_sets: GeneSetLoader) -> None:
        """初始化筛选器.

        Args:
            gene_sets: 已加载的感官基因集。
        """
        self.gene_sets = gene_sets

    def filter_variants(self, variants: List[Variant]) -> List[Variant]:
        """筛选落在感官基因集上的变异.

        支持多基因归属：
        1. 优先使用 gene_symbols（多基因列表）——用于基因重叠区域
        2. fallback 到逗号分隔的 gene_symbol

        Args:
            variants: 已注释的变异列表（含 gene_symbol/gene_symbols）。

        Returns:
            仅包含感官基因变异的子列表。
        """
        filtered = []
        for variant in variants:
            # 优先使用 gene_symbols（多基因列表）
            genes = []
            if variant.gene_symbols:
                genes = variant.gene_symbols
            elif variant.gene_symbol:
                genes = [g.strip() for g in variant.gene_symbol.split(",")]

            seen = set()
            for gene in genes:
                if gene in seen:
                    continue
                seen.add(gene)
                if self.gene_sets.is_sensory_gene(gene):
                    # 为每个感官基因创建一个 variant 副本
                    v = variant.model_copy(update={"gene_symbol": gene})
                    filtered.append(v)

        logger.info(
            "Sensory filter: %d / %d variants retained (%d unique variants x genes)",
            len(filtered),
            len(variants),
            len(set(v.vcf_id for v in filtered)),
        )
        return filtered

    def group_by_subsystem(self, variants: List[Variant]) -> Dict[str, List[Variant]]:
        """按感官子系统分组变异.

        Args:
            variants: 已筛选的感官基因变异列表。

        Returns:
            子系统 -> 变异列表的字典。
        """
        groups: Dict[str, List[Variant]] = {}
        for variant in variants:
            subsystem = self.gene_sets.get_subsystem(variant.gene_symbol)
            if subsystem:
                groups.setdefault(subsystem, []).append(variant)
        return groups

    def group_by_gene(self, variants: List[Variant]) -> Dict[str, List[Variant]]:
        """按基因符号分组变异.

        Args:
            variants: 已筛选的感官基因变异列表。

        Returns:
            基因符号 -> 变异列表的字典。
        """
        groups: Dict[str, List[Variant]] = {}
        for variant in variants:
            gene = variant.gene_symbol
            if gene:
                groups.setdefault(gene, []).append(variant)
        return groups
