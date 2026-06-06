"""感官基因集加载器.

从 YAML/JSON 数据文件加载五大感官子系统基因集，构建内存索引字典。
"""

import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Set

import yaml

from src.logger import get_logger

logger = get_logger(__name__)


class GeneSetLoader:
    """感官基因集加载器.

    从 data/sensory_gene_sets.yaml 加载基因集，构建:
    - gene_index: gene_symbol -> subsystem
    - subsystem_index: subsystem -> set(gene_symbols)
    """

    def __init__(self, data_dir: Optional[str] = None) -> None:
        """初始化并加载基因集.

        Args:
            data_dir: 数据文件目录，默认查找 ../data/。
        """
        if data_dir is None:
            # 默认位于 src/ 同级目录的 data/
            src_dir = Path(__file__).resolve().parent.parent
            data_dir = src_dir.parent / "data"
        else:
            data_dir = Path(data_dir)

        self.data_dir = Path(data_dir)
        self.gene_index: Dict[str, str] = {}
        self.subsystem_index: Dict[str, Set[str]] = {}
        self.gene_details: Dict[str, Dict[str, str]] = {}

        self._load_gene_sets()

    def _load_gene_sets(self) -> None:
        """加载感官基因集 YAML 文件."""
        gene_sets_path = self.data_dir / "sensory_gene_sets.yaml"
        if not gene_sets_path.exists():
            logger.error("sensory_gene_sets.yaml not found at %s", gene_sets_path)
            return

        with open(gene_sets_path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)

        if not data or "subsystems" not in data:
            logger.error("Invalid sensory_gene_sets.yaml format")
            return

        for subsystem_entry in data["subsystems"]:
            subsystem_name = subsystem_entry.get("name", "")
            genes = subsystem_entry.get("genes", [])
            self.subsystem_index[subsystem_name] = set()
            for gene_entry in genes:
                if isinstance(gene_entry, str):
                    gene_symbol = gene_entry
                    details: Dict[str, str] = {}
                elif isinstance(gene_entry, dict):
                    gene_symbol = gene_entry.get("symbol", "")
                    details = {k: str(v) for k, v in gene_entry.items() if k != "symbol"}
                else:
                    continue

                if gene_symbol:
                    self.gene_index[gene_symbol] = subsystem_name
                    self.subsystem_index[subsystem_name].add(gene_symbol)
                    self.gene_details[gene_symbol] = details

        total_genes = len(self.gene_index)
        logger.info(
            "Loaded gene sets: %d genes across %d subsystems",
            total_genes,
            len(self.subsystem_index),
        )

    def is_sensory_gene(self, gene_symbol: str) -> bool:
        """判断基因是否在感官基因集中."""
        return gene_symbol in self.gene_index

    def get_subsystem(self, gene_symbol: str) -> str:
        """获取基因所属感官子系统.

        Returns:
            子系统名称，若不在基因集中返回空字符串。
        """
        return self.gene_index.get(gene_symbol, "")

    def get_genes_for_subsystem(self, subsystem: str) -> List[str]:
        """获取指定子系统的全部基因列表."""
        return sorted(self.subsystem_index.get(subsystem, set()))

    def get_all_genes(self) -> Set[str]:
        """获取全部感官基因符号集合."""
        return set(self.gene_index.keys())

    def get_gene_detail(self, gene_symbol: str, key: str, default: str = "") -> str:
        """获取基因的额外属性."""
        return self.gene_details.get(gene_symbol, {}).get(key, default)

    def get_gene_function(self, gene_symbol: str) -> str:
        """获取基因在感官系统中的功能描述（中文）."""
        return self.get_gene_detail(gene_symbol, "function", "")
