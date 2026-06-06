"""BED 文件坐标到基因名的映射器.

当 VEP 注释不可用时，使用 BED 文件中的基因组坐标区间
将变异映射到对应的基因名。

使用 bisect 实现 O(log n) 区间查找（替代原来的 O(n) 线性扫描）。
"""

import bisect
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from src.logger import get_logger

logger = get_logger(__name__)


class BedMapper:
    """BED 文件坐标映射器."""

    def __init__(self, bed_path: Optional[str] = None) -> None:
        """初始化 BED 映射器.

        Args:
            bed_path: BED 文件路径，默认查找 assets/data/sensory_gene_regions.bed。
        """
        # chrom -> (sorted_starts, ends, genes)
        # sorted_starts: 有序的起始位置列表
        # ends, genes: 与 sorted_starts 一一对应的结束位置和基因名
        self.regions: Dict[str, Tuple[List[int], List[int], List[str]]] = {}
        self._load_bed(bed_path)

    def _load_bed(self, bed_path: Optional[str]) -> None:
        """加载 BED 文件并构建 bisect 索引."""
        if bed_path is None:
            candidates = [
                Path(__file__).resolve().parent.parent.parent / "assets" / "data" / "sensory_gene_regions.bed",
                Path(__file__).resolve().parent.parent.parent / "data" / "sensory_gene_regions.bed",
            ]
            for candidate in candidates:
                if candidate.exists():
                    bed_path = str(candidate)
                    break

        if bed_path is None or not Path(bed_path).exists():
            logger.warning("BED file not found, coordinate mapping disabled")
            return

        # 临时存储：chrom -> List[(start, end, gene)]
        raw_regions: Dict[str, List[Tuple[int, int, str]]] = {}

        with open(bed_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 4:
                    chrom, start, end, gene = parts[0], int(parts[1]), int(parts[2]), parts[3]
                    raw_regions.setdefault(chrom, []).append((start, end, gene))

        # 按起始位置排序并构建 bisect 索引
        for chrom, intervals in raw_regions.items():
            intervals.sort(key=lambda x: x[0])
            starts = [iv[0] for iv in intervals]
            ends = [iv[1] for iv in intervals]
            genes = [iv[2] for iv in intervals]
            self.regions[chrom] = (starts, ends, genes)

        total = sum(len(v[0]) for v in self.regions.values())
        logger.info("Loaded BED regions: %d intervals across %d chromosomes (bisect indexed)", total, len(self.regions))

    def lookup(self, chrom: str, pos: int) -> Optional[str]:
        """根据染色体和位置查找基因名（O(log n) bisect）.

        Args:
            chrom: 染色体名（支持 chr 前缀）。
            pos: 基因组位置（1-based）。

        Returns:
            基因名，如果不在任何区间内返回 None。
        """
        # 统一染色体命名
        chrom_key = chrom.replace("chr", "") if chrom.startswith("chr") else chrom
        chrom_key_with_chr = f"chr{chrom_key}" if not chrom.startswith("chr") else chrom

        for key in (chrom, chrom_key, chrom_key_with_chr):
            if key in self.regions:
                starts, ends, genes = self.regions[key]
                # bisect_right 找到第一个 start > pos 的位置
                # 候选区间在该位置之前
                idx = bisect.bisect_right(starts, pos) - 1
                if idx >= 0 and pos <= ends[idx]:
                    return genes[idx]
        return None
