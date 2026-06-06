"""VCF 流式解析器.

使用 pysam.VariantFile 流式读取 VCF，提取 GT/DP/GQ/QUAL/AD 等核心字段。
"""

import re
from pathlib import Path
from typing import Iterator, List, Optional

from src.exceptions import VcfError
from src.logger import get_logger
from src.models import Sex, Variant

try:
    import pysam
except ImportError:
    pysam = None  # type: ignore

logger = get_logger(__name__)

# MT 染色体名标准化
_MT_CHROM_NAMES = {"MT", "chrM", "chrMT", "M"}


class VcfParser:
    """VCF 流式解析器.

    支持 .vcf 和 .vcf.gz 格式，流式读取以控制内存占用。
    """

    def __init__(self, vcf_path: str, sex: Sex) -> None:
        """初始化解析器.

        Args:
            vcf_path: VCF 文件路径。
            sex: 样本性别 M/F。
        """
        self.vcf_path = Path(vcf_path)
        self.sex = sex
        self.sample_name = ""

        if pysam is None:
            raise VcfError("pysam is not installed, cannot parse VCF.")

        if not self.vcf_path.exists():
            raise VcfError(f"VCF file not found: {vcf_path}")

    def iter_variants(self) -> Iterator[Variant]:
        """流式迭代 VCF 中的变异记录.

        Yields:
            Variant 对象。
        """
        try:
            with pysam.VariantFile(str(self.vcf_path)) as vcf:
                samples = list(vcf.header.samples)
                if not samples:
                    logger.warning("No samples found in VCF header")
                    return

                self.sample_name = samples[0]
                logger.info("Parsing VCF for sample: %s", self.sample_name)

                for record in vcf:
                    try:
                        variant = self._record_to_variant(record)
                        if variant:
                            yield variant
                    except Exception as exc:
                        chrom = record.chrom if hasattr(record, "chrom") else "?"
                        pos = record.pos if hasattr(record, "pos") else "?"
                        logger.warning(
                            "Skipping malformed record at %s:%s: %s",
                            chrom,
                            pos,
                            exc,
                        )
        except Exception as exc:
            raise VcfError(f"Failed to open or read VCF: {exc}") from exc

    def _record_to_variant(self, record) -> Optional[Variant]:
        """将 pysam VariantRecord 转换为 Variant 模型.

        Args:
            record: pysam VariantRecord。

        Returns:
            Variant 对象，若样本无基因型则返回 None。
        """
        sample = self.sample_name
        call = record.samples.get(sample)
        if call is None:
            return None

        gt = self._format_gt(call.get("GT"))
        if not gt or gt == "./." or gt == ".|." or gt == ".":
            return None

        chrom = self._normalize_chrom(record.chrom)
        pos = int(record.pos)
        ref = str(record.ref) if record.ref else ""

        # 遍历所有 ALT
        alts = record.alts if record.alts else []
        if not alts:
            return None

        alt = str(alts[0])
        # 跳过复杂 SV / 多等位基因（v0.1.0 仅处理 SNV/indel）
        if len(alts) > 1:
            return None

        dp = int(call.get("DP", 0)) if call.get("DP") is not None else 0
        qual = float(record.qual) if record.qual is not None else 0.0
        filter_status = self._format_filter(record.filter)

        return Variant(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            gt=gt,
            dp=dp,
            qual=qual,
            filter_status=filter_status,
        )

    @staticmethod
    def _normalize_chrom(chrom: str) -> str:
        """标准化染色体名称."""
        chrom = str(chrom)
        if chrom.upper() in _MT_CHROM_NAMES:
            return "MT"
        # 去掉 chr 前缀
        if chrom.lower().startswith("chr"):
            chrom = chrom[3:]
        return chrom

    @staticmethod
    def _format_gt(gt) -> str:
        """格式化基因型为字符串."""
        if gt is None:
            return ""
        if isinstance(gt, tuple):
            return "/".join(str(x) if x is not None else "." for x in gt)
        return str(gt)

    @staticmethod
    def _format_filter(filter_val) -> str:
        """格式化 FILTER 字段."""
        if filter_val is None:
            return "PASS"
        if isinstance(filter_val, str):
            return filter_val
        if hasattr(filter_val, "__iter__"):
            return ";".join(str(f) for f in filter_val)
        return str(filter_val)

    def fetch_region(
        self, chrom: str, start: int, end: int
    ) -> List[Variant]:
        """获取指定区域的变异列表.

        Args:
            chrom: 染色体。
            start: 起始位置（0-based，pysam 风格）。
            end: 结束位置。

        Returns:
            该区域内的 Variant 列表。
        """
        variants: List[Variant] = []
        try:
            with pysam.VariantFile(str(self.vcf_path)) as vcf:
                for record in vcf.fetch(chrom, start, end):
                    variant = self._record_to_variant(record)
                    if variant:
                        variants.append(variant)
        except Exception as exc:
            logger.warning("fetch_region failed for %s:%d-%d: %s", chrom, start, end, exc)
        return variants
