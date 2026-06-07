"""VCF 性别推断工具.

通过分析 chrX 基因型模式和 chrY 存在性推断样本性别：
- 男性：chrX 以纯合为主（hemizygous 编码为 1/1），chrY 有变异
- 女性：chrX 杂合比例正常（~30-40%），chrY 无变异
"""

import subprocess
from typing import Optional


def detect_sex(vcf_path: str, sample_size: int = 50000) -> Optional[str]:
    """从 VCF 的 chrX/chrY 基因型模式推断样本性别.

    Args:
        vcf_path: VCF 文件路径。
        sample_size: 用于推断的 chrX 变异采样数量。

    Returns:
        "M"（男性）、"F"（女性），或 None（无法推断）。
    """
    try:
        # 1. 检查 chrY 是否有变异（男性标志）
        chr_y_count = _count_chrom_variants(vcf_path, "chrY")
        if chr_y_count is None:
            return None  # VCF 不包含 chrY

        # 2. 采样 chrX 的基因型分布
        chr_x_gts = _sample_chrx_genotypes(vcf_path, sample_size)
        if not chr_x_gts:
            return None  # chrX 无数据

        # 3. 判断逻辑
        hom_count = sum(1 for gt in chr_x_gts if gt == "1/1")
        het_count = sum(1 for gt in chr_x_gts if gt == "0/1" or gt == "1|0" or gt == "0|1")
        total = len(chr_x_gts)
        het_ratio = het_count / total if total > 0 else 0

        # 男性 chrX：绝大部分为 hemizygous (编码为 1/1)，杂合 < 15%
        # 女性 chrX：杂合比例约 30-50%
        if chr_y_count > 10 and het_ratio < 0.15:
            return "M"
        elif chr_y_count < 5 and het_ratio > 0.20:
            return "F"
        elif chr_y_count > 10:
            # 有 chrY 但 chrX 数据不足以判断 → 默认为男性
            return "M"

        return None
    except Exception:
        return None


def _count_chrom_variants(vcf_path: str, chrom: str) -> Optional[int]:
    """统计指定染色体上的变异数量."""
    try:
        result = subprocess.run(
            ["bcftools", "view", "-H", vcf_path, chrom],
            capture_output=True, text=True, timeout=60,
        )
        if result.returncode != 0:
            # Try without chrom filter (some BCFtools versions handle region args differently)
            result = subprocess.run(
                ["bcftools", "view", "-H", vcf_path],
                capture_output=True, text=True, timeout=60,
            )
            if result.returncode != 0:
                return None
            count = sum(1 for line in result.stdout.splitlines()
                       if line.strip() and line.startswith(chrom))
            return count
        return sum(1 for line in result.stdout.splitlines() if line.strip())
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return None


def _sample_chrx_genotypes(vcf_path: str, max_variants: int) -> list:
    """采样 chrX 变异并提取基因型 (使用 bcftools view)."""
    try:
        # Use bcftools view -H (more reliable than query)
        result = subprocess.run(
            ["bcftools", "view", "-H", vcf_path, "chrX"],
            capture_output=True, text=True, timeout=60,
        )
        if result.returncode != 0:
            # Fallback: unfiltered view + grep
            return _sample_chrx_genotypes_fallback(vcf_path, max_variants)

        gts = []
        for line in result.stdout.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 10:
                continue
            gt_field = parts[9].split(":")[0]
            if gt_field and gt_field != ".":
                gt = gt_field.replace("|", "/")
                if gt in ("0/0", "0/1", "1/1"):
                    gts.append(gt)
            if len(gts) >= max_variants:
                break
        return gts
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return _sample_chrx_genotypes_fallback(vcf_path, max_variants)


def _sample_chrx_genotypes_fallback(vcf_path: str, max_variants: int) -> list:
    """回退方案：用 bcftools view 全量 + Python grep chrX."""
    import gzip

    gts = []
    try:
        result = subprocess.run(
            ["bcftools", "view", "-H", vcf_path],
            capture_output=True, text=True, timeout=60,
        )
        if result.returncode != 0:
            return []

        for line in result.stdout.splitlines():
            line = line.strip()
            if not line or not line.startswith("chrX"):
                continue
            parts = line.split("\t")
            if len(parts) < 10:
                continue
            gt_field = parts[9].split(":")[0]
            if gt_field and gt_field != ".":
                gt = gt_field.replace("|", "/")
                if gt in ("0/0", "0/1", "1/1"):
                    gts.append(gt)
            if len(gts) >= max_variants:
                return gts
    except Exception:
        pass
    return gts
