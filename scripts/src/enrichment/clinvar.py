"""ClinVar 本地查询客户端.

仅使用本地 ClinVar VCF 离线查询，不依赖 NCBI API。
VCF 路径: /Users/zhaorongli/.workbuddy/data/clinvar/clinvar.vcf.gz
"""

import subprocess
from typing import Any, Dict, List, Optional
from pathlib import Path


_LOCAL_CLINVAR_VCF = Path("/Users/zhaorongli/.workbuddy/data/clinvar/clinvar.vcf.gz")


def _bcftools_query(args: List[str]) -> str:
    """Run bcftools query and return stdout."""
    cmd = ["bcftools"] + args
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"bcftools failed: {result.stderr}")
    return result.stdout


def query_clinvar_by_gene(gene_symbol: str) -> Optional[Dict[str, Any]]:
    """Query local ClinVar VCF for all records matching a gene symbol.

    Uses bcftools to filter INFO/GENEINFO field for the given gene.
    Returns structured dict, or None if local VCF not found / no records.

    Args:
        gene_symbol: HGNC gene symbol (e.g. "MYO7A").

    Returns:
        Dict with keys: found, gene, count, records, source.
    """
    if not _LOCAL_CLINVAR_VCF.exists():
        return None

    fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%CLNSIG\t%CLNREVSTAT\t%CLNDN\t%GENEINFO\t%ALLELEID\t%CLNHGVS\t%MC\t%ORIGIN\t%RS\n"

    # Filter by GENEINFO containing the gene symbol followed by colon
    # GENEINFO format: SYMBOL:ENTREZ_ID|SYMBOL2:ENTREZ_ID2|...
    filter_expr = f'INFO/GENEINFO ~ "{gene_symbol}:"'

    try:
        stdout = _bcftools_query([
            "query", "-f", fmt,
            "-i", filter_expr,
            str(_LOCAL_CLINVAR_VCF),
        ])
    except RuntimeError:
        return None

    lines = [l for l in stdout.strip().split("\n") if l.strip()]
    if not lines:
        return None

    records = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 14:
            continue
        (
            chrom, pos, vid, ref, alt, clnsig, clnrevstat, clndn,
            geneinfo, alleleid, clnhgvs, mc, origin, rs
        ) = parts

        records.append({
            "uid": vid,
            "title": clnhgvs if clnhgvs != "." else f"{chrom}:{pos} {ref}>{alt}",
            "clinical_significance": {
                "description": clnsig if clnsig != "." else "Unknown",
                "review_status": clnrevstat if clnrevstat != "." else "Unknown",
            },
            "gds": clndn if clndn != "." else "",
            "position": f"{chrom}:{pos}",
            "ref": ref,
            "alt": alt,
            "geneinfo": geneinfo,
            "rs": rs if rs != "." else None,
        })

    return {
        "found": True,
        "gene": gene_symbol,
        "count": len(records),
        "records": records,
        "source": "local_vcf",
    }


def query_clinvar_by_variant(chrom: str, pos: int, ref: str, alt: str) -> Optional[Dict[str, Any]]:
    """Query local ClinVar VCF for a specific variant (position + ref/alt match).

    Returns the first matching record, or None if not found.
    """
    if not _LOCAL_CLINVAR_VCF.exists():
        return None

    fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%CLNSIG\t%CLNREVSTAT\t%CLNDN\t%GENEINFO\t%ALLELEID\t%CLNHGVS\t%MC\t%ORIGIN\t%RS\n"

    chrom_std = chrom.replace("chr", "").replace("CHR", "")

    try:
        stdout = _bcftools_query([
            "query", "-f", fmt,
            "-r", f"{chrom_std}:{pos}-{pos}",
            str(_LOCAL_CLINVAR_VCF),
        ])
    except RuntimeError:
        return None

    lines = [l for l in stdout.strip().split("\n") if l.strip()]
    if not lines:
        return None

    for line in lines:
        parts = line.split("\t")
        if len(parts) < 14:
            continue
        (
            c_chrom, c_pos, vid, c_ref, c_alt, clnsig, clnrevstat, clndn,
            geneinfo, alleleid, clnhgvs, mc, origin, rs
        ) = parts
        # Multi-allelic support
        alts = c_alt.split(",")
        if c_ref == ref and alt in alts:
            return {
                "found": True,
                "chrom": c_chrom,
                "pos": int(c_pos),
                "clnsig": clnsig if clnsig != "." else "Unknown",
                "clnrevstat": clnrevstat if clnrevstat != "." else "Unknown",
                "clndn": clndn if clndn != "." else "",
                "rs": rs if rs != "." else None,
                "source": "local_vcf",
            }

    return None


class ClinVarClient:
    """ClinVar 本地查询客户端（仅离线 VCF，无 API 依赖）.

    不使用 AsyncApiClient（不需要 HTTP session），直接调用 bcftools。
    """

    def __init__(self, cache=None, rate_limit: int = 3, timeout: int = 30) -> None:
        """初始化 ClinVar 客户端.
        
        cache, rate_limit, timeout 保留以兼容原有 API 签名，实际不使用。
        """
        if not _LOCAL_CLINVAR_VCF.exists():
            raise FileNotFoundError(
                f"ClinVar VCF not found: {_LOCAL_CLINVAR_VCF}. "
                f"Please download it first."
            )

    async def query(self, gene_symbol: str) -> Dict[str, Any]:
        """查询 ClinVar 记录（同步 bcftools，包装为 async）"""
        return await self.query_gene(gene_symbol)

    async def query_gene(self, gene_symbol: str) -> Dict[str, Any]:
        """查询基因相关的 ClinVar 记录.

        Args:
            gene_symbol: 基因符号 (HGNC symbol)。

        Returns:
            ClinVar 数据字典。found=False 时表示本地 VCF 中无该基因的记录。
        """
        result = query_clinvar_by_gene(gene_symbol)
        if result is None:
            return {"found": False, "gene": gene_symbol, "source": "local_vcf"}
        return result

    async def close(self) -> None:
        """No-op: 没有 HTTP session 需要关闭."""
        pass
