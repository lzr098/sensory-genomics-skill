"""ClinVar 查询客户端.

通过本地 ClinVar VCF (优先) 或 NCBI E-utilities (fallback) 查询 ClinVar 记录。

本地 VCF 路径: /Users/zhaorongli/.workbuddy/data/clinvar/clinvar.vcf.gz
"""

import subprocess
from typing import Any, Dict, List, Optional
from pathlib import Path

from src.enrichment.cache import CacheManager
from src.enrichment.client_base import AsyncApiClient


_LOCAL_CLINVAR_VCF = Path("/Users/zhaorongli/.workbuddy/data/clinvar/clinvar.vcf.gz")


def _bcftools_query(args: List[str]) -> str:
    """Run bcftools query and return stdout."""
    cmd = ["bcftools"] + args
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"bcftools failed: {result.stderr}")
    return result.stdout


def _query_clinvar_by_gene(gene_symbol: str) -> Optional[Dict[str, Any]]:
    """Query local ClinVar VCF for all records matching a gene symbol.

    Uses bcftools to filter INFO/GENEINFO field for the given gene.
    Returns structured dict compatible with NCBI API response format.
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

    Returns the first matching record, or None if not found / no local VCF.
    """
    if not _LOCAL_CLINVAR_VCF.exists():
        return None

    fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%CLNSIG\t%CLNREVSTAT\t%CLNDN\t%GENEINFO\t%ALLELEID\t%CLNHGVS\t%MC\t%ORIGIN\t%RS\n"

    try:
        stdout = _bcftools_query([
            "query", "-f", fmt,
            "-r", f"{chrom}:{pos}-{pos}",
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


class ClinVarClient(AsyncApiClient):
    """ClinVar E-utilities 客户端（带本地 VCF 优先查询）."""

    def __init__(self, cache: CacheManager, rate_limit: int = 3, timeout: int = 30) -> None:
        """初始化 ClinVar 客户端.

        Args:
            cache: 缓存管理器。
            rate_limit: 每秒请求限制（NCBI 限制较严，默认 3 req/s）。
            timeout: HTTP 超时。
        """
        super().__init__(
            api_name="clinvar",
            base_url="https://eutils.ncbi.nlm.nih.gov",
            cache=cache,
            rate_limit=rate_limit,
            timeout=timeout,
        )

    async def _fetch(self, key: str) -> Dict[str, Any]:
        """查询 ClinVar 记录 (NCBI API fallback).

        Args:
            key: 基因符号或 rsID。

        Returns:
            ClinVar 数据字典。
        """
        session = await self._get_session()

        # Step 1: esearch 获取 ID 列表
        search_url = (
            f"{self.base_url}/entrez/eutils/esearch.fcgi"
            f"?db=clinvar&term={key}[Gene]&retmode=json&retmax=5"
        )
        async with session.get(search_url) as response:
            response.raise_for_status()
            search_data = await response.json()

        id_list = search_data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return {"found": False, "gene": key, "source": "ncbi_api"}

        # Step 2: esummary 获取摘要
        ids = ",".join(id_list)
        summary_url = (
            f"{self.base_url}/entrez/eutils/esummary.fcgi"
            f"?db=clinvar&id={ids}&retmode=json"
        )
        async with session.get(summary_url) as response:
            response.raise_for_status()
            summary_data = await response.json()

        records = []
        result = summary_data.get("result", {})
        for uid in id_list:
            item = result.get(str(uid), {})
            records.append({
                "uid": uid,
                "title": item.get("title", ""),
                "clinical_significance": item.get("clinical_significance", {}),
                "gds": item.get("gds", ""),
            })

        return {
            "found": True,
            "gene": key,
            "count": len(records),
            "records": records,
            "source": "ncbi_api",
        }

    async def query_gene(self, gene_symbol: str) -> Dict[str, Any]:
        """公开接口：查询 ClinVar 记录.

        优先使用本地 VCF 查询（离线、零延迟），本地未找到时 fallback 到 NCBI API。

        Args:
            gene_symbol: 基因符号。

        Returns:
            ClinVar 数据字典。
        """
        # Try local VCF first
        local_result = _query_clinvar_by_gene(gene_symbol)
        if local_result:
            return local_result

        # Fallback to NCBI API
        return await self.query(gene_symbol)
