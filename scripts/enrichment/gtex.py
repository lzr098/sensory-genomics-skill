"""GTEx 查询客户端.

v0.10.16: 本地 GTEx SQLite DB 优先，零外部 API 调用。
v0.10.16b: 本地 miss 时自动从 GTEx API 拉取并缓存（lazy loading）。
查询基因在感官组织中的表达水平。
"""

import asyncio
import sys
from typing import Any, Dict
from pathlib import Path

from src.enrichment.cache import CacheManager
from src.enrichment.client_base import AsyncApiClient


def _query_gtex_local(gene_symbol: str) -> Dict[str, float]:
    """查询本地 GTEx SQLite DB，返回 {tissue: median_tpm}."""
    try:
        gtex_local_path = str(Path.home() / ".workbuddy" / "scripts")
        if gtex_local_path not in sys.path:
            sys.path.insert(0, gtex_local_path)
        from gtex_local import query_gtex_local
        return query_gtex_local(gene_symbol, tissues=None)
    except Exception:
        return {}


def _fetch_and_cache_gtex(gene_symbol: str) -> bool:
    """从 GTEx API 拉取并写入本地 SQLite DB（sync，用于 asyncio.to_thread）."""
    try:
        gtex_local_path = str(Path.home() / ".workbuddy" / "scripts")
        if gtex_local_path not in sys.path:
            sys.path.insert(0, gtex_local_path)
        from gtex_local import query_gtex_api_median, insert_gtex_api_results
        results = query_gtex_api_median(gene_symbol)
        if results:
            insert_gtex_api_results(None, gene_symbol, results)
            return True
    except Exception:
        pass
    return False


class GTExClient(AsyncApiClient):
    """GTEx 本地 DB + lazy loading 客户端.

    v0.10.16b: 本地优先，miss 时自动在线补数据。
    """

    def __init__(self, cache: CacheManager, rate_limit: int = 10, timeout: int = 30) -> None:
        super().__init__(
            api_name="gtex",
            base_url="local_db",  # dummy, not used
            cache=cache,
            rate_limit=rate_limit,
            timeout=timeout,
        )

    async def _fetch(self, key: str) -> Dict[str, Any]:
        """查询本地 GTEx SQLite DB，miss 时自动在线拉取缓存."""
        local_data = _query_gtex_local(key)

        # --- v0.10.16b: LAZY LOAD ---
        if not local_data:
            try:
                fetched = await asyncio.to_thread(_fetch_and_cache_gtex, key)
                if fetched:
                    local_data = _query_gtex_local(key)
            except Exception:
                pass

        if not local_data:
            return {"found": False, "gene": key}

        # 筛选感官相关组织
        sensory_tissues = {
            "Brain", "Nerve", "Skin", "Tongue", "Eye", "Ear",
            "Whole Blood",
        }

        tissue_expressions = []
        for tissue, tpm in local_data.items():
            if any(st.lower() in tissue.lower() for st in sensory_tissues):
                tissue_expressions.append({
                    "tissue": tissue,
                    "tissue_site": tissue,
                    "median_tpm": tpm,
                    "unit": "TPM",
                })

        return {
            "found": True,
            "gene": key,
            "tissue_expressions": tissue_expressions,
            "source": "gtex_local_db",
        }

    async def query_gene(self, gene_symbol: str) -> Dict[str, Any]:
        """公开接口：查询基因在感官组织中的表达."""
        return await self.query(gene_symbol)
