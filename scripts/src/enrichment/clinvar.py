"""ClinVar 查询客户端.

通过 NCBI E-utilities 查询 ClinVar 记录。
"""

from typing import Any, Dict

from src.enrichment.cache import CacheManager
from src.enrichment.client_base import AsyncApiClient


class ClinVarClient(AsyncApiClient):
    """ClinVar E-utilities 客户端."""

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
        """查询 ClinVar 记录.

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
            return {"found": False, "gene": key}

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
        }

    async def query_gene(self, gene_symbol: str) -> Dict[str, Any]:
        """公开接口：查询 ClinVar 记录.

        Args:
            gene_symbol: 基因符号。

        Returns:
            ClinVar 数据字典。
        """
        return await self.query(gene_symbol)
