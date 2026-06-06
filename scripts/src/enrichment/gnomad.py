"""gnomAD 查询客户端.

查询基因或变异的人群等位基因频率。
"""

from typing import Any, Dict

from src.enrichment.cache import CacheManager
from src.enrichment.client_base import AsyncApiClient


class GnomADClient(AsyncApiClient):
    """gnomAD GraphQL API 客户端."""

    def __init__(self, cache: CacheManager, rate_limit: int = 10, timeout: int = 30) -> None:
        """初始化 gnomAD 客户端.

        Args:
            cache: 缓存管理器。
            rate_limit: 每秒请求限制。
            timeout: HTTP 超时。
        """
        super().__init__(
            api_name="gnomad",
            base_url="https://gnomad.broadinstitute.org",
            cache=cache,
            rate_limit=rate_limit,
            timeout=timeout,
        )

    async def _fetch(self, key: str) -> Dict[str, Any]:
        """查询 gnomAD 基因信息.

        Args:
            key: 基因符号。

        Returns:
            gnomAD 数据字典。
        """
        session = await self._get_session()

        # 使用 gnomAD GraphQL API
        query = """
        query GeneInfo($geneSymbol: String!) {
            gene(gene_symbol: $geneSymbol, reference_genome: GRCh38) {
                gene_id
                gene_symbol
                pLI
                oe_lof
                oe_lof_upper
                gnomad_constraint {
                    pLI
                    oe_lof
                    oe_lof_upper
                }
            }
        }
        """
        variables = {"geneSymbol": key}

        async with session.post(
            f"{self.base_url}/api/",
            json={"query": query, "variables": variables},
            headers={"Content-Type": "application/json"},
        ) as response:
            response.raise_for_status()
            data = await response.json()

        if data.get("errors"):
            return {
                "found": False,
                "gene": key,
                "errors": data["errors"],
            }

        gene_data = data.get("data", {}).get("gene")
        if not gene_data:
            return {"found": False, "gene": key}

        return {
            "found": True,
            "gene": key,
            "gene_id": gene_data.get("gene_id"),
            "pLI": gene_data.get("pLI"),
            "oe_lof": gene_data.get("oe_lof"),
            "oe_lof_upper": gene_data.get("oe_lof_upper"),
            "constraint": gene_data.get("gnomad_constraint"),
        }

    async def query_gene(self, gene_symbol: str) -> Dict[str, Any]:
        """公开接口：查询基因约束信息.

        Args:
            gene_symbol: 基因符号。

        Returns:
            gnomAD 数据字典。
        """
        return await self.query(gene_symbol)
