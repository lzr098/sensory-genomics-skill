"""GTEx 查询客户端.

查询基因在感官组织中的表达水平。
"""

from typing import Any, Dict

from src.enrichment.cache import CacheManager
from src.enrichment.client_base import AsyncApiClient


class GTExClient(AsyncApiClient):
    """GTEx REST API 客户端."""

    def __init__(self, cache: CacheManager, rate_limit: int = 10, timeout: int = 30) -> None:
        """初始化 GTEx 客户端.

        Args:
            cache: 缓存管理器。
            rate_limit: 每秒请求限制。
            timeout: HTTP 超时。
        """
        super().__init__(
            api_name="gtex",
            base_url="https://gtexportal.org/api/v2",
            cache=cache,
            rate_limit=rate_limit,
            timeout=timeout,
        )

    async def _fetch(self, key: str) -> Dict[str, Any]:
        """查询 GTEx 基因表达.

        Args:
            key: 基因符号或 GENCODE ID。

        Returns:
            GTEx 表达数据字典。
        """
        session = await self._get_session()

        # 查询基因在各组织中的表达
        url = (
            f"{self.base_url}/expression/geneExpression"
            f"?gencodeId={key}&datasetId=gtex_v8"
        )

        async with session.get(url) as response:
            response.raise_for_status()
            data = await response.json()

        gene_expression = data.get("geneExpression", [])
        if not gene_expression:
            # 尝试用 gene symbol
            url = (
                f"{self.base_url}/expression/geneExpression"
                f"?geneId={key}&datasetId=gtex_v8"
            )
            async with session.get(url) as response:
                response.raise_for_status()
                data = await response.json()
                gene_expression = data.get("geneExpression", [])

        if not gene_expression:
            return {"found": False, "gene": key}

        # 筛选感官相关组织
        sensory_tissues = {
            "Brain", "Nerve", "Skin", "Tongue", "Eye", "Ear",
            "Whole Blood",
        }

        tissue_expressions = []
        for expr in gene_expression:
            tissue = expr.get("tissueSiteDetailId", "")
            if any(st.lower() in tissue.lower() for st in sensory_tissues):
                tissue_expressions.append({
                    "tissue": tissue,
                    "tissue_site": expr.get("tissueSiteDetail", ""),
                    "median_tpm": expr.get("median", 0),
                    "unit": expr.get("unit", "TPM"),
                })

        return {
            "found": True,
            "gene": key,
            "tissue_expressions": tissue_expressions,
        }

    async def query_gene(self, gene_symbol: str) -> Dict[str, Any]:
        """公开接口：查询基因在感官组织中的表达.

        Args:
            gene_symbol: 基因符号。

        Returns:
            GTEx 表达数据字典。
        """
        return await self.query(gene_symbol)
