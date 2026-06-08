"""GTEx 查询客户端.

查询基因在感官组织中的表达水平。
使用本地 gene→ENSG 映射表 (gene_ensg_map.json)，零外部 API 依赖用于 ID 转换。

trust_env=False: 不走系统代理。
"""

import json
from typing import Any, Dict, Optional
from pathlib import Path

from src.enrichment.cache import CacheManager
from src.enrichment.client_base import AsyncApiClient


# 本地 ENSG 映射表路径
_ENSG_MAP_PATH = Path(__file__).resolve().parent.parent.parent.parent / "assets" / "data" / "gene_ensg_map.json"


def _load_ensg_map() -> Dict[str, str]:
    """加载 gene_symbol → ENSG ID 映射表."""
    if _ENSG_MAP_PATH.exists():
        with open(_ENSG_MAP_PATH) as f:
            return json.load(f)
    return {}


class GTExClient(AsyncApiClient):
    """GTEx REST API 客户端.

    使用本地 ENSG 映射表将 gene symbol 转换为 gencodeId，
    避免外部 API 依赖用于 ID 转换。

    trust_env=False: 不走系统代理。
    """

    _API_URLS = ["https://gtexportal.org/api/v2"]

    def __init__(self, cache: CacheManager, rate_limit: int = 10, timeout: int = 30) -> None:
        super().__init__(
            api_name="gtex",
            base_url="https://gtexportal.org/api/v2",
            cache=cache,
            rate_limit=rate_limit,
            timeout=timeout,
            trust_env=False,
        )
        self._ensg_map = _load_ensg_map()

    async def _fetch(self, key: str) -> Dict[str, Any]:
        """查询 GTEx 基因表达.

        策略:
            1. 本地 ENSG 映射 → gencodeId 查询
            2. geneSymbol 查询 (fallback)
            3. geneId 查询 (fallback)
        """
        session = await self._get_session()

        # 策略 1: 本地 ENSG 映射
        ensg = self._ensg_map.get(key)
        if ensg:
            try:
                result = await self._try_query(session, "gencodeId", ensg, key)
                if result.get("found"):
                    return result
            except Exception:
                pass

        # 策略 2: geneSymbol
        try:
            result = await self._try_query(session, "geneSymbol", key, key)
            if result.get("found"):
                return result
        except Exception:
            pass

        # 策略 3: geneId
        try:
            result = await self._try_query(session, "geneId", key, key)
            if result.get("found"):
                return result
        except Exception:
            pass

        return {"found": False, "gene": key}

    async def _try_query(
        self, session, param_name: str, param_value: str, key: str
    ) -> Dict[str, Any]:
        """尝试一种参数组合查询."""
        for api_url in self._API_URLS:
            url = (
                f"{api_url}/expression/geneExpression"
                f"?{param_name}={param_value}&datasetId=gtex_v10"
            )
            async with session.get(url) as response:
                response.raise_for_status()
                data = await response.json()

            gene_expression = data.get("geneExpression", [])
            if gene_expression:
                sensory_tissues = {
                    "Brain", "Nerve", "Skin", "Tongue", "Eye", "Ear",
                    "Whole Blood", "Pituitary", "Spinal cord",
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
                    "source": f"gtex_v10/{param_name}",
                }
        return {"found": False, "gene": key}

    async def query_gene(self, gene_symbol: str) -> Dict[str, Any]:
        """查询基因在感官组织中的表达."""
        return await self.query(gene_symbol)
