"""UniProt 查询客户端.

查询蛋白功能、结构域、拓扑位置等信息。
"""

from typing import Any, Dict

from src.enrichment.cache import CacheManager
from src.enrichment.client_base import AsyncApiClient


class UniProtClient(AsyncApiClient):
    """UniProt REST API 客户端."""

    def __init__(self, cache: CacheManager, rate_limit: int = 10, timeout: int = 30) -> None:
        """初始化 UniProt 客户端.

        Args:
            cache: 缓存管理器。
            rate_limit: 每秒请求限制。
            timeout: HTTP 超时。
        """
        super().__init__(
            api_name="uniprot",
            base_url="https://rest.uniprot.org",
            cache=cache,
            rate_limit=rate_limit,
            timeout=timeout,
        )

    async def _fetch(self, key: str) -> Dict[str, Any]:
        """查询 UniProt 蛋白信息.

        Args:
            key: 基因符号或 UniProt ID。

        Returns:
            蛋白信息字典。
        """
        session = await self._get_session()
        url = f"{self.base_url}/uniprotkb/search?query=gene:{key}+AND+organism_id:9606&format=json&size=1"

        async with session.get(url) as response:
            response.raise_for_status()
            data = await response.json()

        results = data.get("results", [])
        if not results:
            return {"found": False, "gene": key}

        entry = results[0]
        protein_info = {
            "found": True,
            "gene": key,
            "accession": entry.get("primaryAccession", ""),
            "protein_name": entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", ""),
            "functions": [],
            "domains": [],
            "topology": [],
            "subcellular_locations": [],
        }

        # 提取功能注释
        comments = entry.get("comments", [])
        for comment in comments:
            comment_type = comment.get("commentType", "")
            if comment_type == "FUNCTION":
                texts = comment.get("texts", [])
                for text in texts:
                    protein_info["functions"].append(text.get("value", ""))
            elif comment_type == "SUBCELLULAR LOCATION":
                locations = comment.get("subcellularLocations", [])
                for loc in locations:
                    loc_val = loc.get("location", {}).get("value", "")
                    if loc_val:
                        protein_info["subcellular_locations"].append(loc_val)
            elif comment_type == "TOPOLOGY":
                topo = comment.get("topology", {}).get("value", "")
                if topo:
                    protein_info["topology"].append(topo)

        # 提取结构域
        features = entry.get("features", [])
        for feature in features:
            feature_type = feature.get("type", "")
            if feature_type in ("Domain", "Repeat", "Transmembrane", "Topological domain"):
                domain_info = {
                    "type": feature_type,
                    "description": feature.get("description", ""),
                    "location": feature.get("location", {}),
                }
                protein_info["domains"].append(domain_info)

        return protein_info

    async def query_protein(self, gene_symbol: str) -> Dict[str, Any]:
        """公开接口：查询蛋白信息.

        Args:
            gene_symbol: 基因符号。

        Returns:
            蛋白信息字典。
        """
        return await self.query(gene_symbol)
