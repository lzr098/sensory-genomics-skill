"""gnomAD 查询客户端.

双模式:
    基因模式 (query_gene): pLI, oe_lof 约束指标
    变异模式 (query_variant): 人群等位基因频率 (AF, hom count, EAS AF)

GraphQL API, trust_env=False.
"""

from typing import Any, Dict, Optional

from src.enrichment.cache import CacheManager
from src.enrichment.client_base import AsyncApiClient


class GnomADClient(AsyncApiClient):
    """gnomAD API 客户端 (GraphQL).

    trust_env=False: 不走系统代理（gnomAD 境外直连可能有 Clash 干扰）。
    """

    def __init__(self, cache: CacheManager, rate_limit: int = 10, timeout: int = 30) -> None:
        """初始化 gnomAD 客户端."""
        super().__init__(
            api_name="gnomad",
            base_url="https://gnomad.broadinstitute.org",
            cache=cache,
            rate_limit=rate_limit,
            timeout=timeout,
            trust_env=False,
        )

    def _variant_id(self, chrom: str, pos: int, ref: str, alt: str) -> str:
        """构建 gnomAD variant_id: {chrom}-{pos}-{ref}-{alt}."""
        return f"{chrom}-{pos}-{ref}-{alt}"

    async def _fetch(self, key: str) -> Dict[str, Any]:
        """默认: 基因约束查询."""
        return await self._query_gene_graphql(key)

    # ── 基因模式 ──

    async def _query_gene_graphql(self, gene_symbol: str) -> Dict[str, Any]:
        """GraphQL 基因约束查询."""
        session = await self._get_session()
        query = """
        query GeneConstraint($geneSymbol: String!) {
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
        variables = {"geneSymbol": gene_symbol}

        try:
            async with session.post(
                f"{self.base_url}/api/",
                json={"query": query, "variables": variables},
                headers={"Content-Type": "application/json"},
            ) as response:
                response.raise_for_status()
                data = await response.json()
        except Exception:
            return {"found": False, "gene": gene_symbol}

        if data.get("errors"):
            return {"found": False, "gene": gene_symbol, "errors": data["errors"]}

        gene_data = data.get("data", {}).get("gene")
        if not gene_data:
            return {"found": False, "gene": gene_symbol}

        return {
            "found": True,
            "gene": gene_symbol,
            "gene_id": gene_data.get("gene_id"),
            "pLI": gene_data.get("pLI") or gene_data.get("gnomad_constraint", {}).get("pLI"),
            "oe_lof": gene_data.get("oe_lof") or gene_data.get("gnomad_constraint", {}).get("oe_lof"),
            "oe_lof_upper": gene_data.get("oe_lof_upper") or gene_data.get("gnomad_constraint", {}).get("oe_lof_upper"),
            "constraint": gene_data.get("gnomad_constraint"),
            "query_type": "gene",
        }

    # ── 变异模式 ──

    async def _query_variant_graphql(
        self, chrom: str, pos: int, ref: str, alt: str,
        dataset: str = "gnomad_r4",
    ) -> Dict[str, Any]:
        """GraphQL 变异性查询（频率数据）."""
        session = await self._get_session()
        vid = self._variant_id(chrom, pos, ref, alt)

        query = """
        query VariantFreq($variantId: String!, $datasetId: DatasetId!) {
            variant(variantId: $variantId, dataset: $datasetId) {
                variantId
                reference_genome
                rsids
                genome {
                    ac
                    an
                    af
                    homozygote_count
                }
                exome {
                    ac
                    an
                    af
                    homozygote_count
                }
                populations {
                    id
                    ac
                    an
                    af
                    homozygote_count
                }
            }
        }
        """
        variables = {"variantId": vid, "datasetId": dataset}

        try:
            async with session.post(
                f"{self.base_url}/api/",
                json={"query": query, "variables": variables},
                headers={"Content-Type": "application/json"},
            ) as response:
                response.raise_for_status()
                data = await response.json()
        except Exception:
            return {"found": False, "variantId": vid}

        if data.get("errors"):
            return {"found": False, "variantId": vid, "errors": data["errors"]}

        variant_data = data.get("data", {}).get("variant")
        if not variant_data:
            return {"found": False, "variantId": vid}

        # 提取 EAS 人群频率
        pop_data = {}
        for pop in variant_data.get("populations", []):
            pop_id = pop.get("id", "")
            pop_data[pop_id] = {
                "ac": pop.get("ac"),
                "an": pop.get("an"),
                "af": pop.get("af"),
                "homozygote_count": pop.get("homozygote_count"),
            }

        genome = variant_data.get("genome") or {}
        exome = variant_data.get("exome") or {}
        eas = pop_data.get("eas", pop_data.get("EAS", {}))

        return {
            "found": True,
            "variantId": vid,
            "rsids": variant_data.get("rsids", []),
            "genome_af": genome.get("af"),
            "genome_ac": genome.get("ac"),
            "genome_an": genome.get("an"),
            "genome_hom": genome.get("homozygote_count"),
            "exome_af": exome.get("af"),
            "exome_ac": exome.get("ac"),
            "exome_an": exome.get("an"),
            "exome_hom": exome.get("homozygote_count"),
            "eas_af": eas.get("af"),
            "eas_ac": eas.get("ac"),
            "eas_an": eas.get("an"),
            "eas_hom": eas.get("homozygote_count"),
            "all_populations": pop_data,
            "query_type": "variant",
        }

    # ── 公开接口 ──

    async def query_gene(self, gene_symbol: str) -> Dict[str, Any]:
        """基因模式：查询约束指标 (pLI, oe_lof)."""
        return await self.query(gene_symbol)

    async def query_variant(
        self, chrom: str, pos: int, ref: str, alt: str,
        dataset: str = "gnomad_r4",
    ) -> Dict[str, Any]:
        """变异模式：查询人群等位基因频率.

        Args:
            chrom: 染色体 (如 "1", "chr1" 均接受).
            pos: 1-based 位置.
            ref: 参考等位基因.
            alt: 替代等位基因.
            dataset: gnomad_r4 (GRCh38) 或 gnomad_r2_1 (GRCh37).

        Returns:
            AF, hom count, EAS AF 等频率字典.
        """
        # 去 chr 前缀
        chrom = chrom.replace("chr", "").replace("CHR", "")
        return await self._query_variant_graphql(chrom, pos, ref, alt, dataset)
