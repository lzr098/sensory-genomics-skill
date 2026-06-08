"""Local gnomAD client — queries precompute VEP SQLite DB.

零网络依赖，直接查询 sensory_variants_vep.sqlite 中的 gnomad_af 字段。
"""

import sqlite3
from typing import Any, Dict, Optional
from pathlib import Path


_DEFAULT_DB = Path(__file__).resolve().parent.parent.parent.parent / "assets" / "data" / "sensory_variants_vep.sqlite"


class LocalGnomADClient:
    """本地 gnomAD 变异性查询客户端.

    从 VEP precompute DB 读取 gnomad_af (人群等位基因频率)。
    零网络依赖，作为 GnomADClient.query_variant() 的本地替代。
    """

    def __init__(self, db_path: Optional[str] = None):
        self.db_path = str(db_path or _DEFAULT_DB)
        if not Path(self.db_path).exists():
            raise FileNotFoundError(f"Precompute DB not found: {self.db_path}")

    def query_variant(
        self, chrom: str, pos: int, ref: str, alt: str
    ) -> Dict[str, Any]:
        """查询单个变异的 gnomAD 频率.

        Args:
            chrom: 染色体 (如 "1", "chr1" 均接受).
            pos: 1-based 位置.
            ref: 参考等位基因.
            alt: 替代等位基因.

        Returns:
            {"found": True/False, "gnomad_af": float, ...} 或包含 "warning" 的降级结果.
        """
        chrom = str(chrom).replace("chr", "").replace("CHR", "")
        with sqlite3.connect(self.db_path) as conn:
            row = conn.execute(
                "SELECT gnomad_af FROM variant_vep WHERE chrom=? AND pos=? AND ref=? AND alt=?",
                (chrom, pos, ref, alt),
            ).fetchone()

        if row is None:
            return {
                "found": False,
                "variantId": f"{chrom}-{pos}-{ref}-{alt}",
                "warning": "NOT_IN_PRECOMPUTE_DB",
            }

        return {
            "found": True,
            "variantId": f"{chrom}-{pos}-{ref}-{alt}",
            "gnomad_af": row[0],
            "source": "local_precompute_db",
        }

    async def query_variant_async(
        self, chrom: str, pos: int, ref: str, alt: str
    ) -> Dict[str, Any]:
        """异步包装 (兼容 async enrichment pipeline)."""
        return self.query_variant(chrom, pos, ref, alt)

    async def close(self) -> None:
        """No-op."""
        pass
