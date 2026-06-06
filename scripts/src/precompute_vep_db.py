"""预计算感官基因区域内常见变异的 VEP 注释.

从 gnomAD 查询每个基因区域的变异列表，过滤 AF>阈值后做 VEP 注释，
构建本地 SQLite 数据库。运行时优先查库，未命中再调 VEP API.

Usage:
    PYTHONPATH=scripts/src python3 scripts/src/precompute_vep_db.py \
        --bed assets/data/sensory_genes_exons.bed \
        --output assets/data/sensory_variants_vep.sqlite \
        --af-threshold 0.001
"""

import argparse
import asyncio
import json
import os
import sqlite3
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import aiohttp

from src.logger import get_logger

logger = get_logger(__name__)

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
VEP_API = "https://rest.ensembl.org"
BATCH_SIZE = 200
GNOMAD_RATE_LIMIT = 2    # gnomAD API 更严格，低并发避免 429
VEP_RATE_LIMIT = 5       # VEP REST API 并发上限


class VariantPrecomputer:
    """gnomAD + VEP 预计算器."""

    def __init__(
        self,
        bed_path: str,
        output_db: str,
        af_threshold: float = 0.001,
        gnomad_rate_limit: int = GNOMAD_RATE_LIMIT,
        vep_rate_limit: int = VEP_RATE_LIMIT,
    ) -> None:
        self.bed_path = Path(bed_path)
        self.output_db = Path(output_db)
        self.af_threshold = af_threshold
        self.gnomad_sem = asyncio.Semaphore(gnomad_rate_limit)
        self.vep_sem = asyncio.Semaphore(vep_rate_limit)
        self.session: Optional[aiohttp.ClientSession] = None
        self.proxy = os.environ.get("HTTP_PROXY") or os.environ.get("http_proxy")

    async def _get_session(self) -> aiohttp.ClientSession:
        if self.session is None or self.session.closed:
            kwargs: Dict[str, Any] = {
                "timeout": aiohttp.ClientTimeout(total=60),
                "trust_env": True,
            }
            if self.proxy:
                kwargs["proxy"] = self.proxy
            self.session = aiohttp.ClientSession(**kwargs)
        return self.session

    async def close(self) -> None:
        if self.session and not self.session.closed:
            await self.session.close()

    def _init_db(self) -> None:
        """初始化 SQLite 数据库."""
        self.output_db.parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(str(self.output_db))
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS variant_vep (
                id TEXT PRIMARY KEY,
                chrom TEXT NOT NULL,
                pos INTEGER NOT NULL,
                ref TEXT NOT NULL,
                alt TEXT NOT NULL,
                gene_symbol TEXT,
                consequence TEXT,
                hgvsc TEXT,
                hgvsp TEXT,
                canonical INTEGER,
                gnomad_af REAL,
                vep_json TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
            """
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_variant ON variant_vep(chrom, pos, ref, alt)"
        )
        conn.commit()
        conn.close()

    def _load_bed_regions(self) -> List[Tuple[str, int, int]]:
        """加载 BED 文件并合并相邻区间."""
        regions: List[Tuple[str, int, int]] = []
        with open(self.bed_path, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    regions.append((chrom, start, end))

        # 按染色体排序并合并重叠/相邻区间（<1kb 视为可合并）
        regions.sort(key=lambda r: (r[0], r[1]))
        merged: List[Tuple[str, int, int]] = []
        for chrom, start, end in regions:
            if merged and merged[-1][0] == chrom and start <= merged[-1][2] + 1000:
                merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
            else:
                merged.append((chrom, start, end))

        logger.info("Loaded %d BED intervals, merged to %d regions", len(regions), len(merged))
        return merged

    async def _query_gnomad_region(
        self, chrom: str, start: int, end: int, max_retries: int = 3
    ) -> List[Dict[str, Any]]:
        """查询 gnomAD 区域内的变异列表.

        内置指数退避重试，处理 429 Rate Limit。
        """
        query = (
            f'query {{ region(chrom: "{chrom}", start: {start}, stop: {end}, '
            f'reference_genome: GRCh38) {{ '
            f'variants(dataset: gnomad_r4) {{ variant_id exome {{ af }} genome {{ af }} }} '
            f'}} }}'
        )

        for attempt in range(max_retries + 1):
            async with self.gnomad_sem:
                try:
                    session = await self._get_session()
                    async with session.post(
                        GNOMAD_API,
                        headers={"Content-Type": "application/json"},
                        json={"query": query},
                    ) as resp:
                        if resp.status == 200:
                            data = await resp.json()
                            variants = data.get("data", {}).get("region", {}).get("variants", [])
                            return variants
                        elif resp.status == 429:
                            logger.warning(
                                "gnomAD 429 for %s:%d-%d (attempt %d/%d)",
                                chrom, start, end, attempt + 1, max_retries + 1,
                            )
                        else:
                            logger.warning("gnomAD HTTP %d for %s:%d-%d", resp.status, chrom, start, end)
                            return []
                except Exception as exc:
                    logger.error("gnomAD query error for %s:%d-%d: %s", chrom, start, end, exc)

            if attempt < max_retries:
                wait = 2 ** attempt + 1  # 2, 5, 9 秒退避
                logger.info("Retrying gnomAD %s:%d-%d in %ds...", chrom, start, end, wait)
                await asyncio.sleep(wait)

        return []

    def _parse_variant_id(self, variant_id: str) -> Optional[Tuple[str, int, str, str]]:
        """解析 gnomAD variant_id: chrom-pos-ref-alt."""
        try:
            parts = variant_id.split("-")
            if len(parts) >= 4:
                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[2]
                alt = parts[3]
                return chrom, pos, ref, alt
        except Exception:
            pass
        return None

    def _filter_by_af(self, variants: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """按 AF 阈值过滤变异."""
        filtered = []
        for v in variants:
            exome_af = (v.get("exome") or {}).get("af") or 0
            genome_af = (v.get("genome") or {}).get("af") or 0
            max_af = max(exome_af, genome_af)
            if max_af >= self.af_threshold:
                v["max_af"] = max_af
                filtered.append(v)
        return filtered

    async def _vep_annotate_batch(
        self, variants: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """对一批变异做 VEP 注释."""
        if not variants:
            return []

        # 构建 VEP POST payload
        payload = []
        for v in variants:
            parsed = self._parse_variant_id(v["variant_id"])
            if parsed:
                chrom, pos, ref, alt = parsed
                payload.append(f"{chrom} {pos} . {ref} {alt} . . .")

        if not payload:
            return []

        url = f"{VEP_API}/vep/human/region"
        headers = {"Content-Type": "application/json"}
        data = {
            "variants": payload,
            "canonical": "1",
            "hgvs": "1",
        }

        async with self.vep_sem:
            try:
                session = await self._get_session()
                async with session.post(url, headers=headers, json=data) as resp:
                    if resp.status == 200:
                        results = await resp.json()
                        return self._merge_vep_results(variants, results)
                    logger.warning("VEP batch HTTP %d", resp.status)
                    return []
            except Exception as exc:
                logger.error("VEP batch error: %s", exc)
                return []

    def _merge_vep_results(
        self, gnomad_variants: List[Dict[str, Any]], vep_results: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """合并 gnomAD 变异和 VEP 结果."""
        merged = []
        for gv, vr in zip(gnomad_variants, vep_results):
            parsed = self._parse_variant_id(gv["variant_id"])
            if not parsed:
                continue

            chrom, pos, ref, alt = parsed
            record = {
                "id": f"{chrom}:{pos}:{ref}:{alt}",
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gnomad_af": gv.get("max_af", 0),
                "vep_json": json.dumps(vr),
            }

            # 提取最佳转录本信息
            tc = vr.get("transcript_consequences", [])
            if tc:
                # 选 canonical + protein_coding
                best = None
                best_score = -1
                for t in tc:
                    score = 0
                    if t.get("canonical") == "1":
                        score += 10
                    if t.get("biotype") == "protein_coding":
                        score += 5
                    if score > best_score:
                        best_score = score
                        best = t
                if best:
                    record["gene_symbol"] = best.get("gene_symbol", "")
                    record["consequence"] = ",".join(best.get("consequence_terms", []))
                    record["hgvsc"] = best.get("hgvsc", "")
                    record["hgvsp"] = best.get("hgvsp", "")
                    record["canonical"] = 1 if best.get("canonical") == "1" else 0

            merged.append(record)
        return merged

    def _save_to_db(self, records: List[Dict[str, Any]]) -> None:
        """保存记录到 SQLite."""
        if not records:
            return

        conn = sqlite3.connect(str(self.output_db))
        cursor = conn.cursor()
        for r in records:
            cursor.execute(
                """
                INSERT OR REPLACE INTO variant_vep
                (id, chrom, pos, ref, alt, gene_symbol, consequence, hgvsc, hgvsp, canonical, gnomad_af, vep_json)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    r["id"],
                    r["chrom"],
                    r["pos"],
                    r["ref"],
                    r["alt"],
                    r.get("gene_symbol", ""),
                    r.get("consequence", ""),
                    r.get("hgvsc", ""),
                    r.get("hgvsp", ""),
                    r.get("canonical", 0),
                    r.get("gnomad_af", 0),
                    r.get("vep_json", "{}"),
                ),
            )
        conn.commit()
        conn.close()
        logger.info("Saved %d records to database", len(records))

    async def precompute(self) -> None:
        """执行预计算."""
        self._init_db()
        regions = self._load_bed_regions()

        total_gnomad = 0
        total_filtered = 0
        total_vep = 0

        for i, (chrom, start, end) in enumerate(regions):
            logger.info("Processing region %d/%d: %s:%d-%d", i + 1, len(regions), chrom, start, end)

            # 1. 查询 gnomAD
            gnomad_variants = await self._query_gnomad_region(chrom, start, end)
            total_gnomad += len(gnomad_variants)

            # 2. AF 过滤
            filtered = self._filter_by_af(gnomad_variants)
            total_filtered += len(filtered)
            logger.info("  gnomAD: %d, AF>=%.4f: %d", len(gnomad_variants), self.af_threshold, len(filtered))

            if not filtered:
                continue

            # 3. VEP 批量注释
            for j in range(0, len(filtered), BATCH_SIZE):
                batch = filtered[j : j + BATCH_SIZE]
                vep_results = await self._vep_annotate_batch(batch)
                if vep_results:
                    self._save_to_db(vep_results)
                    total_vep += len(vep_results)

        logger.info(
            "Precompute complete. gnomAD total: %d, AF filtered: %d, VEP annotated: %d",
            total_gnomad,
            total_filtered,
            total_vep,
        )


async def main() -> None:
    parser = argparse.ArgumentParser(description="Precompute VEP annotations for common variants in sensory gene regions")
    parser.add_argument("--bed", type=str, required=True, help="Path to exon BED file")
    parser.add_argument("--output", type=str, default="assets/data/sensory_variants_vep.sqlite", help="Output SQLite DB path")
    parser.add_argument("--af-threshold", type=float, default=0.001, help="Minimum gnomAD AF to include")
    args = parser.parse_args()

    pc = VariantPrecomputer(
        bed_path=args.bed,
        output_db=args.output,
        af_threshold=args.af_threshold,
    )
    try:
        await pc.precompute()
    finally:
        await pc.close()


if __name__ == "__main__":
    asyncio.run(main())
