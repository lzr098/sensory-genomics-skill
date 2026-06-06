"""感官基因坐标预计算脚本.

离线生成精确的基因组坐标文件，用于 pipeline Phase 0 快速过滤：
1. 从 Ensembl REST API 获取所有感官基因的 canonical 转录本坐标
2. 生成 exon/CDS/UTR BED 文件（比基因全区域更精确）
3. 生成 gene -> canonical_transcript_id 映射 JSON
4. 生成基因功能域坐标（用于后续功能域注释）

Usage:
    python3 -m scripts.precompute --output-dir assets/data/

输出文件:
    - sensory_genes_exons.bed      : 精确外显子坐标
    - sensory_genes_cds.bed        : 精确 CDS 坐标
    - gene_transcript_map.json     : gene -> canonical_transcript_id 映射

注意: Ensembl overlap API 不返回 UTR 坐标，如需 UTR 需通过其他 endpoint 获取。
"""

import argparse
import asyncio
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import aiohttp

from src.gene_sets.loader import GeneSetLoader
from src.logger import get_logger

logger = get_logger(__name__)

ENSEMBL_API = "https://rest.ensembl.org"
BATCH_SIZE = 100  # Ensembl POST lookup batch size
RATE_LIMIT = 15   # requests per second


class EnsemblPrecomputer:
    """Ensembl 坐标预计算器."""

    def __init__(self, output_dir: Path, rate_limit: int = RATE_LIMIT, proxy: Optional[str] = None) -> None:
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.semaphore = asyncio.Semaphore(rate_limit)
        self.session: Optional[aiohttp.ClientSession] = None
        self.proxy = proxy or os.environ.get("HTTP_PROXY") or os.environ.get("http_proxy")
        self._cache: Dict[str, Any] = {}  # 内存缓存减少重复请求

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

    async def _post_batch_lookup(self, symbols: List[str]) -> Dict[str, Any]:
        """批量查询基因符号 -> Ensembl ID."""
        url = f"{ENSEMBL_API}/lookup/symbol/homo_sapiens"
        headers = {"Content-Type": "application/json"}
        payload = {"symbols": symbols}

        async with self.semaphore:
            try:
                session = await self._get_session()
                async with session.post(url, headers=headers, json=payload) as resp:
                    if resp.status == 200:
                        return await resp.json()
                    logger.warning("Lookup batch HTTP %d for %d symbols", resp.status, len(symbols))
                    return {}
            except Exception as exc:
                logger.error("Lookup batch error: %s", exc)
                return {}

    async def _get_transcripts(self, gene_id: str, max_retries: int = 2) -> List[Dict[str, Any]]:
        """获取基因的所有转录本信息（含 exon、CDS 坐标）.

        内置指数退避重试，应对 API 偶发超时。
        """
        if gene_id in self._cache:
            return self._cache[gene_id]

        url = (
            f"{ENSEMBL_API}/overlap/id/{gene_id}"
            f"?feature=transcript;feature=exon;feature=cds"
            f";content-type=application/json"
        )

        for attempt in range(max_retries + 1):
            async with self.semaphore:
                try:
                    session = await self._get_session()
                    async with session.get(url) as resp:
                        if resp.status == 200:
                            data = await resp.json()
                            self._cache[gene_id] = data
                            return data
                        logger.warning("Transcript fetch HTTP %d for %s", resp.status, gene_id)
                except asyncio.TimeoutError:
                    logger.warning("Transcript fetch timeout for %s (attempt %d)", gene_id, attempt + 1)
                except Exception as exc:
                    logger.error("Transcript fetch error for %s: %s", gene_id, exc)

            if attempt < max_retries:
                wait = 2 ** attempt
                logger.info("Retrying %s in %ds...", gene_id, wait)
                await asyncio.sleep(wait)

        return []

    async def precompute(self, gene_symbols: List[str]) -> None:
        """执行预计算."""
        logger.info("Starting precompute for %d sensory genes", len(gene_symbols))

        # Stage 1: 批量查询基因符号 -> Ensembl gene ID
        gene_id_map: Dict[str, str] = {}  # symbol -> ensembl_gene_id
        canonical_map: Dict[str, str] = {}  # symbol -> canonical_transcript_id

        for i in range(0, len(gene_symbols), BATCH_SIZE):
            batch = gene_symbols[i:i + BATCH_SIZE]
            result = await self._post_batch_lookup(batch)
            for symbol, info in result.items():
                if isinstance(info, dict) and "id" in info:
                    gene_id_map[symbol] = info["id"]
                    canonical_map[symbol] = info.get("canonical_transcript", "")
            logger.info("Lookup batch %d/%d done", i // BATCH_SIZE + 1, (len(gene_symbols) - 1) // BATCH_SIZE + 1)

        logger.info("Resolved %d/%d genes to Ensembl IDs", len(gene_id_map), len(gene_symbols))

        # Stage 2: 并发获取每个基因的转录本/外显子坐标（分批并发，避免压垮 API）
        exon_bed: List[Tuple[str, int, int, str, str]] = []  # chrom, start, end, gene, transcript
        cds_bed: List[Tuple[str, int, int, str, str]] = []
        failed_genes: List[str] = []

        # 分批并发：每批 20 个基因，避免一次性 300+ 请求压垮代理和 Ensembl API
        CONCURRENT_BATCH = 20
        gene_ids = [(symbol, gene_id_map[symbol]) for symbol in gene_symbols if symbol in gene_id_map]
        logger.info("Fetching transcripts for %d genes in batches of %d...", len(gene_ids), CONCURRENT_BATCH)

        for batch_start in range(0, len(gene_ids), CONCURRENT_BATCH):
            batch = gene_ids[batch_start:batch_start + CONCURRENT_BATCH]
            transcript_tasks = [self._get_transcripts(gid) for _, gid in batch]
            transcript_results = await asyncio.gather(*transcript_tasks, return_exceptions=True)

            for (symbol, gene_id), features in zip(batch, transcript_results):
                if isinstance(features, Exception):
                    logger.debug("Transcript fetch failed for %s: %s", symbol, features)
                    failed_genes.append(symbol)
                    continue
                if not features:
                    failed_genes.append(symbol)
                    continue

                # 收集该基因的所有坐标
                for feat in features:
                    feat_type = feat.get("feature_type", "")
                    chrom = feat.get("seq_region_name", "")
                    start = feat.get("start", 0)
                    end = feat.get("end", 0)
                    transcript_id = feat.get("Parent", feat.get("transcript_id", ""))

                    if not chrom or start <= 0 or end <= 0:
                        continue

                    if feat_type == "exon":
                        exon_bed.append((chrom, start, end, symbol, transcript_id))
                    elif feat_type == "cds":
                        cds_bed.append((chrom, start, end, symbol, transcript_id))

            logger.info("Transcript batch %d/%d done (%d genes)", 
                       batch_start // CONCURRENT_BATCH + 1,
                       (len(gene_ids) - 1) // CONCURRENT_BATCH + 1,
                       len(batch))

        logger.info(
            "Fetched coordinates: %d exons, %d CDS. Failed: %d genes",
            len(exon_bed), len(cds_bed), len(failed_genes),
        )

        # Stage 3: 写输出文件
        self._write_bed(exon_bed, "sensory_genes_exons.bed", cols=5)
        self._write_bed(cds_bed, "sensory_genes_cds.bed", cols=5)

        # gene -> canonical transcript mapping
        with open(self.output_dir / "gene_transcript_map.json", "w") as f:
            json.dump(canonical_map, f, indent=2)

        # failed genes log
        if failed_genes:
            with open(self.output_dir / "precompute_failed_genes.txt", "w") as f:
                f.write("\n".join(failed_genes))

        logger.info("Precompute complete. Files written to %s", self.output_dir)

    def _write_bed(
        self,
        records: List[Tuple],
        filename: str,
        cols: int = 5,
    ) -> None:
        """写 BED 文件."""
        if not records:
            logger.warning("No records for %s", filename)
            return

        # 去重并排序
        seen: Set[str] = set()
        unique_records = []
        for rec in records:
            key = "\t".join(str(x) for x in rec)
            if key not in seen:
                seen.add(key)
                unique_records.append(rec)

        # 按 chrom, start 排序
        unique_records.sort(key=lambda r: (r[0], r[1]))

        output_path = self.output_dir / filename
        with open(output_path, "w") as f:
            for rec in unique_records:
                f.write("\t".join(str(x) for x in rec) + "\n")

        logger.info("Wrote %s: %d unique intervals", filename, len(unique_records))


async def main() -> None:
    parser = argparse.ArgumentParser(description="Precompute sensory gene coordinates from Ensembl")
    parser.add_argument(
        "--output-dir",
        type=str,
        default="assets/data",
        help="Output directory for precomputed files",
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default=None,
        help="Directory containing sensory_gene_sets.yaml",
    )
    parser.add_argument(
        "--subsystems",
        nargs="+",
        default=["vision", "hearing", "olfaction", "taste", "somatosensation"],
        help="Subsystems to precompute",
    )
    args = parser.parse_args()

    # 加载基因集
    loader = GeneSetLoader(data_dir=args.data_dir)
    gene_symbols: List[str] = []
    for subsystem in args.subsystems:
        gene_symbols.extend(loader.get_genes_for_subsystem(subsystem))
    gene_symbols = sorted(set(gene_symbols))

    logger.info("Precomputing coordinates for %d unique genes", len(gene_symbols))

    pc = EnsemblPrecomputer(output_dir=Path(args.output_dir))
    try:
        await pc.precompute(gene_symbols)
    finally:
        await pc.close()


if __name__ == "__main__":
    asyncio.run(main())
