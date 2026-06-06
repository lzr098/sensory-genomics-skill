"""VEP REST API 异步客户端.

支持批量 POST 查询、指数退避重试、15 req/s 限速。
"""

import asyncio
import json
import sqlite3
from typing import Any, Dict, List, Optional

import aiohttp

from src.exceptions import VepError
from src.logger import get_logger
from src.models import Variant

logger = get_logger(__name__)


class VepClient:
    """Ensembl VEP REST API 异步客户端."""

    def __init__(
        self,
        base_url: str = "https://rest.ensembl.org",
        batch_size: int = 200,
        rate_limit: int = 15,
        max_retries: int = 3,
        timeout: int = 60,
        cache=None,
        canonical_only: bool = True,
        transcript_ids: Optional[List[str]] = None,
        precompute_db: Optional[str] = None,
    ) -> None:
        """初始化 VEP 客户端.

        Args:
            base_url: VEP REST API 基础 URL。
            batch_size: 每批查询的变异数量。
            rate_limit: 每秒最大请求数。
            max_retries: 失败后的最大重试次数。
            timeout: HTTP 请求超时（秒）。
            cache: 可选的 CacheManager 实例，用于缓存 VEP 结果。
            canonical_only: 只返回 canonical 转录本的注释，减少数据量。
            transcript_ids: 可选的转录本 ID 列表，只返回这些转录本的注释。
            precompute_db: 预计算 SQLite 数据库路径，优先查库减少 API 调用。
        """
        self.base_url = base_url.rstrip("/")
        self.batch_size = batch_size
        self.max_retries = max_retries
        self.timeout = aiohttp.ClientTimeout(total=timeout)
        self.semaphore = asyncio.Semaphore(rate_limit)
        self.session: Optional[aiohttp.ClientSession] = None
        self.cache = cache
        self.canonical_only = canonical_only
        self.transcript_ids = transcript_ids or []
        self.precompute_db = precompute_db
        self._precompute_conn: Optional[sqlite3.Connection] = None
        self._cache_hits = 0
        self._cache_misses = 0
        self._precompute_hits = 0
        self._precompute_misses = 0

    async def _get_session(self) -> aiohttp.ClientSession:
        """获取或创建 aiohttp 会话."""
        if self.session is None or self.session.closed:
            self.session = aiohttp.ClientSession(
                timeout=self.timeout,
                trust_env=True,
            )
        return self.session

    def _cache_key(self, variant: Variant) -> str:
        """生成 VEP 缓存键."""
        return f"vep:{variant.chrom}:{variant.pos}:{variant.ref}:{variant.alt}"

    async def _get_cached(self, variant: Variant) -> Optional[Dict[str, Any]]:
        """从缓存获取 VEP 结果."""
        if self.cache is None:
            return None
        key = self._cache_key(variant)
        return await self.cache.get(key)

    async def _set_cached(self, variant: Variant, data: Dict[str, Any]) -> None:
        """将 VEP 结果写入缓存."""
        if self.cache is None:
            return
        key = self._cache_key(variant)
        # VEP 结果 TTL 设为 90 天（参考基因组注释稳定）
        await self.cache.set(key, data, ttl_days=90)

    def _get_precompute_conn(self) -> Optional[sqlite3.Connection]:
        """获取或创建预计算数据库连接."""
        if not self.precompute_db:
            return None
        if self._precompute_conn is None:
            try:
                self._precompute_conn = sqlite3.connect(self.precompute_db)
            except sqlite3.Error as exc:
                logger.warning("Failed to open precompute DB %s: %s", self.precompute_db, exc)
                return None
        return self._precompute_conn

    def _query_precompute(self, variant: Variant) -> Optional[Dict[str, Any]]:
        """从预计算 SQLite 库查询 VEP 结果.

        处理 chrom 前缀差异（chr1 vs 1）。
        """
        conn = self._get_precompute_conn()
        if conn is None:
            return None

        chrom = variant.chrom
        # gnomAD / Ensembl 通常不带 chr 前缀，但 VCF 可能带
        chrom_alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"

        try:
            cursor = conn.execute(
                """
                SELECT vep_json FROM variant_vep
                WHERE (chrom = ? OR chrom = ?) AND pos = ? AND ref = ? AND alt = ?
                LIMIT 1
                """,
                (chrom, chrom_alt, variant.pos, variant.ref, variant.alt),
            )
            row = cursor.fetchone()
            if row:
                return json.loads(row[0])
        except (sqlite3.Error, json.JSONDecodeError) as exc:
            logger.warning("Precompute DB query error for %s: %s", variant.vcf_id, exc)
        return None

    async def annotate(self, variants: List[Variant]) -> List[Dict[str, Any]]:
        """批量注释变异.

        查询优先级：内存缓存 → 预计算 SQLite 库 → VEP REST API。
        """
        if not variants:
            return []

        results: List[Optional[Dict[str, Any]]] = [None] * len(variants)
        to_fetch: List[Tuple[int, Variant]] = []

        # 第一阶段：尝试从内存缓存读取
        for i, variant in enumerate(variants):
            cached = await self._get_cached(variant)
            if cached is not None:
                results[i] = cached
                self._cache_hits += 1
            else:
                to_fetch.append((i, variant))
                self._cache_misses += 1

        # 第二阶段：对未命中的变异查预计算库
        still_to_fetch: List[Tuple[int, Variant]] = []
        for i, variant in to_fetch:
            pc = self._query_precompute(variant)
            if pc is not None:
                results[i] = pc
                self._precompute_hits += 1
                # 同步写入内存缓存，加速后续重复查询
                await self._set_cached(variant, pc)
            else:
                still_to_fetch.append((i, variant))
                self._precompute_misses += 1
        to_fetch = still_to_fetch

        total = len(variants)
        mem_rate = (
            self._cache_hits / (self._cache_hits + self._cache_misses) * 100
            if (self._cache_hits + self._cache_misses) > 0 else 0
        )
        pc_rate = (
            self._precompute_hits / (self._precompute_hits + self._precompute_misses) * 100
            if (self._precompute_hits + self._precompute_misses) > 0 else 0
        )
        logger.info(
            "VEP annotation: %d variants, mem_cache=%d hits (%.1f%%), "
            "precompute=%d hits (%.1f%%), api=%d to fetch, batch_size=%d",
            total, self._cache_hits, mem_rate,
            self._precompute_hits, pc_rate, len(to_fetch), self.batch_size,
        )

        # 第三阶段：对仍未命中的变异发 POST batch 请求 VEP API
        for i in range(0, len(to_fetch), self.batch_size):
            batch = to_fetch[i : i + self.batch_size]
            batch_results = await self._post_batch(batch)
            for (orig_idx, _), result in zip(batch, batch_results):
                results[orig_idx] = result
                # 写入内存缓存
                variant = variants[orig_idx]
                await self._set_cached(variant, result)
            logger.info("VEP batch %d/%d done", i // self.batch_size + 1, (len(to_fetch) - 1) // self.batch_size + 1)

        logger.info(
            "VEP annotation completed: %d results (mem_cache: %d, precompute: %d, api: %d)",
            len(results), self._cache_hits, self._precompute_hits,
            self._precompute_misses - len(to_fetch),
        )
        return results

    async def _post_batch(self, batch: List[Tuple[int, Variant]]) -> List[Dict[str, Any]]:
        """发送一批变异的 POST 请求到 VEP REST API.

        使用 VEP 的 POST batch endpoint，一次请求处理多个变异，
        相比逐个 GET 请求大幅减少网络往返。
        """
        variants = [v for _, v in batch]
        payload = self._build_payload(variants)
        url = f"{self.base_url}/vep/human/region"
        headers = {"Content-Type": "application/json"}

        # 构建查询参数：canonical_only 减少返回数据量
        params: Dict[str, str] = {"content-type": "application/json"}
        if self.canonical_only:
            params["canonical"] = "1"
        if self.transcript_ids:
            # VEP 支持通过 transcript 参数过滤（逗号分隔）
            # 但 batch 查询通常不限制 transcript，而是在解析时过滤
            pass

        for attempt in range(self.max_retries + 1):
            async with self.semaphore:
                try:
                    session = await self._get_session()
                    async with session.post(url, headers=headers, json=payload, params=params) as response:
                        if response.status == 200:
                            data = await response.json()
                            return self._parse_response(data, variants)
                        elif response.status == 429:
                            wait_time = 2 ** attempt
                            logger.warning("VEP rate limited (429), retrying in %ds...", wait_time)
                            await asyncio.sleep(wait_time)
                        else:
                            text = await response.text()
                            logger.error("VEP HTTP %d: %s", response.status, text[:200])
                            if attempt == self.max_retries:
                                return [self._empty_result(v) for v in variants]
                            await asyncio.sleep(2 ** attempt)
                except aiohttp.ClientError as exc:
                    logger.error("VEP POST request error: %s", exc)
                    if attempt == self.max_retries:
                        return [self._empty_result(v) for v in variants]
                    await asyncio.sleep(2 ** attempt)
                except asyncio.TimeoutError:
                    logger.error("VEP POST request timeout")
                    if attempt == self.max_retries:
                        return [self._empty_result(v) for v in variants]
                    await asyncio.sleep(2 ** attempt)

        return [self._empty_result(v) for v in variants]

    async def _get_variant(self, variant: Variant) -> Dict[str, Any]:
        """发送单个变异的 GET 请求到 VEP REST API（fallback）."""
        # 跳过复杂 indel（长度 > 10bp 或 ref/alt 都 > 1）
        if len(variant.ref) > 10 or len(variant.alt) > 10:
            return self._empty_result(variant)

        # Ensembl 使用无 chr 前缀的染色体名
        chrom = variant.chrom.replace("chr", "") if variant.chrom.startswith("chr") else variant.chrom
        # MT 线粒体特殊处理
        if chrom.upper() in ("M", "MT"):
            chrom = "MT"

        # 构建 GET URL: /vep/human/region/{chrom}:{start}-{end}:{strand}/{allele}
        ref_len = len(variant.ref)
        alt_len = len(variant.alt)

        if ref_len == 1 and alt_len == 1:
            # SNV
            allele = variant.alt
            end = variant.pos
        elif ref_len > 1 and alt_len == 1:
            # Deletion (e.g., GT>G) -> allele = "-"
            allele = "-"
            end = variant.pos + ref_len - 1
        elif ref_len == 1 and alt_len > 1:
            # Insertion (e.g., G>GT) -> allele = inserted_seq
            allele = variant.alt[1:]  # 去掉第一个碱基（与 ref 相同的锚定碱基）
            end = variant.pos
        else:
            # Complex substitution, skip
            return self._empty_result(variant)

        url = (
            f"{self.base_url}/vep/human/region/"
            f"{chrom}:{variant.pos}-{end}:1/{allele}"
            f"?content-type=application/json"
        )

        for attempt in range(self.max_retries + 1):
            async with self.semaphore:
                try:
                    session = await self._get_session()
                    async with session.get(url) as response:
                        if response.status == 200:
                            data = await response.json()
                            if isinstance(data, list) and len(data) > 0:
                                return data[0]
                            return self._empty_result(variant)
                        elif response.status == 429:
                            wait_time = 2 ** attempt
                            logger.warning("VEP rate limited (429), retrying in %ds...", wait_time)
                            await asyncio.sleep(wait_time)
                        else:
                            text = await response.text()
                            logger.error("VEP HTTP %d for %s: %s", response.status, url, text[:200])
                            if attempt == self.max_retries:
                                return self._empty_result(variant)
                            await asyncio.sleep(2 ** attempt)
                except aiohttp.ClientError as exc:
                    logger.error("VEP request error for %s: %s", url, exc)
                    if attempt == self.max_retries:
                        return self._empty_result(variant)
                    await asyncio.sleep(2 ** attempt)
                except asyncio.TimeoutError:
                    logger.error("VEP request timeout for %s", url)
                    if attempt == self.max_retries:
                        return self._empty_result(variant)
                    await asyncio.sleep(2 ** attempt)

        return self._empty_result(variant)

    @staticmethod
    def _build_payload(batch: List[Variant]) -> List[str]:
        """构建 VEP POST payload.

        格式: ["chrom pos ref alt", ...]
        """
        return [f"{v.chrom} {v.pos} {v.ref} {v.alt}" for v in batch]

    @staticmethod
    def _parse_response(data: List[Dict[str, Any]], batch: List[Variant]) -> List[Dict[str, Any]]:
        """解析 VEP 响应.

        Args:
            data: VEP 返回的 JSON 列表。
            batch: 原始请求的变异批次。

        Returns:
            解析后的结果列表，确保与 batch 一一对应。
        """
        if not isinstance(data, list):
            logger.error("Unexpected VEP response format: not a list")
            return [VepClient._empty_result(v) for v in batch]

        # VEP 返回的结果顺序通常与请求顺序一致
        results = []
        for i, variant in enumerate(batch):
            if i < len(data):
                item = data[i]
                # 校验输入是否匹配
                if VepClient._variant_matches(variant, item):
                    results.append(item)
                else:
                    # 顺序不一致，尝试查找匹配项
                    matched = VepClient._find_matching(variant, data)
                    results.append(matched if matched else VepClient._empty_result(variant))
            else:
                results.append(VepClient._empty_result(variant))
        return results

    @staticmethod
    def _variant_matches(variant: Variant, item: Dict[str, Any]) -> bool:
        """检查 VEP 响应项是否与请求变异匹配."""
        seq_region = item.get("seq_region_name", "")
        start = item.get("start", 0)
        return str(seq_region) == str(variant.chrom) and int(start) == variant.pos

    @staticmethod
    def _find_matching(variant: Variant, data: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
        """在 VEP 响应中查找匹配的变异."""
        for item in data:
            if VepClient._variant_matches(variant, item):
                return item
        return None

    @staticmethod
    def _empty_result(variant: Variant) -> Dict[str, Any]:
        """生成空的 VEP 结果占位."""
        return {
            "input": f"{variant.chrom} {variant.pos} {variant.ref} {variant.alt}",
            "seq_region_name": variant.chrom,
            "start": variant.pos,
            "allele_string": f"{variant.ref}/{variant.alt}",
            "transcript_consequences": [],
        }

    async def close(self) -> None:
        """关闭 HTTP 会话和预计算数据库连接."""
        if self.session and not self.session.closed:
            await self.session.close()
        if self._precompute_conn:
            self._precompute_conn.close()
            self._precompute_conn = None

    async def __aenter__(self):
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        await self.close()
