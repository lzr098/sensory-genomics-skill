"""AsyncApiClient 抽象基类.

统一外部 API 的重试、限速、缓存逻辑。
"""

import asyncio
from abc import ABC, abstractmethod
from typing import Any, Dict, Optional

import aiohttp

from src.enrichment.cache import CacheManager
from src.exceptions import ApiError
from src.logger import get_logger

logger = get_logger(__name__)


class AsyncApiClient(ABC):
    """异步 API 客户端抽象基类.

    内置功能：
    - SQLite 缓存查询（Cache-Aside 模式）
    - 指数退避重试
    - 基于 Semaphore 的速率限制
    - 失败时返回带 error 标记的降级结果
    """

    def __init__(
        self,
        api_name: str,
        base_url: str,
        cache: CacheManager,
        rate_limit: int = 10,
        max_retries: int = 3,
        timeout: int = 30,
    ) -> None:
        """初始化 API 客户端.

        Args:
            api_name: API 名称（用于缓存键前缀和日志）。
            base_url: API 基础 URL。
            cache: 缓存管理器实例。
            rate_limit: 每秒最大请求数。
            max_retries: 最大重试次数。
            timeout: HTTP 超时（秒）。
        """
        self.api_name = api_name
        self.base_url = base_url.rstrip("/")
        self.cache = cache
        self.max_retries = max_retries
        self.timeout = aiohttp.ClientTimeout(total=timeout)
        self.semaphore = asyncio.Semaphore(rate_limit)
        self.session: Optional[aiohttp.ClientSession] = None

    async def _get_session(self) -> aiohttp.ClientSession:
        """获取或创建 aiohttp 会话."""
        if self.session is None or self.session.closed:
            self.session = aiohttp.ClientSession(timeout=self.timeout, trust_env=True)
        return self.session

    async def query(self, key: str) -> Dict[str, Any]:
        """查询 API，优先读取缓存.

        Args:
            key: 查询键。

        Returns:
            查询结果字典，失败时返回 {"error": "..."} 降级结果。
        """
        cache_key = f"{self.api_name}:{key}"

        # 1. 检查缓存
        cached = await self.cache.get(cache_key)
        if cached is not None:
            logger.debug("Cache hit for %s", cache_key)
            return cached

        # 2. 限速 + 查询
        async with self.semaphore:
            for attempt in range(self.max_retries + 1):
                try:
                    result = await self._fetch(key)
                    await self.cache.set(cache_key, result)
                    return result
                except aiohttp.ClientResponseError as exc:
                    if exc.status == 429:
                        wait_time = 2 ** attempt
                        logger.warning(
                            "%s rate limited (429), retrying in %ds...",
                            self.api_name,
                            wait_time,
                        )
                        await asyncio.sleep(wait_time)
                    else:
                        logger.error("%s HTTP %d: %s", self.api_name, exc.status, exc.message)
                        if attempt == self.max_retries:
                            return {
                                "error": f"HTTP {exc.status}",
                                "api": self.api_name,
                                "key": key,
                            }
                        await asyncio.sleep(2 ** attempt)
                except aiohttp.ClientError as exc:
                    logger.error("%s request error: %s", self.api_name, exc)
                    if attempt == self.max_retries:
                        return {
                            "error": f"Request failed: {exc}",
                            "api": self.api_name,
                            "key": key,
                        }
                    await asyncio.sleep(2 ** attempt)
                except asyncio.TimeoutError:
                    logger.error("%s request timeout", self.api_name)
                    if attempt == self.max_retries:
                        return {
                            "error": "Timeout",
                            "api": self.api_name,
                            "key": key,
                        }
                    await asyncio.sleep(2 ** attempt)

        return {
            "error": "Max retries exceeded",
            "api": self.api_name,
            "key": key,
        }

    @abstractmethod
    async def _fetch(self, key: str) -> Dict[str, Any]:
        """实际执行 API 查询.

        Args:
            key: 查询键。

        Returns:
            API 响应字典。

        Raises:
            aiohttp.ClientError: HTTP 请求失败。
        """
        raise NotImplementedError

    async def close(self) -> None:
        """关闭 HTTP 会话."""
        if self.session and not self.session.closed:
            await self.session.close()

    async def __aenter__(self):
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        await self.close()
