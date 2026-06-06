"""测试 SQLite 异步缓存管理器.

覆盖 CacheManager 的 get/set/TTL/过期清理/清空，及并发安全。
"""

import asyncio
from pathlib import Path

import pytest

from src.enrichment.cache import CacheManager


@pytest.fixture
def cache(tmp_path: Path) -> CacheManager:
    db_path = tmp_path / "test_cache.sqlite"
    return CacheManager(db_path=str(db_path), default_ttl_days=1)


@pytest.mark.asyncio
class TestCacheManager:
    """CacheManager 异步测试."""

    async def test_set_and_get(self, cache: CacheManager) -> None:
        await cache.set("key1", {"data": "value1"})
        result = await cache.get("key1")
        assert result == {"data": "value1"}

    async def test_get_missing(self, cache: CacheManager) -> None:
        result = await cache.get("nonexistent")
        assert result is None

    async def test_update_existing(self, cache: CacheManager) -> None:
        await cache.set("key1", {"data": "old"})
        await cache.set("key1", {"data": "new"})
        result = await cache.get("key1")
        assert result == {"data": "new"}

    async def test_custom_ttl(self, cache: CacheManager) -> None:
        # TTL = 0 天意味着立即过期（取决于时间精度，可能刚好过期）
        await cache.set("key1", {"data": "value"}, ttl_days=0)
        # 由于 ttl_days=0  expires_at == created_at, 立即过期
        result = await cache.get("key1")
        assert result is None

    async def test_purge_expired(self, cache: CacheManager) -> None:
        await cache.set("key1", {"data": "value"}, ttl_days=-1)
        count = await cache.purge_expired()
        assert count >= 1
        result = await cache.get("key1")
        assert result is None

    async def test_clear(self, cache: CacheManager) -> None:
        await cache.set("key1", {"data": "value1"})
        await cache.set("key2", {"data": "value2"})
        await cache.clear()
        assert await cache.get("key1") is None
        assert await cache.get("key2") is None

    async def test_non_serializable_value(self, cache: CacheManager) -> None:
        """测试非 JSON 原生类型的序列化（如 datetime）."""
        from datetime import datetime

        await cache.set("key1", {"time": datetime(2024, 1, 1, 12, 0, 0)})
        result = await cache.get("key1")
        assert result is not None
        assert "2024-01-01" in result["time"]

    async def test_unicode_value(self, cache: CacheManager) -> None:
        await cache.set("key1", {"text": "中文测试"})
        result = await cache.get("key1")
        assert result == {"text": "中文测试"}

    async def test_concurrent_access(self, cache: CacheManager) -> None:
        """并发读写不应抛出异常."""

        async def writer(n: int) -> None:
            for i in range(10):
                await cache.set(f"key_{n}_{i}", {"n": n, "i": i})

        async def reader(n: int) -> None:
            for i in range(10):
                await cache.get(f"key_{n}_{i}")

        tasks = [writer(i) for i in range(5)] + [reader(i) for i in range(5)]
        await asyncio.gather(*tasks)
        # 只要没抛异常就算通过
        assert True
