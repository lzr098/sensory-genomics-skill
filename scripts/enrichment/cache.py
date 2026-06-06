"""SQLite 异步缓存管理器.

支持 TTL 过期、异步读写、定期清理过期条目。
"""

import json
import os
import sqlite3
from datetime import datetime, timedelta, timezone
from typing import Any, Dict, Optional

import aiosqlite

from src.logger import get_logger

logger = get_logger(__name__)


class CacheManager:
    """异步 SQLite 缓存管理器.

    缓存键格式: {api_name}:{query_key}
    默认 TTL: 30 天
    """

    def __init__(self, db_path: str, default_ttl_days: int = 30) -> None:
        """初始化缓存管理器.

        Args:
            db_path: SQLite 数据库文件路径。
            default_ttl_days: 默认缓存有效期（天）。
        """
        self.db_path = os.path.expanduser(db_path)
        self.default_ttl_days = default_ttl_days
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)
        self._initialized = False

    async def _ensure_table(self) -> None:
        """确保缓存表存在（幂等）."""
        if self._initialized:
            return
        async with aiosqlite.connect(self.db_path) as db:
            await db.execute(
                """
                CREATE TABLE IF NOT EXISTS cache (
                    key TEXT PRIMARY KEY,
                    value TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    expires_at TEXT NOT NULL
                )
                """
            )
            await db.execute(
                "CREATE INDEX IF NOT EXISTS idx_expires ON cache(expires_at)"
            )
            await db.commit()
        self._initialized = True

    async def get(self, key: str) -> Optional[Dict[str, Any]]:
        """从缓存获取数据.

        Args:
            key: 缓存键。

        Returns:
            缓存的字典数据，若不存在或已过期则返回 None。
        """
        await self._ensure_table()
        now = datetime.now(timezone.utc).isoformat()
        async with aiosqlite.connect(self.db_path) as db:
            async with db.execute(
                "SELECT value FROM cache WHERE key = ? AND expires_at > ?",
                (key, now),
            ) as cursor:
                row = await cursor.fetchone()
                if row:
                    try:
                        return json.loads(row[0])
                    except json.JSONDecodeError:
                        logger.warning("Cache corruption for key %s", key)
                        return None
        return None

    async def set(self, key: str, value: Dict[str, Any], ttl_days: Optional[int] = None) -> None:
        """写入缓存.

        Args:
            key: 缓存键。
            value: 要缓存的字典数据。
            ttl_days: 自定义有效期（天），默认使用初始化时的默认值。
        """
        await self._ensure_table()
        ttl = ttl_days if ttl_days is not None else self.default_ttl_days
        now = datetime.now(timezone.utc)
        expires = now + timedelta(days=ttl)
        value_str = json.dumps(value, ensure_ascii=False, default=str)
        async with aiosqlite.connect(self.db_path) as db:
            await db.execute(
                """
                INSERT INTO cache (key, value, created_at, expires_at)
                VALUES (?, ?, ?, ?)
                ON CONFLICT(key) DO UPDATE SET
                    value = excluded.value,
                    created_at = excluded.created_at,
                    expires_at = excluded.expires_at
                """,
                (key, value_str, now.isoformat(), expires.isoformat()),
            )
            await db.commit()

    async def purge_expired(self) -> int:
        """清理过期缓存条目.

        Returns:
            清理的条目数量。
        """
        await self._ensure_table()
        now = datetime.now(timezone.utc).isoformat()
        async with aiosqlite.connect(self.db_path) as db:
            async with db.execute(
                "DELETE FROM cache WHERE expires_at <= ?", (now,)
            ) as cursor:
                await db.commit()
                return cursor.rowcount if cursor.rowcount is not None else 0

    async def clear(self) -> None:
        """清空全部缓存."""
        await self._ensure_table()
        async with aiosqlite.connect(self.db_path) as db:
            await db.execute("DELETE FROM cache")
            await db.commit()
