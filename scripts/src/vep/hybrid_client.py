"""Hybrid VEP 客户端：本地 Docker 优先，自动 fallback 到 REST API.

自动检测本地 Docker VEP 环境，可用时走离线注释，不可用或无 cache 时
无缝 fallback 到 VEP REST API。对外接口与 VepClient / LocalVepClient 完全兼容。
"""

import asyncio
import os
import subprocess
from typing import Any, Dict, List, Optional

from src.logger import get_logger
from src.models import Variant
from src.vep.client import VepClient
from src.vep.local_client import LocalVepClient

logger = get_logger(__name__)


class HybridVepClient:
    """混合 VEP 客户端：本地 Docker 优先，REST API fallback.

    使用策略：
    1. 初始化时检测 Docker + VEP cache 是否就绪
    2. annotate() 优先调用 LocalVepClient（离线、无 rate limit）
    3. 本地失败（Docker 未启动、cache 缺失等）→ 自动 fallback 到 VepClient
    4. 用户可通过 config 强制指定来源
    """

    def __init__(
        self,
        *,
        local_config: Optional[Dict[str, Any]] = None,
        rest_config: Optional[Dict[str, Any]] = None,
        force_source: Optional[str] = None,
    ) -> None:
        """初始化混合 VEP 客户端.

        Args:
            local_config: LocalVepClient 配置 dict，如 {"cache_dir": "..."}。
            rest_config: VepClient 配置 dict，如 {"base_url": "...", "cache": ...}。
            force_source: 强制指定来源。"local" | "rest" | None（自动检测）。
        """
        self.force_source = force_source
        self.local_config = local_config or {}
        self.rest_config = rest_config or {}

        self._local: Optional[LocalVepClient] = None
        self._rest: Optional[VepClient] = None
        self._source: str = "unknown"

    # ------------------------------------------------------------------
    # 环境检测
    # ------------------------------------------------------------------

    @staticmethod
    def _docker_available() -> bool:
        """检测 Docker 是否可用."""
        try:
            result = subprocess.run(
                ["docker", "info"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False

    @staticmethod
    def _docker_image_present(image: str = "ensemblorg/ensembl-vep:latest") -> bool:
        """检测 VEP Docker 镜像是否已拉取."""
        try:
            result = subprocess.run(
                ["docker", "image", "inspect", image],
                capture_output=True,
                text=True,
                timeout=10,
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False

    @staticmethod
    def _cache_ready(cache_dir: str) -> bool:
        """检测 VEP cache 目录是否包含有效数据."""
        expanded = os.path.expanduser(cache_dir)
        # 检查是否存在物种目录和 FASTA
        species_dir = os.path.join(expanded, "homo_sapiens")
        if not os.path.isdir(species_dir):
            return False
        # 检查是否有 chromosome 数据文件（至少 chr1）
        for root, _, files in os.walk(species_dir):
            for f in files:
                if f.endswith(".gz"):
                    return True
        return False

    def _detect_source(self) -> str:
        """检测最佳 VEP 来源.

        Returns:
            "local" 或 "rest"。
        """
        if self.force_source:
            logger.info("HybridVepClient: force_source=%s", self.force_source)
            return self.force_source

        cache_dir = self.local_config.get("cache_dir", "~/.workbuddy/tools/vep/cache")

        if not self._docker_available():
            logger.info("HybridVepClient: Docker unavailable → using REST API")
            return "rest"
        if not self._docker_image_present():
            logger.info("HybridVepClient: VEP Docker image missing → using REST API")
            return "rest"
        if not self._cache_ready(cache_dir):
            logger.info(
                "HybridVepClient: cache not ready at %s → using REST API",
                cache_dir,
            )
            return "rest"

        logger.info("HybridVepClient: Docker + image + cache ready → using local")
        return "local"

    # ------------------------------------------------------------------
    # Lazy init
    # ------------------------------------------------------------------

    def _ensure_local(self) -> LocalVepClient:
        if self._local is None:
            self._local = LocalVepClient(**self.local_config)
        return self._local

    def _ensure_rest(self) -> VepClient:
        if self._rest is None:
            self._rest = VepClient(**self.rest_config)
        return self._rest

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    async def annotate(self, variants: List[Variant]) -> List[Dict[str, Any]]:
        """注释变异列表：优先本地，fallback REST.

        Args:
            variants: 待注释的变异列表。

        Returns:
            与 VepClient.annotate() 格式兼容的注释结果列表。
        """
        if not variants:
            return []

        source = self._detect_source()
        self._source = source

        if source == "local":
            try:
                local = self._ensure_local()
                results = await local.annotate(variants)
                logger.info(
                    "HybridVepClient: annotated %d variants via LOCAL",
                    len(results),
                )
                return results
            except Exception as exc:
                logger.warning(
                    "HybridVepClient: local VEP failed (%s), falling back to REST",
                    exc,
                )
                # fallback 到 REST
                rest = self._ensure_rest()
                results = await rest.annotate(variants)
                logger.info(
                    "HybridVepClient: annotated %d variants via REST (fallback)",
                    len(results),
                )
                return results
        else:
            rest = self._ensure_rest()
            results = await rest.annotate(variants)
            logger.info(
                "HybridVepClient: annotated %d variants via REST",
                len(results),
            )
            return results

    @property
    def active_source(self) -> str:
        """返回当前实际使用的来源（local / rest / unknown）."""
        return self._source

    async def close(self) -> None:
        """关闭所有底层客户端."""
        if self._local is not None:
            await self._local.close()
            self._local = None
        if self._rest is not None:
            await self._rest.close()
            self._rest = None
