"""本地 Docker VEP 客户端.

通过 Docker 调用 ensemblorg/ensembl-vep 镜像进行离线注释。
输出格式与 REST API VepClient 兼容（JSON Lines）。
"""

import asyncio
import json
import os
import subprocess
import tempfile
from typing import Any, Dict, List, Optional

from src.logger import get_logger
from src.models import Variant

logger = get_logger(__name__)


class LocalVepClient:
    """Docker 离线 VEP 客户端.

    将 Variant 列表写入临时 VCF，通过 Docker 挂载本地 cache 运行 VEP，
    解析 JSON Lines 输出为与 REST API 兼容的 dict 列表。
    """

    DOCKER_IMAGE = "ensemblorg/ensembl-vep:latest"
    CONTAINER_CACHE_DIR = "/data/vep_cache"

    def __init__(
        self,
        cache_dir: str = "~/.workbuddy/tools/vep/cache",
        species: str = "homo_sapiens",
        assembly: str = "GRCh38",
        cache_version: int = 115,
    ) -> None:
        """初始化本地 VEP 客户端.

        Args:
            cache_dir: 本地 VEP cache 目录（包含 homo_sapiens/115_GRCh38/）。
            species: 物种名称。
            assembly: 基因组组装版本。
            cache_version: VEP cache 版本号。
        """
        self.cache_dir = os.path.expanduser(cache_dir)
        self.species = species
        self.assembly = assembly
        self.cache_version = cache_version
        self._tmp_files: List[str] = []

    async def annotate(self, variants: List[Variant]) -> List[Dict[str, Any]]:
        """注释变异列表（异步接口，与 VepClient 兼容）.

        Args:
            variants: 待注释的变异列表。

        Returns:
            与 VepClient.annotate() 格式兼容的注释结果列表。
            每个元素对应一个变异的 VEP JSON dict。
        """
        if not variants:
            return []

        # Docker VEP 是 CPU/IO 密集型同步操作，放到线程池避免阻塞事件循环
        loop = asyncio.get_event_loop()
        return await loop.run_in_executor(None, self._annotate_sync, variants)

    def _annotate_sync(self, variants: List[Variant]) -> List[Dict[str, Any]]:
        """同步执行 Docker VEP 注释."""
        # 1. 写入临时 VCF
        vcf_path = self._write_vcf(variants)
        self._tmp_files.append(vcf_path)

        # 2. 创建临时输出目录
        out_dir = tempfile.mkdtemp(prefix="vep_local_")
        self._tmp_files.append(out_dir)

        output_file = os.path.join(out_dir, "vep_output.json")

        # 3. 构造 Docker 命令
        cmd = [
            "docker", "run", "--rm",
            "-v", f"{self.cache_dir}:{self.CONTAINER_CACHE_DIR}:ro",
            "-v", f"{os.path.dirname(vcf_path)}:/data/input:ro",
            "-v", f"{out_dir}:/data/output",
            self.DOCKER_IMAGE,
            "vep",
            "--cache", "--offline",
            "--dir_cache", self.CONTAINER_CACHE_DIR,
            "--input_file", f"/data/input/{os.path.basename(vcf_path)}",
            "--output_file", "/data/output/vep_output.json",
            "--json",
            "--force_overwrite",
            "--species", self.species,
            "--assembly", self.assembly,
            "--cache_version", str(self.cache_version),
            "--symbol",
            "--canonical",
            "--biotype",
            "--numbers",
            "--everything",
        ]

        logger.info("Running local VEP for %d variants via Docker", len(variants))

        # 4. 执行 Docker
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600,
            )
            if result.returncode != 0:
                logger.error("VEP Docker failed: %s", result.stderr)
                raise RuntimeError(f"VEP Docker failed: {result.stderr}")
        except subprocess.TimeoutExpired:
            logger.error("VEP Docker timed out after 600s")
            raise RuntimeError("VEP Docker timed out")

        # 5. 解析 JSON Lines 输出
        if not os.path.exists(output_file):
            logger.error("VEP output file not found: %s", output_file)
            raise RuntimeError("VEP output file not found")

        results = self._parse_json_lines(output_file)
        logger.info("Local VEP annotated %d variants", len(results))
        return results

    def _write_vcf(self, variants: List[Variant]) -> str:
        """将 Variant 列表写入临时 VCF 文件.

        Args:
            variants: 变异列表。

        Returns:
            临时 VCF 文件路径。
        """
        fd, path = tempfile.mkstemp(suffix=".vcf")
        try:
            with os.fdopen(fd, "w") as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                for v in variants:
                    f.write(
                        f"{v.chrom}\t{v.pos}\t.\t{v.ref}\t{v.alt}\t.\tPASS\t.\n"
                    )
        except Exception:
            os.close(fd)
            raise
        return path

    @staticmethod
    def _parse_json_lines(path: str) -> List[Dict[str, Any]]:
        """解析 VEP JSON Lines 输出.

        VEP --json 输出每行一个 JSON 对象。

        Args:
            path: JSON Lines 文件路径。

        Returns:
            JSON dict 列表。
        """
        results = []
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    data = json.loads(line)
                    results.append(data)
                except json.JSONDecodeError as exc:
                    logger.warning("Failed to parse VEP JSON line: %s", exc)
        return results

    async def close(self) -> None:
        """清理临时文件（异步接口，与 VepClient 兼容）."""
        for path in self._tmp_files:
            try:
                if os.path.isfile(path):
                    os.unlink(path)
                elif os.path.isdir(path):
                    import shutil
                    shutil.rmtree(path)
            except OSError as exc:
                logger.warning("Failed to clean up temp file %s: %s", path, exc)
        self._tmp_files.clear()
