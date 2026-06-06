"""Local VEP client using Docker (offline).

Replaces VepClient REST API calls with local Ensembl VEP via Docker.
Requires:
    - Docker Desktop running
    - ensemblorg/ensembl-vep image pulled
    - GRCh38 cache downloaded to ~/.workbuddy/tools/vep/cache/

Usage:
    client = LocalVepClient(cache_dir="~/.workbuddy/tools/vep/cache")
    results = await client.annotate(variants)
"""

import asyncio
import json
import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional

from src.models import Variant

logger = logging.getLogger(__name__)


class LocalVepClient:
    """Offline VEP annotation using local Docker VEP.

    Translates variants to VCF format, mounts cache volume, runs VEP,
    parses JSON output back into the same format as VepClient.
    """

    def __init__(
        self,
        cache_dir: str = "~/.workbuddy/tools/vep/cache",
        docker_image: str = "ensemblorg/ensembl-vep",
        assembly: str = "GRCh38",
        fork: int = 4,
    ) -> None:
        self.cache_dir = Path(cache_dir).expanduser()
        self.docker_image = docker_image
        self.assembly = assembly
        self.fork = fork
        self._checked = False

    def _check_prerequisites(self) -> None:
        """Verify Docker and cache are available."""
        if self._checked:
            return

        # Check docker
        result = subprocess.run(
            ["docker", "info"],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError("Docker is not running. Please start Docker Desktop.")

        # Check image
        result = subprocess.run(
            ["docker", "images", "-q", self.docker_image],
            capture_output=True,
            text=True,
        )
        if not result.stdout.strip():
            raise RuntimeError(
                f"Docker image {self.docker_image} not found. "
                "Run: docker pull ensemblorg/ensembl-vep"
            )

        # Check cache
        cache_path = self.cache_dir / "homo_sapiens"
        if not cache_path.exists():
            logger.warning(
                "VEP cache not found at %s. "
                "Run: docker run -v %s:/data %s INSTALL.pl -a cf -s homo_sapiens -y GRCh38 -d /data",
                cache_path,
                self.cache_dir,
                self.docker_image,
            )
            raise RuntimeError(f"VEP cache not found at {cache_path}")

        self._checked = True
        logger.info("Local VEP ready: cache=%s, image=%s", self.cache_dir, self.docker_image)

    def _variants_to_vcf(self, variants: List[Variant]) -> str:
        """Convert variants to VCF-like string for VEP input.

        VEP region format: chrom start end ref/alt strand
        """
        lines = []
        for v in variants:
            chrom = v.chrom.replace("chr", "")
            # VEP region input: chrom  start  end  allele/strand
            # For SNV: start=end=pos
            # Format: "1  230710048  230710048  A/G  +"
            lines.append(f"{chrom}\t{v.pos}\t{v.pos}\t{v.ref}/{v.alt}\t+")
        return "\n".join(lines) + "\n"

    def _parse_vep_json(self, raw: str) -> List[Dict[str, Any]]:
        """Parse VEP JSON output into list of response dicts.

        Matches the format returned by Ensembl VEP REST API.
        """
        results: List[Dict[str, Any]] = []
        for line in raw.strip().split("\n"):
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
                results.append(record)
            except json.JSONDecodeError:
                logger.warning("Failed to parse VEP JSON line: %s", line[:200])
        return results

    async def annotate(self, variants: List[Variant]) -> List[Dict[str, Any]]:
        """Annotate variants using local Docker VEP.

        Returns list of VEP JSON records, one per input variant.
        """
        if not variants:
            return []

        self._check_prerequisites()

        # Write input to temp file
        vcf_content = self._variants_to_vcf(variants)

        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(vcf_content)
            input_path = Path(f.name)

        output_path = input_path.with_suffix(".json")

        try:
            # Build docker command
            cmd = [
                "docker", "run", "--rm",
                "-v", f"{self.cache_dir}:/data/vep_cache:ro",
                "-v", f"{input_path.parent}:/data/input:ro",
                self.docker_image,
                "vep",
                "--cache", "--offline",
                "--dir_cache", "/data/vep_cache",
                "--assembly", self.assembly,
                "--input_file", f"/data/input/{input_path.name}",
                "--output_file", f"/data/input/{output_path.name}",
                "--format", "region",
                "--json",
                "--canonical",
                "--hgvs",
                "--symbol",
                "--protein",
                "--sift", "b",
                "--polyphen", "b",
                "--variant_class",
                "--fork", str(self.fork),
            ]

            logger.info("Running local VEP for %d variants...", len(variants))

            # Run in executor since subprocess is blocking
            loop = asyncio.get_event_loop()
            proc = await loop.run_in_executor(
                None,
                lambda: subprocess.run(cmd, capture_output=True, text=True),
            )

            if proc.returncode != 0:
                logger.error("VEP failed: %s", proc.stderr)
                raise RuntimeError(f"VEP annotation failed: {proc.stderr[:500]}")

            # Read output
            if not output_path.exists():
                logger.error("VEP output file not created")
                raise RuntimeError("VEP did not produce output file")

            raw_output = output_path.read_text()
            results = self._parse_vep_json(raw_output)

            logger.info("Local VEP completed: %d variants → %d results", len(variants), len(results))
            return results

        finally:
            # Cleanup temp files
            input_path.unlink(missing_ok=True)
            output_path.unlink(missing_ok=True)

    async def close(self) -> None:
        """No-op for local client (no persistent connections)."""
        pass
