#!/usr/bin/env python3
"""
GRCh38 FASTA Enrichment for Sensory Genomics
v0.1.0 - 2026-06-07

Provides genomic sequence context for sensory system variant reports.
No network required (offline capable).

Depends on shared module: ~/.workbuddy/scripts/grch38_fasta_local.py
"""

import sys
from pathlib import Path
from typing import Optional

_FASTA_MODULE = Path("/Users/zhaorongli/.workbuddy/scripts/grch38_fasta_local.py")


def _import_fasta():
    if not _FASTA_MODULE.exists():
        return None
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location("grch38_fasta_local", _FASTA_MODULE)
        if spec is None or spec.loader is None:
            return None
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod
    except Exception:
        return None


def get_variant_context(chrom: str, pos: int, ref: str, alt: str, flank: int = 30) -> Optional[dict]:
    """Get genomic context around a sensory system variant.

    Args:
        chrom: Chromosome (e.g. "1" or "chr1")
        pos: 1-based position
        ref: Reference allele
        alt: Alternate allele
        flank: Bases to include on each side

    Returns:
        Dict with upstream, ref, alt, downstream, full_context, or None.
    """
    fasta = _import_fasta()
    if fasta is None:
        return None
    try:
        return fasta.get_flanking_sequence(chrom, pos, ref, alt, flank=flank)
    except Exception:
        return None


def verify_ref(chrom: str, pos: int, ref: str) -> tuple[bool, Optional[str]]:
    """Verify if REF allele matches GRCh38 reference genome.

    Returns:
        (is_correct, actual_ref_or_error)
    """
    fasta = _import_fasta()
    if fasta is None:
        return True, None  # FASTA unavailable, assume correct
    try:
        ok, actual = fasta.verify_ref(chrom, pos, ref)
        return ok, actual if not ok else None
    except Exception as e:
        return True, str(e)


def annotate_variant(variant: dict, flank: int = 30) -> None:
    """In-place add genomic context to a variant dict.

    Args:
        variant: Dict with keys chrom, pos, ref, alt
        flank: Flanking bases
    """
    if not all(k in variant for k in ("chrom", "pos", "ref", "alt")):
        return

    ctx = get_variant_context(variant["chrom"], variant["pos"], variant["ref"], variant["alt"], flank=flank)
    if ctx:
        variant["genomic_context"] = ctx

    ok, actual = verify_ref(variant["chrom"], variant["pos"], variant["ref"])
    if not ok:
        variant["ref_warning"] = f"REF mismatch: VCF='{variant['ref']}', genome='{actual}'"
