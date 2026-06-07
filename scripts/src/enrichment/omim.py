"""OMIM gene-disease enrichment for sensory-genomics.

Wraps the shared ~/.workbuddy/scripts/omim_local.py for:
- Gene-disease association lookup
- Inheritance pattern retrieval
- Sensory system phenotype search

Uses the local OMIM SQLite database (no network required).
"""

from pathlib import Path
from typing import Any

_SHARED_OMIM = Path.home() / ".workbuddy/scripts/omim_local.py"

_omim = None
_available = False


def _load():
    global _omim, _available
    if _omim is not None:
        return _omim
    if not _SHARED_OMIM.exists():
        _available = False
        return None
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location("omim_local", _SHARED_OMIM)
        if spec is None or spec.loader is None:
            _available = False
            return None
        _omim = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_omim)
        _available = True
        return _omim
    except Exception:
        _available = False
        return None


def is_available() -> bool:
    """Check if OMIM database is available."""
    _load()
    return _available


def get_gene_disease_info(gene: str) -> dict[str, Any] | None:
    """Get OMIM gene-disease annotation for a single gene."""
    omim = _load()
    if omim is None:
        return None
    try:
        return omim.get_gene_phenotype(gene)
    except Exception:
        return None


def is_mendelian_gene(gene: str) -> bool:
    """Check if gene has Mendelian disease association."""
    omim = _load()
    if omim is None:
        return False
    try:
        return omim.is_mendelian_gene(gene)
    except Exception:
        return False


def search_sensory_phenotypes(keywords: list[str], limit: int = 50) -> list[dict[str, Any]]:
    """Search OMIM for sensory system-related phenotypes.

    Args:
        keywords: List of sensory-related keywords (e.g. ['retinal', 'hearing', 'olfactory'])
        limit: Max results per keyword

    Returns:
        Combined list of matching OMIM entries
    """
    omim = _load()
    if omim is None:
        return []

    all_results = []
    seen = set()
    for kw in keywords:
        try:
            results = omim.search_disease(kw, limit=limit)
            for r in results:
                mim = r.get("mim_number")
                if mim and mim not in seen:
                    seen.add(mim)
                    all_results.append(r)
        except Exception:
            continue

    return all_results
