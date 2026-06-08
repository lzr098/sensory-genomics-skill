"""OR 基因分级展示逻辑（v0.2.0 — enrichment-aware）.

Tier A: Known-ligand OR w/ HOM LoF/missense → main report §2.1 detail
Tier B: Known-ligand OR w/ HET protein-affecting → main report §2.1 compact table
Tier C: Known-ligand OR w/ syn/UTR only → appendix one-liner

Unknown-ligand OR (gnomAD AF gate):
    Layer 1 (hide): ALL variants AF > 30% → hidden
    Layer 2 (appendix): at least 1 variant AF 5-30%
    Layer 3 (main report): at least 1 variant AF < 5% OR ClinVar P/LP
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from src.assessment.engine import ImpactEngine
from src.gene_sets.loader import GeneSetLoader
from src.logger import get_logger
from src.models import ImpactAssessment, ORTierResult, Sex, Variant

logger = get_logger(__name__)

# Protein-affecting consequences
PROTEIN_CSQ = {
    "frameshift_variant", "stop_gained",
    "splice_acceptor_variant", "splice_donor_variant",
    "missense_variant", "inframe_deletion", "inframe_insertion",
}
SEVERE_CSQ = {
    "frameshift_variant", "stop_gained",
    "splice_acceptor_variant", "splice_donor_variant",
}


class ORTierClassifier:
    """OR 基因分级分类器（v0.2.0 — enrichment-aware）."""

    def __init__(self, data_path: Optional[str] = None) -> None:
        if data_path is None:
            src_dir = Path(__file__).resolve().parent.parent
            data_path = src_dir.parent / "data" / "or_ligands.json"
        else:
            data_path = Path(data_path)

        self.ligand_map: Dict[str, Dict[str, str]] = {}
        self._load_data(data_path)

    def _load_data(self, path: Path) -> None:
        if not path.exists():
            logger.warning("or_ligands.json not found, using empty defaults")
            return
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        for entry in data.get("ligands", []):
            gene = entry.get("gene_symbol", "")
            if gene:
                self.ligand_map[gene] = entry
        logger.info("Loaded OR ligand map: %d entries", len(self.ligand_map))

    # ── Public API: basic classification (no enrichment) ──

    def classify(
        self, or_variants: List[Variant], sex: Sex, gene_sets: GeneSetLoader
    ) -> List[ORTierResult]:
        """Basic OR classification (pre-enrichment, backward compatible)."""
        return self._classify_core(or_variants, sex, gene_sets, {})

    # ── Public API: enrichment-aware classification ──

    def classify_with_enrichment(
        self,
        or_variants: List[Variant],
        sex: Sex,
        gene_sets: GeneSetLoader,
        enrichment_data: Dict[str, Dict[str, Any]],
    ) -> List[ORTierResult]:
        """OR classification using gnomAD AF + ClinVar enrichment data.

        Args:
            or_variants: OR gene variants.
            sex: Sample sex.
            gene_sets: Gene set loader.
            enrichment_data: {gene: {uniprot/gnomad/clinvar/gtex/gnomad_variants}}.

        Returns:
            ORTierResult list sorted by tier (A → B → C).
        """
        return self._classify_core(or_variants, sex, gene_sets, enrichment_data)

    # ── Core classification logic ──

    def _classify_core(
        self,
        or_variants: List[Variant],
        sex: Sex,
        gene_sets: GeneSetLoader,
        enrichment: Dict[str, Dict[str, Any]],
    ) -> List[ORTierResult]:
        impact_engine = ImpactEngine()

        # Group by gene
        gene_variants: Dict[str, List[Variant]] = {}
        for variant in or_variants:
            gene = variant.gene_symbol
            if not gene or not gene.startswith("OR"):
                continue
            gene_variants.setdefault(gene, []).append(variant)

        results = []
        for gene, variants in gene_variants.items():
            ligand_info = self.ligand_map.get(gene, {})
            known_ligand = ligand_info.get("ligand_zh")
            gene_enrich = enrichment.get(gene, {})

            best_tier = "C"
            best_variant = variants[0]
            best_assessment = impact_engine.assess(best_variant, sex, gene_sets)

            for variant in variants:
                assessment = impact_engine.assess(variant, sex, gene_sets)
                tier = self._determine_tier_enriched(
                    variant, assessment, known_ligand, gene_enrich
                )
                tier_rank = {"A": 0, "B": 1, "C": 2}
                if tier_rank.get(tier, 99) < tier_rank.get(best_tier, 99):
                    best_tier = tier
                    best_variant = variant
                    best_assessment = assessment

            result = ORTierResult(
                tier=best_tier,
                gene_symbol=gene,
                known_ligand_zh=known_ligand,
                odor_description_zh=ligand_info.get("odor_description_zh"),
                variant=best_variant,
                assessment=best_assessment,
            )
            results.append(result)

        tier_order = {"A": 0, "B": 1, "C": 2}
        results.sort(key=lambda x: tier_order.get(x.tier, 99))
        return results

    # ── Tier determination (enrichment-aware) ──

    def _determine_tier_enriched(
        self,
        variant: Variant,
        assessment: ImpactAssessment,
        known_ligand: Optional[str],
        enrichment: Dict[str, Any],
    ) -> str:
        """Determine OR gene tier using enrichment data (gnomAD AF, ClinVar)."""
        consequence = variant.consequence or ""
        is_protein_affecting = consequence in PROTEIN_CSQ
        is_severe = consequence in SEVERE_CSQ

        # Get gnomAD AF for this specific variant
        gnomad_af = self._get_variant_gnomad_af(
            variant.chrom, variant.pos, variant.ref, variant.alt, enrichment
        )

        # ── Known-ligand OR → always shown, tiered by variant severity ──
        if known_ligand:
            if variant.is_homozygous and is_protein_affecting:
                return "A"  # HOM LoF/missense → detail display
            if not variant.is_homozygous and is_protein_affecting:
                # HET protein-affecting: check rarity
                if gnomad_af is not None and gnomad_af < 0.01 and is_severe:
                    return "A"  # Rare HET LoF → upgrade to detail
                return "B"  # Common HET → compact table
            if not is_protein_affecting:
                return "C"  # syn/UTR only → appendix
            return "C"

        # ── Unknown-ligand OR → gnomAD AF gate ──
        if variant.is_homozygous and is_protein_affecting:
            if gnomad_af is None:
                # Not in precompute DB → treat as unknown, keep in appendix
                return "C"
            if gnomad_af < 0.05:
                # Rare or moderately rare → promote to main
                if gnomad_af < 0.01:
                    return "B"  # Very rare, worth mentioning
                return "B"  # AF 1-5%, still worth noting
            # AF >= 5% → appendix only
            return "C"

        # Default: not shown
        return "C"

    def _get_variant_gnomad_af(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        enrichment: Dict[str, Any],
    ) -> Optional[float]:
        """Extract gnomAD AF for a specific variant from enrichment data."""
        variants = enrichment.get("gnomad_variants", [])
        for gv in variants:
            vi = gv.get("variant", {})
            if vi.get("pos") == pos and vi.get("ref") == ref and vi.get("alt") == alt:
                result = gv.get("result", {})
                if result.get("found"):
                    return result.get("gnomad_af")
        return None

    # ── Classification utilities (used by markdown generator for grouped display) ──

    @staticmethod
    def classify_known_ligand_gene(
        gene: str,
        variants: List[Variant],
        ligand_info: Dict[str, Any],
        enrichment: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Classify a single known-ligand OR gene for tiered display.

        Returns:
            Dict with keys: gene, tier, ligand_zh, odor_zh, variants_by_tier, enrichment_summary
        """
        tier_a_variants = []  # HOM LoF/missense, rare HET LoF
        tier_b_variants = []  # Common HET protein-affecting
        tier_c_variants = []  # syn/UTR only

        for v in variants:
            csq = v.consequence or ""
            is_protein = csq in PROTEIN_CSQ
            is_hom = v.is_homozygous

            if is_hom and is_protein:
                tier_a_variants.append(v)
            elif not is_hom and is_protein:
                tier_b_variants.append(v)
            else:
                tier_c_variants.append(v)

        return {
            "gene": gene,
            "tier": "A" if tier_a_variants else ("B" if tier_b_variants else "C"),
            "ligand_zh": ligand_info.get("ligand_zh", ""),
            "odor_zh": ligand_info.get("odor_description_zh", ""),
            "tier_a_variants": tier_a_variants,
            "tier_b_variants": tier_b_variants,
            "tier_c_variants": tier_c_variants,
            "enrichment": enrichment,
        }

    @staticmethod
    def classify_unknown_or_layer(
        gene: str,
        variants: List[Variant],
        enrichment: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Determine which report layer an unknown-ligand OR gene falls into.

        Returns:
            Dict with layer (1=hide, 2=appendix, 3=main), variants, max_af, min_af
        """
        gnmd_variants = enrichment.get("gnomad_variants", [])
        af_values = []
        for gv in gnmd_variants:
            result = gv.get("result", {})
            af = result.get("gnomad_af")
            if af is not None:
                af_values.append(af)

        if not af_values:
            # No gnomAD data → treat conservatively (layer 2)
            return {"gene": gene, "layer": 2, "reason": "no_gnomad_data", "variants": variants}

        max_af = max(af_values)
        min_af = min(af_values)

        if min_af > 0.30:
            return {"gene": gene, "layer": 1, "reason": "all_common", "max_af": max_af, "variants": variants}
        elif min_af < 0.05:
            return {"gene": gene, "layer": 3, "reason": "has_rare", "min_af": min_af, "variants": variants}
        else:
            return {"gene": gene, "layer": 2, "reason": "moderate", "min_af": min_af, "max_af": max_af, "variants": variants}
