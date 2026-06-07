---
name: sensory-genomics
description: Analyze WGS/WES VCF files to generate genetic reports on an individual's five sensory systems (vision, hearing, olfaction, taste, somatosensation). Trigger when the user uploads a VCF file and asks about sensory abilities, taste perception, color vision, hearing genetics, olfactory receptors, or pain/temperature sensitivity. Also trigger for requests involving TAS2R38 haplotype analysis, mitochondrial deafness variants (MT-RNR1), OR gene classification, or sensory-related genetic variant interpretation.
agent_created: true
---

# Sensory Genomics Analysis Skill

## Overview

Analyze individual genomic variants from WGS/WES VCF files to assess the functional impact on five sensory systems. Output structured Markdown reports and JSON data with a 5-level impact assessment scale (complete loss → no impact). This skill performs no medical diagnosis — all outputs use descriptive language with mandatory ethical disclaimers.

## When to Use This Skill

Use this skill when the user:
- Uploads or references a VCF file and asks about sensory genetics
- Wants to know about bitter taste perception (TAS2R38 haplotypes)
- Asks about color vision genetics or red-green color weakness risk
- Wants hearing-related genetic variant analysis (including mitochondrial deafness)
- Asks about olfactory receptor (OR) gene variants or "what smells they might not detect"
- Requests pain sensitivity, temperature perception, or mechanosensation genetic analysis
- Mentions specific genes like GJB2, OPN1LW, SCN9A, PIEZO2, TRPV1, etc.

## Core Capabilities

### 1. VCF Processing Pipeline

To process a VCF file through the full pipeline:

1. Validate the VCF path and sample sex (M/F)
2. Parse the VCF using streaming read (supports `.vcf` and `.vcf.gz`)
3. Apply quality prefilter: QUAL ≥ 30, DP ≥ 10, FILTER = PASS / `.` / `""`
4. **Strict-filter** (optional): use precomputed exon BED to exclude intronic variants before VEP
5. Annotate variants via **tiered VEP lookup** (auto-detect by default):
   - **L1** — Memory cache (SQLite, 90-day TTL)
   - **L2** — Precomputed variant VEP database (`assets/data/sensory_variants_vep.sqlite`, 19,861 records)
   - **L3** — Local Docker VEP (`auto`/`local_docker` mode): offline annotation with full GRCh38 cache (~24 GB)
   - **L4** — VEP REST API POST batch (`rest_api` mode): 200 variants/request, fallback for novel variants

   Default `vep.source: "auto"` uses `HybridVepClient` — automatically detects Docker + cache availability, prefers local, seamlessly falls back to REST API if unavailable.
6. Transcript arbitration: score by canonical(+10) + protein_coding(+5) + consequence severity(+0-40) to resolve gene overlaps (e.g. NLRP3/OR2B11)
7. Filter to sensory gene sets with exact gene-symbol matching
8. Run functional impact assessment
9. Apply specialized logic modules
10. Enrich with UniProt / gnomAD / ClinVar / GTEx (cached, display-only)
11. Generate Markdown report + JSON output

Use `scripts/main.py` as the pipeline entry point. Refer to `references/config.yaml` for default parameters.

### 2. Five-Level Impact Assessment

Assess each coding variant across three dimensions:

**Protein Impact:**
- frameshift / stop-gained → "complete_loss"
- splice site → "severe_damage"
- missense in critical domain → "moderate_change" or "severe_damage"
- missense in flexible loop → "minor_change"
- synonymous / UTR / intron → "none"

**Gene Certainty:**
- high: single-function genes with strong genetic evidence (GJB2, OPN1LW, SCN9A)
- medium: genes with partial evidence (TAS2R38, PIEZO2)
- low: speculative function or high polymorphism (most OR genes, TRPV1 missense)

**Inheritance Pattern Match:**
- recessive: requires homozygous LoF or compound heterozygote
- X-linked recessive: male hemizygous or female homozygous
- dominant: heterozygous sufficient
- mitochondrial: heteroplasmy threshold > 50-90%

Combine into final level: complete loss / significant / partial / minor / none.

### 3. Specialized Logic Modules

**TAS2R38 Haplotype Analyzer**
- Located in `scripts/specialized/tas2r38.py`
- Input: 3 SNP genotypes (rs713598, rs1726866, rs10246939)
- Output: diplotype (PAV/PAV, PAV/AVI, AVI/AVI, etc.) + bitter taste phenotype
- Reference data: `assets/data/tas2r38_snps.json`

**Mitochondrial Deafness Annotator**
- Located in `scripts/specialized/mitochondrial.py`
- Detects known MT-RNR1 / MT-TS1 variants (m.1555A>G, etc.)
- Outputs drug warning: "Avoid aminoglycoside antibiotics"
- Reference data: `assets/data/mitochondrial_variants.json`

**OR Gene Tiered Display**
- Located in `scripts/specialized/or_tiers.py`
- Tier A: Known ligand ORs with homozygous LoF → highlight in main report
- Tier B: Functionally studied ORs with homozygous LoF → table list
- Tier C: Other functional ORs with homozygous LoF → appendix
- Reference data: `assets/data/or_ligands.json`

### 4. Precomputed Reference Data

**Gene Coordinates** (`scripts/precompute.py`)
- Fetched from Ensembl REST API for all 376 sensory genes across 5 subsystems
- 367/376 genes successfully resolved (OR13C2, OR2T27 failed)
- Outputs:
  - `assets/data/sensory_genes_exons.bed` — 6,705 exon intervals (0-based half-open)
  - `assets/data/sensory_genes_cds.bed` — 5,758 CDS intervals
  - `assets/data/gene_transcript_map.json` — gene → canonical transcript ID (229 genes)

**Variant VEP Database** (`scripts/precompute_vep_db.py`)
- Queries gnomAD r4 GraphQL API per BED region, filters by AF ≥ 0.1%
- Annotates retained variants via VEP POST batch
- Stored in `assets/data/sensory_variants_vep.sqlite` (19,861 records, 80 MB)
- At runtime, looked up before VEP API calls → ~75-80% hit rate for common variants

### 5. Report Generation

Generate two output formats simultaneously:

**Markdown Report** (`scripts/report/markdown_generator.py`)
- Uses Jinja2 templates from `assets/templates/`
- Sections: disclaimer → executive summary → 5 sensory chapters → methodology → appendix
- Each gene card contains: variant info, impact assessment, rationale, inheritance, protein topology, reference data

**JSON Output** (`scripts/report/json_generator.py`)
- Machine-readable structured data
- Contains all variants, assessments, metadata

## File Organization

```
sensory-genomics/
├── SKILL.md                          # This file
├── scripts/                          # Python source code
│   ├── main.py                       # Pipeline entry point
│   ├── models.py                     # Pydantic data models
│   ├── config_loader.py              # Configuration loading
│   ├── precompute.py                 # Offline gene coordinate fetcher (Ensembl API)
│   ├── precompute_vep_db.py          # Offline variant VEP precomputer (gnomAD + VEP)
│   ├── vcf/                          # VCF parser and prefilter
│   ├── vep/                          # VEP REST API client (3-tier lookup)
│   ├── gene_sets/                    # Sensory gene set loader/filter (bisect BED)
│   ├── assessment/                   # Impact assessment engine (cached inheritance)
│   ├── specialized/                  # TAS2R38, mitochondrial, OR tiers
│   ├── enrichment/                   # API clients + SQLite cache
│   └── report/                       # Markdown/JSON generators
├── assets/
│   ├── data/                         # Gene sets, SNP definitions, variant DBs, precomputed BED/SQLite
│   │   ├── sensory_gene_sets.yaml
│   │   ├── sensory_genes_exons.bed   # Precomputed exon coordinates (6,705 intervals)
│   │   ├── sensory_genes_cds.bed     # Precomputed CDS coordinates (5,758 intervals)
│   │   ├── gene_transcript_map.json  # Gene → canonical transcript ID
│   │   ├── sensory_variants_vep.sqlite  # Precomputed VEP annotations (19,861 records)
│   │   └── ...
│   └── templates/                    # Jinja2 report templates
└── references/
    ├── prd.md                        # Product requirements document
    ├── architecture.md               # System architecture design
    ├── requirements.txt              # Python dependencies
    └── config.yaml                   # Default runtime configuration
```

## Dependencies

Install dependencies from `references/requirements.txt`:

```
pysam>=0.22.0
aiohttp>=3.9.0
aiosqlite>=0.20.0
pydantic>=2.5.0
Jinja2>=3.1.0
PyYAML>=6.0.1
```

### Local Docker VEP (Optional, Zero-Config)

With `vep.source: "auto"` (default), the pipeline **automatically detects** local Docker VEP and uses it when available. No manual config needed.

To prepare local VEP:

1. **Install Docker Desktop** (macOS/Windows/Linux)
2. **Pull VEP image**: `docker pull ensemblorg/ensembl-vep:latest`
3. **Download GRCh38 cache** (~24 GB):
   ```bash
   curl -C - -O \
     "http://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz"
   tar -xzf homo_sapiens_vep_115_GRCh38.tar.gz
   ```
   Cache path: `~/.workbuddy/tools/vep/cache/homo_sapiens/115_GRCh38/`

Force a specific source via `config.yaml`:
- `"auto"` — auto-detect local, fallback REST (default)
- `"local_docker"` — always use local Docker
- `"rest_api"` — always use REST API

## Usage Example

```python
import asyncio
from scripts.main import SensoryPipeline
from scripts.config_loader import ConfigLoader

async def analyze():
    config = ConfigLoader().load()
    pipeline = SensoryPipeline(config)
    report = await pipeline.run(
        vcf_path="sample.vcf.gz",
        sex="M",
        subsystems=["vision", "hearing", "olfaction", "taste", "somatosensation"]
    )
    report.save_markdown("output/report.md")
    report.save_json("output/report.json")

asyncio.run(analyze())
```

## Ethical Constraints

- Never use diagnostic language ("致病", "患者", "诊断", "pathogenic", "patient", "diagnosis")
- OR gene LoF must be labeled as "specific olfactory receptor loss" not disease
- SCN9A / PIEZO2 positive findings must include "建议遗传咨询" / "genetic counseling recommended"
- Mandatory disclaimer must appear on the first page of every report
- Do not make deterministic predictions for minors
- Do not causally link TAS2R38 variants to eating behavior — only report bitter perception ability

## Performance Targets

| Metric | Target | Notes |
|--------|--------|-------|
| Analysis time (precompute DB hit) | ~5-7 min | 75-80% variants resolved locally |
| Analysis time (cold cache) | ~8-12 min | All VEP via REST API |
| VEP API batches (precompute) | ~30-40 | vs ~150 without precompute DB |
| VEP API batches (cold) | ~150 | Full REST API for 21k variants |
| Memory usage | < 2 GB | |
| Report size | 50–300 KB | |

**Optimization highlights:** BED bisect O(log n) filtering, per-gene inheritance/certainty caching, Pydantic v2 `model_copy()`, flattened enrichment API concurrency, logger handler deduplication.

## Known Limitations

- v0.1.0 analyzes SNVs/indels only — no CNV or SV detection
- OPN1LW/MW CNV rearrangements (primary cause of red-green color weakness) are not detected
- OR gene ligand mapping is incomplete — most OR losses cannot be linked to specific odors
- Taste perception involves polygenic and environmental factors not fully captured
- Pain/temperature sensitivity varies greatly due to psychological and developmental factors
- Precompute DB only covers common variants (AF ≥ 0.1% in gnomAD r4); novel/private variants still require VEP API
- 329 of 20,190 filtered variants failed VEP annotation during precompute (1.6%, mostly API timeouts)
- Ensembl overlap API does not return UTR coordinates; UTR variants rely on VEP runtime annotation

## Runtime Troubleshooting

### macOS pysam Code-Signing Issue

On macOS, the `pysam` package in the default WorkBuddy venv may fail to import due to Team ID code-signing conflicts ("The process has forked and you cannot use this CoreFoundation functionality safely"). Do **not** attempt to reinstall pysam in the default venv — codesign, xattr, and source compilation have all been tried and do not resolve the issue.

**Workaround**: Use the system Python 3.14 runtime instead:
```bash
/Library/Frameworks/Python.framework/Versions/3.14/bin/python3 -m pip install aiosqlite pydantic Jinja2 PyYAML
export PYTHONPATH="/path/to/sensory-genomics/scripts:$PYTHONPATH"
/Library/Frameworks/Python.framework/Versions/3.14/bin/python3 -m src.main --vcf ...
```

### Reference Genome Version Mismatch (GRCh37 vs GRCh38)

The `assets/data/key_trait_snps.yaml` file defines trait-associated SNPs for the Key SNP Inference module. Historically these coordinates were stored in **GRCh37** but modern WGS/WES VCF files are almost always **GRCh38**. A coordinate mismatch causes all SNP lookups to miss, resulting in 100% REF/REF fallback inference — which silently produces incorrect phenotype predictions.

**Prevention**:
- All SNP coordinates in `key_trait_snps.yaml` must be verified as GRCh38 before deployment. The included SNPs were batch-updated via Ensembl API on 2026-06-07 (commit `48cf785`).
- When adding new SNPs, query Ensembl `/variation/human/{rsid}` and use the `GRCh38` mapping (`seq_region_name`, `start`, `allele_string`).

### Phenotype Map Key Normalization

When a SNP's reference allele changes (e.g., due to GRCh37→GRCh38 strand flips or allele reorientation), the genotype keys in `phenotype_map` may no longer match the VCF genotype string. Additionally, heterozygous genotypes like `CT` and `TC` are semantically identical but fail exact string matching.

**Fix**: The `key_snps.py` module now normalizes genotype keys by sorting alleles alphabetically before phenotype lookup:
```python
sorted_pmap = { "".join(sorted(k)): v for k, v in pmap.items() }
pheno = sorted_pmap.get("".join(sorted(gt)), {})
```

### Multi-Allelic SNP Coverage

Some trait SNPs are multi-allelic (e.g., rs1229984 / ADH1B has three alleles in GRCh38). The original `phenotype_map` only covered bi-allelic combinations. When a VCF reports a genotype involving a secondary alt allele (e.g., `CT` where `T` is the primary ref and `C` is alt[1]), it falls through to "unknown phenotype".

**Fix**: Ensure `phenotype_map` covers all combinatorial genotypes for multi-allelic SNPs after coordinate updates. After the 2026-06-07 fix, all 46 key SNPs resolve to defined phenotypes with zero "unknown" entries.
