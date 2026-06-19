# Sensory Genomics &nbsp;·&nbsp; 五感基因组分析

[![version](https://img.shields.io/badge/version-0.1.0-blue)](https://github.com/lzr098/sensory-genomics-skill)
[![python](https://img.shields.io/badge/python-3.10%2B-green)](https://www.python.org/)
[![license](https://img.shields.io/badge/license-MIT-lightgrey)](LICENSE)

> **Genome-to-senses report.** Analyze WGS/WES VCF files for genetic variants affecting the five sensory systems: vision, hearing, olfaction, taste, and somatosensation. Outputs structured Markdown + JSON reports with five-level impact assessment.

> **从基因组到感官报告。** 分析 WGS/WES VCF 中影响视觉、听觉、嗅觉、味觉、体感五大系统的遗传变异，输出五级影响评估的结构化报告。

---

## What It Does

- Processes germline WGS/WES VCF files through a 14-stage pipeline
- **3-tier VEP annotation**: memory cache → precomputed variant DB (19,861 records) → local Docker VEP → REST API
- **376 sensory genes** across 5 subsystems (vision, hearing, olfaction, taste, somatosensation)
- **5-level impact assessment** (complete loss → significant → partial → minor → none)
- **Specialized modules**: TAS2R38 bitter taste haplotypes, mitochondrial deafness (MT-RNR1), OR gene tiered classification
- **Key trait SNP inference**: 50+ trait-associated SNPs with REF/REF fallback for absent positions
- **Personal trait prediction**: 8 traits (eye/hair/skin color, caffeine/alcohol/lactose metabolism, hair texture, earwax type)
- **Ethics**: Mandatory disclaimers, no medical diagnosis language

## 做了什么

- 14 阶段管线处理 germline WGS/WES VCF
- **三级 VEP 注释**：内存缓存 → 预计算变异数据库 (19,861 条) → 本地 Docker VEP → REST API
- **376 个感官基因** 覆盖 5 个子系统
- **五级影响评估**（完全丧失 → 显著 → 部分 → 轻微 → 无）
- **专项模块**：TAS2R38 苦味单倍型、线粒体耳聋 (MT-RNR1)、嗅觉受体基因分级
- **关键性状 SNP 推断**：50+ SNPs，缺失位点自动 REF/REF 推断
- **个人性状预测**：8 个性状（眼/发/肤色、咖啡因/酒精/乳糖代谢、发质、耳垢类型）
- **伦理约束**：强制免责声明，禁止医学诊断语言

---

## Quick Start · 快速开始

```bash
# Basic usage (auto-detect sex, strict filter)
python scripts/main.py --vcf sample.vcf.gz --auto-sex --strict-filter

# Specify sex explicitly
python scripts/main.py --vcf sample.vcf.gz --sex M --strict-filter

# Skip API enrichment for faster run
python scripts/main.py --vcf sample.vcf.gz --auto-sex --no-reference-info
```

---

## Five-Level Impact Assessment · 五级影响评估

| Level | Protein Impact | Gene Certainty | Inheritance |
|---|---|---|---|
| **Complete Loss** | frameshift / stop-gained in single-function gene | High certainty gene | Homozygous / compound het |
| **Significant** | splice site / missense in critical domain | High/medium | Matches pattern |
| **Partial** | missense in flexible region | Medium | Heterozygous AD |
| **Minor** | missense with high population AF | Low | Mismatched pattern |
| **None** | synonymous / UTR / intron | Any | Any |

---

## Specialized Modules · 专项模块

### TAS2R38 Bitter Taste · 苦味感知

Analyzes 3 SNPs (rs713598, rs1726866, rs10246939) → diplotype → bitter taste phenotype:
- **PAV/PAV** → supertaster
- **PAV/AVI** → medium taster
- **AVI/AVI** → non-taster

### Mitochondrial Deafness · 线粒体耳聋

Detects MT-RNR1 / MT-TS1 variants (m.1555A>G, etc.) → outputs aminoglycoside antibiotic warning.

### OR Gene Classification · 嗅觉受体分级

- **Tier A** — Known ligand ORs with homozygous LoF → main report
- **Tier B** — Functionally studied ORs with homozygous LoF → table
- **Tier C** — Other functional ORs with homozygous LoF → appendix

---

## Pipeline Stages · 管线阶段

```
1. VCF validation & sex detection
2. VCF parsing (streaming, .vcf / .vcf.gz)
3. Quality prefilter (QUAL ≥ 30, DP ≥ 10)
4. Strict-filter (BED exon filtering, optional)
5. Tiered VEP annotation (3-level lookup)
6. Transcript arbitration
7. Sensory gene set filter (376 genes, 5 subsystems)
8. Functional impact assessment
9. Specialized modules (TAS2R38, mito, OR)
10. Key trait SNP inference
11. API enrichment (UniProt / gnomAD / ClinVar / GTEx)
12. Personal trait prediction (8 traits)
13. Cross-check & downgrade
14. Report generation (Markdown + JSON)
```

---

## Precomputed Reference Data · 预计算参考数据

| Asset | Size | Records |
|---|---|---|
| Sensory genes exons BED | — | 6,705 intervals |
| Sensory genes CDS BED | — | 5,758 intervals |
| Gene → transcript map | — | 229 genes |
| Variant VEP SQLite DB | 80 MB | 19,861 records |

---

## Performance · 性能

| Metric | Target |
|---|---|
| Analysis (precompute DB hit) | ~5-7 min |
| Analysis (cold cache) | ~8-12 min |
| VEP API batches (precompute) | ~30-40 |
| Memory usage | < 2 GB |
| Report size | 50-300 KB |

---

## Requirements · 运行环境

- Python 3.10+
- System Python 3.14 on macOS (pysam code-signing constraint)
- `pysam`, `aiohttp`, `aiosqlite`, `pydantic`, `Jinja2`, `PyYAML`
- Docker (optional, for local VEP offline annotation)
- VEP Docker cache: `~/.workbuddy/tools/vep/cache/homo_sapiens/115_GRCh38/`

---

## Ethical Constraints · 伦理约束

- ❌ No diagnostic language ("pathogenic", "diagnosis", "patient")
- ✅ OR gene LoF → "specific olfactory receptor loss" (not disease)
- ✅ SCN9A / PIEZO2 findings → genetic counseling recommended
- ✅ Mandatory disclaimer on first page of every report

---

## Known Limitations · 已知限制

- SNVs/indels only — no CNV or SV detection
- OPN1LW/MW CNV rearrangements (red-green color weakness) not detected
- OR gene ligand mapping incomplete
- Pain/temperature sensitivity subject to psychological factors
- Precompute DB only covers common variants (AF ≥ 0.1%)

---

## Related Skills · 相关技能

| Skill | Repo | Purpose |
|---|---|---|
| **GPA** | [lzr098/dgra-genomic-risk](https://github.com/lzr098/dgra-genomic-risk) | Whole-genome phenotype association |
| **variant-impact** | [lzr098/variant-impact](https://github.com/lzr098/variant-impact) | Single variant ACMG classification |
| **disease-risk-query** | [lzr098/Disease-Risk-Query](https://github.com/lzr098/Disease-Risk-Query) | Disease-specific risk analysis |

---

## License

MIT

---

**Maintainer**: [@lzr098](https://github.com/lzr098)
