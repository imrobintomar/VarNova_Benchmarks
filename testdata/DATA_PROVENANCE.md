---
title: Data provenance for this repo's testdata/ and download scripts
---

Everything anyone needs to fully reproduce every tier in this repo, without
re-hosting any third-party database here. Files actually committed in this
repo are small (test/gold-standard VCFs); everything else is fetched on
demand by the scripts listed below.

## Committed in this repo (small, no download needed)

| Path | Source | Size |
|---|---|---|
| `testdata/clinvar_gold/clinvar_gold.vcf` | Derived from NCBI ClinVar `variant_summary.txt`, filtered to Pathogenic/Likely_pathogenic, high-confidence review status, GRCh38 | ~3MB |
| `testdata/clinvar_gold/clinvar_gold_truth.tsv` | Same filter, columns: Chrom/Pos/Ref/Alt/Gene/Transcript/HGVSc/HGVSp/ClinicalSignificance/ReviewStatus | ~12MB |

## Fetched on demand (large — not committed)

| What | Script | Approx size | Source |
|---|---|---|---|
| Reference DBs (refGene, gnomAD 4.1 exome, ClinVar, avsnp151, dbNSFP 4.7a) | `scripts/download_reference_dbs.sh DB_DIR` | ~90GB | ANNOVAR `--downdb` (annovar.openbioinformatics.org mirror) |
| GIAB HG001-HG007 v4.2.1 GRCh38 benchmark VCFs | `scripts/giab_benchmark.sh` | ~700MB-1.6GB each | NIST GIAB FTP (ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/) — verify exact subpaths against the live listing, NIST periodically revises release layouts |
| Synthetic test VCFs (any size, 1K-100M variants) | `scripts/generate_test_vcf.py --variants N` | generated locally | N/A — synthetic, seeded from real known ClinVar/1000G variants |
| VarNova binary | GitHub Releases | ~1.4MB compressed | This repo's own Releases page |

## Why nothing large is committed here

GitHub repos aren't a practical home for hundreds of GB of third-party reference
data, and re-hosting gnomAD/ClinVar/dbSNP/GIAB ourselves duplicates data that
already has a canonical, versioned, authoritative source. Every script above
downloads directly from that canonical source so anyone re-running a benchmark
gets the same (or a clearly-dated newer) version NCBI/NIST actually publishes.
