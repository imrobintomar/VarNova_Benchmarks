#!/usr/bin/env bash
# =============================================================================
# VarNova Benchmark Configuration
# Edit these paths to match your system before running the benchmark
# =============================================================================

# ── Input ─────────────────────────────────────────────────────────────────────
# VCF file to benchmark (exome, hg38)
INPUT_VCF="testdata/test_exome_hg38.vcf"
# Set to number of variants in VCF (used for speed calculations)
VARIANTS="auto"   # or set manually, e.g. 89052

# ── Database directory (ANNOVAR flat .txt files) ──────────────────────────────
DB="/path/to/humandb"               # e.g. /home/user/annovar_humandb

# ── Binary database cache (optional, for max speed) ──────────────────────────
# Convert DBs with: varnova convert $DB/hg38_gnomad41_exome.txt --out-dir $VNIDX_CACHE
VNIDX_CACHE="$HOME/varnova_db"     # set to "" to disable

# ── Tool paths ────────────────────────────────────────────────────────────────
VARNOVA="varnova"                                        # if in PATH
ANNOVAR_DIR="/path/to/annovar"                           # directory with table_annovar.pl
SNPEFF_JAR="/path/to/snpEff/snpEff.jar"
SNPSIFT_JAR="/path/to/snpEff/SnpSift.jar"

# ── VEP cache ─────────────────────────────────────────────────────────────────
VEP_CACHE="$HOME/.vep"             # VEP cache directory
VEP_PLUGIN_DIR=""                  # leave empty if no plugins

# ── Optional: dbNSFP for VEP/SnpSift plugins ─────────────────────────────────
DBNSFP_BGZ=""                      # path to dbNSFP4.7a_combined.gz (optional)
GNOMAD_BGZ=""                      # path to gnomad.exome.v4.1.vcf.bgz (optional)
CLINVAR_BGZ=""                     # path to clinvar.vcf.bgz (optional)

# ── Output ────────────────────────────────────────────────────────────────────
OUT_DIR="results/$(date +%Y%m%d_%H%M%S)"
THREADS=$(nproc)

# ── Which tools to benchmark ──────────────────────────────────────────────────
RUN_VARNOVA=true
RUN_ANNOVAR=true
RUN_VEP=true
RUN_SNPEFF=true
