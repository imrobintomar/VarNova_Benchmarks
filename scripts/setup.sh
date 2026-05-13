#!/usr/bin/env bash
# =============================================================================
# VarNova Benchmark Setup Script
# Installs all required tools and downloads test data
# =============================================================================

set -euo pipefail

BENCHMARK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB_DIR="${DB_DIR:-$HOME/varnova_benchmark_dbs}"
TOOLS_DIR="${TOOLS_DIR:-$HOME/varnova_tools}"
VARNOVA_VERSION="${VARNOVA_VERSION:-latest}"

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

mkdir -p "$DB_DIR" "$TOOLS_DIR"

# ── 1. Install VarNova ────────────────────────────────────────────────────────
log "Installing VarNova..."
if ! command -v varnova &>/dev/null; then
    ARCH=$(uname -m)
    OS=$(uname -s | tr '[:upper:]' '[:lower:]')
    URL="https://github.com/imrobintomar/varnova/releases/latest/download/varnova-${OS}-${ARCH}.tar.gz"
    wget -q --show-progress -O /tmp/varnova.tar.gz "$URL" || \
        die "Cannot download VarNova binary. Build from source: https://github.com/imrobintomar/varnova"
    tar -xzf /tmp/varnova.tar.gz -C ~/.local/bin/ && chmod +x ~/.local/bin/varnova
    export PATH="$HOME/.local/bin:$PATH"
fi
varnova --version && log "VarNova: OK"

# ── 2. Install ANNOVAR ────────────────────────────────────────────────────────
log "Checking ANNOVAR..."
if [[ ! -f "$TOOLS_DIR/annovar/table_annovar.pl" ]]; then
    log "ANNOVAR requires free registration at https://annovar.openbioinformatics.org/"
    log "Download annovar.latest.tar.gz, place in $TOOLS_DIR/ and re-run setup"
    log "Skipping ANNOVAR for now (VarNova-only benchmark will still work)"
fi

# ── 3. Install VEP ───────────────────────────────────────────────────────────
log "Checking VEP..."
if ! conda run -n vep vep --help &>/dev/null 2>&1; then
    log "Installing VEP via conda..."
    conda create -n vep -c bioconda -c conda-forge ensembl-vep -y 2>/dev/null || true
fi

# ── 4. Install SnpEff ────────────────────────────────────────────────────────
log "Checking SnpEff..."
if [[ ! -f "$TOOLS_DIR/snpEff/snpEff.jar" ]]; then
    log "Downloading SnpEff..."
    wget -q --show-progress -O /tmp/snpeff.zip \
        "https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip"
    unzip -q /tmp/snpeff.zip -d "$TOOLS_DIR/"
    java -jar "$TOOLS_DIR/snpEff/snpEff.jar" download hg38 -noLog 2>/dev/null
fi

# ── 5. Download test VCF (1000 Genomes sample, ~500 variants) ────────────────
log "Downloading test VCF..."
TEST_VCF="$BENCHMARK_DIR/testdata/test_exome_hg38.vcf"
if [[ ! -f "$TEST_VCF" ]]; then
    # Small public test VCF from 1000G project
    wget -q --show-progress -O "$TEST_VCF.gz" \
        "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/exome_alignment/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20121211.bam.vcf.gz" \
        2>/dev/null || \
    # Fallback: create minimal test VCF
    python3 "$BENCHMARK_DIR/scripts/generate_test_vcf.py" --output "$TEST_VCF" --variants 1000
    [[ -f "${TEST_VCF}.gz" ]] && gunzip "${TEST_VCF}.gz"
fi
log "Test VCF: $(wc -l < "$TEST_VCF") lines"

# ── 6. Download databases ─────────────────────────────────────────────────────
log "Database download options:"
echo ""
echo "  Option A (quick): Use pre-built .vnidx databases (VarNova only)"
echo "    bash $BENCHMARK_DIR/scripts/download_vnidx_dbs.sh $DB_DIR"
echo ""
echo "  Option B (full): Download all ANNOVAR databases (~60 GB)"
echo "    perl annovar/annotate_variation.pl --downdb refGene $DB_DIR --buildver hg38"
echo "    perl annovar/annotate_variation.pl --downdb gnomad41_exome $DB_DIR --buildver hg38"
echo "    perl annovar/annotate_variation.pl --downdb clinvar_20240730 $DB_DIR --buildver hg38"
echo "    perl annovar/annotate_variation.pl --downdb avsnp151 $DB_DIR --buildver hg38"
echo ""

log "Setup complete. Edit config.sh then run: bash scripts/run_benchmark.sh"
