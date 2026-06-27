#!/usr/bin/env bash
# Downloads every ANNOVAR-format reference database needed across all
# benchmark tiers in this repo. ~90GB total — opt-in, not run automatically
# by setup.sh, since a download this large shouldn't happen without the user
# explicitly asking for it.
#
# Required for:
#   Tier 1/2 (gene + full pipeline)  -> refGene, gnomad41_exome, clinvar, avsnp151
#   Tier 5   (ClinVar gold-standard) -> refGene, gnomad41_exome, dbnsfp47a
#
# Usage: bash scripts/download_reference_dbs.sh /path/to/humandb
set -euo pipefail

DB="${1:?Usage: bash scripts/download_reference_dbs.sh /path/to/humandb}"
ANNOVAR_DIR="${ANNOVAR_DIR:-$HOME/varnova_tools/annovar}"

if [[ ! -f "$ANNOVAR_DIR/annotate_variation.pl" ]]; then
  echo "ANNOVAR not found at $ANNOVAR_DIR/annotate_variation.pl" >&2
  echo "ANNOVAR requires free registration: https://annovar.openbioinformatics.org/" >&2
  echo "Download annovar.latest.tar.gz, extract to $ANNOVAR_DIR, and re-run." >&2
  exit 1
fi

mkdir -p "$DB"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

download() {
  local name="$1"
  if ls "$DB"/hg38_"${name}"* &>/dev/null; then
    log "$name already present in $DB — skipping"
    return
  fi
  log "Downloading $name (this can take a while)..."
  perl "$ANNOVAR_DIR/annotate_variation.pl" --downdb -webfrom annovar "$name" "$DB" --buildver hg38
}

log "== refGene (~50MB) =="
download refGene

log "== gnomAD 4.1 exome (~18GB) =="
download gnomad41_exome

log "== ClinVar, latest available (~1GB) =="
latest_clinvar=$(curl -s "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/" 2>/dev/null \
  | grep -oP 'clinvar_\d{8}' | sort -u | tail -1 || echo "clinvar_20240730")
log "Using $latest_clinvar (verify this is current before relying on it)"
download "$latest_clinvar"

log "== avsnp151 / dbSNP (~29GB) =="
download avsnp151

log "== dbNSFP 4.7a — needed for Tier 5 PP3/BP4/PS1/PM5 computational evidence (~25GB) =="
download dbnsfp47a

log "All reference databases ready in $DB"
log "Point config.sh's DB= variable at this directory before running any tier."
