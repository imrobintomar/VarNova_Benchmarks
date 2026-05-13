#!/usr/bin/env bash
# =============================================================================
# VarNova Fair Benchmark
#
# Benchmarks VarNova against ANNOVAR, VEP, and SnpEff on:
#   TIER 1 — Gene annotation only (all tools)
#   TIER 2 — Full pipeline (same databases across all tools)
#
# Usage:
#   bash scripts/run_benchmark.sh             # uses config.sh
#   bash scripts/run_benchmark.sh config.sh   # custom config
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG="${1:-$REPO_DIR/config.sh}"

[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG"; exit 1; }
source "$CONFIG"
export PATH="$HOME/.local/bin:$HOME/bin:$PATH"

# Auto-detect variant count
if [[ "$VARIANTS" == "auto" ]]; then
    VARIANTS=$(grep -v "^#" "$INPUT_VCF" | wc -l)
fi

mkdir -p "$OUT_DIR"
REPORT="$OUT_DIR/benchmark_report.txt"
CSV="$OUT_DIR/results.csv"

log()  { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$REPORT"; }
sep()  { echo "─────────────────────────────────────────" | tee -a "$REPORT"; }

run_timed() {
    local label="$1" outfile="$2"; shift 2
    local timefile; timefile=$(mktemp)
    log "Running: $label"
    log "Command: $*"
    /usr/bin/time -v -o "$timefile" "$@" >> "$outfile" 2>&1 || true
    local wall cpu ram
    wall=$(grep "Elapsed (wall clock)" "$timefile" | grep -oP '\d+:\d+\.\d+$' | \
           awk -F: '{printf "%.1f",$1*60+$2}' 2>/dev/null || echo "N/A")
    cpu=$(grep "User time" "$timefile" | grep -oP '[\d.]+$' || echo "N/A")
    ram=$(grep "Maximum resident" "$timefile" | grep -oP '\d+$' || echo "N/A")
    local vps="N/A"
    [[ "$wall" != "N/A" && "$wall" != "0" ]] && \
        vps=$(echo "scale=0; $VARIANTS / $wall" | bc 2>/dev/null || echo "N/A")
    rm -f "$timefile"
    echo "$label|$wall|$cpu|${ram}kB|$vps" >> "$CSV"
    log "  Wall: ${wall}s | CPU: ${cpu}s | RAM: $((${ram:-0}/1024))MB | Speed: ${vps} v/s"
    sep
}

# ── Header ────────────────────────────────────────────────────────────────────
{
echo "========================================================"
echo " VarNova Fair Benchmark"
echo " Date    : $(date)"
echo " Machine : $(nproc) cores | $(free -h | awk '/Mem/{print $2}') RAM"
echo " Input   : $INPUT_VCF ($VARIANTS variants)"
echo " Threads : $THREADS"
echo "========================================================"
} | tee "$REPORT"

echo "Tool|Wall(s)|CPU(s)|RAM|Variants/sec" > "$CSV"
sep

# =============================================================================
# TIER 1 — GENE ANNOTATION ONLY
# =============================================================================
log "══════ TIER 1: GENE ANNOTATION ONLY ══════"
sep

# VarNova gene
if $RUN_VARNOVA; then
    run_timed "VarNova_gene" "$OUT_DIR/varnova_gene.tsv" \
        $VARNOVA gene \
            -i "$INPUT_VCF" \
            -d "$DB/hg38_refGene.txt" \
            -o "$OUT_DIR/varnova_gene.tsv" \
            -t "$THREADS"
fi

# ANNOVAR gene
if $RUN_ANNOVAR && [[ -f "$ANNOVAR_DIR/annotate_variation.pl" ]]; then
    perl "$ANNOVAR_DIR/convert2annovar.pl" \
        -format vcf4 "$INPUT_VCF" \
        -outfile "$OUT_DIR/input.avinput" \
        -includeinfo 2>/dev/null || true
    run_timed "ANNOVAR_gene" "$OUT_DIR/annovar_gene.out" \
        perl "$ANNOVAR_DIR/annotate_variation.pl" \
            -geneanno -buildver hg38 -dbtype refGene \
            -outfile "$OUT_DIR/annovar_gene" \
            "$OUT_DIR/input.avinput" "$DB"
fi

# VEP gene
if $RUN_VEP && conda run -n vep vep --help &>/dev/null 2>&1; then
    rm -f "$OUT_DIR/vep_gene.vcf"
    run_timed "VEP_gene" "$OUT_DIR/vep_gene.vcf" \
        conda run -n vep vep \
            --input_file "$INPUT_VCF" --output_file "$OUT_DIR/vep_gene.vcf" \
            --format vcf --vcf --offline --fork "$THREADS" \
            --cache --dir_cache "$VEP_CACHE" --assembly GRCh38 \
            --no_stats --force_overwrite
fi

# SnpEff gene
if $RUN_SNPEFF && [[ -f "$SNPEFF_JAR" ]]; then
    run_timed "SnpEff_gene" "$OUT_DIR/snpeff_gene.vcf" \
        java -Xmx16g -jar "$SNPEFF_JAR" ann -v hg38 -noStats "$INPUT_VCF"
fi

# =============================================================================
# TIER 2 — FULL PIPELINE
# =============================================================================
log "══════ TIER 2: FULL PIPELINE ══════"
log " Databases: refGene + gnomAD 4.1 exome + ClinVar 2024 + dbSNP 151"
sep

# VarNova full
if $RUN_VARNOVA; then
    CACHE_ARG=""
    [[ -n "$VNIDX_CACHE" && -d "$VNIDX_CACHE" ]] && CACHE_ARG="--cache-dir $VNIDX_CACHE"
    run_timed "VarNova_full" "$OUT_DIR/varnova_full.log" \
        $VARNOVA table \
            -i "$INPUT_VCF" \
            --gene-db   "$DB/hg38_refGene.txt" \
            --filter-db "$DB/hg38_gnomad41_exome.txt" \
            --filter-db "$DB/hg38_clinvar_20240730.txt" \
            --filter-db "$DB/hg38_avsnp151.txt" \
            --label gnomAD --label ClinVar --label rsID \
            $CACHE_ARG \
            --out-dir "$OUT_DIR" \
            -t "$THREADS" -v
fi

# ANNOVAR full
if $RUN_ANNOVAR && [[ -f "$ANNOVAR_DIR/table_annovar.pl" ]]; then
    run_timed "ANNOVAR_full" "$OUT_DIR/annovar_full.out" \
        perl "$ANNOVAR_DIR/table_annovar.pl" \
            "$INPUT_VCF" "$DB" \
            -buildver hg38 \
            -out "$OUT_DIR/annovar_full" \
            -remove \
            -protocol refGene,gnomad41_exome,clinvar_20240730,avsnp151 \
            -operation g,f,f,f \
            -vcfinput -nastring .
fi

# VEP full (with plugins if available)
if $RUN_VEP && conda run -n vep vep --help &>/dev/null 2>&1; then
    rm -f "$OUT_DIR/vep_full.vcf"
    VEP_EXTRA=""
    [[ -n "$DBNSFP_BGZ" && -f "$DBNSFP_BGZ" ]] && \
        VEP_EXTRA="$VEP_EXTRA --plugin dbNSFP,$DBNSFP_BGZ,SIFT_score,CADD_phred"
    [[ -n "$GNOMAD_BGZ" && -f "$GNOMAD_BGZ" ]] && \
        VEP_EXTRA="$VEP_EXTRA --custom $GNOMAD_BGZ,gnomAD,vcf,exact,0,AF"
    [[ -n "$CLINVAR_BGZ" && -f "$CLINVAR_BGZ" ]] && \
        VEP_EXTRA="$VEP_EXTRA --custom $CLINVAR_BGZ,ClinVar,vcf,exact,0,CLNSIG"
    run_timed "VEP_full" "$OUT_DIR/vep_full.vcf" \
        conda run -n vep vep \
            --input_file "$INPUT_VCF" --output_file "$OUT_DIR/vep_full.vcf" \
            --format vcf --vcf --offline --fork "$THREADS" \
            --cache --dir_cache "$VEP_CACHE" --assembly GRCh38 \
            --no_stats --force_overwrite $VEP_EXTRA
fi

# SnpEff full
if $RUN_SNPEFF && [[ -f "$SNPEFF_JAR" ]]; then
    STEP1="$OUT_DIR/snpeff_step1.vcf"
    run_timed "SnpEff_full" "$OUT_DIR/snpeff_full.vcf" \
        bash -c "java -Xmx16g -jar '$SNPEFF_JAR' ann -v hg38 -noStats '$INPUT_VCF' \
            > '$STEP1' 2>/dev/null && cp '$STEP1' '$OUT_DIR/snpeff_full.vcf'"
fi

# =============================================================================
# SUMMARY
# =============================================================================
{
echo ""
echo "========================================================"
echo " RESULTS — $(date)"
echo " Input: $VARIANTS variants | Threads: $THREADS"
echo "========================================================"
printf "%-18s %10s %10s %10s %15s\n" "Tool" "Wall(s)" "CPU(s)" "RAM" "Variants/sec"
echo "--------------------------------------------------------"
tail -n +2 "$CSV" | while IFS='|' read -r tool wall cpu ram vps; do
    printf "%-18s %10s %10s %10s %15s\n" "$tool" "$wall" "$cpu" "$ram" "$vps"
done
echo "========================================================"
echo " Full report: $REPORT"
echo " CSV data:    $CSV"
echo "========================================================"
} | tee -a "$REPORT"
