#!/usr/bin/env bash
# =============================================================================
# VarNova Fair Benchmark
#
# TIER 1 — Gene annotation only (all 4 tools)
#   VarNova | ANNOVAR | VEP | SnpEff
#
# TIER 2 — Full pipeline, identical databases (all 4 tools)
#   Databases: refGene + gnomAD 4.1 exome + ClinVar 2024 + dbNSFP 4.7a
#   VarNova | ANNOVAR | VEP + plugins | SnpEff + SnpSift
#
# Metrics per tool: wall time, CPU time, peak RAM, variants/sec, output count
# =============================================================================

set -uo pipefail
export PATH="$HOME/bin:$PATH"

# ── Paths ─────────────────────────────────────────────────────────────────────

INPUT_VCF="/media/drprabudh/m3/varnova/testVCF/SRR36790057.vcf"
DB="/media/drprabudh/m1/annovar/hg38_humandb"
ANNOVAR_DIR="/media/drprabudh/m1/annovar"
SNPEFF_JAR="/media/drprabudh/m1/annovar/snpEff/snpEff.jar"
VNIDX_CACHE="/home/drprabudh/varnova_db"
SNPSIFT_JAR="/media/drprabudh/m1/annovar/snpEff/SnpSift.jar"
DBNSFP_BGZ="/media/drprabudh/m1/annovar/dbNSFP4.7a/dbNSFP4.7a_combined.gz"
GNOMAD_BGZ="/media/drprabudh/m1/annovar/gnomad_vcf/gnomad.exome.v4.1.vcf.bgz"
CLINVAR_BGZ="/media/drprabudh/m1/vep_cache/clinvar.vcf.bgz"
VEP_CACHE="/media/drprabudh/m1/vep_cache"
VEP_PLUGIN_DIR="/home/drprabudh/miniconda3/envs/vep/share/ensembl-vep-115.2-1"
OUT_DIR="/media/drprabudh/m3/varnova/benchmark/results"
ACCURACY_DIR="/media/drprabudh/m3/varnova/benchmark/accuracy"
THREADS=48
VARIANTS=89052

mkdir -p "$OUT_DIR" "$ACCURACY_DIR"
REPORT="$OUT_DIR/benchmark_report.txt"
CSV="$OUT_DIR/benchmark_results.csv"

# ── Helpers ───────────────────────────────────────────────────────────────────

log()  { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$REPORT"; }
sep()  { echo "─────────────────────────────────────────────────" | tee -a "$REPORT"; }

run_timed() {
    local label="$1"
    local outfile="$2"
    shift 2
    local timefile; timefile=$(mktemp)

    log "Running: $label"
    log "Command: $*"

    /usr/bin/time -v -o "$timefile" "$@" >> "$outfile" 2>&1 || true

    local wall_sec cpu_sec ram_kb
    wall_sec=$(grep "Elapsed (wall clock)" "$timefile" | grep -oP '\d+:\d+\.\d+$' | \
               awk -F: '{printf "%.1f", $1*60+$2}' 2>/dev/null || echo "N/A")
    cpu_sec=$(grep "User time"            "$timefile" | grep -oP '[\d.]+$' || echo "N/A")
    ram_kb=$( grep "Maximum resident"     "$timefile" | grep -oP '\d+$'    || echo "N/A")

    local ram_mb="N/A"
    [[ "$ram_kb" != "N/A" ]] && ram_mb=$(echo "$ram_kb / 1024" | bc)

    local vps="N/A"
    if [[ "$wall_sec" != "N/A" && "$wall_sec" != "0" ]]; then
        vps=$(echo "scale=0; $VARIANTS / $wall_sec" | bc 2>/dev/null || echo "N/A")
    fi

    rm -f "$timefile"
    echo "$label|$wall_sec|$cpu_sec|${ram_mb} MB|$vps" >> "$CSV"
    log "  Wall time : ${wall_sec}s | CPU: ${cpu_sec}s | RAM: ${ram_mb} MB | Speed: ${vps} v/s"
    sep
}

count_variants() {
    local f="$1"
    if [[ -f "$f" ]]; then
        # TSV: count non-header lines
        if [[ "$f" == *.tsv ]] || [[ "$f" == *.txt ]]; then
            tail -n +2 "$f" | wc -l
        # VCF: count non-# lines
        else
            grep -v "^#" "$f" | wc -l
        fi
    else
        echo "FILE_MISSING"
    fi
}

# ── Header ────────────────────────────────────────────────────────────────────

{
echo "============================================================"
echo " VarNova Fair Benchmark — $(date)"
echo " Machine : $(nproc) cores | $(free -h | awk '/Mem/{print $2}') RAM"
echo " Input   : $INPUT_VCF ($VARIANTS variants)"
echo " Threads : $THREADS"
echo ""
echo " TIER 1 — Gene annotation only  (all 4 tools)"
echo " TIER 2 — Full pipeline         (all 4 tools, same databases)"
echo "          Gene + gnomAD 4.1 exome + ClinVar 2024 + dbNSFP 4.7a"
echo "============================================================"
} | tee "$REPORT"

echo "Tool|Tier|WallTime(s)|CPUTime(s)|PeakRAM|Variants/sec|OutputCount" > "$CSV"

sep

# =============================================================================
# TIER 1 — GENE ANNOTATION ONLY
# =============================================================================

log "══════════════════════════════════════════════"
log " TIER 1: GENE ANNOTATION ONLY"
log "══════════════════════════════════════════════"
sep

# ── VarNova gene ──────────────────────────────────────────────────────────────
log "T1-1: VarNova — Gene only"
run_timed "VarNova_gene" "$OUT_DIR/varnova_gene.tsv" \
    varnova gene \
        -i "$INPUT_VCF" \
        -d "$DB/hg38_refGene.txt" \
        -o "$OUT_DIR/varnova_gene.tsv" \
        -t $THREADS

# ── ANNOVAR gene ──────────────────────────────────────────────────────────────
log "T1-2: ANNOVAR — Gene only"
perl "$ANNOVAR_DIR/convert2annovar.pl" \
    -format vcf4 "$INPUT_VCF" \
    -outfile "$OUT_DIR/annovar_input.avinput" \
    -includeinfo 2>/dev/null || true

run_timed "ANNOVAR_gene" "$OUT_DIR/annovar_gene.out" \
    perl "$ANNOVAR_DIR/annotate_variation.pl" \
        -geneanno -buildver hg38 -dbtype refGene \
        -outfile "$OUT_DIR/annovar_gene" \
        "$OUT_DIR/annovar_input.avinput" "$DB"

# ── VEP gene only ─────────────────────────────────────────────────────────────
log "T1-3: VEP 115.2 — Gene only"
rm -f "$OUT_DIR/vep_gene_only.vcf"
run_timed "VEP_gene" "$OUT_DIR/vep_gene_only.vcf" \
    conda run -n vep vep \
        --input_file  "$INPUT_VCF" \
        --output_file "$OUT_DIR/vep_gene_only.vcf" \
        --format vcf --vcf --offline \
        --fork $THREADS \
        --cache --dir_cache "$VEP_CACHE" \
        --assembly GRCh38 \
        --no_stats --force_overwrite

# ── SnpEff gene only ──────────────────────────────────────────────────────────
log "T1-4: SnpEff 5.2a — Gene only"
run_timed "SnpEff_gene" "$OUT_DIR/snpeff_gene.vcf" \
    java -Xmx16g -jar "$SNPEFF_JAR" ann -v hg38 -noStats "$INPUT_VCF"

# =============================================================================
# TIER 2 — FULL PIPELINE (identical databases: Gene + gnomAD + ClinVar + dbNSFP)
# =============================================================================

log "══════════════════════════════════════════════"
log " TIER 2: FULL PIPELINE"
log " Databases: refGene + gnomAD 4.1 exome + ClinVar 2024 + avsnp151"
log "══════════════════════════════════════════════"
sep

# ── VarNova full ──────────────────────────────────────────────────────────────
log "T2-1: VarNova"
run_timed "VarNova_full" "$OUT_DIR/varnova_full.log" \
    varnova table \
        -i "$INPUT_VCF" \
        --gene-db   "$DB/hg38_refGene.txt" \
        --filter-db "$DB/hg38_gnomad41_exome.txt" \
        --filter-db "$DB/hg38_clinvar_20240416.txt" \
        --filter-db "$DB/hg38_avsnp151.txt" \
        --label gnomAD --label ClinVar --label rsID \
        --cache-dir "$VNIDX_CACHE" \
        --out-dir "$OUT_DIR" \
        -t $THREADS -v

VARNOVA_FULL_COUNT=$(count_variants "$OUT_DIR/SRR36790057.annotated.tsv")
log "  Output variants: $VARNOVA_FULL_COUNT"

# ── ANNOVAR full ──────────────────────────────────────────────────────────────
log "T2-2: ANNOVAR"
run_timed "ANNOVAR_full" "$OUT_DIR/annovar_full.out" \
    perl "$ANNOVAR_DIR/table_annovar.pl" \
        "$INPUT_VCF" "$DB" \
        -buildver hg38 \
        -out "$OUT_DIR/annovar_full" \
        -remove \
        -protocol refGene,gnomad41_exome,clinvar_20240416,avsnp151 \
        -operation g,f,f,f \
        -vcfinput -nastring .

ANNOVAR_FULL_OUT=$(ls "$OUT_DIR"/annovar_full.hg38_multianno.txt 2>/dev/null | head -1)
ANNOVAR_FULL_COUNT=$(count_variants "$ANNOVAR_FULL_OUT")
log "  Output variants: $ANNOVAR_FULL_COUNT"

# ── VEP full ──────────────────────────────────────────────────────────────────
log "T2-3: VEP 115.2"
rm -f "$OUT_DIR/vep_full.vcf"
run_timed "VEP_full" "$OUT_DIR/vep_full.vcf" \
    conda run -n vep vep \
        --input_file   "$INPUT_VCF" \
        --output_file  "$OUT_DIR/vep_full.vcf" \
        --format vcf --vcf --offline \
        --fork         $THREADS \
        --cache --dir_cache "$VEP_CACHE" \
        --dir_plugins  "$VEP_PLUGIN_DIR" \
        --assembly     GRCh38 \
        --no_stats --force_overwrite \
        --plugin dbNSFP,"$DBNSFP_BGZ",SIFT_score,Polyphen2_HDIV_score,CADD_phred,REVEL_score \
        --custom "$GNOMAD_BGZ",gnomAD_exome,vcf,exact,0,AF,AF_afr,AF_amr,AF_eas,AF_nfe,AF_sas \
        --custom "$CLINVAR_BGZ",ClinVar,vcf,exact,0,CLNSIG,CLNDN

VEP_FULL_COUNT=$(count_variants "$OUT_DIR/vep_full.vcf")
log "  Output variants: $VEP_FULL_COUNT"

# ── SnpEff full (SnpEff + SnpSift) ───────────────────────────────────────────
log "T2-4: SnpEff 5.2a + SnpSift"
SNPEFF_STEP1="$OUT_DIR/snpeff_step1.vcf"
SNPEFF_STEP2="$OUT_DIR/snpeff_step2.vcf"
SNPEFF_STEP3="$OUT_DIR/snpeff_step3.vcf"

run_timed "SnpEff_full" "$OUT_DIR/snpeff_full.vcf" \
    bash -c "
        java -Xmx16g -jar '$SNPEFF_JAR' ann -v hg38 -noStats '$INPUT_VCF' \
            > '$SNPEFF_STEP1' 2>/dev/null && \
        java -Xmx8g -jar '$SNPSIFT_JAR' DbNsfp \
            -db '$DBNSFP_BGZ' \
            -f SIFT_score,Polyphen2_HDIV_score,CADD_phred,REVEL_score \
            '$SNPEFF_STEP1' > '$SNPEFF_STEP2' 2>/dev/null && \
        java -Xmx8g -jar '$SNPSIFT_JAR' annotate \
            '$GNOMAD_BGZ' '$SNPEFF_STEP2' > '$SNPEFF_STEP3' 2>/dev/null && \
        java -Xmx8g -jar '$SNPSIFT_JAR' annotate \
            '$CLINVAR_BGZ' '$SNPEFF_STEP3'
    "

rm -f "$SNPEFF_STEP1" "$SNPEFF_STEP2" "$SNPEFF_STEP3"
SNPEFF_FULL_COUNT=$(count_variants "$OUT_DIR/snpeff_full.vcf")
log "  Output variants: $SNPEFF_FULL_COUNT"

# =============================================================================
# SUMMARY TABLE
# =============================================================================

{
echo ""
echo "============================================================"
echo " TIER 1 — GENE ANNOTATION ONLY"
printf "%-15s %10s %10s %10s %15s\n" "Tool" "Wall(s)" "CPU(s)" "RAM" "Variants/sec"
echo "------------------------------------------------------------"
grep "_gene|" "$CSV" | while IFS='|' read -r tool tier wall cpu ram vps count; do
    printf "%-15s %10s %10s %10s %15s\n" "$tool" "$wall" "$cpu" "$ram" "$vps"
done
echo ""
echo " TIER 2 (Gene + gnomAD + ClinVar + dbNSFP)"
printf "%-15s %10s %10s %10s %15s %12s\n" "Tool" "Wall(s)" "CPU(s)" "RAM" "Variants/sec" "OutputCount"
echo "------------------------------------------------------------"
for tool_label in VarNova_full ANNOVAR_full VEP_full SnpEff_full; do
    line=$(grep "^${tool_label}|" "$CSV" | tail -1)
    IFS='|' read -r tool wall cpu ram vps <<< "$line"
    case "$tool_label" in
        VarNova_full) cnt=$VARNOVA_FULL_COUNT ;;
        ANNOVAR_full) cnt=$ANNOVAR_FULL_COUNT ;;
        VEP_full)     cnt=$VEP_FULL_COUNT ;;
        SnpEff_full)  cnt=$SNPEFF_FULL_COUNT ;;
        *)            cnt="N/A" ;;
    esac
    printf "%-15s %10s %10s %10s %15s %12s\n" "$tool" "$wall" "$cpu" "$ram" "$vps" "$cnt"
done
echo "============================================================"
} | tee -a "$REPORT"

# =============================================================================
# PHASE 2 — ACCURACY ANALYSIS
# =============================================================================

log ""
log "============================================================"
log " PHASE 2: ACCURACY ANALYSIS"
log "============================================================"

log "Step 1: Building gold standard from ClinVar..."
python3 "$ACCURACY_DIR/build_goldstandard.py" 2>&1 | tee -a "$REPORT"

log "Step 2: Computing accuracy metrics and generating plots..."
python3 "$ACCURACY_DIR/compute_accuracy.py" 2>&1 | tee -a "$REPORT"

log ""
log "============================================================"
log " BENCHMARK COMPLETE"
log " Report : $REPORT"
log " CSV    : $CSV"
log " Plots  : $ACCURACY_DIR/plots/"
log "============================================================"
