#!/usr/bin/env bash
# Tier 3 — Scalability curve
#
# Two independent scaling axes:
#   (a) thread scaling at fixed variant count (1, 2, 4, 8, 16, 32, 48 threads)
#   (b) variant-count scaling at fixed (max) thread count (1M, 10M, 50M, 100M)
# Both use synthetic VCFs (scripts/generate_test_vcf.py) so they're runnable
# without any real patient data, scaled to whatever size your disk/RAM allow.
#
# Uses `/usr/bin/time -f` with a fixed-field format string (not the verbose
# `-v` output) so each run's wall time, CPU time, and peak RSS land in known
# columns — no fragile log-grepping.
set -euo pipefail

source "${1:-config.sh}"

OUT="results"
mkdir -p "$OUT" testdata
TIME_FMT='%e|%U|%M'   # wall(s)|user-cpu(s)|max-RSS(KB)

echo "== (a) Thread scaling, fixed variant count (1M synthetic variants) =="
THREAD_VCF="testdata/scalability_thread_test.vcf"
[[ -f "$THREAD_VCF" ]] || python3 scripts/generate_test_vcf.py --output "$THREAD_VCF" --variants 1000000
n_variants=$(grep -vc '^#' "$THREAD_VCF")

echo "Threads|WallTime(s)|CPUTime(s)|PeakRAM_KB|Variants/sec|Speedup" > "$OUT/scalability_results.csv"
baseline_wall=""
for t in 1 2 4 8 16 32 48; do
  [[ "$t" -gt "$(nproc)" ]] && continue
  read -r wall cpu rss < <(/usr/bin/time -f "$TIME_FMT" "$VARNOVA" gene -i "$THREAD_VCF" \
    -d "$DB/hg38_refGene.txt" -o "$OUT/scalability_t${t}.tsv" -t "$t" 2>&1 \
    | tail -1 | tr '|' ' ')
  vps=$(awk "BEGIN{printf \"%.0f\", $n_variants/$wall}")
  [[ -z "$baseline_wall" ]] && baseline_wall="$wall"
  speedup=$(awk "BEGIN{printf \"%.2f\", $baseline_wall/$wall}")
  echo "$t|$wall|$cpu|$rss|$vps|${speedup}x" >> "$OUT/scalability_results.csv"
done

echo "== (b) Variant-count scaling, fixed thread count =="
echo "Tool|Variants|WallTime(s)|CPUTime(s)|PeakRAM_KB|Variants/sec" > "$OUT/cohort_scale_results.csv"
for n in 1000000 10000000 50000000 100000000; do
  vcf="testdata/cohort_${n}.vcf"
  [[ -f "$vcf" ]] || python3 scripts/generate_test_vcf.py --output "$vcf" --variants "$n"
  read -r wall cpu rss < <(/usr/bin/time -f "$TIME_FMT" "$VARNOVA" gene -i "$vcf" \
    -d "$DB/hg38_refGene.txt" -o "$OUT/cohort_${n}_varnova.tsv" -t "${THREADS:-$(nproc)}" 2>&1 \
    | tail -1 | tr '|' ' ')
  vps=$(awk "BEGIN{printf \"%.0f\", $n/$wall}")
  echo "VarNova_cohort_${n}|${n}|${wall}|${cpu}|${rss}|${vps}" >> "$OUT/cohort_scale_results.csv"
done

echo "Done. NOTE: 100M variants needs ~22GB peak RAM and ~25 minutes wall time"
echo "      on the reference 48-core/62GB machine in the published results —"
echo "      size your hardware accordingly before running the full curve."
