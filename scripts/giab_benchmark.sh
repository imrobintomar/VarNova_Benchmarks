#!/usr/bin/env bash
# Tier 4 — GIAB multi-sample stability
#
# Downloads the NIST Genome in a Bottle (GIAB) v4.2.1 GRCh38 small-variant
# benchmark VCFs for HG001-HG007 and re-runs the VarNova-vs-ANNOVAR
# concordance comparison on each, independently of the single SRR sample
# used in Tier 1/2 — demonstrates the concordance numbers aren't a one-sample
# artifact.
#
# Source: NIST GIAB FTP (ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/)
# NOTE: verify these exact subpaths/filenames against the live FTP listing
# before relying on this script — GIAB periodically revises release
# versions, and the per-sample directory layout (trio grouping) below
# reflects NIST's structure as of this writing, not a guarantee it is
# still current.
set -euo pipefail

source "${1:-config.sh}"

GIAB_ROOT="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release"
OUT="results/giab"
DATA="testdata/giab"
mkdir -p "$OUT" "$DATA"

declare -A SAMPLE_PATH=(
  [HG001]="NA12878_HG001/NISTv4.2.1/GRCh38"
  [HG002]="AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"
  [HG003]="AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38"
  [HG004]="AshkenazimTrio/HG004_NA24143_mother/NISTv4.2.1/GRCh38"
  [HG005]="ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38"
  [HG006]="ChineseTrio/HG006_NA24694_father/NISTv4.2.1/GRCh38"
  [HG007]="ChineseTrio/HG007_NA24695_mother/NISTv4.2.1/GRCh38"
)

RESULTS_CSV="$OUT/giab_multisample_results.csv"
echo "Sample|Variants|VarNova_s|ANNOVAR_s|Func_pct|Gene_pct|ExonicFunc_pct|AAChange_pct" > "$RESULTS_CSV"

for sample in "${!SAMPLE_PATH[@]}"; do
  vcf="$DATA/${sample}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
  if [[ ! -f "$vcf" ]]; then
    url="$GIAB_ROOT/${SAMPLE_PATH[$sample]}/${sample}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    echo "== Downloading $sample =="
    curl -fsSL "$url" -o "$vcf" || { echo "Failed to fetch $url — check the path against the live FTP listing" >&2; continue; }
    curl -fsSL "$url.tbi" -o "$vcf.tbi" || true
  fi

  echo "== Annotating $sample with VarNova =="
  t0=$(date +%s.%N)
  "$VARNOVA" gene -i "$vcf" -d "$DB/hg38_refGene.txt" -o "$OUT/${sample}_varnova.tsv" -t "${THREADS:-$(nproc)}"
  t1=$(date +%s.%N)
  varnova_s=$(echo "$t1 - $t0" | bc)

  echo "== Annotating $sample with ANNOVAR =="
  t0=$(date +%s.%N)
  perl "$ANNOVAR_DIR/table_annovar.pl" "$vcf" "$DB" \
    -buildver hg38 -out "$OUT/${sample}_annovar" -remove \
    -protocol refGene -operation g -nastring . -vcfinput
  t1=$(date +%s.%N)
  annovar_s=$(echo "$t1 - $t0" | bc)

  python3 scripts/compute_accuracy.py \
    --varnova "$OUT/${sample}_varnova.tsv" \
    --annovar "$OUT/${sample}_annovar.hg38_multianno.txt" \
    --out-dir "$OUT/${sample}_accuracy" > "$OUT/${sample}_compare.log"

  n_variants=$(zcat -f "$vcf" | grep -vc '^#' || echo "NA")
  func_pct=$(grep "Func.refGene:" "$OUT/${sample}_accuracy/concordance_report.txt" | grep -oP '[\d.]+(?=% concordant)')
  gene_pct=$(grep "Gene.refGene:" "$OUT/${sample}_accuracy/concordance_report.txt" | grep -oP '[\d.]+(?=% concordant)')
  exonic_pct=$(grep "ExonicFunc.refGene:" "$OUT/${sample}_accuracy/concordance_report.txt" | grep -oP '[\d.]+(?=% concordant)')

  echo "${sample}|${n_variants}|${varnova_s}|${annovar_s}|${func_pct}|${gene_pct}|${exonic_pct}|NA" >> "$RESULTS_CSV"
done

echo "Done — see $RESULTS_CSV"
