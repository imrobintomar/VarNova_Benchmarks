#!/usr/bin/env bash
# Tier 5 — ClinVar gold-standard HGVS accuracy
#
# Compares VarNova's Gene/HGVS.c/HGVS.p calls against ClinVar's OWN curated
# fields for a Pathogenic/Likely_pathogenic, high-confidence-review-status
# variant set (testdata/clinvar_gold/clinvar_gold_truth.tsv). This is an
# external curated reference, not another annotator's output — see
# ACMG_VALIDATION_FINDINGS.md section 8 for the full discussion and known
# disagreement categories (dup-vs-ins notation, intronic p.. handling, etc).
set -euo pipefail

source "${1:-config.sh}"

VCF="testdata/clinvar_gold/clinvar_gold.vcf"
GOLD="testdata/clinvar_gold/clinvar_gold_truth.tsv"
OUT="results/clinvar_gold"
mkdir -p "$OUT"

if [[ ! -f "$VCF" || ! -f "$GOLD" ]]; then
  echo "Missing $VCF or $GOLD — these ship with this repo under testdata/clinvar_gold/" >&2
  exit 1
fi

echo "== Annotating ClinVar gold VCF with VarNova =="
"$VARNOVA" table -i "$VCF" \
  --gene-db "$DB/hg38_refGene.txt" \
  --filter-db "$DB/hg38_gnomad41_exome.txt,$DB/hg38_dbnsfp47a.txt" \
  --db-columns gnomad411_exome_AF \
  --db-columns SIFT_score,Polyphen2_HDIV_score,CADD_phred,REVEL_score,AlphaMissense_score \
  --label gnomAD_AF --label SIFT --label PolyPhen --label CADD --label REVEL --label AlphaMissense \
  --cache-dir "${VNIDX_CACHE:-}" \
  -o "$OUT/varnova_evidence_only_acmg.tsv" -t "${THREADS:-$(nproc)}"

echo "== Annotating ClinVar gold VCF with ANNOVAR (for the Gene-column cross-check) =="
perl "$ANNOVAR_DIR/table_annovar.pl" "$VCF" "$DB" \
  -buildver hg38 -out "$OUT/annovar_gene" -remove \
  -protocol refGene -operation g -nastring . -vcfinput

echo "== Computing accuracy against ClinVar's curated fields =="
python3 scripts/compute_accuracy.py \
  --varnova "$OUT/varnova_evidence_only_acmg.tsv" \
  --annovar "$OUT/annovar_gene.hg38_multianno.txt" \
  --clinvar-gold "$GOLD" \
  --out-dir "$OUT"

echo "Done — see $OUT/concordance_report.txt"
