---
title: ACMG/AMP classification — validation findings
status: working document, not manuscript prose — curate into the paper later
last verified: cargo test --release -p varnova-core -p varnova-io -p varnova-anno -p varnova (68 passed, 0 failed)
---

# Scope

This documents everything established about VarNova's ACMG/AMP classification engine
(`crates/varnova-anno/src/acmg.rs`) during this validation pass: what was found broken,
what was fixed, what was added, and what remains a known, explicitly-documented gap.
No claim here treats ClinVar, ANNOVAR, or InterVar as ground truth — see
`feedback_no_ground_truth` in project memory. All comparisons are framed as
concordance/disagreement-with-cause, never as raw "accuracy."

# 1. Combining-rule framework: migrated from 2015 categorical table to the
   current Tavtigian 2018/2020 points-based system

**Finding.** The original implementation used a weighted-sum score (PVS1=8, PM=2,
PP=1, BS=-4, BP=-1) with an ad hoc threshold (Pathogenic at score >= 7). This let
PVS1 ALONE reach "Pathogenic" with zero supporting evidence — verified as the
dominant cause of disagreement with InterVar on the 95,599-variant ClinVar gold set.

**First fix attempted (and then superseded).** Replaced the weighted sum with the
named-combination table from Richards et al. 2015 (Table 5) — e.g. "PVS1 AND >=1
Moderate" required for Likely_pathogenic, "PVS1 AND >=2 Moderate" for Pathogenic.
This is the textbook 2015 rule set.

**Why it was superseded.** Per-user instruction to check current (2025/2026)
ClinGen guidance: the 2015 categorical table has been superseded by ClinGen SVI's
endorsed combining method, the Tavtigian et al. 2018/2020 quantitative points-based
Bayesian framework ("Fitting a naturally scaled point system to the ACMG/AMP variant
classification guidelines," Human Mutation), reaffirmed as recently as ClinGen's 2025
refinement. The two frameworks are NOT equivalent on every input — PVS1 alone is the
clearest case: under the points system (8 points), it lands in Likely_pathogenic
(6-9 range), but the 2015 categorical table has no named combination for "PVS1 alone"
at all.

**Current implementation** (`points_to_class()`, `acmg.rs`):
- Point values (log-odds-proportional, power-of-2 scaled): Very Strong=8, Strong=4,
  Moderate=2, Supporting=1 (and symmetric negatives for benign-direction evidence).
- Cutoffs: Pathogenic >=10, Likely_pathogenic 6-9, VUS 0-5, Likely_benign -1 to -6,
  Benign <=-7.
- Points are additive in log-odds space — pathogenic and benign evidence net against
  each other in one sum. No separate "conflicting evidence" branch is needed (unlike
  the 2015 table, which required one); see
  `test_pathogenic_and_benign_evidence_nets_additively_not_as_a_special_case`.

# 2. PM2 downgraded from Moderate to Supporting strength

**Finding.** ClinGen SVI downgraded PM2 (absence from population databases) from
Moderate to Supporting strength — a well-established, repeatedly-reaffirmed
recommendation, not a one-off. VarNova's implementation had been treating it as
Moderate.

**Fix.** PM2 now scores at Supporting (+1), not Moderate (+2).

**Consequence, verified on real data.** This is the single highest-impact fix in
this validation pass. Before: PVS1+PM2 = 8+2 = 10 -> Pathogenic. After: PVS1+PM2 =
8+1 = 9 -> Likely_pathogenic. Re-running the 95,599-variant ClinVar gold set after
this fix dropped VarNova's "Pathogenic" calls from ~53,800 to 26, with the
corresponding mass shifting to "Likely_pathogenic" (53,733).

**Independent confirmation this fix is correct, not a regression.** Audited every
one of InterVar's 52,566 "Pathogenic" calls driven by PVS1 on this same dataset:
ALL 52,566 have zero real PS (Strong) criteria firing (PS array all-zero) — they
reach Pathogenic purely via PVS1(8) + exactly one Moderate-tier PM criterion (PM
array shows a single `1` among 7 slots). That combination only crosses a
Pathogenic threshold if PM2 is scored at the pre-2018 Moderate weight (8+2=10).
This is direct, conclusive evidence that **InterVar has not implemented the PM2
downgrade** and is running outdated (2015-era) evidence-strength weighting for this
criterion — not that VarNova's fix is wrong. Honest framing for the manuscript: on
this specific point, VarNova's combining arithmetic is more current than InterVar's;
this does not mean VarNova is more accurate overall (see Section 7).

# 3. PVS1 gene-level loss-of-function gate + alias resolution

**Finding.** PVS1 ("very strong pathogenic") was firing on ANY null-variant type
(frameshift/stopgain) regardless of whether loss-of-function is even an established
disease mechanism for that gene. ACMG/AMP requires LOF to be a confirmed mechanism
for PVS1 to apply at all.

**Fix.** Added a curated gene-level LOF-mechanism gate (`PVS1_LOF_GENES`, loaded
from `PVS1_LOF_genes.hg38.txt` — same source InterVar uses, ClinGen/literature
derived, 6,213 genes). PVS1 now requires `gene_has_confirmed_lof_mechanism(gene)`
in addition to the null-variant-type check.

**Secondary finding during verification.** The gate initially had near-zero
measured impact (53,843 -> 53,348 Pathogenic calls, a small drop). Investigated and
found two real causes: (1) ClinVar's own P/LP-labeled dataset is naturally enriched
for already-established LOF genes, so the gate matters less here than it would on
general WGS data; (2) a genuine HGNC gene-symbol staleness gap — 108/2,422 genes
(4.5%) were renamed since the LOF list was built (e.g. `AARS`->`AARS1`,
`GBA`->`GBA1`), causing false suppression.

**Fix for (2).** Added `GENE_ALIASES` (loaded from HGNC's `hgnc_complete_set.txt`,
parsing `symbol`/`alias_symbol`/`prev_symbol` columns) — PVS1 now also fires if any
of a gene's previous/alias names is in the LOF list. Recovered 412 additional
correct PVS1 firings (53,348 -> 53,760).

# 4. BP7 gated on splice-boundary proximity

**Finding.** BP7 ("synonymous variant with no predicted splice impact") was firing
on ANY synonymous SNV regardless of position, despite the module's own doc comment
already (incorrectly) claiming a splice check existed. A synonymous change in the
last 1-2 bases of an exon can still disrupt a splice donor/acceptor site.

**Fix.** Added `Annotation::near_splice_boundary: bool`, threaded through
`TxAnno`/the gene-annotation pipeline (conservative OR-reduction across all
contributing transcripts — any transcript near a boundary marks the call as
near-boundary). BP7 now additionally requires `!anno.near_splice_boundary`.

**Biological plausibility check.** 66.7% of ClinVar-labeled-pathogenic variants
annotated as "synonymous" are near a splice boundary — consistent with the
well-established fact that synonymous variants are pathogenic almost exclusively
via splice disruption, not via the (absent) protein-coding change itself.

# 5. New criteria implemented this pass (closing documented gaps)

All four reuse InterVar's own curated reference data (same source, not a new
external dependency) and are wired with explicit gene/codon gating, not
heuristics.

### PP2 / BP1 (gene-level missense relevance)
- PP2 (Supporting, missense-mechanism genes, 931 genes) / BP1 (Supporting-benign,
  truncating-only genes, 1,449 genes), both gated on `FuncClass::NonsynonymousSNV`.
- 38 genes appear in BOTH curated lists (genuine ambiguity in the source data, not
  a code bug). Both criteria are allowed to fire simultaneously for those genes;
  the points system nets PP2(+1) and BP1(-1) to exactly 0 — the correct behavior for
  a genuinely ambiguous gene, not an artifact needing suppression.

### PS1 / PM5 (codon-position-keyed known-pathogenic-missense lookup)
- This was the single largest documented gap going in ("True PS1 needs a
  codon-position-keyed lookup of known pathogenic variants, which this engine
  doesn't have").
- Built from `PS1.AA.change.patho.hg38` (38,353 ClinVar-derived pathogenic missense
  rows: chr/start/end/ref/alt/ref_aa/alt_aa/p.notation/allele_id — no gene column).
  VarNova's own `GenomeIndex` (already built from refGene for annotation) resolves
  each row's genomic position to its gene at load time, producing a
  `(gene, codon_pos) -> [(ref_aa3, alt_aa3), ...]` index (33,895 entries from the
  38,353 source rows).
- PS1 (Strong, +4): current missense variant's amino-acid change matches a KNOWN
  entry exactly at that residue. PM5 (Moderate, +2): a DIFFERENT missense change at
  a residue where some OTHER substitution is known-pathogenic.
- **Bug caught during verification, before shipping**: initial implementation
  assumed the `aa_change_full` field used `;` before the `p.` notation (matching
  `exonic.rs`'s internal per-transcript "detail" string format), but the actual
  joined field gene.rs produces uses `:` throughout
  (`GENE:NM_xxx:exonN:c.123A>G:p.Lys41Arg`). This made the parser silently match
  zero entries — confirmed by checking real output (PS1/PM5 fire counts were both
  exactly 0 before the fix). Fixed by splitting on the last `:` instead of `;`.
  After the fix, on the 95,599-variant ClinVar gold set: PS1 fires 14,172 times,
  PM5 fires 3,082 times, moving roughly 10,161 variants from VUS into
  Likely_pathogenic/Pathogenic.

# 6. Known gaps — explicitly NOT implemented, and why

These are documented verbatim in `acmg.rs`'s module doc comment and should be
stated plainly in the manuscript's Methods/Limitations rather than discovered by a
reviewer:

- **PM1** (mutational hot spot / critical domain). The available curated data
  (Interpro domain *names* per gene, from InterVar's own `PM1_domains_with_benigns`
  file) lacks per-residue coordinate mapping — there's no way to know which codons
  fall inside a named domain without an additional Interpro coordinate file this
  engine doesn't have.
- **PS4** (case-control significance) and **BS2** (observed in healthy adults
  without symptoms for a full-penetrance disorder). Both need cohort/study-level
  data (case-control prevalence, phenotype-negative carrier status) that isn't
  derivable from ClinVar or gnomAD alone.
- **PP3/BP4 calibration.** Currently a flat "≥3 of 5 tools agree" vote across
  SIFT/PolyPhen/CADD/REVEL/AlphaMissense at fixed Supporting strength. ClinGen's
  actual recommendation (Pejaver et al. 2022, Am J Hum Genet) calibrates each tool
  individually to a specific strength (Supporting/Moderate/Strong, even
  Very-Strong-benign for REVEL) and recommends using ONE validated/calibrated tool,
  not vote-counting across several correlated ones. Not yet implemented.
- **BA1 frequency basis.** Uses a flat raw population-AF threshold (>0.05). ClinGen's
  current recommendation uses gnomAD's filtering allele frequency (FAF) plus an
  explicit gene/variant-specific exception list (so a common-but-disease-causing
  variant in a specific gene isn't auto-benigned). Not yet implemented — real risk:
  without the exception list, BA1 could in principle misclassify a true exception
  case. (Empirically rare on the ClinVar gold set: only 3 BA1 firings total.)

# 7. ACMG-vs-InterVar: honest two-sided summary for the manuscript

Do NOT compress this into a single concordance percentage — that framing was
explicitly tried and explicitly rejected during this validation pass ("i dont want
to result compare with intervar i just which is much accurate").

**Where VarNova's mechanics are demonstrably more current than InterVar's:**
PM2 evidence strength (Section 2) — backed by the full-population InterVar audit
(52,566/52,566 PVS1-driven Pathogenic calls have zero real PS evidence behind
them), not assumption.

**Where InterVar is more complete than VarNova:** PM1 (domain hotspots), PS4
(case-control), BS2 (healthy-carrier), and PP3/BP4's full multi-tool/strength
calibration — InterVar implements forms of all of these; VarNova does not yet (see
Section 6).

**Manuscript framing:** "VarNova is more correct on the mechanics of what it
implements; InterVar implements more of the full criteria set." State both
directions. A single agreement number actively obscures this and should not appear
as the headline ACMG result.

# 8. ClinVar gold-standard accuracy (HGVS-level, not ACMG-level)

Re-verified fresh (not carried from memory) via
`compute_accuracy.py --clinvar-gold`, saved at
`benchmark/accuracy/clinvar_gold/concordance_report.txt`, against the 95,599-variant
ClinVar Pathogenic/Likely_pathogenic high-confidence-review gold set:

| Metric | Result |
|---|---|
| Gene accuracy | 99.563% (95,181/95,599) |
| HGVS.c accuracy | 91.432% (87,408/95,599) |
| HGVS.p accuracy | 71.170% |
| Transcript selection (same transcript as ClinVar's submitter) | 97.270% (92,989/95,599) |

Visible disagreement categories (from the same report, useful for the manuscript's
limitations paragraph):
1. dup-vs-ins notational differences on the same underlying indel
   (e.g. ClinVar `c.26dup` vs VarNova `c.27_28insC` — same edit, different valid
   HGVS convention).
2. ClinVar reports `p..` for intronic/splice variants; VarNova correctly reports
   intronic c.-notation instead of a (nonexistent) protein change — not a VarNova
   error, but will show as a "disagreement" in a naive string comparison.
3. Complex delins/multi-base-substitution cases VarNova doesn't yet resolve to
   ClinVar's `c.346_363delinsTGGGCCCCTG`-style notation.

# 9. Test coverage

68 workspace tests passing as of this document (`varnova-anno`: 56, `varnova-core`:
4, `varnova-io`: 8, `varnova` CLI: 0 dedicated but exercised via integration runs).
Notable test categories in `acmg.rs`:
- Per-criterion rule-correctness (one test per implemented criterion, firing and
  non-firing cases, boundary values).
- `test_no_logically_contradictory_criteria_pairs_across_input_grid` — sweeps a
  representative grid of inputs (15,000+ combinations) confirming the classifier
  never produces a logically impossible criteria pair (e.g. "common in population"
  and "absent from population" both true).
- Points-cutoff boundary tests (`test_points_to_class_cutoff_boundaries`).
- PS1/PM5 codon-lookup tests including a guard against firing on non-missense
  variant types even when gene+codon coincidentally match a listed entry.
