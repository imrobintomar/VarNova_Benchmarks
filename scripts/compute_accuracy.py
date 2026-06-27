#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
import sys
from collections import Counter, defaultdict

CONCORDANCE_THRESHOLD = 99.0

# Ensembl VEP consequence terms -> ANNOVAR ExonicFunc.refGene vocabulary.
# Only the terms relevant to coding-exon classification are mapped; anything
# else (UTR, intronic, regulatory, etc.) is compared at the Func.refGene
# (broad category) level instead, since VEP's broad-category vocabulary
# differs from ANNOVAR's by design and isn't a 1:1 mapping.
VEP_TO_EXONIC_FUNC = {
    "synonymous_variant": "synonymous SNV",
    "missense_variant": "nonsynonymous SNV",
    "stop_gained": "stopgain",
    "stop_lost": "stoploss",
    "start_lost": "startloss",
    "frameshift_variant": None,  # ambiguous ins/del without further parsing
    "inframe_insertion": "nonframeshift insertion",
    "inframe_deletion": "nonframeshift deletion",
}

VEP_SEVERITY_ORDER = [
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "missense_variant", "inframe_insertion", "inframe_deletion",
    "splice_donor_variant", "splice_acceptor_variant",
    "synonymous_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
    "intron_variant", "upstream_gene_variant", "downstream_gene_variant",
    "intergenic_variant",
]


def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def normalize_chrom(c):
    c = c[3:] if c.lower().startswith("chr") else c
    return "M" if c.upper() == "MT" else c


def load_tsv(path):
    """Load a Chr/Start/End/Ref/Alt + annotation TSV into a list of dicts."""
    rows = []
    with open_maybe_gz(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if "Chr" in row and row["Chr"]:
                row["Chr"] = normalize_chrom(row["Chr"])
            rows.append(row)
    return rows


def normalize_cdna_snv(cdna):
    """ANNOVAR and VarNova use DIFFERENT c. notation for the same SNV:
    ANNOVAR:  c.{ref}{pos}{alt}   e.g. c.A1470G
    VarNova:  c.{pos}{ref}>{alt}  e.g. c.1470A>G  (standard HGVS)
    Comparing these as raw strings — as an earlier version of this script
    did — makes every single matching SNV look like a mismatch (verified:
    dropped a real benchmark run's AAChange concordance to 1.7%, when the
    underlying calls actually agreed). Normalize both to (pos, ref, alt)
    before comparing. Returns None for anything that isn't a plain SNV
    (indels keep their raw string — both tools' indel notation already
    starts with the position, so raw comparison is less broken there, just
    not bulletproof for differing dup/ins/del conventions)."""
    m = re.match(r"c\.(\d+)([ACGT])>([ACGT])$", cdna)
    if m:
        return (m.group(1), m.group(2), m.group(3))
    m = re.match(r"c\.([ACGT])(\d+)([ACGT])$", cdna)
    if m:
        return (m.group(2), m.group(1), m.group(3))
    return None


def normalize_dup_suffix(cdna):
    """HGVS dup notation optionally includes the duplicated base(s) for a
    single-position dup (c.667dupC) — ClinVar's Name field always omits it
    (c.667dup), VarNova includes it. Both describe the identical change;
    strip the optional trailing base letters so they compare equal instead
    of registering as a spurious mismatch (verified: 3,195 of the ClinVar
    gold-standard HGVS.c mismatches were this exact cosmetic difference,
    not a real positional or calling disagreement)."""
    return re.sub(r"dup[ACGT]+$", "dup", cdna)


def parse_aachange_set(field):
    """AAChange.refGene: GENE:NM_xxx:exonN:c.XXX:p.YYY, comma-joined per
    transcript. Returns a set of (transcript_without_version, normalized_c.)
    pairs — cDNA-level identity is what actually matters; comparing the full
    string (or including p.) would spuriously fail on cosmetic differences
    (exon numbering conventions, protein 1-letter vs 3-letter style,
    transcript version suffixes, and ANNOVAR-vs-HGVS c. notation order) that
    don't reflect a real disagreement."""
    out = set()
    if not field or field == ".":
        return out
    for entry in field.split(","):
        parts = entry.split(":")
        if len(parts) < 4:
            continue
        tx = parts[1].split(".")[0]  # drop version suffix (NM_007294.4 -> NM_007294)
        cdna = parts[3]
        if not cdna.startswith("c."):
            continue
        snv = normalize_cdna_snv(cdna)
        out.add((tx, snv) if snv else (tx, normalize_dup_suffix(cdna)))
    return out


_AA_1TO3 = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'Q': 'Gln',
    'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys',
    'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp',
    'Y': 'Tyr', 'V': 'Val', 'X': 'Ter', '*': 'Ter',
}


def normalize_hgvsp(p):
    """Normalize HGVS.p across tool conventions before comparing: ANNOVAR
    uses 1-letter AA codes + 'X' for stop + 'fs*NN' frameshift offset (e.g.
    p.S194Gfs*60); VarNova/ClinVar use 3-letter codes + 'Ter' + bare 'fs'
    (e.g. p.Ser194fs). Returns a comparable (pos, aa_from, aa_to_or_'fs')
    tuple for simple substitution/stop/frameshift changes — the vast
    majority of cases — or None for anything more complex (range/delins),
    which falls back to raw-string comparison instead."""
    if not p or not p.startswith("p."):
        return None
    body = p[2:]
    m = re.match(r'^([A-Za-z]{1,3})(\d+)[A-Za-z*]*fs', body)
    if m:
        aa, pos = m.group(1), m.group(2)
        aa3 = _AA_1TO3.get(aa.upper(), aa) if len(aa) == 1 else aa
        return (pos, aa3, "fs")
    m = re.match(r'^([A-Za-z]{1,3})(\d+)([A-Za-z*]{1,3})$', body)
    if m:
        aa1, pos, aa2 = m.group(1), m.group(2), m.group(3)
        aa1_3 = _AA_1TO3.get(aa1.upper(), aa1) if len(aa1) == 1 else aa1
        if aa2 == "*" or aa2.upper() == "X":
            aa2_3 = "Ter"
        else:
            aa2_3 = _AA_1TO3.get(aa2.upper(), aa2) if len(aa2) == 1 else aa2
        return (pos, aa1_3, aa2_3)
    return None


def parse_hgvsp_set(field):
    """Same transcript-level pairing as parse_aachange_set, but keeps the
    protein-level p. notation instead of the cDNA change, normalized via
    normalize_hgvsp() so 1-letter/3-letter and X/Ter/fs*N convention
    differences across tools don't register as spurious mismatches."""
    out = set()
    if not field or field == ".":
        return out
    for entry in field.split(","):
        parts = entry.split(":")
        if len(parts) < 5:
            continue
        tx = parts[1].split(".")[0]
        p = parts[4]
        if p.startswith("p."):
            norm = normalize_hgvsp(p)
            out.add((tx, norm if norm else p))
    return out


def clinvar_gold_join_key(chrom, pos, ref, alt):
    """ClinVar's variant_summary.txt gives raw VCF-anchor-style Ref/Alt
    (e.g. G->GGGGCC), but VarNova's own output strips the anchor base for
    simple indels (Ref=-, Alt=GGGCC) the same way ANNOVAR's avinput
    conversion does — verified against real output for both an insertion
    (chr1:1040717 G>GGGGCC -> pos unchanged, Ref=-, Alt=GGGCC) and a
    deletion (chr1:1041678 CGCT...>C -> pos+1, Ref=GCT..., Alt=-). Complex
    delins/substitutions where ref/alt don't share a simple single-base
    anchor pass through unchanged — verified VarNova does the same."""
    pos = int(pos)
    if len(ref) == 1 and len(alt) > 1 and alt.startswith(ref):
        return (chrom, pos, "-", alt[1:])
    if len(alt) == 1 and len(ref) > 1 and ref.startswith(alt):
        return (chrom, pos + 1, ref[1:], "-")
    return (chrom, pos, ref, alt)


def load_clinvar_gold(path):
    """Load the ClinVar gold-standard truth table built from NCBI's
    variant_summary.txt (GRCh38, Pathogenic/Likely_pathogenic, restricted to
    practice-guideline/expert-panel/multiple-submitters-no-conflicts review
    status) — keyed using the same anchor-stripped convention VarNova's own
    output uses, so lookups join directly without re-deriving it per query."""
    gold = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            chrom = normalize_chrom(row["Chrom"])
            key = clinvar_gold_join_key(chrom, row["Pos"], row["Ref"], row["Alt"])
            gold[key] = row
    return gold


def compare_clinvar_gold(varnova_rows, gold_path, report):
    """Section 3 — accuracy against ClinVar's OWN curated Gene/HGVS.c/
    HGVS.p/Transcript fields, not another annotation tool's output. This is
    the strongest accuracy claim available here: an external curated
    reference rather than a tool-vs-tool concordance number (see
    feedback_no_ground_truth — still not infallible: ClinVar submitters can
    choose a different valid transcript than VarNova does, which is why
    transcript selection is reported separately from cDNA/protein accuracy
    rather than folded into one pass/fail number)."""
    report.append("")
    report.append("=" * 70)
    report.append("Section 3 — Accuracy against ClinVar gold-standard "
                   "(Pathogenic/Likely_pathogenic, high-confidence review status)")
    report.append("=" * 70)

    gold = load_clinvar_gold(gold_path)
    report.append(f"Gold-standard variants: {len(gold)}")

    n = 0
    gene_ok = hgvsc_ok = hgvsp_ok = hgvsp_n = transcript_ok = 0
    gene_ex, hgvsc_ex, hgvsp_ex = [], [], []

    for v in varnova_rows:
        try:
            key = clinvar_gold_join_key(v.get("Chr", ""), v.get("Start", "0"),
                                         v.get("Ref", ""), v.get("Alt", ""))
        except (ValueError, TypeError):
            continue
        g = gold.get(key)
        if g is None:
            continue
        n += 1

        locus = f"{v.get('Chr')}:{v.get('Start')} {v.get('Ref')}>{v.get('Alt')}"
        varnova_aachange = v.get("AAChange.refGene", ".")

        gold_gene = g["Gene"]
        if gold_gene in parse_gene_set(v.get("Gene.refGene", ".")):
            gene_ok += 1
        elif len(gene_ex) < 20:
            gene_ex.append(f"  {locus}  clinvar='{gold_gene}' varnova='{v.get('Gene.refGene', '.')}'")

        gold_tx = g["Transcript"].split(".")[0]
        gold_snv = normalize_cdna_snv(g["HGVSc"])
        gold_cdna_key = (gold_tx, gold_snv if gold_snv else normalize_dup_suffix(g["HGVSc"]))
        if gold_cdna_key in parse_aachange_set(varnova_aachange):
            hgvsc_ok += 1
        elif len(hgvsc_ex) < 20:
            hgvsc_ex.append(f"  {locus}  clinvar='{gold_tx}:{g['HGVSc']}' varnova='{varnova_aachange}'")

        if g["HGVSp"] != "p.":
            hgvsp_n += 1
            gold_p_norm = normalize_hgvsp(g["HGVSp"])
            gold_p_key = (gold_tx, gold_p_norm if gold_p_norm else g["HGVSp"])
            if gold_p_key in parse_hgvsp_set(varnova_aachange):
                hgvsp_ok += 1
            elif len(hgvsp_ex) < 20:
                hgvsp_ex.append(f"  {locus}  clinvar='{gold_tx}:{g['HGVSp']}' varnova='{varnova_aachange}'")

        # Transcript selection: did VarNova's call use the SAME transcript
        # ClinVar's submitter referenced — independent of whether the change
        # itself matched (a different valid isoform isn't necessarily wrong).
        varnova_txs = {tx for tx, _ in parse_aachange_set(varnova_aachange)}
        if gold_tx in varnova_txs:
            transcript_ok += 1

    report.append(f"Rows matched by locus (Chr:Pos:Ref:Alt, anchor-normalized): {n}")
    if n:
        report.append(f"  Gene accuracy:        {100.0*gene_ok/n:.3f}% ({gene_ok}/{n})")
        report.append(f"  HGVS.c accuracy:      {100.0*hgvsc_ok/n:.3f}% ({hgvsc_ok}/{n})")
        if hgvsp_n:
            report.append(f"  HGVS.p accuracy:      {100.0*hgvsp_ok/hgvsp_n:.3f}% ({hgvsp_ok}/{hgvsp_n})")
        report.append(f"  Transcript selection: {100.0*transcript_ok/n:.3f}% ({transcript_ok}/{n}) "
                       f"— uses the SAME transcript ClinVar's submitter referenced")
    if gene_ex:
        report.append("-- Gene disagreements (up to 20) --")
        report.extend(gene_ex)
    if hgvsc_ex:
        report.append("-- HGVS.c disagreements (up to 20) --")
        report.extend(hgvsc_ex)
    if hgvsp_ex:
        report.append("-- HGVS.p disagreements (up to 20) --")
        report.extend(hgvsp_ex)


def parse_gene_set(field):
    """Gene.refGene: ANNOVAR joins genuinely distinct overlapping/flanking
    genes with ';' (e.g. "MIR200A;MIR200B"), VarNova joins them with ','
    (e.g. "MIR200A,MIR200B") — same convention ANNOVAR itself also uses for
    multi-gene readthrough names within one entry (e.g. "PMF1,PMF1-BGLAP").
    Comparing the raw string spuriously fails on this separator difference
    alone even when the gene set is identical, so split on both delimiters
    and compare as sets."""
    if not field or field in (".", "NONE"):
        return set()
    return {g for g in re.split(r"[,;]", field) if g}


def compare_aachange(varnova_rows, annovar_rows, report):
    report.append("")
    n = min(len(varnova_rows), len(annovar_rows))
    exact, overlap, both_empty, compared = 0, 0, 0, 0

    for i in range(n):
        v, a = varnova_rows[i], annovar_rows[i]
        if (v.get("Chr") != a.get("Chr") or v.get("Start") != a.get("Start")
                or v.get("Ref") != a.get("Ref") or v.get("Alt") != a.get("Alt")):
            continue
        v_field = v.get("AAChange.refGene", ".")
        a_field = a.get("AAChange.refGene", ".")
        v_set = parse_aachange_set(v_field)
        a_set = parse_aachange_set(a_field)
        if not v_set and not a_set:
            both_empty += 1
            continue  # neither tool called a coding change here — not informative
        compared += 1
        if v_field == a_field:
            exact += 1
        if v_set & a_set:
            overlap += 1

    exact_pct = 100.0 * exact / compared if compared else 0.0
    overlap_pct = 100.0 * overlap / compared if compared else 0.0
    report.append(
        f"AAChange.refGene (cDNA-level, {compared} rows where either tool called "
        f"a coding change; {both_empty} rows where neither did, excluded):"
    )
    report.append(f"  exact full-string match:        {exact_pct:.3f}% ({exact}/{compared})")
    report.append(
        f"  at least one shared (tx, c.) pair: {overlap_pct:.3f}% ({overlap}/{compared}) "
        f"— the more meaningful number; exact-string match is needlessly strict "
        f"(transcript ordering, exon-numbering convention differ even when the "
        f"underlying call is identical)"
    )
    return overlap_pct


def precision_recall_f1(y_true, y_pred, labels):
    """Macro-averaged precision/recall/F1 over the given label set."""
    tp = Counter()
    fp = Counter()
    fn = Counter()
    for t, p in zip(y_true, y_pred):
        if t == p:
            tp[t] += 1
        else:
            fp[p] += 1
            fn[t] += 1

    precisions, recalls, f1s = [], [], []
    for label in labels:
        tpc, fpc, fnc = tp[label], fp[label], fn[label]
        prec = tpc / (tpc + fpc) if (tpc + fpc) else 0.0
        rec = tpc / (tpc + fnc) if (tpc + fnc) else 0.0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
        precisions.append(prec)
        recalls.append(rec)
        f1s.append(f1)
    n = len(labels) or 1
    return sum(precisions) / n, sum(recalls) / n, sum(f1s) / n


def compare_annovar(varnova_rows, annovar_rows, report):
    report.append("=" * 70)
    report.append("VarNova vs ANNOVAR — direct row-aligned comparison")
    report.append("=" * 70)

    n = min(len(varnova_rows), len(annovar_rows))
    if len(varnova_rows) != len(annovar_rows):
        report.append(
            f"WARNING: row count mismatch — VarNova={len(varnova_rows)} "
            f"ANNOVAR={len(annovar_rows)}. Comparing first {n} rows only; "
            f"this likely means the two runs used different inputs/expansion."
        )

    mismatches = defaultdict(list)
    columns = ["Func.refGene", "Gene.refGene", "ExonicFunc.refGene"]
    matches = {c: 0 for c in columns}

    exonic_true, exonic_pred = [], []

    for i in range(n):
        v, a = varnova_rows[i], annovar_rows[i]
        # sanity: confirm the rows actually describe the same variant
        same_locus = (
            v.get("Chr") == a.get("Chr")
            and v.get("Start") == a.get("Start")
            and v.get("Ref") == a.get("Ref")
            and v.get("Alt") == a.get("Alt")
        )
        if not same_locus:
            mismatches["LOCUS_MISALIGNED"].append(i)
            continue

        for c in columns:
            vv, av = v.get(c, "."), a.get(c, ".")
            if c == "Gene.refGene":
                is_match = parse_gene_set(vv) == parse_gene_set(av)
            else:
                is_match = vv == av
            if is_match:
                matches[c] += 1
            elif len(mismatches[c]) < 50:
                mismatches[c].append(
                    f"  row {i}: {v.get('Chr')}:{v.get('Start')} "
                    f"{v.get('Ref')}>{v.get('Alt')}  varnova='{vv}' annovar='{av}'"
                )

        exonic_true.append(a.get("ExonicFunc.refGene", "."))
        exonic_pred.append(v.get("ExonicFunc.refGene", "."))

    report.append(f"Rows compared: {n}")
    if mismatches["LOCUS_MISALIGNED"]:
        report.append(
            f"LOCUS_MISALIGNED rows (skipped from column comparison): "
            f"{len(mismatches['LOCUS_MISALIGNED'])}"
        )

    headline_pct = {}
    for c in columns:
        pct = 100.0 * matches[c] / n if n else 0.0
        headline_pct[c] = pct
        report.append(f"{c}: {pct:.3f}% concordant ({matches[c]}/{n})")

    exonic_labels = sorted(set(exonic_true) | set(exonic_pred))
    prec, rec, f1 = precision_recall_f1(exonic_true, exonic_pred, exonic_labels)
    report.append(
        f"ExonicFunc.refGene — macro precision={prec:.4f} "
        f"recall={rec:.4f} f1={f1:.4f} (ANNOVAR as pairwise reference, "
        f"not ground truth — see majority_vote.py for a 3-tool comparison)"
    )

    aachange_pct = compare_aachange(varnova_rows, annovar_rows, report)
    headline_pct["AAChange.refGene"] = aachange_pct

    report.append("")
    report.append("Sample disagreements (up to 50 per column):")
    for c in columns:
        if mismatches[c]:
            report.append(f"-- {c} --")
            report.extend(mismatches[c])

    return headline_pct


def parse_vep_csq_header(vcf_path):
    with open_maybe_gz(vcf_path) as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                fmt = line.split("Format:")[-1].strip().rstrip('">')
                return fmt.split("|")
            if not line.startswith("#"):
                break
    return []


def most_severe_csq(csq_field, fields):
    """Pick the single most clinically severe CSQ block for a multi-transcript variant."""
    best = None
    best_rank = len(VEP_SEVERITY_ORDER)
    for block in csq_field.split(","):
        vals = block.split("|")
        rec = dict(zip(fields, vals))
        cons_terms = rec.get("Consequence", "").split("&")
        for term in cons_terms:
            if term in VEP_SEVERITY_ORDER:
                rank = VEP_SEVERITY_ORDER.index(term)
                if rank < best_rank:
                    best_rank = rank
                    best = rec
    return best


def compare_vep(varnova_rows, vep_path, report):
    report.append("")
    report.append("=" * 70)
    report.append("VarNova vs VEP — keyed comparison (most-severe consequence)")
    report.append("=" * 70)

    fields = parse_vep_csq_header(vep_path)
    if not fields:
        report.append("SKIPPED: could not parse CSQ format header from VEP VCF")
        return

    vep_by_key = {}
    with open_maybe_gz(vep_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _id, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            chrom = normalize_chrom(chrom)
            info = cols[7]
            csq_field = None
            for kv in info.split(";"):
                if kv.startswith("CSQ="):
                    csq_field = kv[4:]
                    break
            if csq_field is None:
                continue
            rec = most_severe_csq(csq_field, fields)
            if rec is None:
                continue
            for alt_allele in alt.split(","):
                vep_by_key[(chrom, pos, ref, alt_allele)] = rec

    exonic_true, exonic_pred = [], []
    n_matched, n_unmapped, n_no_key = 0, 0, 0

    for v in varnova_rows:
        chrom = v.get("Chr", "")
        key = (chrom, v.get("Start"), v.get("Ref"), v.get("Alt"))
        rec = vep_by_key.get(key)
        if rec is None:
            n_no_key += 1
            continue
        cons_terms = rec.get("Consequence", "").split("&")
        mapped = None
        for term in cons_terms:
            if term in VEP_TO_EXONIC_FUNC:
                mapped = VEP_TO_EXONIC_FUNC[term]
                break
        if mapped is None:
            n_unmapped += 1
            continue
        n_matched += 1
        exonic_true.append(mapped)
        exonic_pred.append(v.get("ExonicFunc.refGene", "."))

    report.append(f"VarNova rows: {len(varnova_rows)}")
    report.append(f"  matched by Chr:Pos:Ref:Alt + mappable consequence: {n_matched}")
    report.append(f"  no VEP record at this exact locus (expected — VEP doesn't "
                   f"expand multi-allelics the same way): {n_no_key}")
    report.append(f"  VEP consequence term not in our coding-exon mapping table "
                   f"(expected for non-exonic/non-coding classes): {n_unmapped}")

    if n_matched:
        agree = sum(1 for t, p in zip(exonic_true, exonic_pred) if t == p)
        pct = 100.0 * agree / n_matched
        report.append(f"ExonicFunc concordance on matched coding variants: "
                       f"{pct:.3f}% ({agree}/{n_matched})")
        labels = sorted(set(exonic_true) | set(exonic_pred))
        prec, rec, f1 = precision_recall_f1(exonic_true, exonic_pred, labels)
        report.append(f"  macro precision={prec:.4f} recall={rec:.4f} f1={f1:.4f}")


# ExonicFunc.refGene categories plausible for a damaging/pathogenic coding
# variant. Note: a pathogenic variant CAN legitimately land outside this set
# (e.g. a synonymous change that disrupts splicing, or a deep-intronic
# variant), so this is a plausibility signal, not a correctness requirement.
DAMAGING_EXONIC_FUNCS = {
    "nonsynonymous SNV", "stopgain", "stoploss", "startloss",
    "frameshift insertion", "frameshift deletion",
    "nonframeshift insertion", "nonframeshift deletion",
    "frameshift substitution", "nonframeshift substitution",
}


def compare_clinvar(full_annotated_path, report):
    """Section 2 — clinical validation. Uses the full `varnova table` output
    (gene-db + ClinVar filter-db), which carries VarNova's own independently
    -derived Gene.refGene/ExonicFunc.refGene (from refGene + codon logic)
    alongside a ClinVar significance column populated by a *separate* filter
    -DB lookup. Checking the two against each other is not circular: ClinVar
    significance doesn't feed into how Gene.refGene/ExonicFunc.refGene are
    computed (only the ACMG_Classification column's ClinVar-override path
    does that, which isn't used here)."""
    report.append("")
    report.append("=" * 70)
    report.append("Section 2 — Clinical validation against ClinVar Pathogenic/Likely_pathogenic")
    report.append("=" * 70)

    rows = load_tsv(full_annotated_path)
    if not rows or "ClinVar" not in rows[0]:
        report.append("SKIPPED: no 'ClinVar' column found in this file — was it run with "
                       "--filter-db <clinvar.txt>?")
        return

    n_total, n_gene_ok, n_damaging = 0, 0, 0
    examples_no_gene = []
    examples_benign_looking = []

    for r in rows:
        sig = r.get("ClinVar", ".")
        if "Pathogenic" not in sig or "Benign" in sig:
            continue  # keep Pathogenic / Likely_pathogenic / Pathogenic/Likely_pathogenic
        n_total += 1
        gene = r.get("Gene.refGene", "NONE")
        exonic = r.get("ExonicFunc.refGene", ".")
        gene_ok = gene not in ("NONE", ".", "")
        if gene_ok:
            n_gene_ok += 1
        elif len(examples_no_gene) < 10:
            examples_no_gene.append(
                f"  {r.get('Chr')}:{r.get('Start')} {r.get('Ref')}>{r.get('Alt')} "
                f"ClinVar='{sig}' Gene.refGene='{gene}'"
            )
        if exonic in DAMAGING_EXONIC_FUNCS:
            n_damaging += 1
        elif len(examples_benign_looking) < 10:
            examples_benign_looking.append(
                f"  {r.get('Chr')}:{r.get('Start')} {r.get('Ref')}>{r.get('Alt')} "
                f"gene={gene} ClinVar='{sig}' ExonicFunc.refGene='{exonic}' "
                f"(not necessarily wrong — could be a splice-disrupting synonymous "
                f"change or a non-exonic regulatory variant)"
            )

    report.append(f"ClinVar Pathogenic/Likely_pathogenic variants in this run: {n_total}")
    if n_total:
        report.append(f"  assigned a real gene (Gene.refGene != NONE): "
                       f"{100.0*n_gene_ok/n_total:.2f}% ({n_gene_ok}/{n_total})")
        report.append(f"  exonic consequence plausible for damage "
                       f"({sorted(DAMAGING_EXONIC_FUNCS)}): "
                       f"{100.0*n_damaging/n_total:.2f}% ({n_damaging}/{n_total})")
    if examples_no_gene:
        report.append("-- Pathogenic/Likely_pathogenic with no assigned gene --")
        report.extend(examples_no_gene)
    if examples_benign_looking:
        report.append("-- Pathogenic/Likely_pathogenic with a non-damaging-looking ExonicFunc --")
        report.extend(examples_benign_looking)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--varnova", required=True)
    ap.add_argument("--annovar", required=True)
    ap.add_argument("--vep")
    ap.add_argument("--full-annotated",
                     help="optional — output of `varnova table` with a ClinVar filter-db, "
                          "enables the ClinVar clinical-validation section")
    ap.add_argument("--clinvar-gold",
                     help="optional — clinvar_gold_truth.tsv (Chrom/Pos/Ref/Alt/Gene/Transcript/"
                          "HGVSc/HGVSp), enables the ClinVar gold-standard accuracy section")
    ap.add_argument("--out-dir", default=".")
    ap.add_argument("--threshold", type=float, default=CONCORDANCE_THRESHOLD)
    args = ap.parse_args()

    report = []
    varnova_rows = load_tsv(args.varnova)
    annovar_rows = load_tsv(args.annovar)

    headline_pct = compare_annovar(varnova_rows, annovar_rows, report)

    if args.vep:
        compare_vep(varnova_rows, args.vep, report)

    if args.full_annotated:
        compare_clinvar(args.full_annotated, report)

    if args.clinvar_gold:
        compare_clinvar_gold(varnova_rows, args.clinvar_gold, report)

    text = "\n".join(report)
    print(text)

    out_path = f"{args.out_dir}/concordance_report.txt"
    with open(out_path, "w") as f:
        f.write(text + "\n")

    gate_cols = ["Func.refGene", "ExonicFunc.refGene"]
    failing = [c for c in gate_cols if headline_pct.get(c, 0.0) < args.threshold]
    if failing:
        print(f"\nFAIL: concordance below {args.threshold}% for: {', '.join(failing)}",
              file=sys.stderr)
        sys.exit(1)
    print(f"\nPASS: all gated columns >= {args.threshold}% concordance with ANNOVAR")


if __name__ == "__main__":
    main()
