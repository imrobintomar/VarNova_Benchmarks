# VarNova Benchmarks

Reproducible benchmarks comparing **VarNova** against ANNOVAR, VEP, and SnpEff for genomic variant annotation.

> **Paper:** VarNova: a high-performance Rust-based genomic variant annotator *(in preparation)*  
> **Download VarNova binary:** [Releases →](https://github.com/imrobintomar/VarNova_Benchmarks/releases/latest)

---

## Results (48-core server, 62 GB RAM, NVMe + HDD)

### Tier 1 — Gene Annotation Only

| Tool | Wall Time | Speed | vs VarNova |
|------|-----------|-------|-----------|
| **VarNova 0.1** | **0.8 s** | **111,315 v/s** | — |
| ANNOVAR 2020 | 7.3 s | 12,198 v/s | 9.1× slower |
| SnpEff 5.2a | 53.6 s | 1,661 v/s | 67× slower |
| VEP 115.2 | 69.5 s | 1,281 v/s | 87× slower |

### Tier 2 — Full Pipeline (gene + gnomAD + ClinVar + dbSNP)

| Tool | Wall Time | Speed | Output Variants | vs VarNova |
|------|-----------|-------|----------------|-----------|
| **VarNova 0.1** | **18.9 s** | **4,711 v/s** | **90,188** | — |
| ANNOVAR 2020 | 263.8 s | 337 v/s | 90,188 | **14× slower** |
| VEP 115.2 | 202.4 s | 439 v/s | 89,052 | **10.7× slower** |
| SnpEff 5.2a | 1,619 s | 54 v/s | — | **86× slower** |

**Input:** 89,052 variants (GATK HaplotypeCaller, exome, hg38)  
**Databases:** gnomAD 4.1 exome (18 GB) + ClinVar 2024 (1 GB) + dbSNP 151 (29 GB)  
VarNova used `.vnidx` binary databases on NVMe. ANNOVAR/VEP/SnpEff used standard formats on HDD.

---

### Tier 3 — Scalability

**(a) Thread scaling** (1M synthetic variants, reference 48-core machine):

| Threads | Variants/sec | Speedup | Efficiency |
|---|---|---|---|
| 1 | 6,142 | 1.00× | 100.0% |
| 4 | 23,752 | 3.86× | 96.5% |
| 16 | 55,248 | 8.99× | 56.1% |
| 48 | 72,992 | 11.88× | 24.7% |

Efficiency drops sharply past 16 threads — shown here in full, not hidden, since
it's a real property of the current implementation, not a benchmarking artifact.
Full curve: [`results/scalability_results.csv`](results/scalability_results.csv).

**(b) Variant-count scaling** (48 threads, gene-only annotation):

| Variants | Wall Time | Variants/sec | Peak RAM |
|---|---|---|---|
| 1,000,000 | 19.9 s | 50,251 | 680 MB |
| 10,000,000 | 136.5 s | 73,260 | 2.4 GB |
| 50,000,000 | 694.3 s | 72,014 | 11.3 GB |
| 100,000,000 | 1,379.8 s (23 min) | 72,474 | 21.8 GB |

For comparison, on the same 1M-variant input: ANNOVAR 7,751 v/s, SnpEff 14,471 v/s,
VEP 1,880 v/s. Full data: [`results/cohort_scale_results.csv`](results/cohort_scale_results.csv).

Run it: `bash scripts/scalability_benchmark.sh my_config.sh`

---

### Tier 4 — GIAB multi-sample stability

Re-runs the Tier 1 concordance comparison independently on all seven NIST Genome
in a Bottle reference samples (HG001–HG007), to show the numbers above aren't a
one-sample artifact:

| Sample | Variants | Func % | Gene % | ExonicFunc % | AAChange % |
|---|---|---|---|---|---|
| HG001 | 3,893,341 | 99.802 | 98.587 | 99.979 | 98.810 |
| HG002 | 4,518,208 | 99.432 | 98.191 | 99.604 | 98.586 |
| HG003 | 4,000,097 | 99.809 | 98.472 | 99.980 | 98.598 |
| HG004 | 4,031,346 | 99.799 | 98.459 | 99.977 | 98.630 |
| HG005 | 3,856,856 | 99.800 | 98.559 | 99.973 | 98.785 |
| HG006 | 3,839,315 | 99.810 | 98.593 | 99.973 | 98.758 |
| HG007 | 3,859,704 | 99.801 | 98.565 | 99.985 | 98.884 |

Full data: [`results/giab_multisample_results.csv`](results/giab_multisample_results.csv).
Run it: `bash scripts/giab_benchmark.sh my_config.sh` — downloads the official NIST GIAB
v4.2.1 GRCh38 benchmark VCFs automatically (verify the FTP subpaths in the script
against the live listing first; NIST periodically revises release layouts).

---

### Tier 5 — ClinVar gold-standard HGVS accuracy

Compares VarNova's Gene/HGVS.c/HGVS.p calls against ClinVar's **own curated fields**
(not another annotator's output) for 95,599 Pathogenic/Likely_pathogenic,
high-confidence-review-status variants:

| Metric | Result |
|---|---|
| Gene accuracy | 99.563% (95,181/95,599) |
| HGVS.c accuracy | 91.432% (87,408/95,599) |
| HGVS.p accuracy | 71.170% |
| Transcript selection (same transcript ClinVar's submitter used) | 97.270% (92,989/95,599) |

This is the strongest external-reference validation in this repo — see
[`ACMG_VALIDATION_FINDINGS.md`](ACMG_VALIDATION_FINDINGS.md) §8 for the known
disagreement categories (dup-vs-ins notation, intronic `p..` handling, complex
delins cases) before citing this number without context.

Run it: `bash scripts/clinvar_gold_benchmark.sh my_config.sh` (gold set ships in
`testdata/clinvar_gold/` — no download needed).

---

### Tier 6 — ACMG/AMP classification validation

VarNova implements ACMG/AMP classification using the Tavtigian et al. 2018/2020
points-based combining method (ClinGen SVI's current standard), not the static
2015 categorical table. Validation here is **not** a single agreement percentage
against another tool — see [`ACMG_VALIDATION_FINDINGS.md`](ACMG_VALIDATION_FINDINGS.md)
for the full writeup, including:
- Per-criterion rule-correctness and internal logical-consistency tests (68 passing).
- A two-sided comparison against InterVar: where VarNova's combining arithmetic is
  demonstrably more current (PM2 evidence strength — backed by an audit showing
  100% of InterVar's PVS1-driven Pathogenic calls rely on pre-2018 Moderate-tier
  PM2 weighting), and where InterVar implements more of the full criteria set
  (PM1, PS4, BS2) than VarNova currently does.

---

## Run It Yourself

### Requirements

| Tool | Version | Required for |
|------|---------|-------------|
| VarNova | ≥ 0.1.0 | All benchmarks |
| ANNOVAR | 2020+ | ANNOVAR comparison |
| VEP | ≥ 110 | VEP comparison |
| SnpEff | ≥ 5.0 | SnpEff comparison |
| GNU time | any | Timing |
| Python 3 | ≥ 3.8 | Test VCF generation |

### Quick start (VarNova only, 5 minutes)

```bash
# 1. Clone benchmark repo
git clone https://github.com/imrobintomar/VarNova_Benchmarks.git
cd VarNova_Benchmarks

# 2. Download VarNova binary (Linux x86_64)
wget https://github.com/imrobintomar/VarNova_Benchmarks/releases/latest/download/varnova-linux-x86_64.tar.gz
tar -xzf varnova-linux-x86_64.tar.gz
mv varnova ~/.local/bin/ && chmod +x ~/.local/bin/varnova
export PATH="$HOME/.local/bin:$PATH"
varnova --version

# 3. Generate test VCF (1000 synthetic variants — no real data needed)
python3 scripts/generate_test_vcf.py --output testdata/test.vcf --variants 1000

# 4. Edit config.sh — set DB= to your humandb/ directory
cp config.sh my_config.sh
# nano my_config.sh

# 5. Run (VarNova only if ANNOVAR/VEP/SnpEff not installed)
bash scripts/run_benchmark.sh my_config.sh
```

### Full benchmark (all 4 tools)

```bash
# Install all tools
bash scripts/setup.sh

# Download annotation databases (~60 GB)
# Edit config.sh first, then:
perl annovar/annotate_variation.pl --downdb refGene humandb/ --buildver hg38
perl annovar/annotate_variation.pl --downdb gnomad41_exome humandb/ --buildver hg38
perl annovar/annotate_variation.pl --downdb clinvar_20240730 humandb/ --buildver hg38
perl annovar/annotate_variation.pl --downdb avsnp151 humandb/ --buildver hg38

# (Optional) Convert to VarNova binary format for maximum speed
mkdir -p ~/varnova_db
varnova convert humandb/hg38_gnomad41_exome.txt --out-dir ~/varnova_db/
varnova convert humandb/hg38_clinvar_20240730.txt --out-dir ~/varnova_db/
varnova convert humandb/hg38_avsnp151.txt --out-dir ~/varnova_db/

# Run full benchmark
bash scripts/run_benchmark.sh
```

### Expected output

```
results/YYYYMMDD_HHMMSS/
├── benchmark_report.txt    ← full timing log
├── results.csv             ← machine-readable results
├── varnova_gene.tsv        ← VarNova gene annotation
├── *.annotated.tsv         ← VarNova full pipeline output
├── annovar_gene.*
├── annovar_full.*
├── vep_*.vcf
└── snpeff_*.vcf
```

---

## What VarNova output includes

Every variant automatically gets standard ANNOVAR-compatible columns **plus** gene-disease context:

```
Func.refGene | Gene.refGene | ExonicFunc.refGene | AAChange.refGene
GenCC_Classification | ClinGen_Validity | OMIM_Inheritance | OMIM_Gene_MIM
HPO_Disease_Count | Disease_Names | PanelApp_Panels
```

Example for a BRCA1 variant:
```
exonic | BRCA1 | nonsynonymous SNV | BRCA1:NM_007294:exon10:c.981A>T:p.Lys327Asn
Definitive | Definitive | AD | 113705 | 84 | Breast-ovarian cancer;... | Hereditary BRCA(Green)
```

---

## Reproducing the Published Results

The exact hardware and commands used:

- **Machine:** 48-core CPU, 62 GB RAM, 7.3 TB HDD (82 MB/s) + 3.7 TB NVMe
- **OS:** Ubuntu 22.04, kernel 6.8.0
- **VarNova databases:** `.vnidx` on NVMe (`~/varnova_db/`)
- **ANNOVAR/VEP/SnpEff databases:** standard format on HDD

```bash
# VarNova (as benchmarked)
varnova table \
  -i sample_89052variants_hg38.vcf \
  --gene-db humandb/hg38_refGene.txt \
  --filter-db hg38_gnomad41_exome.txt,hg38_clinvar_20240730.txt,hg38_avsnp151.txt \
  --cache-dir ~/varnova_db/ \
  --out-dir results/ \
  --threads 48 -v
```

---

## Differences from Other Benchmarks

This benchmark is designed to be **fair and reproducible**:

1. **Same databases** — all tools use identical data (gnomAD 4.1 + ClinVar 2024 + dbSNP 151)
2. **Same input VCF** — identical variants for all tools
3. **Same hardware** — all tools run on the same machine in the same session
4. **Cold cache** — first-run timing (OS page cache flushed between runs)
5. **Output verification** — variant counts verified across tools (VarNova = ANNOVAR = 90,188)

---

## Contributing Your Results

Run the benchmark on your hardware and share results:

1. Fork this repo
2. Run `bash scripts/run_benchmark.sh`
3. Copy `results/*/results.csv` → `community_results/<machine_name>.csv`
4. Add specs to `community_results/README.md`
5. Open a PR

---

## Citation

```bibtex
@article{varnova2026,
  title   = {VarNova: a high-performance genomic variant annotator},
  author  = {Robin Tomar},
  year    = {2026},
  note    = {manuscript in preparation}
}
```

---

## License

MIT — see [LICENSE](LICENSE)
