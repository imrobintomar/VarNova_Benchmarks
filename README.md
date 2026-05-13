# VarNova Benchmarks

Reproducible benchmarks comparing **VarNova** against ANNOVAR, VEP, and SnpEff for genomic variant annotation.

> **Paper:** VarNova: a high-performance Rust-based genomic variant annotator *(in preparation)*  
> **VarNova:** https://github.com/imrobintomar/varnova

---

## Results (on 48-core server, 62 GB RAM, NVMe + HDD)

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
VarNova used `.vnidx` binary databases on NVMe. ANNOVAR used ANNOVAR flat `.txt` files on HDD.

---

## Run It Yourself

### Requirements

| Tool | Version | Required for |
|------|---------|-------------|
| VarNova | ≥ 0.1.0 | All benchmarks |
| ANNOVAR | 2020+ | ANNOVAR comparison |
| VEP | ≥ 110 | VEP comparison |
| SnpEff | ≥ 5.0 | SnpEff comparison |
| GNU time | any | Timing measurements |
| bcftools | any | (optional) VCF sort |

### Quick start (VarNova only)

```bash
# Clone benchmark repo
git clone https://github.com/imrobintomar/VarNova_Benchmarks.git
cd VarNova_Benchmarks

# Install VarNova
wget https://github.com/imrobintomar/varnova/releases/latest/download/varnova-linux-x86_64.tar.gz
tar -xzf varnova-linux-x86_64.tar.gz
mv varnova ~/.local/bin/

# Generate a test VCF (1000 synthetic variants)
python3 scripts/generate_test_vcf.py --output testdata/test.vcf --variants 1000

# Edit config to set your database paths
cp config.sh my_config.sh
# Edit my_config.sh: set DB=, ANNOVAR_DIR=, etc.

# Run benchmark
bash scripts/run_benchmark.sh my_config.sh
```

### Full benchmark (all 4 tools)

```bash
# 1. Setup all tools (installs VEP via conda, checks SnpEff)
bash scripts/setup.sh

# 2. Download databases (required for ANNOVAR/VEP/SnpEff comparison)
#    Skip if you only want to benchmark VarNova
perl annovar/annotate_variation.pl --downdb refGene humandb/ --buildver hg38
perl annovar/annotate_variation.pl --downdb gnomad41_exome humandb/ --buildver hg38
perl annovar/annotate_variation.pl --downdb clinvar_20240730 humandb/ --buildver hg38
perl annovar/annotate_variation.pl --downdb avsnp151 humandb/ --buildver hg38

# 3. (Optional) Convert to VarNova binary format for max speed
mkdir -p ~/varnova_db
varnova convert humandb/hg38_gnomad41_exome.txt --out-dir ~/varnova_db/
varnova convert humandb/hg38_clinvar_20240730.txt --out-dir ~/varnova_db/
varnova convert humandb/hg38_avsnp151.txt --out-dir ~/varnova_db/

# 4. Edit config.sh with your paths
# 5. Run
bash scripts/run_benchmark.sh
```

### Expected output

```
results/YYYYMMDD_HHMMSS/
├── benchmark_report.txt    ← full log with timing
├── results.csv             ← machine-readable results
├── varnova_gene.tsv        ← VarNova gene annotation output
├── varnova_full/           ← VarNova full pipeline output
│   └── *.annotated.tsv
├── annovar_gene.*          ← ANNOVAR gene output
├── annovar_full.*          ← ANNOVAR full output
├── vep_gene.vcf            ← VEP gene output
├── vep_full.vcf            ← VEP full output
└── snpeff_*.vcf            ← SnpEff output
```

---

## Reproducing the Published Results

The exact benchmark used in the paper:

- **Machine:** 48-core AMD EPYC, 62 GB RAM, 7.3 TB HDD (82 MB/s) + 3.7 TB NVMe
- **Input:** `testdata/SRR36790057_exome_hg38.vcf` (89,052 GATK variants)
- **VarNova databases:** `.vnidx` format on NVMe (`~/varnova_db/`)
- **Other tools:** databases on HDD (standard setup)
- **OS:** Ubuntu 22.04, kernel 6.8.0

```bash
# Exact command used for VarNova (as in paper)
varnova table \
  -i SRR36790057_exome_hg38.vcf \
  --gene-db humandb/hg38_refGene.txt \
  --filter-db hg38_gnomad41_exome.txt,hg38_clinvar_20240730.txt,hg38_avsnp151.txt \
  --cache-dir ~/varnova_db/ \
  --out-dir results/ \
  --threads 48 -v
```

---

## What VarNova output includes (beyond speed)

Every annotated variant automatically gets:

```
Func.refGene | Gene.refGene | ExonicFunc.refGene | AAChange.refGene
GenCC_Classification | ClinGen_Validity | OMIM_Inheritance | OMIM_Gene_MIM
HPO_Disease_Count | Disease_Names | PanelApp_Panels
```

Example for a BRCA1 variant:
```
exonic | BRCA1 | nonsynonymous SNV | BRCA1:NM_007294:exon10:c.981A>T:p.Lys327Asn
Definitive | Definitive | AD | 113705 | 84 | Breast-ovarian cancer, familial 1;... | Hereditary BRCA(Green)
```

---

## Differences from Other Benchmarks

This benchmark is designed to be fair:

1. **Same databases** — all tools use gnomAD 4.1 exome + ClinVar 2024 + dbSNP 151
2. **Same input** — identical VCF file
3. **Same hardware** — all tools on same machine in same session
4. **Cold cache** — first run timing (no warm OS page cache advantage)
5. **Output verification** — variant counts compared across tools

---

## Contributing

If you run this benchmark on your hardware, please open a PR with your results:

1. Fork this repo
2. Run `bash scripts/run_benchmark.sh`
3. Copy `results/*/results.csv` to `community_results/your_machine_name.csv`
4. Add machine specs to `community_results/README.md`
5. Open a PR

---

## Citation

```bibtex
@article{varnova2026,
  title={VarNova: a high-performance Rust-based genomic variant annotator},
  author={Robin Tomar},
  journal={},
  year={2026},
  note={manuscript in preparation}
}
```

---

## License

MIT — see [LICENSE](LICENSE)
