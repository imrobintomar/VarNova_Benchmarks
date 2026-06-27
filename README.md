# VarNova Benchmarks

Reproducible benchmarks comparing **VarNova** against ANNOVAR, VEP, and SnpEff for genomic variant annotation.

> **Paper:** VarNova: a high-performance Rust-based genomic variant annotator *(in preparation)*  
> **Download VarNova binary:** [Releases →](https://github.com/imrobintomar/VarNova_Benchmarks/releases/latest)

---

## Results (48-core server, 62 GB RAM, NVMe + HDD)

### Tier 1 - Gene Annotation Only

| Tool | Wall Time | Speed | vs VarNova |
|------|-----------|-------|-----------|
| **VarNova 0.1** | **0.8 s** | **111,315 v/s** | - |
| ANNOVAR 2020 | 7.3 s | 12,198 v/s | 9.1× slower |
| SnpEff 5.2a | 53.6 s | 1,661 v/s | 67× slower |
| VEP 115.2 | 69.5 s | 1,281 v/s | 87× slower |

### Tier 2 - Full Pipeline (gene + gnomAD + ClinVar + dbSNP)

| Tool | Wall Time | Speed | Output Variants | vs VarNova |
|------|-----------|-------|----------------|-----------|
| **VarNova 0.1** | **18.9 s** | **4,711 v/s** | **90,188** | - |
| ANNOVAR 2020 | 263.8 s | 337 v/s | 90,188 | **14× slower** |
| VEP 115.2 | 202.4 s | 439 v/s | 89,052 | **10.7× slower** |
| SnpEff 5.2a | 1,619 s | 54 v/s | - | **86× slower** |

**Input:** 89,052 variants (GATK HaplotypeCaller, exome, hg38)  
**Databases:** gnomAD 4.1 exome (18 GB) + ClinVar 2024 (1 GB) + dbSNP 151 (29 GB)  
VarNova used `.vnidx` binary databases on NVMe. ANNOVAR/VEP/SnpEff used standard formats on HDD.

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

# 3. Generate test VCF (1000 synthetic variants - no real data needed)
python3 scripts/generate_test_vcf.py --output testdata/test.vcf --variants 1000

# 4. Edit config.sh - set DB = to your humandb/ directory
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

1. **Same databases** - all tools use identical data (gnomAD 4.1 + ClinVar 2024 + dbSNP 151)
2. **Same input VCF** - identical variants for all tools
3. **Same hardware** - all tools run on the same machine in the same session
4. **Cold cache** - first-run timing (OS page cache flushed between runs)
5. **Output verification** - variant counts verified across tools (VarNova = ANNOVAR = 90,188)

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

MIT - see [LICENSE](LICENSE)
