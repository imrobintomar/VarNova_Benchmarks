#!/usr/bin/env python3
"""
Generate a minimal test VCF with known variants for benchmarking.
Uses real variants from ClinVar/1000G to ensure realistic annotation results.
"""
import argparse, random

HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr17,length=83257441>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""

# Real known variants (ClinVar pathogenic + common population variants)
KNOWN_VARIANTS = [
    # BRCA1 missense
    ("chr17", 43045733, "rs28897696", "A", "T", "BRCA1 pathogenic"),
    # BRCA2 nonsense
    ("chr13", 32338271, "rs80359550", "A", "T", "BRCA2 pathogenic"),
    # TP53 hotspot
    ("chr17", 7674220, "rs28934578", "C", "T", "TP53 Li-Fraumeni"),
    # CFTR common variant
    ("chr7", 117559590, "rs113993960", "ATT", "A", "CFTR deletion"),
    # Wilson disease (ATP7B)
    ("chr13", 52517374, "rs76151636", "G", "A", "ATP7B Wilson disease"),
    # Common synonymous
    ("chr1", 69270, "rs2691305", "A", "G", "OR4F5 synonymous"),
    ("chr1", 69511, "rs2691276", "A", "G", "OR4F5 missense"),
]

def generate_vcf(output, n_variants=1000, seed=42):
    random.seed(seed)
    with open(output, 'w') as f:
        f.write(HEADER)
        written = 0

        # Write known variants first
        for chrom, pos, rsid, ref, alt, note in KNOWN_VARIANTS:
            f.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t100\tPASS\t.\tGT\t0/1\n")
            written += 1

        # Fill remaining with synthetic variants spread across genome
        chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
        bases = "ACGT"
        while written < n_variants:
            chrom = random.choice(chroms)
            pos   = random.randint(100000, 100000000)
            ref   = random.choice(bases)
            alt   = random.choice([b for b in bases if b != ref])
            qual  = random.randint(30, 200)
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\tPASS\t.\tGT\t0/1\n")
            written += 1

    print(f"Generated {written} variants → {output}")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--output", required=True)
    p.add_argument("--variants", type=int, default=1000)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()
    generate_vcf(args.output, args.variants, args.seed)
