#!/usr/bin/env python3
"""
Step 10: Genome-Wide Unbiased Candidate Inventory & Exclusion Audit Builder.

Extracts all variants from WGS_EX2312012_HGWCNDSX7.vcf.gz across the entire genome
without pre-filtering or panel restrictions.

Outputs (Protected local workspace):
- runs/track1_model2_genomewide/model2_all_variants.parquet
- runs/track1_model2_genomewide/model2_exclusions.parquet
"""

import argparse
import json
import sys
import time
from pathlib import Path
import pandas as pd
import numpy as np
from cyvcf2 import VCF

def normalize_contig(chrom: str) -> str:
    chrom_str = str(chrom).strip()
    if chrom_str.startswith("chr"):
        chrom_str = chrom_str[3:]
    return chrom_str

def parse_genotype(gt_tuple):
    # gt_tuple in cyvcf2 is array of alleles: e.g. [0, 1, False]
    if len(gt_tuple) < 2:
        return "./."
    a1, a2 = gt_tuple[0], gt_tuple[1]
    if a1 < 0 or a2 < 0:
        return "./."
    phased = "|" if (len(gt_tuple) >= 3 and gt_tuple[2]) else "/"
    return f"{a1}{phased}{a2}"

def build_inventory(vcf_path: Path, output_dir: Path):
    print(f"Opening source VCF: {vcf_path}")
    vcf = VCF(str(vcf_path))
    samples = vcf.samples
    sample_name = samples[0] if samples else "WGS_EX2312012"
    print(f"Target proband sample in callset: {sample_name}")

    all_records = []
    exclusion_records = []

    start_time = time.time()
    count = 0

    for record in vcf:
        count += 1
        if count % 100000 == 0:
            print(f"Processed {count:,} variants... ({time.time() - start_time:.1f}s)")

        chrom = normalize_contig(record.CHROM)
        pos = int(record.POS)
        ref = str(record.REF).strip().upper()
        filter_status = str(record.FILTER) if record.FILTER else "PASS"

        # Format fields
        dps = record.format("DP")
        dp = int(dps[0][0]) if (dps is not None and len(dps) > 0 and len(dps[0]) > 0) else record.INFO.get("DP", 0)

        gqs = record.format("GQ")
        gq = float(gqs[0][0]) if (gqs is not None and len(gqs) > 0 and len(gqs[0]) > 0) else 0.0

        ads = record.format("AD")
        gt_types = record.gt_types[0] if record.gt_types is not None else 3 # 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN

        # Genotype string
        gt_arr = record.genotypes[0] if record.genotypes else [-1, -1, False]
        gt_str = parse_genotype(gt_arr)

        for i, alt in enumerate(record.ALT):
            alt_str = str(alt).strip().upper()
            var_id = f"{chrom}:{pos}:{ref}:{alt_str}"

            # Calculate Allele Balance (AB)
            ab = np.nan
            if ads is not None and len(ads) > 0 and len(ads[0]) > i + 1:
                ref_ad = float(ads[0][0])
                alt_ad = float(ads[0][i + 1])
                total_ad = ref_ad + alt_ad
                if total_ad > 0:
                    ab = alt_ad / total_ad

            # Exclusion checks
            included = True
            exclusion_reason = "RETAINED"

            if gt_types == 0: # HOM_REF
                included = False
                exclusion_reason = "HOM_REF"
            elif gt_types == 3: # UNKNOWN
                included = False
                exclusion_reason = "MISSING_GENOTYPE"
            elif dp < 5:
                included = False
                exclusion_reason = "LOW_DEPTH_DP_LT_5"

            rec = {
                "source_variant_id": var_id,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt_str,
                "filter": filter_status,
                "qual": float(record.QUAL) if record.QUAL is not None else 0.0,
                "gt": gt_str,
                "gt_type": int(gt_types),
                "dp": int(dp),
                "gq": float(gq),
                "allele_balance": float(ab) if not np.isnan(ab) else None,
                "variant_allele_fraction": float(ab) if not np.isnan(ab) else None,
                "included": included,
                "exclusion_reason": exclusion_reason,
            }

            all_records.append(rec)
            if not included:
                exclusion_records.append(rec)

    print(f"Completed extraction of {len(all_records):,} variant alleles from {count:,} VCF records.")

    df_all = pd.DataFrame(all_records)
    df_exclusions = pd.DataFrame(exclusion_records)

    output_dir.mkdir(parents=True, exist_ok=True)
    all_parquet_path = output_dir / "model2_all_variants.parquet"
    exclusions_parquet_path = output_dir / "model2_exclusions.parquet"

    df_all.to_parquet(all_parquet_path, index=False)
    df_exclusions.to_parquet(exclusions_parquet_path, index=False)

    print(f"Saved complete inventory: {all_parquet_path} ({all_parquet_path.stat().st_size / 1e6:.1f} MB)")
    print(f"Saved exclusion audit: {exclusions_parquet_path} ({exclusions_parquet_path.stat().st_size / 1e6:.1f} MB)")

    # Print summary counts
    retained_df = df_all[df_all["included"]]
    print("\nSummary of Extracted Callset:")
    print(f"  Total Alleles: {len(df_all):,}")
    print(f"  Retained Candidates: {len(retained_df):,}")
    print(f"  Excluded Alleles: {len(df_exclusions):,}")
    print("\nExclusion Breakdown:")
    print(df_all["exclusion_reason"].value_counts().to_string())

def main():
    parser = argparse.ArgumentParser(description="Build Model 2 Genome-Wide Variant Inventory & Audit")
    parser.add_argument("--vcf", default="data/raw/WGS_EX2312012_HGWCNDSX7.vcf.gz", help="Path to raw VCF")
    parser.add_argument("--outdir", default="runs/track1_model2_genomewide", help="Output directory")
    args = parser.parse_args()

    vcf_path = Path(args.vcf)
    outdir = Path(args.outdir)

    if not vcf_path.exists():
        print(f"ERROR: Raw VCF not found at {vcf_path}", file=sys.stderr)
        sys.exit(1)

    build_inventory(vcf_path, outdir)

if __name__ == "__main__":
    main()
