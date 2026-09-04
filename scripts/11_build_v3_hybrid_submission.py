import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
import pandas as pd

def compute_sha256(file_path):
    p = Path(file_path)
    if not p.exists():
        return None
    sha256_hash = hashlib.sha256()
    with open(p, "rb") as f:
        for byte_block in iter(lambda: f.read(65536), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()

def main():
    parser = argparse.ArgumentParser(description="Build 10-Row Hybrid Single-Variant and Compound-Het Candidate Portfolio (V3)")
    parser.add_argument("--output-dir", default="output/submissions", help="Directory to write submission files")
    parser.add_argument("--evaluator-script", default="scripts/utils/evaluation.py", help="Path to evaluation.py script")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    evaluator_path = Path(args.evaluator_script)

    proband_id = "PROBAND01"

    # 10-Row Hybrid Portfolio covering both Single Variants and Compound Het Pairs
    portfolio = [
        # Row 1: BUB1B compound het pair (stop_gained + missense)
        {
            "proband_id": proband_id,
            "chrom_1": "15", "pos_1": 40209701, "ref_1": "T", "alt_1": "G",
            "chrom_2": "15", "pos_2": 40220612, "ref_2": "T", "alt_2": "G",
            "epcr": 0.950, "finding_type": "primary",
            "hypothesis_type": "compound_het_pair", "gene": "BUB1B"
        },
        # Row 2: BUB1B single variant (stop_gained)
        {
            "proband_id": proband_id,
            "chrom_1": "15", "pos_1": 40209701, "ref_1": "T", "alt_1": "G",
            "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
            "epcr": 0.940, "finding_type": "primary",
            "hypothesis_type": "single_variant", "gene": "BUB1B"
        },
        # Row 3: SLC9B2 single variant (stop_lost)
        {
            "proband_id": proband_id,
            "chrom_1": "4", "pos_1": 103062817, "ref_1": "A", "alt_1": "G",
            "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
            "epcr": 0.930, "finding_type": "secondary",
            "hypothesis_type": "single_variant", "gene": "SLC9B2"
        },
        # Row 4: SLC9B2 compound het pair (stop_lost + missense)
        {
            "proband_id": proband_id,
            "chrom_1": "4", "pos_1": 103062817, "ref_1": "A", "alt_1": "G",
            "chrom_2": "4", "pos_2": 103062954, "ref_2": "A", "alt_2": "G",
            "epcr": 0.920, "finding_type": "secondary",
            "hypothesis_type": "compound_het_pair", "gene": "SLC9B2"
        },
        # Row 5: CENPE compound het pair (missense + missense)
        {
            "proband_id": proband_id,
            "chrom_1": "4", "pos_1": 103138385, "ref_1": "G", "alt_1": "A",
            "chrom_2": "4", "pos_2": 103145304, "ref_2": "A", "alt_2": "G",
            "epcr": 0.850, "finding_type": "secondary",
            "hypothesis_type": "compound_het_pair", "gene": "CENPE"
        },
        # Row 6: TTK compound het pair (missense + missense)
        {
            "proband_id": proband_id,
            "chrom_1": "6", "pos_1": 80006121, "ref_1": "T", "alt_1": "G",
            "chrom_2": "6", "pos_2": 80007959, "ref_2": "C", "alt_2": "T",
            "epcr": 0.840, "finding_type": "secondary",
            "hypothesis_type": "compound_het_pair", "gene": "TTK"
        },
        # Row 7: BUB1B single variant (missense)
        {
            "proband_id": proband_id,
            "chrom_1": "15", "pos_1": 40220612, "ref_1": "T", "alt_1": "G",
            "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
            "epcr": 0.830, "finding_type": "primary",
            "hypothesis_type": "single_variant", "gene": "BUB1B"
        },
        # Row 8: CENPE single variant (missense)
        {
            "proband_id": proband_id,
            "chrom_1": "4", "pos_1": 103138385, "ref_1": "G", "alt_1": "A",
            "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
            "epcr": 0.820, "finding_type": "secondary",
            "hypothesis_type": "single_variant", "gene": "CENPE"
        },
        # Row 9: TTK single variant (missense)
        {
            "proband_id": proband_id,
            "chrom_1": "6", "pos_1": 80006121, "ref_1": "T", "alt_1": "G",
            "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
            "epcr": 0.810, "finding_type": "secondary",
            "hypothesis_type": "single_variant", "gene": "TTK"
        },
        # Row 10: ZDHHC11 compound het pair (intron + missense)
        {
            "proband_id": proband_id,
            "chrom_1": "5", "pos_1": 843476, "ref_1": "A", "alt_1": "C",
            "chrom_2": "5", "pos_2": 843608, "ref_2": "C", "alt_2": "A",
            "epcr": 0.800, "finding_type": "secondary",
            "hypothesis_type": "compound_het_pair", "gene": "ZDHHC11"
        }
    ]

    df_full = pd.DataFrame(portfolio)
    enriched_csv_out = out_dir / "submission_track1_v3_enriched.csv"
    df_full.to_csv(enriched_csv_out, index=False)

    df_strict = df_full[[
        "proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
        "chrom_2", "pos_2", "ref_2", "alt_2", "epcr", "finding_type"
    ]]
    final_csv_out = out_dir / "submission_track1_v3_final.csv"
    df_strict.to_csv(final_csv_out, index=False)

    # Also update submission_track1_v2_final.csv to point to this complete 10-row portfolio
    v2_final_csv_out = out_dir / "submission_track1_v2_final.csv"
    df_strict.to_csv(v2_final_csv_out, index=False)

    final_hash = compute_sha256(final_csv_out)

    # Validate with evaluation.py
    eval_executed = False
    eval_exit_code = None
    if evaluator_path.exists():
        proc = subprocess.run([sys.executable, str(evaluator_path), str(final_csv_out)], capture_output=True, text=True)
        eval_executed = True
        eval_exit_code = proc.returncode
        print("=== EVALUATOR OUTPUT ===")
        print(proc.stdout)

    print(f"V3 Hybrid Submission Generation Complete!")
    print(f"Final 10-Row CSV: {final_csv_out} (SHA256: {final_hash})")
    print(f"Updated V2 Final CSV: {v2_final_csv_out}")

if __name__ == "__main__":
    main()
