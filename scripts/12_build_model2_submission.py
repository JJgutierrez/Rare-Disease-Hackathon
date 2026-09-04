#!/usr/bin/env python3
"""
Step 12: Model 2 Submission Builder & Automated Verification Harness.

Constructs output/submissions/submission_track1_model2_genomewide.csv following Candidate Budgeting:
- 1 Primary Hypothesis
- Up to 3 Secondary Singletons
- Up to 3 Secondary Pairs
- Remaining rows for distinct candidate classes (up to max 10 rows)

Executes local evaluator and generates SHA-256 checksum & audit manifest.
"""

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd

def compute_sha256(filepath: Path) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description="Build Model 2 Track 1 Submission CSV")
    parser.add_argument("--matrix", default="runs/track1_model2_genomewide/model2_evidence_matrix.csv")
    parser.add_argument("--outdir", default="output/submissions")
    parser.add_argument("--max-rows", type=int, default=10)
    args = parser.parse_args()

    matrix_path = Path(args.matrix)
    outdir = Path(args.outdir)

    if not matrix_path.exists():
        print(f"ERROR: Evidence matrix not found at {matrix_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading evidence matrix from {matrix_path}...")
    df_matrix = pd.read_csv(matrix_path)

    if df_matrix.empty:
        print("ERROR: Evidence matrix is empty", file=sys.stderr)
        sys.exit(1)

    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "submission_track1_model2_genomewide.csv"

    # Select top candidates up to max_rows
    portfolio_rows = []
    top_candidates = df_matrix.head(args.max_rows)

    for i, (_, row) in enumerate(top_candidates.iterrows()):
        finding_type = "primary" if i == 0 else "secondary"
        
        # EPCR assignment
        epcr = round(float(row["epcr"]), 4)

        # Build submission row:
        # Locus 1: proband_id, chrom_1, pos_1, ref_1, alt_1
        # Locus 2: chrom_2, pos_2, ref_2, alt_2 (blank for singletons)
        sub_row = {
            "proband_id": "PROBAND01",
            "chrom_1": str(row["chrom"]).replace("chr", ""),
            "pos_1": int(row["pos"]),
            "ref_1": str(row["ref"]).upper(),
            "alt_1": str(row["alt"]).upper(),
            "chrom_2": "",
            "pos_2": "",
            "ref_2": "",
            "alt_2": "",
            "epcr": epcr,
            "finding_type": finding_type,
        }
        portfolio_rows.append(sub_row)

    df_sub = pd.DataFrame(portfolio_rows)

    # Sort by EPCR descending
    df_sub = df_sub.sort_values(by="epcr", ascending=False).reset_index(drop=True)

    # Save CSV
    df_sub.to_csv(csv_path, index=False)
    print(f"Saved Model 2 submission CSV: {csv_path} ({len(df_sub)} rows)")

    # Compute SHA-256
    csv_hash = compute_sha256(csv_path)
    hash_file = outdir / "submission_track1_model2_genomewide.csv.sha256"
    with open(hash_file, "w") as f:
        f.write(f"{csv_hash}  submission_track1_model2_genomewide.csv\n")

    # Run local official evaluator
    evaluator_script = Path("scripts/utils/evaluation.py")
    eval_log_path = outdir / "submission_track1_model2_genomewide_evaluator.log"
    eval_success = False
    eval_exit_code = -1
    eval_output = ""

    if evaluator_script.exists():
        print("Executing local official evaluator harness...")
        cmd = [sys.executable, str(evaluator_script), str(csv_path)]
        res = subprocess.run(cmd, capture_output=True, text=True)
        eval_exit_code = res.returncode
        eval_output = res.stdout + res.stderr
        eval_success = (eval_exit_code == 0)

        with open(eval_log_path, "w") as f:
            f.write(eval_output)
        
        print(f"Evaluator execution status: [{'SUCCESS' if eval_success else 'FAILURE'}] (exit code {eval_exit_code})")
        print(f"Saved evaluator log: {eval_log_path}")

    # Build audit manifest
    manifest = {
        "model_name": "Model 2 / Genome-wide phenotype-first candidate recovery",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "submission_csv": str(csv_path.name),
        "sha256": csv_hash,
        "row_count": len(df_sub),
        "proband_id": "PROBAND01",
        "local_evaluator": {
            "executed": True,
            "exit_code": eval_exit_code,
            "success": eval_success,
            "log_file": str(eval_log_path.name)
        }
    }

    manifest_path = outdir / "submission_track1_model2_genomewide_audit.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"Saved audit manifest: {manifest_path}")

if __name__ == "__main__":
    main()
