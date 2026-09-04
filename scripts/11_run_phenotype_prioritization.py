#!/usr/bin/env python3
"""
Step 11: Phenotype-First Ranking & Vectorized Implementation-Safe Evidence Score Builder.

Prioritizes genome-wide variant candidates from model2_all_variants.parquet using HPO terms
derived from Challenge_Clinical_Phenotype_1.docx.

Vectorized for ultra-fast processing across millions of candidates.
"""

import argparse
import json
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# MVA and phenotype-associated high priority genes
GENE_PHENOTYPE_SCORES = {
    "BUB1B": 1.00,  # MVA1 (Mosaic Variegated Aneuploidy Syndrome 1)
    "CEP57": 0.95,  # MVA2 (Mosaic Variegated Aneuploidy Syndrome 2)
    "TRIP13": 0.95, # MVA3 (Mosaic Variegated Aneuploidy Syndrome 3)
    "MAD2L1": 0.85,
    "CDC20": 0.85,
    "BUB1": 0.80,
    "BUBR1": 0.80,
    "TTK": 0.75,
    "MPS1": 0.75,
    "PLK1": 0.75,
    "AURKA": 0.70,
    "AURKB": 0.70,
    "CENPE": 0.70,
    "CENPF": 0.70,
    "KIF11": 0.70,
    "ASPM": 0.65,
    "WDR62": 0.65,
    "MCPH1": 0.65,
    "WT1": 0.60,
}

CONSEQUENCE_WEIGHTS = {
    "stop_gained": 1.0,
    "frameshift_variant": 1.0,
    "splice_acceptor_variant": 0.95,
    "splice_donor_variant": 0.95,
    "start_lost": 0.90,
    "stop_lost": 0.85,
    "inframe_insertion": 0.75,
    "inframe_deletion": 0.75,
    "missense_variant": 0.60,
    "protein_altering_variant": 0.55,
    "splice_region_variant": 0.50,
    "synonymous_variant": 0.10,
    "intron_variant": 0.05,
    "5_prime_utr_variant": 0.10,
    "3_prime_utr_variant": 0.10,
    "upstream_gene_variant": 0.05,
    "downstream_gene_variant": 0.05,
}

def main():
    parser = argparse.ArgumentParser(description="Phenotype-First Ranking & Evidence Matrix Builder")
    parser.add_argument("--candidates", default="runs/track1_model2_genomewide/model2_all_variants.parquet")
    parser.add_argument("--outdir", default="runs/track1_model2_genomewide")
    args = parser.parse_args()

    candidates_path = Path(args.candidates)
    outdir = Path(args.outdir)

    if not candidates_path.exists():
        print(f"ERROR: Candidates parquet file not found at {candidates_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading genome-wide candidate inventory from {candidates_path}...")
    df = pd.read_parquet(candidates_path)

    # Filter to retained candidates
    retained = df[df["included"]].copy()
    print(f"Retained candidates for evidence scoring: {len(retained):,}")

    # Vectorized P: Phenotype Score
    gene_col = retained.get("gene", pd.Series([""] * len(retained), index=retained.index)).astype(str).str.upper()
    P = gene_col.map(GENE_PHENOTYPE_SCORES).fillna(0.40)

    # Vectorized F: Functional Consequence
    csq_col = retained.get("consequence", pd.Series([""] * len(retained), index=retained.index)).astype(str).str.lower()
    F = csq_col.map(CONSEQUENCE_WEIGHTS).fillna(0.30)

    # Vectorized I: Genotype
    gt_type = retained.get("gt_type", pd.Series([1] * len(retained), index=retained.index))
    I = np.where(gt_type == 2, 1.0, np.where(gt_type == 1, 0.70, 0.10))

    # Vectorized Q: Quality
    dp = retained.get("dp", pd.Series([30] * len(retained), index=retained.index)).astype(float)
    gq = retained.get("gq", pd.Series([99.0] * len(retained), index=retained.index)).astype(float)
    ab = retained.get("allele_balance", pd.Series([0.5] * len(retained), index=retained.index)).fillna(0.5).astype(float)

    dp_score = np.clip(dp / 30.0, 0.0, 1.0)
    gq_score = np.clip(gq / 99.0, 0.0, 1.0)
    ab_score = 1.0 - np.abs(0.5 - ab) * 2.0
    Q = (dp_score + gq_score + ab_score) / 3.0

    # Rarity (R) & Phase (H)
    R = 0.90
    w_p, w_f, w_i, w_q, w_r, w_h = 0.30, 0.25, 0.15, 0.15, 0.10, 0.05

    phase_status = retained.get("phase_status", pd.Series(["unresolved"] * len(retained), index=retained.index))
    H = np.where(phase_status == "confirmed_trans", 1.0, 0.0)

    # Raw evidence calculation
    raw_evidence = (w_p * P) + (w_f * F) + (w_i * I) + (w_q * Q) + (w_r * R) + (w_h * H)
    
    # Binary phase exclusion gate: confirmed cis pairs excluded (score -1)
    raw_evidence = np.where(phase_status == "confirmed_cis", -1.0, raw_evidence)

    retained["raw_evidence"] = raw_evidence

    # Filter out binary phase excluded
    eligible = retained[retained["raw_evidence"] >= 0.0].copy()
    eligible = eligible.sort_values(by="raw_evidence", ascending=False)

    # Map raw evidence to monotonically decreasing epcr in (0, 1]
    if not eligible.empty:
        max_raw = eligible["raw_evidence"].max()
        min_raw = eligible["raw_evidence"].min()
        if max_raw > min_raw:
            eligible["epcr"] = 0.50 + 0.49 * ((eligible["raw_evidence"] - min_raw) / (max_raw - min_raw))
        else:
            eligible["epcr"] = 0.90

    outdir.mkdir(parents=True, exist_ok=True)
    matrix_path = outdir / "model2_evidence_matrix.csv"
    eligible.to_csv(matrix_path, index=False)
    print(f"Saved evidence matrix: {matrix_path} ({len(eligible):,} eligible candidates)")

    # Print top 15 ranked candidates
    top_cols = ["chrom", "pos", "ref", "alt", "gt", "dp", "gq", "raw_evidence", "epcr"]
    print("\nTop 15 Evidence-Ranked Genome-Wide Candidates:")
    print(eligible[top_cols].head(15).to_string())

if __name__ == "__main__":
    main()
