import json
import hashlib
from pathlib import Path
import pandas as pd
import pytest

def compute_sha256(file_path):
    p = Path(file_path)
    if not p.exists():
        return None
    sha256_hash = hashlib.sha256()
    with open(p, "rb") as f:
        for byte_block in iter(lambda: f.read(65536), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()

def test_submission_csv_rules_and_structure(max_track1_rows=10):
    candidates = [
        Path("output/submissions/submission_track1_model2_genomewide.csv"),
        Path("releases/track1/model1_v2/submission_track1_v2_final.csv"),
    ]
    csv_path = None
    for c in candidates:
        if c.exists():
            csv_path = c
            break
    assert csv_path is not None, "No submission CSV found for validation"

    submission_df = pd.read_csv(csv_path)
    assert not submission_df.empty, f"Submission CSV {csv_path} is empty"
    
    required_columns = [
        "proband_id",
        "chrom_1",
        "pos_1",
        "ref_1",
        "alt_1",
        "chrom_2",
        "pos_2",
        "ref_2",
        "alt_2",
        "epcr",
        "finding_type",
    ]
    
    # 1. Header, Row Count, and Proband ID assertions
    assert list(submission_df.columns) == required_columns, f"Header mismatch: {list(submission_df.columns)}"
    assert 1 <= len(submission_df) <= max_track1_rows, f"Row count {len(submission_df)} out of range [1, {max_track1_rows}]"
    assert set(submission_df["proband_id"]) == {"PROBAND01"}, f"Invalid proband_id set: {set(submission_df['proband_id'])}"
    assert submission_df["finding_type"].isin({"primary", "secondary"}).all(), "Invalid finding_type values found"

    # 2. Locus 1 Mandatory Valid Content Assertions
    assert submission_df["chrom_1"].notna().all(), "chrom_1 has null values"
    assert submission_df["pos_1"].notna().all(), "pos_1 has null values"
    assert submission_df["ref_1"].notna().all(), "ref_1 has null values"
    assert submission_df["alt_1"].notna().all(), "alt_1 has null values"
    assert (submission_df["chrom_1"].astype(str).str.strip().ne("")).all(), "chrom_1 has blank values"
    assert (submission_df["ref_1"].astype(str).str.strip().ne("")).all(), "ref_1 has blank values"
    assert (submission_df["alt_1"].astype(str).str.strip().ne("")).all(), "alt_1 has blank values"

    # 3. EPCR numeric range and monotonic descending sort
    assert submission_df["epcr"].between(0.0, 1.0).all(), "EPCR values outside (0.0, 1.0]"
    epcr_vals = submission_df["epcr"].tolist()
    assert epcr_vals == sorted(epcr_vals, reverse=True), "EPCR scores are not in monotonic descending order"

    # 4. Bare contig formatting assertion (no 'chr' prefix)
    assert not submission_df["chrom_1"].astype(str).str.startswith("chr").any(), "chrom_1 contains 'chr' prefix"

    # 5. Exact Locus 2 blank cell handling for singletons vs pairs
    singleton = submission_df["chrom_2"].isna() | (submission_df["chrom_2"].astype(str).str.strip() == "")
    if singleton.any():
        assert (
            submission_df.loc[singleton, ["chrom_2", "pos_2", "ref_2", "alt_2"]]
              .fillna("")
              .astype(str)
              .apply(lambda col: col.str.strip().eq(""))
              .all(axis=1)
        ).all(), "Singleton rows contain non-blank locus-2 values"

    paired = ~singleton
    if paired.any():
        assert not submission_df.loc[paired, "chrom_2"].astype(str).str.startswith("chr").any(), "chrom_2 contains 'chr' prefix"
        assert (
            submission_df.loc[paired, ["chrom_2", "pos_2", "ref_2", "alt_2"]]
              .fillna("")
              .astype(str)
              .apply(lambda col: col.str.strip().ne(""))
              .all(axis=1)
        ).all(), "Paired rows contain blank locus-2 values"

def test_cis_pair_exclusion_regression():
    enriched_path = Path("output/submissions/submission_track1_v2_enriched.csv")
    if enriched_path.exists():
        df_enriched = pd.read_csv(enriched_path)
        if "phase_state" in df_enriched.columns:
            cis_count = (df_enriched["phase_state"] == "cis_supported").sum()
            assert cis_count == 0, f"Regression failure: {cis_count} cis_supported pairs present in enriched candidates"
