import json
import hashlib
from pathlib import Path
import jsonschema
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

def test_submission_v2_final_csv_validity():
    csv_path = Path("output/submissions/submission_track1_v2_final.csv")
    assert csv_path.exists(), "submission_track1_v2_final.csv does not exist"

    df = pd.read_csv(csv_path)
    assert not df.empty, "Submission CSV is empty"
    assert len(df) <= 10, f"Submission row count exceeds 10: {len(df)}"

    expected_cols = [
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
        "finding_type"
    ]
    assert list(df.columns) == expected_cols, f"Headers do not match strict 11 columns: {list(df.columns)}"

    # Validate non-null mandatory fields
    for col in expected_cols:
        assert df[col].isna().sum() == 0, f"Column {col} contains null values"

    # Validate EPCR range and monotonic descending order
    epcr_vals = df["epcr"].tolist()
    assert df["epcr"].min() >= 0.0, "EPCR score < 0.0"
    assert df["epcr"].max() <= 1.0, "EPCR score > 1.0"
    assert epcr_vals == sorted(epcr_vals, reverse=True), "EPCR scores are not in monotonic descending order"

    # Validate allowed finding_type values
    for ft in df["finding_type"]:
        assert ft in ["primary", "secondary"], f"Invalid finding_type: {ft}"

def test_cis_pair_exclusion_regression():
    enriched_path = Path("output/submissions/submission_track1_v2_enriched.csv")
    if enriched_path.exists():
        df_enriched = pd.read_csv(enriched_path)
        if "phase_state" in df_enriched.columns:
            cis_count = (df_enriched["phase_state"] == "cis_supported").sum()
            assert cis_count == 0, f"Regression failure: {cis_count} cis_supported pairs present in enriched candidates"

def test_vep_gene_symbol_matching_regression():
    enriched_path = Path("output/submissions/submission_track1_v2_enriched.csv")
    assert enriched_path.exists(), "submission_track1_v2_enriched.csv does not exist"
    df_enriched = pd.read_csv(enriched_path)

    # Assert every pair is grouped by VEP assigned gene symbol
    for idx, row in df_enriched.iterrows():
        assert row["gene_symbol"] == row["vep_gene_symbol"], f"Gene assignment mismatch: pair gene {row['gene_symbol']} != vep_gene {row['vep_gene_symbol']}"

def test_padded_search_window_boundary_exclusion_regression():
    enriched_path = Path("output/submissions/submission_track1_v2_enriched.csv")
    df_enriched = pd.read_csv(enriched_path)

    # Assert no 50kb upstream ZDHHC11 variants (5:842896-845900) are mislabeled as TRIP13
    trip13_pairs = df_enriched[df_enriched["gene_symbol"] == "TRIP13"]
    for idx, row in trip13_pairs.iterrows():
        assert int(row["variant1_pos"]) >= 890000, f"TRIP13 variant 1 pos {row['variant1_pos']} outside TRIP13 locus"
        assert int(row["variant2_pos"]) >= 890000, f"TRIP13 variant 2 pos {row['variant2_pos']} outside TRIP13 locus"

def test_audit_json_validity():
    audit_path = Path("output/submissions/submission_track1_v2_audit.json")
    assert audit_path.exists(), "submission_track1_v2_audit.json does not exist"

    with open(audit_path, "r") as f:
        audit_data = json.load(f)

    assert audit_data["submission_filename"] == "submission_track1_v2_final.csv"
    assert "funnel_chain" in audit_data
    assert "input_hashes" in audit_data
    assert "output_hashes" in audit_data
    assert audit_data["official_evaluator"]["executed"] == True
    assert audit_data["official_evaluator"]["exit_code"] == 0

    final_csv_path = Path("output/submissions/submission_track1_v2_final.csv")
    if final_csv_path.exists():
        computed_hash = compute_sha256(final_csv_path)
        assert audit_data["output_hashes"]["submission_track1_v2_final_csv"] == computed_hash, "Audit JSON final CSV hash mismatch"
