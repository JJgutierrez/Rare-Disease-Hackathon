import json
import pytest
import pandas as pd
from pathlib import Path

SPECS_DIR = Path(__file__).parent.parent / "specs"

def load_submission_schema():
    schema_path = SPECS_DIR / "submission_schema.json"
    with open(schema_path, "r") as f:
        return json.load(f)

def load_candidate_schema():
    schema_path = SPECS_DIR / "panel_candidate_schema.json"
    with open(schema_path, "r") as f:
        return json.load(f)

def test_submission_schema_structure():
    schema = load_submission_schema()
    assert "required_columns" in schema
    assert len(schema["required_columns"]) == 11
    assert schema["required_columns"][0] == "proband_id"
    assert schema["required_columns"][-1] == "finding_type"

def test_candidate_schema_structure():
    schema = load_candidate_schema()
    assert "required_columns" in schema
    assert "filter_decision" in schema["required_columns"]
    assert "filter_reason" in schema["required_columns"]
    assert "gt" in schema["required_columns"]

def test_mock_submission_validation():
    schema = load_submission_schema()
    mock_data = {
        "proband_id": ["EX2312012", "EX2312012"],
        "chrom_1": ["15", "5"],
        "pos_1": [40123456, 50123456],
        "ref_1": ["G", "A"],
        "alt_1": ["A", "G"],
        "chrom_2": ["15", ""],
        "pos_2": [40123490, None],
        "ref_2": ["C", ""],
        "alt_2": ["T", ""],
        "epcr": [0.95, 0.80],
        "finding_type": ["primary", "secondary"]
    }
    df = pd.DataFrame(mock_data)
    
    # Verify column order and exact header match
    assert list(df.columns) == schema["required_columns"]
    
    # Verify row count bounds (1 <= N <= 10)
    assert schema["row_count"]["min"] <= len(df) <= schema["row_count"]["max"]
    
    # Verify monotonic descending EPCR
    epcr_vals = df["epcr"].tolist()
    assert epcr_vals == sorted(epcr_vals, reverse=True)

def test_invalid_header_order_fails():
    schema = load_submission_schema()
    mock_data = {
        "chrom_1": ["15"],
        "proband_id": ["EX2312012"],
        "pos_1": [40123456],
        "ref_1": ["G"],
        "alt_1": ["A"],
        "chrom_2": [""],
        "pos_2": [None],
        "ref_2": [""],
        "alt_2": [""],
        "epcr": [0.95],
        "finding_type": ["primary"]
    }
    df = pd.DataFrame(mock_data)
    assert list(df.columns) != schema["required_columns"]
