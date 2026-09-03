import pytest
import re
import pandas as pd

ALLOWED_CONTIGS = set([str(i) for i in range(1, 23)] + ["X", "Y", "MT"])
IUPAC_BASES = set(["A", "C", "G", "T"])

def is_bare_contig(contig_str):
    if not contig_str or pd.isna(contig_str):
        return True
    return str(contig_str) in ALLOWED_CONTIGS

def is_valid_dna(base_str):
    if not base_str or pd.isna(base_str):
        return True
    return str(base_str).upper() in IUPAC_BASES

def test_contig_bare_formatting():
    valid_contigs = ["15", "5", "X", "Y"]
    invalid_contigs = ["chr15", "chr5", "chrX", "CHR1"]
    
    for c in valid_contigs:
        assert is_bare_contig(c)
        
    for c in invalid_contigs:
        assert not is_bare_contig(c)

def test_allele_base_formatting():
    valid_alleles = ["A", "C", "G", "T"]
    invalid_alleles = ["N", "AC", "chr", "1"]
    
    for a in valid_alleles:
        assert is_valid_dna(a)
        
    for a in invalid_alleles:
        assert not is_valid_dna(a)

def test_epcr_tie_breaking_sort():
    # Test sorting by epcr DESC, then proband_id, chrom, pos
    records = [
        {"proband_id": "EX2312012", "chrom_1": "5", "pos_1": 100, "epcr": 0.80},
        {"proband_id": "EX2312012", "chrom_1": "15", "pos_1": 200, "epcr": 0.95},
        {"proband_id": "EX2312012", "chrom_1": "15", "pos_1": 100, "epcr": 0.95},
    ]
    df = pd.DataFrame(records)
    
    sorted_df = df.sort_values(
        by=["epcr", "proband_id", "chrom_1", "pos_1"],
        ascending=[False, True, True, True]
    )
    
    pos_order = sorted_df["pos_1"].tolist()
    assert pos_order == [100, 200, 100]  # epcr=0.95 pos 100 before epcr=0.95 pos 200, then epcr=0.80
