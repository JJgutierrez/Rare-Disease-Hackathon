import importlib.util
import sys
from pathlib import Path
import pytest

root_dir = Path(__file__).resolve().parent.parent
if str(root_dir) not in sys.path:
    sys.path.insert(0, str(root_dir))

def load_module_from_path(module_name, file_path):
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = mod
    spec.loader.exec_module(mod)
    return mod

build_mini_ref = load_module_from_path("script_04", root_dir / "scripts" / "04_build_mini_ref.py")
phase_read_pairs = load_module_from_path("script_07", root_dir / "scripts" / "07_phase_read_pairs.py")
annotate_candidates = load_module_from_path("script_08", root_dir / "scripts" / "08_annotate_candidates.py")

merge_intervals = build_mini_ref.merge_intervals
compute_candidate_id = phase_read_pairs.compute_candidate_id
select_best_transcript = annotate_candidates.select_best_transcript

def test_interval_merging():
    mock_regions = [
        {"contig": "15", "search_start": 40100000, "search_end": 40200000, "gene_symbol": "BUB1B"},
        {"contig": "15", "search_start": 40150000, "search_end": 40250000, "gene_symbol": "BUB1B_EXT"},
        {"contig": "5", "search_start": 800000, "search_end": 900000, "gene_symbol": "TRIP13"}
    ]
    slices = merge_intervals(mock_regions)
    assert len(slices) == 2
    chr15_slice = [s for s in slices if s["chrom"] == "15"][0]
    assert chr15_slice["search_start"] == 40100000
    assert chr15_slice["search_end"] == 40250000
    assert "BUB1B" in chr15_slice["genes"]
    assert "BUB1B_EXT" in chr15_slice["genes"]

def test_candidate_id_reproducibility():
    cid1 = compute_candidate_id("WGS_EX2312012", "15", 40123456, "G", "A", 1)
    cid2 = compute_candidate_id("WGS_EX2312012", "15", 40123456, "G", "A", 1)
    cid3 = compute_candidate_id("WGS_EX2312012", "15", 40123456, "G", "T", 2)
    assert cid1 == cid2
    assert cid1 != cid3

def test_transcript_selection_precedence():
    mock_txs = [
        {"transcript_id": "ENST001", "biotype": "protein_coding", "consequence_terms": ["intron_variant"]},
        {"transcript_id": "ENST002", "mane_select": "NM_001", "biotype": "protein_coding", "consequence_terms": ["missense_variant"]},
    ]
    best_tx, reason = select_best_transcript(mock_txs)
    assert best_tx["transcript_id"] == "ENST002"
    assert reason == "MANE_Select"
