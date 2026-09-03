import argparse
import hashlib
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
import requests

ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
CORE_MVA_GENES = ["BUB1B", "CEP57", "TRIP13"]

def fetch_canonical_transcript_exons(gene_symbol, cache_dir=Path("cache/ensembl")):
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"lookup_{gene_symbol}.json"

    if cache_file.exists():
        with open(cache_file, "r") as f:
            data = json.load(f)
    else:
        url = f"{ENSEMBL_LOOKUP_URL}/{gene_symbol}?expand=1;content-type=application/json"
        headers = {"Accept": "application/json"}
        
        data = None
        for attempt in range(5):
            try:
                response = requests.get(url, headers=headers, timeout=60)
                response.raise_for_status()
                data = response.json()
                with open(cache_file, "w") as f:
                    json.dump(data, f, indent=2)
                break
            except Exception as e:
                print(f"WARNING: Ensembl lookup for {gene_symbol} failed (attempt {attempt+1}/5): {e}", file=sys.stderr)
                time.sleep(2 ** attempt)

        if not data:
            return None, "fetch_failed", []

    transcripts = data.get("Transcript", [])
    selected_tx = None
    selection_reason = "none"

    for tx in transcripts:
        if tx.get("is_mane_select"):
            selected_tx = tx
            selection_reason = "MANE_Select"
            break

    if not selected_tx:
        for tx in transcripts:
            if tx.get("is_canonical") and tx.get("biotype") == "protein_coding":
                selected_tx = tx
                selection_reason = "Canonical_protein_coding"
                break

    if not selected_tx:
        for tx in transcripts:
            if tx.get("is_canonical"):
                selected_tx = tx
                selection_reason = "Canonical_fallback"
                break

    if not selected_tx and transcripts:
        selected_tx = transcripts[0]
        selection_reason = "First_transcript_fallback"

    if not selected_tx:
        return None, selection_reason, []

    exons = selected_tx.get("Exon", [])
    exons_sorted = sorted(exons, key=lambda x: int(x["start"]))
    
    parsed_exons = []
    for rank, ex in enumerate(exons_sorted, start=1):
        parsed_exons.append({
            "exon_rank": rank,
            "exon_id": ex.get("id"),
            "chrom": str(ex.get("seq_region_name")).replace("chr", ""),
            "start": int(ex.get("start")),
            "end": int(ex.get("end")),
            "strand": int(ex.get("strand", 1))
        })

    return selected_tx.get("id"), selection_reason, parsed_exons

def main():
    parser = argparse.ArgumentParser(description="Exon Depth & CNV Coverage Anomaly Screen")
    parser.add_argument("--bam", default="output/mini_panel.bam", help="Path to mini_panel.bam")
    parser.add_argument("--mini-manifest", default="work/mini_panel_ref_manifest.json", help="Path to mini_panel_ref_manifest.json")
    parser.add_argument("--output-csv", default="output/exon_coverage_profiles.csv", help="Path to write exon_coverage_profiles.csv")
    parser.add_argument("--output-json", default="output/cnv_coverage_anomaly_screen.json", help="Path to write cnv_coverage_anomaly_screen.json")
    args = parser.parse_args()

    bam_path = Path(args.bam)
    manifest_path = Path(args.mini_manifest)
    csv_out = Path(args.output_csv)
    json_out = Path(args.output_json)

    if not bam_path.exists():
        print(f"ERROR: BAM file not found: {bam_path}", file=sys.stderr)
        sys.exit(1)

    if not manifest_path.exists():
        print(f"ERROR: Mini reference manifest not found: {manifest_path}", file=sys.stderr)
        sys.exit(1)

    with open(manifest_path, "r") as f:
        mini_manifest = json.load(f)

    gene_to_mini = mini_manifest.get("gene_to_mini_contig", {})
    bam = pysam.AlignmentFile(str(bam_path), "rb")

    all_exon_profiles = []
    gene_depths = {}

    for gene in CORE_MVA_GENES:
        print(f"Analyzing exon coverage for {gene}...")
        tx_id, selection_reason, exons = fetch_canonical_transcript_exons(gene)
        if not exons:
            print(f"WARNING: No exons found for {gene}", file=sys.stderr)
            continue

        mini_slice = gene_to_mini.get(gene)
        if not mini_slice:
            print(f"WARNING: No mini reference slice mapping for {gene}", file=sys.stderr)
            continue

        mini_contig = mini_slice["mini_contig_id"]
        slice_start = mini_slice["search_start"]

        gene_pos_depths = []

        for ex in exons:
            ex_start = ex["start"]
            ex_end = ex["end"]
            
            mini_start = ex_start - slice_start + 1
            mini_end = ex_end - slice_start + 1

            try:
                cov_tuple = bam.count_coverage(mini_contig, start=mini_start - 1, stop=mini_end, quality_threshold=20)
                pos_depths = [cov_tuple[0][i] + cov_tuple[1][i] + cov_tuple[2][i] + cov_tuple[3][i] for i in range(len(cov_tuple[0]))]
                gene_pos_depths.extend(pos_depths)
            except Exception as e:
                print(f"WARNING: count_coverage failed for {gene} exon {ex['exon_rank']}: {e}", file=sys.stderr)
                pos_depths = []

            if not pos_depths:
                mean_depth = 0.0
                median_depth = 0.0
                b10 = 0.0
                b20 = 0.0
            else:
                mean_depth = float(np.mean(pos_depths))
                median_depth = float(np.median(pos_depths))
                b10 = float(np.mean([1 if d >= 10 else 0 for d in pos_depths]))
                b20 = float(np.mean([1 if d >= 20 else 0 for d in pos_depths]))

            all_exon_profiles.append({
                "gene_symbol": gene,
                "transcript_id": tx_id,
                "transcript_selection": selection_reason,
                "exon_rank": ex["exon_rank"],
                "exon_id": ex["exon_id"],
                "chrom": ex["chrom"],
                "genomic_start": ex_start,
                "genomic_end": ex_end,
                "mini_contig": mini_contig,
                "mini_start": mini_start,
                "mini_end": mini_end,
                "length_bp": ex_end - ex_start + 1,
                "mean_depth": round(mean_depth, 2),
                "median_depth": round(median_depth, 2),
                "breadth_at_10x": round(b10, 4),
                "breadth_at_20x": round(b20, 4),
                "mapq_min_threshold": 20
            })

        gene_depths[gene] = float(np.median(gene_pos_depths)) if gene_pos_depths else 1.0

    panel_median_depth = float(np.median(list(gene_depths.values()))) if gene_depths else 1.0

    results_list = []

    for prof in all_exon_profiles:
        g_symbol = prof["gene_symbol"]
        g_med = gene_depths.get(g_symbol, 1.0)
        ex_med = prof["median_depth"]

        ratio_gene = (ex_med + 0.1) / (g_med + 0.1)
        ratio_panel = (ex_med + 0.1) / (panel_median_depth + 0.1)

        log2_gene = round(math.log2(ratio_gene), 3)
        log2_panel = round(math.log2(ratio_panel), 3)

        prof["gene_median_depth"] = round(g_med, 2)
        prof["panel_median_depth"] = round(panel_median_depth, 2)
        prof["log2_ratio_to_gene_median"] = log2_gene
        prof["log2_ratio_to_panel_median"] = log2_panel

        if prof["mean_depth"] < 5.0:
            status = "no_call"
        elif log2_gene <= -0.8 and log2_panel <= -0.8:
            status = "possible_exon_depth_reduction"
        else:
            status = "normal_coverage"

        prof["screen_status"] = status
        prof["coverage_source"] = "mini_reference_realignment"
        prof["reference_cohort"] = "none_single_sample"
        prof["cnv_interpretation"] = "screening_only"

        results_list.append(prof)

    df_exons = pd.DataFrame(results_list)
    csv_out.parent.mkdir(parents=True, exist_ok=True)
    df_exons.to_csv(csv_out, index=False)

    flagged_reductions = [p for p in results_list if p["screen_status"] == "possible_exon_depth_reduction"]

    json_payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "bam_file": str(bam_path),
        "coverage_source": "mini_reference_realignment",
        "reference_cohort": "none_single_sample",
        "cnv_interpretation": "screening_only",
        "panel_median_depth": round(panel_median_depth, 2),
        "gene_median_depths": {k: round(v, 2) for k, v in gene_depths.items()},
        "total_exons_evaluated": len(results_list),
        "flagged_depth_reductions_count": len(flagged_reductions),
        "flagged_reductions": flagged_reductions,
        "all_exon_profiles": results_list
    }

    json_out.parent.mkdir(parents=True, exist_ok=True)
    with open(json_out, "w") as f:
        json.dump(json_payload, f, indent=2)

    print(f"Exon depth screen complete. Evaluated {len(results_list)} exons across core MVA genes.")
    print(f"Flagged depth reductions: {len(flagged_reductions)}")
    print(f"Profiles saved to {csv_out} and {json_out}")

if __name__ == "__main__":
    main()
