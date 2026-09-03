import argparse
import hashlib
import json
import time
import sys
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
import requests

VEP_POST_URL = "https://rest.ensembl.org/vep/homo_sapiens/region"

CONSEQUENCE_SEVERITY_RANK = {
    "transcript_ablation": 1,
    "splice_acceptor_variant": 2,
    "splice_donor_variant": 3,
    "stop_gained": 4,
    "frameshift_variant": 5,
    "stop_lost": 6,
    "start_lost": 7,
    "transcript_amplification": 8,
    "inframe_insertion": 9,
    "inframe_deletion": 10,
    "missense_variant": 11,
    "protein_altering_variant": 12,
    "splice_region_variant": 13,
    "incomplete_terminal_codon_variant": 14,
    "stop_retained_variant": 15,
    "synonymous_variant": 16,
    "coding_sequence_variant": 17,
    "mature_miRNA_variant": 18,
    "5_prime_UTR_variant": 19,
    "3_prime_UTR_variant": 20,
    "non_coding_transcript_exon_variant": 21,
    "intron_variant": 22,
    "NMD_transcript_variant": 23,
    "non_coding_transcript_variant": 24,
    "upstream_gene_variant": 25,
    "downstream_gene_variant": 26,
    "tfbs_ablation": 27,
    "tfbs_amplification": 28,
    "tf_binding_site_variant": 29,
    "regulatory_region_ablation": 30,
    "regulatory_region_amplification": 31,
    "feature_elongation": 32,
    "feature_truncation": 33,
    "regulatory_region_variant": 34,
    "intergenic_variant": 35
}

def post_vep_batch(variants_chunk, max_retries=5):
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }
    payload = {"variants": variants_chunk}
    
    for attempt in range(max_retries):
        try:
            resp = requests.post(VEP_POST_URL, headers=headers, json=payload, timeout=60)
            if resp.status_code == 429:  # Rate limited
                retry_after = int(resp.headers.get("Retry-After", 2))
                print(f"Rate limited (429). Retrying after {retry_after} seconds...")
                time.sleep(retry_after)
                continue
            resp.raise_for_status()
            return resp.json()
        except Exception as e:
            if attempt == max_retries - 1:
                print(f"ERROR: VEP POST request failed after {max_retries} attempts: {e}", file=sys.stderr)
                raise
            time.sleep(2 ** attempt)

def select_best_transcript(transcripts_list):
    if not transcripts_list:
        return None, "none"

    # 1. MANE Select
    for tx in transcripts_list:
        if tx.get("mane_select") or tx.get("mane"):
            return tx, "MANE_Select"

    # 2. MANE Plus Clinical
    for tx in transcripts_list:
        if tx.get("mane_plus_clinical"):
            return tx, "MANE_Plus_Clinical"

    # 3. Canonical protein coding
    for tx in transcripts_list:
        if tx.get("is_canonical") and tx.get("biotype") == "protein_coding":
            return tx, "Canonical_protein_coding"

    # 4. Any protein coding with most severe consequence
    pc_txs = [tx for tx in transcripts_list if tx.get("biotype") == "protein_coding"]
    if pc_txs:
        pc_txs.sort(key=lambda t: min([CONSEQUENCE_SEVERITY_RANK.get(c, 99) for c in t.get("consequence_terms", ["intergenic_variant"])]))
        return pc_txs[0], "Most_severe_protein_coding"

    # 5. Fallback most severe
    transcripts_list.sort(key=lambda t: min([CONSEQUENCE_SEVERITY_RANK.get(c, 99) for c in t.get("consequence_terms", ["intergenic_variant"])]))
    return transcripts_list[0], "Most_severe_transcript_fallback"

def main():
    parser = argparse.ArgumentParser(description="Adaptive Batched VEP & Predictor Candidate Annotation")
    parser.add_argument("--candidates-parquet", default="output/panel_candidates_raw.parquet", help="Path to panel_candidates_raw.parquet")
    parser.add_argument("--output-parquet", default="output/panel_candidates_annotated.parquet", help="Path to write panel_candidates_annotated.parquet")
    parser.add_argument("--output-csv", default="output/panel_candidates_annotated.csv", help="Path to write panel_candidates_annotated.csv")
    parser.add_argument("--output-tx-parquet", default="output/panel_candidate_transcripts.parquet", help="Path to write panel_candidate_transcripts.parquet")
    parser.add_argument("--manifest", default="output/annotation_manifest.json", help="Path to write annotation_manifest.json")
    parser.add_argument("--batch-size", type=int, default=50, help="VEP REST batch size")
    args = parser.parse_args()

    parquet_in = Path(args.candidates_parquet)
    parquet_out = Path(args.output_parquet)
    csv_out = Path(args.output_csv)
    tx_parquet_out = Path(args.output_tx_parquet)
    manifest_out = Path(args.manifest)

    if not parquet_in.exists():
        print(f"ERROR: Candidate parquet not found: {parquet_in}", file=sys.stderr)
        sys.exit(1)

    df_candidates = pd.read_parquet(parquet_in)
    if df_candidates.empty:
        print("WARNING: Candidate table is empty. Exiting annotation.", file=sys.stderr)
        sys.exit(0)

    print(f"Annotating {len(df_candidates)} candidate alleles via Ensembl VEP REST API...")

    # Format VEP region queries: e.g. "15 40123456 40123456 G/A 1"
    vep_variant_strings = []
    cand_id_map = {}

    for idx, row in df_candidates.iterrows():
        chrom = str(row["chrom"]).replace("chr", "")
        pos = int(row["pos"])
        ref = str(row["ref"])
        alt = str(row["alt"])
        cand_id = row.get("candidate_id") if pd.notna(row.get("candidate_id")) else f"{chrom}_{pos}_{ref}_{alt}"

        # VEP region string format: chrom start end allele_string strand
        # Handle indels vs SNVs
        if len(ref) == 1 and len(alt) == 1:
            start = pos
            end = pos
            allele_str = f"{ref}/{alt}"
        elif len(ref) > len(alt):  # Deletion
            start = pos + 1
            end = pos + len(ref) - 1
            allele_str = "-/" + alt if alt != "" else "-"
        else:  # Insertion
            start = pos
            end = pos + 1
            allele_str = f"-/{alt[1:]}"

        v_str = f"{chrom} {start} {end} {allele_str} 1"
        vep_variant_strings.append(v_str)
        cand_id_map[v_str] = cand_id

    # Adaptive batching
    batch_size = args.batch_size
    all_vep_responses = []

    cache_dir = Path("cache/vep")
    cache_dir.mkdir(parents=True, exist_ok=True)

    for i in range(0, len(vep_variant_strings), batch_size):
        chunk = vep_variant_strings[i:i + batch_size]
        print(f"Annotating batch {i // batch_size + 1}/{(len(vep_variant_strings) + batch_size - 1) // batch_size} ({len(chunk)} variants)...")
        
        # Check cache
        chunk_hash = hashlib.sha256(json.dumps(chunk).encode("utf-8")).hexdigest()
        cache_file = cache_dir / f"vep_batch_{chunk_hash[:12]}.json"

        if cache_file.exists():
            with open(cache_file, "r") as f:
                res_data = json.load(f)
        else:
            res_data = post_vep_batch(chunk)
            with open(cache_file, "w") as f:
                json.dump(res_data, f, indent=2)

        if isinstance(res_data, list):
            all_vep_responses.extend(res_data)

    # Process VEP responses
    candidate_annotations = {}
    all_transcript_records = []

    for item in all_vep_responses:
        input_str = item.get("input")
        cand_id = cand_id_map.get(input_str)
        most_severe = item.get("most_severe_consequence")
        v_id = item.get("id")

        tx_consequences = item.get("transcript_consequences", [])
        
        # Collect all transcripts
        for tx in tx_consequences:
            all_transcript_records.append({
                "candidate_id": cand_id,
                "input_string": input_str,
                "transcript_id": tx.get("transcript_id"),
                "gene_symbol": tx.get("gene_symbol"),
                "gene_id": tx.get("gene_id"),
                "biotype": tx.get("biotype"),
                "consequence_terms": ",".join(tx.get("consequence_terms", [])),
                "impact": tx.get("impact"),
                "hgvsc": tx.get("hgvsc"),
                "hgvsp": tx.get("hgvsp"),
                "is_canonical": bool(tx.get("is_canonical")),
                "is_mane_select": bool(tx.get("mane_select")),
                "polyphen_score": tx.get("polyphen_score"),
                "sift_score": tx.get("sift_score")
            })

        best_tx, tx_reason = select_best_transcript(tx_consequences)

        # Extract gnomAD AF from colocated variants
        gnomad_af = None
        gnomad_source = "VEP_colocated"
        colocated = item.get("colocated_variants", [])
        for cv in colocated:
            frequencies = cv.get("frequencies", {})
            for allele, pop_freqs in frequencies.items():
                if "gnomad" in pop_freqs:
                    gnomad_af = float(pop_freqs["gnomad"])
                    break
                elif "gnomad_e" in pop_freqs:
                    gnomad_af = float(pop_freqs["gnomad_e"])
                    break

        candidate_annotations[cand_id] = {
            "ensembl_variant_id": v_id,
            "vep_gene_symbol": best_tx.get("gene_symbol") if best_tx else None,
            "vep_most_severe_consequence": most_severe,
            "vep_consequence_terms": ",".join(best_tx.get("consequence_terms", [])) if best_tx and isinstance(best_tx.get("consequence_terms"), list) else (best_tx.get("consequence_terms") if best_tx else None),
            "vep_selected_transcript_id": best_tx.get("transcript_id") if best_tx else None,
            "vep_selected_transcript_reason": tx_reason,
            "vep_all_transcript_count": len(tx_consequences),
            "vep_hgvsc": best_tx.get("hgvsc") if best_tx else None,
            "vep_hgvsp": best_tx.get("hgvsp") if best_tx else None,
            "vep_impact": best_tx.get("impact") if best_tx else None,
            "gnomad_af": gnomad_af,
            "gnomad_af_source": gnomad_source if gnomad_af is not None else "unavailable",
            "cadd_phred": None,
            "cadd_source": "unavailable",
            "alphamissense_score": None,
            "alphamissense_source": "unavailable",
            "spliceai_max_ds": None,
            "spliceai_source": "unavailable",
            "annotation_status": "success",
            "annotation_failure_reason": None
        }

    # Merge annotations into df_candidates
    ann_rows = []
    for idx, row in df_candidates.iterrows():
        r_chrom = str(row["chrom"]).replace("chr", "")
        r_pos = int(row["pos"])
        r_ref = str(row["ref"])
        r_alt = str(row["alt"])
        cand_id = row.get("candidate_id") if ("candidate_id" in row and pd.notna(row.get("candidate_id"))) else f"{r_chrom}_{r_pos}_{r_ref}_{r_alt}"
        ann = candidate_annotations.get(cand_id, {
            "annotation_status": "unannotated",
            "annotation_failure_reason": "not_returned_by_vep"
        })
        row_dict = row.to_dict()
        row_dict["panel_region_gene"] = row_dict.get("gene_symbol")
        merged_row = {**row_dict, **ann}
        ann_rows.append(merged_row)

    df_annotated = pd.DataFrame(ann_rows)
    df_transcripts = pd.DataFrame(all_transcript_records)

    # Save output Parquet files
    parquet_out.parent.mkdir(parents=True, exist_ok=True)
    df_annotated.to_parquet(parquet_out, index=False)

    csv_out.parent.mkdir(parents=True, exist_ok=True)
    df_annotated.to_csv(csv_out, index=False)

    tx_parquet_out.parent.mkdir(parents=True, exist_ok=True)
    df_transcripts.to_parquet(tx_parquet_out, index=False)

    manifest_payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "total_candidates_annotated": len(df_annotated),
        "total_transcripts_parsed": len(df_transcripts),
        "vep_endpoint": VEP_POST_URL,
        "batch_size_used": batch_size,
        "output_parquet": str(parquet_out),
        "output_csv": str(csv_out),
        "output_transcripts_parquet": str(tx_parquet_out)
    }

    manifest_out.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_out, "w") as f:
        json.dump(manifest_payload, f, indent=2)

    print(f"Annotation complete! Annotated {len(df_annotated)} candidates across {len(df_transcripts)} transcript records.")
    print(f"Annotated Parquet: {parquet_out}")
    print(f"Transcripts Parquet: {tx_parquet_out}")

if __name__ == "__main__":
    main()
