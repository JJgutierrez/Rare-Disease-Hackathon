import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from cyvcf2 import VCF, Writer
import pandas as pd
import pysam

def compute_candidate_id(sample_id, chrom, pos, ref, alt, alt_index):
    canonical_str = f"{sample_id}|GRCh38|{chrom}|{pos}|{ref}|{alt}|{alt_index}"
    return hashlib.sha256(canonical_str.encode("utf-8")).hexdigest()

def main():
    parser = argparse.ArgumentParser(description="Coordinate-Aware WhatsHap Phasing & Evidence Interpretation")
    parser.add_argument("--candidates-parquet", default="output/panel_candidates_raw.parquet", help="Path to panel_candidates_raw.parquet")
    parser.add_argument("--proband-resolution", default="output/proband_resolution.json", help="Path to proband_resolution.json")
    parser.add_argument("--mini-manifest", default="work/mini_panel_ref_manifest.json", help="Path to mini_panel_ref_manifest.json")
    parser.add_argument("--fasta", default="work/mini_panel_ref.fasta", help="Path to mini_panel_ref.fasta")
    parser.add_argument("--bam", default="output/mini_panel.bam", help="Path to mini_panel.bam")
    parser.add_argument("--config", default="specs/pipeline_config.json", help="Path to specs/pipeline_config.json")
    parser.add_argument("--output-json", default="output/phasing_evidence.json", help="Path to write phasing_evidence.json")
    parser.add_argument("--manifest", default="output/phasing_manifest.json", help="Path to write phasing_manifest.json")
    args = parser.parse_args()

    parquet_path = Path(args.candidates_parquet)
    proband_res_path = Path(args.proband_resolution)
    manifest_path = Path(args.mini_manifest)
    fasta_path = Path(args.fasta)
    bam_path = Path(args.bam)
    config_path = Path(args.config)

    for p in [parquet_path, proband_res_path, manifest_path, fasta_path, bam_path, config_path]:
        if not p.exists():
            print(f"ERROR: Required file not found: {p}", file=sys.stderr)
            sys.exit(1)

    with open(proband_res_path, "r") as f:
        proband_res = json.load(f)
    sample_id = proband_res.get("vcf_sample_id", "WGS_EX2312012")

    with open(manifest_path, "r") as f:
        mini_manifest = json.load(f)
    gene_to_mini = mini_manifest.get("gene_to_mini_contig", {})

    with open(config_path, "r") as f:
        config = json.load(f)

    buffer_bp = config.get("mini_reference", {}).get("phasing_boundary_buffer_bp", 1000)

    # 1. Sample ID Validation Gate
    bam_file = pysam.AlignmentFile(str(bam_path), "rb")
    bam_samples = set()
    for rg in bam_file.header.get("RG", []):
        if "SM" in rg:
            bam_samples.add(rg["SM"])

    if sample_id not in bam_samples:
        print(f"ERROR: Sample ID '{sample_id}' not found in BAM @RG SM headers: {bam_samples}", file=sys.stderr)
        sys.exit(1)

    df_cand = pd.read_parquet(parquet_path)
    if df_cand.empty:
        print("WARNING: Candidate dataset is empty. Exiting phasing step.", file=sys.stderr)
        sys.exit(0)

    # Compute candidate_id
    df_cand["candidate_id"] = df_cand.apply(
        lambda r: compute_candidate_id(sample_id, r["chrom"], r["pos"], r["ref"], r["alt"], r["alt_index"]), axis=1
    )

    # Fetch mini reference sequences for REF verification
    mini_fasta = pysam.FastaFile(str(fasta_path))

    grch38_vcf_records = []
    mini_ref_vcf_records = []
    exclusions = []

    for idx, row in df_cand.iterrows():
        cand_id = row["candidate_id"]
        chrom = str(row["chrom"])
        pos = int(row["pos"])
        ref = str(row["ref"])
        alt = str(row["alt"])
        gene = str(row["gene_symbol"])
        is_snv = bool(row["is_snv"])

        if not is_snv or len(ref) != 1 or len(alt) != 1:
            exclusions.append({
                "candidate_id": cand_id,
                "source_record_key": row.get("source_record_key"),
                "exclusion_reason": "non_snv_variant",
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt
            })
            continue

        mini_slice = gene_to_mini.get(gene)
        if not mini_slice:
            exclusions.append({
                "candidate_id": cand_id,
                "source_record_key": row.get("source_record_key"),
                "exclusion_reason": "unmapped_gene_slice",
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt
            })
            continue

        mini_contig = mini_slice["mini_contig_id"]
        slice_start = mini_slice["search_start"]
        slice_end = mini_slice["search_end"]

        # Coordinate translation: mini_pos = pos - slice_start + 1
        mini_pos = pos - slice_start + 1

        if mini_pos < 1 or mini_pos > (slice_end - slice_start + 1):
            exclusions.append({
                "candidate_id": cand_id,
                "source_record_key": row.get("source_record_key"),
                "exclusion_reason": "translated_pos_out_of_bounds",
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt
            })
            continue

        # REF allele validation against mini FASTA
        ref_seq_in_fasta = mini_fasta.fetch(mini_contig, mini_pos - 1, mini_pos).upper()
        if ref_seq_in_fasta != ref.upper():
            exclusions.append({
                "candidate_id": cand_id,
                "source_record_key": row.get("source_record_key"),
                "exclusion_reason": f"ref_mismatch_in_fasta: {ref} != {ref_seq_in_fasta}",
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt
            })
            continue

        # Boundary buffer check
        dist_start = mini_pos
        dist_end = (slice_end - slice_start + 1) - mini_pos
        boundary_status = "edge_limited" if (dist_start < buffer_bp or dist_end < buffer_bp) else "eligible"

        gt_val = row["gt"]

        mini_ref_vcf_records.append({
            "candidate_id": cand_id,
            "mini_contig": mini_contig,
            "mini_pos": mini_pos,
            "ref": ref,
            "alt": alt,
            "gt": gt_val,
            "boundary_status": boundary_status,
            "orig_chrom": chrom,
            "orig_pos": pos,
            "gene": gene
        })

    # Save exclusions
    df_exclusions = pd.DataFrame(exclusions)
    excl_path = Path("output/phasing_input_exclusions.parquet")
    excl_path.parent.mkdir(parents=True, exist_ok=True)
    df_exclusions.to_parquet(excl_path, index=False)

    print(f"Prepared {len(mini_ref_vcf_records)} eligible candidates for WhatsHap phasing ({len(exclusions)} excluded).")

    if not mini_ref_vcf_records:
        print("WARNING: No eligible candidates for WhatsHap phasing.", file=sys.stderr)
        sys.exit(0)

    # Sort mini_ref_vcf_records by mini_contig, mini_pos
    df_mini_vcf = pd.DataFrame(mini_ref_vcf_records).sort_values(by=["mini_contig", "mini_pos"])

    # Build mini_ref VCF string and write BGZF
    vcf_header_lines = [
        "##fileformat=VCFv4.2",
        f"##reference=file://{fasta_path.absolute()}",
        f"##SAMPLE=<ID={sample_id}>"
    ]
    for s in mini_manifest.get("slices", []):
        vcf_header_lines.append(f"##contig=<ID={s['mini_contig_id']},length={s['length_bp']}>")

    vcf_header_lines.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    vcf_header_lines.append("##INFO=<ID=CAND_ID,Number=1,Type=String,Description=\"Candidate ID\">")
    vcf_header_lines.append("##INFO=<ID=ORIG_CHROM,Number=1,Type=String,Description=\"Original Chromosome\">")
    vcf_header_lines.append("##INFO=<ID=ORIG_POS,Number=1,Type=Integer,Description=\"Original Position\">")
    vcf_header_lines.append(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_id}")

    mini_vcf_uncomp = Path("work/panel_candidates_raw.mini_ref.vcf")
    with open(mini_vcf_uncomp, "w") as f:
        f.write("\n".join(vcf_header_lines) + "\n")
        for idx, r in df_mini_vcf.iterrows():
            info_field = f"CAND_ID={r['candidate_id']};ORIG_CHROM={r['orig_chrom']};ORIG_POS={r['orig_pos']}"
            line = f"{r['mini_contig']}\t{r['mini_pos']}\t{r['candidate_id']}\t{r['ref']}\t{r['alt']}\t.\tPASS\t{info_field}\tGT\t{r['gt']}"
            f.write(line + "\n")

    mini_vcf_gz = Path("output/panel_candidates_raw.mini_ref.vcf.gz")
    pysam.tabix_compress(str(mini_vcf_uncomp), str(mini_vcf_gz), force=True)
    pysam.tabix_index(str(mini_vcf_gz), preset="vcf", force=True)
    mini_vcf_uncomp.unlink()

    print(f"Created tabix-indexed VCF: {mini_vcf_gz}")

    # Invoke WhatsHap phase
    phased_vcf_gz = Path("output/phased_candidates.mini_ref.vcf.gz")
    whatshap_cmd = [
        sys.executable, "-m", "whatshap", "phase",
        "--reference", str(fasta_path),
        "--output", str(phased_vcf_gz),
        str(mini_vcf_gz),
        str(bam_path)
    ]
    print(f"Executing WhatsHap: {' '.join(whatshap_cmd)}")
    subprocess.run(whatshap_cmd, check=True)

    # Index phased VCF
    pysam.tabix_index(str(phased_vcf_gz), preset="vcf", force=True)
    print(f"Indexed phased VCF: {phased_vcf_gz}.tbi")

    # Read phased VCF and parse phase blocks
    phased_vcf = VCF(str(phased_vcf_gz))
    phased_candidates = {}

    for var in phased_vcf:
        cand_id = var.ID
        gt_arr = var.genotypes[0]
        is_phased = bool(gt_arr[2])
        
        # Check PS tag
        ps_val = None
        try:
            ps_raw = var.format("PS")
            if ps_raw is not None:
                ps_val = str(ps_raw[0][0])
        except Exception:
            pass

        phased_candidates[cand_id] = {
            "is_phased": is_phased,
            "phase_set": ps_val,
            "phased_gt": f"{gt_arr[0]}|{gt_arr[1]}" if is_phased else f"{gt_arr[0]}/{gt_arr[1]}"
        }

    # Evaluate Candidate Pairs for Biallelic Phasing Evidence
    df_mini_indexed = df_mini_vcf.set_index("candidate_id")
    candidate_pairs = []

    # Group candidate alleles by gene
    by_gene = df_mini_vcf.groupby("gene")

    for gene_name, group in by_gene:
        if len(group) < 2:
            continue
        g_list = group.to_dict("records")
        for i in range(len(g_list)):
            for j in range(i + 1, len(g_list)):
                c1 = g_list[i]
                c2 = g_list[j]

                dist = abs(c1["orig_pos"] - c2["orig_pos"])
                pair_id = f"{c1['orig_chrom']}:{c1['orig_pos']}:{c1['ref']}:{c1['alt']}__{c2['orig_chrom']}:{c2['orig_pos']}:{c2['ref']}:{c2['alt']}"

                ph1 = phased_candidates.get(c1["candidate_id"], {})
                ph2 = phased_candidates.get(c2["candidate_id"], {})

                is_p1 = ph1.get("is_phased", False)
                is_p2 = ph2.get("is_phased", False)
                ps1 = ph1.get("phase_set")
                ps2 = ph2.get("phase_set")

                edge_limited = (c1["boundary_status"] == "edge_limited" or c2["boundary_status"] == "edge_limited")

                evidence_state = "phase_ambiguous"
                phase_score = None
                phase_conclusion = "unknown"

                if is_p1 and is_p2 and ps1 and ps2 and ps1 == ps2:
                    gt1 = ph1.get("phased_gt", "")
                    gt2 = ph2.get("phased_gt", "")
                    if gt1 in ["0|1", "1|0"] and gt2 in ["0|1", "1|0"]:
                        if gt1 != gt2:
                            evidence_state = "trans_supported"
                            phase_score = 1.0
                            phase_conclusion = "in_trans"
                        else:
                            evidence_state = "cis_supported"
                            phase_score = 0.0
                            phase_conclusion = "in_cis"
                elif edge_limited:
                    evidence_state = "edge_limited_ambiguous"
                    phase_score = None
                    phase_conclusion = "unknown"

                candidate_pairs.append({
                    "pair_id": pair_id,
                    "gene_symbol": gene_name,
                    "candidate_id_1": c1["candidate_id"],
                    "candidate_id_2": c2["candidate_id"],
                    "chrom": c1["orig_chrom"],
                    "pos_1": c1["orig_pos"],
                    "pos_2": c2["orig_pos"],
                    "distance_bp": dist,
                    "evidence_state": evidence_state,
                    "phase_score": phase_score,
                    "phase_conclusion": phase_conclusion,
                    "phase_set_1": ps1,
                    "phase_set_2": ps2,
                    "boundary_status_1": c1["boundary_status"],
                    "boundary_status_2": c2["boundary_status"],
                    "sample_id": sample_id
                })

    # Save phasing evidence JSON
    output_json_path = Path(args.output_json)
    output_json_path.parent.mkdir(parents=True, exist_ok=True)

    evidence_payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "sample_id": sample_id,
        "total_candidate_pairs_evaluated": len(candidate_pairs),
        "trans_supported_pairs": len([p for p in candidate_pairs if p["evidence_state"] == "trans_supported"]),
        "cis_supported_pairs": len([p for p in candidate_pairs if p["evidence_state"] == "cis_supported"]),
        "candidate_pairs": candidate_pairs
    }

    with open(output_json_path, "w") as f:
        json.dump(evidence_payload, f, indent=2)

    # Manifest
    manifest_out = Path(args.manifest)
    manifest_out.parent.mkdir(parents=True, exist_ok=True)
    manifest_payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "sample_id": sample_id,
        "eligible_candidates_phased": len(mini_ref_vcf_records),
        "exclusions_count": len(exclusions),
        "mini_ref_vcf": str(mini_vcf_gz),
        "phased_mini_ref_vcf": str(phased_vcf_gz),
        "phasing_evidence_json": str(output_json_path),
        "exclusions_parquet": str(excl_path)
    }
    with open(manifest_out, "w") as f:
        json.dump(manifest_payload, f, indent=2)

    print(f"WhatsHap phasing complete. Evaluated {len(candidate_pairs)} candidate pairs.")
    print(f"Trans-supported: {evidence_payload['trans_supported_pairs']}, Cis-supported: {evidence_payload['cis_supported_pairs']}")
    print(f"Phasing evidence saved to {output_json_path}")

if __name__ == "__main__":
    main()
