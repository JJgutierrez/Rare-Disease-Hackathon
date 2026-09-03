import argparse
import hashlib
import json
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path
from cyvcf2 import VCF
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

def determine_variant_class(ref, alt):
    if len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT":
        return "SNV", True, False, False
    elif alt.startswith("<") or alt.endswith(">"):
        return "SYMBOLIC", False, False, True
    else:
        return "INDEL", False, True, False

def main():
    parser = argparse.ArgumentParser(description="Sample-Aware & Multiallelic-Safe VCF Slicing")
    parser.add_argument("--vcf", required=True, help="Path to input VCF")
    parser.add_argument("--preflight", required=True, help="Path to output/vcf_preflight.json")
    parser.add_argument("--proband-resolution", required=True, help="Path to output/proband_resolution.json")
    parser.add_argument("--panel", required=True, help="Path to output/target_panel_regions.json")
    parser.add_argument("--config", required=True, help="Path to specs/pipeline_config.json")
    parser.add_argument("--output-parquet", required=True, help="Path to write panel_candidates_raw.parquet")
    parser.add_argument("--output-csv", required=True, help="Path to write panel_candidates_raw.csv")
    parser.add_argument("--manifest", required=True, help="Path to write run_manifest.json")
    args = parser.parse_args()

    vcf_path = Path(args.vcf)
    preflight_path = Path(args.preflight)
    proband_res_path = Path(args.proband_resolution)
    panel_path = Path(args.panel)
    config_path = Path(args.config)
    parquet_out = Path(args.output_parquet)
    csv_out = Path(args.output_csv)
    manifest_out = Path(args.manifest)

    for p in [vcf_path, preflight_path, proband_res_path, panel_path, config_path]:
        if not p.exists():
            print(f"ERROR: Required file not found: {p}", file=sys.stderr)
            sys.exit(1)

    with open(preflight_path, "r") as f:
        preflight = json.load(f)
    if preflight.get("preflight_status") not in ["pass", "warning"]:
        print(f"ERROR: Preflight status failed: {preflight.get('preflight_status')}", file=sys.stderr)
        sys.exit(1)

    with open(proband_res_path, "r") as f:
        proband_res = json.load(f)
    if proband_res.get("status") != "pass":
        print(f"ERROR: Proband resolution status failed: {proband_res.get('status')}", file=sys.stderr)
        sys.exit(1)

    proband_sample_id = proband_res.get("vcf_sample_id")

    with open(panel_path, "r") as f:
        panel_data = json.load(f)
    regions = panel_data.get("regions", [])

    with open(config_path, "r") as f:
        config = json.load(f)

    min_gq = config.get("min_gq", 20)
    min_dp = config.get("min_dp", 10)
    accept_filter_values = set(config.get("accept_filter_values", ["PASS"]))
    accept_unfiltered_dot = config.get("accept_unfiltered_dot", True)

    vcf = VCF(str(vcf_path))
    samples = vcf.samples
    if proband_sample_id not in samples:
        print(f"ERROR: Proband sample ID '{proband_sample_id}' not in VCF samples: {samples}", file=sys.stderr)
        sys.exit(1)

    proband_idx = samples.index(proband_sample_id)
    run_id = f"run_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}_{uuid.uuid4().hex[:8]}"

    # Trace counters
    counts = {
        "total_records_processed": 0,
        "region_overlap_records": 0,
        "pass_filtered_records": 0,
        "proband_called_records": 0,
        "proband_non_ref_records": 0,
        "gq_passed_or_flagged_records": 0,
        "dp_passed_or_flagged_records": 0,
        "final_candidate_alleles": 0
    }

    candidates = []

    for reg in regions:
        contig = reg["contig"]
        start = reg["search_start"]
        end = reg["search_end"]
        gene_symbol = reg["gene_symbol"]
        tier = reg["tier"]

        query_str = f"{contig}:{start}-{end}"
        try:
            region_iter = vcf(query_str)
        except Exception as e:
            print(f"WARNING: Region query failed for {query_str}: {e}", file=sys.stderr)
            continue

        for variant in region_iter:
            counts["total_records_processed"] += 1
            counts["region_overlap_records"] += 1

            # 1. Filter evaluation
            v_filter = variant.FILTER
            if v_filter is None:
                is_pass = accept_unfiltered_dot
                v_filter_str = "."
            elif v_filter == "PASS" or v_filter in accept_filter_values:
                is_pass = True
                v_filter_str = "PASS"
            else:
                is_pass = False
                v_filter_str = str(v_filter)

            if not is_pass:
                continue

            counts["pass_filtered_records"] += 1

            # 2. Genotype extraction
            gt_arr = variant.genotypes[proband_idx]  # e.g. [0, 1, False]
            allele1, allele2 = gt_arr[0], gt_arr[1]
            is_phased = bool(gt_arr[2])

            if allele1 < 0 and allele2 < 0:
                continue  # Uncalled

            counts["proband_called_records"] += 1

            # Check non-reference
            proband_allele_indices = [idx for idx in [allele1, allele2] if idx > 0]
            if not proband_allele_indices:
                continue  # Hom ref (0/0)

            counts["proband_non_ref_records"] += 1

            # 3. Quality extraction using cyvcf2 vector properties
            gq_val = None
            gq_missing = True
            try:
                gq_raw = float(variant.gt_quals[proband_idx])
                if gq_raw >= 0:
                    gq_val = gq_raw
                    gq_missing = False
            except Exception:
                pass

            dp_val = None
            dp_missing = True
            try:
                ref_dp = int(variant.gt_ref_depths[proband_idx])
                alt_dp = int(variant.gt_alt_depths[proband_idx])
                dp_val = ref_dp + alt_dp
                dp_missing = False
            except Exception:
                pass

            # GQ Policy Evaluation
            gq_policy_applied = "threshold_pass"
            if gq_missing:
                if config["missing_gq_policy"] == "retain_with_flag":
                    gq_policy_applied = "missing_retain_with_flag"
                else:
                    continue
            elif gq_val < min_gq:
                if config["low_gq_policy"] == "exclude":
                    continue
                else:
                    gq_policy_applied = "low_gq_retained"

            counts["gq_passed_or_flagged_records"] += 1

            # DP Policy Evaluation
            dp_policy_applied = "threshold_pass"
            if dp_missing:
                if config["missing_dp_policy"] == "retain_with_flag":
                    dp_policy_applied = "missing_retain_with_flag"
                else:
                    continue
            elif dp_val < min_dp:
                if config["low_dp_policy"] == "exclude":
                    continue
                else:
                    dp_policy_applied = "low_dp_retained"

            counts["dp_passed_or_flagged_records"] += 1

            # PID / PGT extraction
            pid_val = None
            try:
                pid_raw = variant.format("PID")
                if pid_raw is not None and len(pid_raw) > proband_idx:
                    pid_val = str(pid_raw[proband_idx][0])
            except Exception:
                pass

            pgt_val = None
            try:
                pgt_raw = variant.format("PGT")
                if pgt_raw is not None and len(pgt_raw) > proband_idx:
                    pgt_val = str(pgt_raw[proband_idx][0])
            except Exception:
                pass

            ref_base = variant.REF
            alt_alleles = variant.ALT
            pos = variant.POS
            chrom = str(variant.CHROM)
            qual = float(variant.QUAL) if variant.QUAL is not None else None

            # Source record key
            alts_comma = ",".join(alt_alleles)
            source_record_key = f"{chrom}:{pos}:{ref_base}:{alts_comma}"
            gt_str = f"{allele1}/{allele2}" if not is_phased else f"{allele1}|{allele2}"

            # 4. Multiallelic Normalization: emit 1 row per ALT allele represented in proband GT
            unique_proband_alt_indices = sorted(list(set(proband_allele_indices)))

            for alt_idx in unique_proband_alt_indices:
                if alt_idx - 1 < len(alt_alleles):
                    alt_base = alt_alleles[alt_idx - 1]
                else:
                    alt_base = "N"

                var_class, is_snv, is_indel, is_symbolic = determine_variant_class(ref_base, alt_base)

                filter_reasons = [v_filter_str, f"GT_{gt_str}"]
                if gq_missing:
                    filter_reasons.append("GQ_MISSING")
                else:
                    filter_reasons.append(f"GQ_{gq_val}")
                if dp_missing:
                    filter_reasons.append("DP_MISSING")
                else:
                    filter_reasons.append(f"DP_{dp_val}")

                filter_reason_str = ";".join(filter_reasons)

                candidates.append({
                    "chrom": chrom,
                    "pos": int(pos),
                    "ref": ref_base,
                    "alt": alt_base,
                    "alt_index": int(alt_idx),
                    "variant_class": var_class,
                    "is_snv": is_snv,
                    "is_indel": is_indel,
                    "is_symbolic": is_symbolic,
                    "gene_symbol": gene_symbol,
                    "tier": tier,
                    "filter": v_filter_str,
                    "is_filtered_pass": is_pass,
                    "qual": qual,
                    "sample_id": proband_sample_id,
                    "gt": gt_str,
                    "proband_allele_indices": str(proband_allele_indices),
                    "is_non_reference": True,
                    "gq": gq_val,
                    "gq_missing": gq_missing,
                    "dp": dp_val,
                    "dp_missing": dp_missing,
                    "pid": pid_val,
                    "pgt": pgt_val,
                    "is_phased": is_phased,
                    "phase_set": pid_val,
                    "filter_decision": "retained",
                    "filter_reason": filter_reason_str,
                    "gq_policy_applied": gq_policy_applied,
                    "dp_policy_applied": dp_policy_applied,
                    "source_record_key": source_record_key,
                    "source_vcf": str(vcf_path),
                    "vcf_sha256": preflight.get("vcf_sha256"),
                    "panel_sha256": panel_data.get("panel_sha256"),
                    "reference_build": panel_data.get("reference_build", "GRCh38"),
                    "run_id": run_id
                })

                counts["final_candidate_alleles"] += 1

    df_candidates = pd.DataFrame(candidates)
    if not df_candidates.empty:
        # Deduplicate by chrom, pos, ref, alt, sample_id
        df_candidates = df_candidates.drop_duplicates(subset=["chrom", "pos", "ref", "alt", "sample_id"])

    # Output Parquet
    parquet_out.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pandas(df_candidates)
    pq.write_table(table, parquet_out)

    # Output CSV inspection
    csv_out.parent.mkdir(parents=True, exist_ok=True)
    df_candidates.to_csv(csv_out, index=False)

    # Output Manifest
    manifest_payload = {
        "run_id": run_id,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_vcf": str(vcf_path),
        "vcf_sha256": preflight.get("vcf_sha256"),
        "panel_sha256": panel_data.get("panel_sha256"),
        "proband_sample_id": proband_sample_id,
        "counts": counts,
        "output_parquet": str(parquet_out),
        "output_csv": str(csv_out),
        "parquet_candidate_count": len(df_candidates)
    }

    manifest_out.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_out, "w") as f:
        json.dump(manifest_payload, f, indent=2)

    print(f"VCF Slicing complete. Retained candidate alleles: {len(df_candidates)}")
    print(f"Parquet saved to: {parquet_out}")
    print(f"Run manifest saved to: {manifest_out}")

if __name__ == "__main__":
    main()
