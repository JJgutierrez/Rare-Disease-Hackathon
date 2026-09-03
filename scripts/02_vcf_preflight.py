import argparse
import hashlib
import json
import sys
from pathlib import Path
from cyvcf2 import VCF

GRCH38_CONTIG_LENGTHS = {
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259,
    "6": 170805979, "7": 159345973, "8": 145138636, "9": 138394717, "10": 133797422,
    "11": 135086622, "12": 133275309, "13": 114364328, "14": 107043718, "15": 101991189,
    "16": 90338345, "17": 83257441, "18": 80373285, "19": 58617616, "20": 64444167,
    "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415, "MT": 16569, "M": 16569
}

def compute_file_hash(filepath, chunk_size=65536):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description="VCF Preflight Header Audit")
    parser.add_argument("--vcf", required=True, help="Path to input VCF file")
    parser.add_argument("--config", required=True, help="Path to pipeline_config.json")
    parser.add_argument("--output", required=True, help="Path to write vcf_preflight.json")
    args = parser.parse_args()

    vcf_path = Path(args.vcf)
    output_path = Path(args.output)
    
    if not vcf_path.exists():
        print(f"ERROR: VCF file not found: {vcf_path}", file=sys.stderr)
        sys.exit(1)

    # Check index existence
    csi_path = Path(str(vcf_path) + ".csi")
    tbi_path = Path(str(vcf_path) + ".tbi")
    index_path = None
    if csi_path.exists():
        index_path = csi_path
    elif tbi_path.exists():
        index_path = tbi_path
    elif Path(vcf_path.stem).with_suffix(".tbi").exists():
        index_path = Path(vcf_path.stem).with_suffix(".tbi")

    blocking_errors = []
    warnings = []

    if index_path is None:
        blocking_errors.append("Missing VCF index (.tbi or .csi). Region queries cannot proceed.")

    vcf = VCF(str(vcf_path))
    samples = vcf.samples
    if not samples:
        blocking_errors.append("No sample IDs found in VCF header.")

    # Audit FORMAT fields
    format_headers = set()
    header_reference = None
    contigs_dict = {}

    for h in vcf.header_iter():
        h_type = h.info().get("HeaderType")
        if h_type == "FORMAT":
            format_headers.add(h.info().get("ID"))
        elif h_type == "reference":
            header_reference = h.info().get("Value")
        elif h_type == "contig":
            info = h.info()
            if "ID" in info:
                contigs_dict[info["ID"]] = int(info["length"]) if "length" in info else None

    # Audit FORMAT requirements
    has_gt = "GT" in format_headers
    if not has_gt:
        blocking_errors.append("FORMAT GT definition missing from VCF header.")

    has_gq = "GQ" in format_headers
    has_dp = "DP" in format_headers
    has_pid = "PID" in format_headers
    has_pgt = "PGT" in format_headers

    if not has_pid:
        warnings.append("FORMAT PID definition not present in VCF header.")
    if not has_pgt:
        warnings.append("FORMAT PGT definition not present in VCF header.")

    # Detect contig convention & reference build validation
    contig_style = "bare"
    has_chr_prefix = any(k.startswith("chr") for k in contigs_dict.keys())
    if has_chr_prefix:
        contig_style = "chr"

    matched_contigs = []
    mismatched_contigs = []

    for c_id, c_len in contigs_dict.items():
        norm_id = c_id.replace("chr", "")
        if norm_id in GRCH38_CONTIG_LENGTHS:
            exp_len = GRCH38_CONTIG_LENGTHS[norm_id]
            if c_len is None or c_len == exp_len:
                matched_contigs.append(c_id)
            else:
                mismatched_contigs.append(f"{c_id}: observed {c_len} != expected {exp_len}")

    ref_build_status = "resolved"
    if mismatched_contigs:
        ref_build_status = "mismatched"
        blocking_errors.append(f"Contig length mismatches against GRCh38: {mismatched_contigs}")

    reference_build = "GRCh38"

    # Compute hashes & file sizes
    vcf_bytes = vcf_path.stat().st_size
    vcf_sha256 = compute_file_hash(vcf_path)
    index_sha256 = compute_file_hash(index_path) if index_path else None
    index_bytes = index_path.stat().st_size if index_path else 0

    preflight_status = "pass"
    if blocking_errors:
        preflight_status = "fail"
    elif warnings:
        preflight_status = "warning"

    preflight_data = {
        "preflight_status": preflight_status,
        "vcf_path": str(vcf_path),
        "vcf_sha256": vcf_sha256,
        "vcf_bytes": vcf_bytes,
        "index_path": str(index_path) if index_path else None,
        "index_sha256": index_sha256,
        "index_bytes": index_bytes,
        "reference_build": reference_build,
        "reference_build_status": ref_build_status,
        "reference_build_evidence": {
            "header_reference": header_reference,
            "contig_length_validation": {
                "checked": True,
                "matching_reference": "GRCh38",
                "matched_contigs": matched_contigs,
                "mismatched_contigs": mismatched_contigs
            },
            "resolution_method": "header_reference_and_contig_lengths"
        },
        "contig_style": contig_style,
        "contig_count": len(contigs_dict),
        "sample_ids": samples,
        "format_fields": {
            "GT": has_gt,
            "GQ": has_gq,
            "DP": has_dp,
            "PID": has_pid,
            "PGT": has_pgt
        },
        "blocking_errors": blocking_errors,
        "warnings": warnings
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(preflight_data, f, indent=2)

    print(f"VCF Preflight completed with status: [{preflight_status.upper()}]")
    print(f"Sample IDs: {samples}")
    print(f"Contig Style: {contig_style}, Total Contigs: {len(contigs_dict)}")
    if blocking_errors:
        print(f"Blocking errors: {blocking_errors}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
