import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
import docx

def parse_docx_content(doc_path):
    doc = docx.Document(doc_path)
    content_blocks = []
    
    # Parse paragraphs
    for i, p in enumerate(doc.paragraphs):
        text = p.text.strip()
        if text:
            content_blocks.append({"type": "paragraph", "index": i, "text": text})
            
    # Parse tables
    for t_idx, table in enumerate(doc.tables):
        for r_idx, row in enumerate(table.rows):
            row_text = [cell.text.strip() for cell in row.cells if cell.text.strip()]
            if row_text:
                content_blocks.append({
                    "type": "table_row",
                    "table_index": t_idx,
                    "row_index": r_idx,
                    "text": " | ".join(row_text)
                })
                
    return content_blocks

def main():
    parser = argparse.ArgumentParser(description="Clinical Phenotype Parsing & Proband Resolution")
    parser.add_argument("--input", required=True, help="Path to Challenge_Clinical_Phenotype_1.docx")
    parser.add_argument("--preflight", required=True, help="Path to output/vcf_preflight.json")
    parser.add_argument("--output-clinical", required=True, help="Path to write clinical_summary.json")
    parser.add_argument("--output-resolution", required=True, help="Path to write proband_resolution.json")
    args = parser.parse_args()

    input_path = Path(args.input)
    preflight_path = Path(args.preflight)
    output_clinical_path = Path(args.output_clinical)
    output_res_path = Path(args.output_resolution)

    if not input_path.exists():
        print(f"ERROR: Input clinical doc not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    if not preflight_path.exists():
        print(f"ERROR: Preflight file not found: {preflight_path}", file=sys.stderr)
        sys.exit(1)

    with open(preflight_path, "r") as f:
        preflight = json.load(f)

    vcf_samples = preflight.get("sample_ids", [])

    content_blocks = parse_docx_content(input_path)
    full_text = "\n".join([b["text"] for b in content_blocks])

    # Extract Proband ID (looking for patterns like EX2312012 or WGS_EX2312012)
    proband_id_match = re.search(r"\b(EX\d+|WGS_EX\d+|\w*EX\d+\w*)\b", full_text, re.IGNORECASE)
    extracted_proband_id = proband_id_match.group(1) if proband_id_match else None

    # Clinical features search
    features = {
        "microcephaly": bool(re.search(r"microcephaly", full_text, re.IGNORECASE)),
        "growth_retardation": bool(re.search(r"growth|retardation|IUGR", full_text, re.IGNORECASE)),
        "wilms_tumor": bool(re.search(r"wilms|tumor|neoplasm", full_text, re.IGNORECASE)),
        "aneuploidy_mosaicism": bool(re.search(r"aneuploidy|mosaic|karyotype|trisomy", full_text, re.IGNORECASE))
    }

    clinical_summary = {
        "source_document": str(input_path),
        "extracted_at_utc": datetime.now(timezone.utc).isoformat(),
        "extracted_proband_id": extracted_proband_id,
        "blocks_parsed": len(content_blocks),
        "clinical_features": features,
        "content_blocks": content_blocks
    }

    # Proband Resolution Logic
    normalized_vcf_samples = {s.upper().replace("-", "_"): s for s in vcf_samples}
    resolved_sample_id = None
    resolution_method = "none"
    res_status = "fail"

    if extracted_proband_id:
        norm_id = extracted_proband_id.upper().replace("-", "_")
        if norm_id in normalized_vcf_samples:
            resolved_sample_id = normalized_vcf_samples[norm_id]
            resolution_method = "exact_normalized_match"
            res_status = "pass"
        else:
            for n_sample, orig_sample in normalized_vcf_samples.items():
                if norm_id in n_sample or n_sample in norm_id:
                    resolved_sample_id = orig_sample
                    resolution_method = "substring_normalized_match"
                    res_status = "pass"
                    break

    if res_status != "pass" and len(vcf_samples) == 1:
        resolved_sample_id = vcf_samples[0]
        resolution_method = "single_vcf_sample_fallback"
        res_status = "pass"

    proband_resolution = {
        "status": res_status,
        "clinical_proband_id": extracted_proband_id,
        "vcf_sample_id": resolved_sample_id,
        "resolution_method": resolution_method,
        "normalization_rules": [
            "trim_whitespace",
            "uppercase",
            "replace_hyphen_with_underscore_for_comparison"
        ],
        "candidate_sample_ids": vcf_samples,
        "resolved_at_utc": datetime.now(timezone.utc).isoformat()
    }

    output_clinical_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_clinical_path, "w") as f:
        json.dump(clinical_summary, f, indent=2)

    output_res_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_res_path, "w") as f:
        json.dump(proband_resolution, f, indent=2)

    print(f"Clinical parsing completed. Extracted Proband ID: {extracted_proband_id}")
    print(f"Proband Resolution Status: [{res_status.upper()}], Resolved VCF Sample ID: {resolved_sample_id} ({resolution_method})")

if __name__ == "__main__":
    main()
