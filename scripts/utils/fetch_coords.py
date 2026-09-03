import argparse
import hashlib
import json
import sys
import requests
from datetime import datetime, timezone
from pathlib import Path

ENSEMBL_POST_URL = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"

CORE_MVA_GENES = ["BUB1B", "TRIP13", "CEP57"]
EXTENDED_SAC_GENES = ["MAD2L1", "CDC20", "BUB1", "BUB3", "TTK", "CENPE", "PLK1"]

def fetch_ensembl_batch(gene_symbols):
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }
    payload = {"symbols": gene_symbols}
    response = requests.post(ENSEMBL_POST_URL, headers=headers, json=payload, timeout=30)
    response.raise_for_status()
    return response.json()

def main():
    parser = argparse.ArgumentParser(description="Batch Ensembl Gene Lookup with Full Provenance")
    parser.add_argument("--preflight", required=True, help="Path to output/vcf_preflight.json")
    parser.add_argument("--config", required=True, help="Path to specs/pipeline_config.json")
    parser.add_argument("--output", required=True, help="Path to write output/target_panel_regions.json")
    args = parser.parse_args()

    preflight_path = Path(args.preflight)
    config_path = Path(args.config)
    output_path = Path(args.output)

    if not preflight_path.exists():
        print(f"ERROR: Preflight file not found: {preflight_path}", file=sys.stderr)
        sys.exit(1)

    with open(preflight_path, "r") as f:
        preflight = json.load(f)

    if preflight.get("preflight_status") not in ["pass", "warning"]:
        print(f"ERROR: Blocking preflight status: {preflight.get('preflight_status')}", file=sys.stderr)
        sys.exit(1)

    ref_build = preflight.get("reference_build", "GRCh38")
    contig_style = preflight.get("contig_style", "bare")

    with open(config_path, "r") as f:
        config = json.load(f)
    padding_bp = config.get("panel_padding_bp", 50000)

    all_symbols = CORE_MVA_GENES + EXTENDED_SAC_GENES
    print(f"Performing Ensembl batch lookup for {len(all_symbols)} gene symbols...")
    raw_response = fetch_ensembl_batch(all_symbols)

    # Save cached raw response in cache/ensembl
    cache_dir = Path("cache/ensembl")
    cache_dir.mkdir(parents=True, exist_ok=True)
    timestamp_str = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    cache_file = cache_dir / f"ensembl_batch_{timestamp_str}.json"
    with open(cache_file, "w") as f:
        json.dump(raw_response, f, indent=2)

    regions = []
    for symbol in all_symbols:
        gene_info = raw_response.get(symbol)
        if not gene_info or "seq_region_name" not in gene_info:
            print(f"WARNING: Symbol {symbol} not resolved by Ensembl REST API", file=sys.stderr)
            continue

        raw_seq_name = str(gene_info["seq_region_name"])
        if contig_style == "bare":
            contig = raw_seq_name.replace("chr", "")
        else:
            contig = f"chr{raw_seq_name.replace('chr', '')}"

        gene_start = int(gene_info["start"])
        gene_end = int(gene_info["end"])
        search_start = max(1, gene_start - padding_bp)
        search_end = gene_end + padding_bp
        strand = int(gene_info.get("strand", 1))
        ensembl_id = gene_info.get("id")

        tier = "core_mva" if symbol in CORE_MVA_GENES else "mitotic_checkpoint_extended"

        regions.append({
            "gene_symbol": symbol,
            "tier": tier,
            "ensembl_gene_id": ensembl_id,
            "contig": contig,
            "gene_start": gene_start,
            "gene_end": gene_end,
            "search_start": search_start,
            "search_end": search_end,
            "strand": strand
        })

    result_payload = {
        "reference_build": ref_build,
        "coordinate_source": {
            "provider": "Ensembl REST API",
            "endpoint": ENSEMBL_POST_URL,
            "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            "cache_file": str(cache_file)
        },
        "padding_bp": padding_bp,
        "contig_style": contig_style,
        "total_regions": len(regions),
        "regions": regions
    }

    # Add content hash
    content_str = json.dumps(result_payload, sort_keys=True)
    panel_sha256 = hashlib.sha256(content_str.encode("utf-8")).hexdigest()
    result_payload["panel_sha256"] = panel_sha256

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(result_payload, f, indent=2)

    print(f"Target panel regions written to {output_path} (Total regions: {len(regions)}, Hash: {panel_sha256[:10]}...)")

if __name__ == "__main__":
    main()
