import argparse
import hashlib
import json
import subprocess
import sys
import requests
from datetime import datetime, timezone
from pathlib import Path

ENSEMBL_SEQ_URL = "https://rest.ensembl.org/sequence/region/homo_sapiens"

def fetch_region_fasta(chrom, start, end):
    url = f"{ENSEMBL_SEQ_URL}/{chrom}:{start}..{end}:1?content-type=application/json"
    headers = {"Accept": "application/json"}
    response = requests.get(url, headers=headers, timeout=60)
    response.raise_for_status()
    data = response.json()
    return data.get("seq", "")

def merge_intervals(regions):
    by_contig = {}
    for r in regions:
        c = str(r["contig"])
        if c not in by_contig:
            by_contig[c] = []
        by_contig[c].append(r)

    slices = []
    slice_counter = 1

    for c in sorted(by_contig.keys()):
        regs = sorted(by_contig[c], key=lambda x: x["search_start"])
        current_start = None
        current_end = None
        current_genes = []

        for r in regs:
            s_start = r["search_start"]
            s_end = r["search_end"]
            gene = r["gene_symbol"]

            if current_start is None:
                current_start = s_start
                current_end = s_end
                current_genes = [gene]
            elif s_start <= current_end:  # Overlaps or touches
                current_end = max(current_end, s_end)
                if gene not in current_genes:
                    current_genes.append(gene)
            else:
                # Flush previous slice
                mini_id = f"MVA_{slice_counter:03d}"
                slices.append({
                    "mini_contig_id": mini_id,
                    "chrom": c,
                    "search_start": current_start,
                    "search_end": current_end,
                    "genes": current_genes,
                    "length_bp": current_end - current_start + 1
                })
                slice_counter += 1
                current_start = s_start
                current_end = s_end
                current_genes = [gene]

        if current_start is not None:
            mini_id = f"MVA_{slice_counter:03d}"
            slices.append({
                "mini_contig_id": mini_id,
                "chrom": c,
                "search_start": current_start,
                "search_end": current_end,
                "genes": current_genes,
                "length_bp": current_end - current_start + 1
            })
            slice_counter += 1

    return slices

def main():
    parser = argparse.ArgumentParser(description="Merged Mini-Reference FASTA Generation")
    parser.add_argument("--panel", required=True, help="Path to output/target_panel_regions.json")
    parser.add_argument("--work-dir", default="work", help="Path to work directory")
    args = parser.parse_args()

    panel_path = Path(args.panel)
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    if not panel_path.exists():
        print(f"ERROR: Panel file not found: {panel_path}", file=sys.stderr)
        sys.exit(1)

    with open(panel_path, "r") as f:
        panel_data = json.load(f)

    regions = panel_data.get("regions", [])
    merged_slices = merge_intervals(regions)
    print(f"Merged {len(regions)} panel regions into {len(merged_slices)} non-overlapping mini-reference slices.")

    fasta_path = work_dir / "mini_panel_ref.fasta"
    fasta_lines = []

    for s in merged_slices:
        mini_id = s["mini_contig_id"]
        chrom = s["chrom"]
        start = s["search_start"]
        end = s["search_end"]
        genes_str = ",".join(s["genes"])
        print(f"Fetching FASTA sequence for {mini_id} ({chrom}:{start}-{end}, genes: {genes_str})...")
        seq = fetch_region_fasta(chrom, start, end)

        if len(seq) != s["length_bp"]:
            print(f"WARNING: Sequence length {len(seq)} != expected {s['length_bp']}", file=sys.stderr)

        header = f">{mini_id} GRCh38|{chrom}:{start}-{end}|{genes_str}"
        fasta_lines.append(header)
        for i in range(0, len(seq), 80):
            fasta_lines.append(seq[i:i+80])

    fasta_content = "\n".join(fasta_lines) + "\n"
    with open(fasta_path, "w") as f:
        f.write(fasta_content)

    fasta_sha256 = hashlib.sha256(fasta_content.encode("utf-8")).hexdigest()
    print(f"Wrote mini-reference FASTA to {fasta_path} (SHA256: {fasta_sha256[:10]}...)")

    print("Generating samtools faidx...")
    subprocess.run(["samtools", "faidx", str(fasta_path)], check=True)

    print("Generating samtools dict...")
    dict_path = work_dir / "mini_panel_ref.dict"
    subprocess.run(["samtools", "dict", "-o", str(dict_path), str(fasta_path)], check=True)

    print("Generating minimap2 index...")
    mmi_path = work_dir / "mini_panel_ref.mmi"
    subprocess.run(["minimap2", "-d", str(mmi_path), str(fasta_path)], check=True)

    gene_to_mini = {}
    for s in merged_slices:
        for g in s["genes"]:
            gene_to_mini[g] = s

    manifest_payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "reference_build": panel_data.get("reference_build", "GRCh38"),
        "panel_sha256": panel_data.get("panel_sha256"),
        "fasta_path": str(fasta_path),
        "fasta_sha256": fasta_sha256,
        "dict_path": str(dict_path),
        "mmi_path": str(mmi_path),
        "total_slices": len(merged_slices),
        "slices": merged_slices,
        "gene_to_mini_contig": gene_to_mini,
        "translation_formula": "mini_pos = genomic_pos - search_start + 1"
    }

    manifest_path = work_dir / "mini_panel_ref_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest_payload, f, indent=2)

    print(f"Mini-reference manifest saved to {manifest_path}")

if __name__ == "__main__":
    main()
