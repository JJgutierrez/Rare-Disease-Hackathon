import argparse
import glob
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Dynamic FASTQ Lane Alignment to Mini-BAM")
    parser.add_argument("--raw-dir", default="data/raw", help="Path to raw FASTQ directory")
    parser.add_argument("--proband-resolution", default="output/proband_resolution.json", help="Path to proband_resolution.json")
    parser.add_argument("--mmi", default="work/mini_panel_ref.mmi", help="Path to mini_panel_ref.mmi")
    parser.add_argument("--work-dir", default="work", help="Path to work directory")
    parser.add_argument("--output-bam", default="output/mini_panel.bam", help="Path to write mini_panel.bam")
    parser.add_argument("--manifest", default="output/alignment_manifest.json", help="Path to write alignment_manifest.json")
    args = parser.parse_args()

    raw_dir = Path(args.raw_dir)
    proband_res_path = Path(args.proband_resolution)
    mmi_path = Path(args.mmi)
    work_dir = Path(args.work_dir)
    output_bam_path = Path(args.output_bam)
    manifest_path = Path(args.manifest)

    if not proband_res_path.exists():
        print(f"ERROR: Proband resolution file not found: {proband_res_path}", file=sys.stderr)
        sys.exit(1)

    with open(proband_res_path, "r") as f:
        proband_res = json.load(f)

    sample_id = proband_res.get("vcf_sample_id", "WGS_EX2312012")
    print(f"Aligning FASTQs for resolved sample ID: '{sample_id}'")

    if not mmi_path.exists():
        print(f"ERROR: MMI index file not found: {mmi_path}", file=sys.stderr)
        sys.exit(1)

    # Discover lane pairs
    r1_files = sorted(glob.glob(str(raw_dir / "*_R1_*.fastq.gz")))
    lane_pairs = []

    for r1 in r1_files:
        r2 = r1.replace("_R1_", "_R2_")
        if not Path(r2).exists():
            print(f"WARNING: R2 mate missing for R1: {r1}", file=sys.stderr)
            continue

        lane_match = re.search(r"_(L\d+)_R1_", r1)
        lane_id = lane_match.group(1) if lane_match else f"L{len(lane_pairs)+1:03d}"

        lane_pairs.append({
            "lane_id": lane_id,
            "r1": r1,
            "r2": r2
        })

    if not lane_pairs:
        print(f"ERROR: No valid FASTQ lane pairs found in {raw_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Discovered {len(lane_pairs)} FASTQ lane pairs.")

    lane_bams = []

    for pair in lane_pairs:
        lane_id = pair["lane_id"]
        rg_header = f"@RG\\tID:{sample_id}.{lane_id}\\tSM:{sample_id}\\tLB:{sample_id}\\tPL:ILLUMINA\\tPU:{sample_id}.{lane_id}"
        lane_bam = work_dir / f"lane_{lane_id}.bam"

        align_cmd = f"minimap2 -ax sr -t 8 -R '{rg_header}' '{mmi_path}' '{pair['r1']}' '{pair['r2']}' | samtools sort -@ 4 -o '{lane_bam}'"
        print(f"Aligning lane {lane_id} to {lane_bam}...")
        subprocess.run(align_cmd, shell=True, check=True)
        lane_bams.append(str(lane_bam))

    # Merge lane BAMs
    output_bam_path.parent.mkdir(parents=True, exist_ok=True)
    if len(lane_bams) == 1:
        subprocess.run(["mv", lane_bams[0], str(output_bam_path)], check=True)
    else:
        print(f"Merging {len(lane_bams)} lane BAMs into {output_bam_path}...")
        merge_cmd = ["samtools", "merge", "-f", str(output_bam_path)] + lane_bams
        subprocess.run(merge_cmd, check=True)

    # Index merged BAM
    print("Indexing merged mini-BAM...")
    subprocess.run(["samtools", "index", str(output_bam_path)], check=True)

    # QC checks
    print("Running samtools quickcheck...")
    subprocess.run(["samtools", "quickcheck", "-v", str(output_bam_path)], check=True)

    print("Running samtools flagstat...")
    flagstat_proc = subprocess.run(["samtools", "flagstat", str(output_bam_path)], capture_output=True, text=True, check=True)
    flagstat_out = flagstat_proc.stdout

    print("Running samtools idxstats...")
    idxstats_proc = subprocess.run(["samtools", "idxstats", str(output_bam_path)], capture_output=True, text=True, check=True)
    idxstats_out = idxstats_proc.stdout

    manifest_payload = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "sample_id": sample_id,
        "mmi_path": str(mmi_path),
        "output_bam": str(output_bam_path),
        "output_bai": f"{output_bam_path}.bai",
        "lane_pairs_processed": len(lane_pairs),
        "lane_pairs": lane_pairs,
        "flagstat": flagstat_out,
        "idxstats": idxstats_out
    }

    with open(manifest_path, "w") as f:
        json.dump(manifest_payload, f, indent=2)

    print(f"Mini-BAM alignment complete! Manifest written to {manifest_path}")

if __name__ == "__main__":
    main()
