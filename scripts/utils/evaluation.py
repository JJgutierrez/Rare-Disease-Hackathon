"""
Track 1 (Variant Prediction) automated scoring.

Adapted from the CAGI6 Rare Genomes Project Challenge assessment methodology
(Stenton et al., 2024), simplified for a single proband ("N-of-1") with a
clinically validated compound-heterozygous answer key.

Submission format (CSV):
    proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr

    - One row per proposed causal variant (or compound-het variant pair).
    - chrom_2/pos_2/ref_2/alt_2 are left blank for single-variant proposals.
    - epcr = Estimated Probability of Causal Relationship, in (0, 1].
    - Up to 10 rows per proband. Rows do NOT need to be pre-sorted by the
      submitter; this script sorts by epcr descending before scoring.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

Variant = tuple[str, int, str, str]  # (chrom, pos, ref, alt)

RANK_POINT_TIERS = [
    (1, 100),
    (3, 50),
    (5, 25),
    (10, 10),
]


@dataclass
class SubmissionRow:
    variants: frozenset[Variant]  # 1 element = single variant, 2 = compound het
    epcr: float
    rank: int  # 1-indexed position after sorting by epcr descending
    finding_type: str = "primary"


@dataclass
class ScoreResult:
    proband_id: str
    full_match_rank: int | None
    partial_match_rank: int | None
    rank_points: float
    f_max: float
    f_max_threshold: float | None
    n_predictions_at_f_max: int


def _parse_variant(chrom: str, pos: str, ref: str, alt: str) -> Variant | None:
    if not chrom or not pos or pd_isna(chrom) or pd_isna(pos):
        return None
    return (str(chrom).strip(), int(pos), str(ref).strip().upper(), str(alt).strip().upper())

def pd_isna(val):
    return val is None or str(val).strip() in ("", "nan", "None", "null")


def load_submission(csv_path: str) -> dict[str, list[SubmissionRow]]:
    """Load a Track 1 submission CSV and group rows by proband, sorted by EPCR desc."""
    by_proband: dict[str, list[tuple[frozenset[Variant], float, str]]] = {}

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pid = row["proband_id"].strip()
            v1 = _parse_variant(row["chrom_1"], row["pos_1"], row["ref_1"], row["alt_1"])
            v2 = _parse_variant(row.get("chrom_2", ""), row.get("pos_2", ""),
                                 row.get("ref_2", ""), row.get("alt_2", ""))
            if v1 is None:
                raise ValueError(f"Row missing primary variant for proband {pid}")
            variants = frozenset([v1, v2]) if v2 else frozenset([v1])
            epcr = float(row["epcr"])
            if not (0 < epcr <= 1):
                raise ValueError(f"EPCR {epcr} out of range (0,1] for proband {pid}")
            finding_type = (row.get("finding_type") or "primary").strip().lower()
            if finding_type not in ("primary", "secondary"):
                raise ValueError(
                    f"finding_type must be 'primary' or 'secondary' for proband {pid}, "
                    f"got '{finding_type}'"
                )
            by_proband.setdefault(pid, []).append((variants, epcr, finding_type))

    result: dict[str, list[SubmissionRow]] = {}
    for pid, rows in by_proband.items():
        if len(rows) > 10:
            raise ValueError(f"Proband {pid} has {len(rows)} rows; max is 10")
        rows_sorted = sorted(enumerate(rows), key=lambda x: (-x[1][1], x[0]))
        result[pid] = [
            SubmissionRow(variants=v, epcr=e, rank=i + 1, finding_type=ft)
            for i, (_, (v, e, ft)) in enumerate(rows_sorted)
        ]
    return result


def _rank_to_points(rank: int) -> int:
    for max_rank, points in RANK_POINT_TIERS:
        if rank <= max_rank:
            return points
    return 0


def score_proband(
    proband_id: str,
    submission_rows: list[SubmissionRow],
    true_variants: frozenset[Variant],
) -> ScoreResult:
    is_compound_het = len(true_variants) == 2

    full_match_rank = None
    partial_match_rank = None

    for row in submission_rows:
        if row.variants == true_variants:
            full_match_rank = row.rank
            break

    if is_compound_het and full_match_rank is None:
        for row in submission_rows:
            if row.variants & true_variants:
                partial_match_rank = row.rank
                break

    if full_match_rank is not None:
        rank_points = float(_rank_to_points(full_match_rank))
    elif partial_match_rank is not None:
        rank_points = 0.5 * _rank_to_points(partial_match_rank)
    else:
        rank_points = 0.0

    thresholds = sorted({row.epcr for row in submission_rows}, reverse=True)
    best_f = 0.0
    best_threshold = None
    best_n = 0

    for t in thresholds:
        predicted_variants: set[Variant] = set()
        n_rows_at_threshold = 0
        for row in submission_rows:
            if row.epcr >= t:
                predicted_variants |= row.variants
                n_rows_at_threshold += 1

        tp = len(predicted_variants & true_variants)
        fp = len(predicted_variants - true_variants)
        fn = len(true_variants - predicted_variants)

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0

        if f > best_f:
            best_f = f
            best_threshold = t
            best_n = n_rows_at_threshold

    return ScoreResult(
        proband_id=proband_id,
        full_match_rank=full_match_rank,
        partial_match_rank=partial_match_rank,
        rank_points=rank_points,
        f_max=best_f,
        f_max_threshold=best_threshold,
        n_predictions_at_f_max=best_n,
    )

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python evaluation.py <submission_csv>", file=sys.stderr)
        sys.exit(1)

    csv_path = sys.argv[1]
    try:
        loaded = load_submission(csv_path)
        print(f"VALIDATION SUCCESS: Loaded {len(loaded)} proband submission(s) from {csv_path}")
        for pid, rows in loaded.items():
            print(f"Proband {pid}: {len(rows)} valid row(s) loaded.")
    except Exception as e:
        print(f"VALIDATION FAILURE: {e}", file=sys.stderr)
        sys.exit(2)
