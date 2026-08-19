#!/usr/bin/env python3
"""Promote completed newfamily batch-verifier results into compact certificates.

Raw search data remains under ``artifacts/local/elliptic-curves/newfamily/``.
Compact theorem-strength lower-bound records belong under
``artifacts/generated-results/elliptic-curves/``.

Examples::

    # One specialization
    python3 elliptic-curves/cas/newfamily/pin_batch_rank_result.py \
      --input-json artifacts/local/elliptic-curves/newfamily/batch_exact_rank_gain_hits_v3.json \
      --parameter 83/6 \
      --output artifacts/generated-results/elliptic-curves/newfamily_rank14_t83_6_v2.json

    # Entire completed batch
    python3 elliptic-curves/cas/newfamily/pin_batch_rank_result.py \
      --input-json artifacts/local/elliptic-curves/newfamily/batch_exact_rank_gain_hits_v3.json \
      --all \
      --output artifacts/generated-results/elliptic-curves/newfamily_rank_gain_batch_v1.json
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

ROOTS = [-47, -43, -31, 30, 45, 46]


def compact_record(row):
    return {
        "parameter": row.get("parameter"),
        "numerator": row.get("numerator"),
        "denominator": row.get("denominator"),
        "known_subgroup_rank": row.get("known_subgroup_rank"),
        "known_specialization_rank_drop": row.get("known_specialization_rank_drop"),
        "known_specialization_rank_drop_by": row.get("known_specialization_rank_drop_by"),
        "exact_rank_gain_over_known": row.get("exact_rank_gain_over_known"),
        "processed_subgroup_rank": row.get("processed_subgroup_rank"),
        "proved_rank_at_least": row.get("proved_rank_at_least"),
        "root_number": row.get("root_number"),
        "minimal_discriminant_bits": row.get("minimal_discriminant_bits"),
        "point_bits_min": row.get("point_bits_min"),
        "point_bits_median": row.get("point_bits_median"),
        "point_bits_max": row.get("point_bits_max"),
        "unique_hit_count": row.get("unique_hit_count"),
        "independent_hit_count_observed": row.get("independent_hit_count_observed"),
        "hit_rank_increases": row.get("hit_rank_increases"),
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input-json", required=True)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--parameter")
    mode.add_argument("--all", action="store_true")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    rows = json.loads(Path(args.input_json).read_text())
    if not isinstance(rows, list):
        raise SystemExit("input JSON must be a list of batch records")

    completed = [r for r in rows if r.get("status") == "completed"]

    if args.parameter is not None:
        matches = [r for r in completed if r.get("parameter") == args.parameter]
        if len(matches) != 1:
            raise SystemExit(
                f"expected one completed record for {args.parameter}, found {len(matches)}"
            )
        row = matches[0]
        certificate = {
            "schema": "newfamily_rank_lower_bound_batch_v1",
            "family": "newfamily six-root quartic",
            "roots": ROOTS,
            "source_batch": args.input_json,
            **compact_record(row),
            "claim_boundary": "exact subgroup rank lower bound only; no full-rank upper bound",
        }
        message = (
            f"PINNED T={certificate['parameter']} "
            f"rank>={certificate['processed_subgroup_rank']}"
        )
    else:
        if not completed:
            raise SystemExit("input contains no completed specializations")

        records = [compact_record(r) for r in completed]
        records.sort(
            key=lambda r: (
                -(r.get("processed_subgroup_rank") or -1),
                -(r.get("exact_rank_gain_over_known") or -1),
                r.get("parameter") or "",
            )
        )

        rank_distribution = Counter(
            int(r["processed_subgroup_rank"])
            for r in records
            if r.get("processed_subgroup_rank") is not None
        )
        baseline_distribution = Counter(
            int(r["known_subgroup_rank"])
            for r in records
            if r.get("known_subgroup_rank") is not None
        )
        gain_distribution = Counter(
            int(r["exact_rank_gain_over_known"])
            for r in records
            if r.get("exact_rank_gain_over_known") is not None
        )

        maximum = max(
            (r.get("processed_subgroup_rank") or -1 for r in records),
            default=-1,
        )
        strongest = [
            r["parameter"] for r in records
            if r.get("processed_subgroup_rank") == maximum
        ]

        certificate = {
            "schema": "newfamily_rank_gain_batch_v1",
            "family": "newfamily six-root quartic",
            "roots": ROOTS,
            "source_batch": args.input_json,
            "completed_specializations": len(records),
            "maximum_processed_subgroup_rank": maximum,
            "strongest_parameters": strongest,
            "processed_rank_distribution": {
                str(k): v for k, v in sorted(rank_distribution.items())
            },
            "known_baseline_rank_distribution": {
                str(k): v for k, v in sorted(baseline_distribution.items())
            },
            "exact_gain_distribution": {
                str(k): v for k, v in sorted(gain_distribution.items())
            },
            "records": records,
            "claim_boundary": (
                "each completed record is an exact processed-subgroup rank lower bound; "
                "no full-rank upper bound or saturation claim is implied"
            ),
        }
        message = (
            f"PINNED batch completed={len(records)} "
            f"max_rank>={maximum} strongest={','.join(strongest)}"
        )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(f"{message} -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
