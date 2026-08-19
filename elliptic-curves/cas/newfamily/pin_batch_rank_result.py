#!/usr/bin/env python3
"""Extract one completed specialization from a batch verifier JSON.

The output is intended for artifacts/generated-results/elliptic-curves/.
Raw search data remains under artifacts/local/elliptic-curves/newfamily/.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input-json", required=True)
    p.add_argument("--parameter", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    rows = json.loads(Path(args.input_json).read_text())
    matches = [r for r in rows if r.get("parameter") == args.parameter]
    if len(matches) != 1:
        raise SystemExit(f"expected one record for {args.parameter}, found {len(matches)}")

    row = matches[0]
    if row.get("status") != "completed":
        raise SystemExit(f"record for {args.parameter} is not completed")

    certificate = {
        "schema": "newfamily_rank_lower_bound_batch_v1",
        "family": "newfamily six-root quartic",
        "roots": [-47, -43, -31, 30, 45, 46],
        "source_batch": args.input_json,
        "parameter": row.get("parameter"),
        "numerator": row.get("numerator"),
        "denominator": row.get("denominator"),
        "known_subgroup_rank": row.get("known_subgroup_rank"),
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
        "claim_boundary": "exact subgroup rank lower bound only; no full-rank upper bound",
    }

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(
        f"PINNED T={certificate['parameter']} "
        f"rank>={certificate['processed_subgroup_rank']} -> {out}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
