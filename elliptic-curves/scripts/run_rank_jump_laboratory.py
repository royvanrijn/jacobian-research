#!/usr/bin/env python3
"""Replay the provenance-checked rank-jump retrieval laboratory."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
sys.path.insert(0, str(ELLIPTIC_ROOT))

from ecsearch.rank_jump_benchmark import evaluate_lab_manifest  # noqa: E402


DEFAULT_INPUT = ELLIPTIC_ROOT / "data" / "rank_jump_laboratory_v1.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument(
        "--output",
        type=Path,
        help="optional JSON result; local experiments belong under artifacts/local/",
    )
    args = parser.parse_args()

    manifest = json.loads(args.input.read_text(encoding="utf-8"))
    result = evaluate_lab_manifest(manifest, ROOT)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

    summary = result["summary"]
    r17 = next(
        family
        for family in result["families"]
        if family["id"] == "published-r17-compact-t"
    )
    first_run = r17["ranking_runs"][0]
    ranks = [
        observation["population_rank"]
        for observation in first_run["positive_observations"]
    ]
    print(
        "RANK_JUMP_LAB_PASS "
        f"families={summary['family_count']} "
        f"ranked_families={summary['ranked_family_count']} "
        f"certified_positives={summary['certified_positive_count']} "
        f"ranked_positives={summary['ranked_positive_count']} "
        f"r17_population={first_run['population_count']} "
        f"r17_control_ranks={','.join(map(str, ranks))}"
    )
    if args.output is not None:
        print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
