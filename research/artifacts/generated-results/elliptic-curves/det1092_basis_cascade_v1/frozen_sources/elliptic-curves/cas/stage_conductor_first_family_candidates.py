#!/usr/bin/env python3
"""Validate staged Mestre/Fermigier candidates and emit strict work queues."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from conductor_first_pipeline import pareto_frontier, rank_first_order, work_queues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.input.read_text())
    records = source["candidates"]
    payload = {
        "schema": "elliptic-curves.conductor-first-family-queues.v1",
        "source": str(args.input),
        "queues": work_queues(records),
        "rank_first_order": [record["id"] for record in rank_first_order(records)],
        "pareto_frontier": [record["id"] for record in pareto_frontier(records)],
        "claim_boundary": (
            "queues enforce computation order; only point_recovery status=certified "
            "enters the exact rank/conductor frontier"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"WROTE {args.output}")
    for stage, ids in payload["queues"].items():
        print(f"{stage}: {len(ids)}")


if __name__ == "__main__":
    main()
