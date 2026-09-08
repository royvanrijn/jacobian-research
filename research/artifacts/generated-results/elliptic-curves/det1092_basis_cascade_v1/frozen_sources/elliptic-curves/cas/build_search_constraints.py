#!/usr/bin/env python3
"""Import theorem-bound search constraints; --check is read-only and cheap."""

import argparse
import json
from pathlib import Path

from research_runtime.pruning import known_constraints
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--output",type=Path,default=ROOT/"elliptic-curves/data/search_constraints_v1.json")
    args=parser.parse_args()
    payload=known_constraints(ROOT)
    if args.check:
        if json.loads(args.output.read_text())!=payload:raise SystemExit("stale search constraints; rebuild from the authority")
    else:checkpoint(args.output,payload)
    print(f"SEARCH_CONSTRAINTS|count={len(payload['constraints'])}|status=PASS")


if __name__=="__main__":main()
