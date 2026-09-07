#!/usr/bin/env python3
"""Cheap structural/observability audits, not a new rank search.

Never mutates a frozen campaign. 'family' consumes the retained fixed-field
anchor; 'visibility' is an explicitly retrospective oracle; 'mask' writes
separate detector and oracle files. No mode invokes Sage/PARI or enumeration.
"""
import argparse
import json
from pathlib import Path

from fixed_cubic_geometry import generic_geometry, alternating_rank_distribution
from search_observability import point_visibility, masked_control


def load(path):
    return json.loads(path.read_text())


def write(path, result):
    text = json.dumps(result, indent=2, sort_keys=True)+"\n"
    if path is None:
        print(text, end="")
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("x") as stream:  # Do not overwrite a retained experiment.
            stream.write(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)
    family = sub.add_parser("family")
    family.add_argument("--input", type=Path, required=True)
    family.add_argument("--output", type=Path)
    visibility = sub.add_parser("visibility")
    visibility.add_argument("--record", type=Path, required=True)
    visibility.add_argument("--oracle", type=Path, required=True)
    visibility.add_argument("--output", type=Path)
    masking = sub.add_parser("mask")
    masking.add_argument("--input", type=Path, required=True)
    masking.add_argument("--withhold", required=True, help="comma-separated zero-based indices")
    masking.add_argument("--search-input", type=Path, required=True)
    masking.add_argument("--oracle", type=Path, required=True)
    null = sub.add_parser("ct-null")
    null.add_argument("--dimension", type=int, required=True)
    null.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.mode == "family":
        source = load(args.input)
        coefficients = source["anchor"]["base_polynomial_ascending"]
        if len(coefficients) != 4 or coefficients[2:] != ["0", "1"]:
            raise ValueError("expected the pinned depressed monic cubic format")
        result = generic_geometry(coefficients[1], coefficients[0])
        result["anchor_id"] = source["anchor"]["id"]
        write(args.output, result)
    elif args.mode == "visibility":
        record = load(args.record)
        record = record.get("record", record)  # Shared search checkpoint envelope.
        oracle = load(args.oracle)
        points = oracle["withheld_points"] if isinstance(oracle, dict) else oracle
        write(args.output, {"schema": "elliptic-curves.point-observability.v1",
            "retrospective_only": True, "points": [point_visibility(record, p) for p in points],
            "claim_boundary": "Pointwise coverage of known points, not recovery of their entire rational span or rank-jump incidence."})
    elif args.mode == "mask":
        if args.search_input.resolve() == args.oracle.resolve() or args.search_input.exists() or args.oracle.exists():
            raise ValueError("two different new output files required")
        source = load(args.input)
        search, oracle = masked_control(source["curve"], source["points"], source["metric_gram"],
                                        [int(i) for i in args.withhold.split(",")])
        write(args.search_input, search)
        write(args.oracle, oracle)
    else:
        distribution = alternating_rank_distribution(args.dimension)
        write(args.output, {"model": "uniform alternating matrices over F2; illustrative, not an elliptic-curve sampling law",
            "dimension": args.dimension,
            "probabilities_by_rank": {str(k): str(v) for k, v in distribution.items()},
            "maximum_rank_probability": str(distribution[max(distribution)])})


if __name__ == "__main__":
    main()
