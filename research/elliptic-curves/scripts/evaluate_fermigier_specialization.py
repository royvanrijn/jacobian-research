#!/usr/bin/env python3
"""Evaluate Fermigier's baseline sections and optionally search its quartic."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.fermigier_rank import (  # noqa: E402
    certify_fermigier_rank_sections,
    parse_ratpoints_output,
    search_quartic_points_with_gp,
    search_quartic_points_with_ratpoints,
    section_and_point_cloud_differences,
    specialize_fermigier_rank_sections,
    write_json_exclusively,
)
from ecsearch.rank_certification import select_independent_subset  # noqa: E402


def rational_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def parse_rational(value: str) -> Fraction:
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("adapter_parameter", type=parse_rational)
    parser.add_argument("--quartic-height", type=int)
    parser.add_argument(
        "--ratpoints-output",
        type=Path,
        help="replay an existing quiet, abscissa-only ratpoints output file",
    )
    parser.add_argument(
        "--search-engine", choices=("gp", "ratpoints"), default="gp"
    )
    parser.add_argument("--denominator-bound", type=int)
    parser.add_argument("--ratpoints-executable", default="ratpoints")
    parser.add_argument(
        "--certify-searched-subset",
        action="store_true",
        help="select a mod-ell independent subset of baseline plus searched points",
    )
    parser.add_argument("--relation-prime", type=int, default=5)
    parser.add_argument("--maximum-reduction-prime", type=int, default=2000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    specialization = specialize_fermigier_rank_sections(args.adapter_parameter)
    certificate = certify_fermigier_rank_sections(
        specialization,
        maximum_reduction_prime=args.maximum_reduction_prime,
    )
    result: dict[str, object] = {
        "schema": "elliptic-curves.fermigier-specialization-evaluation.v1",
        "adapter_parameter": rational_text(specialization.adapter_parameter),
        "literal_shift": rational_text(specialization.quartic_model.shift),
        "baseline": {
            "quartic_point_count": len(specialization.quartic_points),
            "independent_section_differences": len(
                specialization.section_differences
            ),
            "certificate": certificate.to_json_object(),
        },
        "limitations": [
            "the certificate is a rank lower bound, not an upper bound or saturation proof",
            "a bounded point search is not complete; only the selected modular subset is certified independent",
        ],
    }
    if args.quartic_height is not None and args.ratpoints_output is not None:
        parser.error("choose --quartic-height or --ratpoints-output, not both")
    if (
        args.certify_searched_subset
        and args.quartic_height is None
        and args.ratpoints_output is None
    ):
        parser.error(
            "--certify-searched-subset requires a bounded search or point file"
        )
    if args.quartic_height is not None or args.ratpoints_output is not None:
        if args.ratpoints_output is not None:
            found = parse_ratpoints_output(
                specialization.quartic_model,
                args.ratpoints_output.read_text(),
            )
            engine = "ratpoints output replay"
        elif args.search_engine == "gp":
            if args.denominator_bound is not None:
                parser.error("--denominator-bound is only supported by ratpoints")
            found = search_quartic_points_with_gp(
                specialization.quartic_model, args.quartic_height
            )
            engine = "PARI/GP hyperellratpoints"
        else:
            found = search_quartic_points_with_ratpoints(
                specialization.quartic_model,
                args.quartic_height,
                denominator_bound=args.denominator_bound,
                executable=args.ratpoints_executable,
            )
            engine = "ratpoints"
        result["bounded_quartic_search"] = {
            "engine": engine,
            "height_bound": args.quartic_height,
            "denominator_bound": args.denominator_bound,
            "replayed_output": (
                str(args.ratpoints_output)
                if args.ratpoints_output is not None
                else None
            ),
            "point_count_including_ordinate_signs": len(found),
            "points": [
                [rational_text(x_coordinate), rational_text(y_coordinate)]
                for x_coordinate, y_coordinate in found
            ],
        }
        if args.certify_searched_subset:
            cloud = section_and_point_cloud_differences(specialization, found)
            indices, cloud_certificate = select_independent_subset(
                specialization.canonical_model,
                cloud,
                relation_prime=args.relation_prime,
                maximum_reduction_prime=args.maximum_reduction_prime,
            )
            result["certified_point_cloud_subset"] = {
                "relation_prime": args.relation_prime,
                "deduplicated_difference_count": len(cloud),
                "selected_indices": list(indices),
                "selected_count": len(indices),
                "all_twelve_baseline_differences_selected": set(range(12)).issubset(
                    indices
                ),
                "certificate": cloud_certificate.to_json_object(),
                "interpretation": (
                    "an exact rank lower bound for this specialization; modular "
                    "nonselection does not prove dependence"
                ),
            }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        write_json_exclusively(args.output, rendered)
        print(f"WROTE {args.output}")


if __name__ == "__main__":
    main()
