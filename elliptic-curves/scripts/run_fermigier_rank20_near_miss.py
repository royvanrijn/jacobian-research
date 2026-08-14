#!/usr/bin/env python3
"""Generate the pinned low-conductor Fermigier rank-20 near miss."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.fermigier_near_miss import (  # noqa: E402
    FERMIGIER_RANK20_DENOMINATOR_BOUND,
    FERMIGIER_RANK20_PARAMETER,
    FERMIGIER_RANK20_SEARCH_HEIGHT,
    build_fermigier_rank20_manifest,
)
from ecsearch.fermigier_rank import (  # noqa: E402
    search_quartic_points_with_ratpoints,
    specialize_fermigier_rank_sections,
    write_json_exclusively,
)
from ecsearch.fermigier_near_miss import canonical_ratpoints_output  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ratpoints-output", type=Path)
    parser.add_argument("--ratpoints-executable", default="ratpoints")
    parser.add_argument("--maximum-reduction-prime", type=int, default=2000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.ratpoints_output is None:
        specialization = specialize_fermigier_rank_sections(
            FERMIGIER_RANK20_PARAMETER
        )
        points = search_quartic_points_with_ratpoints(
            specialization.quartic_model,
            FERMIGIER_RANK20_SEARCH_HEIGHT,
            denominator_bound=FERMIGIER_RANK20_DENOMINATOR_BOUND,
            executable=args.ratpoints_executable,
        )
        abscissas = tuple(
            point[0]
            for index, point in enumerate(points)
            if index == 0 or point[0] != points[index - 1][0]
        )
        raw_output = canonical_ratpoints_output(abscissas)
    else:
        raw_output = args.ratpoints_output.read_text()
    manifest = build_fermigier_rank20_manifest(
        raw_output,
        maximum_reduction_prime=args.maximum_reduction_prime,
    )
    print(
        "FERMIGIER_RANK20_NEAR_MISS "
        f"rank_lower_bound={manifest['point_cloud']['selected_count']} "
        f"conductor={manifest['global_curve']['conductor']}"
    )
    if args.output is not None:
        write_json_exclusively(
            args.output, json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        )
        print(f"WROTE {args.output}")


if __name__ == "__main__":
    main()
