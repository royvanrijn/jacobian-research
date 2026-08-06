#!/usr/bin/env python3
"""Verify all dense off-diagonal two-row rank-two quartic charts.

There are ten row pairs and six orbits under simultaneous reversal.  On
each pair, retain all eight off-diagonal positions in the two rows.  The
factor gauge is fixed by U=(e_r,e_s), and the minor in columns r,s is
nonzero on the dense coefficient torus, so every point has exact rank
two.  This checker proves over QQ that moments through mu_8 generate the
unit ideal on every orbit representative.

Coordinate boundaries and non-coordinate U charts are outside the scope.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
import json
from pathlib import Path

from explore_two_pair_sic_bidegree44_two_row_off_diagonal import (
    ROW_PAIR_REPRESENTATIVES,
    classify,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree44_two_row_off_diagonal.json"
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    if arguments.workers < 1 or arguments.timeout < 1:
        raise ValueError("--workers and --timeout must be positive")

    covered_pairs: set[tuple[int, int]] = set()
    for rows in ROW_PAIR_REPRESENTATIVES:
        covered_pairs.add(rows)
        covered_pairs.add((4 - rows[1], 4 - rows[0]))
    if covered_pairs != set(combinations(range(5), 2)):
        raise AssertionError("row-pair reversal orbits do not cover all ten pairs")

    records: list[dict[str, object]] = []
    with ProcessPoolExecutor(max_workers=arguments.workers) as executor:
        futures = [
            executor.submit(classify, rows, 8, arguments.timeout)
            for rows in ROW_PAIR_REPRESENTATIVES
        ]
        for future in as_completed(futures):
            record = future.result()
            result = record["classification_through_requested_order"]
            if (
                result["returncode"] != 0
                or result["status"] != "unit_ideal"
                or record["component_rational_parametrization"] is not None
            ):
                raise AssertionError(f"nonunit off-diagonal row chart: {record}")
            if record["residual_coordinate_count"] != 6:
                raise AssertionError("unexpected quotient coordinate count")
            records.append(record)
    records.sort(key=lambda record: record["row_pair"])
    if len(records) != 6:
        raise AssertionError("not every row-pair orbit was classified")

    artifact = {
        "format": "two-pair-sic-bidegree44-two-row-off-diagonal-v1",
        "field": "characteristic zero",
        "factor_chart": {
            "normal_form": "U=(e_r,e_s), C=U*B",
            "internal_GL2_gauge_removed": True,
            "dense_support": (
                "all positions in rows r,s except the two diagonal positions"
            ),
            "orbit_normalization": (
                "overall scaling and diagonal torus set two distinct-weight "
                "coefficients to one"
            ),
            "residual_coordinates": 6,
            "exact_rank_two_minor": "det C[rows r,s; columns r,s]=-c_rs*c_sr",
        },
        "row_pair_count": 10,
        "reversal_orbit_count": 6,
        "moment_orders": [1, 8],
        "exact_unit_system_count": 6,
        "component_conclusion": (
            "the dense coefficient-torus moment scheme through mu_8 is "
            "empty over QQ on every row-pair orbit"
        ),
        "recurrence_gate": (
            "no exact component survives through mu_8, so no recurrence or "
            "mixed-tail search is required on these dense charts"
        ),
        "global_conclusion": (
            "all ten dense off-diagonal two-row exact-rank-two charts contain "
            "no all-order pure-moment point"
        ),
        "scope": (
            "complete dense-torus classification of this ten-chart coordinate "
            "family; coordinate boundaries and generic factor charts remain open"
        ),
        "orbits": records,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")

    print("PASS six reversal orbits cover all ten row pairs")
    print("PASS every dense chart has exact coefficient rank two")
    print("PASS all six QQ coefficient-torus ideals are units through mu_8")


if __name__ == "__main__":
    main()
