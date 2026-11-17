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
    support_for,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree44_two_row_off_diagonal.json"
)


def row_pair_orbits() -> tuple[tuple[int, int], ...]:
    covered_pairs: set[tuple[int, int]] = set()
    for rows in ROW_PAIR_REPRESENTATIVES:
        covered_pairs.add(rows)
        covered_pairs.add((4 - rows[1], 4 - rows[0]))
    if covered_pairs != set(combinations(range(5), 2)):
        raise AssertionError("row-pair reversal orbits do not cover all ten pairs")
    return ROW_PAIR_REPRESENTATIVES


def validate_existing_artifact(path: Path) -> None:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    if artifact.get("format") != "two-pair-sic-bidegree44-two-row-off-diagonal-v1":
        raise AssertionError("unexpected dense off-diagonal artifact format")

    representatives = row_pair_orbits()
    records = artifact.get("orbits")
    if not isinstance(records, list):
        raise AssertionError("dense off-diagonal artifact has no orbits list")
    stored_pairs = [tuple(record.get("row_pair", ())) for record in records]
    if stored_pairs != list(representatives):
        raise AssertionError(
            "stored row pairs are not the exact ordered reversal representatives"
        )
    if len(stored_pairs) != len(set(stored_pairs)):
        raise AssertionError("stored dense chart census contains duplicate row pairs")
    if artifact.get("row_pair_count") != 10:
        raise AssertionError("unexpected declared row-pair count")
    if artifact.get("reversal_orbit_count") != len(representatives) or len(
        representatives
    ) != 6:
        raise AssertionError("unexpected declared reversal-orbit count")
    if artifact.get("exact_unit_system_count") != len(representatives):
        raise AssertionError("stored unit-system total is not the exact census")
    if artifact.get("moment_orders") != [1, 8]:
        raise AssertionError("dense theorem artifact must record moments 1 through 8")

    for record, rows in zip(records, representatives, strict=True):
        reversal_partner = [4 - rows[1], 4 - rows[0]]
        if record.get("reversal_partner") != reversal_partner:
            raise AssertionError(f"wrong reversal partner for row pair {rows}")
        if record.get("support") != [list(position) for position in support_for(rows)]:
            raise AssertionError(f"wrong dense support for row pair {rows}")
        if record.get("residual_coordinate_count") != 6:
            raise AssertionError(f"wrong quotient dimension for row pair {rows}")
        profiles = record.get("moment_profiles")
        if not isinstance(profiles, list) or [
            profile.get("order") for profile in profiles
        ] != list(range(1, 9)):
            raise AssertionError(f"incomplete moment profile for row pair {rows}")
        result = record.get("classification_through_requested_order")
        if (
            not isinstance(result, dict)
            or result.get("returncode") != 0
            or result.get("status") != "unit_ideal"
        ):
            raise AssertionError(f"stored dense chart is not a unit for {rows}")
        if record.get("component_rational_parametrization") is not None:
            raise AssertionError(f"spurious surviving component for row pair {rows}")

    print("PASS six exact row-pair keys cover all ten reversal-related charts")
    print("PASS every stored dense system is a QQ unit through mu_8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--audit-existing-only",
        action="store_true",
        help=(
            "validate exact stored row-pair coverage and unit outcomes without "
            "rerunning msolve"
        ),
    )
    arguments = parser.parse_args()
    if arguments.workers < 1 or arguments.timeout < 1:
        raise ValueError("--workers and --timeout must be positive")
    if arguments.audit_existing_only:
        validate_existing_artifact(arguments.output)
        return

    representatives = row_pair_orbits()

    records: list[dict[str, object]] = []
    with ProcessPoolExecutor(max_workers=arguments.workers) as executor:
        futures = [
            executor.submit(classify, rows, 8, arguments.timeout)
            for rows in representatives
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
