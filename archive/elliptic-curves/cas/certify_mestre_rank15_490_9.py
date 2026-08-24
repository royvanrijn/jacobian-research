#!/usr/bin/env python3
"""Freeze the exact rank-15 certificate found at T=490/9.

The input is the pinned, leakage-controlled two-family frontier artifact.  This
script performs no point search: it replays the existing H=10^6 numerical
subset, verifies every point on the exact short Weierstrass model, and proves
their independence with the finite-reduction mod-3 certificate machinery.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import json
import os
from pathlib import Path
import sys
from typing import Any

from search_mestre_rank14_pair_rational_frontier import family_coefficients
from search_mestre_root_tuple_scale import (
    point_digest,
    point_on_short_curve,
    sha256_file,
)
from search_mestre_root_tuple_scale_max100 import stable_json_digest
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

FRONTIER_FILENAME = "elliptic_mestre_rank14_pair_rational_frontier.json"
EXPECTED_FRONTIER_SHA256 = (
    "87e2d278cc1ee0653d1a4f871c1e34ed3d03babe1c1cd2ffe6712b7608efaee7"
)
EXPECTED_FRONTIER_RESULT_SHA256 = (
    "d33c0cf0a5e2364bd49e18363e0a1f3ca51512fdf60be7052dd05f1cbfa9d610"
)
EXPECTED_FRONTIER_SCRIPT_SHA256 = (
    "2f6251c67e2eb3cee2eca37d7e866913e9d5de73d30e3bfcb253641454d40d5f"
)
FAMILY_INDEX = 1
ROOTS = (0, 7, 121, 128, 183, 194)
PARAMETER = Q(490, 9)
STAGE = "H1000000"
CERTIFICATE_PRIME_BOUND = 499


def rational_string(value: Fraction) -> str:
    value = Q(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def exclusive_write(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[2]
    parser.add_argument(
        "--frontier",
        type=Path,
        default=root / "artifacts/generated-results" / FRONTIER_FILENAME,
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts/generated-results"
            / "elliptic_mestre_rank15_490_9.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.output.exists():
        raise SystemExit("refusing to overwrite the rank-15 certificate artifact")
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    frontier_script = script_path.with_name(
        "search_mestre_rank14_pair_rational_frontier.py"
    )
    if sha256_file(args.frontier) != EXPECTED_FRONTIER_SHA256:
        raise SystemExit("the pinned rational-frontier artifact changed")
    if sha256_file(frontier_script) != EXPECTED_FRONTIER_SCRIPT_SHA256:
        raise SystemExit("the pinned rational-frontier script changed")
    frontier = json.loads(args.frontier.read_text())
    if frontier["result_sha256"] != EXPECTED_FRONTIER_RESULT_SHA256:
        raise AssertionError("the rational-frontier result digest changed")
    if tuple(frontier["scope"]["included_families"][FAMILY_INDEX]) != ROOTS:
        raise AssertionError("the certified family changed")

    record = next(
        row
        for row in frontier["selected_records"]
        if row["family_index"] == FAMILY_INDEX
        and row["numerator"] == PARAMETER.numerator
        and row["denominator"] == PARAMETER.denominator
    )
    stage = record["point_stages"][STAGE]
    if stage["status"] != "completed" or stage["stable_numerical_rank"] != 15:
        raise AssertionError("the frozen H=10^6 leader changed")
    points = tuple(
        (Q(row["x"]), Q(row["y"])) for row in stage["numerical_subset"]
    )
    coefficients = family_coefficients(FAMILY_INDEX, PARAMETER)
    if len(points) != 15 or not all(
        point_on_short_curve(coefficients, point) for point in points
    ):
        raise AssertionError("the proposed exact rank-15 point set changed")
    if point_digest(points) != (
        "fcbc9e86472490b3e5db4980d8d579c42657c41e3685be61953e1b95f5f9ed8e"
    ):
        raise AssertionError("the proposed exact rank-15 point digest changed")

    certificate = mod3_independence_certificate(
        coefficients, points, prime_bound=CERTIFICATE_PRIME_BOUND
    )
    if certificate["certified_algebraic_rank_lower_bound"] != 15:
        raise AssertionError("finite reduction did not certify all 15 points")
    conductor = record["conductor_phase"]
    if not conductor["status"].startswith("completed exact"):
        raise AssertionError("the exact conductor record changed")

    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "certified exact algebraic rank lower bound 15",
        "curve": {
            "roots": list(ROOTS),
            "parameter": rational_string(PARAMETER),
            "sign_equivalent_parameter": rational_string(-PARAMETER),
            "weierstrass_coefficients": [
                rational_string(value) for value in coefficients
            ],
            "conductor": conductor["conductor"],
            "log_conductor": conductor["log_conductor"],
            "minimal_model": conductor["minimal_model"],
            "minimal_discriminant": conductor["minimal_discriminant"],
            "root_number": conductor["root_number"],
            "strict_log_conductor_target": "182.72",
            "below_strict_log_conductor_target_numerically": conductor[
                "below_strict_log_conductor_target_numerically"
            ],
        },
        "point_source": {
            "frontier_stage": STAGE,
            "height_bound": stage["height_bound"],
            "mapping_truncated": stage["mapping_truncated"],
            "point_count": len(points),
            "point_sha256": point_digest(points),
            "exact_curve_membership_replayed": True,
            "selection_statement": (
                "the finite certificate replays the fixed H=10^6 subset and "
                "performs no new parameter or point search"
            ),
        },
        "points": [
            {"x": rational_string(x_value), "y": rational_string(y_value)}
            for x_value, y_value in points
        ],
        "finite_reduction_certificate": certificate,
        "claim": {
            "certified_algebraic_rank_lower_bound": 15,
            "does_not_claim_exact_mordell_weil_rank": True,
            "does_not_hit_rank21_or_rank30_target": True,
        },
        "provenance": {
            "script_path": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "frontier_path": str(args.frontier.relative_to(root)),
            "frontier_sha256": EXPECTED_FRONTIER_SHA256,
            "frontier_result_sha256": EXPECTED_FRONTIER_RESULT_SHA256,
            "frontier_script_sha256": EXPECTED_FRONTIER_SCRIPT_SHA256,
            "certificate_prime_bound": CERTIFICATE_PRIME_BOUND,
            "external_process_calls": 0,
        },
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    artifact["result_sha256"] = stable_json_digest(
        {
            "curve": artifact["curve"],
            "point_source": artifact["point_source"],
            "points": artifact["points"],
            "certificate": certificate,
            "claim": artifact["claim"],
        }
    )
    exclusive_write(args.output, artifact)
    print(
        f"certified rank>={certificate['certified_algebraic_rank_lower_bound']} "
        f"at T={rational_string(PARAMETER)} output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
