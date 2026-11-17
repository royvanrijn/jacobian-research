#!/usr/bin/env python3
"""Certify all complex points on the two six-support survivor systems."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import combinations
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "artifacts" / "generated-results"
SIX_CENSUS = (
    GENERATED
    / "two_pair_sic_bidegree33_sparse_six_support_screen.json"
)
SEVEN_CENSUS = (
    GENERATED
    / "two_pair_sic_bidegree33_sparse_support7_screen.json"
)
OUTPUT = (
    GENERATED
    / "two_pair_sic_bidegree33_sparse_survivor_rur.json"
)
TARGETS = (
    (
        3568,
        ((0, 1), (0, 3), (1, 0), (2, 1), (2, 3), (3, 2)),
        (-1, -2, -1, -1),
    ),
    (
        4106,
        ((0, 1), (1, 0), (1, 2), (2, 3), (3, 0), (3, 2)),
        (2, 1, -1, -1),
    ),
)
POSITIONS = tuple((row, column) for row in range(4) for column in range(4))


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--audit-census-only",
        action="store_true",
        help=(
            "validate exact size-six/seven support coverage and stored "
            "outcome counts without rerunning the two RURs"
        ),
    )
    return parser.parse_args()


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalized_result(text: str) -> str:
    return " ".join(text.split())


def expected_mixed_supports(
    support_size: int,
) -> tuple[tuple[tuple[int, int], ...], ...]:
    return tuple(
        support
        for support in combinations(POSITIONS, support_size)
        if min(row - column for row, column in support) < 0
        and max(row - column for row, column in support) > 0
    )


def validate_complete_census(
    payload: dict[str, object],
    support_size: int,
    expected_counts: dict[str, int],
) -> None:
    expected = expected_mixed_supports(support_size)
    recorded_support_size = payload.get("support_size")
    if recorded_support_size is None:
        # The first complete size-six artifact predates the explicit field.
        # Accept that immutable schema only for its exact calculation name;
        # current and future merger outputs always carry ``support_size``.
        if not (
            support_size == 6
            and payload.get("calculation")
            == "two_pair_sic_bidegree33_sparse_six_support_screen"
        ):
            raise AssertionError("missing support_size outside legacy size six")
        recorded_support_size = 6
    if (
        recorded_support_size != support_size
        or not payload["complete_mixed_enumeration"]
        or payload["start"] != 0
        or payload["stop"] != len(expected)
        or payload["mixed_support_count"] != len(expected)
    ):
        raise AssertionError(f"unexpected size-{support_size} census scope")
    records = payload["records"]
    if not isinstance(records, list) or len(records) != len(expected):
        raise AssertionError(f"incomplete size-{support_size} record list")
    actual = tuple(
        tuple(tuple(map(int, position)) for position in record["support"])
        for record in records
    )
    if actual != expected:
        raise AssertionError(
            f"size-{support_size} records do not match the exact ordered domain"
        )
    actual_counts = Counter(str(record["status"]) for record in records)
    if dict(actual_counts) != expected_counts or payload["counts"] != expected_counts:
        raise AssertionError(
            f"unexpected size-{support_size} census counts: {actual_counts}"
        )


def run_target(
    start: int,
    support: tuple[tuple[int, int], ...],
    residual: tuple[int, ...],
    directory: Path,
) -> dict[str, object]:
    path = directory / f"support-{start}.json"
    subprocess.run(
        [
            sys.executable,
            "scripts/research_two_pair_sic_bidegree33_sparse_six_counterexample.py",
            "--support-size",
            "6",
            "--through",
            "12",
            "--start",
            str(start),
            "--limit",
            "1",
            "--timeout",
            "60",
            "--rational-parametrization",
            "--output",
            str(path),
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(path.read_text(encoding="utf-8"))
    record = payload["records"][0]
    if tuple(map(tuple, record["support"])) != support:
        raise AssertionError(f"unexpected support at mixed index {start}")
    if record["status"] != "msolve_nonempty":
        raise AssertionError(f"unexpected msolve status at {start}")
    result = normalized_result(record["msolve"]["result_head"])
    if "5, 1, ['h', 'z0', 'z1', 'z2', 'z3']" not in result:
        raise AssertionError(
            f"the RUR at {start} does not have degree one in five variables"
        )
    exact_box = "[[1 / 2, 1 / 2], " + ", ".join(
        f"[{value}, {value}]" for value in residual
    ) + "]"
    if exact_box not in result:
        raise AssertionError(f"unexpected unique complex point at {start}")

    matrix = sp.zeros(4)
    for coefficient, (i, j) in zip(
        (1, 1, *residual),
        support,
        strict=True,
    ):
        matrix[i, j] = coefficient
    if matrix.det() != 1 or matrix.rank() != 4:
        raise AssertionError("a six-support survivor is not full rank")
    return {
        "mixed_support_index": start,
        "support": [list(position) for position in support],
        "normalized_residual_coordinates": list(residual),
        "coefficient_matrix": [
            [int(matrix[i, j]) for j in range(4)]
            for i in range(4)
        ],
        "coefficient_determinant": 1,
        "coefficient_rank": 4,
        "complex_solution_count": 1,
        "rur_degree": 1,
        "rur": record["msolve"]["result_head"],
    }


def main() -> None:
    options = parse_arguments()
    six = json.loads(SIX_CENSUS.read_text(encoding="utf-8"))
    seven = json.loads(SEVEN_CENSUS.read_text(encoding="utf-8"))
    six_counts = {
        "excluded_on_coefficient_torus": 7586,
        "msolve_nonempty": 2,
    }
    seven_counts = {
        "excluded_on_coefficient_torus": 11200,
    }
    validate_complete_census(six, 6, six_counts)
    validate_complete_census(seven, 7, seven_counts)
    survivor_supports = {
        tuple(tuple(map(int, position)) for position in record["support"])
        for record in six["records"]
        if record["status"] == "msolve_nonempty"
    }
    expected_survivors = {support for _, support, _ in TARGETS}
    if survivor_supports != expected_survivors:
        raise AssertionError("unexpected size-six survivor supports")
    mixed_six = expected_mixed_supports(6)
    for index, support, _ in TARGETS:
        if mixed_six[index] != support:
            raise AssertionError(
                f"stored size-six target index {index} no longer matches its support"
            )
    if options.audit_census_only:
        print("PASS all 7588 mixed size-six supports occur exactly once")
        print("PASS all 11200 mixed size-seven supports occur exactly once")
        print("PASS the two size-six survivor indices and supports agree")
        return

    with tempfile.TemporaryDirectory(
        prefix="sic33-sparse-six-rur-"
    ) as directory:
        points = [
            run_target(start, support, residual, Path(directory))
            for start, support, residual in TARGETS
        ]

    payload = {
        "calculation": "two_pair_sic_bidegree33_sparse_survivor_rur",
        "status": "proved",
        "field": "characteristic zero",
        "scope": (
            "full complex rational-univariate certification of the two "
            "nonunit coefficient-torus systems in the complete size-six "
            "census; each RUR has degree one"
        ),
        "size_six_census": {
            "mixed_supports": 7588,
            "excluded_supports": 7586,
            "survivor_systems": 2,
            "sha256": digest(SIX_CENSUS),
        },
        "size_seven_census": {
            "mixed_supports": 11200,
            "excluded_supports": 11200,
            "sha256": digest(SEVEN_CENSUS),
        },
        "points": points,
        "conclusion": (
            "the only algebraic coefficient-torus points left by the "
            "size-six census are the two displayed full-rank "
            "Weyl/torus copies of the Rodrigues orbit"
        ),
    }
    OUTPUT.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("PASS both exceptional size-six systems have degree-one RURs")
    print("PASS both unique complex points have determinant one")
    print("PASS all size-seven coefficient tori are excluded")


if __name__ == "__main__":
    main()
