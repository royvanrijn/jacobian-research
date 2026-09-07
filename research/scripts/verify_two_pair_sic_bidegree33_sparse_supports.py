#!/usr/bin/env python3
"""Replay the exact support-at-most-five bidegree-(3,3) exclusion."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
from itertools import combinations
import json
from pathlib import Path
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT_DIRECTORY = ROOT / "artifacts" / "generated-results"
THREE = (
    ARTIFACT_DIRECTORY
    / "two_pair_sic_bidegree33_sparse_three_support_screen.json"
)
FOUR = (
    ARTIFACT_DIRECTORY
    / "two_pair_sic_bidegree33_sparse_four_support_screen.json"
)
FIVE = (
    ARTIFACT_DIRECTORY
    / "two_pair_sic_bidegree33_sparse_five_support_screen.json"
)
OUTPUT = (
    ARTIFACT_DIRECTORY
    / "two_pair_sic_bidegree33_sparse_supports_certificate.json"
)
POSITIONS = tuple((row, column) for row in range(4) for column in range(4))


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--audit-existing-only",
        action="store_true",
        help=(
            "validate exact stored support keys and outcome counts without "
            "rerunning the three elimination screens"
        ),
    )
    return parser.parse_args()


def run(command: list[str]) -> None:
    subprocess.run(command, cwd=ROOT, check=True)


def load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def support_key(record: dict[str, object]) -> tuple[tuple[int, int], ...]:
    return tuple(
        tuple(map(int, position)) for position in record["support"]
    )


def exact_mixed_supports(
    support_size: int,
) -> tuple[tuple[tuple[int, int], ...], ...]:
    return tuple(
        support
        for support in combinations(POSITIONS, support_size)
        if min(row - column for row, column in support) < 0
        and max(row - column for row, column in support) > 0
    )


def validate_record_domain(
    payload: dict[str, object],
    expected: tuple[tuple[tuple[int, int], ...], ...],
    expected_counts: dict[str, int],
    label: str,
) -> None:
    records = payload["records"]
    if not isinstance(records, list) or len(records) != len(expected):
        raise AssertionError(f"{label} record count does not match its domain")
    actual = tuple(support_key(record) for record in records)
    if actual != expected:
        raise AssertionError(
            f"{label} records do not match the exact ordered support domain"
        )
    actual_counts = Counter(str(record["status"]) for record in records)
    if dict(actual_counts) != expected_counts or payload["counts"] != expected_counts:
        raise AssertionError(f"{label} stored outcome counts disagree with records")


def validate_artifacts() -> tuple[dict[str, object], ...]:
    three = load(THREE)
    four = load(FOUR)
    five = load(FIVE)
    expected_three = tuple(combinations(POSITIONS, 3))
    expected_four = exact_mixed_supports(4)
    expected_five = exact_mixed_supports(5)
    if three["support_count"] != len(expected_three):
        raise AssertionError("unexpected three-support count")
    if not three["all_mixed_supports_excluded"]:
        raise AssertionError("a mixed three-support survived")
    validate_record_domain(
        three,
        expected_three,
        {
            "covered_by_exact_diagonal_slice": 4,
            "excluded_on_coefficient_torus": 324,
            "reduces_to_at_most_two_diagonal_entries": 192,
            "strict_half_space_nullcone": 40,
        },
        "size-three census",
    )
    if (
        four["mixed_support_count"] != len(expected_four)
        or not four["complete_mixed_enumeration"]
    ):
        raise AssertionError("unexpected mixed four-support domain")
    validate_record_domain(
        four,
        expected_four,
        {"excluded_on_coefficient_torus": 1401},
        "size-four census",
    )
    if (
        five["mixed_support_count"] != len(expected_five)
        or not five["complete_mixed_enumeration"]
    ):
        raise AssertionError("unexpected mixed five-support domain")
    validate_record_domain(
        five,
        expected_five,
        {"excluded_on_coefficient_torus": 3864},
        "size-five census",
    )
    return three, four, five


def main() -> None:
    options = parse_arguments()
    started = time.monotonic()
    if not options.audit_existing_only:
        run(
            [
                sys.executable,
                "scripts/research_two_pair_sic_bidegree33_sparse_counterexample.py",
                "--through",
                "12",
                "--output",
                str(THREE),
            ]
        )
        run(
            [
                sys.executable,
                "scripts/research_two_pair_sic_bidegree33_sparse_four_counterexample.py",
                "--through",
                "12",
                "--output",
                str(FOUR),
            ]
        )
        run(
            [
                sys.executable,
                "scripts/research_two_pair_sic_bidegree33_sparse_five_counterexample.py",
                "--backend",
                "msolve",
                "--through",
                "12",
                "--timeout",
                "60",
                "--output",
                str(FIVE),
            ]
        )

    three, four, five = validate_artifacts()
    if options.audit_existing_only:
        print("PASS all 560 size-three supports occur exactly once")
        print("PASS all 1401 mixed size-four supports occur exactly once")
        print("PASS all 3864 mixed size-five supports occur exactly once")
        return

    payload = {
        "calculation": (
            "two_pair_sic_bidegree33_sparse_supports_certificate"
        ),
        "status": "proved",
        "scope": (
            "exact characteristic-zero exclusion through mu12 of every "
            "mixed-sign coefficient support of size at most five; "
            "nonmixed supports reduce to the exact diagonal theorem and "
            "a strict Hilbert--Mumford half-space"
        ),
        "orders": [1, 12],
        "support_counts": {
            "all_size_three": 560,
            "mixed_size_three": int(
                three["counts"]["excluded_on_coefficient_torus"]
            ),
            "mixed_size_four": 1401,
            "mixed_size_five": 3864,
        },
        "formula_checks": {
            "size_four": four["independent_formula_checks"],
            "size_five": five["independent_formula_check"],
        },
        "artifacts": {
            path.name: f"sha256:{digest(path)}"
            for path in (THREE, FOUR, FIVE)
        },
        "seconds": round(time.monotonic() - started, 6),
    }
    OUTPUT.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("PASS all 560 three-entry supports are classified")
    print("PASS all 1,401 mixed four-entry supports are excluded")
    print("PASS all 3,864 mixed five-entry supports are excluded")
    print("PASS every bidegree-(3,3) support of size at most five is safe")


if __name__ == "__main__":
    main()
