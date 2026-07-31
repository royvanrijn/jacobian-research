#!/usr/bin/env python3
"""Replay the exact support-at-most-five bidegree-(3,3) exclusion."""

from __future__ import annotations

import hashlib
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


def run(command: list[str]) -> None:
    subprocess.run(command, cwd=ROOT, check=True)


def load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    started = time.monotonic()
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

    three = load(THREE)
    four = load(FOUR)
    five = load(FIVE)
    if three["support_count"] != 560:
        raise AssertionError("unexpected three-support count")
    if not three["all_mixed_supports_excluded"]:
        raise AssertionError("a mixed three-support survived")
    if four["mixed_support_count"] != 1401:
        raise AssertionError("unexpected mixed four-support count")
    if not four["complete_mixed_enumeration"]:
        raise AssertionError("the mixed four-support enumeration is incomplete")
    if five["mixed_support_count"] != 3864:
        raise AssertionError("unexpected mixed five-support count")
    if not five["complete_mixed_enumeration"]:
        raise AssertionError("the mixed five-support enumeration is incomplete")
    if not all(
        record["status"] == "excluded_on_coefficient_torus"
        for record in four["records"]
    ):
        raise AssertionError("a mixed four-support was not excluded")
    if not all(
        record["status"] == "excluded_on_coefficient_torus"
        for record in five["records"]
    ):
        raise AssertionError("a mixed five-support was not excluded")

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
