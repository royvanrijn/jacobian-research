#!/usr/bin/env python3
"""Dependency-free replay of the mixed BCH leading-coefficient recurrence."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


BERNOULLI_EVEN = {
    2: Fraction(1, 6),
    4: Fraction(-1, 30),
    6: Fraction(1, 42),
    8: Fraction(-1, 30),
    10: Fraction(5, 66),
    12: Fraction(-691, 2730),
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "certificate",
        nargs="?",
        type=Path,
        default=Path("artifacts/generated-results/lr_mixed_bch_classes.json"),
    )
    arguments = parser.parse_args()
    data = json.loads(arguments.certificate.read_text(encoding="utf-8"))
    rows = data["computed_regression"]

    previous = None
    for row in rows:
        k = row["k"]
        leading = Fraction(row["residue_leading_coefficient"])
        assert leading != 0
        assert row["residue_degree_u"] == 4 * k + 11
        if previous is None:
            assert leading == Fraction(14438891520, 2401)
        else:
            assert leading == -73440 * (k + 2) * (2 * k + 5) * previous

        bernoulli = BERNOULLI_EVEN[2 * k]
        scalar = (
            -bernoulli
            * math.comb(2 * k, k)
            / math.factorial(2 * k)
        )
        assert scalar == Fraction(row["balanced_bch_scalar"])
        assert scalar * leading == Fraction(
            row["balanced_bch_leading_coefficient"]
        )
        previous = leading

    assert data["commutation"]["D_B_D_C_bracket_zero"] is True
    assert "s^k*t^k" in data["descent_audit"]["lower_target_amplitudes"]
    print("PASS: exact mixed-word leading recurrence")
    print("PASS: exact Bernoulli/binomial BCH scalars")
    print("PASS: every recorded balanced BCH normal coefficient is nonzero")
    print("PASS: certificate records the lower-jet universality failure")


if __name__ == "__main__":
    main()
