#!/usr/bin/env python3
"""Exact finite controls for the unbounded reciprocal-stability obstruction.

The written note proves the all-parameter theorem. This script derives the
recorded rational matrices directly from factorial exponent columns and
retains exact inertia certificates. No numerical eigenvalue solver is used.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, factorial, prod
from pathlib import Path
import sys

if not __debug__:
    raise SystemExit("REJECT: optimized execution is not an approved replay")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
NOTE = HERE / "SFC_DEAMPLIFICATION_BRIDGES.md"
EXPECTED_NOTE = "df9f19e4418c2abdf47c0f2696a60637dcfb761f93a9b91457aad4b5c671e6b0"
IMPORTED = ROOT / "research/factorial/image-rigidity/SFC_BALANCED_PRIVATE_AMPLIFICATION.md"
EXPECTED_IMPORTED = "f4ebaa9028cb04258cf65be178685650d53e2852644ed0678f1a655d134c065e"


def require(condition: bool, description: str) -> None:
    if not condition:
        raise ValueError(description)


def pin(path: Path, expected: str) -> str:
    actual = sha256(path.read_bytes()).hexdigest()
    require(actual == expected, f"source pin mismatch: {path.name}")
    return actual


def rational(value: Fraction | int) -> dict[str, int]:
    value = Fraction(value)
    return {"numerator": value.numerator, "denominator": value.denominator}


def serial_matrix(matrix: list[list[Fraction]]) -> list[list[dict[str, int]]]:
    return [[rational(value) for value in row] for row in matrix]


def moment(columns: list[tuple[int, ...]], powers: tuple[int, ...]) -> int:
    require(len(columns) == len(powers), "column/power shape")
    require(len({len(col) for col in columns}) == 1, "ragged columns")
    return prod(
        factorial(sum(powers[j] * columns[j][i] for j in range(len(columns))))
        for i in range(len(columns[0]))
    )


def reciprocal_matrix(columns: list[tuple[int, ...]]) -> list[list[Fraction]]:
    """Matrix of Q_2, including both multinomial and symmetric factors."""
    count = len(columns)
    matrix = [[Fraction(0) for _ in columns] for _ in columns]
    for i in range(count):
        powers = tuple(2 if j == i else 0 for j in range(count))
        matrix[i][i] = Fraction(1, moment(columns, powers))
        for j in range(i):
            powers = tuple(int(k == i) + int(k == j) for k in range(count))
            matrix[i][j] = matrix[j][i] = Fraction(1, 4 * moment(columns, powers))
    return matrix


def scale(matrix: list[list[Fraction]], factor: int) -> list[list[Fraction]]:
    return [[entry * factor for entry in row] for row in matrix]


def det3(matrix: list[list[Fraction]]) -> Fraction:
    return sum(
        matrix[0][j]
        * (
            matrix[1][(j + 1) % 3] * matrix[2][(j + 2) % 3]
            - matrix[1][(j + 2) % 3] * matrix[2][(j + 1) % 3]
        )
        for j in range(3)
    )


def pivots(matrix: list[list[Fraction]]) -> list[Fraction]:
    d1 = matrix[0][0]
    d2 = matrix[0][0] * matrix[1][1] - matrix[0][1] ** 2
    d3 = det3(matrix)
    require(d1 != 0 and d2 != 0 and d3 != 0, "singular leading minor")
    return [d1, d2 / d1, d3 / d2]


def bilinear(left: list[int], right: list[int], matrix: list[list[Fraction]]) -> Fraction:
    return sum(left[i] * matrix[i][j] * right[j] for i in range(3) for j in range(3))


def source_reciprocal_controls() -> int:
    columns = [
        (1, 1, 0, 0, 0, 0),
        (0, 0, 1, 1, 0, 0),
        (0, 0, 0, 0, 1, 1),
    ]
    count = 0
    for m in range(1, 6):
        for k1 in range(m + 1):
            for k2 in range(m - k1 + 1):
                powers = (k1, k2, m - k1 - k2)
                multinomial = factorial(m) // prod(factorial(k) for k in powers)
                actual_coefficient = multinomial * moment(columns, powers)
                require(actual_coefficient == factorial(m) * prod(factorial(k) for k in powers), "source moment identity")
                require(Fraction(1, actual_coefficient) == Fraction(multinomial, factorial(m) ** 2), "source reciprocal identity")
                count += 1
    return count


def small_row_control() -> dict:
    columns = [
        (1, 1, 0, 0, 0, 0, 0),
        (0, 0, 1, 1, 0, 0, 1),
        (0, 0, 0, 0, 1, 1, 2),
    ]
    matrix = scale(reciprocal_matrix(columns), 4)
    expected = [
        [Fraction(1), Fraction(1), Fraction(1, 2)],
        [Fraction(1), Fraction(1, 2), Fraction(1, 6)],
        [Fraction(1, 2), Fraction(1, 6), Fraction(1, 24)],
    ]
    require(matrix == expected, "weighted-row matrix mismatch")
    require(det3(matrix) == Fraction(-1, 144), "weighted-row determinant")
    require(pivots(matrix) == [Fraction(1), Fraction(-1, 2), Fraction(1, 72)], "weighted-row inertia pivots")
    x, y = [23, -94, 140], [1, 1, 1]
    xy, xx, yy = bilinear(x, y, matrix), bilinear(x, x, matrix), bilinear(y, y, matrix)
    require((xy, xx, yy) == (0, 273, Fraction(39, 8)), "explicit half-plane zero bilinear identities")
    require(xx / yy == 56, "explicit half-plane zero scale")
    return {
        "columns": columns,
        "scaled_reciprocal_matrix": serial_matrix(matrix),
        "determinant": rational(det3(matrix)),
        "ldl_pivots": [rational(value) for value in pivots(matrix)],
        "zero_real_numerator": x,
        "zero_real_denominator_squared": 56,
        "zero_imaginary_part": y,
        "bilinear_identities": {"xMy": rational(xy), "xMx": rational(xx), "yMy": rational(yy)},
    }


def amplified_control(D: int) -> dict:
    a = 4 * D
    columns = [
        (D, 0, 0, 2 * a, a, 0),
        (0, D, 0, a, 2 * a, 0),
        (0, 0, D, 0, 2 * a, a),
    ]
    profiles = [sorted(value for value in col if value) for col in columns]
    require(profiles[0] == profiles[1] == profiles[2] == [D, a, 2 * a], "balanced profile")
    H = factorial(4 * a) * factorial(2 * a)
    Hij = {
        "12": factorial(3 * a) ** 2,
        "13": factorial(2 * a) * factorial(3 * a) * factorial(a),
        "23": factorial(a) ** 2 * factorial(4 * a),
    }
    R = {ij: Fraction(H, value) for ij, value in Hij.items()}
    require(R["13"] == comb(4 * a, a) and R["23"] == comb(2 * a, a), "ratio cancellations")
    b = Fraction(comb(2 * D, D), 4)
    A, B, C = [b * R[ij] for ij in ("12", "13", "23")]
    matrix = scale(reciprocal_matrix(columns), factorial(2 * D) * H)
    require(matrix == [[1, A, B], [A, 1, C], [B, C, 1]], "direct factorial reciprocal matrix")
    L = Fraction(factorial(a) * factorial(3 * a), factorial(2 * a) ** 2)
    product_L = prod(Fraction(2 * a + i, a + i) for i in range(1, a + 1))
    require(L == product_L == R["13"] / (R["12"] * R["23"]), "L identity")
    require(L >= Fraction(3, 2) ** a > 4 ** D > 2 * b, "uniform L lower bound control")
    require(B > 2 * A * C and A > 1 and C > 1, "determinant comparison hypotheses")
    determinant = det3(matrix)
    require(determinant == 1 + 2 * A * B * C - A * A - B * B - C * C, "determinant expansion")
    require(determinant < 1 - A * A - C * C < 0, "negative determinant")
    diagonal = pivots(matrix)
    require([value > 0 for value in diagonal] == [True, False, True], "inertia must be (2,1)")
    return {
        "D": D,
        "a": a,
        "columns": columns,
        "SFC_theorem_threshold_met": D >= 18,
        "b_D": rational(b),
        "R": {ij: rational(value) for ij, value in R.items()},
        "L_a": rational(L),
        "scaled_reciprocal_matrix": serial_matrix(matrix),
        "determinant": rational(determinant),
        "ldl_pivots": [rational(value) for value in diagonal],
        "inertia": {"positive": 2, "negative": 1, "zero": 0},
    }


def private_only_control(D: int) -> dict:
    columns = [tuple(D if i == j else 0 for i in range(3)) for j in range(3)]
    matrix = scale(reciprocal_matrix(columns), factorial(2 * D))
    b = Fraction(comb(2 * D, D), 4)
    require(matrix == [[1, b, b], [b, 1, b], [b, b, 1]], "private-only matrix")
    diagonal = pivots(matrix)
    require([value > 0 for value in diagonal] == [True, False, False], "private-only quadratic has (1,2) inertia")
    return {"D": D, "ldl_pivots": [rational(value) for value in diagonal], "inertia": [1, 2, 0]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note", type=Path, default=NOTE)
    parser.add_argument("--receipt", type=Path, default=HERE / "sfc_deamplification_reciprocal_verification.json")
    args = parser.parse_args()
    sources = {"note": pin(args.note, EXPECTED_NOTE), "balanced_amplification": pin(IMPORTED, EXPECTED_IMPORTED)}
    receipt = {
        "status": "PASS",
        "scope": "Exact finite algebraic and inertia controls; the written proof, not this finite replay, proves the unbounded family.",
        "source_sha256": sources,
        "source_reciprocal_coefficient_controls": source_reciprocal_controls(),
        "small_weighted_row": small_row_control(),
        "amplified_family": [amplified_control(D) for D in (2, 3, 18, 19, 24)],
        "private_only_controls": [private_only_control(D) for D in (2, 18)],
        "no_SFC_counterexample_claim": True,
    }
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print("PASS: weighted integer row, explicit upper-half-plane zero, five amplified inertia certificates, two private-only controls")
    print(f"Receipt: {args.receipt.relative_to(ROOT) if args.receipt.is_relative_to(ROOT) else args.receipt}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError) as error:
        raise SystemExit(f"REJECT: {error}") from error
