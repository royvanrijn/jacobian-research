#!/usr/bin/env python3
"""Exact replay of the bidegree-(4,4) SIC(2) counterexample.

For the contraction pairs (xi1,z1), (xi2,z2), set

    R = xi1*z1 + xi2*z2
    Z = xi1*z2
    W = 2*xi2*z1
    T = xi1*z1 - xi2*z2
    F = (R+Z) * (R^2*W - (2*R+Z)*T^2/2)
    Q = Z.

The written proof gives, for every m >= 1,

    E_2(F^m) = 0,
    E_2(Q*F^m) = (4*m+2)!*m!/(2*m+1)!!.

This dependency-free checker constructs the formula from R,Z,W,T, checks
the rank-one-cone relation T^2=R^2-2*Z*W, performs exact sparse
contractions through the declared cutoff, and replays the two finite-sum
identities underlying the all-order beta-integral proof.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import comb, factorial
from pathlib import Path


Exponent = tuple[int, int, int, int]  # xi1, xi2, z1, z2
Polynomial = dict[Exponent, Fraction]
CUTOFF = 8
IDENTITY_CUTOFF = 99
SCRIPT_PATH = Path(__file__).resolve()
ROOT = SCRIPT_PATH.parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_image_mathieu_counterexample.json"
)
ZERO = (0, 0, 0, 0)


def monomial(exponent: Exponent, coefficient: int = 1) -> Polynomial:
    return {exponent: Fraction(coefficient)}


def add(*polynomials: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            result[exponent] = result.get(exponent, Fraction(0)) + coefficient
    return {
        exponent: coefficient
        for exponent, coefficient in result.items()
        if coefficient
    }


def scale(coefficient: Fraction | int, polynomial: Polynomial) -> Polynomial:
    scalar = Fraction(coefficient)
    return {
        exponent: scalar * value
        for exponent, value in polynomial.items()
        if scalar * value
    }


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(
                left_exponent[index] + right_exponent[index]
                for index in range(4)
            )
            result[exponent] = (
                result.get(exponent, Fraction(0))
                + left_coefficient * right_coefficient
            )
    return {
        exponent: coefficient
        for exponent, coefficient in result.items()
        if coefficient
    }


def power(polynomial: Polynomial, exponent: int) -> Polynomial:
    result = monomial(ZERO)
    base = polynomial
    value = exponent
    while value:
        if value % 2:
            result = multiply(result, base)
        base = multiply(base, base)
        value //= 2
    return result


def contraction(polynomial: Polynomial) -> Polynomial:
    """Apply E_2 and retain exponents only in z1,z2 slots."""
    result: Polynomial = {}
    for (xi1, xi2, z1, z2), coefficient in polynomial.items():
        if xi1 > z1 or xi2 > z2:
            continue
        residual = (0, 0, z1 - xi1, z2 - xi2)
        value = (
            coefficient
            * Fraction(factorial(z1), factorial(z1 - xi1))
            * Fraction(factorial(z2), factorial(z2 - xi2))
        )
        result[residual] = result.get(residual, Fraction(0)) + value
    return {
        exponent: coefficient
        for exponent, coefficient in result.items()
        if coefficient
    }


def double_factorial_odd(order: int) -> int:
    result = 1
    for value in range(1, order + 1, 2):
        result *= value
    return result


def witness() -> tuple[Polynomial, Polynomial, dict[str, Polynomial]]:
    xi1z1 = monomial((1, 0, 1, 0))
    xi2z2 = monomial((0, 1, 0, 1))
    r = add(xi1z1, xi2z2)
    z = monomial((1, 0, 0, 1))
    w = scale(2, monomial((0, 1, 1, 0)))
    t = add(xi1z1, scale(-1, xi2z2))
    f = multiply(
        add(r, z),
        add(
            multiply(power(r, 2), w),
            scale(
                Fraction(-1, 2),
                multiply(add(scale(2, r), z), power(t, 2)),
            ),
        ),
    )
    return f, z, {"R": r, "Z": z, "W": w, "T": t}


def coefficient_matrix(f: Polynomial) -> list[list[Fraction]]:
    """Coefficient matrix in xi1^i xi2^(4-i), z1^j z2^(4-j)."""
    matrix = [[Fraction(0) for _ in range(5)] for _ in range(5)]
    for (xi1, xi2, z1, z2), coefficient in f.items():
        assert xi1 + xi2 == z1 + z2 == 4
        matrix[xi1][z1] = coefficient
    return matrix


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    work = [row[:] for row in matrix]
    result = Fraction(1)
    for column in range(len(work)):
        pivot = next(
            row
            for row in range(column, len(work))
            if work[row][column]
        )
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        pivot_value = work[column][column]
        result *= pivot_value
        for index in range(column, len(work)):
            work[column][index] /= pivot_value
        for row in range(column + 1, len(work)):
            factor = work[row][column]
            for index in range(column, len(work)):
                work[row][index] -= factor * work[column][index]
    return result


def main() -> None:
    f, q, generators = witness()
    r = generators["R"]
    z = generators["Z"]
    w = generators["W"]
    t = generators["T"]

    # The four bilinears lie on the rank-one quadric.
    assert add(power(t, 2), scale(-1, power(r, 2)), scale(2, multiply(z, w))) == {}
    assert len(f) == 16
    assert all(
        xi1 + xi2 == z1 + z2 == 4
        for xi1, xi2, z1, z2 in f
    )
    assert determinant(coefficient_matrix(f)) == 48

    f_power = monomial(ZERO)
    sparse_term_counts: dict[str, int] = {}
    for order in range(1, CUTOFF + 1):
        f_power = multiply(f_power, f)
        sparse_term_counts[str(order)] = len(f_power)
        assert contraction(f_power) == {}
        expected = Fraction(
            factorial(4 * order + 2) * factorial(order),
            double_factorial_odd(2 * order + 1),
        )
        assert contraction(multiply(q, f_power)) == {ZERO: expected}

    # Phase extraction followed by the uniform Hopf-coordinate integral
    # gives these two sums.  The written proof evaluates them for all m by
    # differentiating x^(m-1)*integral_0^x (1-y^2)^m dy at x=1.
    for order in range(1, IDENTITY_CUTOFF + 1):
        pure = Fraction(0)
        mixed = Fraction(0)
        for index in range(order + 1):
            common = (
                Fraction((-1) ** index * comb(order, index), 2**order)
                / (2 * index + 1)
            )
            pure += common * comb(order + 2 * index, order)
            mixed += common * comb(order + 2 * index, order - 1)
        assert pure == 0
        assert mixed == Fraction(
            factorial(order),
            double_factorial_odd(2 * order + 1),
        )

    artifact = {
        "format": "two-pair-image-mathieu-counterexample-v1",
        "field": "characteristic zero",
        "contraction_pairs": [["xi1", "z1"], ["xi2", "z2"]],
        "bilinears": {
            "R": "xi1*z1+xi2*z2",
            "Z": "xi1*z2",
            "W": "2*xi2*z1",
            "T": "xi1*z1-xi2*z2",
            "relation": "T^2=R^2-2*Z*W",
        },
        "F": "(R+Z)*(R^2*W-(2*R+Z)*T^2/2)",
        "Q": "Z",
        "expanded_F_term_count": len(f),
        "ordinary_total_degree_F": 8,
        "bidegrees": {"F": [4, 4], "Q": [1, 1]},
        "coefficient_matrix_determinant": 48,
        "all_order_identities": {
            "E_2(F^m)": "0",
            "E_2(Q*F^m)": "(4*m+2)!*m!/(2*m+1)!!",
        },
        "exact_sparse_replay_cutoff": CUTOFF,
        "sparse_power_term_counts": sparse_term_counts,
        "finite_sum_identity_cutoff": IDENTITY_CUTOFF,
        "written_proof": (
            "extended-geometry/"
            "TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print(
        "PASS SIC(2): E_2(F^m)=0 and "
        "E_2(Q*F^m)=(4m+2)!m!/(2m+1)!! through m=8"
    )
    print("PASS SIC(2): all-order finite-sum identities through m=99")
    print("PASS SIC(2): F has 16 terms and bidegree (4,4)")
    print("PASS SIC(2): the 5x5 coefficient matrix has determinant 48")
    print(f"PASS SIC(2): wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
