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

For every d >= 4, the propagated witness F_d=R^(d-4)*F satisfies

    E_2(F_d^m) = 0,
    E_2(Q*F_d^m) = (d*m+2)!*m!/(2*m+1)!!.

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
PROPAGATION_DEGREES = range(4, 11)
PROPAGATION_CUTOFF = 4
BOUNDED_RADIAL_DEGREES = range(4, 16)
BOUNDED_RADIAL_CUTOFF = 3
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


def evaluate(polynomial: Polynomial, point: Exponent) -> Fraction:
    result = Fraction(0)
    for exponent, coefficient in polynomial.items():
        term = coefficient
        for value, order in zip(point, exponent):
            term *= value**order
        result += term
    return result


def hopf_profile_witness(
    height: int,
    r: Polynomial,
    z: Polynomial,
    w: Polynomial,
    t: Polynomial,
) -> Polynomial:
    """Return Phi_h for rho_h(u)=(1-u)(1+u)^(h-1)."""
    r_plus_z = add(r, z)
    d = add(
        scale(2, multiply(w, power(r, 2))),
        scale(-2, multiply(power(t, 2), r)),
        scale(-1, multiply(z, power(t, 2))),
    )
    profile_sum: Polynomial = {}
    for index in range(height):
        term = multiply(
            multiply(
                power(r, 4 * (height - 1 - index)),
                power(t, 2 * index),
            ),
            power(r_plus_z, 2 * index),
        )
        profile_sum = add(
            profile_sum,
            scale(comb(height - 1, index), term),
        )
    return multiply(multiply(r_plus_z, d), profile_sum)


def hopf_profile_integral(height: int, order: int) -> Fraction:
    """Integral of (1-v^2)^m(1+v^2)^((h-1)m) on [0,1]."""
    result = Fraction(0)
    for left in range(order + 1):
        for right in range((height - 1) * order + 1):
            result += Fraction(
                (-1) ** left
                * comb(order, left)
                * comb((height - 1) * order, right),
                2 * (left + right) + 1,
            )
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
    radial_zero_point = (1, 1, 1, -1)
    assert evaluate(r, radial_zero_point) == 0
    assert evaluate(f, radial_zero_point) == -2

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

    propagation_term_counts: dict[str, dict[str, int]] = {}
    for degree in PROPAGATION_DEGREES:
        f_degree = multiply(power(r, degree - 4), f)
        assert all(
            xi1 + xi2 == z1 + z2 == degree
            for xi1, xi2, z1, z2 in f_degree
        )
        degree_power = monomial(ZERO)
        propagation_term_counts[str(degree)] = {}
        for order in range(1, PROPAGATION_CUTOFF + 1):
            degree_power = multiply(degree_power, f_degree)
            propagation_term_counts[str(degree)][str(order)] = len(
                degree_power
            )
            assert contraction(degree_power) == {}
            expected = Fraction(
                factorial(degree * order + 2) * factorial(order),
                double_factorial_odd(2 * order + 1),
            )
            assert contraction(multiply(q, degree_power)) == {
                ZERO: expected
            }

    bounded_radial_term_counts: dict[str, dict[str, int]] = {}
    for degree in BOUNDED_RADIAL_DEGREES:
        seed_power, radial_order = divmod(degree, 4)
        f_degree = multiply(power(r, radial_order), power(f, seed_power))
        assert all(
            xi1 + xi2 == z1 + z2 == degree
            for xi1, xi2, z1, z2 in f_degree
        )
        bounded_radial_term_counts[str(degree)] = {}
        degree_power = monomial(ZERO)
        for order in range(1, BOUNDED_RADIAL_CUTOFF + 1):
            degree_power = multiply(degree_power, f_degree)
            bounded_radial_term_counts[str(degree)][str(order)] = len(
                degree_power
            )
            assert contraction(degree_power) == {}
            expected = Fraction(
                factorial(degree * order + 2)
                * factorial(seed_power * order),
                double_factorial_odd(2 * seed_power * order + 1),
            )
            assert contraction(multiply(q, degree_power)) == {
                ZERO: expected
            }

    hopf_profile_term_counts: dict[str, dict[str, int]] = {}
    for height in range(1, 5):
        profile_witness = hopf_profile_witness(height, r, z, w, t)
        degree = 4 * height
        assert all(
            xi1 + xi2 == z1 + z2 == degree
            for xi1, xi2, z1, z2 in profile_witness
        )
        assert evaluate(profile_witness, radial_zero_point) != 0
        if height == 1:
            assert profile_witness == scale(2, f)
        hopf_profile_term_counts[str(height)] = {}
        profile_power = monomial(ZERO)
        for order in range(1, 4):
            profile_power = multiply(profile_power, profile_witness)
            hopf_profile_term_counts[str(height)][str(order)] = len(
                profile_power
            )
            assert contraction(profile_power) == {}
            integral = hopf_profile_integral(height, order)
            assert integral > 0
            expected = factorial(degree * order + 2) * integral
            assert contraction(multiply(q, profile_power)) == {
                ZERO: expected
            }

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
            "F_d": "R^(d-4)*F for every d>=4",
            "E_2(F_d^m)": "0",
            "E_2(Q*F_d^m)": "(d*m+2)!*m!/(2*m+1)!!",
            "G_(r,k)": "R^k*F^r of degree d=4r+k",
            "E_2(G_(r,k)^m)": "0",
            "E_2(Q*G_(r,k)^m)": "(d*m+2)!*(r*m)!/(2*r*m+1)!!",
            "Phi_h_degree": "4*h",
            "E_2(Phi_h^m)": "0",
            "E_2(Q*Phi_h^m)": (
                "(4*h*m+2)!*integral_0^1 "
                "(1-v^2)^m*(1+v^2)^((h-1)*m)dv"
            ),
        },
        "exact_sparse_replay_cutoff": CUTOFF,
        "sparse_power_term_counts": sparse_term_counts,
        "finite_sum_identity_cutoff": IDENTITY_CUTOFF,
        "propagation_replay": {
            "degrees": [
                min(PROPAGATION_DEGREES),
                max(PROPAGATION_DEGREES),
            ],
            "moment_cutoff": PROPAGATION_CUTOFF,
            "sparse_power_term_counts": propagation_term_counts,
        },
        "bounded_radial_order_replay": {
            "degrees": [
                min(BOUNDED_RADIAL_DEGREES),
                max(BOUNDED_RADIAL_DEGREES),
            ],
            "moment_cutoff": BOUNDED_RADIAL_CUTOFF,
            "sparse_power_term_counts": bounded_radial_term_counts,
        },
        "non_power_hopf_profile_replay": {
            "heights": [1, 4],
            "moment_cutoff": 3,
            "sparse_power_term_counts": hopf_profile_term_counts,
        },
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
    print(
        "PASS SIC(2): bounded replay of the written all-order finite-sum "
        "formulas through m=99"
    )
    print(
        "PASS MN_d: F_d=R^(d-4)F has the claimed pure and mixed "
        "contractions for 4<=d<=10 and m<=4"
    )
    print(
        "PASS MN_d: R^(d mod 4)F^floor(d/4) has the claimed moments "
        "for 4<=d<=15 and m<=3"
    )
    print(
        "PASS MN_(4h): primitive non-power Hopf profiles have the "
        "claimed moments for 1<=h<=4 and m<=3"
    )
    print("PASS SIC(2): F has 16 terms and bidegree (4,4)")
    print("PASS SIC(2): the 5x5 coefficient matrix has determinant 48")
    print(f"PASS SIC(2): wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
