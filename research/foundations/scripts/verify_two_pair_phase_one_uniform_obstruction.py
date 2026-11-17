#!/usr/bin/env python3
"""Exact symbolic certificate for the uniform opposite phase-one pair."""

from __future__ import annotations

import json
from math import factorial
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_phase_one_uniform_obstruction.json"
)
h, u, a, b = sp.symbols("h u a b")
PROFILE = {
    -1: (1 - u) / 2,
    0: (1 - 3 * u) / 2,
    1: -3 * u / 2,
    2: -u / 2,
}


def convolve(
    left: dict[int, sp.Expr],
    right: dict[int, sp.Expr],
) -> dict[int, sp.Expr]:
    result: dict[int, sp.Expr] = {}
    for left_phase, left_coefficient in left.items():
        for right_phase, right_coefficient in right.items():
            phase = left_phase + right_phase
            result[phase] = sp.expand(
                result.get(phase, 0)
                + left_coefficient * right_coefficient
            )
    return result


def profile_power(order: int) -> dict[int, sp.Expr]:
    result: dict[int, sp.Expr] = {0: sp.Integer(1)}
    for _ in range(order):
        result = convolve(result, PROFILE)
    return result


def height_integral(polynomial: sp.Expr, base_power: sp.Expr) -> sp.Expr:
    expanded = sp.Poly(sp.expand(polynomial), u)
    return sp.factor(
        sum(
            coefficient / (base_power + 2 * power[0] + 1)
            for power, coefficient in expanded.terms()
        )
    )


def angular_term(
    seed_count: int,
    positive_count: int,
    negative_count: int,
) -> sp.Expr:
    required_seed_phase = negative_count - positive_count
    coefficient = profile_power(seed_count).get(required_seed_phase, 0)
    if coefficient == 0:
        return sp.Integer(0)
    return height_integral(
        coefficient * (1 - u) ** negative_count / 2**negative_count,
        2 * h * (positive_count + negative_count),
    )


def angular_moment(order: int) -> sp.Expr:
    result = 0
    for positive_count in range(order + 1):
        for negative_count in range(order - positive_count + 1):
            seed_count = order - positive_count - negative_count
            multinomial = sp.Rational(
                factorial(order),
                factorial(seed_count)
                * factorial(positive_count)
                * factorial(negative_count),
            )
            result += (
                multinomial
                * a**positive_count
                * b**negative_count
                * angular_term(
                    seed_count,
                    positive_count,
                    negative_count,
                )
            )
    return sp.factor(result)


def coefficient_table(
    polynomial: sp.Expr,
    variable: sp.Symbol,
) -> list[list[str]]:
    poly = sp.Poly(polynomial, variable)
    return [
        [str(power[0]), str(coefficient)]
        for power, coefficient in poly.terms()
    ]


def main() -> None:
    moment_2, moment_3, moment_4 = (
        angular_moment(order) for order in (2, 3, 4)
    )
    a_solution = sp.factor(sp.solve(moment_2, a)[0])
    expected_a_solution = sp.factor(
        3
        * b
        * (2 * h + 1)
        * (4 * h + 1)
        * (4 * h + 3)
        / (
            2
            * (2 * h + 5)
            * (
                4 * b * h**2
                + 8 * b * h
                + 3 * b
                + 16 * h**2
                + 16 * h
                + 3
            )
        )
    )
    assert sp.factor(a_solution - expected_a_solution) == 0

    numerator_3 = sp.factor(
        sp.together(moment_3.subs(a, a_solution)).as_numer_denom()[0]
    )
    numerator_4 = sp.factor(
        sp.together(moment_4.subs(a, a_solution)).as_numer_denom()[0]
    )
    factor_3 = sp.factor_list(numerator_3)[1]
    factor_4 = sp.factor_list(numerator_4)[1]
    quadratic = max(
        (factor for factor, _ in factor_3),
        key=lambda factor: sp.degree(factor, b),
    )
    cubic = max(
        (factor for factor, _ in factor_4),
        key=lambda factor: sp.degree(factor, b),
    )
    assert sp.degree(quadratic, b) == 2
    assert sp.degree(cubic, b) == 3

    elimination_resultant = sp.factor(sp.resultant(quadratic, cubic, b))
    resultant_factorization = sp.factor_list(elimination_resultant)
    degree_31_factor = max(
        (factor for factor, _ in resultant_factorization[1]),
        key=lambda factor: sp.degree(factor, h),
    )
    assert sp.degree(degree_31_factor, h) == 31

    expected_linear_factors = {
        (str(sp.factor(factor)), multiplicity)
        for factor, multiplicity in resultant_factorization[1]
        if sp.degree(factor, h) <= 1
    }
    assert expected_linear_factors == {
        ("4*h + 5", 1),
        ("h", 2),
        ("2*h + 1", 2),
        ("2*h + 3", 2),
        ("2*h + 7", 2),
        ("4*h + 1", 2),
        ("4*h + 7", 2),
        ("4*h + 3", 6),
    }

    modulus = 29
    residues = [
        int(degree_31_factor.subs(h, residue)) % modulus
        for residue in range(modulus)
    ]
    assert all(residues)

    artifact = {
        "format": "two-pair-phase-one-uniform-obstruction-v1",
        "field": "characteristic zero",
        "family": (
            "R^(d-4)F+a*Z*T^(d-1)+b*W*T^(d-1), "
            "d=2h+1 and h>=2"
        ),
        "moments_used": [2, 3, 4],
        "a_after_moment_2": str(a_solution),
        "moment_3_quadratic_in_b": coefficient_table(quadratic, b),
        "moment_4_cubic_in_b": coefficient_table(cubic, b),
        "resultant_content": str(resultant_factorization[0]),
        "resultant_linear_factors": [
            [factor, int(multiplicity)]
            for factor, multiplicity in sorted(expected_linear_factors)
        ],
        "degree_31_factor_coefficients_descending": [
            str(coefficient)
            for coefficient in sp.Poly(degree_31_factor, h).all_coeffs()
        ],
        "modulus": modulus,
        "nonzero_residues": residues,
        "written_source": (
            "extended-geometry/"
            "TWO_PAIR_OPPOSITE_MONOMIAL_OBSTRUCTION.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(
        "PASS uniform phase one: moments 2,3,4 reduce to a "
        "quadratic/cubic pair with nonzero resultant for every h>=2"
    )
    print(
        "PASS resultant: its only nonlinear factor has degree 31 and "
        "no root modulo 29"
    )
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
