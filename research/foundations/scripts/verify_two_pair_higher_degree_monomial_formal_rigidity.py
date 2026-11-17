#!/usr/bin/env python3
"""Exact Hopf-angular formal-rigidity certificates in degrees six and seven."""

from __future__ import annotations

import json
from collections import Counter
from math import comb
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_higher_degree_monomial_formal_rigidity.json"
)
u = sp.symbols("u")
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


PROFILE_POWERS = [{0: sp.Integer(1)}]
for _ in range(24):
    PROFILE_POWERS.append(convolve(PROFILE_POWERS[-1], PROFILE))


def height_integral(polynomial: sp.Expr, base_power: int) -> sp.Expr:
    if polynomial == 0:
        return sp.Integer(0)
    return sp.factor(
        sum(
            coefficient / (base_power + 2 * exponent[0] + 1)
            for exponent, coefficient in sp.Poly(
                sp.expand(polynomial),
                u,
            ).terms()
        )
    )


def correction(
    side: str,
    phase: int,
    degree: int,
) -> tuple[int, int, sp.Expr]:
    height = degree - phase
    if side == "positive":
        return phase, height, sp.Integer(1)
    return -phase, height, (1 - u) ** phase / 2**phase


def linear_term(
    order: int,
    term: tuple[int, int, sp.Expr],
) -> sp.Expr:
    phase, height, coefficient = term
    return height_integral(
        PROFILE_POWERS[order - 1].get(-phase, 0) * coefficient,
        height,
    )


def quadratic_term(
    order: int,
    left: tuple[int, int, sp.Expr],
    right: tuple[int, int, sp.Expr],
) -> sp.Expr:
    left_phase, left_height, left_coefficient = left
    right_phase, right_height, right_coefficient = right
    return height_integral(
        PROFILE_POWERS[order - 2].get(
            -(left_phase + right_phase),
            0,
        )
        * left_coefficient
        * right_coefficient,
        left_height + right_height,
    )


def standard_monomials(
    basis: sp.GroebnerBasis,
    variable_count: int,
    degree_bound: int,
) -> list[tuple[int, ...]]:
    leading_monomials = [
        polynomial.LM(order=basis.order).exponents
        for polynomial in basis.polys
    ]
    result: list[tuple[int, ...]] = []

    def visit(prefix: tuple[int, ...], remaining: int) -> None:
        if len(prefix) == variable_count - 1:
            exponent = prefix + (remaining,)
            if not any(
                all(
                    exponent[index] >= leading[index]
                    for index in range(variable_count)
                )
                for leading in leading_monomials
            ):
                result.append(exponent)
            return
        for value in range(remaining + 1):
            visit(prefix + (value,), remaining - value)

    for total_degree in range(degree_bound + 1):
        visit((), total_degree)
    return result


def degree_certificate(degree: int) -> dict[str, object]:
    even_phases = [
        phase
        for phase in range(1, degree + 1)
        if (degree - phase) % 2 == 0
    ]
    odd_phases = [
        phase
        for phase in range(1, degree + 1)
        if (degree - phase) % 2 == 1
    ]
    even_terms = [
        correction(side, phase, degree)
        for side in ("positive", "negative")
        for phase in even_phases
    ]
    odd_terms = [
        correction(side, phase, degree)
        for side in ("positive", "negative")
        for phase in odd_phases
    ]
    variables = sp.symbols(f"x0:{len(odd_terms)}")
    linear_count = len(even_terms)
    kernel_count = len(odd_terms)
    orders = list(range(2, 2 + linear_count + kernel_count))

    linear_rows: list[list[sp.Expr]] = []
    quadratic_coefficients: list[sp.Expr] = []
    for order in orders:
        linear_rows.append(
            [
                sp.factor(order * linear_term(order, term))
                for term in even_terms
            ]
        )
        quadratic = 0
        for left in range(kernel_count):
            for right in range(left, kernel_count):
                quadratic += (
                    sp.binomial(order, 2)
                    * (1 if left == right else 2)
                    * quadratic_term(
                        order,
                        odd_terms[left],
                        odd_terms[right],
                    )
                    * variables[left]
                    * variables[right]
                )
        quadratic_coefficients.append(sp.factor(quadratic))

    pivot_matrix = sp.Matrix(linear_rows[:linear_count])
    pivot_determinant = sp.factor(pivot_matrix.det())
    assert pivot_determinant
    eliminated_correction = (
        -pivot_matrix.inv()
        * sp.Matrix(quadratic_coefficients[:linear_count])
    )

    projected_obstructions = []
    for index in range(linear_count, linear_count + kernel_count):
        projected = sp.factor(
            quadratic_coefficients[index]
            + (
                sp.Matrix([linear_rows[index]])
                * eliminated_correction
            )[0]
        )
        _, primitive = sp.Poly(projected, *variables).primitive()
        projected_obstructions.append(primitive.as_expr())

    basis = sp.groebner(
        projected_obstructions,
        *variables,
        order="grevlex",
    )
    assert basis.is_zero_dimensional
    for variable in variables:
        assert basis.reduce(variable ** (kernel_count + 1))[1] == 0

    monomials = standard_monomials(
        basis,
        kernel_count,
        kernel_count + 2,
    )
    hilbert_counter = Counter(sum(exponent) for exponent in monomials)
    hilbert_vector = [
        hilbert_counter[degree_index]
        for degree_index in range(max(hilbert_counter) + 1)
    ]
    expected_hilbert_vector = [
        comb(kernel_count, degree_index)
        for degree_index in range(kernel_count + 1)
    ]
    assert hilbert_vector == expected_hilbert_vector
    assert len(monomials) == 2**kernel_count

    return {
        "degree": degree,
        "even_height_phases": even_phases,
        "odd_height_phases": odd_phases,
        "linear_pivot_orders": orders[:linear_count],
        "linear_pivot_determinant": str(pivot_determinant),
        "projected_quadratic_orders": orders[
            linear_count : linear_count + kernel_count
        ],
        "projected_quadratic_obstructions": [
            str(obstruction)
            for obstruction in projected_obstructions
        ],
        "leading_monomial_exponents": [
            polynomial.LM(order=basis.order).exponents
            for polynomial in basis.polys
        ],
        "variable_nilpotence_exponent": kernel_count + 1,
        "hilbert_vector": hilbert_vector,
        "length": len(monomials),
    }


def main() -> None:
    certificates = {
        str(degree): degree_certificate(degree)
        for degree in (6, 7)
    }
    artifact = {
        "format": "two-pair-higher-degree-monomial-formal-rigidity-v1",
        "field": "characteristic zero",
        "degrees": certificates,
        "written_source": (
            "extended-geometry/"
            "TWO_PAIR_OPPOSITE_MONOMIAL_OBSTRUCTION.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(
        "PASS degrees 6 and 7: consecutive linear moment blocks eliminate "
        "all even-height corrections"
    )
    print(
        "PASS projected odd-height obstructions: six quadrics form a "
        "complete intersection with Hilbert vector (1,6,15,20,15,6,1)"
    )
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
