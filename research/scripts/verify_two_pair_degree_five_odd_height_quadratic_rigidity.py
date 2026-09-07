#!/usr/bin/env python3
"""Quadratic rigidity of the degree-five odd-height correction sector."""

from __future__ import annotations

import json
import sys
from collections import Counter
from functools import reduce
from math import gcd
from math import comb
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from verify_two_pair_image_mathieu_counterexample import (  # noqa: E402
    ZERO,
    contraction,
    multiply,
    power,
    witness,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_degree_five_odd_height_quadratic_rigidity.json"
)


def standard_monomials_from_basis(
    basis: sp.GroebnerBasis,
    variable_count: int,
    degree_bound: int,
) -> list[tuple[int, ...]]:
    leading_monomials = [
        polynomial.LM(order=basis.order).exponents
        for polynomial in basis.polys
    ]
    standard_monomials: list[tuple[int, ...]] = []

    def visit(prefix: tuple[int, ...], remaining_degree: int) -> None:
        if len(prefix) == variable_count - 1:
            exponent = prefix + (remaining_degree,)
            if not any(
                all(
                    exponent[index] >= leading[index]
                    for index in range(variable_count)
                )
                for leading in leading_monomials
            ):
                standard_monomials.append(exponent)
            return
        for value in range(remaining_degree + 1):
            visit(prefix + (value,), remaining_degree - value)

    for total_degree in range(degree_bound + 1):
        visit((), total_degree)
    return standard_monomials


def main() -> None:
    f, _, generators = witness()
    r = generators["R"]
    z = generators["Z"]
    w = generators["W"]
    t = generators["T"]
    seed = multiply(r, f)

    a2, a4, b2, b4 = sp.symbols("a2 a4 b2 b4")
    variables = (a2, a4, b2, b4)
    corrections = (
        multiply(power(z, 2), power(t, 3)),
        multiply(power(z, 4), t),
        multiply(power(w, 2), power(t, 3)),
        multiply(power(w, 4), t),
    )
    even_height_corrections = tuple(
        multiply(power(base, phase), power(t, 5 - phase))
        for base in (z, w)
        for phase in (1, 3, 5)
    )
    seed_powers = {0: power(seed, 0)}
    for exponent in range(1, 11):
        seed_powers[exponent] = multiply(seed_powers[exponent - 1], seed)
    correction_products = {
        (left, right): multiply(corrections[left], corrections[right])
        for left in range(4)
        for right in range(left, 4)
    }

    quadratic_obstructions: list[sp.Expr] = []
    quadratic_coefficients: list[sp.Expr] = []
    linear_rows: list[list[sp.Rational]] = []
    for order in range(2, 12):
        linear_row = []
        for correction in even_height_corrections:
            value = order * contraction(
                multiply(seed_powers[order - 1], correction)
            ).get(ZERO, 0)
            linear_row.append(sp.Rational(value.numerator, value.denominator))
        linear_rows.append(linear_row)

        obstruction = 0
        for left in range(4):
            for right in range(left, 4):
                contraction_value = contraction(
                    multiply(
                        seed_powers[order - 2],
                        correction_products[(left, right)],
                    )
                ).get(ZERO, 0)
                symmetry = 1 if left == right else 2
                coefficient = (
                    comb(order, 2) * symmetry * contraction_value
                )
                obstruction += (
                    sp.Rational(coefficient.numerator, coefficient.denominator)
                    * variables[left]
                    * variables[right]
                )
        quadratic_coefficients.append(obstruction)
        if order > 10:
            continue
        _, primitive = sp.Poly(
            obstruction,
            *variables,
        ).primitive()
        quadratic_obstructions.append(sp.factor(primitive.as_expr()))

    basis = sp.groebner(
        quadratic_obstructions,
        *variables,
        order="grevlex",
    )
    assert basis.is_zero_dimensional
    nilpotence = {}
    for variable in variables:
        assert basis.reduce(variable**3)[1] == 0
        nilpotence[str(variable)] = 3

    standard_monomials = standard_monomials_from_basis(basis, 4, 6)
    hilbert_vector_counter = Counter(
        sum(exponent) for exponent in standard_monomials
    )
    hilbert_vector = [
        hilbert_vector_counter[degree]
        for degree in range(max(hilbert_vector_counter) + 1)
    ]
    assert hilbert_vector == [1, 4, 1]
    assert len(standard_monomials) == 6

    quadratic_monomials = (
        a2**2,
        a2 * a4,
        a2 * b2,
        a2 * b4,
        a4**2,
        a4 * b2,
        a4 * b4,
        b2**2,
        b2 * b4,
        b4**2,
    )
    coefficient_matrix = sp.Matrix(
        [
            [
                sp.Poly(obstruction, *variables).coeff_monomial(monomial)
                for monomial in quadratic_monomials
            ]
            for obstruction in quadratic_obstructions
        ]
    )
    assert coefficient_matrix.rank() == 9
    socle_functional = coefficient_matrix.nullspace()[0]
    common_denominator = sp.ilcm(
        *[coefficient.q for coefficient in socle_functional]
    )
    integral_functional = [
        int(coefficient * common_denominator)
        for coefficient in socle_functional
    ]
    content = reduce(gcd, (abs(value) for value in integral_functional))
    integral_functional = [
        value // content for value in integral_functional
    ]
    socle_pairing = sp.zeros(4)
    functional_index = 0
    for row in range(4):
        for column in range(row, 4):
            value = integral_functional[functional_index]
            socle_pairing[row, column] = value
            socle_pairing[column, row] = value
            functional_index += 1
    socle_determinant = int(socle_pairing.det())
    assert socle_determinant % 11 == 9

    pivot_matrix = sp.Matrix(linear_rows[:6])
    pivot_determinant = int(pivot_matrix.det())
    assert pivot_determinant
    eliminated_even_correction = (
        -pivot_matrix.inv() * sp.Matrix(quadratic_coefficients[:6])
    )
    projected_obstructions: list[sp.Expr] = []
    for index in range(6, 10):
        projected = sp.factor(
            quadratic_coefficients[index]
            + (
                sp.Matrix([linear_rows[index]])
                * eliminated_even_correction
            )[0]
        )
        _, primitive = sp.Poly(projected, *variables).primitive()
        projected_obstructions.append(sp.factor(primitive.as_expr()))

    projected_basis = sp.groebner(
        projected_obstructions,
        *variables,
        order="grevlex",
    )
    assert projected_basis.is_zero_dimensional
    projected_nilpotence = {}
    for variable in variables:
        assert projected_basis.reduce(variable**5)[1] == 0
        projected_nilpotence[str(variable)] = 5
    projected_standard_monomials = standard_monomials_from_basis(
        projected_basis,
        4,
        8,
    )
    projected_hilbert_counter = Counter(
        sum(exponent) for exponent in projected_standard_monomials
    )
    projected_hilbert_vector = [
        projected_hilbert_counter[degree]
        for degree in range(max(projected_hilbert_counter) + 1)
    ]
    assert projected_hilbert_vector == [1, 4, 6, 4, 1]
    assert len(projected_standard_monomials) == 16

    artifact = {
        "format": "two-pair-degree-five-odd-height-quadratic-rigidity-v1",
        "field": "characteristic zero",
        "seed": "R*F",
        "correction": (
            "a2*Z^2*T^3+a4*Z^4*T+b2*W^2*T^3+b4*W^4*T"
        ),
        "orders": list(range(2, 11)),
        "primitive_quadratic_obstructions": [
            str(obstruction)
            for obstruction in quadratic_obstructions
        ],
        "groebner_basis": [
            str(polynomial.as_expr())
            for polynomial in basis.polys
        ],
        "variable_nilpotence_in_initial_obstruction_ideal": nilpotence,
        "radical": "(a2,a4,b2,b4)",
        "leading_monomial_exponents": [
            polynomial.LM(order=basis.order).exponents
            for polynomial in basis.polys
        ],
        "standard_monomial_exponents": standard_monomials,
        "hilbert_vector": hilbert_vector,
        "length": len(standard_monomials),
        "socle_degree": 2,
        "socle_functional_on_quadratic_monomials": {
            str(monomial): value
            for monomial, value in zip(
                quadratic_monomials,
                integral_functional,
            )
        },
        "socle_pairing_determinant": str(socle_determinant),
        "socle_pairing_determinant_mod_11": socle_determinant % 11,
        "full_monomial_formal_rigidity": {
            "even_height_variables": [
                "a1",
                "a3",
                "a5",
                "b1",
                "b3",
                "b5",
            ],
            "linear_pivot_orders": [2, 3, 4, 5, 6, 7],
            "linear_pivot_determinant": str(pivot_determinant),
            "projected_quadratic_orders": [8, 9, 10, 11],
            "projected_quadratic_obstructions": [
                str(obstruction)
                for obstruction in projected_obstructions
            ],
            "projected_groebner_basis": [
                str(polynomial.as_expr())
                for polynomial in projected_basis.polys
            ],
            "variable_nilpotence": projected_nilpotence,
            "hilbert_vector": projected_hilbert_vector,
            "length": len(projected_standard_monomials),
        },
        "written_source": (
            "extended-geometry/"
            "TWO_PAIR_OPPOSITE_MONOMIAL_OBSTRUCTION.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(
        "PASS degree-five odd-height sector: quadratic obstructions from "
        "moments 2..10 have radical (a2,a4,b2,b4)"
    )
    print(
        "PASS exact Groebner certificate: every correction variable cube "
        "belongs to the quadratic obstruction ideal"
    )
    print(
        "PASS obstruction algebra: Hilbert vector (1,4,1), length 6, "
        "and nondegenerate socle pairing"
    )
    print(
        "PASS full monomial formal rigidity: moments 2..7 eliminate the "
        "six even-height corrections and projected moments 8..11 form "
        "a length-16 complete intersection"
    )
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
