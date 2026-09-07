#!/usr/bin/env python3
"""Bounded exact search on the fixed-quintic trace descent variety.

For

    A = Q[theta]/(theta^5 + u*theta^3 + v*theta + v),

enumerate small trace-zero elements eta and test the intrinsic square/cube
conditions from FIXED_QUINTIC_MODULI_DOMINANCE.md.  This is an experiment,
not a proof of rational solubility.
"""

from __future__ import annotations

import argparse
import itertools
import math
from fractions import Fraction
from functools import reduce

import sympy as sp


def rational_square_root(value: Fraction) -> Fraction | None:
    if value <= 0:
        return None
    numerator_root = math.isqrt(value.numerator)
    denominator_root = math.isqrt(value.denominator)
    if (
        numerator_root * numerator_root == value.numerator
        and denominator_root * denominator_root == value.denominator
    ):
        return Fraction(numerator_root, denominator_root)
    return None


def integer_cube_root(value: int) -> int | None:
    sign = -1 if value < 0 else 1
    absolute = abs(value)
    candidate, exact = sp.integer_nthroot(absolute, 3)
    if exact:
        return sign * int(candidate)
    return None


def rational_cube_root(value: Fraction) -> Fraction | None:
    numerator_root = integer_cube_root(value.numerator)
    denominator_root = integer_cube_root(value.denominator)
    if numerator_root is None or denominator_root is None:
        return None
    return Fraction(numerator_root, denominator_root)


def evaluate_terms(
    terms: list[tuple[tuple[int, ...], int]], coefficients: tuple[int, ...]
) -> int:
    return sum(
        coefficient
        * math.prod(
            value**exponent
            for value, exponent in zip(coefficients, monomial, strict=True)
        )
        for monomial, coefficient in terms
    )


def polynomial_terms(expression: sp.Expr, variables: tuple[sp.Symbol, ...]):
    polynomial = sp.Poly(expression, *variables, domain=sp.ZZ)
    return [(monomial, int(coefficient)) for monomial, coefficient in polynomial.terms()]


def companion_matrix(u: int, v: int) -> sp.Matrix:
    matrix = sp.zeros(5)
    for column in range(4):
        matrix[column + 1, column] = 1
    matrix[:, 4] = sp.Matrix((-v, -v, 0, -u, 0))
    return matrix


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--u", type=int, required=True)
    parser.add_argument("--v", type=int, required=True)
    parser.add_argument("--bound", type=int, default=12)
    parser.add_argument("--max-results", type=int, default=1)
    args = parser.parse_args()

    matrix = companion_matrix(args.u, args.v)
    variables = sp.symbols("c1:5")
    trace_zero_basis = [
        5 * matrix**power - sp.trace(matrix**power) * sp.eye(5)
        for power in range(1, 5)
    ]
    numerator_operator = sum(
        (
            coefficient * basis
            for coefficient, basis in zip(
                variables, trace_zero_basis, strict=True
            )
        ),
        sp.zeros(5),
    )
    trace2_numerator = sp.expand(sp.trace(numerator_operator**2))
    trace4_numerator = sp.expand(sp.trace(numerator_operator**4))
    trace2_terms = polynomial_terms(trace2_numerator, variables)
    trace4_terms = polynomial_terms(trace4_numerator, variables)

    results = 0
    tested = 0
    coefficient_range = range(-args.bound, args.bound + 1)
    for coefficients in itertools.product(coefficient_range, repeat=4):
        if coefficients == (0, 0, 0, 0):
            continue
        if reduce(math.gcd, coefficients) != 1:
            continue
        if next(value for value in coefficients if value != 0) < 0:
            continue
        tested += 1

        trace2_scaled = evaluate_terms(trace2_terms, coefficients)
        normalization = rational_square_root(Fraction(trace2_scaled, 250))
        if normalization is None:
            continue
        trace4_scaled = evaluate_terms(trace4_terms, coefficients)
        cube_class = Fraction(
            25 * (trace2_scaled**2 - 2 * trace4_scaled),
            8 * trace2_scaled**2,
        )
        pi = rational_cube_root(cube_class)
        if pi in (None, 0):
            continue

        raw_coefficients = [
            Fraction(
                sum(
                    coefficients[index] * int(trace_zero_basis[index][row, 0])
                    for index in range(4)
                ),
                5,
            )
            for row in range(5)
        ]
        # The first column is eta(1), hence its coordinates in
        # 1,theta,...,theta^4.
        normalized_coefficients = [
            coefficient / normalization for coefficient in raw_coefficients
        ]
        normalized_operator = numerator_operator.subs(
            dict(zip(variables, coefficients, strict=True))
        ) / (5 * sp.Rational(normalization.numerator, normalization.denominator))
        characteristic = sp.Poly(
            normalized_operator.charpoly().as_expr(), domain=sp.QQ
        )
        characteristic_coefficients = characteristic.all_coeffs()
        target_b = -characteristic_coefficients[3] / (2 * sp.Rational(pi))
        target_c = -characteristic_coefficients[5] / (2 * sp.Rational(pi) ** 5)
        print(
            {
                "u": args.u,
                "v": args.v,
                "search_coefficients": coefficients,
                "eta_coefficients": [str(value) for value in normalized_coefficients],
                "Pi": str(pi),
                "B": str(target_b),
                "C": str(target_c),
                "characteristic_polynomial": str(characteristic.as_expr()),
                "trace2": "10",
                "trace4": str(Fraction(50) - 16 * pi**3),
            }
        )
        results += 1
        if results >= args.max_results:
            break

    print(
        {
            "tested_primitive_vectors": tested,
            "bound": args.bound,
            "results": results,
        }
    )


if __name__ == "__main__":
    main()
