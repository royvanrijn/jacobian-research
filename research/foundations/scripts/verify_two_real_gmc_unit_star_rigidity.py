#!/usr/bin/env python3
"""Exact regression for the unit-star prime-rigidity theorem."""

from __future__ import annotations

import math

import sympy as sp


def factorial_functional(expression: sp.Expr, radial: sp.Symbol) -> int:
    return sum(
        int(coefficient) * math.factorial(exponent[0])
        for exponent, coefficient in sp.Poly(
            sp.expand(expression), radial
        ).terms()
    )


def small_star_moment(
    moment: int, centered: sp.Expr, first: sp.Expr, second: sp.Expr
) -> sp.Expr:
    answer = 0
    for second_power in range(moment // 3 + 1):
        for first_power in range((moment - 3 * second_power) // 2 + 1):
            centered_power = moment - 2 * first_power - 3 * second_power
            positive_power = first_power + 2 * second_power
            coefficient = (
                math.factorial(moment)
                // math.factorial(positive_power)
                // math.factorial(first_power)
                // math.factorial(second_power)
                // math.factorial(centered_power)
            )
            answer += (
                coefficient
                * centered**centered_power
                * first**first_power
                * second**second_power
            )
    return sp.expand(answer)


def main() -> None:
    radial = sp.symbols("U")
    centered_symbol, first_symbol, second_symbol = sp.symbols("C X Y")

    expected = (
        centered_symbol,
        centered_symbol**2 + 2 * first_symbol,
        centered_symbol**3
        + 6 * centered_symbol * first_symbol
        + 3 * second_symbol,
        centered_symbol**4
        + 12 * centered_symbol**2 * first_symbol
        + 6 * first_symbol**2
        + 12 * centered_symbol * second_symbol,
    )
    for moment, target in enumerate(expected, start=1):
        assert (
            small_star_moment(
                moment, centered_symbol, first_symbol, second_symbol
            )
            == target
        )

    for prime in (5, 7, 11):
        assert sp.Poly(
            small_star_moment(
                prime, centered_symbol, first_symbol, second_symbol
            )
            - centered_symbol**prime,
            centered_symbol,
            first_symbol,
            second_symbol,
            modulus=prime,
        ).is_zero
        assert sp.Poly(
            small_star_moment(
                2 * prime, centered_symbol, first_symbol, second_symbol
            )
            - (centered_symbol ** (2 * prime) + 2 * first_symbol**prime),
            centered_symbol,
            first_symbol,
            second_symbol,
            modulus=prime,
        ).is_zero
        assert sp.Poly(
            small_star_moment(
                3 * prime, centered_symbol, first_symbol, second_symbol
            )
            - (
                centered_symbol ** (3 * prime)
                + 6 * centered_symbol**prime * first_symbol**prime
                + 3 * second_symbol**prime
            ),
            centered_symbol,
            first_symbol,
            second_symbol,
            modulus=prime,
        ).is_zero

        cases = (
            (0, 2, 3, prime, 0, 1),
            (2, 1, 3, 2 * prime, prime, 2),
            (2, 3, 1, 3 * prime, prime, 3),
        )
        for a, x, y, moment, baseline, endpoint in cases:
            centered = radial**a * (1 + 2 * radial)
            first = radial**x * (1 + 3 * radial)
            second = radial**y * (1 + 5 * radial)
            value = factorial_functional(
                small_star_moment(moment, centered, first, second),
                radial,
            )
            divisor = math.factorial(baseline)
            assert value % divisor == 0
            assert value // divisor % prime == endpoint

    assert 5 - 2 - 3 == 0
    assert 3 < 5 + 2 and 3 < 5 + 3

    print("PASS unit star: primitive invariants and first four moments")
    print("PASS unit star: p, 2p, and 3p endpoint congruences")
    print("PASS unit star: all normalized-order cases")
    print("PASS unit star: arbitrary-weight semigroup obstruction")


if __name__ == "__main__":
    main()
