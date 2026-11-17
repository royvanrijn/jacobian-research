#!/usr/bin/env python3
"""Exact checks for complete equal-radial unions in binary GVC.

Fix an order ``d`` and a radial vector ``rho``.  At scale ``N``, sum
*all* operator and polynomial selection states of order ``N*d`` having
common radial vector ``N*rho``.  Before the common radial factorial is
restored, this complete union is

    [X^(N*rho)] lambda(X)^(N*d)
    [Y^(N*rho)] P(Y)^(N*d)

and hence is the constant term of

    (X^(-rho) Y^(-rho) lambda(X)^d P(Y)^d)^N.

The script verifies this identity on an exact two-dimensional example.
It also records the smallest warning against the stronger, unnecessary
claim that achievable color counts must be face-saturated.  For

    C(z) = a(1+z^2) + b*z

at order ``2*N`` and level ``2*N``, the achievable number of selections
from ``a(1+z^2)`` is exactly the even interval ``0,2,...,2*N``.  Thus
there are holes at every scale.  Nevertheless the whole union is one
constant-term sequence.  Its first two rows are

    b^2 + 2*a^2,
    b^4 + 12*a^2*b^2 + 6*a^4,

and have no common nonzero projective solution.

The identities are exact.  They prove that color-count saturation is
not required once a complete radial-vector union is exposed.  They do
not prove that the Hall/jet filtration always exposes such a union.
"""

from __future__ import annotations

from collections import defaultdict
from math import factorial

import sympy as sp


Exponent = tuple[int, int]


def compositions(total: int, length: int):
    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, length - 1):
            yield (first,) + tail


def multinomial(total: int, entries: tuple[int, ...]) -> int:
    answer = factorial(total)
    for entry in entries:
        answer //= factorial(entry)
    return answer


def coefficient_of_power(
    exponents: tuple[Exponent, ...],
    coefficients: tuple[int, ...],
    order: int,
    target: Exponent,
) -> int:
    answer = 0
    for entries in compositions(order, len(exponents)):
        exponent = (
            sum(entry * alpha[0] for entry, alpha in zip(entries, exponents)),
            sum(entry * alpha[1] for entry, alpha in zip(entries, exponents)),
        )
        if exponent == target:
            coefficient = multinomial(order, entries)
            for entry, value in zip(entries, coefficients):
                coefficient *= value**entry
            answer += coefficient
    return answer


def complete_return_sum(
    operator_exponents: tuple[Exponent, ...],
    operator_coefficients: tuple[int, ...],
    polynomial_exponents: tuple[Exponent, ...],
    polynomial_coefficients: tuple[int, ...],
    base_order: int,
    radial: Exponent,
    scale: int,
) -> int:
    order = base_order * scale
    target = (scale * radial[0], scale * radial[1])
    operator_coefficient = coefficient_of_power(
        operator_exponents,
        operator_coefficients,
        order,
        target,
    )
    polynomial_coefficient = coefficient_of_power(
        polynomial_exponents,
        polynomial_coefficients,
        order,
        target,
    )
    return operator_coefficient * polynomial_coefficient


def direct_complete_return_sum(
    operator_exponents: tuple[Exponent, ...],
    operator_coefficients: tuple[int, ...],
    polynomial_exponents: tuple[Exponent, ...],
    polynomial_coefficients: tuple[int, ...],
    base_order: int,
    radial: Exponent,
    scale: int,
) -> int:
    order = base_order * scale
    target = (scale * radial[0], scale * radial[1])
    operator_states = []
    polynomial_states = []

    for entries in compositions(order, len(operator_exponents)):
        exponent = (
            sum(
                entry * alpha[0]
                for entry, alpha in zip(entries, operator_exponents)
            ),
            sum(
                entry * alpha[1]
                for entry, alpha in zip(entries, operator_exponents)
            ),
        )
        if exponent == target:
            value = multinomial(order, entries)
            for entry, coefficient in zip(entries, operator_coefficients):
                value *= coefficient**entry
            operator_states.append(value)

    for entries in compositions(order, len(polynomial_exponents)):
        exponent = (
            sum(
                entry * beta[0]
                for entry, beta in zip(entries, polynomial_exponents)
            ),
            sum(
                entry * beta[1]
                for entry, beta in zip(entries, polynomial_exponents)
            ),
        )
        if exponent == target:
            value = multinomial(order, entries)
            for entry, coefficient in zip(entries, polynomial_coefficients):
                value *= coefficient**entry
            polynomial_states.append(value)

    return sum(
        operator_value * polynomial_value
        for operator_value in operator_states
        for polynomial_value in polynomial_states
    )


def achievable_red_counts(scale: int) -> tuple[int, ...]:
    """Counts for R={0,2}, B={1}, order 2N, target level 2N."""

    order = 2 * scale
    target = 2 * scale
    counts = set()
    for red_zero in range(order + 1):
        for red_two in range(order - red_zero + 1):
            blue_one = order - red_zero - red_two
            level = 2 * red_two + blue_one
            if level == target:
                counts.add(red_zero + red_two)
    return tuple(sorted(counts))


def first_color_count_hole(maximum_scale: int) -> tuple[int, int]:
    for scale in range(1, maximum_scale + 1):
        counts = achievable_red_counts(scale)
        for count in range(counts[0], counts[-1] + 1):
            if count not in counts:
                return scale, count
    raise AssertionError("expected the parity hole")


def verify() -> None:
    operator_exponents = ((2, 0), (1, 1), (0, 2))
    operator_coefficients = (2, -1, 3)
    polynomial_exponents = ((2, 0), (1, 1), (0, 2), (1, 0), (0, 1))
    polynomial_coefficients = (1, 4, -2, 3, -1)
    base_order = 2
    radial = (2, 2)

    complete_rows = {}
    for scale in range(1, 5):
        factored = complete_return_sum(
            operator_exponents,
            operator_coefficients,
            polynomial_exponents,
            polynomial_coefficients,
            base_order,
            radial,
            scale,
        )
        direct = direct_complete_return_sum(
            operator_exponents,
            operator_coefficients,
            polynomial_exponents,
            polynomial_coefficients,
            base_order,
            radial,
            scale,
        )
        assert direct == factored
        complete_rows[scale] = factored

    first_hole = first_color_count_hole(8)
    assert first_hole == (1, 1)
    parity_rows = {}
    for scale in range(1, 9):
        counts = achievable_red_counts(scale)
        expected = tuple(range(0, 2 * scale + 1, 2))
        assert counts == expected
        parity_rows[scale] = counts

    a, b = sp.symbols("a b")
    z = sp.symbols("z")
    color_polynomial = a * (1 + z**2) + b * z
    row_one = sp.expand(color_polynomial**2).coeff(z, 2)
    row_two = sp.expand(color_polynomial**4).coeff(z, 4)
    assert row_one == b**2 + 2 * a**2
    assert row_two == b**4 + 12 * a**2 * b**2 + 6 * a**4
    assert sp.rem(
        sp.Poly(row_two, b),
        sp.Poly(row_one, b),
    ).as_expr() == -14 * a**4

    print(
        "PASS complete equal-radial union factorization through scale 4: "
        f"{complete_rows}"
    )
    print(
        "first persistent color-count saturation failure: "
        "R={0,2}, B={1}, order=2, level=2, missing red count 1"
    )
    print(f"scaled achievable red counts: {parity_rows}")
    print(
        "whole-union rows: "
        f"M1={row_one}, M2={row_two}, "
        "M2 mod M1=-14*a^4"
    )
    print(
        "STATUS: color-count face saturation is false, but a complete "
        "radial-vector union is one Laurent constant-term sequence; "
        "Hall/jet exposure remains unproved in the parked route and is "
        "not needed by Hall-envelope separation"
    )


if __name__ == "__main__":
    verify()
