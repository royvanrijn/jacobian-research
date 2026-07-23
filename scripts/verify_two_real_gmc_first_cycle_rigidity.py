#!/usr/bin/env python3
"""Exact regression for prime-valuation rigidity of {-2,-1,1,2}.

The checker enumerates constant terms in the four circuit invariants,
verifies K_(kp) = K_k^p modulo p in the toric normal form, and checks the
finite initial-form eliminations for every possible valuation tie.
"""

from __future__ import annotations

import math
from collections import defaultdict
from fractions import Fraction

import sympy as sp


X11, X12, X21, X22 = sp.symbols("X11 X12 X21 X22")


def canonical_exponents(
    positive_two: int,
    positive_one: int,
    negative_one: int,
    negative_two: int,
) -> tuple[int, int, int, int]:
    """Reduce a zero-weight count vector by X12*X21=X11^2*X22."""

    difference = (negative_one - positive_one) // 2
    assert (negative_one - positive_one) % 2 == 0
    if difference >= 0:
        exponents = (positive_one, 0, difference, negative_two)
    else:
        exponents = (
            negative_one,
            -difference,
            0,
            negative_two + difference,
        )
    assert min(exponents) >= 0
    return exponents


def moment_dictionary(moment: int) -> dict[tuple[int, int, int, int], int]:
    """Return CT(P^moment) in the canonical toric monomial basis."""

    terms: defaultdict[tuple[int, int, int, int], int] = defaultdict(int)
    for positive_two in range(moment + 1):
        for positive_one in range(moment - positive_two + 1):
            for negative_one in range(
                moment - positive_two - positive_one + 1
            ):
                negative_two = (
                    moment - positive_two - positive_one - negative_one
                )
                weight = (
                    2 * positive_two
                    + positive_one
                    - negative_one
                    - 2 * negative_two
                )
                if weight != 0:
                    continue
                coefficient = (
                    math.factorial(moment)
                    // math.factorial(positive_two)
                    // math.factorial(positive_one)
                    // math.factorial(negative_one)
                    // math.factorial(negative_two)
                )
                terms[
                    canonical_exponents(
                        positive_two,
                        positive_one,
                        negative_one,
                        negative_two,
                    )
                ] += coefficient
    return dict(terms)


def as_expression(
    terms: dict[tuple[int, int, int, int], int],
) -> sp.Expr:
    variables = (X11, X12, X21, X22)
    return sp.expand(
        sum(
            coefficient
            * sp.prod(
                variable**exponent
                for variable, exponent in zip(variables, exponents)
            )
            for exponents, coefficient in terms.items()
        )
    )


def modulo_prime(
    terms: dict[tuple[int, int, int, int], int], prime: int
) -> dict[tuple[int, int, int, int], int]:
    return {
        exponents: coefficient % prime
        for exponents, coefficient in terms.items()
        if coefficient % prime
    }


def frobenius_power(
    terms: dict[tuple[int, int, int, int], int], prime: int
) -> dict[tuple[int, int, int, int], int]:
    return {
        tuple(prime * exponent for exponent in exponents):
        pow(coefficient, prime, prime)
        for exponents, coefficient in terms.items()
        if coefficient % prime
    }


def main() -> None:
    moments = {
        moment: as_expression(moment_dictionary(moment))
        for moment in (2, 3, 4, 6, 12)
    }

    assert moments[2] == 2 * (X11 + X22)
    assert moments[3] == 3 * (X12 + X21)
    assert moments[4] == 6 * (
        X11**2 + 4 * X11 * X22 + X22**2
    )
    assert all(
        sp.expand(
            moments[moment].subs(
                {X12: X21, X21: X12}, simultaneous=True
            )
            - moments[moment]
        )
        == 0
        for moment in (6, 12)
    )

    for prime in (3, 5, 7):
        for base_moment in (2, 3, 4, 6):
            assert modulo_prime(
                moment_dictionary(base_moment * prime), prime
            ) == frobenius_power(
                moment_dictionary(base_moment), prime
            )

    # The toric valuation equation permits only the stated minimum faces.
    adjacent_faces = {
        frozenset((0, 1)),
        frozenset((0, 2)),
        frozenset((3, 1)),
        frozenset((3, 2)),
    }
    for alpha_11 in range(13):
        for alpha_12 in range(13):
            for alpha_21 in range(13):
                for alpha_22 in range(13):
                    if alpha_12 + alpha_21 != (
                        2 * alpha_11 + alpha_22
                    ):
                        continue
                    slopes = (
                        Fraction(alpha_11, 2),
                        Fraction(alpha_12, 3),
                        Fraction(alpha_21, 3),
                        Fraction(alpha_22, 2),
                    )
                    minimum = min(slopes)
                    face = frozenset(
                        index
                        for index, slope in enumerate(slopes)
                        if slope == minimum
                    )
                    assert (
                        len(face) == 1
                        or face in adjacent_faces
                        or len(face) == 4
                    )

    x, y = sp.symbols("x y")
    initial_six = 20 * x**3 + 15 * y**2
    assert moments[6].subs(
        {X11: x, X12: y, X21: 0, X22: 0}
    ) == initial_six
    assert moments[6].subs(
        {X11: 0, X12: y, X21: 0, X22: x}
    ) == initial_six

    initial_twelve_inner = moments[12].subs(
        {X11: x, X12: y, X21: 0, X22: 0}
    )
    initial_twelve_outer = moments[12].subs(
        {X11: 0, X12: y, X21: 0, X22: x}
    )
    tie_equation = 3 * y**2 + 4 * x**3
    assert sp.rem(initial_twelve_inner, tie_equation, y) == (
        -8756 * x**6
    )
    assert sp.rem(initial_twelve_outer, tie_equation, y) == (
        -35156 * x**6
    )

    a, b, c, d = sp.symbols("a b c d")
    assert moments[4].subs(
        {X11: a, X12: 0, X21: 0, X22: -a}
    ) == -12 * a**2
    assert moments[6].subs(
        {X11: 0, X12: b, X21: -b, X22: 0}
    ) == 30 * b**2
    assert sp.expand(b * c - a**2 * d).subs(
        {d: -a, c: -b}
    ) == a**3 - b**2

    print("PASS first cycle: exact toric moments through order twelve")
    print("PASS first cycle: K_(kp) equals K_k^p modulo p")
    print("PASS first cycle: every valuation-tie initial system is empty")


if __name__ == "__main__":
    main()
