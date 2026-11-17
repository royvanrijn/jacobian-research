#!/usr/bin/env python3
"""Bounded exact search for a quartic three-family fiber collision.

The degree-four cancellation inverse is

    F(T)=T-Q^2*T^2/2+2*Pi*Q*T^3/3-Pi^2*T^4/4-R.

Write Q=1, Pi=1/rho, and make the affine generator change
T=rho*(a+d*W).  The weighted tangent-chord condition is the conic

    y^2=-18*a^2+24*a-2.

The rational parameter k through (a,y)=(1,2) gives

    a=(k^2-4*k+6)/(k^2+18).

This script enumerates small k, rho, both chord branches, and four small
constant terms.  It retains boundary-clean weighted presentations,
admissible quadratic-gauge presentations, and quartics irreducible modulo
17.  Since 17 splits in Q(sqrt(-2)), the last test also certifies
irreducibility over the cancellation map's coefficient field.
"""

from __future__ import annotations

from fractions import Fraction
from functools import reduce
from math import gcd, lcm

import sympy as sp


W = sp.symbols("W")


def bounded_fractions(numerator_bound: int, denominator_bound: int) -> list[Fraction]:
    """Return the nonzero reduced fractions in the declared box."""

    return sorted(
        {
            Fraction(numerator, denominator)
            for denominator in range(1, denominator_bound + 1)
            for numerator in range(-numerator_bound, numerator_bound + 1)
            if numerator
        }
    )


def primitive_coefficients(coefficients: tuple[Fraction, ...]) -> tuple[int, ...]:
    """Clear denominators and content, with positive leading coefficient."""

    denominator = 1
    for coefficient in coefficients:
        denominator = lcm(denominator, coefficient.denominator)
    integers = [
        coefficient.numerator * (denominator // coefficient.denominator)
        for coefficient in coefficients
    ]
    content = reduce(gcd, (abs(value) for value in integers if value))
    integers = [value // content for value in integers]
    if integers[0] < 0:
        integers = [-value for value in integers]
    return tuple(integers)


def irreducible_mod_17(coefficients: tuple[int, ...]) -> bool:
    """Test the quartic certificate at a split prime of Q(sqrt(-2))."""

    if coefficients[0] % 17 == 0:
        return False
    polynomial = sum(
        coefficient * W ** (4 - index)
        for index, coefficient in enumerate(coefficients)
    )
    return bool(sp.Poly(polynomial, W, modulus=17).is_irreducible)


def search() -> list[tuple[tuple[int, int], tuple[int, ...], tuple[Fraction, ...]]]:
    """Return all certified candidates, ordered by primitive coefficient size."""

    candidates = []
    seen_seeds: set[tuple[Fraction, Fraction, Fraction]] = set()
    constants = (
        Fraction(1),
        Fraction(-1),
        Fraction(1, 2),
        Fraction(-1, 2),
    )

    for k in bounded_fractions(4, 4):
        a = (k * k - 4 * k + 6) / (k * k + 18)
        y = 2 + k * (a - 1)
        for sign in (-1, 1):
            d = -2 * a + Fraction(4, 3) + sign * y / 3
            if not d:
                continue

            for rho in bounded_fractions(12, 6):
                denominator = 1 - rho * a * (1 - a) ** 2
                if not denominator:
                    continue

                # Coefficients of the linear-normalized polynomial
                # u+W+c2*W^2+c3*W^3+c4*W^4.
                c2 = (
                    -rho
                    * d
                    * (1 - 4 * a + 3 * a * a)
                    / (2 * denominator)
                )
                c3 = -rho * d * d * (a - Fraction(2, 3)) / denominator
                c4 = -rho * d**3 / (4 * denominator)

                # Weighted boundary cleanliness and admissibility.
                if not c2 or not c3 or not c4 or c2 == c4:
                    continue
                if c2 + c3 + c4:
                    raise AssertionError("tangent-chord parameterization failed")
                weighted_c = c2 - c4
                if (2 * c2 + 6 * c3 + 12 * c4) / weighted_c == -2:
                    continue

                seed = (c2, c3, c4)
                if seed in seen_seeds:
                    continue
                seen_seeds.add(seed)

                for constant in constants:
                    primitive = primitive_coefficients(
                        (c4, c3, c2, Fraction(1), constant)
                    )
                    if not irreducible_mod_17(primitive):
                        continue

                    # Recover the selected cancellation target and affine
                    # generator T_old=source_a+source_d*W.
                    source_a = rho * a
                    source_d = rho * d
                    pi = 1 / rho
                    integral_at_a = rho * (
                        a
                        - rho
                        * (
                            a * a / 2
                            - 2 * a**3 / 3
                            + a**4 / 4
                        )
                    )
                    linear_coefficient = (
                        rho * d * (1 - rho * a * (1 - a) ** 2)
                    )
                    target_r = integral_at_a - linear_coefficient * constant
                    metadata = (
                        k,
                        rho,
                        a,
                        d,
                        constant,
                        pi,
                        source_a,
                        source_d,
                        target_r,
                    )
                    score = (
                        max(abs(value) for value in primitive),
                        sum(abs(value) for value in primitive),
                    )
                    candidates.append((score, primitive, metadata))

    return sorted(candidates)


if __name__ == "__main__":
    results = search()
    best = results[0]
    expected = (
        (19, 50),
        (9, -19, 10, -8, -4),
        (
            Fraction(2),
            Fraction(11, 4),
            Fraction(1, 11),
            Fraction(12, 11),
            Fraction(1, 2),
            Fraction(4, 11),
            Fraction(1, 4),
            Fraction(3),
            -Fraction(22481, 23232),
        ),
    )
    assert best == expected
    assert results[1][0] > best[0]
    print(f"PASS: {len(results)} certified collisions in the bounded box")
    print("best primitive polynomial coefficients:", best[1])
    print(
        "k, rho, a, d, constant, Pi, T-origin, T-step, R:",
        best[2],
    )
    print("NOTE: this is bounded presentation minimality, not global minimality")
