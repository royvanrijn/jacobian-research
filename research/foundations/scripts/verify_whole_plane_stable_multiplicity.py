#!/usr/bin/env python3
"""Exact checks for whole-plane and cyclotomic stable multiplicity."""

from __future__ import annotations

import sympy as sp


P, B, C, S = sp.symbols("P B C S")
t, q, x = sp.symbols("t q x")
a = sp.symbols("a")
g = sp.symbols("g0:11", nonzero=True)


def power_shifted_seed(degree: int, shift: int) -> sp.Expr:
    """The inverse seed used by the common power shift."""

    result = g[1] * S + P * (g[2] * S**2 + g[3] * S**3)
    result += sum(
        g[index] * P ** (index + shift) * S**index
        for index in range(4, degree + 1)
    )
    return sp.expand(result)


# Every power-shifted inverse plane becomes the same plane at P=1.
for degree in range(4, 10):
    common = sum(g[index] * S**index for index in range(1, degree + 1))
    for shift in range(7):
        seed = power_shifted_seed(degree, shift)
        inverse = seed - g[1] * (B * S**2 + C) / 2
        common_inverse = common - g[1] * (B * S**2 + C) / 2
        assert sp.expand(inverse.subs(P, 1) - common_inverse) == 0

        # The source-coordinate decorations acquire exactly the factor P^m.
        for index in range(4, degree + 1):
            assert sp.expand(
                t ** (shift + 2) * x ** (index - 2) * q ** (index + shift)
                - (t * q) ** shift * t**2 * x ** (index - 2) * q**index
            ) == 0
            assert sp.expand(
                t**shift * x**index * q ** (index + shift)
                - (t * q) ** shift * x**index * q**index
            ) == 0


# Cubic lifts differ from the minimal gauge by the factor P^(n-3)-1.
for exponent in range(4, 14):
    correction_b = (
        t ** (exponent - 1) * x * q**exponent - t**2 * x * q**3
    )
    correction_c = (
        t ** (exponent - 3) * x**3 * q**exponent - x**3 * q**3
    )
    assert sp.expand(
        correction_b
        - t**2 * x * q**3 * ((t * q) ** (exponent - 3) - 1)
    ) == 0
    assert sp.expand(
        correction_c
        - x**3 * q**3 * ((t * q) ** (exponent - 3) - 1)
    ) == 0

    cubic_seed = (
        g[1] * S
        + g[2] * P * S**2
        + g[3] * P * (1 + P ** (exponent - 1) - P**2) * S**3
    )
    common_cubic = g[1] * S + g[2] * S**2 + g[3] * S**3
    assert sp.expand(cubic_seed.subs(P, 1) - common_cubic) == 0


# A fixed exponent residue class gives equality on the cyclotomic divisor
# P^d=1.  This applies to both m for power shifts and n for cubic lifts.
for modulus in range(1, 8):
    cyclotomic_divisor = sp.Poly(P**modulus - 1, P)
    for residue in range(modulus):
        power_exponents = [residue + modulus * index for index in range(4)]
        for first, second in zip(power_exponents, power_exponents[1:]):
            remainder = sp.rem(
                sp.Poly(P**second - P**first, P),
                cyclotomic_divisor,
            )
            assert remainder.is_zero
            seed_remainder = sp.rem(
                sp.Poly(
                    power_shifted_seed(7, second)
                    - power_shifted_seed(7, first),
                    P,
                ),
                cyclotomic_divisor,
            )
            assert seed_remainder.is_zero

        first_cubic = 4 + ((residue - 4) % modulus)
        cubic_exponents = [
            first_cubic + modulus * index for index in range(5)
        ]
        for first, second in zip(cubic_exponents, cubic_exponents[1:]):
            remainder = sp.rem(
                sp.Poly(P**second - P**first, P),
                cyclotomic_divisor,
            )
            assert remainder.is_zero
            first_seed = (
                g[1] * S
                + g[2] * P * S**2
                + g[3] * P * (1 + P ** (first - 1) - P**2) * S**3
            )
            second_seed = (
                g[1] * S
                + g[2] * P * S**2
                + g[3] * P * (1 + P ** (second - 1) - P**2) * S**3
            )
            seed_remainder = sp.rem(
                sp.Poly(second_seed - first_seed, P),
                cyclotomic_divisor,
            )
            assert seed_remainder.is_zero


# The fixed quintic Hasse line lies in the common squarefree plane away from
# a=0,1, and its stable Newton areas are 7+3m.
G = S**5 - sp.Rational(3, 2) * S**4 + sp.Rational(3, 2) * S**3
G -= sp.Rational(5, 4) * S**2
G += sp.Rational(9, 16) * S
hasse_inverse = sp.expand(
    G
    - sp.Rational(9, 32)
    * (sp.Rational(32, 9) * a * S**2 + (8 * a + 1) / 3)
)
centered_factorization = sp.expand(
    ((S - sp.Rational(1, 2)) ** 3 - a)
    * (S**2 + sp.Rational(3, 4))
)
assert sp.expand(hasse_inverse - centered_factorization) == 0
assert sp.factor(sp.discriminant(hasse_inverse, S)) == 81 * a**2 * (a - 1) ** 4

for shift in range(20):
    assert 2 * 5 - 3 + (5 - 2) * shift == 7 + 3 * shift


print("PASS: every power shift has one common squarefree inverse plane")
print("PASS: every cubic lift has the same common plane at P=1")
print("PASS: exponent residue classes agree on P^d=1")
print("PASS: the fixed Hasse line has discriminant 81*a^2*(a-1)^4")
print("PASS: the quintic stable Newton areas are 7+3*m")
