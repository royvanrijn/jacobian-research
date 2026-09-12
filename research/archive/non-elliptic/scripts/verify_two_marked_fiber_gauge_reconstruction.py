#!/usr/bin/env python3
"""Exact checks for reconstruction from two marked power-shift fibers."""

from __future__ import annotations

from math import lcm

import sympy as sp


P, S, B, C = sp.symbols("P S B C")
u, v, a, b = sp.symbols("u v a b", nonzero=True)
g = sp.symbols("g0:13", nonzero=True)


def inverse_equation(degree: int, shift: int, base: sp.Expr) -> sp.Expr:
    """Power-shift inverse equation specialized to ``P=base``."""

    seed = g[1] * S + base * (g[2] * S**2 + g[3] * S**3)
    seed += sum(
        g[index] * base ** (index + shift) * S**index
        for index in range(4, degree + 1)
    )
    return sp.expand(seed - g[1] * (B * S**2 + C) / 2)


def monic_linear_coefficient(
    degree: int,
    shift: int,
    base: sp.Expr,
) -> sp.Expr:
    equation = inverse_equation(degree, shift, base)
    leading = sp.expand(equation).coeff(S, degree)
    return sp.cancel(sp.expand(equation / leading).coeff(S, 1))


# A marked fiber remembers the monic annihilator of its distinguished root.
# Its linear coefficient is independent of the target values B and C.
for degree in range(4, 12):
    for shift in range(8):
        linear_at_one = monic_linear_coefficient(degree, shift, sp.Integer(1))
        linear_at_c = monic_linear_coefficient(degree, shift, a)
        assert linear_at_one == g[1] / g[degree]
        assert sp.cancel(
            linear_at_c - g[1] / (g[degree] * a ** (degree + shift))
        ) == 0
        assert sp.cancel(linear_at_one / linear_at_c) == a ** (degree + shift)


# The P=1 marked polynomial recovers the normalized seed coefficients.
degree = 8
shift = 5
monic_at_one = sp.expand(
    inverse_equation(degree, shift, 1) / g[degree]
)
linear_at_one = monic_at_one.coeff(S, 1)
recovered = {
    1: linear_at_one,
    2: sp.cancel(monic_at_one.coeff(S, 2) + B * linear_at_one / 2),
}
recovered.update(
    {
        index: monic_at_one.coeff(S, index)
        for index in range(3, degree + 1)
    }
)
for index in range(1, degree + 1):
    assert sp.cancel(recovered[index] - g[index] / g[degree]) == 0


# One marked fiber away from P=1 does not recover both the seed and shift:
# changing m to n and rescaling every higher seed coefficient by c^(m-n)
# leaves the specialized equation literally unchanged.
for degree in range(4, 9):
    for first, second in ((0, 1), (2, 5), (6, 3)):
        original = inverse_equation(degree, first, a)
        alternative = (
            g[1] * S
            + a * (g[2] * S**2 + g[3] * S**3)
            + sum(
                g[index]
                * a ** (first - second)
                * a ** (index + second)
                * S**index
                for index in range(4, degree + 1)
            )
            - g[1] * (B * S**2 + C) / 2
        )
        assert sp.expand(original - alternative) == 0


# The universally available choice c=2 separates every pair of shifts.
for degree in range(4, 13):
    for first in range(10):
        for second in range(10):
            ratio_first = sp.Integer(2) ** (degree + first)
            ratio_second = sp.Integer(2) ** (degree + second)
            assert (ratio_first == ratio_second) == (first == second)
            recovered_shift = sp.factorint(
                int(ratio_first / (sp.Integer(2) ** degree))
            ).get(2, 0)
            assert recovered_shift == first
            recovered_area = 2 * degree - 3 + (degree - 2) * recovered_shift
            assert recovered_area == 2 * degree - 3 + (degree - 2) * first


# Torsion samples are sharp counterexamples.  On any finite collection of
# roots-of-unity planes, shifts separated by the lcm of their orders agree.
for orders in ((1,), (2,), (2, 3), (3, 4, 5), (4, 6, 10)):
    period = lcm(*orders)
    for order in orders:
        cyclotomic = sp.Poly(P**order - 1, P)
        for shift in range(5):
            first = inverse_equation(7, shift, P)
            second = inverse_equation(7, shift + period, P)
            assert sp.rem(sp.Poly(second - first, P), cyclotomic).is_zero
            assert (
                2 * 7 - 3 + (7 - 2) * shift
                != 2 * 7 - 3 + (7 - 2) * (shift + period)
            )


# On a transverse affine line P=a+b*u, the marked monic polynomial has
# linear coefficient g_1/(g_N P^(N+m)).  Its pole order at P=0 is N+m.
line_base = a + b * u
for degree in range(4, 10):
    for shift in range(7):
        line_linear = sp.factor(
            monic_linear_coefficient(degree, shift, line_base)
        )
        numerator, denominator = sp.fraction(line_linear)
        assert sp.expand(numerator - g[1]) == 0
        assert sp.expand(
            denominator - g[degree] * line_base ** (degree + shift)
        ) == 0
        centered_denominator = sp.expand(
            denominator.subs(u, v - a / b)
        )
        assert sp.cancel(
            centered_denominator
            / (g[degree] * b ** (degree + shift) * v ** (degree + shift))
        ) == 1


# Representative squarefree fibers show that the reconstruction data occur
# on a nonempty finite-etale locus at both P=1 and P=2.
quartic = S + S**3 + S**4
for shift_value in range(6):
    for base_value in (1, 2):
        equation = sp.expand(
            quartic
            + (base_value - 1) * S**3
            + (base_value ** (4 + shift_value) - 1) * S**4
        )
        assert sp.discriminant(equation, S) != 0


print("PASS: P=1 recovers the normalized seed from the marked annihilator")
print("PASS: the planes P=1 and P=2 recover the gauge exponent exactly")
print("PASS: torsion plane samples retain infinite stable ambiguity")
print("PASS: a transverse marked line sees m as the pole order N+m")
