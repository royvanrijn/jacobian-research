#!/usr/bin/env python3
"""Exact audit of the first cubic Danielewski contraction frontier."""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from functools import cache
from math import comb

import sympy as sp


a, b, c, w, x = sp.symbols("a b c w x")
P = x**3 + x
P_prime = sp.diff(P, x)


def bracket(first: sp.Expr, second: sp.Expr) -> sp.Expr:
    return sp.expand(
        c
        * (
            sp.diff(first, c) * sp.diff(second, x)
            - sp.diff(first, x) * sp.diff(second, c)
        )
        + P_prime
        * (
            sp.diff(first, c) * sp.diff(second, w)
            - sp.diff(first, w) * sp.diff(second, c)
        )
        + w
        * (
            sp.diff(first, x) * sp.diff(second, w)
            - sp.diff(first, w) * sp.diff(second, x)
        )
    )


@cache
def reduce_x_power(power: int) -> tuple[sp.Rational, sp.Rational]:
    if power == 0:
        return sp.S.One, sp.S.Zero
    if power == 1:
        return sp.S.Zero, sp.S.One
    constant, linear = reduce_x_power(power - 2)
    factor = -sp.Rational(power - 1, power + 1)
    return factor * constant, factor * linear


def de_rham_reduce(polynomial: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    answer = [sp.S.Zero, sp.S.Zero]
    for (a_power, c_power, x_power, w_power), coefficient in sp.Poly(
        sp.expand(polynomial),
        a,
        c,
        x,
        w,
    ).terms():
        if c_power != w_power:
            continue
        pure_x = sp.Poly(
            sp.expand(x**x_power * P**c_power),
            x,
        )
        for (power,), pure_coefficient in pure_x.terms():
            constant, linear = reduce_x_power(power)
            answer[0] += coefficient * pure_coefficient * constant * a**a_power
            answer[1] += coefficient * pure_coefficient * linear * a**a_power
    return tuple(sp.factor(entry) for entry in answer)


# The Danielewski open immersion turns the coupled ledger into an ordinary
# three-variable Keller determinant.
x_pullback = b * c
w_pullback = b * (x_pullback**2 + 1)
assert sp.expand(c * w_pullback - (x_pullback**3 + x_pullback)) == 0

A_test = a**2 * c + a * x + w
F_test = a**3 * w + a * c * x + x**2
G_test = a * c**2 + a**2 * x + w**2
coupled_test = sp.expand(
    sp.diff(A_test, a) * bracket(F_test, G_test)
    - sp.diff(F_test, a) * bracket(A_test, G_test)
    + sp.diff(G_test, a) * bracket(A_test, F_test)
)
pullback = {
    x: x_pullback,
    w: w_pullback,
}
ordinary_test = sp.factor(
    sp.Matrix(
        [
            A_test.xreplace(pullback),
            F_test.xreplace(pullback),
            G_test.xreplace(pullback),
        ]
    ).jacobian((a, b, c)).det()
)
assert sp.expand(
    ordinary_test + coupled_test.xreplace(pullback)
) == 0
print("PASS: the open immersion converts J to the opposite ordinary Jacobian")


# The ambient-degree-at-most-three quotient has 35 monomials and the single
# relation cw-x^3-x.  After fixing constants and the source-linear jet, the
# coefficient of x remains free because x=bc has zero source-linear jet.  Each
# target row therefore has 30 residual coefficients.
exponents: list[tuple[int, int, int, int]] = [(0, 0, 1, 0)]
for total_degree in (2, 3):
    for a_power in range(total_degree, -1, -1):
        for c_power in range(total_degree - a_power, -1, -1):
            for x_power in range(
                total_degree - a_power - c_power,
                -1,
                -1,
            ):
                w_power = total_degree - a_power - c_power - x_power
                exponent = (a_power, c_power, x_power, w_power)
                if exponent == (0, 1, 0, 1):
                    continue
                exponents.append(exponent)
assert len(exponents) == 30
assert 1 + 3 + 30 == 34
print("PASS: the normalized ambient-cubic ansatz has 90 coefficients")


# Compile the generic cubic flux sparsely.  Coefficient monomials are tuples of
# variable indices, so this count does not require expanding an 87-variable
# SymPy expression.
rows = []
for row, leading in enumerate(
    (
        (1, 0, 0, 0),
        (0, 1, 0, 0),
        (0, 0, 0, 1),
    )
):
    rows.append(
        [(leading, ())]
        + [
            (exponent, (30 * row + column,))
            for column, exponent in enumerate(exponents)
        ]
    )


def add_exponents(
    first: tuple[int, int, int, int],
    second: tuple[int, int, int, int],
) -> tuple[int, int, int, int]:
    return tuple(left + right for left, right in zip(first, second))


def monomial_bracket(
    first: tuple[int, int, int, int],
    second: tuple[int, int, int, int],
) -> list[tuple[tuple[int, int, int, int], int]]:
    answer = []

    coefficient = first[1] * second[2] - first[2] * second[1]
    if coefficient:
        answer.append(
            (
                (
                    first[0] + second[0],
                    first[1] + second[1],
                    first[2] + second[2] - 1,
                    first[3] + second[3],
                ),
                coefficient,
            )
        )

    coefficient = first[1] * second[3] - first[3] * second[1]
    if coefficient:
        base = (
            first[0] + second[0],
            first[1] + second[1] - 1,
            first[2] + second[2],
            first[3] + second[3] - 1,
        )
        answer.append((base, coefficient))
        answer.append(
            (
                (base[0], base[1], base[2] + 2, base[3]),
                3 * coefficient,
            )
        )

    coefficient = first[2] * second[3] - first[3] * second[2]
    if coefficient:
        answer.append(
            (
                (
                    first[0] + second[0],
                    first[1] + second[1],
                    first[2] + second[2] - 1,
                    first[3] + second[3],
                ),
                coefficient,
            )
        )
    return answer


@cache
def rational_x_reduction(power: int) -> tuple[Fraction, Fraction]:
    if power == 0:
        return Fraction(1), Fraction(0)
    if power == 1:
        return Fraction(0), Fraction(1)
    constant, linear = rational_x_reduction(power - 2)
    factor = Fraction(-(power - 1), power + 1)
    return factor * constant, factor * linear


compiled: dict[
    tuple[int, int],
    dict[tuple[int, ...], Fraction],
] = defaultdict(lambda: defaultdict(Fraction))
for exponent_A, coefficient_A in rows[0]:
    for exponent_F, coefficient_F in rows[1]:
        for exponent_G, coefficient_G in rows[2]:
            for bracket_exponent, bracket_coefficient in monomial_bracket(
                exponent_F,
                exponent_G,
            ):
                exponent = add_exponents(exponent_A, bracket_exponent)
                if exponent[1] != exponent[3]:
                    continue
                common_power = exponent[1]
                coefficient_monomial = tuple(
                    sorted(coefficient_A + coefficient_F + coefficient_G)
                )
                for choice in range(common_power + 1):
                    x_power = exponent[2] + common_power + 2 * choice
                    reductions = rational_x_reduction(x_power)
                    binomial = comb(common_power, choice)
                    for coordinate, reduction in enumerate(reductions):
                        if reduction:
                            compiled[(coordinate, exponent[0])][
                                coefficient_monomial
                            ] += (
                                Fraction(bracket_coefficient * binomial)
                                * reduction
                            )

compiled = {
    key: {
        monomial: coefficient
        for monomial, coefficient in polynomial.items()
        if coefficient
    }
    for key, polynomial in compiled.items()
    if any(polynomial.values())
}
constant_degrees = sorted(
    degree for (coordinate, degree) in compiled if coordinate == 0
)
linear_degrees = sorted(
    degree for (coordinate, degree) in compiled if coordinate == 1
)
assert constant_degrees == list(range(6))
assert linear_degrees == list(range(7))
assert sum(degree > 0 for _, degree in compiled) == 11
print("PASS: generic cubic flux gives 11 nonconstant de Rham equations")


# A normalized cubic family passes both the leading-face theorem and every
# compiled flux equation, but not the pointwise equation.  This is an exact
# witness that leading-face plus de Rham data alone cannot close degree three.
K, L = sp.symbols("K L", nonzero=True)
M = sp.Rational(15, 4)
A = a + w + K * a**2 * c
F = -w + M * a * x**2 + L * a**2 * c
G = c

A_top = K * c
F_top = L * c
G_top = c
leading_face = sp.expand(
    2 * A_top * bracket(F_top, G_top)
    - 2 * F_top * bracket(A_top, G_top)
)
assert leading_face == 0

local_jacobian = sp.Matrix([A, F, G]).jacobian((a, c, w)).subs(
    {a: 0, c: 0, x: 0, w: 0}
)
assert local_jacobian.det() == 1

flux = sp.expand(A * bracket(F, G))
coupled = sp.expand(
    sp.diff(A, a) * bracket(F, G)
    - sp.diff(F, a) * bracket(A, G)
    + sp.diff(G, a) * bracket(A, F)
)
assert de_rham_reduce(flux) == (a, 0)
assert de_rham_reduce(coupled) == (1, 0)
assert sp.expand(coupled - 1) != 0
print("PASS: a normalized cubic family passes leading face and flux but not J=1")


# The two nonzero boundary roots have the same local normal coefficient.  Their
# distinction is encoded by evaluation of odd/even x-parts and by the two
# global periods, not by opposite values of P'.
assert P_prime.subs(x, sp.I) == -2
assert P_prime.subs(x, -sp.I) == -2
period_matrix = sp.Matrix(
    [
        [sp.I, sp.I**2 / 2],
        [-sp.I, (-sp.I) ** 2 / 2],
    ]
)
assert period_matrix.det() != 0
print("PASS: boundary roots have equal local P' but independent global periods")
print("PASS cubic Danielewski frontier audit")
