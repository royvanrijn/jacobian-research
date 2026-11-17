#!/usr/bin/env python3
"""Verify the exact moment identities for balanced cubic GVC in two variables.

For a binary cubic operator symbol A(u,v) and

    P(x,y) = a*x^3 + b*x^2*y + c*x*y^2 + d*y^3,

the scalar contraction is the apolar pairing of A^m and P^m.  The three
nonzero GL_2-orbits of binary cubics are represented by

    u^3, u^2*v, u*v*(u+v).

The checker derives the first required moments directly and verifies the
branch factorizations used in the proof.  The all-order conclusion is the
degree/annihilating-direction argument in the canonical note.
"""

from __future__ import annotations

from math import factorial

import sympy as sp


x, y, u, v = sp.symbols("x y u v")
a, b, c, d = sp.symbols("a b c d")

P = a * x**3 + b * x**2 * y + c * x * y**2 + d * y**3


def assert_zero(expression: sp.Expr) -> None:
    """Assert an exact polynomial identity."""

    assert sp.expand(expression) == 0


def apolar_moment(symbol: sp.Expr, order: int) -> sp.Expr:
    """Return A(partial)^order(P^order) by coefficientwise apolarity."""

    symbol_power = sp.Poly(sp.expand(symbol**order), u, v)
    polynomial_power = sp.Poly(sp.expand(P**order), x, y)
    value = 0
    for (x_order, y_order), coefficient in symbol_power.terms():
        value += (
            coefficient
            * polynomial_power.coeff_monomial(x**x_order * y**y_order)
            * factorial(x_order)
            * factorial(y_order)
        )
    return sp.expand(value)


# Triple-root orbit.
triple = u**3
assert apolar_moment(triple, 1) == 6 * a
assert apolar_moment(triple, 2) == 720 * a**2


# Double-root orbit.
double = u**2 * v
double_moments = [apolar_moment(double, order) for order in range(1, 4)]
assert double_moments == [
    2 * b,
    48 * (2 * a * c + b**2),
    4320 * (3 * a**2 * d + 6 * a * b * c + b**3),
]
assert sp.factor(double_moments[1].subs(b, 0)) == 96 * a * c
assert sp.factor(double_moments[2].subs(b, 0)) == 12960 * a**2 * d


# Squarefree orbit.
squarefree = u * v * (u + v)
squarefree_moments = [
    apolar_moment(squarefree, order) for order in range(1, 5)
]
assert squarefree_moments[0] == 2 * (b + c)

restricted = [sp.expand(moment.subs(c, -b)) for moment in squarefree_moments]
r2 = -2 * a * b + 3 * a * d - b**2 + 2 * b * d
r4 = (
    14 * a**2 * b**2
    - 42 * a**2 * b * d
    + 45 * a**2 * d**2
    + 12 * a * b**3
    - 40 * a * b**2 * d
    + 42 * a * b * d**2
    + 3 * b**4
    - 12 * b**3 * d
    + 14 * b**2 * d**2
)
assert restricted[1] == 48 * r2
assert sp.expand(restricted[2] - 12960 * a * d * (a + d)) == 0
assert restricted[3] == 414720 * r4

# The third moment gives a=0, d=0, or d=-a.  On each branch the second
# moment has two factors, and the fourth rejects exactly the spurious one.
assert_zero(r2.subs(a, 0) - b * (2 * d - b))
assert_zero(
    r4.subs(a, 0) - b**2 * (3 * b**2 - 12 * b * d + 14 * d**2)
)
assert_zero(r4.subs({a: 0, b: 2 * d}) - 8 * d**4)

assert_zero(r2.subs(d, 0) + b * (2 * a + b))
assert_zero(
    r4.subs(d, 0) - b**2 * (14 * a**2 + 12 * a * b + 3 * b**2)
)
assert_zero(r4.subs({d: 0, b: -2 * a}) - 8 * a**4)

assert_zero(r2.subs(d, -a) + (a + b) * (3 * a + b))
assert_zero(
    r4.subs(d, -a)
    - (3 * a + b) ** 2 * (5 * a**2 + 6 * a * b + 3 * b**2)
)
assert_zero(r4.subs({d: -a, b: -a}) - 8 * a**4)
assert_zero(r4.subs({d: -a, b: -3 * a}))


# The surviving squarefree forms are cubes annihilated by one of the three
# commuting linear factors of the operator.
survivors = (x**3, y**3, (x - y) ** 3)
annihilators = (
    sp.diff(survivors[0], y),
    sp.diff(survivors[1], x),
    sp.diff(survivors[2], x) + sp.diff(survivors[2], y),
)
assert annihilators == (0, 0, 0)

print("PASS: the three binary-cubic symbol orbits have the claimed moments")
print("PASS: moments 1 through 4 leave only one-sided cubic Segre forms")
print("PASS: every squarefree survivor has a constant annihilating direction")
