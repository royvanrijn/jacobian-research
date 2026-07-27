#!/usr/bin/env python3
"""Exact five-row certificate for the fixed universal quintic calculator.

The bounded search is separate.  This checker uses only exact rational
arithmetic and finite-field factorization; it does not call PARI/GP's
Galois-group classifier.
"""

from __future__ import annotations

import math
import os
import warnings
from functools import reduce

os.environ.setdefault("SYMPY_GROUND_TYPES", "python")

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning

warnings.filterwarnings("ignore", category=SymPyDeprecationWarning)


x, y, z, S, T, X, Y = sp.symbols("x y z S T X Y")


def primitive_integer_polynomial(poly: sp.Poly) -> sp.Poly:
    """Clear denominators and content, with positive leading coefficient."""

    _, cleared = sp.Poly(poly, domain=sp.QQ).clear_denoms(convert=True)
    _, primitive = cleared.primitive()
    if primitive.LC() < 0:
        primitive = -primitive
    return sp.Poly(primitive, poly.gens[0], domain=sp.ZZ)


def factor_degrees_mod_prime(poly: sp.Poly, prime: int) -> tuple[int, ...]:
    """Return the squarefree factor-degree partition modulo ``prime``."""

    poly = primitive_integer_polynomial(poly)
    variable = poly.gens[0]
    assert int(poly.LC()) % prime
    assert int(sp.discriminant(poly.as_expr(), variable)) % prime
    factors = sp.factor_list(poly.as_expr(), modulus=prime)[1]
    return tuple(
        sorted(
            (
                int(sp.degree(factor, variable))
                for factor, exponent in factors
                for _ in range(exponent)
            ),
            reverse=True,
        )
    )


def projective_target(
    coordinates: tuple[int, int, int, int],
) -> tuple[sp.Rational, sp.Rational, sp.Rational]:
    """Convert primitive ``[W:P:B:C]`` to the affine target."""

    w, pi, b, c = coordinates
    assert w > 0 and pi
    assert reduce(math.gcd, (abs(value) for value in coordinates)) == 1
    return sp.Rational(pi, w), sp.Rational(b, w), sp.Rational(c, w)


def inverse_polynomial(
    target: tuple[sp.Rational, sp.Rational, sp.Rational],
) -> sp.Poly:
    pi, b, c = target
    return sp.Poly(
        pi**5 * S**5 - 5 * pi * S**3 - 2 * b * S**2 + 4 * S - 2 * c,
        S,
        domain=sp.QQ,
    )


def normalized_polynomial(
    target: tuple[sp.Rational, sp.Rational, sp.Rational],
) -> sp.Poly:
    """Return ``Pi^5 E(Pi^-2 T)``."""

    pi, b, c = target
    return sp.Poly(
        T**5
        - 5 * T**3
        - 2 * pi * b * T**2
        + 4 * pi**3 * T
        - 2 * pi**5 * c,
        T,
        domain=sp.QQ,
    )


def pair_sum_resolvent(a: sp.Rational, b: sp.Rational, c: sp.Rational) -> sp.Poly:
    """Pair-sum resolvent of ``T^5-5T^3+aT^2+bT+c``."""

    return sp.Poly(
        X**10
        - 15 * X**8
        + a * X**7
        + (75 - 3 * b) * X**6
        + (-10 * a - 11 * c) * X**5
        + (-a**2 + 10 * b - 125) * X**4
        + (-4 * a * b + 25 * a + 20 * c) * X**3
        + (5 * a**2 + 7 * a * c - 4 * b**2 + 25 * b) * X**2
        + (-a**3 + 4 * b * c - 25 * c) * X
        - a**2 * b
        - 5 * a * c
        - c**2,
        X,
        domain=sp.QQ,
    )


def dummit_resolvent(
    p: sp.Rational,
    q: sp.Rational,
    r: sp.Rational,
    s: sp.Rational,
) -> sp.Poly:
    """Dummit's F_20 sextic for ``T^5+pT^3+qT^2+rT+s``.

    Formula (2) of D. S. Dummit, "Solving Solvable Quintics",
    Math. Comp. 57 (1991), 387--401.  The author's exact Mathematica
    notebook is the transcription source.
    """

    c4 = 2 * p * q**2 - 6 * p**2 * r + 40 * r**2 - 50 * q * s
    c3 = -(
        2 * q**4
        - 21 * p * q**2 * r
        + 40 * p**2 * r**2
        - 160 * r**3
        + 15 * p**2 * q * s
        + 400 * q * r * s
        - 125 * p * s**2
    )
    c2 = (
        p**2 * q**4
        - 6 * p**3 * q**2 * r
        - 8 * q**4 * r
        + 9 * p**4 * r**2
        + 76 * p * q**2 * r**2
        - 136 * p**2 * r**3
        + 400 * r**4
        - 50 * p * q**3 * s
        + 90 * p**2 * q * r * s
        - 1400 * q * r**2 * s
        + 625 * q**2 * s**2
        + 500 * p * r * s**2
    )
    c1 = -(
        2 * p * q**6
        - 19 * p**2 * q**4 * r
        + 51 * p**3 * q**2 * r**2
        - 3 * q**4 * r**2
        - 32 * p**4 * r**3
        - 76 * p * q**2 * r**3
        + 256 * p**2 * r**4
        - 512 * r**5
        + 31 * p**3 * q**3 * s
        + 58 * q**5 * s
        - 117 * p**4 * q * r * s
        - 105 * p * q**3 * r * s
        - 260 * p**2 * q * r**2 * s
        + 2400 * q * r**3 * s
        + 108 * p**5 * s**2
        + 325 * p**2 * q**2 * s**2
        - 525 * p**3 * r * s**2
        - 2750 * q**2 * r * s**2
        + 500 * p * r**2 * s**2
        - 625 * p * q * s**3
        + 3125 * s**4
    )
    c0 = (
        q**8
        - 13 * p * q**6 * r
        + p**5 * q**2 * r**2
        + 65 * p**2 * q**4 * r**2
        - 4 * p**6 * r**3
        - 128 * p**3 * q**2 * r**3
        + 17 * q**4 * r**3
        + 48 * p**4 * r**4
        - 16 * p * q**2 * r**4
        - 192 * p**2 * r**5
        + 256 * r**6
        - 4 * p**5 * q**3 * s
        - 12 * p**2 * q**5 * s
        + 18 * p**6 * q * r * s
        + 12 * p**3 * q**3 * r * s
        - 124 * q**5 * r * s
        + 196 * p**4 * q * r**2 * s
        + 590 * p * q**3 * r**2 * s
        - 160 * p**2 * q * r**3 * s
        - 1600 * q * r**4 * s
        - 27 * p**7 * s**2
        - 150 * p**4 * q**2 * s**2
        - 125 * p * q**4 * s**2
        - 99 * p**5 * r * s**2
        - 725 * p**2 * q**2 * r * s**2
        + 1200 * p**3 * r**2 * s**2
        + 3250 * q**2 * r**2 * s**2
        - 2000 * p * r**3 * s**2
        - 1250 * p * q * r * s**3
        + 3125 * p**2 * s**4
        - 9375 * r * s**4
    )
    return sp.Poly(
        X**6 + 8 * r * X**5 + c4 * X**4 + c3 * X**3 + c2 * X**2 + c1 * X + c0,
        X,
        domain=sp.QQ,
    )


# The fixed determinant-minus-two Keller map.
t = 1 + x * y
q_source = t**2 * z - sp.Rational(4, 5) * y**2 * (1 + 3 * t)
mapping = (
    t * q_source,
    y
    - sp.Rational(15, 4) * x * q_source
    + sp.Rational(5, 4) * t**2 * x**3 * q_source**5,
    x * (5 - 3 * t)
    + sp.Rational(5, 4) * x**3 * z
    - sp.Rational(3, 4) * (x * q_source) ** 5,
)
assert sp.factor(sp.Matrix(mapping).jacobian((x, y, z)).det()) == -2


coordinates = {
    "S5": (1, -1, -1, -1),
    "A5": (5, 5, 0, -2),
    "C5": (10, 10, 0, -7),
    "D5": (10, 4, -21, 20),
    "F20": (10, 5, 15, 4),
}
targets = {name: projective_target(value) for name, value in coordinates.items()}
inverse = {name: inverse_polynomial(value) for name, value in targets.items()}
normalized = {name: normalized_polynomial(value) for name, value in targets.items()}
primitive_inverse = {
    name: primitive_integer_polynomial(value) for name, value in inverse.items()
}
primitive_normalized = {
    name: primitive_integer_polynomial(value) for name, value in normalized.items()
}

assert max(max(abs(value) for value in row) for row in coordinates.values()) == 21
assert {name: max(abs(value) for value in row) for name, row in coordinates.items()} == {
    "S5": 1,
    "A5": 5,
    "C5": 10,
    "D5": 21,
    "F20": 15,
}
assert {
    name: max(abs(int(value)) for value in poly.all_coeffs())
    for name, poly in primitive_inverse.items()
} == {
    "S5": 5,
    "A5": 25,
    "C5": 25,
    "D5": 13125,
    "F20": 640,
}


# S_5: transitivity, a transposition, and odd discriminant square class.
assert normalized["S5"].as_expr() == T**5 - 5 * T**3 - 2 * T**2 - 4 * T - 2
assert factor_degrees_mod_prime(normalized["S5"], 5) == (5,)
assert factor_degrees_mod_prime(normalized["S5"], 43) == (2, 1, 1, 1)
assert sp.discriminant(normalized["S5"].as_expr(), T) == -16 * 3 * 61813


# A_5: transitivity, containment in A_5, and a 3-cycle.
assert primitive_normalized["A5"].as_expr() == 5 * T**5 - 25 * T**3 + 20 * T + 4
assert factor_degrees_mod_prime(normalized["A5"], 3) == (5,)
assert sp.discriminant(normalized["A5"].as_expr(), T) == 232**2
assert factor_degrees_mod_prime(normalized["A5"], 23) == (3, 1, 1)


# C_5: irreducibility and a nontrivial explicit order-five automorphism.
assert primitive_normalized["C5"].as_expr() == 5 * T**5 - 25 * T**3 + 20 * T + 7
assert factor_degrees_mod_prime(normalized["C5"], 2) == (5,)
c5_automorphism = (
    45 * T**4 - 20 * T**3 - 240 * T**2 + 78 * T + 174
) / 43
c5_equation = normalized["C5"].as_expr()
assert sp.rem(
    c5_equation.subs(T, c5_automorphism),
    c5_equation,
    T,
) == 0
c5_iterate = T
for exponent in range(1, 6):
    c5_iterate = sp.cancel(
        sp.rem(
            sp.together(c5_iterate.subs(T, c5_automorphism)),
            c5_equation,
            T,
        )
    )
    if exponent < 5:
        assert sp.expand(c5_iterate - T) != 0
assert sp.expand(c5_iterate - T) == 0


# D_5: square discriminant, two pair orbits, and a reflection.
assert primitive_inverse["D5"].as_expr() == (
    32 * S**5 - 6250 * S**3 + 13125 * S**2 + 12500 * S - 12500
)
d5_discriminant = sp.discriminant(normalized["D5"].as_expr(), T)
assert d5_discriminant == sp.Rational(2040064, 78125) ** 2
d5_resolvent = pair_sum_resolvent(
    sp.Rational(42, 25),
    sp.Rational(32, 125),
    -sp.Rational(128, 3125),
)
d5_equation = normalized["D5"].as_expr()
assert sp.expand(
    sp.resultant(
        d5_equation.subs(T, Y),
        d5_equation.subs(T, X - Y),
        Y,
    )
    - 2**5
    * d5_equation.subs(T, X / 2)
    * d5_resolvent.as_expr() ** 2
) == 0
d5_factors = (
    sp.Poly(
        3125 * X**5 - 28750 * X**3 - 2000 * X**2 + 62725 * X + 4784,
        X,
        domain=sp.ZZ,
    ),
    sp.Poly(
        3125 * X**5 - 18125 * X**3 + 7250 * X**2 + 2500 * X - 776,
        X,
        domain=sp.ZZ,
    ),
)
assert sp.expand(
    d5_resolvent.as_expr()
    - d5_factors[0].as_expr() * d5_factors[1].as_expr() / 9765625
) == 0
assert all(factor_degrees_mod_prime(factor, 3) == (5,) for factor in d5_factors)
assert factor_degrees_mod_prime(normalized["D5"], 11) == (2, 2, 1)


# F_20: transitivity, Dummit's solvability sextic, and odd square class.
assert primitive_inverse["F20"].as_expr() == (
    5 * S**5 - 400 * S**3 - 480 * S**2 + 640 * S - 128
)
assert factor_degrees_mod_prime(normalized["F20"], 29) == (5,)
f20_discriminant = sp.discriminant(normalized["F20"].as_expr(), T)
assert f20_discriminant == 5 * sp.Rational(2127, 320) ** 2
f20_resolvent = dummit_resolvent(
    sp.Integer(-5),
    -sp.Rational(3, 2),
    sp.Rational(1, 2),
    -sp.Rational(1, 40),
)
f20_resolvent_integer = primitive_integer_polynomial(f20_resolvent)
f20_resolvent_cofactor = (
    20480 * X**5
    - 51200 * X**4
    - 1497600 * X**3
    + 1947840 * X**2
    + 24055920 * X
    + 24084531
)
assert sp.expand(
    f20_resolvent_integer.as_expr() - (2 * X + 13) * f20_resolvent_cofactor
) == 0
assert f20_resolvent.eval(-sp.Rational(13, 2)) == 0


print("PASS: the fixed quintic map has determinant -2")
print("PASS: all five transitive quintic groups occur below projective height 22")
print("PASS: the row heights are S5=1, A5=5, C5=10, F20=15, D5=21")
print("PASS: every row uses at most three exact integer or modular checks")
print("PASS: no numerical root approximation or Galois-group oracle is used")
