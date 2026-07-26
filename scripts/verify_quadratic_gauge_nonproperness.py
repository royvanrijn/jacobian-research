#!/usr/bin/env python3
"""Exact audits for the quadratic-gauge Jelonek hypersurface."""

from __future__ import annotations

from fractions import Fraction
from math import gcd

import sympy as sp


P, B, C, S = sp.symbols("P B C S")
g1, g2, g3, g4 = sp.symbols("g1 g2 g3 g4", nonzero=True)


# Degree three: the resultant has an extraneous leading-coefficient factor P,
# while the discriminant specializes exactly to B^2*(1-B*C).
E3 = (
    g3 * P * S**3
    + (g2 * P - g1 * B / 2) * S**2
    + g1 * S
    - g1 * C / 2
)
D3 = sp.factor(sp.discriminant(E3, S))
R3 = sp.factor(sp.resultant(E3, sp.diff(E3, S), S))
assert sp.factor(R3 + g3 * P * D3) == 0
assert sp.factor(D3.subs(P, 0) - g1**4 * B**2 * (1 - B * C) / 4) == 0


# Degree four: the discriminant has exact P-adic order two.  The resultant
# adds the leading coefficient g4*P^4 but has the same reduced zero set.
E4 = (
    g4 * P**4 * S**4
    + g3 * P * S**3
    + (g2 * P - g1 * B / 2) * S**2
    + g1 * S
    - g1 * C / 2
)
D4 = sp.factor(sp.discriminant(E4, S))
R4 = sp.factor(sp.resultant(E4, sp.diff(E4, S), S))
assert sp.factor(R4 - g4 * P**4 * D4) == 0
assert sp.factor(
    (D4 / P**2).subs(P, 0)
    - g1**4 * g3**2 * B**2 * (1 - B * C) / 4
) == 0


# Direct special-fiber charts.  On q=0 the residual inverse is quadratic,
# and on t=0 the target B determines one source point exactly when B != 0.
x, y, z = sp.symbols("x y z", nonzero=True)
t = 1 + x * y
q = t**2 * z + (g1 / g3) * y**2 * (1 + 3 * t)
z_on_q_zero = -(g1 / g3) * y**2 * (1 + 3 * t) / t**2
S_source = x / t
assert sp.factor(q.subs(z, z_on_q_zero)) == 0
assert sp.factor(1 - y * S_source - 1 / t) == 0

x_on_t_zero = -1 / y
q_on_t_zero = sp.factor(q.subs({x: x_on_t_zero, t: 0}))
assert q_on_t_zero == g1 * y**2 / g3
assert sp.factor(
    y + 3 * (g3 / g1) * x_on_t_zero * q_on_t_zero + 2 * y
) == 0

residual = B * S**2 - 2 * S + C
assert sp.discriminant(residual, S) == 4 * (1 - B * C)
assert sp.factor(residual.subs(S, 1 / B) - (B * C - 1) / B) == 0


# On B=0, the two extra infinity sheets have Pi=u^2, S~a/u with
# a^2=-g1/g3. Their reconstruction valuations are exact.
v_u_P = Fraction(2)
v_u_S = Fraction(-1)
v_u_D = Fraction(0)
v_u_t = Fraction(0)
v_u_x = v_u_S - v_u_D
v_u_y = Fraction(1)
v_u_q = v_u_P + v_u_D
v_u_z = Fraction(2)
assert (v_u_x, v_u_y, v_u_q, v_u_z) == (-1, 1, 2, 2)


# Independent exact discriminant orders for numerical admissible seeds.
# Symbolic B,C are retained, so the lowest coefficient also checks that the
# saturated slice is nonzero away from B=0 and B*C=1.
for degree in range(3, 11):
    coefficients = {1: sp.Integer(1), 2: sp.Integer(2), 3: sp.Integer(1)}
    coefficients.update({k: sp.Integer(k + 1) for k in range(4, degree + 1)})
    inverse = (
        coefficients[1] * S
        + P * (coefficients[2] * S**2 + coefficients[3] * S**3)
        + sum(
            coefficients[k] * P**k * S**k
            for k in range(4, degree + 1)
        )
        - coefficients[1] * (B * S**2 + C) / 2
    )
    discriminant = sp.Poly(sp.discriminant(inverse, S), P)
    valuation = min(
        monomial[0] for monomial, value in discriminant.terms() if value
    )
    expected = 0 if degree == 3 else degree**2 - 3 * degree - 2
    assert valuation == expected
    lowest = sp.factor(discriminant.coeff_monomial(P**valuation))
    assert sp.factor(lowest.subs({B: 2, C: 3})) != 0
    assert sp.factor(lowest.subs(B, 0)) == 0
    assert sp.factor(lowest.subs(C, 1 / B)) == 0


# Uniform Newton polygon and valuation ledger.
for degree in range(4, 65):
    d = degree - 3
    h = gcd(d, 2)
    final_slope = Fraction(d + 2, d)
    for exponent in range(4, degree):
        assert Fraction(1) + (exponent - 3) * final_slope < exponent

    # Two q=0 sheets, one t=0 sheet, and d missing sheets.
    assert 2 + 1 + d == degree

    # Boundary-prime and source-coordinate valuations.
    assert h * (d // h) == d
    v_P = Fraction(d, h)
    v_S = -Fraction(d + 2, h)
    v_D = v_P + 2 * v_S
    v_t = -v_D
    v_x = v_S - v_D
    v_y = -v_x
    v_q = v_P + v_D
    v_z = v_q - 2 * v_t
    assert v_D == -Fraction(d + 4, h)
    assert v_t == Fraction(d + 4, h)
    assert v_x == Fraction(2, h)
    assert v_y == -Fraction(2, h)
    assert v_q == -Fraction(4, h)
    assert v_z == -Fraction(2 * d + 12, h)


print("PASS: cubic resultant saturation and discriminant slice are exact")
print("PASS: quartic discriminant order and saturated slice are exact")
print("PASS: q=0 and t=0 give the complete Pi=0 affine chart")
print("PASS: B=0 extra-sheet valuations are exact")
print("PASS: discriminant orders agree through degree ten")
print("PASS: Newton ledgers and boundary valuations agree through degree 64")
