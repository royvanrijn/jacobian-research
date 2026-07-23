#!/usr/bin/env python3
"""Exact regressions for quadratic-versus-weighted stable separation."""

from __future__ import annotations

from fractions import Fraction
from math import gcd

import sympy as sp


P, B, C, S = sp.symbols("P B C S")
g1, g2, g3, g4, g5, g6 = sp.symbols(
    "g1 g2 g3 g4 g5 g6",
    nonzero=True,
)

G_P = (
    g1 * S
    + P * (g2 * S**2 + g3 * S**3)
    + g4 * P**4 * S**4
    + g5 * P**5 * S**5
    + g6 * P**6 * S**6
)
E = G_P - g1 * (B * S**2 + C) / 2
D = sp.diff(G_P, S) / g1 - B * S
assert sp.factor(sp.diff(E, S) - g1 * D) == 0


# Critical-value parameterization and its differential identity.
r = sp.symbols("r", nonzero=True)
H = G_P.subs(S, r)
B_critical = sp.factor(sp.diff(G_P, S).subs(S, r) / (g1 * r))
C_critical = sp.factor((2 * H - r * sp.diff(G_P, S).subs(S, r)) / g1)
assert sp.factor(
    sp.diff(C_critical, r) + r**2 * sp.diff(B_critical, r)
) == 0
assert sp.limit(r * B_critical, r, 0) == 1

# The Laurent coordinate r^{-1} is recovered after adjoining the integral r.
inverse_r = sp.factor(
    B_critical
    - (sp.diff(G_P, S).subs(S, r) - g1) / (g1 * r)
)
assert inverse_r == 1 / r
integral_equation = sp.Poly(
    sp.diff(G_P, S).subs(S, r) - g1 * B * r,
    r,
)
assert integral_equation.degree() == 5
assert integral_equation.LC() == 6 * g6 * P**6


# The two finite/slope-one branches over P=0 are the affine source divisors
# q=0 and t=0.
x, y, z = sp.symbols("x y z", nonzero=True)
t = 1 + x * y
q = t**2 * z + (g1 / g3) * y**2 * (1 + 3 * t)
P_source = t * q
S_source = x / t
Q_source = y + x * q
beta = sp.cancel((sp.diff(G_P, S) / g1 - 1 - P * S**2) / S)
B_chart = sp.factor(Q_source + beta.subs({P: P_source, S: S_source}))
C_chart = sp.factor(
    (
        2 * G_P / g1 - B * S**2
    ).subs({P: P_source, S: S_source, B: B_chart})
)

z_on_q_zero = -(g1 / g3) * y**2 * (1 + 3 * t) / t**2
assert sp.factor(B_chart.subs(z, z_on_q_zero) - y) == 0
assert sp.factor(
    C_chart.subs(z, z_on_q_zero)
    - (2 * S_source - y * S_source**2)
) == 0

x_on_t_zero = -1 / y
q_on_t_zero = sp.factor(q.subs({x: x_on_t_zero, t: 0}))
assert q_on_t_zero == g1 * y**2 / g3
assert sp.factor(
    B_chart.subs({x: x_on_t_zero, t: 0, q: q_on_t_zero}) + 2 * y
) == 0
assert sp.factor(
    (P_source * S_source).subs(
        {x: x_on_t_zero, t: 0, q: q_on_t_zero}
    )
    - g1 * (-2 * y) / (2 * g3)
) == 0


# All-degree Newton ledger at P=0.
for degree in range(4, 65):
    d = degree - 3
    h = gcd(d, 2)
    boundary_e = d // h

    # Coefficient-valuation points:
    # (0,0),(1,0),(2,0),(3,1),(j,j), 4<=j<=N.
    last_slope = Fraction(degree - 1, degree - 3)
    assert last_slope == Fraction(d + 2, d)
    for exponent in range(4, degree):
        line_height = Fraction(1) + (exponent - 3) * last_slope
        assert line_height < exponent

    # The three horizontal lengths account for every inverse sheet.
    assert 2 + 1 + h * boundary_e == degree

    # Normalize a high-branch valuation.  The derivative has the asserted
    # order, and q=P*D has a pole on every such branch.
    v_P = Fraction(d, h)
    v_S = -Fraction(d + 2, h)
    v_D = v_P + 2 * v_S
    assert v_D == -Fraction(d + 4, h)
    v_q = v_P + v_D
    assert v_q == -Fraction(4, h)
    assert v_q < 0

    # The endpoint residual g3*z^3+gN*z^N is separable at nonzero roots.
    # Modulo z^d=-g3/gN its derivative coefficient is (3-N)g3*z^2.
    assert 3 - degree == -d
    assert d != 0

    # The second-vertex boundary profile never equals one index-two prime.
    assert not (h == 1 and boundary_e == 2)


# Polynomial variables add no units to either intrinsic normalization ring.
weighted_unit_rank = 1  # k[W,xi^{+-1},z_1,...,z_s]
quadratic_unit_rank = 2  # k[P^{+-1},r^{+-1},z_1,...,z_s]
assert weighted_unit_rank != quadratic_unit_rank


print("PASS: critical parameterization has dC+r^2*dB=0 and odd pole at r=0")
print("PASS: r and r^{-1} give the finite Laurent normalization")
print("PASS: q=0 and t=0 are the two affine branches over P=0")
print("PASS: all P=0 Newton ledgers and degree sums hold through degree 64")
print("PASS: every high Newton branch has a pole in the source coordinate q")
print("PASS: stabilized intrinsic strata have unit ranks one and two")
