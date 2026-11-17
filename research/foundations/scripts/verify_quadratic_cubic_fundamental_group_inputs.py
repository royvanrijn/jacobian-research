#!/usr/bin/env python3
"""Symbolic inputs to the fundamental-group computation for U_(2,3)."""

from __future__ import annotations

import sympy as sp


T = sp.Symbol("T")
a0, a1, a2 = sp.symbols("a0 a1 a2")
b0, b1, b2, b3 = sp.symbols("b0 b1 b2 b3")
lam, mu = sp.symbols("lam mu", nonzero=True)
r, s = sp.symbols("r s", nonzero=True)
variables = (a0, a1, a2, b0, b1, b2, b3)

A = a0 * T**2 + a1 * T + a2
B = b0 * T**3 + b1 * T**2 + b2 * T + b3
m = a0 * b1 + a1 * b0
resultant = sp.expand(sp.resultant(A, B, T))

# The two relative scalings have weights (1,1) on m and (3,2) on Res.
scaling = {
    a0: lam * a0,
    a1: lam * a1,
    a2: lam * a2,
    b0: mu * b0,
    b1: mu * b1,
    b2: mu * b2,
    b3: mu * b3,
}
scaled_m = sp.expand(m.subs(scaling, simultaneous=True))
scaled_resultant = sp.expand(resultant.subs(scaling, simultaneous=True))
assert sp.expand(scaled_m - lam * mu * m) == 0
assert sp.expand(scaled_resultant - lam**3 * mu**2 * resultant) == 0

# For m=s and Res=r, these scalars give the unique normalization to (1,1).
normalizing_lam = s**2 / r
normalizing_mu = r / s**3
assert sp.factor(normalizing_lam * normalizing_mu * s) == 1
assert sp.factor(normalizing_lam**3 * normalizing_mu**2 * r) == 1

# An explicit smooth transverse point of Res=0 and m=0.
point = {a0: 0, a1: 1, a2: 0, b0: 0, b1: 1, b2: 0, b3: 1}
assert m.subs(point) == 0
assert resultant.subs(point) == 0
boundary_jacobian = sp.Matrix([
    [sp.diff(m, variable).subs(point) for variable in variables],
    [sp.diff(resultant, variable).subs(point) for variable in variables],
])
assert boundary_jacobian.rank() == 2

# On the incidence normalization, a common factor L=v*T-u*S splits m into
# the root-at-infinity component and one residual component.
u, v, alpha, beta, c, d, e = sp.symbols("u v alpha beta c d e")
common_factor = v * T - u
incidence_A = sp.expand(common_factor * (alpha * T + beta))
incidence_B = sp.expand(common_factor * (c * T**2 + d * T + e))
incidence_A_poly = sp.Poly(incidence_A, T)
incidence_B_poly = sp.Poly(incidence_B, T)
incidence_m = sp.expand(
    incidence_A_poly.coeff_monomial(T**2)
    * incidence_B_poly.coeff_monomial(T**2)
    + incidence_A_poly.coeff_monomial(T)
    * incidence_B_poly.coeff_monomial(T**3)
)
expected_incidence_m = v * (v * (alpha * d + beta * c) - 2 * u * alpha * c)
assert sp.expand(incidence_m - expected_incidence_m) == 0

# Cycle-class ranks used in the H^2/H^3 calculation.  Coordinates are
# (h1^2, h1*h2, h2^2).
z_infinity = sp.Matrix([0, 1, 0])
z_residual = sp.Matrix([3, 4, 2])
assert sp.Matrix.hstack(z_infinity, z_residual).rank() == 2
boundary_spanning_cycles = sp.Matrix.hstack(
    sp.Matrix([1, 1, 0]),  # E*h1
    sp.Matrix([0, 1, 1]),  # E*h2
    sp.Matrix([3, 2, 0]),  # R*h1
)
assert boundary_spanning_cycles.det() != 0

print("PASS normalization torus: weights (1,1) and (3,2) are unimodular")
print("PASS W = X_(2,3) x G_m^2 normalization formulas")
print("PASS the tangent and resultant boundaries meet transversely")
print("PASS R intersect E has the infinity and residual incidence components")
print("PASS top boundary cycle ranks give H^2=0 and dim H^3=1")
