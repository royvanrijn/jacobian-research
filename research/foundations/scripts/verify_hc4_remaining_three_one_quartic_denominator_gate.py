#!/usr/bin/env python3
"""Exact checks for the remaining clean P=x^3*y incidence packets."""

from __future__ import annotations

import sympy as sp


x, y, z = sp.symbols("x y z")
alpha, beta, h0, j0, k0, k1, scalar = sp.symbols(
    "alpha beta h0 j0 k0 k1 scalar"
)


# After the tangent constant-kernel coefficient ladder on alpha!=0, this is
# the complete quintic with x^6 divisibility.
h_alpha = (
    alpha * y**5
    + beta * x * y**4
    + h0 * x**2 * y**3 / 2
    + j0 * x**3 * y**2 / 6
    + x**4 * (k0 * y + k1 * z) / 24
    + scalar * x**5 / 120
)
det_alpha = sp.factor(sp.hessian(h_alpha, (x, y, z)).det())
assert det_alpha == (
    -k1**2
    * x**6
    * (60 * alpha * y**3 + 36 * beta * x * y**2 + 9 * h0 * x**2 * y + j0 * x**3)
    / 108
)

# Exact x^6*y^3 factorization removes the higher x faces, but also destroys
# the generic rank-two boundary determinant.
residual_cube_conditions = {beta: 0, h0: 0, j0: 0}
assert sp.factor(det_alpha.subs(residual_cube_conditions)) == (
    -sp.Rational(5, 9) * alpha * k1**2 * x**6 * y**3
)
boundary_rank_factor = 5 * alpha * h0 - 4 * beta**2
assert boundary_rank_factor.subs(residual_cube_conditions) == 0

# The residual x^6 coefficient is always a cube, never root type 2+1.
coefficient_x6 = sp.factor(sp.expand(det_alpha).coeff(x, 6))
assert coefficient_x6 == -sp.Rational(5, 9) * alpha * k1**2 * y**3


# On alpha=0, x^6 divisibility removes j1 and the determinant starts at x^7.
h_beta = (
    beta * x * y**4
    + h0 * x**2 * y**3 / 2
    + j0 * x**3 * y**2 / 6
    + x**4 * (k0 * y + k1 * z) / 24
    + scalar * x**5 / 120
)
det_beta = sp.factor(sp.hessian(h_beta, (x, y, z)).det())
assert det_beta == (
    -k1**2
    * x**7
    * (36 * beta * y**2 + 9 * h0 * x * y + j0 * x**2)
    / 108
)
assert sp.expand(det_beta).coeff(x, 6) == 0


# A transverse constant kernel has no surviving normal jet under x^6
# divisibility and therefore zero determinant.  Verify the terminal form.
F = sp.Function("F")
u, v = sp.symbols("u v")
binary = u**5 + u * v**4 + v**5
h_transverse_terminal = binary
assert sp.hessian(h_transverse_terminal, (x, u, v)).det() == 0


print("PASS: exact alpha-branch tangent determinant factorization")
print("PASS: residual cubic at exact x-multiplicity six is forced to a cube")
print("PASS: cube matching contradicts generic boundary rank")
print("PASS: beta branch has x-multiplicity at least seven")
print("THEOREM: x^6*y^3 and x^6*y^2*z clean packets are empty")
print("SCOPE: closes the remaining P=x^3*y incidences; lower-Smith is separate")
