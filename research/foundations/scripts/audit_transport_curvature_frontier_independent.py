#!/usr/bin/env python3
"""Independent direct-expansion audit for the transport-curvature frontier.

Unlike ``verify_transport_curvature_frontier.py``, this checker does not use
or reconstruct the block inverse and determinant-lemma proof.  It builds
ordinary four-variable Hessians from explicit polynomials and compares their
determinants with the closed formulas in the note.
"""
from __future__ import annotations

import sympy as sp

x, y, t, m = sp.symbols("x y t m")
variables = (x, y, t, m)


def direct_det(poly: sp.Expr) -> sp.Expr:
    return sp.factor(sp.hessian(sp.expand(poly), variables).det())


def dominant_rhs(a: sp.Expr, b: sp.Expr, degree: int, scalar: int) -> sp.Expr:
    c = sp.Matrix([a, b])
    jac = c.jacobian((x, y))
    radial = sp.simplify(jac.inv() * c)
    linear = sp.expand(a * t + b * m)
    source_hessian = t * sp.hessian(a, (x, y)) + m * sp.hessian(b, (x, y))
    quadratic = sp.simplify((radial.T * source_hessian * radial)[0])
    return sp.factor(
        (scalar * degree) ** 4
        * linear ** (4 * degree - 5)
        * jac.det() ** 2
        * ((2 * degree - 1) * linear - (degree - 1) * quadratic)
    )


# Direct test of (2.1), with a nonlinear triangular coefficient map, d=3.
a1 = x
b1 = y + x**2
phi1 = 2 * (a1 * t + b1 * m) ** 3
assert sp.factor(direct_det(phi1) - dominant_rhs(a1, b1, 3, 2)) == 0

# A second, unrelated nonlinear coefficient map, d=2.
a2 = x + y**2
b2 = y
phi2 = -3 * (a2 * t + b2 * m) ** 2
assert sp.factor(direct_det(phi2) - dominant_rhs(a2, b2, 2, -3)) == 0

# Direct one-latitude test of (5.6) on the nonlinear latitude h=xy.
h = x * y
A = sp.Integer(1)
B = h
mu = sp.Integer(1)
degree = 2
linear = A * t + B * m
phi3 = mu * linear**degree
bordered = sp.factor(
    sp.diff(h, x) ** 2 * sp.diff(h, y, 2)
    - 2 * sp.diff(h, x) * sp.diff(h, y) * sp.diff(h, x, y)
    + sp.diff(h, y) ** 2 * sp.diff(h, x, 2)
)
wronskian = sp.Integer(1)
linear_h = m
rhs3 = sp.factor(
    -degree**3
    * (degree - 1)
    * mu**3
    * bordered
    * wronskian**2
    * linear ** (4 * degree - 5)
    * (degree * mu * linear_h)
)
assert sp.factor(direct_det(phi3) - rhs3) == 0
assert bordered == -2 * x * y

# An actual affine latitude must give zero determinant, despite nontrivial
# amplitude and a moving projective line over that latitude.
ell = x + 2 * y
phi4 = (1 + ell**2) * ((1 + ell) * t + (1 + ell**3) * m) ** 3
assert direct_det(phi4) == 0

# Rank collapse alone is not enough: this profile has coefficient algebra
# k[xy] but a nonzero nonlinear-latitude curvature.
phi5 = x**2 * (y * t + x * y**2 * m) ** 2
expected5 = 32 * x**8 * y**8 * (t + m * x * y) ** 3 * (t + 2 * m * x * y)
assert sp.factor(direct_det(phi5) - expected5) == 0

# Direct fixed instances of the two Cohn-profile axis obstructions.
phi6 = ((1 + x * y) * t + x**2 * m) ** 2
axis6 = sp.Poly(sp.expand(direct_det(phi6).subs({y: 0, m: 0})), t, x)
assert axis6.coeff_monomial(t**4 * x**4) == 128

phi7 = (x**2 * t + (1 + x * y) * m) ** 2
axis7 = sp.factor(direct_det(phi7).subs(t, 0))
assert axis7 != 0

print("PASS independent direct Hessian expansion: dominant formula, d=3")
print("PASS independent direct Hessian expansion: dominant formula, d=2")
print("PASS independent direct Hessian expansion: nonlinear one-latitude formula")
print("PASS affine linear latitude has zero four-variable Hessian determinant")
print("PASS rank collapse without affine latitude retains a nonzero determinant")
print("PASS direct fixed Cohn-profile axis obstructions")
