#!/usr/bin/env python3
"""Exact checks for the three-puncture double-incidence core."""

import math

import sympy as sp


s, t, u, v, xi = sp.symbols("s t u v xi")
ell = t - s
D0 = 1 - s * u
D1 = 1 - t * v


def diagonal_primitive(a: int, b: int) -> sp.Expr:
    """Return the explicit polynomial in equation (1)."""

    total = 0
    for i in range(a + 1):
        for j in range(b + 1):
            for q in range(j + 1):
                total += (
                    math.comb(a, i)
                    * math.comb(b, j)
                    * math.comb(j, q)
                    * (-u) ** i
                    * (-v) ** j
                    * ell ** (j - q)
                    * s ** (i + q + 1)
                    / sp.Integer(i + q + 1)
                )
    return sp.expand(total)


for a in range(1, 5):
    for b in range(1, 5):
        R = diagonal_primitive(a, b)

        integrand = sp.expand(
            (1 - xi * u) ** a * (1 - (xi + ell) * v) ** b
        )
        integral = sp.integrate(integrand, (xi, 0, s))
        assert sp.expand(R - integral) == 0

        diagonal_derivative = sp.diff(R, s) + sp.diff(R, t)
        assert sp.factor(diagonal_derivative - D0**a * D1**b) == 0

        outputs = sp.Matrix([t - s + 1, u, v, R])
        jacobian = sp.factor(outputs.jacobian((s, t, u, v)).det())
        assert sp.factor(jacobian + D0**a * D1**b) == 0

# The selected complete intersection is smooth and has the asserted
# parametrization.  A nonzero 3-by-3 minor suffices for the rank check.
L = t - s + 1
constraint_jacobian = sp.Matrix([D0, D1, L]).jacobian((s, t, u, v))
minor = sp.factor(constraint_jacobian[:, (1, 2, 3)].det())
minor_on_curve = sp.factor(
    minor.subs({t: s - 1, u: 1 / s, v: 1 / (s - 1)})
)
assert minor_on_curve == s * (s - 1)

curve_substitution = {t: s - 1, u: 1 / s, v: 1 / (s - 1)}
assert sp.factor(D0.subs(curve_substitution)) == 0
assert sp.factor(D1.subs(curve_substitution)) == 0
assert sp.factor(L.subs(curve_substitution)) == 0

# The valuations of s and s-1 at (0,1,infinity) are independent and span
# the degree-zero puncture lattice.
valuation_matrix = sp.Matrix([[1, 0, -1], [0, 1, -1]])
assert valuation_matrix.rank() == 2
assert all(sum(valuation_matrix.row(i)) == 0 for i in range(2))

# Universal affine-linear completion identity.
r, z = sp.symbols("r z")
A0 = sp.Function("A0")(r)
B0 = sp.Function("B0")(r)
aa = sp.Function("aa")(r)
bb = sp.Function("bb")(r)
A = A0 + aa * z
B = B0 + bb * z
affine_linear_jacobian = sp.expand(
    sp.Matrix([A, B]).jacobian((r, z)).det()
)
expected = sp.expand(
    sp.diff(A0, r) * bb
    - aa * sp.diff(B0, r)
    + z * (sp.diff(aa, r) * bb - aa * sp.diff(bb, r))
)
assert sp.simplify(affine_linear_jacobian - expected) == 0

print(
    "PASS: double-incidence determinant, three-puncture normalization, "
    "and affine-linear completion gate"
)
