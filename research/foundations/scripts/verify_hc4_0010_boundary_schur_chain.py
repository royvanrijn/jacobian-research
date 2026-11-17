#!/usr/bin/env python3
"""Boundary Schur chain for the unresolved PC(2) chart 0010.

Use adapted chart coordinates

    u0=a,  u1=b,  u2=c,  u3=b+d.

The caustic image is the plane a=d=0.  Write the two-jet of a potential V
on that plane as

    f_bb, f_bc, f_cc,
    g_b=V_ab, g_c=V_ac,
    h_b=V_bd, h_c=V_cd,
    A=V_aa, B=V_ad, C=V_dd.

Restrict det d(m+grad(V)(q)) to X=0.  It is quadratic in the collapsed
fiber coordinate W.  This script proves that its W^2, W^1, and W^0
coefficients are triangular in C, B, and A, with one common pivot

    L = 5*b*f_cc - (f_bb*f_cc-f_bc^2).

On L != 0 the constant-determinant boundary equation therefore determines
C, then B, then A rationally.  This is not an obstruction by itself: the
remaining gates are polynomial divisibility, compatibility with higher
normal jets, and the separate characteristic branch L=0.
"""

from __future__ import annotations

import runpy

import sympy as sp


graph = runpy.run_path("scripts/search_hc4_graph_polarizations.py")
h_variables = graph["h_variables"]
Q = list(graph["position_coordinates_h"])
M = list(graph["momentum_coordinates_h"])
X, Y, W, D = h_variables

mask = (0, 0, 1, 0)
q = sp.Matrix([M[index] if mask[index] else Q[index] for index in range(4)])
m = sp.Matrix(
    [-Q[index] if mask[index] else M[index] for index in range(4)]
)
jq_boundary = q.jacobian(h_variables).subs(X, 0)
jm_boundary = m.jacobian(h_variables).subs(X, 0)

A, g_b, g_c, B, f_bb, f_bc, h_b, f_cc, h_c, C = sp.symbols(
    "A g_b g_c B f_bb f_bc h_b f_cc h_c C"
)

hessian_adapted = sp.Matrix(
    [
        [A, g_b, g_c, B],
        [g_b, f_bb, f_bc, h_b],
        [g_c, f_bc, f_cc, h_c],
        [B, h_b, h_c, C],
    ]
)

# u=(a,b,c,b+d).  Hessians transform by H_z=P^t H_u P, so recover H_u.
coordinate_matrix = sp.Matrix(
    [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 1, 0, 1],
    ]
)
hessian_original = (
    coordinate_matrix.inv().T
    * hessian_adapted
    * coordinate_matrix.inv()
)

boundary_determinant = sp.Poly(
    sp.expand(
        (
            jm_boundary + hessian_original * jq_boundary
        ).det(method="berkowitz")
    ),
    W,
)
assert boundary_determinant.degree() == 2
coefficient_2 = sp.factor(boundary_determinant.nth(2))
coefficient_1 = sp.factor(boundary_determinant.nth(1))
coefficient_0 = sp.factor(boundary_determinant.nth(0))

tangential_hessian_determinant = f_bb * f_cc - f_bc**2
L = 5 * Y * f_cc - tangential_hessian_determinant

boundary_three_hessian = sp.Matrix(
    [
        [f_bb - 5 * Y, f_bc, h_b - Y / 2],
        [f_bc, f_cc, h_c],
        [h_b - Y / 2, h_c, C + 13 * Y / 64],
    ]
)
assert sp.factor(
    coefficient_2
    + sp.Rational(32, 9) * boundary_three_hessian.det()
) == 0

pivot_C = sp.factor(sp.diff(coefficient_2, C))
pivot_B = sp.factor(sp.diff(coefficient_1, B))
pivot_A = sp.factor(sp.diff(coefficient_0, A))
assert sp.factor(pivot_C - sp.Rational(32, 9) * L) == 0
assert sp.factor(pivot_B - sp.Rational(8, 3) * L) == 0
assert sp.factor(pivot_A - sp.Rational(1, 2) * L) == 0

kappa = sp.symbols("kappa")
forced_C = sp.factor(-coefficient_2.subs(C, 0) / pivot_C)
forced_B = sp.factor(
    -coefficient_1.subs({C: forced_C, B: 0}, simultaneous=True)
    / pivot_B
)
forced_A = sp.factor(
    (
        kappa
        - coefficient_0.subs(
            {C: forced_C, B: forced_B, A: 0}, simultaneous=True
        )
    )
    / pivot_A
)

assert sp.factor(coefficient_2.subs(C, forced_C)) == 0
assert sp.factor(
    coefficient_1.subs({C: forced_C, B: forced_B}, simultaneous=True)
) == 0
assert sp.factor(
    coefficient_0.subs(
        {C: forced_C, B: forced_B, A: forced_A},
        simultaneous=True,
    )
    - kappa
) == 0

def main() -> None:
    print("PASS: the chart 0010 boundary determinant has W-degree two")
    print("PASS: its W^2 coefficient is a three-variable Hessian determinant")
    print("PASS: W^2, W^1, W^0 solve C, B, A with one common pivot L")
    print("SCOPE: L != 0 gives rational forced data, not a polynomial solution")
    print("OPEN: test divisibility of the forced data and the branch L=0")


if __name__ == "__main__":
    main()
