#!/usr/bin/env python3
"""Boundary Schur chain for the PC(2) chart 1000.

Use adapted chart coordinates

    u0=a,  u1=b,  u2=c,  u3=b+d.

At X=0 the independent projection has image (a,b,0,b), where

    D = -a/2 - W^2/3 - 3*W*b^2 - 29*b^4/4.

Write f=V|_{c=d=0}, g=V_c|_{c=d=0}, and h=V_d|_{c=d=0}.
The remaining normal Hessian entries are V_cc, V_cd, and V_dd.

After the fixed-image substitution, the determinant is quadratic in W.
Its W^2, W^1, and W^0 coefficients are triangular in V_dd, V_cd, and
V_cc, with common pivot

    L = 5*b*f_aa - (f_aa*f_bb-f_ab^2).

Thus L != 0 determines the three normal Hessian entries rationally.  This
is a reduction, not an obstruction: polynomial divisibility and higher
normal-jet compatibility remain to be checked on that branch.
"""

from __future__ import annotations

import runpy

import sympy as sp


graph = runpy.run_path("scripts/search_hc4_graph_polarizations.py")
h_variables = graph["h_variables"]
Q = list(graph["position_coordinates_h"])
M = list(graph["momentum_coordinates_h"])
X, Y, W, D = h_variables

mask = (1, 0, 0, 0)
q = sp.Matrix([M[index] if mask[index] else Q[index] for index in range(4)])
m = sp.Matrix(
    [-Q[index] if mask[index] else M[index] for index in range(4)]
)

a, b = sp.symbols("a b")
fixed_image_D = (
    -a / 2 - W**2 / 3 - 3 * W * b**2 - sp.Rational(29, 4) * b**4
)
boundary_substitution = {X: 0, Y: b, D: fixed_image_D}
assert [
    sp.factor(coordinate.subs(boundary_substitution)) for coordinate in q
] == [a, b, 0, b]

jq_boundary = q.jacobian(h_variables).subs(boundary_substitution)
jm_boundary = m.jacobian(h_variables).subs(boundary_substitution)

(
    f_aa,
    f_ab,
    g_a,
    h_a,
    f_bb,
    g_b,
    h_b,
    V_cc,
    V_cd,
    V_dd,
) = sp.symbols("f_aa f_ab g_a h_a f_bb g_b h_b V_cc V_cd V_dd")

hessian_adapted = sp.Matrix(
    [
        [f_aa, f_ab, g_a, h_a],
        [f_ab, f_bb, g_b, h_b],
        [g_a, g_b, V_cc, V_cd],
        [h_a, h_b, V_cd, V_dd],
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

tangential_hessian_determinant = f_aa * f_bb - f_ab**2
L = 5 * b * f_aa - tangential_hessian_determinant

boundary_three_hessian = sp.Matrix(
    [
        [f_aa, f_ab, h_a],
        [f_ab, f_bb - 5 * b, h_b - b / 2],
        [h_a, h_b - b / 2, V_dd - 5 * b / 64],
    ]
)
assert sp.factor(
    coefficient_2
    + sp.Rational(32, 9) * boundary_three_hessian.det()
) == 0

pivot_V_dd = sp.factor(sp.diff(coefficient_2, V_dd))
pivot_V_cd = sp.factor(sp.diff(coefficient_1, V_cd))
pivot_V_cc = sp.factor(sp.diff(coefficient_0, V_cc))
assert sp.factor(pivot_V_dd - sp.Rational(32, 9) * L) == 0
assert sp.factor(pivot_V_cd - sp.Rational(16, 3) * L) == 0
assert sp.factor(pivot_V_cc - 2 * L) == 0

kappa = sp.symbols("kappa")
forced_V_dd = sp.factor(
    -coefficient_2.subs(V_dd, 0) / pivot_V_dd
)
forced_V_cd = sp.factor(
    -coefficient_1.subs(
        {V_dd: forced_V_dd, V_cd: 0}, simultaneous=True
    )
    / pivot_V_cd
)
forced_V_cc = sp.factor(
    (
        kappa
        - coefficient_0.subs(
            {
                V_dd: forced_V_dd,
                V_cd: forced_V_cd,
                V_cc: 0,
            },
            simultaneous=True,
        )
    )
    / pivot_V_cc
)

assert sp.factor(coefficient_2.subs(V_dd, forced_V_dd)) == 0
assert sp.factor(
    coefficient_1.subs(
        {V_dd: forced_V_dd, V_cd: forced_V_cd}, simultaneous=True
    )
) == 0
assert sp.factor(
    coefficient_0.subs(
        {
            V_dd: forced_V_dd,
            V_cd: forced_V_cd,
            V_cc: forced_V_cc,
        },
        simultaneous=True,
    )
    - kappa
) == 0


def main() -> None:
    print("PASS: the chart 1000 fixed-image boundary is (a,b,0,b)")
    print("PASS: the boundary determinant has W-degree two")
    print("PASS: its W^2 coefficient is a three-variable Hessian determinant")
    print("PASS: W^2, W^1, W^0 solve V_dd, V_cd, V_cc with pivot L")
    print("SCOPE: L != 0 gives rational forced data, not a polynomial solution")
    print("OPEN: test divisibility on L != 0 and classify the branch L=0")


if __name__ == "__main__":
    main()
