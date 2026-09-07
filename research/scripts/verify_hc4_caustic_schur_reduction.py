#!/usr/bin/env python3
"""Exact caustic data for the four unresolved singular PC(2) charts.

For a graph chart (q,m), put

    A = dq/dh,  C = dm/dh,  Delta = det(A),  S = C*A^{-1}.

The Lagrangian identity makes S symmetric.  In charts 0010, 0111, 1000,
and 1101, Delta has a simple zero at the graph origin and A has rank three.
Consequently the polar part of S is rank one.  This script verifies the
exact local normal data and the stronger global X-residues in charts 0010
and 1000.

It also checks an all-degree collision-incidence obstruction: in charts
1100 and 1101 the two symmetric collision points have identical q-values
and distinct m-values, so no shear m+grad(V)(q), for any polynomial V, can
retain that pair.
"""

from __future__ import annotations

import runpy

import sympy as sp


graph = runpy.run_path("scripts/search_hc4_graph_polarizations.py")
h = graph["h_variables"]
Q = list(graph["position_coordinates_h"])
M = list(graph["momentum_coordinates_h"])
position_collision_values = graph["position_collision_values"]
momentum_collision_values = graph["momentum_collision_values"]

X, Y, W, D = h
origin = {variable: 0 for variable in h}


def chart_data(
    chart: str,
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix, sp.Expr]:
    mask = tuple(int(bit) for bit in chart)
    q = sp.Matrix([M[index] if mask[index] else Q[index] for index in range(4)])
    m = sp.Matrix(
        [-Q[index] if mask[index] else M[index] for index in range(4)]
    )
    A = q.jacobian(h)
    C = m.jacobian(h)
    return q, m, A, C, sp.factor(A.det())


def collision_chart_values(
    chart: str, point_index: int
) -> tuple[sp.Matrix, sp.Matrix]:
    mask = tuple(int(bit) for bit in chart)
    q = sp.Matrix(
        [
            momentum_collision_values[point_index][index]
            if mask[index]
            else position_collision_values[point_index][index]
            for index in range(4)
        ]
    )
    m = sp.Matrix(
        [
            -position_collision_values[point_index][index]
            if mask[index]
            else momentum_collision_values[point_index][index]
            for index in range(4)
        ]
    )
    return q, m


local_data = {
    "0010": (
        (-sp.Rational(8, 3), 0, 0, 0),
        -sp.Rational(1, 2)
        * sp.Matrix([0, 1, 0, -1])
        * sp.Matrix([0, 1, 0, -1]).T,
    ),
    "0111": (
        (0, sp.Rational(5, 2), 0, 0),
        sp.Rational(1, 2)
        * sp.Matrix([0, 1, 0, 1])
        * sp.Matrix([0, 1, 0, 1]).T,
    ),
    "1000": (
        (-sp.Rational(32, 3), 0, 0, 0),
        -2
        * sp.Matrix([0, 1, 0, -1])
        * sp.Matrix([0, 1, 0, -1]).T,
    ),
    "1101": (
        (0, 10, 0, 0),
        2
        * sp.Matrix([0, 1, 0, 1])
        * sp.Matrix([0, 1, 0, 1]).T,
    ),
}

for chart, (expected_gradient, expected_polar_numerator) in local_data.items():
    _, _, A, C, delta = chart_data(chart)
    A0 = A.subs(origin)
    C0 = C.subs(origin)
    polar_numerator = (C0 * A0.adjugate()).applyfunc(sp.factor)
    delta_gradient = tuple(
        sp.factor(sp.diff(delta, variable).subs(origin)) for variable in h
    )

    assert delta.subs(origin) == 0
    assert delta_gradient == expected_gradient
    assert A0.rank() == 3
    assert polar_numerator == expected_polar_numerator
    assert polar_numerator == polar_numerator.T
    assert polar_numerator.rank() == 1


def x_residue(chart: str) -> sp.Matrix:
    _, _, A, C, delta = chart_data(chart)
    numerator_at_x_zero = (C * A.adjugate()).subs(X, 0)
    normal_factor = sp.cancel(delta / X).subs(X, 0)
    return numerator_at_x_zero.applyfunc(
        lambda entry: sp.factor(sp.cancel(entry / normal_factor))
    )


r = 2 * W + 9 * Y**2
v_0010 = sp.Matrix([4 * r / 3, 1, 0, -1])
v_1000 = sp.Matrix([0, 1, 2 * r / 3, -1])

residue_0010 = x_residue("0010")
residue_1000 = x_residue("1000")
assert (
    residue_0010 - sp.Rational(3, 16) * v_0010 * v_0010.T
).applyfunc(sp.factor) == sp.zeros(4)
assert (
    residue_1000 - sp.Rational(3, 16) * v_1000 * v_1000.T
).applyfunc(sp.factor) == sp.zeros(4)
assert residue_0010.rank() == 1
assert residue_1000.rank() == 1


for chart in ("1100", "1101"):
    q_plus, m_plus = collision_chart_values(chart, 1)
    q_minus, m_minus = collision_chart_values(chart, 2)
    assert q_plus == q_minus
    assert m_plus != m_minus


print("PASS: the four caustic charts have smooth rank-three projection drops")
print("PASS: their local Schur polar numerators are symmetric of rank one")
print("PASS: charts 0010 and 1000 have the asserted global rank-one X-residues")
print(
    "PASS: no single shear in charts 1100 or 1101 can retain the "
    "symmetric collision pair"
)
print(
    "SCOPE: the rank-one caustic equations do not yet rule out arbitrary "
    "polynomial Hessian data"
)
