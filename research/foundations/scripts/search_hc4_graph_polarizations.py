#!/usr/bin/env python3
"""Exact first search for an HC(4) graph polarization of the PC(2) map.

The graph of the symplectic map G=(R,T,D,S) is Lagrangian in

    (A^4_source x A^4_target, -omega_source + omega_target).

In ambient Darboux pairs

    (x,-p), (q,-z), (R,D), (T,S),

each of the 2^4 coordinate polarizations chooses either the first coordinate
as an independent variable and the second as its gradient coordinate, or
the second as independent and the negative of the first as complementary.

For every polarization this script checks exact necessary conditions for an
HC(4) counterexample:

* constant nonzero Jacobian of the independent projection;
* constant nonzero Hessian determinant, computed as det(dB)/det(dA);
* symmetry of dB/dA on the fraction field; and
* survival of a collision from the known complete three-point fiber.

Polynomial inversion and integration are only needed for candidates passing
these gates.
"""

from __future__ import annotations

from itertools import combinations, product

import sympy as sp


x, q, p, z = sp.symbols("x q p z")
source_variables = (x, q, p, z)


def poisson(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    return sp.expand(
        sp.diff(f, p) * sp.diff(g, x)
        - sp.diff(f, x) * sp.diff(g, p)
        + sp.diff(f, z) * sp.diff(g, q)
        - sp.diff(f, q) * sp.diff(g, z)
    )


X = x
Q = q
Z = 3 * x**2 * p + (2 - 6 * x * q) * z
E = (1 + 3 * x * q) * p / 2 - 3 * q**2 * z
W = sp.expand(Z - 9 * Q**2)
Y = sp.expand(Q - X * W / 3)
U = sp.expand(1 + X * Y)

R = sp.expand(2 * X - 3 * X**2 * Y - X**3 * W)
S = sp.expand((U**3 * W + Y**2 * U * (4 + 3 * X * Y)) / 2)
T = sp.expand(Y + 3 * X * U**2 * W + 3 * X * Y**2 * (4 + 3 * X * Y))
K = -(
    10 * W**3 * X**2
    + 90 * W**2 * X * Y
    + 20 * W**2
    - 18 * W * X**3 * Y**5
    - 90 * W * X**2 * Y**4
    - 180 * W * X * Y**3
    + 90 * W * Y**2
    - 54 * X**2 * Y**6
    - 234 * X * Y**5
    - 375 * Y**4
) / 60
D = sp.expand(E + K)

position_coordinates = (x, q, R, T)
momentum_coordinates = (-p, -z, D, S)
pair_names = (("x", "-p"), ("q", "-z"), ("R", "D"), ("T", "S"))

collision_points = (
    (sp.Rational(0), sp.Rational(0), sp.Rational(1, 24), sp.Rational(-1, 8)),
    (sp.Rational(1), sp.Rational(2, 3), sp.Rational(247, 96), sp.Rational(-89, 64)),
    (sp.Rational(-1), sp.Rational(-2, 3), sp.Rational(247, 96), sp.Rational(-89, 64)),
)

# The Darboux trivialization H=(X,Y,W,D) is a polynomial automorphism.  It is
# much faster to compute projection determinants in these coordinates than
# after expanding D back into (x,q,p,z); the determinant quotient det(dB)/det(dA)
# is unchanged.  These formulas are the exact inverse of H.
Xh, Yh, Wh, Dh = sp.symbols("Xh Yh Wh Dh")
Qh = Yh + Xh * Wh / 3
Zh = Wh + 9 * Qh**2
Uh = 1 + Xh * Yh
Rh = 2 * Xh - 3 * Xh**2 * Yh - Xh**3 * Wh
Sh = (Uh**3 * Wh + Yh**2 * Uh * (4 + 3 * Xh * Yh)) / 2
Th = Yh + 3 * Xh * Uh**2 * Wh + 3 * Xh * Yh**2 * (4 + 3 * Xh * Yh)
Kh = -(
    10 * Wh**3 * Xh**2
    + 90 * Wh**2 * Xh * Yh
    + 20 * Wh**2
    - 18 * Wh * Xh**3 * Yh**5
    - 90 * Wh * Xh**2 * Yh**4
    - 180 * Wh * Xh * Yh**3
    + 90 * Wh * Yh**2
    - 54 * Xh**2 * Yh**6
    - 234 * Xh * Yh**5
    - 375 * Yh**4
) / 60
Eh = Dh - Kh
ph = -6 * Eh * Xh * Qh + 2 * Eh + 3 * Zh * Qh**2
zh = -3 * Eh * Xh**2 + 3 * Zh * Xh * Qh / 2 + Zh / 2
h_variables = (Xh, Yh, Wh, Dh)
position_coordinates_h = (Xh, Qh, Rh, Th)
momentum_coordinates_h = (-ph, -zh, Dh, Sh)

# Values of the four ambient Darboux pairs on the certified collision fiber.
# These are inexpensive exact data already checked by the completion script.
position_collision_values = (
    (sp.Rational(0), sp.Rational(0), sp.Rational(0), sp.Rational(0)),
    (sp.Rational(1), sp.Rational(2, 3), sp.Rational(0), sp.Rational(0)),
    (sp.Rational(-1), sp.Rational(-2, 3), sp.Rational(0), sp.Rational(0)),
)
momentum_collision_values = (
    (
        sp.Rational(-1, 24),
        sp.Rational(1, 8),
        sp.Rational(0),
        sp.Rational(-1, 8),
    ),
    (
        sp.Rational(-247, 96),
        sp.Rational(89, 64),
        sp.Rational(0),
        sp.Rational(-1, 8),
    ),
    (
        sp.Rational(-247, 96),
        sp.Rational(89, 64),
        sp.Rational(0),
        sp.Rational(-1, 8),
    ),
)


def is_nonzero_constant(expr: sp.Expr) -> bool:
    return bool(expr != 0 and not expr.free_symbols)


def at_point(expr: sp.Expr, point: tuple[sp.Rational, ...]) -> sp.Expr:
    return sp.factor(expr.subs(dict(zip(source_variables, point, strict=True))))


def main() -> None:
    # The graph must be Lagrangian before any polarization test is meaningful.
    graph_form = sp.zeros(4, 4)
    for a, b in zip(position_coordinates, momentum_coordinates, strict=True):
        da = sp.Matrix([sp.diff(a, variable) for variable in source_variables])
        db = sp.Matrix([sp.diff(b, variable) for variable in source_variables])
        graph_form += da * db.T - db * da.T
    assert graph_form.applyfunc(sp.expand) == sp.zeros(4, 4)

    rows: list[dict[str, object]] = []
    survivors: list[dict[str, object]] = []
    for mask in product((0, 1), repeat=4):
        independent = []
        complementary = []
        independent_names = []
        complementary_names = []
        for bit, a, b, names in zip(
            mask,
            position_coordinates_h,
            momentum_coordinates_h,
            pair_names,
            strict=True,
        ):
            if bit == 0:
                independent.append(a)
                complementary.append(b)
                independent_names.append(names[0])
                complementary_names.append(names[1])
            else:
                independent.append(b)
                complementary.append(-a)
                independent_names.append(names[1])
                complementary_names.append(f"-{names[0]}")

        complementary_values = []
        for point_index in range(len(collision_points)):
            values = []
            for pair_index, bit in enumerate(mask):
                if bit == 0:
                    values.append(momentum_collision_values[point_index][pair_index])
                else:
                    values.append(-position_collision_values[point_index][pair_index])
            complementary_values.append(tuple(values))
        collision_pairs = [
            pair
            for pair in combinations(range(len(collision_points)), 2)
            if complementary_values[pair[0]] == complementary_values[pair[1]]
        ]

        # Collision survival is necessary and eliminates 12 cases before any
        # high-degree determinant is formed.
        det_a: sp.Expr | str = "not needed"
        det_b: sp.Expr | str = "not needed"
        projection_keller = False
        hessian_constant = False
        symmetric = False
        if collision_pairs:
            jac_a = sp.Matrix(independent).jacobian(h_variables)
            jac_b = sp.Matrix(complementary).jacobian(h_variables)
            det_a = sp.factor(jac_a.det(method="berkowitz"))
            projection_keller = is_nonzero_constant(det_a)
            if projection_keller:
                det_b = sp.factor(jac_b.det(method="berkowitz"))
                hessian_constant = is_nonzero_constant(sp.cancel(det_b / det_a))
                relative_jacobian = (jac_b * jac_a.inv()).applyfunc(sp.cancel)
                symmetric = all(
                    sp.cancel(relative_jacobian[i, j] - relative_jacobian[j, i]) == 0
                    for i in range(4)
                    for j in range(i)
                )

        row = {
            "mask": "".join(map(str, mask)),
            "independent": tuple(independent_names),
            "complementary": tuple(complementary_names),
            "det_a": det_a,
            "det_b": det_b,
            "projection_keller": projection_keller,
            "hessian_constant": hessian_constant,
            "symmetric": symmetric,
            "collision_pairs": collision_pairs,
        }
        rows.append(row)
        if projection_keller and hessian_constant and symmetric and collision_pairs:
            survivors.append(row)

    print("mask  independent       det(dA)  det(dB)  symmetric  collisions")
    for row in rows:
        print(
            f"{row['mask']}  {','.join(row['independent']):<18} "
            f"{str(row['det_a']):<9} {str(row['det_b']):<9} "
            f"{str(row['symmetric']):<9} {row['collision_pairs']}"
        )

    print()
    print(f"coordinate polarizations checked: {len(rows)}")
    print(f"full-gate survivors: {len(survivors)}")
    for row in survivors:
        print(row)


if __name__ == "__main__":
    main()
