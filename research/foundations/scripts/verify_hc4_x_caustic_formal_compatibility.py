#!/usr/bin/env python3
"""Formal compatibility on the generic X-caustic branches 0010 and 1000.

On the fixed-image caustic plane c=d=0, the boundary Schur equations solve

    V_dd, V_cd, V_cc

rationally in the tangential Cauchy two-jet whenever

    L = 5*b*f_aa - (f_aa*f_bb-f_ab^2)

is nonzero.  Differentiate those forced entries in a and b.  Together with
the derivatives of f=V, g=V_c, and h=V_d, this determines every symmetric
third derivative having at least one tangential index.  The only new third
normal derivatives are

    n0=V_ccc, n1=V_ccd, n2=V_cdd, n3=V_ddd.

Expand the full graph determinant to first order in the source-normal
coordinate X while holding the boundary source parameters fixed.  After the
order-zero Schur solution is substituted and denominators are cleared, the
X coefficient has W-degree three.  Its W^3,W^2,W^1,W^0 equations are
triangular in n3,n2,n1,n0, with nonzero pivots proportional to L^4.

The boundary linearization has a stronger binary-form description.  If

    tau=2*W+9*b^2,  v=(2,4*tau/3),

then for a normal Hessian variation Q,

    d(det)_H(Q) = L*Q(v,v)/2.

At prolongation order r, the new pure-normal tensor T has order r+2 and
the principal symbol is L*T(v,...,v)/2.  Since the second component of v
is affine-linear and nonconstant in W, this maps binary forms of degree
r+2 isomorphically to W-polynomials of degree at most r+2.  Consequently
every formal normal order is recursively solvable after L is inverted.

Both X-caustic charts have this binary-form symbol, with first normal
velocity component 1 in chart 0010 and 2 in chart 1000.  The checker
verifies both boundary identities.  In chart 1000 it additionally verifies
the complete first prolongation and the five-dimensional second principal
symbol.  It also extracts and factors the first necessary divisibility
condition for polynomial extendability.  Formal solvability does not imply
polynomial extendability.
"""

from __future__ import annotations

from itertools import permutations
import runpy

import sympy as sp


boundary = runpy.run_path("scripts/verify_hc4_1000_boundary_schur_chain.py")
boundary_0010 = runpy.run_path("scripts/verify_hc4_0010_boundary_schur_chain.py")
graph = runpy.run_path("scripts/search_hc4_graph_polarizations.py")

X, Y, W, D = graph["h_variables"]
h_variables = graph["h_variables"]
q = boundary["q"]
m = boundary["m"]
a = boundary["a"]
b = boundary["b"]

f_aa = boundary["f_aa"]
f_ab = boundary["f_ab"]
f_bb = boundary["f_bb"]
g_a = boundary["g_a"]
g_b = boundary["g_b"]
h_a = boundary["h_a"]
h_b = boundary["h_b"]
kappa = boundary["kappa"]
L = boundary["L"]

V_cc = boundary["forced_V_cc"]
V_cd = boundary["forced_V_cd"]
V_dd = boundary["forced_V_dd"]

tau = 2 * W + 9 * b**2
normal_velocity_pair = (sp.Integer(2), 4 * tau / 3)
boundary_determinant = boundary["boundary_determinant"].as_expr()
normal_linearization = (
    sp.diff(boundary_determinant, boundary["V_cc"]),
    sp.diff(boundary_determinant, boundary["V_cd"]),
    sp.diff(boundary_determinant, boundary["V_dd"]),
)
expected_linearization = (
    L * normal_velocity_pair[0] ** 2 / 2,
    L * normal_velocity_pair[0] * normal_velocity_pair[1],
    L * normal_velocity_pair[1] ** 2 / 2,
)
assert all(
    sp.factor(got - expected) == 0
    for got, expected in zip(normal_linearization, expected_linearization, strict=True)
)

tau_0010 = 2 * W + 9 * boundary_0010["Y"] ** 2
normal_velocity_0010 = (sp.Integer(1), 4 * tau_0010 / 3)
boundary_determinant_0010 = boundary_0010["boundary_determinant"].as_expr()
normal_linearization_0010 = (
    sp.diff(boundary_determinant_0010, boundary_0010["A"]),
    sp.diff(boundary_determinant_0010, boundary_0010["B"]),
    sp.diff(boundary_determinant_0010, boundary_0010["C"]),
)
expected_linearization_0010 = (
    boundary_0010["L"] * normal_velocity_0010[0] ** 2 / 2,
    boundary_0010["L"] * normal_velocity_0010[0] * normal_velocity_0010[1],
    boundary_0010["L"] * normal_velocity_0010[1] ** 2 / 2,
)
assert all(
    sp.factor(got - expected) == 0
    for got, expected in zip(
        normal_linearization_0010,
        expected_linearization_0010,
        strict=True,
    )
)

(
    f_aaa,
    f_aab,
    f_abb,
    f_bbb,
    g_aa,
    g_ab,
    g_bb,
    h_aa,
    h_ab,
    h_bb,
) = sp.symbols("f_aaa f_aab f_abb f_bbb " "g_aa g_ab g_bb h_aa h_ab h_bb")

tangential_a = {
    a: 1,
    b: 0,
    f_aa: f_aaa,
    f_ab: f_aab,
    f_bb: f_abb,
    g_a: g_aa,
    g_b: g_ab,
    h_a: h_aa,
    h_b: h_ab,
    kappa: 0,
}
tangential_b = {
    a: 0,
    b: 1,
    f_aa: f_aab,
    f_ab: f_abb,
    f_bb: f_bbb,
    g_a: g_ab,
    g_b: g_bb,
    h_a: h_ab,
    h_b: h_bb,
    kappa: 0,
}


def total_derivative(
    expression: sp.Expr, derivation: dict[sp.Symbol, sp.Expr]
) -> sp.Expr:
    return sp.together(
        sum(
            sp.diff(expression, variable) * image
            for variable, image in derivation.items()
        )
    )


V_cc_a = total_derivative(V_cc, tangential_a)
V_cc_b = total_derivative(V_cc, tangential_b)
V_cd_a = total_derivative(V_cd, tangential_a)
V_cd_b = total_derivative(V_cd, tangential_b)
V_dd_a = total_derivative(V_dd, tangential_a)
V_dd_b = total_derivative(V_dd, tangential_b)

n0, n1, n2, n3 = sp.symbols("n0 n1 n2 n3")

hessian_adapted = sp.Matrix(
    [
        [f_aa, f_ab, g_a, h_a],
        [f_ab, f_bb, g_b, h_b],
        [g_a, g_b, V_cc, V_cd],
        [h_a, h_b, V_cd, V_dd],
    ]
)

third_derivatives: dict[tuple[int, int, int], sp.Expr] = {}


def set_symmetric(indices: tuple[int, int, int], value: sp.Expr) -> None:
    for permuted in set(permutations(indices)):
        third_derivatives[permuted] = value


for indices, value in (
    ((0, 0, 0), f_aaa),
    ((0, 0, 1), f_aab),
    ((0, 1, 1), f_abb),
    ((1, 1, 1), f_bbb),
    ((0, 0, 2), g_aa),
    ((0, 1, 2), g_ab),
    ((1, 1, 2), g_bb),
    ((0, 0, 3), h_aa),
    ((0, 1, 3), h_ab),
    ((1, 1, 3), h_bb),
    ((0, 2, 2), V_cc_a),
    ((1, 2, 2), V_cc_b),
    ((0, 2, 3), V_cd_a),
    ((1, 2, 3), V_cd_b),
    ((0, 3, 3), V_dd_a),
    ((1, 3, 3), V_dd_b),
    ((2, 2, 2), n0),
    ((2, 2, 3), n1),
    ((2, 3, 3), n2),
    ((3, 3, 3), n3),
):
    set_symmetric(indices, value)

fixed_image_substitution = {
    Y: b,
    D: boundary["fixed_image_D"],
}
adapted_q = sp.Matrix([q[0], q[1], q[2], q[3] - q[1]])
normal_velocity = [
    sp.factor(sp.diff(coordinate, X).subs(fixed_image_substitution).subs(X, 0))
    for coordinate in adapted_q
]

hessian_velocity_adapted = sp.Matrix(
    4,
    4,
    lambda row, column: sum(
        third_derivatives[(row, column, index)] * normal_velocity[index]
        for index in range(4)
    ),
)

coordinate_matrix = boundary["coordinate_matrix"]
hessian_original = coordinate_matrix.inv().T * hessian_adapted * coordinate_matrix.inv()
hessian_velocity_original = (
    coordinate_matrix.inv().T * hessian_velocity_adapted * coordinate_matrix.inv()
)

jq = q.jacobian(h_variables).subs(fixed_image_substitution)
jm = m.jacobian(h_variables).subs(fixed_image_substitution)
jacobian_zero = (jm + hessian_original * jq).subs(X, 0)
jacobian_one = sp.diff(jm + hessian_original * jq, X).subs(
    X, 0
) + hessian_velocity_original * jq.subs(X, 0)

# The derivative of det(J0+X*J1) at X=0 is the cofactor pairing.
determinant_one = sum(
    jacobian_zero.cofactor(row, column) * jacobian_one[row, column]
    for row in range(4)
    for column in range(4)
)
cleared_numerator = sp.together(determinant_one).as_numer_denom()[0]
prolongation = sp.Poly(sp.expand(cleared_numerator), W)
assert prolongation.degree() == 3

coefficients = [prolongation.nth(power) for power in range(4)]

# Triangularity: W^3 sees only n3; W^2 adds n2; W^1 adds n1; W^0 adds n0.
new_normal_derivatives = (n0, n1, n2, n3)
allowed_by_power = {
    3: {n3},
    2: {n2, n3},
    1: {n1, n2, n3},
    0: {n0, n1, n2, n3},
}
for power in range(4):
    occurring = {
        variable
        for variable in new_normal_derivatives
        if sp.diff(prolongation.nth(power), variable) != 0
    }
    assert occurring == allowed_by_power[power]

expected_pivots = {
    (3, n3): 819200 * L**4,
    (2, n2): 1843200 * L**4,
    (1, n1): 1382400 * L**4,
    (0, n0): 345600 * L**4,
}
for (power, variable), expected in expected_pivots.items():
    pivot = sp.factor(sp.diff(prolongation.nth(power), variable))
    assert sp.factor(pivot - expected) == 0

# The leading equation is
#
#     256*L^3*n3 + top_remainder = 0.
#
# Thus a polynomial n3 requires a third-order divisibility condition.  Its
# first residue on L=0 has a compact exact factorization.
top_remainder = sp.cancel(
    prolongation.nth(3).subs(n3, 0) / (3200 * L)
)
assert sp.together(top_remainder).as_numer_denom()[1] == 1
assert (
    sp.factor(
        prolongation.nth(3)
        - 3200 * L * (256 * L**3 * n3 + top_remainder)
    )
    == 0
)

top_remainder_on_L = sp.together(
    top_remainder.subs(f_bb, 5 * b + f_ab**2 / f_aa)
)
residue_numerator, residue_denominator = (
    top_remainder_on_L.as_numer_denom()
)
first_residue_factor = (
    b * f_aa - 2 * f_aa * h_b + 2 * f_ab * h_a
)
second_residue_factor = (
    f_aa**3 * f_bbb
    - 5 * f_aa**3
    - 3 * f_aa**2 * f_ab * f_abb
    + 3 * f_aa * f_aab * f_ab**2
    - f_aaa * f_ab**3
)
assert sp.factor(residue_denominator - f_aa**3) == 0
assert (
    sp.factor(
        residue_numerator
        - 64 * first_residue_factor**3 * second_residue_factor
    )
    == 0
)

# The next principal symbol can be isolated without constructing the
# inhomogeneous fourth-jet prolongation.  Only the five pure-normal fourth
# derivatives contribute to this symbol.
r0, r1, r2, r3, r4 = sp.symbols("r0 r1 r2 r3 r4")
fourth_normal_derivatives = (r0, r1, r2, r3, r4)
binary_fourth_evaluation = sum(
    sp.binomial(4, d_indices)
    * fourth_normal_derivatives[d_indices]
    * normal_velocity_pair[0] ** (4 - d_indices)
    * normal_velocity_pair[1] ** d_indices
    for d_indices in range(5)
)

fourth_tensor: dict[tuple[int, int, int, int], sp.Expr] = {}
for d_indices, value in enumerate(fourth_normal_derivatives):
    indices = (2,) * (4 - d_indices) + (3,) * d_indices
    for permuted in set(permutations(indices)):
        fourth_tensor[permuted] = value

fourth_hessian_velocity_adapted = sp.zeros(4)
for row in (2, 3):
    for column in (2, 3):
        fourth_hessian_velocity_adapted[row, column] = sp.Rational(1, 2) * sum(
            fourth_tensor[(row, column, first, second)]
            * normal_velocity[first]
            * normal_velocity[second]
            for first in (2, 3)
            for second in (2, 3)
        )
fourth_hessian_velocity_original = (
    coordinate_matrix.inv().T
    * fourth_hessian_velocity_adapted
    * coordinate_matrix.inv()
)
second_symbol_matrix = fourth_hessian_velocity_original * jq.subs(X, 0)
second_symbol = sum(
    jacobian_zero.cofactor(row, column) * second_symbol_matrix[row, column]
    for row in range(4)
    for column in range(4)
)
assert sp.factor(second_symbol - L * binary_fourth_evaluation / 4) == 0

second_symbol_polynomial = sp.Poly(sp.expand(81 * second_symbol), W)
assert second_symbol_polynomial.degree() == 4
expected_second_pivots = {
    (4, r4): 1024 * L,
    (3, r3): 3072 * L,
    (2, r2): 3456 * L,
    (1, r1): 1728 * L,
    (0, r0): 324 * L,
}
for (power, variable), expected in expected_second_pivots.items():
    pivot = sp.factor(sp.diff(second_symbol_polynomial.nth(power), variable))
    assert sp.factor(pivot - expected) == 0


def main() -> None:
    print("PASS: charts 0010 and 1000 have boundary symbol L*Q(v,v)/2")
    print("PASS: the first chart-1000 prolongation has W-degree three")
    print("PASS: mixed third derivatives are represented by one symmetric tensor")
    print("PASS: W^3,W^2,W^1,W^0 solve n3,n2,n1,n0 triangularly")
    print("PASS: the four cleared pivots are nonzero multiples of L^4")
    print("PASS: polynomial V_ddd requires L^3 divisibility")
    print("PASS: the first residue on L=0 factors into two explicit branches")
    print("PASS: the second principal symbol has five nonzero L-pivots")
    print("RESULT: every formal normal order is solvable after inverting L")
    print("SCOPE: formal solvability does not imply polynomiality")


if __name__ == "__main__":
    main()
