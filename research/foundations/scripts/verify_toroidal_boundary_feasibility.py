#!/usr/bin/env python3
"""Exact regressions for the colored-fan boundary feasibility front end."""

from __future__ import annotations

import json
from dataclasses import replace
from fractions import Fraction

import sympy as sp

from boundary_package_compiler import (
    BoundaryColor,
    ColoredDivisorSpanProblem,
    IntegralVariable,
    LinearConstraint,
    PackageStatus,
    SheetColor,
    ToroidalAffineScreen,
    ToroidalBoundaryDatum,
    ToroidalFanDatum,
    TropicalFeasibilityStatus,
    TropicalRay,
    ValuationFeasibilityProblem,
    a4_three_puncture_package,
    a4_toroidal_ledger_datum,
    audit_toroidal_boundary,
    compile_boundary_package,
    colored_proportionality_witnesses,
    d5_toroidal_ledger_datum,
    davenport_toroidal_ledger_datum,
    f20_base_factor_mask_datum,
    f20_toroidal_ledger_datum,
    smith_invariant_factors,
)


def invariant(audit, name):
    return dict(audit.invariants)[name]


def obstruction_codes(audit):
    return {
        diagnostic.code
        for diagnostic in audit.diagnostics
        if diagnostic.obstruction
    }


a4 = audit_toroidal_boundary(a4_toroidal_ledger_datum())
assert a4.status is TropicalFeasibilityStatus.FEASIBLE
assert invariant(a4, "valuation_identities")[
    "pure_target_log_keller_balance"
] == {
    "actual_orders": (0, 0, 0),
    "expected_orders": (0, 0, 0),
    "balanced": True,
}
assert all(
    cone["smooth"] for cone in invariant(a4, "fan")["maximal_cones"]
)
assert not obstruction_codes(a4)

compiled_a4 = compile_boundary_package(
    replace(
        a4_three_puncture_package(),
        name="A4 with toroidal front end",
        toroidal_boundary=a4_toroidal_ledger_datum(),
    )
)
assert compiled_a4.status is PackageStatus.UNKNOWN
assert dict(compiled_a4.invariants)["toroidal_boundary"]["status"] == (
    "feasible"
)
print("PASS: the A4 W/K/L ledger is one exact valuation-matrix identity")


d5 = audit_toroidal_boundary(d5_toroidal_ledger_datum())
assert d5.status is TropicalFeasibilityStatus.FEASIBLE
d5_fan = invariant(d5, "fan")
assert tuple(
    profile["smith_diagonal"] for profile in d5_fan["maximal_cones"]
) == ((1, 1), (1, 1), (1, 1))

# A quadratic a^2-lambda*u has monomial order min(2*p,q) on ray (p,q).
# This derives the exceptional rows rather than merely rereading them.
e11_ray = d5_fan["rays"]["first_blowup"]["vector"]
e12_ray = d5_fan["rays"]["parabolic_scale"]["vector"]
quadratic_orders = tuple(
    min(2 * ray[0], ray[1]) for ray in (e11_ray, e12_ray)
)
assert quadratic_orders == (1, 2)
derived_exceptional_rows = tuple(
    (2 * order, 5 * order, order, order, order)
    for order in quadratic_orders
)
assert derived_exceptional_rows == (
    (2, 5, 1, 1, 1),
    (4, 10, 2, 2, 2),
)

# In Q(sqrt(5)), store x+y*sqrt(5) as (x,y).  The three residues after
# the (1,2) extraction are 1/4, (3-sqrt(5))/2, (3+sqrt(5))/2.
parabolic_residues = {
    (Fraction(1, 4), Fraction(0)),
    (Fraction(3, 2), Fraction(-1, 2)),
    (Fraction(3, 2), Fraction(1, 2)),
}
assert len(parabolic_residues) == 3
d5_problem = invariant(d5, "feasibility_problems")[
    "branch_supported_log_balance_m_le_4"
]
assert d5_problem["search_size"] == 5 * 9**3
assert d5_problem["model_count"] == 4
assert d5_problem["minimal_models"] == (
    {"m": 1, "s_C": 1, "s_plus": 1, "s_minus": 1},
)
assert invariant(d5, "valuation_matrix")["matrix"][-2:] == (
    (2, 5, 1, 1, 1),
    (4, 10, 2, 2, 2),
)
print("PASS: the D5 Newton fan adds both exceptional valuation rows")
print("PASS: the D5 bounded system has one primitive Pareto model")


davenport = audit_toroidal_boundary(davenport_toroidal_ledger_datum())
assert davenport.status is TropicalFeasibilityStatus.FEASIBLE
cox_block = invariant(davenport, "unimodular_blocks")[
    "cox_lattice_completion"
]
assert cox_block["determinant"] == 1
assert cox_block["unimodular"] is True
assert "new divisor L(Y)" in invariant(davenport, "nonlinear_residue")[0]
print("PASS: the Davenport three-column ledger is integrally unimodular")
print("PASS: the new L(Y) divisor remains typed as nonlinear residue")


f20 = audit_toroidal_boundary(f20_toroidal_ledger_datum())
assert f20.status is TropicalFeasibilityStatus.FEASIBLE
assert all(
    cone["smooth"] for cone in invariant(f20, "fan")["maximal_cones"]
)
f20_matrix = invariant(f20, "valuation_matrix")
expected_f20_divisor_rows = (
    (1, 0, 0, 0),
    (4, 0, 0, 3),
    (0, 1, 0, 1),
    (0, 1, 0, 1),
    (0, 1, 0, 0),
    (0, 1, 0, 0),
    (0, 1, 0, 0),
    (0, 2, 0, 1),
    (0, 2, 0, 1),
    (0, 2, 0, 1),
    (0, 2, 0, 1),
    (0, 2, 0, 0),
    (0, 0, 2, 1),
    (0, 0, 2, 1),
    (0, 0, 1, 0),
    (0, 0, 10, 4),
    (0, 0, 20, 8),
    (0, 0, 5, 2),
    (0, 0, 10, 4),
    (0, 0, 10, 4),
    (0, 0, 10, 4),
    (0, 0, 10, 4),
    (0, 0, 10, 4),
    (0, 0, 10, 4),
    (0, 0, 10, 4),
    (4, 4, 4, 7),
    (1, 1, 1, 0),
    (4, 4, 4, 7),
    (1, 1, 1, 0),
    (2, 2, 1, 3),
    (2, 2, 1, 3),
    (2, 2, 1, 3),
    (2, 2, 1, 3),
    (2, 2, 1, 0),
    (2, 2, 1, 3),
    (2, 2, 1, 3),
    (2, 2, 1, 3),
    (2, 2, 1, 3),
    (2, 2, 1, 0),
) + (
    (0, 2, 2, 2),
    (0, 1, 1, 1),
    (0, 2, 2, 1),
    (0, 2, 2, 2),
    (0, 2, 2, 2),
    (0, 2, 2, 2),
    (0, 2, 2, 1),
    (0, 2, 2, 1),
) * 3
assert f20_matrix["matrix"] == tuple(
    row + row[:3] + (0,) for row in expected_f20_divisor_rows
)
assert f20_matrix["columns"][-4:] == (
    "mask_d",
    "mask_q",
    "mask_r",
    "q_selector_w_minus_1",
)

# Lecacheux's primary formula has t^2*d in the quartic coefficient.  The
# missing square in JLY Theorem 2.3.6 destroys the F20 constraint away from
# t=0,1, so the exact discriminant is also a transcription guard.
s, t, X = sp.symbols("s t X")
d = s**2 + 4
q = 4 * s**2 * t**2 + 4 * s**2 * t + 8 * s * t + 6 * s - 8 * t - 5
r = (
    16 * s**2 * t**3
    + 4 * s**2 * t**2
    - 76 * s * t
    - 16 * s
    + 64 * t**3
    + 16 * t**2
    - 164 * t
    - 199
)
P = sp.expand(
    X**5
    + (t**2 * d - 2 * s - sp.Rational(17, 4)) * X**4
    + (3 * t * d + d + sp.Rational(13, 2) * s + 1) * X**3
    - (t * d + sp.Rational(11, 2) * s - 8) * X**2
    + (s - 6) * X
    + 1
)
P_X = sp.diff(P, X)
P_t = sp.diff(P, t)
assert sp.factor(sp.discriminant(P, X)) == d**3 * q**2 * r**2 / 256
assert sp.factor(P_t - X**2 * d * (2 * X**2 * t + 3 * X - 1)) == 0

# Over d=0 there is one simple color and one tame index-four color.
d_profile = (X - sp.Rational(1, 4)) * (X - 1 - s / 2) ** 4
assert sp.rem(P - d_profile, d, s) == 0

coefficient_field = sp.QQ.frac_field(s)


def residue_in_t(expression, divisor):
    divisor_poly = sp.Poly(divisor, t, domain=coefficient_field)
    return sp.Poly(expression, t, domain=coefficient_field).rem(
        divisor_poly
    ).as_expr()


def vanishes_on_divisor(expression, divisor):
    return all(
        residue_in_t(coefficient, divisor) == 0
        for coefficient in sp.Poly(expression, X).all_coeffs()
    )


# Along q=0 the polynomial has one generic double root.  Both P_t and P_X
# vanish there, so this is an unramified transverse crossing, not an e=2
# inertia color.  The nonzero quadratic tangent discriminant separates its
# two geometric normalization branches, and P_X has order one on each.
a_q = (2 * s**2 * t + 2 * s**2 + 3 * s - 4) / (2 * (s - 1))
assert residue_in_t(P.subs(X, a_q), q) == 0
assert residue_in_t(P_X.subs(X, a_q), q) == 0
assert residue_in_t(P_t.subs(X, a_q), q) == 0
q_P_xx = residue_in_t(sp.diff(P, X, 2).subs(X, a_q), q)
q_P_xt = residue_in_t(sp.diff(P, X, t).subs(X, a_q), q)
q_P_tt = residue_in_t(sp.diff(P, t, 2).subs(X, a_q), q)
assert q_P_xx != 0
assert residue_in_t(q_P_xt**2 - q_P_xx * q_P_tt, q) != 0

subresultants = sp.subresultants(P, P_X, X)
assert tuple(sp.degree(item, X) for item in subresultants) == (5, 4, 3, 2, 1, 0)
assert not vanishes_on_divisor(subresultants[-2], q)
assert vanishes_on_divisor(subresultants[-1], q)

# Along r=0 the generic gcd has degree two.  The rational sample below is
# off the other factors and shows that P_t is a unit at both repeated roots;
# hence the two geometric colors have e=2 and derivative order one.
assert not vanishes_on_divisor(subresultants[-3], r)
assert vanishes_on_divisor(subresultants[-2], r)
assert vanishes_on_divisor(subresultants[-1], r)
r_sample = sp.Poly(P.subs({s: 2, t: sp.Rational(7, 4)}), X)
r_repeated = sp.Poly(X**2 + 8 * X - 2, X)
assert r_sample == sp.Poly((4 * X + 1) * r_repeated.as_expr() ** 2 / 4, X)
assert sp.gcd(
    r_repeated,
    sp.Poly(P_t.subs({s: 2, t: sp.Rational(7, 4)}), X),
).degree() == 0

# Each generic target color has total root-cover degree five.  The sums of
# derivative orders recover the exponents (3,2,2) of the three discriminant
# factors without another family-specific ledger.
matrix = tuple(row[:4] for row in f20_matrix["matrix"])
assert sum(row[0] for row in matrix[0:2]) == 5
assert sum(row[1] for row in matrix[2:7]) == 5
assert sum(row[2] for row in matrix[12:15]) == 5
assert (
    sum(row[3] for row in matrix[0:2]),
    sum(row[3] for row in matrix[2:7]),
    sum(row[3] for row in matrix[12:15]),
) == (3, 2, 2)
print("PASS: the corrected Lecacheux discriminant factors as d^3*q^2*r^2/256")
print("PASS: the F20 root cover has exact (4,1), crossing, and (2,2,1) colors")

# The base-boundary incidence atlas is finite.  The q-curve has one ordinary
# node.  After the tangent change w=(s-11)-50*(t+1/2), the r-curve is
# quadratic in w with discriminant v^5 times a unit, hence has an A4
# ramphoid cusp.  These assertions identify the singular blowup centers
# without pretending that their root-cover exceptional colors are known.
u, v, w = sp.symbols("u v w")
q_singular_basis = sp.groebner(
    (q, sp.diff(q, s), sp.diff(q, t)), t, s, domain=sp.QQ
)
assert sp.expand(q_singular_basis.polys[0].as_expr() - (t + sp.Rational(1, 2))) == 0
assert sp.expand(q_singular_basis.polys[1].as_expr() - (s - 1)) == 0
q_node = sp.expand(q.subs({s: 1 + u, t: -sp.Rational(1, 2) + v}))
assert q_node == 4 * u**2 * v**2 - u**2 + 8 * u * v**2 + 8 * u * v + 4 * v**2
assert sp.discriminant(-u**2 + 8 * u * v + 4 * v**2, u) == 80 * v**2

# Blow up the node in the chart u=epsilon, v=z*epsilon.  The special root
# polynomial is (X+1)*(X^2-3X+1)^2.  For either geometric root a of the
# quadratic, X=a+epsilon*Y has a separable quadratic slope polynomial.
# Hence the exceptional divisor has four unramified derivative-order-one
# colors and the simple X=-1 color.
epsilon, z, Y, a = sp.symbols("epsilon z Y a")
node_root_factor = a**2 - 3 * a + 1
node_deformation = sp.expand(
    P.subs(
        {
            s: 1 + epsilon,
            t: -sp.Rational(1, 2) + z * epsilon,
            X: a + epsilon * Y,
        }
    )
)


def residue_at_node_root(expression):
    return sp.rem(expression, node_root_factor, a)


assert residue_at_node_root(sp.Poly(node_deformation, epsilon).nth(0)) == 0
assert residue_at_node_root(sp.Poly(node_deformation, epsilon).nth(1)) == 0
node_slope_polynomial = residue_at_node_root(
    sp.Poly(node_deformation, epsilon).nth(2)
)
assert sp.Poly(node_slope_polynomial, Y).degree() == 2
assert sp.expand(
    sp.Poly(node_slope_polynomial, Y).LC() - 5 * (a + 1)
) == 0
node_slope_discriminant = residue_at_node_root(
    sp.discriminant(node_slope_polynomial, Y)
)
assert sp.expand(
    node_slope_discriminant
    + sp.Rational(25, 4)
    * (
        796 * a * z**2
        - 188 * a * z
        + 11 * a
        - 304 * z**2
        + 72 * z
        - 4
    )
) == 0
assert node_slope_discriminant != 0

node_derivative = sp.expand(
    P_X.subs(
        {
            s: 1 + epsilon,
            t: -sp.Rational(1, 2) + z * epsilon,
            X: a + epsilon * Y,
        }
    )
)
node_derivative_lead = residue_at_node_root(
    sp.Poly(node_derivative, epsilon).nth(1)
)
assert sp.expand(node_derivative_lead - sp.diff(node_slope_polynomial, Y)) == 0
assert sp.expand(
    P.subs({s: 1, t: -sp.Rational(1, 2)})
    - (X + 1) * (X**2 - 3 * X + 1) ** 2
) == 0
assert P_X.subs({s: 1, t: -sp.Rational(1, 2), X: -1}) != 0
node_q_pullback = sp.expand(
    q.subs(
        {
            s: 1 + epsilon,
            t: -sp.Rational(1, 2) + z * epsilon,
        }
    )
)
assert sp.Poly(node_q_pullback, epsilon).nth(0) == 0
assert sp.Poly(node_q_pullback, epsilon).nth(1) == 0
assert sp.Poly(node_q_pullback, epsilon).nth(2) == 4 * z**2 + 8 * z - 1
assert sum(row[1] for row in matrix[7:12]) == 10
assert sum(row[3] for row in matrix[7:12]) == 4
print("PASS: the q-node blowup adds five exact unramified exceptional colors")

r_singular_basis = sp.groebner(
    (r, sp.diff(r, s), sp.diff(r, t)), t, s, domain=sp.QQ
)
assert sp.expand(r_singular_basis.polys[1].as_expr() - (s - 11) ** 4) == 0
assert sp.expand(
    r_singular_basis.polys[0].as_expr().subs(s, 11)
    - (t + sp.Rational(1, 2))
) == 0
r_cusp = sp.expand(
    r.subs(
        {
            s: 11 + w + 50 * v,
            t: -sp.Rational(1, 2) + v,
        }
    )
)
r_cusp_in_w = sp.Poly(r_cusp, w)
r_cusp_A, r_cusp_B, r_cusp_C = r_cusp_in_w.all_coeffs()
assert r_cusp_A.subs(v, 0) == -1
assert sp.factor(r_cusp_B**2 - 4 * r_cusp_A * r_cusp_C) == (
    -2048 * v**5 * (2 * v - 5)
)

# The standard embedded resolution of the A4 cusp has four exceptional
# valuations.  The following substitutions are their generic charts in the
# tangent coordinates v=t+1/2 and w=(s-11)-50*v.  Newton residuals determine
# the root-cover colors without a Galois-group oracle.
tau, x = sp.symbols("tau x")


def lowest_coefficient(expression, variable):
    polynomial = sp.Poly(sp.expand(expression), variable)
    orders = tuple(
        order
        for order in range(polynomial.degree() + 1)
        if polynomial.nth(order) != 0
    )
    order = min(orders)
    return order, sp.factor(polynomial.nth(order))


def scaled_root_residual(source_s, source_t, base_power, root_power):
    shifted = sp.expand(
        P.subs({s: source_s, t: source_t, X: x - 1})
    )
    scaled = shifted.subs(
        {epsilon: tau**base_power, x: tau**root_power * Y}
    )
    order, residual = lowest_coefficient(scaled, tau)
    scaled_derivative = sp.expand(
        P_X.subs(
            {
                s: source_s,
                t: source_t,
                X: x - 1,
            }
        ).subs(
            {epsilon: tau**base_power, x: tau**root_power * Y}
        )
    )
    derivative_order, derivative_residual = lowest_coefficient(
        scaled_derivative, tau
    )
    assert derivative_order == order - root_power
    assert sp.expand(
        derivative_residual - sp.diff(residual, Y)
    ) == 0
    return order, residual


cusp_t = -sp.Rational(1, 2) + epsilon
cusp_E1_s = 11 + epsilon * (z + 50)
E1_r_order, E1_r_lead = lowest_coefficient(
    r.subs({s: cusp_E1_s, t: cusp_t}), epsilon
)
assert (E1_r_order, E1_r_lead) == (2, -z**2)
E1_root_order, E1_residual = scaled_root_residual(
    cusp_E1_s, cusp_t, 5, 1
)
assert E1_root_order == 5
assert sp.expand(E1_residual - (Y**5 + sp.Rational(25, 2) * z)) == 0
assert sp.factor(sp.discriminant(E1_residual, Y)) == (
    sp.Rational(1220703125, 16) * z**4
)

cusp_E2_s = 11 + 50 * epsilon + epsilon**2 * z
E2_r_order, E2_r_lead = lowest_coefficient(
    r.subs({s: cusp_E2_s, t: cusp_t}), epsilon
)
assert (E2_r_order, E2_r_lead) == (4, -(z - 180) ** 2)
E2_root_order, E2_residual = scaled_root_residual(
    cusp_E2_s, cusp_t, 5, 2
)
assert E2_root_order == 10
assert sp.expand(
    E2_residual - (Y**5 + sp.Rational(25, 2) * (z - 180))
) == 0
assert sp.factor(sp.discriminant(E2_residual, Y)) == (
    sp.Rational(1220703125, 16) * (z - 180) ** 4
)

cusp_E3_s = (
    11 + 50 * epsilon + 180 * epsilon**2 + epsilon**3 * z
)
E3_r_order, E3_r_lead = lowest_coefficient(
    r.subs({s: cusp_E3_s, t: cusp_t}), epsilon
)
assert (E3_r_order, E3_r_lead) == (5, 2560)
E3_unramified_order, E3_unramified_residual = scaled_root_residual(
    cusp_E3_s, cusp_t, 1, 1
)
assert E3_unramified_order == 3
assert sp.expand(
    E3_unramified_residual - sp.Rational(25, 2) * (40 * Y + z - 696)
) == 0
E3_ramified_order, E3_ramified_residual = scaled_root_residual(
    cusp_E3_s, cusp_t, 2, 1
)
assert E3_ramified_order == 5
assert sp.expand(
    E3_ramified_residual - Y * (Y**4 - 50 * Y**2 + 500)
) == 0
E3_quartic = sp.cancel(E3_ramified_residual / Y)
assert sp.discriminant(E3_quartic, Y) == 2000000000

cusp_E4_s = (
    11
    + 50 * epsilon**2 * z
    + 180 * epsilon**4 * z**2
    + epsilon**5 * z**2
)
cusp_E4_t = -sp.Rational(1, 2) + epsilon**2 * z
E4_r_order, E4_r_lead = lowest_coefficient(
    r.subs({s: cusp_E4_s, t: cusp_E4_t}), epsilon
)
assert (E4_r_order, E4_r_lead) == (
    10,
    z**4 * (2560 * z - 1),
)
E4_root_order, E4_residual = scaled_root_residual(
    cusp_E4_s, cusp_E4_t, 1, 1
)
assert E4_root_order == 5
assert sp.expand(
    E4_residual
    - (
        Y**5
        - 50 * z * Y**3
        + 500 * z**2 * Y
        + sp.Rational(25, 2) * z**2
    )
) == 0
assert sp.factor(sp.discriminant(E4_residual, Y)) == (
    sp.Rational(1220703125, 16) * z**8 * (2560 * z - 1) ** 2
)

# Source valuations multiply base pullback orders by the ramification index.
# Their derivative sums recover twice the four base r-orders 2,4,5,10.
assert matrix[15] == (0, 0, 10, 4)
assert matrix[16] == (0, 0, 20, 8)
assert matrix[17:20] == (
    (0, 0, 5, 2),
    (0, 0, 10, 4),
    (0, 0, 10, 4),
)
assert matrix[20:25] == ((0, 0, 10, 4),) * 5
assert (
    matrix[15][3],
    matrix[16][3],
    sum(row[3] for row in matrix[17:20]),
    sum(row[3] for row in matrix[20:25]),
) == (4, 8, 10, 20)
print("PASS: four A4-cusp blowups add ten exact F20 exceptional colors")

# Pairwise elimination gives two conjugate triple points where d and q are
# tangent and r is transverse, one additional transverse q-r point, and
# three q-r tangencies.  The factor t=0 in the raw q-r resultant is only a
# leading-coefficient intersection at infinity and is absent from the affine
# Groebner elimination polynomial.
h2 = 16 * t**2 + 24 * t + 13
h3 = 8 * t**3 + 16 * t**2 + 2 * t - 7
assert sp.expand(sp.resultant(d, q, s) - h2**2) == 0
assert sp.expand(sp.resultant(d, r, s) - 3125 * h2) == 0
assert sp.expand(
    sp.resultant(q, r, s) - 64 * t * (t - 2) * h2 * h3**2
) == 0

d_r_basis = sp.groebner((d, r), s, t, domain=sp.QQ)
assert sp.expand(d_r_basis.polys[0].as_expr() - (s - 4 * t - 3)) == 0
assert sp.expand(d_r_basis.polys[1].as_expr() - h2 / 16) == 0
assert d_r_basis.reduce(q)[1] == 0

q_r_basis = sp.groebner((q, r), s, t, domain=sp.QQ)
assert sp.expand(
    q_r_basis.polys[-1].as_expr() - (t - 2) * h2 * h3**2 / 1024
) == 0
assert q.subs({s: sp.Rational(7, 12), t: 2}) == 0
assert r.subs({s: sp.Rational(7, 12), t: 2}) == 0
q_r_jacobian = sp.det(
    sp.Matrix(
        (
            (sp.diff(q, s), sp.diff(q, t)),
            (sp.diff(r, s), sp.diff(r, t)),
        )
    )
)
assert q_r_jacobian.subs({s: sp.Rational(7, 12), t: 2}) != 0
print("PASS: the F20 base incidence has one node, one A4 cusp, and five finite tangencies")

# Resolve either conjugate triple point over Q(i).  On the first exceptional
# divisor the fourfold root has a square first residual; its second residual
# is separable and gives one index-four color with derivative order seven.
# The infinitely-near divisor has a separable quartic cluster and one simple
# sheet.  Complex conjugation supplies the identical integer rows at the
# second center.
imaginary_unit = sp.I
triple_s0 = 2 * imaginary_unit
triple_t0 = -sp.Rational(3, 4) + imaginary_unit / 2
triple_root = 1 + imaginary_unit
triple_simple_root = sp.Rational(1, 4)
assert sp.expand(
    P.subs({s: triple_s0, t: triple_t0})
    - (X - triple_simple_root) * (X - triple_root) ** 4
) == 0

triple_E1_s = triple_s0 + epsilon
triple_E1_t = triple_t0 + z * epsilon
assert lowest_coefficient(
    d.subs({s: triple_E1_s, t: triple_E1_t}), epsilon
)[0] == 1
assert lowest_coefficient(
    q.subs({s: triple_E1_s, t: triple_E1_t}), epsilon
)[0] == 1
assert lowest_coefficient(
    r.subs({s: triple_E1_s, t: triple_E1_t}), epsilon
)[0] == 1
triple_E1_shifted = sp.expand(
    P.subs(
        {
            s: triple_E1_s,
            t: triple_E1_t,
            X: triple_root + x,
        }
    )
)
triple_E1_first = triple_E1_shifted.subs(
    {epsilon: tau**2, x: tau * Y}
)
triple_E1_first_order, triple_E1_first_residual = lowest_coefficient(
    triple_E1_first, tau
)
assert triple_E1_first_order == 4
assert sp.expand(
    triple_E1_first_residual
    - (sp.Rational(3, 4) + imaginary_unit)
    * (Y**2 + imaginary_unit) ** 2
) == 0

b = sp.symbols("b")
triple_first_root_relation = b**2 + imaginary_unit


def lowest_coefficient_mod(expression, variable, modulus, modulus_variable):
    polynomial = sp.Poly(sp.expand(expression), variable)
    for order in range(polynomial.degree() + 1):
        residual = sp.rem(
            polynomial.nth(order), modulus, modulus_variable
        )
        if residual != 0:
            return order, sp.factor(residual)
    raise AssertionError("the reduced expansion is zero")


triple_E1_second = sp.expand(
    P.subs(
        {
            s: triple_s0 + tau**4,
            t: triple_t0 + z * tau**4,
            X: triple_root + tau**2 * (b + tau * Y),
        }
    )
)
triple_E1_second_order, triple_E1_second_residual = lowest_coefficient_mod(
    triple_E1_second,
    tau,
    triple_first_root_relation,
    b,
)
assert triple_E1_second_order == 10
expected_triple_E1_second = (16 + 8 * imaginary_unit) * (
    (sp.Rational(1, 8) - imaginary_unit / 4) * Y**2
    + b * (z + sp.Rational(1, 20) + imaginary_unit / 10)
)
assert sp.expand(
    triple_E1_second_residual - expected_triple_E1_second
) == 0
assert sp.rem(
    sp.discriminant(triple_E1_second_residual, Y),
    triple_first_root_relation,
    b,
) != 0
triple_E1_second_derivative = sp.expand(
    P_X.subs(
        {
            s: triple_s0 + tau**4,
            t: triple_t0 + z * tau**4,
            X: triple_root + tau**2 * (b + tau * Y),
        }
    )
)
triple_E1_derivative_order, triple_E1_derivative_residual = (
    lowest_coefficient_mod(
        triple_E1_second_derivative,
        tau,
        triple_first_root_relation,
        b,
    )
)
assert triple_E1_derivative_order == 7
assert sp.expand(
    triple_E1_derivative_residual
    - sp.diff(triple_E1_second_residual, Y)
) == 0
assert P_X.subs(
    {s: triple_s0, t: triple_t0, X: triple_simple_root}
) != 0

# The second blowup uses t-t0=epsilon and s-s0=z*epsilon^2.
triple_E2_s = triple_s0 + z * epsilon**2
triple_E2_t = triple_t0 + epsilon
assert lowest_coefficient(
    d.subs({s: triple_E2_s, t: triple_E2_t}), epsilon
)[0] == 2
assert lowest_coefficient(
    q.subs({s: triple_E2_s, t: triple_E2_t}), epsilon
)[0] == 2
assert lowest_coefficient(
    r.subs({s: triple_E2_s, t: triple_E2_t}), epsilon
)[0] == 1
triple_E2_shifted = sp.expand(
    P.subs(
        {
            s: triple_E2_s,
            t: triple_E2_t,
            X: triple_root + epsilon * Y,
        }
    )
)
triple_E2_order, triple_E2_residual = lowest_coefficient(
    triple_E2_shifted, epsilon
)
assert triple_E2_order == 4
expected_triple_E2 = (
    (sp.Rational(3, 4) + imaginary_unit) * Y**4
    + z * (-2 + 3 * imaginary_unit / 2) * Y**2
    + z * (16 + 8 * imaginary_unit) * Y
    + z
    * (-sp.Rational(3, 4) - imaginary_unit)
    * (z + sp.Rational(256, 25) + 192 * imaginary_unit / 25)
)
assert sp.expand(triple_E2_residual - expected_triple_E2) == 0
assert sp.factor(
    sp.discriminant(triple_E2_residual, Y),
    extension=imaginary_unit,
) == (
    z**3
    * (-199424 + 15168 * imaginary_unit)
    * (z - sp.Rational(64, 25) - 48 * imaginary_unit / 25) ** 2
)
triple_E2_derivative = sp.expand(
    P_X.subs(
        {
            s: triple_E2_s,
            t: triple_E2_t,
            X: triple_root + epsilon * Y,
        }
    )
)
triple_E2_derivative_order, triple_E2_derivative_residual = (
    lowest_coefficient(triple_E2_derivative, epsilon)
)
assert triple_E2_derivative_order == 3
assert sp.expand(
    triple_E2_derivative_residual - sp.diff(triple_E2_residual, Y)
) == 0

assert matrix[25:29] == (
    (4, 4, 4, 7),
    (1, 1, 1, 0),
    (4, 4, 4, 7),
    (1, 1, 1, 0),
)
assert matrix[29:34] == ((2, 2, 1, 3),) * 4 + ((2, 2, 1, 0),)
assert matrix[34:39] == ((2, 2, 1, 3),) * 4 + ((2, 2, 1, 0),)
assert sum(row[3] for row in matrix[25:27]) == 7
assert sum(row[3] for row in matrix[29:34]) == 12
print("PASS: both triple tangencies add fourteen exact F20 exceptional colors")

# The remaining three tangencies are one cubic orbit.  Work over the residue
# field Q[alpha]/(8*alpha^3+16*alpha^2+2*alpha-7), where
# s0=4*alpha^2-5.  The center polynomial has a triple root A and a double
# root B.  This one calculation therefore certifies all three geometric
# centers without choosing radicals for the cubic.
alpha = sp.symbols("alpha")
qr_center_polynomial = 8 * alpha**3 + 16 * alpha**2 + 2 * alpha - 7
qr_s0 = 4 * alpha**2 - 5
qr_A = -2 - 2 * alpha
qr_B = 8 - 12 * alpha**2 - 8 * alpha


def reduce_at_qr_center(expression):
    return sp.rem(sp.expand(expression), qr_center_polynomial, alpha)


qr_center_factorization_error = sp.Poly(
    P.subs({s: qr_s0, t: alpha})
    - (X - qr_A) ** 3 * (X - qr_B) ** 2,
    X,
)
assert all(
    reduce_at_qr_center(coefficient) == 0
    for coefficient in qr_center_factorization_error.all_coeffs()
)
assert reduce_at_qr_center(d.subs({s: qr_s0, t: alpha})) != 0
assert reduce_at_qr_center(q.subs({s: qr_s0, t: alpha})) == 0
assert reduce_at_qr_center(r.subs({s: qr_s0, t: alpha})) == 0

# On the first blowup use t-alpha=epsilon and
# s-s0=z*epsilon.  Both boundary equations have order one.  At A the first
# Newton residual is Y*(c*Y^2+l): its two nonzero roots form an index-two
# color, while the zero root refines to one unramified color.  The B-cluster
# gives a second index-two color.
qr_E1_s = qr_s0 + z * epsilon
qr_E1_t = alpha + epsilon
qr_E1_d_order, _qr_E1_d_lead = lowest_coefficient_mod(
    d.subs({s: qr_E1_s, t: qr_E1_t}),
    epsilon,
    qr_center_polynomial,
    alpha,
)
qr_E1_q_order, qr_E1_q_lead = lowest_coefficient_mod(
    q.subs({s: qr_E1_s, t: qr_E1_t}),
    epsilon,
    qr_center_polynomial,
    alpha,
)
qr_E1_r_order, qr_E1_r_lead = lowest_coefficient_mod(
    r.subs({s: qr_E1_s, t: qr_E1_t}),
    epsilon,
    qr_center_polynomial,
    alpha,
)
assert (qr_E1_d_order, qr_E1_q_order, qr_E1_r_order) == (0, 1, 1)

qr_A_cubic_coefficient = 48 * alpha**2 + 42 * alpha - 26
qr_E1_A_linear_coefficient = (
    34 * alpha**2 * z
    - 200 * alpha**2
    + 90 * alpha * z
    - 504 * alpha
    + 59 * z
    - 317
)
qr_E1_A_deformation = P.subs(
    {
        s: qr_E1_s,
        t: qr_E1_t,
        X: qr_A + tau * Y,
    }
).subs(epsilon, tau**2)
qr_E1_A_order, qr_E1_A_residual = lowest_coefficient_mod(
    qr_E1_A_deformation,
    tau,
    qr_center_polynomial,
    alpha,
)
assert qr_E1_A_order == 3
assert sp.expand(
    qr_E1_A_residual
    - Y
    * (
        qr_A_cubic_coefficient * Y**2
        + qr_E1_A_linear_coefficient
    )
) == 0
qr_E1_A_discriminant = reduce_at_qr_center(
    sp.discriminant(qr_E1_A_residual, Y)
)
assert qr_E1_A_discriminant != 0
qr_discriminant_A_unit = (
    1072 * alpha**2
    + sp.Rational(43561, 16) * alpha
    + sp.Rational(111205, 64)
)
assert reduce_at_qr_center(
    qr_E1_A_discriminant
    - qr_discriminant_A_unit * qr_E1_q_lead**2 * qr_E1_r_lead
) == 0

qr_E1_A_derivative = P_X.subs(
    {
        s: qr_E1_s,
        t: qr_E1_t,
        X: qr_A + tau * Y,
    }
).subs(epsilon, tau**2)
qr_E1_A_derivative_order, qr_E1_A_derivative_residual = (
    lowest_coefficient_mod(
        qr_E1_A_derivative,
        tau,
        qr_center_polynomial,
        alpha,
    )
)
assert qr_E1_A_derivative_order == 2
assert sp.expand(
    qr_E1_A_derivative_residual - sp.diff(qr_E1_A_residual, Y)
) == 0

qr_E1_A_fine = P.subs(
    {s: qr_E1_s, t: qr_E1_t, X: qr_A + epsilon * Y}
)
qr_E1_A_fine_order, qr_E1_A_fine_residual = lowest_coefficient_mod(
    qr_E1_A_fine,
    epsilon,
    qr_center_polynomial,
    alpha,
)
qr_E1_A_fine_constant = (
    -17 * alpha**2 * z**2
    + 516 * alpha**2
    - 43 * alpha * z**2
    + 1304 * alpha
    - sp.Rational(109, 4) * z**2
    + 828
)
assert qr_E1_A_fine_order == 2
assert sp.expand(
    qr_E1_A_fine_residual
    - (
        qr_E1_A_linear_coefficient * Y
        + qr_E1_A_fine_constant
    )
) == 0
qr_E1_A_fine_derivative = P_X.subs(
    {s: qr_E1_s, t: qr_E1_t, X: qr_A + epsilon * Y}
)
qr_E1_A_fine_derivative_order, qr_E1_A_fine_derivative_residual = (
    lowest_coefficient_mod(
        qr_E1_A_fine_derivative,
        epsilon,
        qr_center_polynomial,
        alpha,
    )
)
assert qr_E1_A_fine_derivative_order == 1
assert sp.expand(
    qr_E1_A_fine_derivative_residual
    - sp.diff(qr_E1_A_fine_residual, Y)
) == 0

qr_B_quadratic_coefficient = 72 * alpha**2 + 36 * alpha - 110
qr_E1_B_constant = (
    2364 * alpha**2 * z
    - 6288 * alpha**2
    + 1786 * alpha * z
    - 4404 * alpha
    - 1653 * z
    + 4210
)
qr_E1_B_deformation = P.subs(
    {
        s: qr_E1_s,
        t: qr_E1_t,
        X: qr_B + tau * Y,
    }
).subs(epsilon, tau**2)
qr_E1_B_order, qr_E1_B_residual = lowest_coefficient_mod(
    qr_E1_B_deformation,
    tau,
    qr_center_polynomial,
    alpha,
)
assert qr_E1_B_order == 2
assert sp.expand(
    qr_E1_B_residual
    + (qr_B_quadratic_coefficient * Y**2 + qr_E1_B_constant) / 2
) == 0
qr_E1_B_discriminant = reduce_at_qr_center(
    sp.discriminant(qr_E1_B_residual, Y)
)
qr_discriminant_B_unit = (
    597 * alpha**2 + 581 * alpha - sp.Rational(975, 2)
)
assert qr_E1_B_discriminant != 0
assert reduce_at_qr_center(
    qr_E1_B_discriminant
    - qr_discriminant_B_unit * qr_E1_r_lead
) == 0
qr_E1_B_derivative = P_X.subs(
    {
        s: qr_E1_s,
        t: qr_E1_t,
        X: qr_B + tau * Y,
    }
).subs(epsilon, tau**2)
qr_E1_B_derivative_order, qr_E1_B_derivative_residual = (
    lowest_coefficient_mod(
        qr_E1_B_derivative,
        tau,
        qr_center_polynomial,
        alpha,
    )
)
assert qr_E1_B_derivative_order == 1
assert sp.expand(
    qr_E1_B_derivative_residual - sp.diff(qr_E1_B_residual, Y)
) == 0

# The common tangent slope is m=-8*alpha^2-4*alpha+10.  In the second
# blowup chart s=s0+m*epsilon+z*epsilon^2, both q and r have order two.
# The A and B residuals are generically separable cubics and quadratics, so
# all five colors are unramified.  Their discriminants factor as q^2*r and r
# respectively, exactly matching the two local root clusters.
qr_tangent_slope = -8 * alpha**2 - 4 * alpha + 10
qr_q_s = reduce_at_qr_center(sp.diff(q, s).subs({s: qr_s0, t: alpha}))
qr_q_t = reduce_at_qr_center(sp.diff(q, t).subs({s: qr_s0, t: alpha}))
qr_r_s = reduce_at_qr_center(sp.diff(r, s).subs({s: qr_s0, t: alpha}))
qr_r_t = reduce_at_qr_center(sp.diff(r, t).subs({s: qr_s0, t: alpha}))
assert reduce_at_qr_center(qr_q_s * qr_tangent_slope + qr_q_t) == 0
assert reduce_at_qr_center(qr_r_s * qr_tangent_slope + qr_r_t) == 0

qr_E2_s = qr_s0 + qr_tangent_slope * epsilon + z * epsilon**2
qr_E2_t = alpha + epsilon
qr_E2_d_order, _qr_E2_d_lead = lowest_coefficient_mod(
    d.subs({s: qr_E2_s, t: qr_E2_t}),
    epsilon,
    qr_center_polynomial,
    alpha,
)
qr_E2_q_order, qr_E2_q_lead = lowest_coefficient_mod(
    q.subs({s: qr_E2_s, t: qr_E2_t}),
    epsilon,
    qr_center_polynomial,
    alpha,
)
qr_E2_r_order, qr_E2_r_lead = lowest_coefficient_mod(
    r.subs({s: qr_E2_s, t: qr_E2_t}),
    epsilon,
    qr_center_polynomial,
    alpha,
)
assert (qr_E2_d_order, qr_E2_q_order, qr_E2_r_order) == (0, 2, 2)

qr_E2_A_deformation = P.subs(
    {s: qr_E2_s, t: qr_E2_t, X: qr_A + epsilon * Y}
)
qr_E2_A_order, qr_E2_A_residual = lowest_coefficient_mod(
    qr_E2_A_deformation,
    epsilon,
    qr_center_polynomial,
    alpha,
)
expected_qr_E2_A = (
    (48 * alpha**2 + 42 * alpha - 26) * Y**3
    + (-40 * alpha**2 - 48 * alpha + 6) * Y**2
    + (
        34 * alpha**2 * z
        + 200 * alpha**2
        + 90 * alpha * z
        + 480 * alpha
        + 59 * z
        + 290
    )
    * Y
    - 188 * alpha**2 * z
    - 1152 * alpha**2
    - 474 * alpha * z
    - 2936 * alpha
    - 300 * z
    - 1880
)
assert qr_E2_A_order == 3
assert sp.expand(qr_E2_A_residual - expected_qr_E2_A) == 0
qr_E2_A_discriminant = reduce_at_qr_center(
    sp.discriminant(qr_E2_A_residual, Y)
)
assert qr_E2_A_discriminant != 0
assert reduce_at_qr_center(
    qr_E2_A_discriminant
    - qr_discriminant_A_unit * qr_E2_q_lead**2 * qr_E2_r_lead
) == 0
qr_E2_A_derivative = P_X.subs(
    {s: qr_E2_s, t: qr_E2_t, X: qr_A + epsilon * Y}
)
qr_E2_A_derivative_order, qr_E2_A_derivative_residual = (
    lowest_coefficient_mod(
        qr_E2_A_derivative,
        epsilon,
        qr_center_polynomial,
        alpha,
    )
)
assert qr_E2_A_derivative_order == 2
assert sp.expand(
    qr_E2_A_derivative_residual - sp.diff(qr_E2_A_residual, Y)
) == 0

qr_E2_B_deformation = P.subs(
    {s: qr_E2_s, t: qr_E2_t, X: qr_B + epsilon * Y}
)
qr_E2_B_order, qr_E2_B_residual = lowest_coefficient_mod(
    qr_E2_B_deformation,
    epsilon,
    qr_center_polynomial,
    alpha,
)
expected_qr_E2_B = -(
    qr_B_quadratic_coefficient * Y**2
    + (432 * alpha**2 + 264 * alpha - 256) * Y
    + 2364 * alpha**2 * z
    - 8016 * alpha**2
    + 1786 * alpha * z
    - 5664 * alpha
    - 1653 * z
    + 5392
) / 2
assert qr_E2_B_order == 2
assert sp.expand(qr_E2_B_residual - expected_qr_E2_B) == 0
qr_E2_B_discriminant = reduce_at_qr_center(
    sp.discriminant(qr_E2_B_residual, Y)
)
assert qr_E2_B_discriminant != 0
assert reduce_at_qr_center(
    qr_E2_B_discriminant
    - qr_discriminant_B_unit * qr_E2_r_lead
) == 0
qr_E2_B_derivative = P_X.subs(
    {s: qr_E2_s, t: qr_E2_t, X: qr_B + epsilon * Y}
)
qr_E2_B_derivative_order, qr_E2_B_derivative_residual = (
    lowest_coefficient_mod(
        qr_E2_B_derivative,
        epsilon,
        qr_center_polynomial,
        alpha,
    )
)
assert qr_E2_B_derivative_order == 1
assert sp.expand(
    qr_E2_B_derivative_residual - sp.diff(qr_E2_B_residual, Y)
) == 0

qr_tangent_block = (
    (0, 2, 2, 2),
    (0, 1, 1, 1),
    (0, 2, 2, 1),
    (0, 2, 2, 2),
    (0, 2, 2, 2),
    (0, 2, 2, 2),
    (0, 2, 2, 1),
    (0, 2, 2, 1),
)
assert matrix[39:63] == qr_tangent_block * 3
assert sum(row[3] for row in qr_tangent_block[:3]) == 4
assert sum(row[3] for row in qr_tangent_block[3:]) == 8
print("PASS: the cubic q-r orbit adds twenty-four exact exceptional colors")

# Globalize the two geometric crossing slopes over q=0.  The normalization
# of the affine nodal q-curve is rational:
#
#   t=(y^2-9)/8,  s=4(y+2)/((y+1)(y+3)).
#
# The node has the two conductor preimages y^2=5.  Pulling P back to this
# normalization exhibits the repeated root X=2/(y+3) globally.
q_parameter, q_selector = sp.symbols("q_parameter q_selector")
q_normal_t = (q_parameter**2 - 9) / 8
q_normal_s = (
    4
    * (q_parameter + 2)
    / ((q_parameter + 1) * (q_parameter + 3))
)
q_inverse_parameter = (
    8 * t * (t + 1) * s + 8 * t + 6
) / (2 * (2 * t + 1))
assert sp.cancel(q.subs({s: q_normal_s, t: q_normal_t})) == 0
assert sp.cancel(
    q_inverse_parameter.subs({s: q_normal_s, t: q_normal_t})
    - q_parameter
) == 0

q_repeated_root = 2 / (q_parameter + 3)
q_linear_factor = X * (q_parameter + 1) ** 2 + 2 * (q_parameter + 3)
q_quadratic_factor = (
    16 * X**2
    + X * q_parameter**4
    - 14 * X * q_parameter**2
    - 16 * X * q_parameter
    - 3 * X
    + 2 * q_parameter**3
    + 10 * q_parameter**2
    + 14 * q_parameter
    + 6
)
q_root_factorization = (
    (X * (q_parameter + 3) - 2) ** 2
    * q_linear_factor
    * q_quadratic_factor
    / (
        16
        * (q_parameter + 1) ** 2
        * (q_parameter + 3) ** 2
    )
)
assert sp.cancel(
    P.subs({s: q_normal_s, t: q_normal_t})
    - q_root_factorization
) == 0

# Holding t fixed and moving transversely in s gives the quadratic normal
# slope polynomial.  Its discriminant has square class
# (y-5)/(y+3), so the two slopes do not split over Q(y).  They instead form
# one connected quadratic cover, ramified at y=5 and y=-3.
q_slope_deformation = sp.cancel(
    P.subs(
        {
            s: q_normal_s + epsilon,
            t: q_normal_t,
            X: q_repeated_root + epsilon * Y,
        }
    )
)
q_slope_expansion = sp.Poly(q_slope_deformation, epsilon)
assert q_slope_expansion.nth(0) == 0
assert q_slope_expansion.nth(1) == 0
q_slope_residual = sp.factor(q_slope_expansion.nth(2))
assert sp.Poly(q_slope_residual, Y).degree() == 2
assert sp.factor(sp.discriminant(q_slope_residual, Y)) == (
    q_parameter**2
    * (q_parameter - 5)
    * (q_parameter**2 + 4 * q_parameter + 5) ** 2
    / (q_parameter + 3) ** 5
)

q_branch_y = (
    5 + 3 * q_selector**2
) / (1 - q_selector**2)
assert sp.cancel(
    ((q_parameter - 5) / (q_parameter + 3)).subs(
        q_parameter, q_branch_y
    )
    - q_selector**2
) == 0

# The node conductor y^2=5 pulls back to the irreducible quartic below.
# Removing the affine points at infinity and the conductor gives the two
# unit lattices
#
# Q[y,1/((y+1)(y+3)(y^2-5))]^*/Q^* = Z^3,
# Q[w,1/((w-1)(w+1)(w^2+3)c(w))]^*/Q^* = Z^4.
q_branch_conductor = q_selector**4 + 10 * q_selector**2 + 5
assert sp.factor(q_branch_conductor) == q_branch_conductor
assert sp.cancel(
    (q_parameter**2 - 5).subs(q_parameter, q_branch_y)
    - 4
    * q_branch_conductor
    / ((q_selector - 1) ** 2 * (q_selector + 1) ** 2)
) == 0
assert sp.cancel(
    (q_parameter + 1).subs(q_parameter, q_branch_y)
    - 2 * (q_selector**2 + 3) / (1 - q_selector**2)
) == 0
assert sp.cancel(
    (q_parameter + 3).subs(q_parameter, q_branch_y)
    - 8 / (1 - q_selector**2)
) == 0

# Rows are the cover-unit basis
# (w-1,w+1,w^2+3,c(w)); columns are pullbacks of
# (y+1,y+3,y^2-5).  The cokernel is free of rank one.  Adjoining w-1,
# rather than merely the anti-invariant ratio (w-1)/(w+1), completes the
# lattice unimodularly.
q_conductor_unit_pullback = (
    (-1, -1, -2),
    (-1, -1, -2),
    (1, 0, 0),
    (0, 0, 1),
)
assert sp.Matrix(q_conductor_unit_pullback).rank() == 3
assert smith_invariant_factors(q_conductor_unit_pullback) == (1, 1, 1)
q_selector_column = sp.Matrix((1, 0, 0, 0))
q_unit_completion = sp.Matrix(q_conductor_unit_pullback).row_join(
    q_selector_column
)
assert q_unit_completion.det() == -1
assert smith_invariant_factors(tuple(map(tuple, q_unit_completion.tolist()))) == (
    1,
    1,
    1,
    1,
)
q_anti_invariant_completion = sp.Matrix(
    q_conductor_unit_pullback
).row_join(sp.Matrix((1, -1, 0, 0)))
assert abs(q_anti_invariant_completion.det()) == 2
assert sp.factor(sp.together(
    (q_selector - 1) * (-q_selector - 1)
    - 8 / (q_branch_y + 3)
)) == 0
print("PASS: the global q crossing is one rational conductor double cover")
print("PASS: one q selector completes the conductor unit lattice unimodularly")

# General colored proportional-row theorem.  If generator rows are a*p and
# b*p, every Laurent monomial in those generators has orders satisfying
# b*t_i=a*t_j.  The helper returns one exact mismatch per proportionality
# class and also detects nonzero targets on zero generator rows.
synthetic_proportionality = colored_proportionality_witnesses(
    ((1, 2), (3, 6), (0, 0)),
    (1, 4, 2),
    ("C_1", "C_2", "C_0"),
)
assert synthetic_proportionality == (
    {
        "kind": "zero_generator_row",
        "colors": ("C_0",),
        "generator_rows": ((0, 0),),
        "target_orders": (2,),
        "mismatch": 2,
    },
    {
        "kind": "proportional_generator_rows",
        "colors": ("C_1", "C_2"),
        "generator_rows": ((1, 2), (3, 6)),
        "row_scales": (1, 3),
        "target_orders": (1, 4),
        "mismatch": -1,
    },
)
synthetic_span_datum = ToroidalBoundaryDatum(
    fan=ToroidalFanDatum(
        lattice_basis=("ord_D",),
        rays=(TropicalRay("D_ray", (1,)),),
        maximal_cones=(("D_ray",),),
        incidence_certificate="synthetic two-color divisor-span regression",
    ),
    boundary_colors=(
        BoundaryColor("D_1", SheetColor.BOUNDARY, "D_ray"),
        BoundaryColor("D_2", SheetColor.BOUNDARY, "D_ray"),
    ),
    valuation_functions=("g",),
    valuation_matrix=((2,), (4,)),
    valuation_certificate="exact synthetic orders",
    divisor_span_problems=(
        ColoredDivisorSpanProblem(
            "torsion_target",
            ("g",),
            (1, 2),
            True,
            "the architecture has only the generator g",
        ),
        ColoredDivisorSpanProblem(
            "integral_target",
            ("g",),
            (4, 8),
            True,
            "the architecture has only the generator g",
        ),
    ),
    nonlinear_residue=(),
)
synthetic_span_audit = audit_toroidal_boundary(synthetic_span_datum)
synthetic_span_profiles = invariant(
    synthetic_span_audit, "divisor_span_problems"
)
assert synthetic_span_profiles["torsion_target"] == {
    "generator_functions": ("g",),
    "matrix_rank": 1,
    "augmented_rank": 1,
    "target_class_order": 2,
    "in_integral_span": False,
    "proportionality_witnesses": (),
}
assert synthetic_span_profiles["integral_target"] == {
    "generator_functions": ("g",),
    "matrix_rank": 1,
    "augmented_rank": 1,
    "target_class_order": 1,
    "in_integral_span": True,
    "proportionality_witnesses": (),
}
assert synthetic_span_audit.status is TropicalFeasibilityStatus.OBSTRUCTED
print("PASS: proportional colored rows give an unbounded divisor-span theorem")
print("PASS: the divisor-span gate distinguishes torsion from membership")

# First Cox/mask feasibility problem.  It is deliberately narrow but
# exhaustive for masks d^a*q^b*r^c*(w-1)^e with nonnegative base exponents:
# the selector parity is reduced by its descended norm, and exponents above
# eight overshoot the derivative target.  The unramified d-color forces
# a=0, while the index-four d-color then demands 4a=3, so no model exists.
f20_base_mask = audit_toroidal_boundary(f20_base_factor_mask_datum())
assert f20_base_mask.status is TropicalFeasibilityStatus.OBSTRUCTED
f20_span_problem = invariant(f20_base_mask, "divisor_span_problems")[
    "uncolored_base_factor_plus_q_selector_span"
]
assert (
    f20_span_problem["matrix_rank"],
    f20_span_problem["augmented_rank"],
    f20_span_problem["target_class_order"],
    f20_span_problem["in_integral_span"],
) == (3, 4, None, False)
assert tuple(
    witness["colors"]
    for witness in f20_span_problem["proportionality_witnesses"]
) == (
    ("d_unramified", "d_ramified_4"),
    ("q_collision_plus", "q_residual_1"),
    ("r_ramified_2_plus", "r_unramified"),
    ("triple_plus_E1_ramified_4", "triple_plus_E1_simple"),
    ("triple_plus_E2_cluster_1", "triple_plus_E2_simple"),
    (
        "qr_tangent_1_E1_A_ramified_2",
        "qr_tangent_1_E1_B_ramified_2",
    ),
)
f20_base_mask_problem = invariant(f20_base_mask, "feasibility_problems")[
    "base_factor_plus_q_selector_derivative_principalization"
]
assert f20_base_mask_problem == {
    "variables": (
        "exponent_d",
        "exponent_q",
        "exponent_r",
        "selector_parity",
    ),
    "search_size": 9**3 * 2,
    "model_count": 0,
    "minimal_models": (),
}
f20_candidate_matrix = invariant(f20_base_mask, "valuation_matrix")["matrix"]
assert f20_candidate_matrix[0][4:] == (1, 0, 0, 0)
assert f20_candidate_matrix[0][3] == 0
assert f20_candidate_matrix[1][4:] == (4, 0, 0, 0)
assert f20_candidate_matrix[1][3] == 3
assert "toroidal.divisor_span_obstruction" in obstruction_codes(f20_base_mask)
assert "toroidal.feasibility_infeasible" in obstruction_codes(f20_base_mask)
print("PASS: six F20 ray classes require genuinely colored Cox generators")
print("PASS: the exhaustive base-factor-plus-selector mask box is empty")

# Stage three is conditional on surviving mask models.  There are none in
# this certified architecture, so no adjugate division or affine-recognition
# claim is made.  Keeping the two report tuples empty is an exact guard
# against accidentally promoting this scoped obstruction to a global no-go.
f20_mask_survivors = f20_base_mask_problem["minimal_models"]
f20_inverse_adjugate_reports = tuple(
    {"model": model, "entrywise_divisible": None}
    for model in f20_mask_survivors
)
f20_affine_recognition_reports = tuple(
    {"model": model, "source_is_affine_space": None}
    for model in f20_mask_survivors
)
assert f20_inverse_adjugate_reports == ()
assert f20_affine_recognition_reports == ()
print("PASS: no scoped survivor reaches adjugate or affine-space testing")


positive_support_problem = ValuationFeasibilityProblem(
    name="positive_support_weights",
    variables=(
        IntegralVariable("w_W", 0, 2),
        IntegralVariable("w_K", 0, 2),
        IntegralVariable("w_L", 0, 2),
    ),
    fixed_function_coefficients=(0, 0, 0),
    variable_function_coefficients=(
        (0, 0, 0),
        (0, 0, 0),
        (0, 0, 0),
    ),
    target_orders=(0, 0, 0),
    constraints=(
        LinearConstraint("positive_W", (1, 0, 0), ">=", 1),
        LinearConstraint("positive_K", (0, 1, 0), ">=", 1),
        LinearConstraint("positive_L", (0, 0, 1), ">=", 1),
    ),
)
positive_support = audit_toroidal_boundary(
    replace(
        a4_toroidal_ledger_datum(),
        feasibility_problems=(positive_support_problem,),
    )
)
assert positive_support.status is TropicalFeasibilityStatus.FEASIBLE
assert invariant(positive_support, "feasibility_problems")[
    "positive_support_weights"
]["minimal_models"] == ({"w_W": 1, "w_K": 1, "w_L": 1},)
print("PASS: integral inequalities return the primitive positive support")


bad_identity = audit_toroidal_boundary(
    replace(
        a4_toroidal_ledger_datum(),
        valuation_matrix=((2, 1, 4), (3, 0, 3), (1, 1, 2)),
    )
)
assert bad_identity.status is TropicalFeasibilityStatus.OBSTRUCTED
assert "toroidal.valuation_identity" in obstruction_codes(bad_identity)
print("PASS: a coefficientwise tropical ledger mismatch is rejected")


torsion_screen = audit_toroidal_boundary(
    replace(
        davenport_toroidal_ledger_datum(),
        valuation_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 4)),
        unimodular_blocks=(),
        affine_screen=ToroidalAffineScreen(
            boundary_colors=("E3", "E6", "J"),
            unit_functions=(
                "pullback_Delta",
                "jacobian_J",
                "primitive_E3_character",
            ),
            normal_core_certificate=(
                "synthetic normal UFD core with complete three-prime boundary"
            ),
        ),
    )
)
assert torsion_screen.status is TropicalFeasibilityStatus.OBSTRUCTED
assert invariant(torsion_screen, "affine_screen")[
    "class_group_torsion"
] == (4,)
assert "toroidal.affine_class_group" in obstruction_codes(torsion_screen)
print("PASS: the shared Smith screen rejects class-group torsion Z/4")


empty_problem = ValuationFeasibilityProblem(
    name="certified_empty_box",
    variables=(IntegralVariable("x", 0, 2),),
    fixed_function_coefficients=(0, 0, 0),
    variable_function_coefficients=((1, 0, 0),),
    target_orders=(1, 1, 1),
    infeasibility_is_obstruction=True,
    exhaustive_scope_certificate="the declared ansatz is exactly x=0,1,2",
)
empty_audit = audit_toroidal_boundary(
    replace(
        a4_toroidal_ledger_datum(),
        feasibility_problems=(empty_problem,),
    )
)
assert empty_audit.status is TropicalFeasibilityStatus.OBSTRUCTED
assert "toroidal.feasibility_infeasible" in obstruction_codes(empty_audit)
print("PASS: only certificate-scoped integral infeasibility is obstructing")


summary = {
    "A4": a4.to_dict(),
    "D5": d5.to_dict(),
    "Davenport": davenport.to_dict(),
    "F20": f20.to_dict(),
}
print(json.dumps(summary, indent=2, sort_keys=True))
print("PASS toroidal boundary feasibility compiler")
