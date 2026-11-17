#!/usr/bin/env python3
"""Exact obstructions for the direct one-variable Schur route to HC(4).

There are two complementary calculations.

First, classify auxiliary coordinates for the PC(2) graph.  In a coordinate
graph chart (q,m), omit one q-coordinate and adjoin a polynomial function of
the graph parameters h=(X,Y,W,D).  Collision transfer first leaves four
charts.  Exact coefficient-ideal and two-variable Jacobian-mate obstructions
then leave only charts 0010 and 0011, in both cases with the fourth
q-coordinate omitted.  Every polynomial slice there is W up to scaling and
a polynomial in the three retained coordinates.

For either surviving chart, the exact Lagrangian primitive gives a canonical
generating family

    F(w,t) = P(w) + A(w) (t-O(w)),

where t=O(w) is the omitted graph equation and A is its complementary
coordinate.  Every quadratic Schur ascent preserving the graph to first
order is

    F_K = F + K(w) (t-O(w))^2.

At the certified symmetric collision, partial Legendre descent at s=A(P+)
would give

    psi_K = P-sO-(A-s)^2/(4K).

The irreducibility of A-s classifies every polynomial choice making psi_K
polynomial: K is a scalar times (A-s)^e for e=0,1,2.  A handful of exact
Hessian evaluations rules out constant nonzero determinant in every case.

Second, test a further linear Schur descent of the Meng--Yang five-variable
potential itself.  Its only constant-direction quadratic restrictions are
directions in the three dual variables.  Exact divisibility equations then
show that no nonzero such direction has a polynomial partial Legendre
transform at any constant dual value.

These are restricted negative results.  They do not rule out nonlinear
auxiliary coordinates, non-coordinate graph embeddings, or higher-degree
generating families.
"""

from __future__ import annotations

from itertools import product
import runpy

import sympy as sp


graph = runpy.run_path("scripts/search_hc4_graph_polarizations.py")
h = tuple(graph["h_variables"])
X, Y, W, D = h
Q = list(graph["position_coordinates_h"])
M = list(graph["momentum_coordinates_h"])


def chart_coordinates(chart: str) -> tuple[sp.Matrix, sp.Matrix]:
    mask = tuple(int(bit) for bit in chart)
    q = sp.Matrix([M[i] if mask[i] else Q[i] for i in range(4)])
    m = sp.Matrix([-Q[i] if mask[i] else M[i] for i in range(4)])
    return q, m


def exact_primitive(q: sp.Matrix, m: sp.Matrix) -> sp.Expr:
    """Integrate the closed graph one-form sum(m_i dq_i) in h-coordinates."""

    one_form = (q.jacobian(h).T * m).applyfunc(sp.expand)
    assert (
        one_form.jacobian(h) - one_form.jacobian(h).T
    ).applyfunc(sp.expand) == sp.zeros(4)

    primitive = sp.integrate(one_form[0], h[0])
    for index, variable in enumerate(h[1:], start=1):
        residual = sp.expand(one_form[index] - sp.diff(primitive, variable))
        assert all(
            sp.diff(residual, h[earlier]) == 0 for earlier in range(index)
        )
        primitive = sp.expand(primitive + sp.integrate(residual, variable))

    assert all(
        sp.expand(sp.diff(primitive, variable) - one_form[index]) == 0
        for index, variable in enumerate(h)
    )
    return primitive


def jacobian_derivation(q: sp.Matrix, omitted: int) -> tuple[sp.Expr, ...]:
    """Return V(h_i)=det d(q_without_omitted,h_i)/dh."""

    retained_jacobian = sp.Matrix(
        [q[index] for index in range(4) if index != omitted]
    ).jacobian(h)
    components = []
    for index in range(4):
        coordinate_row = [0] * 4
        coordinate_row[index] = 1
        components.append(
            sp.factor(
                sp.Matrix(
                    list(retained_jacobian.tolist()) + [coordinate_row]
                ).det(method="domain-ge")
            )
        )
    return tuple(components)


# A literal prescription of a Schur complement as the PC(2) Jacobian cannot
# work: every Hessian Schur complement is symmetric, while this Jacobian is
# not.  The displayed skew entry is already 16/3 at (x,q,p,z)=(0,0,0,1).
source_variables = tuple(graph["source_variables"])
pc2_map = sp.Matrix([graph[name] for name in ("R", "T", "D", "S")])
pc2_jacobian = pc2_map.jacobian(source_variables)
skew_entry = sp.expand(pc2_jacobian[1, 0] - pc2_jacobian[0, 1])
skew_point = dict(zip(source_variables, (0, 0, 0, 1), strict=True))
assert skew_entry.subs(skew_point) == sp.Rational(16, 3)


# Classify all coordinate-chart projections that become polynomial
# coordinates after adjoining a linear auxiliary function of h.
auxiliary_coefficients = sp.symbols("a0:4")
linear_auxiliary = sum(
    coefficient * variable
    for coefficient, variable in zip(auxiliary_coefficients, h, strict=True)
)
surviving_augmented_charts: set[tuple[str, int]] = set()

for mask in product((0, 1), repeat=4):
    chart = "".join(map(str, mask))
    q, _ = chart_coordinates(chart)
    for omitted in range(4):
        retained = [q[index] for index in range(4) if index != omitted]
        determinant = sp.expand(
            sp.Matrix(retained + [linear_auxiliary])
            .jacobian(h)
            .det(method="domain-ge")
        )
        determinant_poly = sp.Poly(determinant, *h)
        nonconstant_rows: list[list[sp.Expr]] = []
        constant_row = sp.zeros(1, 4)
        for monomial, coefficient in determinant_poly.terms():
            row = [
                sp.expand(coefficient).coeff(parameter)
                for parameter in auxiliary_coefficients
            ]
            if sum(monomial) == 0:
                constant_row = sp.Matrix([row])
            else:
                nonconstant_rows.append(row)
        nonconstant_matrix = (
            sp.Matrix(nonconstant_rows)
            if nonconstant_rows
            else sp.zeros(0, 4)
        )
        if any(
            (constant_row * vector)[0] != 0
            for vector in nonconstant_matrix.nullspace()
        ):
            surviving_augmented_charts.add((chart, omitted))

assert surviving_augmented_charts == {("0010", 3), ("0011", 3)}
for chart in ("0010", "0011"):
    q, _ = chart_coordinates(chart)
    determinant = sp.factor(
        sp.Matrix([q[0], q[1], q[2], linear_auxiliary])
        .jacobian(h)
        .det(method="domain-ge")
    )
    assert sp.expand(
        determinant
        - X * auxiliary_coefficients[1] / 3
        + auxiliary_coefficients[2]
    ) == 0


# Upgrade the preceding linear calculation to arbitrary polynomial auxiliary
# coordinates, under the requirement that a pair from the certified
# three-point fiber survive.  Equal generating-family gradients require all
# four complementary chart coordinates to agree at the selected pair.
position_collision_values = graph["position_collision_values"]
momentum_collision_values = graph["momentum_collision_values"]
collision_preserving_charts: set[str] = set()
for mask in product((0, 1), repeat=4):
    chart = "".join(map(str, mask))
    complementary_values = []
    for point_index in range(3):
        complementary_values.append(
            tuple(
                -position_collision_values[point_index][index]
                if mask[index]
                else momentum_collision_values[point_index][index]
                for index in range(4)
            )
        )
    if any(
        complementary_values[left] == complementary_values[right]
        for left in range(3)
        for right in range(left + 1, 3)
    ):
        collision_preserving_charts.add(chart)
assert collision_preserving_charts == {"0000", "0001", "0010", "0011"}


# In charts 0000 and 0001 every omission has a nonunit common factor in the
# derivation components (or the zero derivation), so V(u) cannot be a unit.
for chart in ("0000", "0001"):
    q_coordinates, _ = chart_coordinates(chart)
    for omitted in range(4):
        components = jacobian_derivation(q_coordinates, omitted)
        component_gcd = sp.factor(sp.gcd_list(components))
        assert component_gcd not in (1, -1)


# Chart 0010, omission 0 has a fixed point at the graph origin; omission 2
# has a nonunit common factor.  Chart 0011, omission 0 has the displayed
# fixed point over Q(sqrt(34)); omission 2 again has a nonunit common factor.
origin_h = dict.fromkeys(h, 0)
q_0010, _ = chart_coordinates("0010")
components_0010_0 = jacobian_derivation(q_0010, 0)
assert all(component.subs(origin_h) == 0 for component in components_0010_0)
components_0010_2 = jacobian_derivation(q_0010, 2)
assert sp.factor(sp.gcd_list(components_0010_2)) not in (1, -1)

q_0011, _ = chart_coordinates("0011")
components_0011_0 = jacobian_derivation(q_0011, 0)
sqrt_34 = sp.sqrt(34)
fixed_point_0011 = {
    X: (sqrt_34 - 4) / 6,
    Y: 1,
    W: (84 + 3 * sqrt_34) / 25,
    D: 0,
}
assert all(
    sp.simplify(component.subs(fixed_point_0011)) == 0
    for component in components_0011_0
)
components_0011_2 = jacobian_derivation(q_0011, 2)
assert sp.factor(sp.gcd_list(components_0011_2)) not in (1, -1)


# For omission 1 the retained triples are (X,D,T) and (X,D,S).
# Both omitted-coordinate equations are linear in W:
#
#     f(Y,W)=a(Y)W+b(Y),
#
# with nonconstant a over Q(X,D).  If
# V=-a d/dY+(a'W+b')d/dW and V(u) is a unit, comparison of the top W-degree
# lets one subtract a polynomial in f until u=c(Y); then -a c'(Y) is a unit,
# impossible for nonconstant a.  The assertions below certify the two a's
# and the corresponding Jacobian derivations exactly.
U = 1 + X * Y
for chart, expected_w_coefficient in (
    ("0010", 3 * X * U**2),
    ("0011", U**3 / 2),
):
    q_coordinates, _ = chart_coordinates(chart)
    omitted_equation = q_coordinates[3]
    assert sp.Poly(omitted_equation, W).degree() == 1
    assert sp.expand(
        sp.diff(omitted_equation, W) - expected_w_coefficient
    ) == 0
    components = jacobian_derivation(q_coordinates, 1)
    assert components[0] == 0
    assert components[3] == 0
    assert sp.expand(components[1] + expected_w_coefficient) == 0
    assert sp.Poly(expected_w_coefficient, Y).degree() > 0


# For omission 3 both charts retain (X,Q,D).  In the polynomial coordinates
# (X,Q,D,W), V=-d/dW.  Thus all polynomial solutions of V(u)=c are
# u=-cW+f(X,Q,D), so W is the unique essential auxiliary coordinate modulo
# retained-coordinate gauge.
Q_base = Y + X * W / 3
for chart in ("0010", "0011"):
    q_coordinates, _ = chart_coordinates(chart)
    components = jacobian_derivation(q_coordinates, 3)

    def apply_derivation(expression: sp.Expr) -> sp.Expr:
        return sp.expand(
            sum(
                component * sp.diff(expression, variable)
                for component, variable in zip(components, h, strict=True)
            )
        )

    assert tuple(
        apply_derivation(expression)
        for expression in (X, Q_base, D, W)
    ) == (0, 0, 0, -1)


x, qvar, d, t, u, r = sp.symbols("x q d t u r")
w = (x, qvar, d, u)
coordinate_substitution = {
    X: x,
    Y: qvar - x * u / 3,
    W: u,
    D: d,
}


def determinant_at(
    potential: sp.Expr, point: tuple[int, int, int, int]
) -> sp.Expr:
    hessian = sp.hessian(potential, w)
    return sp.factor(hessian.subs(dict(zip(w, point, strict=True))).det())


schur_data: dict[str, dict[str, object]] = {}
for chart in ("0010", "0011"):
    q_coordinates, m_coordinates = chart_coordinates(chart)
    primitive_h = exact_primitive(q_coordinates, m_coordinates)
    primitive = sp.expand(
        primitive_h.subs(coordinate_substitution, simultaneous=True)
    )
    omitted_coordinate = sp.expand(
        q_coordinates[3].subs(coordinate_substitution, simultaneous=True)
    )
    complementary_coordinate = sp.expand(
        m_coordinates[3].subs(coordinate_substitution, simultaneous=True)
    )
    tau = sp.symbols(f"tau_{chart}")
    generating_family = sp.expand(
        primitive
        + complementary_coordinate * (tau - omitted_coordinate)
    )

    # The generating family has exactly the graph first jet on t=O(w).
    graph_substitution = {tau: omitted_coordinate}
    generating_gradient = sp.Matrix(
        [
            sp.diff(generating_family, variable)
            for variable in (x, qvar, d, tau, u)
        ]
    )
    expected_graph_gradient = sp.Matrix(
        [
            m_coordinates[0].subs(
                coordinate_substitution, simultaneous=True
            ),
            m_coordinates[1].subs(
                coordinate_substitution, simultaneous=True
            ),
            m_coordinates[2].subs(
                coordinate_substitution, simultaneous=True
            ),
            complementary_coordinate,
            0,
        ]
    )
    assert (
        generating_gradient.subs(graph_substitution, simultaneous=True)
        - expected_graph_gradient
    ).applyfunc(sp.expand) == sp.zeros(5, 1)

    collision_dual_value = (
        -sp.Rational(1, 8) if chart == "0010" else sp.Rational(0)
    )
    collision_lifts = (
        (sp.Rational(1), sp.Rational(2, 3), 0, sp.Rational(13, 2)),
        (sp.Rational(-1), sp.Rational(-2, 3), 0, sp.Rational(13, 2)),
    )
    collision_gradients = []
    for lift in collision_lifts:
        lift_w = dict(zip(w, lift, strict=True))
        lift_point = dict(lift_w)
        lift_point[tau] = omitted_coordinate.subs(lift_w)
        collision_gradients.append(
            tuple(
                sp.factor(entry.subs(lift_point))
                for entry in generating_gradient
            )
        )
    assert collision_gradients[0] == collision_gradients[1]
    assert collision_gradients[0][3] == collision_dual_value
    assert collision_gradients[0][4] == 0

    irreducible_factor = sp.Poly(
        complementary_coordinate - collision_dual_value,
        x,
        qvar,
        u,
        domain=sp.QQ,
    )
    assert irreducible_factor.is_irreducible

    schur_data[chart] = {
        "primitive": primitive,
        "omitted": omitted_coordinate,
        "factor": complementary_coordinate - collision_dual_value,
        "dual_value": collision_dual_value,
    }


# If K is polynomial and (A-s)^2/(4K) is polynomial, irreducibility forces
# K=c(A-s)^e, e=0,1,2.  It is therefore enough to test correction powers
# p=2-e in {0,1,2}.
data_0010 = schur_data["0010"]
base_0010 = sp.expand(
    data_0010["primitive"]
    - data_0010["dual_value"] * data_0010["omitted"]
)
factor_0010 = sp.expand(data_0010["factor"])

potential = base_0010  # p=0; subtracting a scalar has no Hessian effect.
assert determinant_at(potential, (0, 0, 0, 0)) == 0
assert determinant_at(potential, (1, 0, 0, 0)) == -sp.Rational(55, 4)

potential = sp.expand(base_0010 - r * factor_0010)  # p=1.
assert determinant_at(potential, (0, 0, 0, 0)) == 0

potential = sp.expand(base_0010 - r * factor_0010**2)  # p=2.
assert determinant_at(potential, (0, 0, 0, 0)) == -2 * r**2
assert determinant_at(potential, (0, 0, 0, 1)) == -10 * r**2


data_0011 = schur_data["0011"]
base_0011 = sp.expand(
    data_0011["primitive"]
    - data_0011["dual_value"] * data_0011["omitted"]
)
factor_0011 = sp.expand(data_0011["factor"])

potential = base_0011  # p=0.
assert determinant_at(potential, (0, 0, 0, 0)) == 1
assert determinant_at(potential, (1, 1, 0, 0)) == -sp.Rational(120632, 3)

potential = sp.expand(base_0011 - r * factor_0011)  # p=1.
det_origin = determinant_at(potential, (0, 0, 0, 0))
det_x = determinant_at(potential, (1, 0, 0, 0))
det_xq = determinant_at(potential, (1, 1, 0, 0))
assert det_origin == 1
assert det_x == 144 * r**2 + 120 * r + 1
assert sp.expand(
    det_xq - (936 * r**2 - 19319 * r - 120632) / 3
) == 0
assert sp.gcd(
    sp.Poly(det_x - det_origin, r, domain=sp.QQ),
    sp.Poly(det_xq - det_origin, r, domain=sp.QQ),
).degree() == 0

potential = sp.expand(base_0011 - r * factor_0011**2)  # p=2.
assert determinant_at(potential, (0, 0, 0, 0)) == 1
assert sp.expand(
    determinant_at(potential, (1, 0, 0, 0)) - (32 * r + 3) / 3
) == 0


# Direct second descent of the Meng--Yang potential.  Its notation is
# Psi=A^2+11A+2B, with the first two variables primal and the last three dual.
x1, x2, y1, y2, y3 = sp.symbols("x1 x2 y1 y2 y3")
meng_variables = (x1, x2, y1, y2, y3)
direction = sp.symbols("v0:5")
unit = 1 + x1 * x2
meng_a = y1 * unit**3 + 3 * x1 * y2 * unit**2 - x1**3 * y3
meng_b = (
    y1 * x2**2 * unit * (4 + 3 * x1 * x2)
    + y2 * (x2 + 3 * x1 * x2**2 * (4 + 3 * x1 * x2))
    + y3 * (2 * x1 - 3 * x1**2 * x2)
)
meng_potential = sp.expand(meng_a**2 + 11 * meng_a + 2 * meng_b)

directional_derivative = meng_potential
for _ in range(3):
    directional_derivative = sp.expand(
        sum(
            coefficient * sp.diff(directional_derivative, variable)
            for coefficient, variable in zip(
                direction, meng_variables, strict=True
            )
        )
    )
third_directional_coefficients = sp.Poly(
    directional_derivative, *meng_variables
).coeffs()
third_directional_basis = sp.groebner(
    third_directional_coefficients, *direction, order="grevlex"
)
third_basis_expressions = {
    sp.expand(polynomial.as_expr())
    for polynomial in third_directional_basis.polys
}
assert direction[0] ** 3 in third_basis_expressions
assert direction[1] ** 3 in third_basis_expressions
assert sp.expand(
    directional_derivative.subs({direction[0]: 0, direction[1]: 0})
) == 0


# In a dual direction (a,b,c), let L=D_v A and K=D_v B.  Polynomiality of
# the partial Legendre transform at a constant dual value sigma requires
# L | (2K-sigma).  The total degrees of L and K agree for every nonzero
# (a,b,c), so the quotient would be a scalar ell.  The coefficient ideal of
# 2K-sigma-ell*L has only the trivial direction.
a, b, c, ell, sigma = sp.symbols("a b c ell sigma")
linear_a = sp.expand(
    sum(
        coefficient * sp.diff(meng_a, variable)
        for coefficient, variable in zip(
            (a, b, c), (y1, y2, y3), strict=True
        )
    )
)
linear_b = sp.expand(
    sum(
        coefficient * sp.diff(meng_b, variable)
        for coefficient, variable in zip(
            (a, b, c), (y1, y2, y3), strict=True
        )
    )
)
linear_a_poly = sp.Poly(linear_a, x1, x2)
linear_b_poly = sp.Poly(linear_b, x1, x2)
assert linear_a_poly.coeff_monomial(x1**3 * x2**3) == a
assert linear_b_poly.coeff_monomial(x1**2 * x2**4) == 3 * a
assert linear_a_poly.coeff_monomial(x1**3 * x2**2) == 3 * b
assert linear_b_poly.coeff_monomial(x1**2 * x2**3) == 9 * b
assert linear_a_poly.coeff_monomial(x1**3) == -c
assert linear_b_poly.coeff_monomial(x1**2 * x2) == -3 * c
divisibility_coefficients = sp.Poly(
    2 * linear_b - sigma - ell * linear_a, x1, x2
).coeffs()
divisibility_basis = sp.groebner(
    divisibility_coefficients,
    ell,
    sigma,
    a,
    b,
    c,
    order="grevlex",
)
assert {
    sp.expand(polynomial.as_expr())
    for polynomial in divisibility_basis.polys
} == {sigma, a, b, c}


print("PASS: the PC(2) Jacobian is not itself a Hessian Schur complement")
print(
    "PASS: collision transfer and all-degree slice obstructions leave only "
    "charts 0010 and 0011 with auxiliary coordinate W modulo gauge"
)
print("PASS: both canonical generating families preserve the rational collision")
print(
    "PASS: every polynomial quadratic-pivot Schur descent in those two "
    "families has nonconstant or zero Hessian determinant"
)
print(
    "PASS: the Meng--Yang potential has no further polynomial linear "
    "quadratic-direction descent"
)
print(
    "SCOPE: non-coordinate graph embeddings and higher-degree critical "
    "equations remain open"
)
