#!/usr/bin/env python3
"""Exact weighted coefficient scheme at the foundational cubic map."""

from itertools import product

import sympy as sp


x, y, z, u, v, tau = sp.symbols("x y z u v tau")

# Weight/degree/linearity constraints give 7+6+3=16 monomials.
source_weights = (1, -1, -2)
output_data = ((-2, 7), (-1, 6), (1, 4))
supports = []
for output_weight, degree_bound in output_data:
    support = []
    for i, j, k in product(range(degree_bound + 1), repeat=3):
        if i + j + k > degree_bound or k > 1:
            continue
        if i - j - 2 * k == output_weight:
            support.append((i, j, k))
    supports.append(tuple(support))
assert tuple(map(len, supports)) == (7, 6, 3)
assert sum(map(len, supports)) == 16

A20, A11, A30, A21, A40, A31 = sp.symbols(
    "A20 A11 A30 A21 A40 A31"
)
B01, B20, B11, B30, B21 = sp.symbols("B01 B20 B11 B30 B21")
C10, C01 = sp.symbols("C10 C01")

A = (
    v
    + A20 * u**2
    + A11 * u * v
    + A30 * u**3
    + A21 * u**2 * v
    + A40 * u**4
    + A31 * u**3 * v
)
B = (
    u
    + B01 * v
    + B20 * u**2
    + B11 * u * v
    + B30 * u**3
    + B21 * u**2 * v
)
C = 2 + C10 * u + C01 * v

reduced_matrix = sp.Matrix(
    [
        (-2 * A, sp.diff(A, u), sp.diff(A, v)),
        (-B, sp.diff(B, u), sp.diff(B, v)),
        (C, sp.diff(C, u), sp.diff(C, v)),
    ]
)
keller_polynomial = sp.Poly(sp.expand(reduced_matrix.det() + 2), u, v)
keller_equations = [coefficient for _monomial, coefficient in keller_polynomial.terms()]
assert len(keller_equations) == 16

# Normalize the two effective diagonal parameters by C10=-3, C01=-1.
normalized_equations = [
    sp.expand(equation.subs({C10: -3, C01: -1}))
    for equation in keller_equations
]
coefficient_variables = (
    A20,
    A11,
    A30,
    A21,
    A40,
    A31,
    B20,
    B11,
    B30,
    B21,
    B01,
)
computed = sp.groebner(
    normalized_equations,
    *coefficient_variables,
    order="grevlex",
    domain=sp.QQ,
)

triangular_relations = (
    (B01 - 3) ** 2,
    12 * A20 - 7 * B01 - 27,
    2 * A11 - B01 - 3,
    4 * A30 - 9 * B01 - 1,
    A21 - B01,
    2 * A40 - 3 * B01 + 3,
    2 * A31 - B01 + 1,
    2 * B20 - 11 * B01 + 9,
    B11 - 3 * B01 + 3,
    B30 - 6 * B01 + 9,
    B21 - 2 * B01 + 3,
)
triangular = sp.groebner(
    triangular_relations,
    *coefficient_variables,
    order="grevlex",
    domain=sp.QQ,
)
assert all(triangular.reduce(poly.as_expr())[1] == 0 for poly in computed.polys)
assert all(computed.reduce(poly.as_expr())[1] == 0 for poly in triangular.polys)

# Extract F + epsilon H.
source_u = x * y
source_v = x**2 * z
foundational_substitution = {
    A20: 4,
    A11: 3,
    A30: 7,
    A21: 3,
    A40: 3,
    A31: 1,
    B01: 3,
    B20: 12,
    B11: 6,
    B30: 9,
    B21: 3,
    C10: -3,
    C01: -1,
}
F = sp.Matrix(
    [
        sp.cancel(x**-2 * A.subs({u: source_u, v: source_v})),
        sp.cancel(x**-1 * B.subs({u: source_u, v: source_v})),
        sp.expand(x * C.subs({u: source_u, v: source_v})),
    ]
).subs(foundational_substitution)

H = sp.Matrix(
    [
        sp.Rational(7, 12) * y**2
        + sp.Rational(1, 2) * x * y * z
        + sp.Rational(9, 4) * x * y**3
        + x**2 * y**2 * z
        + sp.Rational(3, 2) * x**2 * y**4
        + sp.Rational(1, 2) * x**3 * y**3 * z,
        x * z
        + sp.Rational(11, 2) * x * y**2
        + 3 * x**2 * y * z
        + 6 * x**2 * y**3
        + 2 * x**3 * y**2 * z,
        0,
    ]
)

variables = (x, y, z)
deformed_determinant = sp.Poly(
    sp.expand((F + tau * H).jacobian(variables).det() + 2),
    tau,
)
assert deformed_determinant.coeff_monomial(tau) == 0
quadratic_obstruction = sp.factor(deformed_determinant.coeff_monomial(tau**2))
assert quadratic_obstruction != 0
assert deformed_determinant.degree() == 2

# H is not tangent to the affine left-right orbit.
left_parameters = sp.symbols("left0:9")
right_parameters = sp.symbols("right0:9")
left_translation = sp.symbols("left_translation0:3")
right_translation = sp.symbols("right_translation0:3")
left_matrix = sp.Matrix(3, 3, left_parameters)
right_matrix = sp.Matrix(3, 3, right_parameters)
affine_variation = (
    left_matrix * F
    + sp.Matrix(left_translation)
    + F.jacobian(variables)
    * (right_matrix * sp.Matrix(variables) + sp.Matrix(right_translation))
)
parameters = (
    *left_parameters,
    *right_parameters,
    *left_translation,
    *right_translation,
)
linear_equations = []
for component in range(3):
    polynomial = sp.Poly(sp.expand(affine_variation[component] - H[component]), *variables)
    linear_equations.extend(polynomial.coeffs())
matrix, vector = sp.linear_eq_to_matrix(linear_equations, parameters)
assert matrix.rank() == 23
assert matrix.row_join(vector).rank() == 24

# Two exact reduced triangular families on the omitted p=q=0 boundary.
lam, alpha2, alpha3, alpha4 = sp.symbols("lam alpha2 alpha3 alpha4")
boundary_families = (
    (v, u + lam * v, 2),
    (v + alpha2 * u**2 + alpha3 * u**3 + alpha4 * u**4, u, 2),
)
for boundary_A, boundary_B, boundary_C in boundary_families:
    boundary_matrix = sp.Matrix(
        [
            (-2 * boundary_A, sp.diff(boundary_A, u), sp.diff(boundary_A, v)),
            (-boundary_B, sp.diff(boundary_B, u), sp.diff(boundary_B, v)),
            (boundary_C, sp.diff(boundary_C, u), sp.diff(boundary_C, v)),
        ]
    )
    assert sp.expand(boundary_matrix.det()) == -2

# Diagonal gauge normalizes the two one-sided nonconstant-C boundary charts
# to (C10,C01)=(0,-1) and (-3,0).  Both specialized Keller ideals are units.
for one_sided_C in (2 - v, 2 - 3 * u):
    one_sided_matrix = sp.Matrix(
        [
            (-2 * A, sp.diff(A, u), sp.diff(A, v)),
            (-B, sp.diff(B, u), sp.diff(B, v)),
            (
                one_sided_C,
                sp.diff(one_sided_C, u),
                sp.diff(one_sided_C, v),
            ),
        ]
    )
    one_sided_equations = [
        coefficient
        for _monomial, coefficient in sp.Poly(
            sp.expand(one_sided_matrix.det() + 2), u, v
        ).terms()
    ]
    one_sided_ideal = sp.groebner(
        one_sided_equations,
        *coefficient_variables,
        order="grevlex",
        domain=sp.QQ,
    )
    assert len(one_sided_ideal.polys) == 1
    assert one_sided_ideal.polys[0].as_expr() == 1

print("PASS: weights, degree bounds, and z-linearity give exactly 16 monomials")
print("PASS: normalized Keller scheme is Q[epsilon]/(epsilon^2)")
print("PASS: explicit H is first-order Keller and has a nonzero quadratic obstruction")
print("PASS: H is not tangent to the affine left-right orbit")
print("PASS: two reduced triangular automorphism families lie on p=q=0")
print("PASS: both one-sided nonconstant-C boundary charts have unit ideal")
