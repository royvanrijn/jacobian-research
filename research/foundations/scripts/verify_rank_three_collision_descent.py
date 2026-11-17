#!/usr/bin/env python3
"""Exact algebra for the rank-three collision-framed descent audit.

The checker separates three statements:

* ordered distinct pairs in a cubic fiber give the full root frame, and the
  unique projective transformation between two framed root triples satisfies
  the descent cocycle;
* arbitrary quadratic Tschirnhaus changes have the displayed symmetric
  target formulas and exactly two boundary factors;
* the normalized linear--quadratic factorization map transports this
  projective action after target localization, while its denominator-free
  projective stabilizer is only the diagonal torus.

The last statement is deliberately a claim about the canonical
factorization transport.  It does not classify every possible nonlinear
polynomial self-equivalence of the foundational Keller map.
"""

from __future__ import annotations

from itertools import combinations, permutations

import sympy as sp


def assert_zero(expression: sp.Expr) -> None:
    """Assert an exact rational identity."""

    assert sp.cancel(expression) == 0


def permutation_sign(permutation: tuple[int, ...]) -> int:
    """Return the sign of a small permutation."""

    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def mobius_matrix(
    source_roots: tuple[sp.Expr, sp.Expr, sp.Expr],
    target_roots: tuple[sp.Expr, sp.Expr, sp.Expr],
) -> sp.Matrix:
    """Signed-maximal-minor matrix mapping three source roots to targets."""

    interpolation = sp.Matrix(
        [
            [
                source_roots[index],
                1,
                -target_roots[index] * source_roots[index],
                -target_roots[index],
            ]
            for index in range(3)
        ]
    )
    coefficients = []
    for column in range(4):
        remaining = [index for index in range(4) if index != column]
        coefficients.append(
            (-1) ** column * interpolation[:, remaining].det()
        )
    return sp.Matrix(
        [
            [coefficients[0], coefficients[1]],
            [coefficients[2], coefficients[3]],
        ]
    )


# ---------------------------------------------------------------------------
# 1. The cubic off-diagonal sheet is a frame torsor, and framed triples
#    determine a unique PGL_2 transition with an exact cocycle.
# ---------------------------------------------------------------------------

r = sp.symbols("r1:4")
u = sp.symbols("u1:4")
v = sp.symbols("v1:4")

root_vandermonde = (r[0] - r[1]) * (r[0] - r[2]) * (r[1] - r[2])
target_vandermonde = (u[0] - u[1]) * (u[0] - u[2]) * (u[1] - u[2])

root_to_target = mobius_matrix(r, u)
mobius_a, mobius_b, mobius_c, mobius_d = tuple(root_to_target)

for source_root, target_root in zip(r, u):
    assert_zero(
        mobius_a * source_root
        + mobius_b
        - target_root * (mobius_c * source_root + mobius_d)
    )

assert_zero(
    root_to_target.det() - root_vandermonde * target_vandermonde
)

# Simultaneous relabeling multiplies all four maximal minors by one sign.
# Hence the matrix is unchanged in PGL_2 and descends from the S_3 frame
# torsor.
for permutation in permutations(range(3)):
    permuted_matrix = mobius_matrix(
        tuple(r[index] for index in permutation),
        tuple(u[index] for index in permutation),
    )
    expected_sign = permutation_sign(permutation)
    for actual, original in zip(permuted_matrix, root_to_target):
        assert_zero(actual - expected_sign * original)

# The projective transition satisfies the cocycle.  Equality in PGL_2 is
# checked by the six 2-by-2 minors of the two flattened matrices.
target_to_third = mobius_matrix(u, v)
root_to_third = mobius_matrix(r, v)
composite = target_to_third * root_to_target
for first, second in combinations(range(4), 2):
    assert_zero(
        composite[first] * root_to_third[second]
        - composite[second] * root_to_third[first]
    )

# An ordered distinct pair (i,j) in a three-element set determines the
# missing third label, so there are 3*2=6 frames and S_3 acts simply
# transitively.
ordered_pairs = tuple(
    (first, second)
    for first in range(3)
    for second in range(3)
    if first != second
)
assert len(ordered_pairs) == 6
for ordered_pair in ordered_pairs:
    missing = tuple(set(range(3)) - set(ordered_pair))
    assert len(missing) == 1
frames = {
    ordered_pair + tuple(set(range(3)) - set(ordered_pair))
    for ordered_pair in ordered_pairs
}
assert frames == set(permutations(range(3)))


# ---------------------------------------------------------------------------
# 2. Exact quadratic Tschirnhaus target formulas, boundary ledger, and
#    composition.
# ---------------------------------------------------------------------------

q0, q1, q2 = sp.symbols("q0 q1 q2")
p0, p1, p2 = sp.symbols("p0 p1 p2")

e1 = sum(r)
e2 = r[0] * r[1] + r[0] * r[2] + r[1] * r[2]
e3 = r[0] * r[1] * r[2]


def symmetric_triple(
    values: tuple[sp.Expr, sp.Expr, sp.Expr],
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Elementary symmetric functions of three values."""

    return (
        sum(values),
        values[0] * values[1]
        + values[0] * values[2]
        + values[1] * values[2],
        values[0] * values[1] * values[2],
    )


def tschirnhaus_symmetric_transform(
    source_e1: sp.Expr,
    source_e2: sp.Expr,
    source_e3: sp.Expr,
    constant: sp.Expr,
    linear: sp.Expr,
    quadratic: sp.Expr,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Symmetric functions after u=q0+q1*r+q2*r^2."""

    anchored_e1 = (
        linear * source_e1
        + quadratic * (source_e1**2 - 2 * source_e2)
    )
    anchored_e2 = (
        linear**2 * source_e2
        + linear
        * quadratic
        * (source_e1 * source_e2 - 3 * source_e3)
        + quadratic**2
        * (source_e2**2 - 2 * source_e1 * source_e3)
    )
    anchored_e3 = source_e3 * (
        linear**3
        + linear**2 * quadratic * source_e1
        + linear * quadratic**2 * source_e2
        + quadratic**3 * source_e3
    )
    return (
        sp.expand(anchored_e1 + 3 * constant),
        sp.expand(
            anchored_e2
            + 2 * constant * anchored_e1
            + 3 * constant**2
        ),
        sp.expand(
            anchored_e3
            + constant * anchored_e2
            + constant**2 * anchored_e1
            + constant**3
        ),
    )


transformed_roots = tuple(q0 + q1 * root + q2 * root**2 for root in r)
actual_transformed_symmetric = symmetric_triple(transformed_roots)
formula_transformed_symmetric = tschirnhaus_symmetric_transform(
    e1, e2, e3, q0, q1, q2
)
for actual, formula in zip(
    actual_transformed_symmetric,
    formula_transformed_symmetric,
):
    assert_zero(actual - formula)

f1, f2, f3 = formula_transformed_symmetric
pi = 1 / e2
b = -e1 / e2
c = 2 * e3 / e2
transformed_pi = 1 / f2
transformed_b = -f1 / f2
transformed_c = 2 * f3 / f2

S = sp.symbols("S")
normalized_cubic = pi * S**3 + b * S**2 + S - c / 2
transformed_cubic = (
    transformed_pi * S**3
    + transformed_b * S**2
    + S
    - transformed_c / 2
)
for root in r:
    assert_zero(normalized_cubic.subs(S, root))
for root in transformed_roots:
    assert_zero(transformed_cubic.subs(S, root))

# The primitive-element boundary is the Vandermonde multiplier Theta.  The
# second boundary f2=0 is where the transformed relation cannot be normalized
# to have linear coefficient one.
theta = (
    q1**3
    + 2 * q1**2 * q2 * e1
    + q1 * q2**2 * (e1**2 + e2)
    + q2**3 * (e1 * e2 - e3)
)
transformed_vandermonde = (
    (transformed_roots[0] - transformed_roots[1])
    * (transformed_roots[0] - transformed_roots[2])
    * (transformed_roots[1] - transformed_roots[2])
)
assert_zero(transformed_vandermonde - root_vandermonde * theta)

old_discriminant = sp.discriminant(normalized_cubic, S)
assert_zero(old_discriminant - pi**4 * root_vandermonde**2)
assert_zero(
    transformed_pi**4 * transformed_vandermonde**2
    - old_discriminant * (transformed_pi / pi) ** 4 * theta**2
)

# A second quadratic change has the same formula on the already transformed
# elementary symmetric functions.  Comparing with the literal composite on
# roots proves the target cocycle without choosing numerical witnesses.
third_roots = tuple(p0 + p1 * root + p2 * root**2 for root in transformed_roots)
actual_third_symmetric = symmetric_triple(third_roots)
iterated_third_symmetric = tschirnhaus_symmetric_transform(
    f1, f2, f3, p0, p1, p2
)
for actual, formula in zip(actual_third_symmetric, iterated_third_symmetric):
    assert_zero(actual - formula)

# Reduce the literal composite to the basis 1,r,r^2.  These formulas are the
# coefficient-level composition law on the cubic presentation groupoid.
k0 = q0**2 + 2 * q1 * q2 * e3 + q2**2 * e1 * e3
k1 = (
    2 * q0 * q1
    - 2 * q1 * q2 * e2
    + q2**2 * (-e1 * e2 + e3)
)
k2 = (
    2 * q0 * q2
    + q1**2
    + 2 * q1 * q2 * e1
    + q2**2 * (e1**2 - e2)
)
composite_q0 = p0 + p1 * q0 + p2 * k0
composite_q1 = p1 * q1 + p2 * k1
composite_q2 = p1 * q2 + p2 * k2
for root, literal_value in zip(r, third_roots):
    reduced_value = (
        composite_q0 + composite_q1 * root + composite_q2 * root**2
    )
    assert_zero(literal_value - reduced_value)

composite_symmetric = tschirnhaus_symmetric_transform(
    e1,
    e2,
    e3,
    composite_q0,
    composite_q1,
    composite_q2,
)
for iterated, composed in zip(
    iterated_third_symmetric,
    composite_symmetric,
):
    assert_zero(iterated - composed)


# ---------------------------------------------------------------------------
# 3. Canonical target-localized PGL_2 transport of the normalized
#    linear--quadratic factorization map.
# ---------------------------------------------------------------------------

U, V = sp.symbols("U V")
linear_0, linear_1 = sp.symbols("linear_0 linear_1")
quad_0, quad_1, quad_2 = sp.symbols("quad_0 quad_1 quad_2")
matrix_a, matrix_b, matrix_c, matrix_d = sp.symbols(
    "matrix_a matrix_b matrix_c matrix_d"
)

linear_form = linear_0 * U + linear_1 * V
quadratic_form = quad_0 * U**2 + quad_1 * U * V + quad_2 * V**2
cubic_form = sp.expand(linear_form * quadratic_form)


def transform_binary_form(
    form: sp.Expr,
    matrix: sp.Matrix,
) -> sp.Expr:
    """Substitute a GL_2 change of binary variables."""

    return sp.expand(
        form.subs(
            {
                U: matrix[0, 0] * U + matrix[0, 1] * V,
                V: matrix[1, 0] * U + matrix[1, 1] * V,
            },
            simultaneous=True,
        )
    )


def tangent_coefficient(form: sp.Expr) -> sp.Expr:
    """The coefficient functional ell=[U V^2]."""

    return sp.Poly(sp.expand(form), U, V).coeff_monomial(U * V**2)


def linear_quadratic_resultant(
    linear: sp.Expr,
    quadratic: sp.Expr,
) -> sp.Expr:
    """Binary resultant of a linear and a quadratic form."""

    linear_poly = sp.Poly(sp.expand(linear), U, V)
    quadratic_poly = sp.Poly(sp.expand(quadratic), U, V)
    a0 = linear_poly.coeff_monomial(U)
    a1 = linear_poly.coeff_monomial(V)
    q0_coefficient = quadratic_poly.coeff_monomial(U**2)
    q1_coefficient = quadratic_poly.coeff_monomial(U * V)
    q2_coefficient = quadratic_poly.coeff_monomial(V**2)
    return sp.expand(
        a0**2 * q2_coefficient
        - a0 * a1 * q1_coefficient
        + a1**2 * q0_coefficient
    )


projective_matrix = sp.Matrix(
    [[matrix_a, matrix_b], [matrix_c, matrix_d]]
)
matrix_determinant = projective_matrix.det()
transformed_linear = transform_binary_form(linear_form, projective_matrix)
transformed_quadratic = transform_binary_form(
    quadratic_form, projective_matrix
)
transformed_product = transform_binary_form(cubic_form, projective_matrix)
normalizing_denominator = tangent_coefficient(transformed_product)

original_resultant = linear_quadratic_resultant(
    linear_form, quadratic_form
)
raw_transformed_resultant = linear_quadratic_resultant(
    transformed_linear, transformed_quadratic
)
assert_zero(
    raw_transformed_resultant
    - matrix_determinant**2 * original_resultant
)

normalized_linear = (
    normalizing_denominator
    / matrix_determinant**2
    * transformed_linear
)
normalized_quadratic = (
    matrix_determinant**2
    / normalizing_denominator**2
    * transformed_quadratic
)
normalized_product = transformed_product / normalizing_denominator

assert_zero(
    linear_quadratic_resultant(normalized_linear, normalized_quadratic)
    - original_resultant
)
assert_zero(normalized_linear * normalized_quadratic - normalized_product)
assert_zero(tangent_coefficient(normalized_product) - 1)

# The source and target normalization is independent of the chosen GL_2
# scalar lift, so it is intrinsically a PGL_2 transport.
rho = sp.symbols("rho", nonzero=True)
scaled_matrix = rho * projective_matrix
scaled_linear_raw = transform_binary_form(linear_form, scaled_matrix)
scaled_quadratic_raw = transform_binary_form(quadratic_form, scaled_matrix)
scaled_product_raw = transform_binary_form(cubic_form, scaled_matrix)
scaled_determinant = scaled_matrix.det()
scaled_denominator = tangent_coefficient(scaled_product_raw)
scaled_normalized_linear = (
    scaled_denominator / scaled_determinant**2 * scaled_linear_raw
)
scaled_normalized_quadratic = (
    scaled_determinant**2
    / scaled_denominator**2
    * scaled_quadratic_raw
)
assert_zero(scaled_normalized_linear - normalized_linear)
assert_zero(scaled_normalized_quadratic - normalized_quadratic)
assert_zero(
    scaled_product_raw / scaled_denominator - normalized_product
)

# The normalized factorization transport has its own exact cocycle.  With
# the substitution convention F -> F(gx), applying g and then h composes as
# the matrix product g*h.
second_a, second_b, second_c, second_d = sp.symbols(
    "second_a second_b second_c second_d"
)
second_matrix = sp.Matrix(
    [[second_a, second_b], [second_c, second_d]]
)
second_determinant = second_matrix.det()
iterated_product_raw = transform_binary_form(
    normalized_product, second_matrix
)
iterated_denominator = tangent_coefficient(iterated_product_raw)
iterated_product = iterated_product_raw / iterated_denominator
iterated_linear = (
    iterated_denominator
    / second_determinant**2
    * transform_binary_form(normalized_linear, second_matrix)
)
iterated_quadratic = (
    second_determinant**2
    / iterated_denominator**2
    * transform_binary_form(normalized_quadratic, second_matrix)
)

composed_matrix = projective_matrix * second_matrix
composed_determinant = composed_matrix.det()
composed_product_raw = transform_binary_form(cubic_form, composed_matrix)
composed_denominator = tangent_coefficient(composed_product_raw)
composed_product = composed_product_raw / composed_denominator
composed_linear = (
    composed_denominator
    / composed_determinant**2
    * transform_binary_form(linear_form, composed_matrix)
)
composed_quadratic = (
    composed_determinant**2
    / composed_denominator**2
    * transform_binary_form(quadratic_form, composed_matrix)
)
assert_zero(iterated_product - composed_product)
assert_zero(iterated_linear - composed_linear)
assert_zero(iterated_quadratic - composed_quadratic)

# The only constant-denominator projective transformations for the tangent
# coefficient chart are diagonal.  Saturating the three unwanted
# coefficient equations by det(g) gives the exact ideal (b,c).
determinant_inverse = sp.symbols("determinant_inverse")
stabilizer_basis = sp.groebner(
    [
        matrix_a * matrix_b**2,
        matrix_b * (2 * matrix_a * matrix_d + matrix_b * matrix_c),
        matrix_c * matrix_d**2,
        determinant_inverse * matrix_determinant - 1,
    ],
    determinant_inverse,
    matrix_a,
    matrix_b,
    matrix_c,
    matrix_d,
    order="lex",
    domain=sp.QQ,
)
assert stabilizer_basis.reduce(matrix_b)[1] == 0
assert stabilizer_basis.reduce(matrix_c)[1] == 0
assert {
    polynomial.as_expr() for polynomial in stabilizer_basis.polys
} == {
    matrix_a * matrix_d * determinant_inverse - 1,
    matrix_b,
    matrix_c,
}


# ---------------------------------------------------------------------------
# 4. The diagonal stabilizer is exactly the known global cubic Keller torus.
# ---------------------------------------------------------------------------

alpha = sp.symbols("alpha", nonzero=True)
x, y, z = sp.symbols("x y z")
t = 1 + x * y
q = t**2 * z + y**2 * (1 + 3 * t)
target_pi = t * q
target_b = -(y + 3 * x * q) / 2
target_c = x * (5 - 3 * t) - x**3 * z

scaled_source = {
    x: alpha * x,
    y: y / alpha,
    z: z / alpha**2,
}
assert_zero(target_pi.subs(scaled_source, simultaneous=True) - target_pi / alpha**2)
assert_zero(target_b.subs(scaled_source, simultaneous=True) - target_b / alpha)
assert_zero(target_c.subs(scaled_source, simultaneous=True) - alpha * target_c)

independent_pi, independent_b, independent_c = sp.symbols(
    "independent_pi independent_b independent_c"
)
old_relation = (
    independent_pi * S**3
    + independent_b * S**2
    + S
    - independent_c / 2
)
scaled_relation = (
    independent_pi / alpha**2 * S**3
    + independent_b / alpha * S**2
    + S
    - alpha * independent_c / 2
)
assert_zero(scaled_relation.subs(S, alpha * S) - alpha * old_relation)

print("PASS: cubic ordered pairs are full S3 frames")
print("PASS: framed root triples give a descended PGL2 cocycle")
print("PASS: quadratic Tschirnhaus formulas and both boundaries are exact")
print("PASS: normalized factorization transport closes after target localization")
print("PASS: its global projective stabilizer is the cubic scaling torus")
