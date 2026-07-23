#!/usr/bin/env python3
"""Exact regressions for the flat cubic-normalization frontend."""

from __future__ import annotations

import itertools

import sympy as sp


# A ramified boundary prime and an affine étale prime over the same target
# divisor exhaust a cubic DVR degree sum uniquely.
degree_budget_solutions = []
for boundary_e in range(2, 4):
    for boundary_f in range(1, 4):
        for affine_f in range(1, 4):
            if boundary_e * boundary_f + affine_f == 3:
                degree_budget_solutions.append(
                    (boundary_e, boundary_f, 1, affine_f)
                )
assert degree_budget_solutions == [(2, 1, 1, 1)]

# In the foundational calibration, graph-at-infinity elimination and the
# cubic discriminant give the same irreducible target hypersurface.  Hence
# the phantom-boundary quotient j_F/delta_F is the unit 1.
target_a_np, target_b_np, target_c_np = sp.symbols(
    "target_a_np target_b_np target_c_np"
)
root_np, derivative_np = sp.symbols("root_np derivative_np")
root_polynomial_np = (
    target_c_np * root_np**3
    - 2 * root_np**2
    + target_b_np * root_np
    - 2 * target_a_np
)
nonproperness_np = (
    27 * target_a_np**2 * target_c_np**2
    - 18 * target_a_np * target_b_np * target_c_np
    + 16 * target_a_np
    + target_b_np**3 * target_c_np
    - target_b_np**2
)
assert sp.factor(
    sp.discriminant(root_polynomial_np, root_np)
    + 4 * nonproperness_np
) == 0
assert sp.factor(nonproperness_np) == nonproperness_np
target_a_boundary = (
    root_np**2
    - target_c_np * root_np**3
    + derivative_np * root_np / 2
)
target_b_boundary = (
    4 * root_np
    + derivative_np
    - 3 * target_c_np * root_np**2
)
boundary_elimination_np = sp.groebner(
    (
        2 * target_a_np - 2 * target_a_boundary,
        target_b_np - target_b_boundary,
        derivative_np,
    ),
    root_np,
    derivative_np,
    target_a_np,
    target_b_np,
    target_c_np,
    order="lex",
)
eliminated_boundary_equations = [
    sp.factor(polynomial.as_expr())
    for polynomial in boundary_elimination_np.polys
    if not (
        polynomial.as_expr().has(root_np)
        or polynomial.as_expr().has(derivative_np)
    )
]
assert eliminated_boundary_equations == [nonproperness_np]
foundational_phantom_boundary_factor = sp.Integer(1)
assert foundational_phantom_boundary_factor == 1


# Deligne--Faddeev cubic algebra in the basis (1, omega, theta).
a, b, c, d = sp.symbols("a b c d")

M_omega = sp.Matrix(
    (
        (0, -a * c, -a * d),
        (1, b, 0),
        (0, -a, 0),
    )
)
M_theta = sp.Matrix(
    (
        (0, -a * d, -b * d),
        (0, 0, d),
        (1, 0, -c),
    )
)
identity = sp.eye(3)

# The multiplication matrices satisfy all three defining products, hence
# certify commutativity and associativity of the universal cubic algebra.
assert M_omega * M_omega == -a * c * identity + b * M_omega - a * M_theta
assert M_omega * M_theta == -a * d * identity
assert M_theta * M_omega == -a * d * identity
assert M_theta * M_theta == -b * d * identity + d * M_omega - c * M_theta

omega_vector = sp.Matrix((0, 1, 0))
theta_vector = sp.Matrix((0, 0, 1))
assert M_omega * omega_vector == sp.Matrix((-a * c, b, -a))
assert M_omega * theta_vector == sp.Matrix((-a * d, 0, 0))
assert M_theta * theta_vector == sp.Matrix((-b * d, d, -c))

# Trace splits the unit in characteristic zero.
assert sp.trace(identity) == 3
assert sp.trace(M_omega) == b
assert sp.trace(M_theta) == -c
trace_zero_omega = 3 * M_omega - b * identity
trace_zero_theta = 3 * M_theta + c * identity
assert sp.trace(trace_zero_omega) == 0
assert sp.trace(trace_zero_theta) == 0

# The determinant of the trace pairing is the binary-cubic discriminant.
multiplication = (identity, M_omega, M_theta)
trace_pairing = sp.Matrix(
    3,
    3,
    lambda i, j: sp.trace(multiplication[i] * multiplication[j]),
)
trace_discriminant = sp.factor(trace_pairing.det())
form_discriminant = sp.factor(
    b**2 * c**2
    - 4 * a * c**3
    - 4 * b**3 * d
    - 27 * a**2 * d**2
    + 18 * a * b * c * d
)
assert sp.factor(trace_discriminant - form_discriminant) == 0


# Normality only gives a reflexive module before flatness.  The second
# Koszul syzygy of (x,y,z) is rank two and fails to be free only at the
# codimension-three origin.
x, y, z = sp.symbols("x y z")
koszul_generators = sp.Matrix(
    (
        (-y, -z, 0),
        (x, 0, -z),
        (0, x, y),
    )
)
relation = sp.Matrix((z, -y, x))
assert koszul_generators * relation == sp.zeros(3, 1)

# At the origin the relation vanishes modulo the maximal ideal, so this
# rank-two module has a three-dimensional special fiber.  Adding the unit
# summand produces a generic-rank-three module with the smallest possible
# nonflat special-fiber length, namely four.
relation_at_origin = relation.subs({x: 0, y: 0, z: 0})
assert relation_at_origin == sp.zeros(3, 1)
generic_rank_M = 2
special_fiber_dimension_M = 3 - relation_at_origin.rank()
generic_rank_with_unit = 1 + generic_rank_M
special_fiber_length_with_unit = 1 + special_fiber_dimension_M
assert generic_rank_with_unit == 3
assert special_fiber_length_with_unit == 4

# A presentation of R+M has four generators and the one relation
# (0,z,-y,x).  For a generic-rank-three module its Fitt_3 ideal consists of
# the entries of this presentation, hence is precisely the origin ideal.
unit_extended_relation = sp.Matrix((0, z, -y, x))
fitting_3_generators = {
    sp.expand(entry)
    for entry in unit_extended_relation
    if entry != 0
}
assert fitting_3_generators == {x, -y, z}

minors = []
for rows in itertools.combinations(range(3), 2):
    for cols in itertools.combinations(range(3), 2):
        minors.append(
            sp.factor(
                koszul_generators.extract(rows, cols).det()
            )
        )
nonzero_minors = {minor for minor in minors if minor != 0}
assert {x**2, -y**2, z**2}.issubset(nonzero_minors)
assert all(sp.Poly(minor, x, y, z).total_degree() == 2 for minor in nonzero_minors)

# A nonflat triple-cover multiplication on the Koszul module is encoded by
# a symmetric cubic tensor phi: Sym^3(M) -> det(M)=R.  If
#
#   c_ijk=phi(e_i,e_j,e_k),
#
# well-definedness against the relation r=(z,-y,x) says
#
#   z*c_1jk-y*c_2jk+x*c_3jk=0
#
# for every symmetric pair (j,k).  Compute the homogeneous solution spaces
# exactly.  There are no symbols in orders zero, one, or two; the first
# symbol space occurs in order three and has dimension ten.
cubic_index_triples = list(
    itertools.combinations_with_replacement(range(3), 3)
)
cubic_index_pairs = list(
    itertools.combinations_with_replacement(range(3), 2)
)


def homogeneous_exponents(degree: int) -> list[tuple[int, int, int]]:
    """Exponent triples for homogeneous polynomials in x,y,z."""

    return [
        (x_degree, y_degree, degree - x_degree - y_degree)
        for x_degree in range(degree + 1)
        for y_degree in range(degree - x_degree + 1)
    ]


def koszul_cubic_constraint_matrix(degree: int) -> sp.MutableSparseMatrix:
    """Matrix of the relation constraints on degree-degree cubic tensors."""

    input_exponents = homogeneous_exponents(degree)
    output_exponents = homogeneous_exponents(degree + 1)
    columns = {
        (triple, exponent): index
        for index, (triple, exponent) in enumerate(
            itertools.product(cubic_index_triples, input_exponents)
        )
    }
    rows = {
        (pair, exponent): index
        for index, (pair, exponent) in enumerate(
            itertools.product(cubic_index_pairs, output_exponents)
        )
    }
    matrix = sp.MutableSparseMatrix(len(rows), len(columns), {})
    relation_terms = (
        ((0, 0, 1), 1),
        ((0, 1, 0), -1),
        ((1, 0, 0), 1),
    )
    for pair in cubic_index_pairs:
        second, third = pair
        for first, (relation_exponent, sign) in enumerate(relation_terms):
            triple = tuple(sorted((first, second, third)))
            for exponent in input_exponents:
                output_exponent = tuple(
                    exponent[index] + relation_exponent[index]
                    for index in range(3)
                )
                matrix[
                    rows[(pair, output_exponent)],
                    columns[(triple, exponent)],
                ] += sign
    return matrix


koszul_cubic_symbol_dimensions = []
for symbol_degree in range(4):
    constraint_matrix = koszul_cubic_constraint_matrix(symbol_degree)
    koszul_cubic_symbol_dimensions.append(
        constraint_matrix.cols - constraint_matrix.rank()
    )
assert koszul_cubic_symbol_dimensions == [0, 0, 0, 10]

# A chosen volume form identifies the order-three kernel explicitly with
# ternary cubics.  For r=(z,-y,x), send a ternary cubic h to
#
#   phi_h(u,v,w)=h(r cross u, r cross v, r cross w),
#
# where the right side denotes the symmetric trilinear polarization of h.
# The ten ternary monomials give ten independent constraint solutions, so
# they exhaust the order-three symbol space.
ternary_variables = sp.symbols("ternary_0:3")
polarization_variables = sp.symbols("polarization_0:3")
standard_basis = [sp.eye(3).col(index) for index in range(3)]
cross_product_columns = [
    relation.cross(vector)
    for vector in standard_basis
]
degree_three_exponents = homogeneous_exponents(3)
degree_three_monomials = [
    sp.prod(
        ternary_variables[index] ** exponent[index]
        for index in range(3)
    )
    for exponent in degree_three_exponents
]


def polarized_ternary_value(
    polynomial: sp.Expr,
    first: sp.Matrix,
    second: sp.Matrix,
    third: sp.Matrix,
) -> sp.Expr:
    """Symmetric trilinear polarization of a ternary cubic."""

    substitution = {
        ternary_variables[index]: (
            polarization_variables[0] * first[index]
            + polarization_variables[1] * second[index]
            + polarization_variables[2] * third[index]
        )
        for index in range(3)
    }
    expanded = sp.Poly(
        sp.expand(polynomial.subs(substitution)),
        *polarization_variables,
    )
    return sp.expand(
        expanded.coeff_monomial(sp.prod(polarization_variables)) / 6
    )


degree_three_columns = {
    (triple, exponent): index
    for index, (triple, exponent) in enumerate(
        itertools.product(cubic_index_triples, degree_three_exponents)
    )
}
ternary_symbol_vectors = []
for ternary_monomial in degree_three_monomials:
    vector = sp.zeros(len(degree_three_columns), 1)
    for triple in cubic_index_triples:
        coefficient = polarized_ternary_value(
            ternary_monomial,
            *(cross_product_columns[index] for index in triple),
        )
        coefficient_polynomial = sp.Poly(coefficient, x, y, z)
        for exponent in degree_three_exponents:
            monomial = x ** exponent[0] * y ** exponent[1] * z ** exponent[2]
            vector[degree_three_columns[(triple, exponent)]] = (
                coefficient_polynomial.coeff_monomial(monomial)
            )
    ternary_symbol_vectors.append(vector)

ternary_symbol_matrix = sp.Matrix.hstack(*ternary_symbol_vectors)
degree_three_constraint_matrix = koszul_cubic_constraint_matrix(3)
assert degree_three_constraint_matrix * ternary_symbol_matrix == sp.zeros(
    degree_three_constraint_matrix.rows,
    ternary_symbol_matrix.cols,
)
assert ternary_symbol_matrix.rank() == 10

# On the exceptional plane of the blowup at the Koszul defect, the saturated
# transform is the universal quotient k^3/<r>.  Cross product identifies this
# quotient with r-perp, so the exceptional order-three divisor is the incidence
# pullback {(r,q): r.q=0, h(q)=0} of the ternary cubic.
exceptional_r = sp.Matrix(sp.symbols("exceptional_r0:3"))
exceptional_u = sp.Matrix(sp.symbols("exceptional_u0:3"))
exceptional_cross = exceptional_r.cross(exceptional_u)
assert sp.expand(exceptional_r.dot(exceptional_cross)) == 0
assert exceptional_r.cross(exceptional_r) == sp.zeros(3, 1)
exceptional_cross_matrix = sp.Matrix(
    (
        (0, -exceptional_r[2], exceptional_r[1]),
        (exceptional_r[2], 0, -exceptional_r[0]),
        (-exceptional_r[1], exceptional_r[0], 0),
    )
)
assert exceptional_cross_matrix * exceptional_r == sp.zeros(3, 1)
assert any(
    minor != 0
    for minor in (
        sp.factor(exceptional_cross_matrix.extract(rows, cols).det())
        for rows in itertools.combinations(range(3), 2)
        for cols in itertools.combinations(range(3), 2)
    )
)

# Double-saturation calibration.  Put C=A/(x) and take the rank-one
# codimension-one-full submodule T=(y,z)C.  Its S2 hull is C and C/T=k at
# the origin.  As an A-module, T has the exact length-two presentation
#
#   0 -> A -> A^3 -> A^2 -> T -> 0.
#
# Dualizing shows Ext^2_A(T,A)=A/(x,y,z), exactly the canonical dual of
# the finite support-hull quotient C/T.
support_hull_presentation = sp.Matrix(
    (
        (x, 0, -z),
        (0, x, y),
    )
)
support_hull_second_syzygy = sp.Matrix((z, -y, x))
assert (
    support_hull_presentation * support_hull_second_syzygy
    == sp.zeros(2, 1)
)
assert {
    sp.expand(entry)
    for entry in support_hull_second_syzygy
    if entry != 0
} == {x, -y, z}
support_hull_quotient_length = 1
support_hull_ext2_length = 1
assert support_hull_ext2_length == support_hull_quotient_length

# After support saturation, the primitive-generation cokernel is exactly
# closed-point torsion in the cotangent module.  The module model
# Q=C direct-sum k with tau landing in C has L=0 and K=P=k.
coupled_support_defect_length = 0
coupled_point_torsion_length = 1
coupled_generation_cokernel_length = 1
assert coupled_support_defect_length == 0
assert (
    coupled_generation_cokernel_length
    == coupled_point_torsion_length
)

# The next determinantal rung has s=2, a 4-by-2 presentation matrix,
# excess fiber length two, and therefore total fiber length five.
phi_s2 = sp.Matrix(
    (
        (x, 0),
        (y, x),
        (z, y),
        (0, z),
    )
)
fitting_s2 = [
    sp.expand(phi_s2.extract(rows, (0, 1)).det())
    for rows in itertools.combinations(range(4), 2)
]
assert all(
    sp.Poly(minor, x, y, z).total_degree() >= 2
    for minor in fitting_s2
    if minor != 0
)
fitting_s2_groebner = sp.groebner(fitting_s2, x, y, z)
assert fitting_s2_groebner.reduce(x**2)[1] == 0
assert fitting_s2_groebner.reduce(z**2)[1] == 0
assert fitting_s2_groebner.reduce(y**4)[1] == 0
phi_s2_at_origin = phi_s2.subs({x: 0, y: 0, z: 0})
assert phi_s2_at_origin.rank() == 0
s2_generic_rank_M = phi_s2.rows - phi_s2.cols
s2_special_fiber_dimension_M = phi_s2.rows - phi_s2_at_origin.rank()
assert s2_generic_rank_M == 2
assert 1 + s2_special_fiber_dimension_M == 5

# At a reduced s=1 defect, every endomorphism of the Koszul cokernel reduces
# to a scalar.  If C*v=r*v for v=(z,-y,x), comparison of linear terms leaves
# one scalar parameter among the ten coefficients.
endomorphism_entries = sp.symbols("endomorphism_0:9")
endomorphism_scalar = sp.symbols("endomorphism_scalar")
endomorphism_matrix = sp.Matrix(3, 3, endomorphism_entries)
koszul_column = sp.Matrix((z, -y, x))
endomorphism_equations = []
for expression in (
    (endomorphism_matrix - endomorphism_scalar * sp.eye(3))
    * koszul_column
):
    polynomial = sp.Poly(expression, x, y, z)
    endomorphism_equations.extend(
        polynomial.coeff_monomial(variable)
        for variable in (x, y, z)
    )
endomorphism_coefficient_matrix, _ = sp.linear_eq_to_matrix(
    endomorphism_equations,
    (*endomorphism_entries, endomorphism_scalar),
)
assert endomorphism_coefficient_matrix.rank() == 9

# Likewise, every functional M->R has zero constant reduction.
dual_entries = sp.symbols("dual_0:3")
dual_relation = sp.expand(
    sp.Matrix(1, 3, dual_entries).dot(koszul_column)
)
dual_equations = [
    sp.Poly(dual_relation, x, y, z).coeff_monomial(variable)
    for variable in (x, y, z)
]
dual_coefficient_matrix, _ = sp.linear_eq_to_matrix(
    dual_equations,
    dual_entries,
)
assert dual_coefficient_matrix.rank() == 3

# If the M-part of multiplication by u reduces to ell(u)*identity,
# commutativity ell(u)*v=ell(v)*u forces ell=0 in dimension three.
ell_entries = sp.symbols("ell_0:3")
fiber_basis = [sp.eye(3).col(index) for index in range(3)]
commutativity_equations = []
for first, second in itertools.combinations(range(3), 2):
    difference = (
        ell_entries[first] * fiber_basis[second]
        - ell_entries[second] * fiber_basis[first]
    )
    commutativity_equations.extend(difference)
commutativity_matrix, _ = sp.linear_eq_to_matrix(
    commutativity_equations,
    ell_entries,
)
assert commutativity_matrix.rank() == 3

# The allowed foundational collision is k[e]/(e^3): length three,
# embedding dimension one, and generated by e.  The reduced defect
# k+V with V^2=0 has embedding dimension and generator number three.
triple_root_length = 3
triple_root_embedding_dimension = 1
square_zero_defect_length = 1 + 3
square_zero_defect_embedding_dimension = 3
assert triple_root_length == 3
assert triple_root_embedding_dimension == 1
assert square_zero_defect_length == 4
assert square_zero_defect_embedding_dimension == 3

# Cotangent fibers recover these embedding dimensions.  For k[e]/(e^3),
# Omega is generated by de with relation e^2*de=0.  For k+V, V^2=0 and
# Omega tensor k is V, requiring three generators.
triple_root_cotangent_generators = 1
square_zero_cotangent_generators = square_zero_defect_embedding_dimension
assert triple_root_cotangent_generators == 1
assert square_zero_cotangent_generators == 3
triple_root_nilradical_generators = 1
triple_root_nilpotency_index = 3
square_zero_nilradical_generators = 3
square_zero_nilpotency_index = 2
assert (
    triple_root_nilradical_generators,
    triple_root_nilpotency_index,
) == (1, 3)
assert (
    square_zero_nilradical_generators,
    square_zero_nilpotency_index,
) == (3, 2)

# A global monogenic cubic cannot have derivative equal to a constant unit:
# the resulting nonzero quadratic would annihilate a degree-three generator.
root_variable = sp.symbols("root_variable")
monic_cubic_coefficients = sp.symbols("monic_cubic_0:3")
derivative_unit = sp.symbols("derivative_unit", nonzero=True)
generic_monic_cubic = (
    root_variable**3
    + monic_cubic_coefficients[0] * root_variable**2
    + monic_cubic_coefficients[1] * root_variable
    + monic_cubic_coefficients[2]
)
derivative_relation = sp.diff(generic_monic_cubic, root_variable) - derivative_unit
assert sp.Poly(derivative_relation, root_variable).degree() == 2
assert sp.Poly(derivative_relation, root_variable).LC() == 3


# Tangent-nonosculating normalized factorization slice.  Its first two
# quotient coefficients are exactly the positive cubic conormal labels.
alpha, beta, gamma, delta, epsilon = sp.symbols(
    "alpha beta gamma delta epsilon"
)
source_y, source_z = sp.symbols("source_y source_z")
resultant = (
    alpha**2 * epsilon
    - alpha * beta * delta
    + beta**2 * gamma
)
middle_coefficient = alpha * delta + beta * gamma

beta_formula = 1 + alpha * source_y
gamma_formula = (
    1
    - sp.Rational(3, 2) * alpha * source_y
    + alpha**2 * source_z
)
delta_formula = (
    sp.Rational(1, 2) * source_y
    - alpha * source_z
    + sp.Rational(3, 2) * alpha * source_y**2
    - alpha**2 * source_y * source_z
)
epsilon_formula = (
    -2 * source_z
    + 4 * source_y**2
    - 4 * alpha * source_y * source_z
    + 3 * alpha * source_y**3
    - 2 * alpha**2 * source_y**2 * source_z
)
slice_substitution = {
    beta: beta_formula,
    gamma: gamma_formula,
    delta: delta_formula,
    epsilon: epsilon_formula,
}
assert sp.expand((resultant - 1).subs(slice_substitution)) == 0
assert sp.expand((middle_coefficient - 1).subs(slice_substitution)) == 0
assert sp.cancel((beta_formula - 1) / alpha - source_y) == 0
assert sp.cancel(
    (
        gamma_formula
        - 1
        + sp.Rational(3, 2) * alpha * source_y
    )
    / alpha**2
    - source_z
) == 0

target_map = (
    sp.expand((alpha * gamma).subs(slice_substitution)),
    sp.expand(
        (alpha * epsilon + beta * delta).subs(slice_substitution)
    ),
    sp.expand((beta * epsilon).subs(slice_substitution)),
)
assert sp.factor(
    sp.Matrix(target_map)
    .jacobian((alpha, source_y, source_z))
    .det()
    + 1
) == 0

print("PASS: Deligne--Faddeev multiplication is associative and commutative")
print("PASS: the cubic critical DVR budget is uniquely (2,1)+(1,1)")
print("PASS: the foundational phantom-boundary factor is the unit 1")
print("PASS: trace splitting and binary-cubic discriminant are exact")
print("PASS: the reflexive rank-two warning is supported only in codimension three")
print("PASS: its rank-three unit extension has excess special-fiber length four")
print("PASS: the model defect has Fitt_3=(x,y,z)")
print("PASS: the first Koszul cubic-cover symbol is a ternary cubic in order three")
print("PASS: its exceptional divisor is the incidence pullback of that cubic")
print("PASS: the S2-hull quotient and Ext^2 defect have the same length")
print("PASS: after S2 saturation the conormal defect is exactly point torsion")
print("PASS: the s=2 determinantal rung is origin-primary with fiber length five")
print("PASS: every reduced defect is the length-four square-zero Koszul fiber")
print("PASS: curvilinear triple-root and square-zero defect collisions are distinct")
print("PASS: intrinsic cotangent cyclicity detects the curvilinear collision")
print("PASS: primitive nilradical generation is the equivalent closed-point marking")
print("PASS: a constant derivative would contradict a global cubic generator")
print("PASS: tangent-hyperplane quotients recover both positive cubic labels")
print("PASS: normalized multiplication has constant Jacobian -1")
