#!/usr/bin/env python3
"""Exact checks for the nonlinear three-puncture A6 completion frontier."""

from __future__ import annotations

import sympy as sp
from sympy.polys.matrices import DomainMatrix


c, v, r, u, z, w, xi = sp.symbols("c v r u z w xi")
block_variables = (r, u, z, w)

D0 = 1 - r * u
D1 = 1 - (r + c) * v
q = 1 - c * v
R = sp.expand(
    sp.integrate((1 - xi * u) * (1 - (xi + c) * v), (xi, 0, r))
)

# If c=L-1, then differentiation in r fixes L,u,v and is the diagonal
# derivative in the original (s,t)-coordinates.
assert sp.factor(sp.diff(R, r) - D0 * D1) == 0

# Retaining R after replacing u cannot give a Keller completion in any
# number of padded variables: every such determinant lies in (R_r,R_u).
# This rational point is an explicit common zero of those two generators.
retained_witness = {c: 0, v: sp.Rational(3, 2), r: 1, u: 1}
assert sp.diff(R, r).subs(retained_witness) == 0
assert sp.diff(R, u).subs(retained_witness) == 0

# First nonlinear orientation: replace R by A=R+D1*z and expose r.  Any
# padded determinant retaining A and r lies in (A_u,A_z).  Both generators
# vanish on V(r,1-cv), so neither A^5 nor identity padding in A^6 can be
# Keller.
A = sp.expand(R + D1 * z)
rank_drop_substitution = {r: 0, c: 1 / v}
assert sp.factor(sp.diff(A, u).subs(rank_drop_substitution)) == 0
assert sp.factor(sp.diff(A, z).subs(rank_drop_substitution)) == 0

# The minimal polynomial transfer is nevertheless exact.  Its determinant
# is q^3, showing explicitly how the rank-two boundary ledger collapses to
# one residual character.
N = sp.expand(
    c**2 * u * v**2
    - c * r * u * v**2
    - 2 * c * u * v
    - 2 * r**2 * u * v**2
    + r * u * v
    + u
    + 6 * v**2 * z
)
transfer_outputs = sp.Matrix([c, v, A, r, N])
transfer_jacobian = sp.factor(
    transfer_outputs.jacobian((c, v, r, u, z)).det()
)
assert sp.factor(transfer_jacobian - q**3) == 0

# On the selected curve L=0, hence c=-1, u=1/r, v=1/(r-1).
curve_substitution = {c: -1, u: 1 / r, v: 1 / (r - 1)}
assert sp.factor(q.subs(curve_substitution) - r / (r - 1)) == 0

# Symmetric orientation: replacing v by A'=R+D0*z and retaining r has a
# different unavoidable rank-drop divisor.  It meets the selected curve at
# r=3 and therefore adds a fourth puncture rather than preserving the core.
A_symmetric = sp.expand(R + D0 * z)
symmetric_rank_drop = {r: 1 / u, c: -1 / (3 * u)}
assert sp.factor(sp.diff(A_symmetric, v).subs(symmetric_rank_drop)) == 0
assert sp.factor(sp.diff(A_symmetric, z).subs(symmetric_rank_drop)) == 0

M = sp.expand(
    -3 * c * r * u**2 * v
    + 3 * c * u * v
    - 2 * r**2 * u**2 * v
    + r * u * v
    + 6 * u**2 * z
    + v
)
symmetric_outputs = sp.Matrix([c, u, A_symmetric, r, M])
symmetric_jacobian = sp.factor(
    symmetric_outputs.jacobian((c, u, r, v, z)).det()
)
assert sp.factor(symmetric_jacobian - (1 + 3 * c * u)) == 0
assert sp.factor(
    (1 + 3 * c * u).subs(curve_substitution) - (r - 3) / r
) == 0


def cofactor_derivation(first_three: tuple[sp.Expr, ...]) -> list[sp.Expr]:
    """Coefficients of D -> det d(first_three,D)/d(r,u,z,w)."""

    rows = sp.Matrix(first_three).jacobian(block_variables)
    coefficients = []
    for column in range(4):
        other_columns = [index for index in range(4) if index != column]
        cofactor = (-1) ** (3 + column) * rows[:, other_columns].det()
        coefficients.append(sp.expand(cofactor))
    return coefficients


def monomials_through(total_degree: int) -> list[sp.Expr]:
    result = []
    for i in range(total_degree + 1):
        for j in range(total_degree + 1 - i):
            for k in range(total_degree + 1 - i - j):
                for ell in range(total_degree + 1 - i - j - k):
                    result.append(r**i * u**j * z**k * w**ell)
    return result


# Sanity-check the cofactor orientation on the identity block.
identity_derivation = cofactor_derivation((r, u, z))
assert sp.expand(
    sum(
        identity_derivation[index] * sp.diff(w, block_variables[index])
        for index in range(4)
    )
    - 1
) == 0

# A bounded A^6 screen.  It couples both modification variables into the
# primitive, replaces the exposed r-coordinate by five elementary nonlinear
# choices, and tries eight transverse coordinate skeletons.  For each of the
# 80 cases, no fourth coordinate of total block degree <=3 over Q(c,v) has
# determinant one.  This is a finite ansatz result, not a general A^6 no-go.
primitives = (
    sp.expand(R + D1 * z + D0 * w),
    sp.expand(R + D1 * z + D0 * w + z * w),
)
second_coordinates = (
    r + z,
    r + w,
    r + z + w,
    r + u * z,
    r + u * w,
)
third_coordinates = (
    u,
    z,
    w,
    u + z,
    u + w,
    z + w,
    D0 + z,
    D1 + w,
)

tested_skeletons = 0
for primitive in primitives:
    for second in second_coordinates:
        for third in third_coordinates:
            tested_skeletons += 1
            derivation = cofactor_derivation((primitive, second, third))
            found_slice = False
            for degree in range(1, 4):
                monomials = monomials_through(degree)
                coefficients = sp.symbols(f"a0:{len(monomials)}")
                candidate = sum(
                    coefficient * monomial
                    for coefficient, monomial in zip(coefficients, monomials)
                )
                equation = sp.Poly(
                    sp.expand(
                        sum(
                            derivation[index]
                            * sp.diff(candidate, block_variables[index])
                            for index in range(4)
                        )
                        - 1
                    ),
                    *block_variables,
                )
                solution = sp.linsolve(equation.coeffs(), coefficients)
                if solution != sp.EmptySet:
                    found_slice = True
                    break
            assert not found_slice, (primitive, second, third)

assert tested_skeletons == 80


# The proposed affine-coupling moonshot can be screened uniformly when its
# two remaining block outputs are affine-linear.  Write
#
#   Rtilde = R + D1*z + D0*w + z*w*H,
#   Btilde = r + z*P + w*Q + z*w*S,
#
# where P,Q are arbitrary affine functions of r,u over K=Q(c,v).  On the
# slice z=w=0, H and S disappear even if they have arbitrary degree.  The
# gradients of two arbitrary affine transverse outputs C,D enter only
# through their six Pluecker coordinates.  In four variables, the single
# Pluecker quadric is necessary and sufficient for such a two-plane.
p0, p1, p2, q0, q1, q2 = sp.symbols("p0 p1 p2 q0 q1 q2")
P_affine = p0 + p1 * r + p2 * u
Q_affine = q0 + q1 * r + q2 * u

A_zero_gradient = (
    sp.diff(R, r),
    sp.diff(R, u),
    D1,
    D0,
)
B_zero_gradient = (
    1,
    0,
    P_affine,
    Q_affine,
)

p12, p13, p14, p23, p24, p34 = sp.symbols(
    "p12 p13 p14 p23 p24 p34"
)
pluecker_coordinates = {
    (0, 1): p12,
    (0, 2): p13,
    (0, 3): p14,
    (1, 2): p23,
    (1, 3): p24,
    (2, 3): p34,
}

zero_slice_determinant = 0
for first_column in range(4):
    for second_column in range(first_column + 1, 4):
        complementary_columns = [
            column
            for column in range(4)
            if column not in (first_column, second_column)
        ]
        permutation = [
            first_column,
            second_column,
            *complementary_columns,
        ]
        inversions = sum(
            permutation[left] > permutation[right]
            for left in range(4)
            for right in range(left + 1, 4)
        )
        first_minor = (
            A_zero_gradient[first_column]
            * B_zero_gradient[second_column]
            - A_zero_gradient[second_column]
            * B_zero_gradient[first_column]
        )
        zero_slice_determinant += (
            (-1) ** inversions
            * first_minor
            * pluecker_coordinates[tuple(complementary_columns)]
        )

zero_slice_polynomial = sp.Poly(
    sp.expand(zero_slice_determinant - 1),
    r,
    u,
)
pluecker_relation = p12 * p34 - p13 * p24 + p14 * p23
affine_completion_equations = [
    *zero_slice_polynomial.coeffs(),
    pluecker_relation,
]
affine_completion_unknowns = (
    p0,
    p1,
    p2,
    q0,
    q1,
    q2,
    p12,
    p13,
    p14,
    p23,
    p24,
    p34,
)
assert len(affine_completion_equations) == 12
affine_completion_basis = sp.groebner(
    affine_completion_equations,
    *affine_completion_unknowns,
    order="grevlex",
    domain=sp.QQ.frac_field(c, v),
)
assert len(affine_completion_basis.polys) == 1
assert affine_completion_basis.contains(sp.Integer(1))

# The next zero-slice screen lets P,Q remain arbitrary affine functions and
# lets the fourth output be a completely general polynomial of block degree
# at most two.  It tests the same eight transverse skeletons used above.
# H,S again disappear on z=w=0.  Every resulting coefficient ideal is the
# unit ideal over Q(c,v).  This is an exact eight-row screen, not a
# classification of arbitrary quadratic third outputs.
quadratic_monomials = [
    monomial
    for monomial in monomials_through(2)
    if monomial != 1
]
quadratic_coefficients = sp.symbols(
    f"d0:{len(quadratic_monomials)}"
)
quadratic_fourth_output = sum(
    coefficient * monomial
    for coefficient, monomial in zip(
        quadratic_coefficients,
        quadratic_monomials,
        strict=True,
    )
)
quadratic_fourth_gradient = tuple(
    sp.diff(quadratic_fourth_output, variable).subs({z: 0, w: 0})
    for variable in block_variables
)

quadratic_zero_slice_skeletons = 0
for third_coordinate in third_coordinates:
    quadratic_zero_slice_skeletons += 1
    third_gradient = tuple(
        sp.diff(third_coordinate, variable).subs({z: 0, w: 0})
        for variable in block_variables
    )
    determinant_equation = sp.expand(
        sp.Matrix(
            [
                A_zero_gradient,
                B_zero_gradient,
                third_gradient,
                quadratic_fourth_gradient,
            ]
        ).det(method="domain-ge")
        - 1
    )
    coefficient_equations = sp.Poly(
        determinant_equation,
        r,
        u,
    ).coeffs()
    coefficient_basis = sp.groebner(
        coefficient_equations,
        p0,
        p1,
        p2,
        q0,
        q1,
        q2,
        *quadratic_coefficients,
        order="grevlex",
        domain=sp.QQ.frac_field(c, v),
    )
    assert len(coefficient_basis.polys) == 1
    assert coefficient_basis.contains(sp.Integer(1))

assert quadratic_zero_slice_skeletons == 8

# The affine directions with nonzero r coefficient can be normalized to
#
#   C = r + g*u + a*z + b*w.
#
# Their general quadratic fourth-output equation is linear in the fourteen
# quadratic coefficients.  Exact row reduction over successive
# rational-function fields gives a complete exceptional-pivot tree:
#
#   generic -> p1=0 -> q1=0 -> p2=0 -> q2=0 -> g=0
#           -> p0=a.
#
# The coefficient matrix and augmented matrix have ranks 8 and 9 on each of
# the first four opens.  Thereafter their ranks are 6 and 7, including the
# terminal branch.  Thus every branch is inconsistent.
g, a, b = sp.symbols("g a b")


def quadratic_zero_slice_system(
    third_gradient: tuple[sp.Expr | int, ...],
) -> tuple[sp.Matrix, sp.Matrix]:
    """Linear system for a general quadratic fourth output on z=w=0."""

    equation = sp.expand(
        sp.Matrix(
            [
                A_zero_gradient,
                B_zero_gradient,
                third_gradient,
                quadratic_fourth_gradient,
            ]
        ).det(method="domain-ge")
        - 1
    )
    coefficient_equations = sp.Poly(equation, r, u).coeffs()
    return sp.linear_eq_to_matrix(
        coefficient_equations,
        quadratic_coefficients,
    )


nonzero_r_affine_matrix, nonzero_r_affine_rhs = (
    quadratic_zero_slice_system((1, g, a, b))
)


def rational_rref_audit(
    matrix: sp.Matrix,
    substitutions: dict[sp.Symbol, sp.Expr | int],
    parameters: tuple[sp.Symbol, ...],
) -> tuple[int, set[sp.Expr]]:
    """Return rank and visible pivot denominators over the remaining field."""

    remaining_parameters = tuple(
        parameter
        for parameter in parameters
        if parameter not in substitutions
    )
    coefficient_field = sp.QQ.frac_field(c, v, *remaining_parameters)
    reduced, pivots = DomainMatrix.from_Matrix(
        matrix.subs(substitutions),
        fmt="sparse",
    ).convert_to(coefficient_field).rref()
    denominators: set[sp.Expr] = set()
    for row in reduced.rep.values():
        for value in row.values():
            if str(value.denom) != "1":
                denominators.add(sp.factor(value.denom.as_expr()))
    return len(pivots), denominators


def verify_pivot_tree(
    coefficient_matrix: sp.Matrix,
    rhs: sp.Matrix,
    parameters: tuple[sp.Symbol, ...],
    tree: tuple[
        tuple[dict[sp.Symbol, sp.Expr | int], int, int, set[sp.Expr]],
        ...,
    ],
) -> None:
    """Check every rank and every exceptional denominator in a pivot tree."""

    for substitutions, expected_rank, expected_augmented_rank, expected_poles in (
        tree
    ):
        actual_rank, actual_poles = rational_rref_audit(
            coefficient_matrix,
            substitutions,
            parameters,
        )
        augmented_rank, augmented_poles = rational_rref_audit(
            coefficient_matrix.row_join(rhs),
            substitutions,
            parameters,
        )
        assert actual_rank == expected_rank
        assert augmented_rank == expected_augmented_rank
        assert actual_poles == expected_poles
        assert augmented_poles == expected_poles


nonzero_r_parameters = (g, a, b, p0, p1, p2, q0, q1, q2)
nonzero_r_pivot_tree = (
    ({}, 8, 9, {p1}),
    ({p1: 0}, 8, 9, {q1}),
    ({p1: 0, q1: 0}, 8, 9, {p2}),
    ({p1: 0, q1: 0, p2: 0}, 8, 9, {q2}),
    ({p1: 0, q1: 0, p2: 0, q2: 0}, 6, 7, {g}),
    ({p1: 0, q1: 0, p2: 0, q2: 0, g: 0}, 6, 7, {a - p0}),
    (
        {p1: 0, q1: 0, p2: 0, q2: 0, g: 0, p0: a},
        6,
        7,
        set(),
    ),
)
verify_pivot_tree(
    nonzero_r_affine_matrix,
    nonzero_r_affine_rhs,
    nonzero_r_parameters,
    nonzero_r_pivot_tree,
)

# It remains to cover affine gradients with zero r coefficient.  Projective
# normalization gives three charts, according to the first nonzero entry:
#
#   (0,1,a,b),  (0,0,1,b),  (0,0,0,1).
#
# Their shorter exact pivot trees close the entire r-free boundary.
r_free_u_matrix, r_free_u_rhs = quadratic_zero_slice_system((0, 1, a, b))
r_free_u_parameters = (a, b, p0, p1, p2, q0, q1, q2)
r_free_u_pivot_tree = (
    ({}, 8, 9, {p1}),
    ({p1: 0}, 8, 9, {q1}),
    ({p1: 0, q1: 0}, 8, 9, set()),
)
verify_pivot_tree(
    r_free_u_matrix,
    r_free_u_rhs,
    r_free_u_parameters,
    r_free_u_pivot_tree,
)

r_free_z_matrix, r_free_z_rhs = quadratic_zero_slice_system((0, 0, 1, b))
r_free_z_parameters = (b, p0, p1, p2, q0, q1, q2)
r_free_z_pivot_tree = (
    ({}, 7, 8, {b * p1 - q1}),
    ({q1: b * p1}, 7, 8, set()),
)
verify_pivot_tree(
    r_free_z_matrix,
    r_free_z_rhs,
    r_free_z_parameters,
    r_free_z_pivot_tree,
)

r_free_w_matrix, r_free_w_rhs = quadratic_zero_slice_system((0, 0, 0, 1))
r_free_w_parameters = (p0, p1, p2, q0, q1, q2)
r_free_w_pivot_tree = (
    ({}, 7, 8, {p1}),
    ({p1: 0}, 7, 8, set()),
)
verify_pivot_tree(
    r_free_w_matrix,
    r_free_w_rhs,
    r_free_w_parameters,
    r_free_w_pivot_tree,
)

# The first simultaneous-quadratic chart retains B=r, so P=Q=0 on the
# zero slice, and lets C,D both be completely general of block degree at
# most two.  Solving linearly for D gives another complete pivot tree in
# the visible coefficients of C.  Pure z,w quadratics do not occur because
# their gradients vanish on z=w=0.
quadratic_third_coefficients = sp.symbols(
    f"e0:{len(quadratic_monomials)}"
)
quadratic_third_output = sum(
    coefficient * monomial
    for coefficient, monomial in zip(
        quadratic_third_coefficients,
        quadratic_monomials,
        strict=True,
    )
)
quadratic_third_gradient = tuple(
    sp.diff(quadratic_third_output, variable).subs({z: 0, w: 0})
    for variable in block_variables
)
exposed_r_equation = sp.expand(
    sp.Matrix(
        [
            A_zero_gradient,
            (1, 0, 0, 0),
            quadratic_third_gradient,
            quadratic_fourth_gradient,
        ]
    ).det(method="domain-ge")
    - 1
)
exposed_r_matrix, exposed_r_rhs = sp.linear_eq_to_matrix(
    sp.Poly(exposed_r_equation, r, u).coeffs(),
    quadratic_coefficients,
)
exposed_r_visible_indices = (0, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13)
exposed_r_parameters = tuple(
    quadratic_third_coefficients[index]
    for index in exposed_r_visible_indices
)
e = quadratic_third_coefficients
exposed_r_quadratic_tree = (
    ({}, 8, 9, {e[12]}),
    ({e[12]: 0}, 8, 9, {e[11]}),
    ({e[12]: 0, e[11]: 0}, 8, 9, {e[10]}),
    ({e[12]: 0, e[11]: 0, e[10]: 0}, 8, 9, {e[8]}),
    ({e[12]: 0, e[11]: 0, e[10]: 0, e[8]: 0}, 8, 9, {e[7]}),
    (
        {e[12]: 0, e[11]: 0, e[10]: 0, e[8]: 0, e[7]: 0},
        8,
        9,
        {e[6]},
    ),
    (
        {
            e[12]: 0,
            e[11]: 0,
            e[10]: 0,
            e[8]: 0,
            e[7]: 0,
            e[6]: 0,
        },
        6,
        7,
        {e[5]},
    ),
    (
        {
            e[12]: 0,
            e[11]: 0,
            e[10]: 0,
            e[8]: 0,
            e[7]: 0,
            e[6]: 0,
            e[5]: 0,
        },
        6,
        7,
        {e[2]},
    ),
    (
        {
            e[12]: 0,
            e[11]: 0,
            e[10]: 0,
            e[8]: 0,
            e[7]: 0,
            e[6]: 0,
            e[5]: 0,
            e[2]: 0,
        },
        6,
        7,
        set(),
    ),
)
verify_pivot_tree(
    exposed_r_matrix,
    exposed_r_rhs,
    exposed_r_parameters,
    exposed_r_quadratic_tree,
)

# Degree three is the sharp zero-slice threshold for the exposed-r chart.
# With C=w, the displayed rational cubic D has determinant one on z=w=0.
# Clearing q^3 gives a polynomial numerator, but its full H=0 determinant
# is q^3+6*r*v^2*w rather than a unit.  In fact the cofactor derivation for
# any fourth output vanishes where A_u=A_z=0.
exposed_r_primitive = sp.expand(R + D1 * z + D0 * w)
exposed_r_cubic_numerator = sp.expand(
    u * (-q**2 - v * r * q + 2 * v**2 * r**2) - 6 * v**2 * z
)
exposed_r_rational_cubic = exposed_r_cubic_numerator / q**3
exposed_r_cubic_outputs = sp.Matrix(
    [exposed_r_primitive, r, w, exposed_r_cubic_numerator]
)
exposed_r_cubic_determinant = sp.factor(
    exposed_r_cubic_outputs.jacobian(block_variables).det()
)
assert sp.factor(
    exposed_r_cubic_determinant - (q**3 + 6 * r * v**2 * w)
) == 0
assert sp.factor(
    sp.Matrix(
        [exposed_r_primitive, r, w, exposed_r_rational_cubic]
    ).jacobian(block_variables).det().subs({z: 0, w: 0})
    - 1
) == 0

exposed_r_derivation = cofactor_derivation((exposed_r_primitive, r, w))
assert tuple(
    sp.factor(coefficient.subs({r: 0, z: 0, w: 0}))
    for coefficient in exposed_r_derivation
) == (0, -q, 0, 0)
rank_drop_r = q / v
rank_drop_w = sp.factor(sp.diff(R, u).subs(r, rank_drop_r) / rank_drop_r)
assert all(
    sp.factor(
        coefficient.subs(
            {
                r: rank_drop_r,
                w: rank_drop_w,
            }
        )
    )
    == 0
    for coefficient in exposed_r_derivation
)


print(
    "PASS: retained-primitive and one-sided nonlinear rank-drop gates; "
    "exact rank-one transfer determinants"
)
print(
    "PASS: 80 coupled A^6 coordinate skeletons have no degree-<=3 "
    "polynomial Jacobian slice over Q(c,v)"
)
print(
    "PASS: the proposed affine P,Q coupling has no completion with two "
    "affine transverse outputs over Q(c,v)"
)
print(
    "PASS: all 8 transverse skeletons reject a general degree-<=2 fourth "
    "output with arbitrary affine P,Q"
)
print(
    "PASS: every nonconstant affine C rejects a general degree-<=2 fourth "
    "output on four complete projective pivot trees"
)
print(
    "PASS: exposed B=r rejects two general quadratics; its first rational "
    "cubic zero-slice survivor has polynomial numerator determinant "
    "q^3+6*r*v^2*w and a full rank-drop locus"
)
