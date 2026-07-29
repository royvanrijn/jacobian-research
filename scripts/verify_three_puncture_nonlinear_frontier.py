#!/usr/bin/env python3
"""Exact checks for the nonlinear three-puncture A6 completion frontier."""

from __future__ import annotations

import sympy as sp


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
