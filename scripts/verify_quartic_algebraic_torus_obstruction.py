#!/usr/bin/env python3
"""Exact algebra for the intrinsic quartic torus obstruction.

The conceptual faithfulness argument is proved in
``cancellation/NO_ALGEBRAIC_TORUS_EQUIVARIANCE.md``. This checker verifies
the specialized Fitting generator, its scheme-theoretic stabilizer, the
zero tangent algebra, the explicit residual mu_5 left-right symmetries, and
the rational support-minimal representative in the quartic quadratic-gauge
normal form.
"""

from __future__ import annotations

from functools import reduce
from itertools import permutations
from math import gcd

import sympy as sp


x, y, z = sp.symbols("x y z")
P, r = sp.symbols("P r", nonzero=True)
alpha, beta = sp.symbols("alpha beta", nonzero=True)
epsilon, u, v = sp.symbols("epsilon u v")
rho = sp.symbols("rho", nonzero=True)


# After division by g_1 and removal of the quadratic coefficient, the seed
# G(S)=2S-S^2-2S^3+S^4 has a_3=-1 and a_4=1/2.
a3 = -sp.Integer(1)
a4 = sp.Rational(1, 2)
h = r + a3 * P * r**3 + a4 * P**4 * r**4
B = sp.diff(h, r) / r
C = 2 * h - r * sp.diff(h, r)
J = sp.factor(r**2 * sp.diff(B, r))
expected_J = -1 - 3 * P * r**2 + 4 * P**4 * r**3

assert sp.factor(C + r**2 * B - 2 * h) == 0
assert sp.factor(sp.diff(C, r) + r**2 * sp.diff(B, r)) == 0
assert sp.factor(J - expected_J) == 0


# Enumerate the complete affine GL_2(Z)-stabilizer of the Newton support.
# Exponents are row vectors.  For a permutation pi, solve
#     support[i] * M + shift = support[pi[i]]
# over Q, and retain exactly the integral unimodular matrices.
support = (
    sp.Matrix([[0, 0]]),
    sp.Matrix([[1, 2]]),
    sp.Matrix([[4, 3]]),
)
edge_matrix = sp.Matrix.vstack(
    support[1] - support[0],
    support[2] - support[0],
)
lattice_candidates: list[
    tuple[tuple[int, int, int], sp.Matrix, sp.Matrix]
] = []
nonlattice_candidates: list[
    tuple[tuple[int, int, int], sp.Matrix, sp.Matrix]
] = []

for permutation in permutations(range(3)):
    permuted_edges = sp.Matrix.vstack(
        support[permutation[1]] - support[permutation[0]],
        support[permutation[2]] - support[permutation[0]],
    )
    exponent_matrix = edge_matrix.inv() * permuted_edges
    support_shift = (
        support[permutation[0]] - support[0] * exponent_matrix
    )
    candidate = (permutation, exponent_matrix, support_shift)
    is_integral = all(entry.q == 1 for entry in exponent_matrix)
    is_unimodular = abs(exponent_matrix.det()) == 1
    if is_integral and is_unimodular:
        lattice_candidates.append(candidate)
    else:
        nonlattice_candidates.append(candidate)

identity_matrix = sp.eye(2)
reflection_matrix = sp.Matrix([[-2, -1], [3, 2]])
assert [
    (permutation, matrix, shift)
    for permutation, matrix, shift in lattice_candidates
] == [
    ((0, 1, 2), identity_matrix, sp.zeros(1, 2)),
    ((0, 2, 1), reflection_matrix, sp.zeros(1, 2)),
]
assert nonlattice_candidates == [
    (
        (1, 0, 2),
        sp.Rational(1, 5) * sp.Matrix([[9, 8], [-7, -9]]),
        sp.Matrix([[1, 2]]),
    ),
    (
        (1, 2, 0),
        sp.Rational(1, 5) * sp.Matrix([[-11, -7], [13, 6]]),
        sp.Matrix([[1, 2]]),
    ),
    (
        (2, 0, 1),
        sp.Rational(1, 5) * sp.Matrix([[6, 7], [-13, -11]]),
        sp.Matrix([[4, 3]]),
    ),
    (
        (2, 1, 0),
        sp.Rational(1, 5) * sp.Matrix([[1, -3], [-8, -1]]),
        sp.Matrix([[4, 3]]),
    ),
]
assert all(
    any(entry.q == 5 for entry in matrix)
    for _, matrix, _ in nonlattice_candidates
)
assert reflection_matrix.det() == -1
assert reflection_matrix**2 == identity_matrix


# The nonidentity lattice matrix really preserves the bare Fitting divisor:
# it swaps the two nonconstant terms.  It is rejected only when the
# intrinsic base character P is required to map to a scalar multiple of P.
reflection_alpha = -sp.Rational(4, 3)
reflection_beta = -sp.Rational(3, 4)
reflected_J = sp.factor(
    J.subs(
        {
            P: reflection_beta * P**-2 * r**-1,
            r: reflection_alpha * P**3 * r**2,
        },
        simultaneous=True,
    )
)
assert sp.factor(reflected_J - J) == 0
assert reflection_matrix[0, :] != sp.Matrix([[1, 0]])


# An intrinsic decorated-stratum automorphism has P -> beta P and
# r -> alpha r. Preservation of J gives these two character equations.
transformed_J = sp.expand(J.subs({P: beta * P, r: alpha * r}))
assert transformed_J.coeff(P, 1).coeff(r, 2) == -3 * beta * alpha**2
assert transformed_J.coeff(P, 4).coeff(r, 3) == 4 * beta**4 * alpha**3

stabilizer_equations = (
    beta * alpha**2 - 1,
    beta**4 * alpha**3 - 1,
)
assert sp.factor(
    stabilizer_equations[1].subs(beta, alpha**-2)
    + (alpha**5 - 1) / alpha**5
) == 0


# The dual-number tangent equations are v+2u=0 and 4v+3u=0.
linearized = tuple(
    sp.expand(
        equation.subs(
            {
                alpha: 1 + epsilon * u,
                beta: 1 + epsilon * v,
            }
        )
    ).coeff(epsilon, 1)
    for equation in stabilizer_equations
)
assert linearized == (2 * u + v, 3 * u + 4 * v)
tangent_matrix = sp.Matrix([[2, 1], [3, 4]])
assert tangent_matrix.det() == 5
assert tangent_matrix.rank() == 2


# The finite stabilizer is realized on the displayed determinant-one map.
local_t = 1 + x * y
local_q = local_t**2 * z - y**2 * (1 + 3 * local_t)
mapping = (
    -sp.Rational(1, 2) * local_t * local_q,
    y
    - 3 * x * local_q
    - local_t * local_q
    + 2 * local_t**2 * x**2 * local_q**4,
    x * (5 - 3 * local_t) + x**3 * z - (x * local_q) ** 4,
)

source_action = {x: rho * x, y: y / rho, z: z / rho**2}
source_transformed = tuple(
    sp.together(component.subs(source_action, simultaneous=True))
    for component in mapping
)
target_transformed = (
    mapping[0] / rho**2,
    mapping[1] / rho
    - 2 * (rho**-1 - rho**-2) * mapping[0],
    rho * mapping[2],
)


def reduce_mod_mu5(expression: sp.Expr) -> sp.Expr:
    """Reduce a Laurent expression modulo rho^5-1 exactly."""

    numerator = sp.together(expression).as_numer_denom()[0]
    polynomial = sp.Poly(sp.expand(numerator), rho)
    return sp.factor(polynomial.rem(sp.Poly(rho**5 - 1, rho)).as_expr())


assert all(
    reduce_mod_mu5(got - expected) == 0
    for got, expected in zip(source_transformed, target_transformed)
)

# The target matrices form a representation, including the shear caused by
# the retained quadratic seed coefficient.
eta = sp.symbols("eta", nonzero=True)


def target_matrix(parameter: sp.Expr) -> sp.Matrix:
    return sp.Matrix(
        [
            [parameter**-2, 0, 0],
            [
                -2 * (parameter**-1 - parameter**-2),
                parameter**-1,
                0,
            ],
            [0, 0, parameter],
        ]
    )


assert all(
    sp.factor(entry) == 0
    for entry in (
        target_matrix(rho) * target_matrix(eta) - target_matrix(rho * eta)
    )
)


# A rational-root seed with g_2=0 gives the support-minimal representative
# inside the displayed quartic quadratic-gauge normal form.
sparse_q = local_t**2 * z - sp.Rational(4, 7) * y**2 * (1 + 3 * local_t)
sparse_mapping = (
    -sp.Rational(1, 2) * local_t * sparse_q,
    y
    - sp.Rational(21, 4) * x * sparse_q
    + 3 * local_t**2 * x**2 * sparse_q**4,
    x * (5 - 3 * local_t)
    + sp.Rational(7, 4) * x**3 * z
    - sp.Rational(3, 2) * (x * sparse_q) ** 4,
)
sparse_polynomials = tuple(
    sp.Poly(sp.expand(component), x, y, z) for component in sparse_mapping
)
assert tuple(len(poly.terms()) for poly in sparse_polynomials) == (7, 51, 38)
assert tuple(poly.total_degree() for poly in sparse_polynomials) == (7, 26, 24)
assert sp.factor(sp.det(sp.Matrix(sparse_mapping).jacobian((x, y, z)))) == 1

sparse_collision = (
    (0, 0, 1),
    (-sp.Rational(4, 5), sp.Rational(9, 4), -sp.Rational(265, 32)),
    (sp.Rational(1, 2), -sp.Rational(3, 2), 100),
    (
        sp.Rational(3, 10),
        -sp.Rational(29, 6),
        -sp.Rational(24820, 729),
    ),
)
assert len(set(sparse_collision)) == 4
assert all(
    tuple(
        sp.factor(component.subs(dict(zip((x, y, z), point))))
        for component in sparse_mapping
    )
    == (-sp.Rational(1, 2), 0, 0)
    for point in sparse_collision
)

# Every symbolic coefficient is one Laurent monomial in the nonzero
# parameters a_3,a_4, with exactly seven additional coefficients containing
# a_2.  Thus a_2=0 is the unique support drop in this normal form.
a2, generic_a3, generic_a4 = sp.symbols(
    "a2 generic_a3 generic_a4", nonzero=True
)
generic_q = (
    local_t**2 * z
    + y**2 * (1 + 3 * local_t) / generic_a3
)
generic_mapping = (
    -sp.Rational(1, 2) * local_t * generic_q,
    y
    + 3 * generic_a3 * x * generic_q
    + 2 * a2 * local_t * generic_q
    + 4 * generic_a4 * local_t**2 * x**2 * generic_q**4,
    x * (5 - 3 * local_t)
    - generic_a3 * x**3 * z
    - 2 * generic_a4 * (x * generic_q) ** 4,
)
generic_coefficients = tuple(
    tuple(
        coefficient
        for _, coefficient in sp.Poly(
            sp.expand(component), x, y, z
        ).terms()
    )
    for component in generic_mapping
)
assert tuple(len(coefficients) for coefficients in generic_coefficients) == (
    7,
    58,
    38,
)


def is_laurent_monomial(expression: sp.Expr) -> bool:
    """Return whether expression is one nonzero Laurent parameter monomial."""

    numerator, denominator = sp.cancel(expression).as_numer_denom()
    parameters = (a2, generic_a3, generic_a4)
    return (
        len(sp.Poly(numerator, *parameters).terms()) == 1
        and len(sp.Poly(denominator, *parameters).terms()) == 1
    )


assert all(
    is_laurent_monomial(coefficient)
    for coefficients in generic_coefficients
    for coefficient in coefficients
)
a2_coefficients = tuple(
    coefficient
    for coefficients in generic_coefficients
    for coefficient in coefficients
    if coefficient.has(a2)
)
assert len(a2_coefficients) == 7
assert all(
    sp.factor(coefficient / a2).subs(a2, 0) != 0
    for coefficient in a2_coefficients
)
assert all(
    coefficient.subs(a2, 0) != 0
    for coefficients in generic_coefficients
    for coefficient in coefficients
    if not coefficient.has(a2)
)

# With g_2=0 the residual mu_5 action is diagonal on both sides.
sparse_source_transformed = tuple(
    sp.together(component.subs(source_action, simultaneous=True))
    for component in sparse_mapping
)
sparse_target_transformed = (
    sparse_mapping[0] / rho**2,
    sparse_mapping[1] / rho,
    rho * sparse_mapping[2],
)
assert all(
    reduce_mod_mu5(got - expected) == 0
    for got, expected in zip(
        sparse_source_transformed, sparse_target_transformed
    )
)

# The sparse representative has its own small affine-linear certificate.
# The column order is A(9), source translation(3), B(9), and target
# translation(3), with both matrices in row-major order.  The 24 selected
# coefficient rows have 46 nonzero primitive-integer entries.
sparse_a_symbols = sp.symbols(
    "sparse_a11 sparse_a12 sparse_a13 "
    "sparse_a21 sparse_a22 sparse_a23 "
    "sparse_a31 sparse_a32 sparse_a33"
)
sparse_source_constants = sp.symbols("sparse_source_constant_1:4")
sparse_b_symbols = sp.symbols(
    "sparse_b11 sparse_b12 sparse_b13 "
    "sparse_b21 sparse_b22 sparse_b23 "
    "sparse_b31 sparse_b32 sparse_b33"
)
sparse_target_constants = sp.symbols("sparse_target_constant_1:4")
sparse_affine_unknowns = (
    sparse_a_symbols
    + sparse_source_constants
    + sparse_b_symbols
    + sparse_target_constants
)
sparse_A = sp.Matrix(3, 3, sparse_a_symbols)
sparse_B = sp.Matrix(3, 3, sparse_b_symbols)
sparse_vector = sp.Matrix(sparse_mapping)
sparse_affine_identity = (
    sparse_B * sparse_vector
    + sp.Matrix(sparse_target_constants)
    - sparse_vector.jacobian((x, y, z))
    * (
        sparse_A * sp.Matrix((x, y, z))
        + sp.Matrix(sparse_source_constants)
    )
)
sparse_affine_polynomials = tuple(
    sp.Poly(sp.expand(component), x, y, z)
    for component in sparse_affine_identity
)
sparse_affine_labels = (
    (1, (12, 10, 4)),
    (1, (12, 8, 4)),
    (1, (4, 3, 0)),
    (1, (4, 2, 1)),
    (1, (3, 4, 0)),
    (1, (3, 3, 1)),
    (1, (2, 2, 0)),
    (1, (3, 2, 2)),
    (2, (12, 9, 4)),
    (1, (2, 4, 1)),
    (1, (2, 4, 0)),
    (1, (2, 3, 2)),
    (1, (1, 4, 0)),
    (1, (2, 2, 1)),
    (1, (0, 0, 0)),
    (2, (12, 10, 4)),
    (2, (12, 8, 4)),
    (2, (3, 3, 1)),
    (2, (3, 2, 1)),
    (2, (0, 0, 0)),
    (3, (12, 10, 4)),
    (3, (12, 8, 4)),
    (3, (3, 3, 1)),
    (3, (0, 0, 0)),
)


def primitive_integer_row(row: list[sp.Expr]) -> list[int]:
    """Clear denominators and common factors with positive leading entry."""

    denominator_lcm = sp.ilcm(*(sp.denom(entry) for entry in row))
    integer_row = [int(denominator_lcm * entry) for entry in row]
    common_factor = reduce(
        gcd, (abs(entry) for entry in integer_row if entry)
    )
    integer_row = [entry // common_factor for entry in integer_row]
    first_nonzero = next(entry for entry in integer_row if entry)
    if first_nonzero < 0:
        integer_row = [-entry for entry in integer_row]
    return integer_row


sparse_affine_minor = sp.Matrix(
    [
        primitive_integer_row(
            [
                sparse_affine_polynomials[component - 1]
                .coeff_monomial(monomial)
                .coeff(unknown)
                for unknown in sparse_affine_unknowns
            ]
        )
        for component, monomial in sparse_affine_labels
    ]
)
assert sparse_affine_minor.shape == (24, 24)
assert sum(entry != 0 for entry in sparse_affine_minor) == 46
assert sparse_affine_minor.det() == 10


# The rational stable-moduli scaling (alpha,beta)=(1/4,12/5) improves
# coefficient and collision height without changing support.
balanced_q = (
    local_t**2 * z - sp.Rational(3, 35) * y**2 * (1 + 3 * local_t)
)
balanced_mapping = (
    -sp.Rational(1, 2) * local_t * balanced_q,
    y
    - 35 * x * balanced_q
    + sp.Rational(625, 108) * local_t**2 * x**2 * balanced_q**4,
    x * (5 - 3 * local_t)
    + sp.Rational(35, 3) * x**3 * z
    - sp.Rational(625, 216) * (x * balanced_q) ** 4,
)
balanced_polynomials = tuple(
    sp.Poly(sp.expand(component), x, y, z)
    for component in balanced_mapping
)
assert tuple(len(poly.terms()) for poly in balanced_polynomials) == (7, 51, 38)
assert tuple(poly.total_degree() for poly in balanced_polynomials) == (7, 26, 24)
assert sp.factor(sp.det(sp.Matrix(balanced_mapping).jacobian((x, y, z)))) == 1

balanced_collision = (
    (0, 0, sp.Rational(12, 5)),
    (-sp.Rational(1, 5), 9, -sp.Rational(159, 8)),
    (sp.Rational(1, 8), -6, 240),
    (
        sp.Rational(3, 40),
        -sp.Rational(58, 3),
        -sp.Rational(19856, 243),
    ),
)
assert all(
    tuple(
        sp.factor(component.subs(dict(zip((x, y, z), point))))
        for component in balanced_mapping
    )
    == (-sp.Rational(6, 5), 0, 0)
    for point in balanced_collision
)


def rational_height(value: sp.Expr) -> int:
    """Naive reduced height max(abs(numerator), denominator)."""

    numerator, denominator = sp.Rational(value).as_numer_denom()
    return max(abs(int(numerator)), int(denominator))


def expanded_coefficient_height(polynomials: tuple[sp.Poly, ...]) -> int:
    return max(
        rational_height(coefficient)
        for polynomial in polynomials
        for _, coefficient in polynomial.terms()
    )


def point_height(points: tuple[tuple[sp.Expr, ...], ...]) -> int:
    return max(rational_height(value) for point in points for value in point)


assert expanded_coefficient_height(sparse_polynomials) == 2248704
assert expanded_coefficient_height(balanced_polynomials) == 21875
assert point_height(sparse_collision) == 24820
assert point_height(balanced_collision) == 19856

print("PASS: J(P,r) = -1 - 3*P*r^2 + 4*P^4*r^3")
print("PASS: exactly two Newton-support lattice matrices survive")
print("PASS: the sole involution is rejected by the intrinsic base character")
print("PASS: the decorated stabilizer is beta=alpha^-2, alpha^5=1")
print("PASS: its tangent matrix has determinant 5 and zero nullspace")
print("PASS: the residual mu_5 is realized by explicit linear LR symmetries")
print("PASS: sparse quartic has support (7,51,38) and degree (7,26,24)")
print("PASS: sparse quartic has determinant one and a rational four-point fiber")
print("PASS: g_2=0 gives the unique seven-term support drop in this normal form")
print("PASS: sparse affine-linear 24-by-24 minor has determinant 10")
print("PASS: balanced LR scaling lowers coefficient height 2248704 -> 21875")
print("PASS: balanced LR scaling lowers collision height 24820 -> 19856")
