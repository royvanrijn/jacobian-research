#!/usr/bin/env python3
"""Exact certificate for a Keller counterexample with no linear torus symmetry."""

from __future__ import annotations

import sympy as sp


x, y, z = sp.symbols("x y z")
source_variables = (x, y, z)

# The root-engineered quadratic-gauge map for
# G(S) = S(S-1)(S+1)(S-2) = 2S - S^2 - 2S^3 + S^4.
t = 1 + x * y
q = t**2 * z - y**2 * (1 + 3 * t)
F = sp.Matrix(
    [
        -sp.Rational(1, 2) * t * q,
        y - 3 * x * q - t * q + 2 * t**2 * x**2 * q**4,
        x * (5 - 3 * t) + x**3 * z - (x * q) ** 4,
    ]
).applyfunc(sp.expand)

# Keller certificate.
determinant = sp.factor(F.jacobian(source_variables).det())
assert determinant == 1

# Collision certificate.
collision_points = (
    (sp.Rational(0), sp.Rational(1), sp.Rational(5)),
    (-sp.Rational(1), sp.Rational(2), -sp.Rational(9)),
    (sp.Rational(1, 3), -sp.Rational(4), -sp.Rational(27)),
    (sp.Rational(2, 3), -sp.Rational(1), sp.Rational(45)),
)
for point in collision_points:
    assert F.subs(dict(zip(source_variables, point))) == sp.Matrix(
        [-sp.Rational(1, 2), 0, 0]
    )

# Coefficientwise linear-symmetry system
#
#     B F(x) = JF(x) A x.
#
# Its 18 columns correspond, in row-major order, to the entries of A and B.
a_symbols = sp.symbols("a11 a12 a13 a21 a22 a23 a31 a32 a33")
b_symbols = sp.symbols("b11 b12 b13 b21 b22 b23 b31 b32 b33")
A = sp.Matrix(3, 3, a_symbols)
B = sp.Matrix(3, 3, b_symbols)
unknowns = a_symbols + b_symbols

# A stronger affine-linear check allows constant terms in both vector fields.
# It excludes torus actions that become linear only after translating source
# and target origins.
source_constants = sp.symbols("source_constant_1:4")
target_constants = sp.symbols("target_constant_1:4")
affine_unknowns = (
    a_symbols + source_constants + b_symbols + target_constants
)
affine_identity = (
    B * F
    + sp.Matrix(target_constants)
    - F.jacobian(source_variables)
    * (A * sp.Matrix(source_variables) + sp.Matrix(source_constants))
)
affine_rows = []
for polynomial in affine_identity:
    poly = sp.Poly(sp.expand(polynomial), source_variables)
    for _, coefficient in poly.terms():
        affine_rows.append(
            [
                sp.expand(coefficient).coeff(unknown)
                for unknown in affine_unknowns
            ]
        )
affine_coefficient_matrix = sp.Matrix(affine_rows)
assert affine_coefficient_matrix.cols == 24
assert affine_coefficient_matrix.nullspace() == []

certificate_labels = (
    (1, (12, 10, 4)),
    (1, (12, 8, 4)),
    (1, (4, 3, 0)),
    (1, (4, 2, 1)),
    (1, (3, 4, 0)),
    (1, (3, 3, 1)),
    (1, (3, 2, 2)),
    (1, (2, 4, 1)),
    (1, (2, 4, 0)),
    (1, (2, 3, 2)),
    (1, (2, 2, 1)),
    (2, (12, 10, 4)),
    (2, (12, 8, 4)),
    (2, (3, 3, 1)),
    (2, (3, 2, 1)),
    (3, (12, 10, 4)),
    (3, (12, 8, 4)),
    (3, (3, 3, 1)),
)


def primitive_integer_row(row: list[sp.Expr]) -> list[int]:
    """Clear denominators and common factors, with positive first entry."""

    denominator_lcm = sp.ilcm(*(sp.denom(entry) for entry in row))
    integer_row = [int(denominator_lcm * entry) for entry in row]
    common_factor = sp.igcd(*(abs(entry) for entry in integer_row))
    integer_row = [entry // common_factor for entry in integer_row]
    first_nonzero = next(entry for entry in integer_row if entry)
    if first_nonzero < 0:
        integer_row = [-entry for entry in integer_row]
    return integer_row


def selected_symmetry_minor(mapping: sp.Matrix) -> sp.Matrix:
    """Extract the fixed 18 coefficient rows from B F - JF A x."""

    identity = (
        B * mapping
        - mapping.jacobian(source_variables)
        * A
        * sp.Matrix(source_variables)
    )
    polynomials = [
        sp.Poly(sp.expand(polynomial), source_variables)
        for polynomial in identity
    ]
    return sp.Matrix(
        [
            [
                sp.diff(
                    polynomials[component - 1].coeff_monomial(monomial),
                    unknown,
                )
                for unknown in unknowns
            ]
            for component, monomial in certificate_labels
        ]
    )


selected_matrix = selected_symmetry_minor(F)
assert selected_matrix.shape == (18, 18)
certificate_matrix = sp.Matrix(
    [primitive_integer_row(row) for row in selected_matrix.tolist()]
)
certificate_determinant = sp.factor(certificate_matrix.det())
assert certificate_determinant == -5

# The same fixed rows settle every admissible quartic seed at once.
g1, g2, g3, g4 = sp.symbols("g1 g2 g3 g4", nonzero=True)
generic_q = t**2 * z + (g1 / g3) * y**2 * (1 + 3 * t)
generic_F = sp.Matrix(
    [
        -sp.Rational(1, 2) * t * generic_q,
        y
        + 3 * g3 * x * generic_q / g1
        + 2 * g2 * t * generic_q / g1
        + 4 * g4 * t**2 * x**2 * generic_q**4 / g1,
        x * (5 - 3 * t)
        - g3 * x**3 * z / g1
        - 2 * g4 * (x * generic_q) ** 4 / g1,
    ]
).applyfunc(sp.expand)
generic_minor_determinant = sp.factor(selected_symmetry_minor(generic_F).det())
assert generic_minor_determinant == (
    sp.Rational(10935, 4) * g4**6 / g1**6
)

print("PASS: det(JF) =", determinant)
print("PASS: four rational source points map to (-1/2, 0, 0):")
for point in collision_points:
    print(" ", point)
print("PASS: the selected 18 coefficient rows have full rank")
print(
    "PASS: the affine-linear coefficient matrix has shape",
    affine_coefficient_matrix.shape,
    "and nullspace zero",
)
print("PASS: an 18-by-18 certificate minor has determinant", certificate_determinant)
print("PASS: the universal quartic minor is", generic_minor_determinant)
print("certificate rows (component, exponent triple):")
for label, row in zip(
    certificate_labels,
    certificate_matrix.tolist(),
):
    print(label, row)
