#!/usr/bin/env python3
"""Exact checks for the absolute D5 affine-modification frontier.

The all-degree affine-linear mask theorem is proved in the accompanying
note.  This checker verifies its maximal-minor determinant ledger for two
and three auxiliary coordinates, together with the D5 Cox geometry.
"""

from __future__ import annotations

import sympy as sp


a, u, v = sp.symbols("a u v")
sqrt5 = sp.sqrt(5)
alpha = (3 + sqrt5) / 2
beta = (3 - sqrt5) / 2

P = a**5 - 5 * a**3 * u + 5 * a * u**2
J = sp.diff(P, a)
C = a**2 - 4 * u
R_plus = a**2 - alpha * u
R_minus = a**2 - beta * u
Q = R_plus * R_minus

assert sp.expand(alpha + beta - 3) == 0
assert sp.expand(alpha * beta - 1) == 0
assert sp.expand(J - 5 * Q) == 0
assert sp.expand(P**2 - 4 * u**5 - C * Q**2) == 0
assert sp.factor(sp.discriminant(P - v, a)) == 5**5 * (4 * u**5 - v**2) ** 2

# In the source boundary basis (C,R_plus,R_minus), the target branch pulls
# back with vector (1,2,2), whereas J has vector (0,1,1).
branch_vector = sp.Matrix([1, 2, 2])
derivative_vector = sp.Matrix([0, 1, 1])
assert branch_vector.rank() == 1
assert branch_vector.cross(derivative_vector) != sp.zeros(3, 1)
print("PASS: the D5 derivative and branch-pullback ledgers are exact")

# The two ramification colors have intersection multiplicity two.
R_minus_on_R_plus = sp.expand(R_minus.subs(u, a**2 / alpha))
assert sp.expand(
    R_minus_on_R_plus - (1 - beta / alpha) * a**2
) == 0
assert sp.expand(1 - beta / alpha) != 0
print("PASS: the two D5 ramification colors are tangent to order two")

# The total product fill xy=Q has a unique singular point.  Its singular
# equations reduce to the maximal ideal (x,y,a,u).
x, y = sp.symbols("x y")
product_singular_ideal = sp.groebner(
    [x, y, sp.diff(Q, a), sp.diff(Q, u)],
    x,
    y,
    u,
    a,
    extension=sqrt5,
)
nilpotence_exponents = {x: 1, y: 1, u: 2, a: 3}
for variable, exponent in nilpotence_exponents.items():
    assert product_singular_ideal.reduce(variable**exponent)[1] == 0
assert all(
    polynomial.as_expr().subs({x: 0, y: 0, a: 0, u: 0}) == 0
    for polynomial in product_singular_ideal.polys
)
print("PASS: the total-product Danielewski fill has only the singular origin")

# Filling either R color alone eliminates u and gives a polynomial ring.
xp, yp, xm, ym = sp.symbols("xp yp xm ym")
u_from_plus = (a**2 - xp * yp) / alpha
assert sp.expand((xp * yp - R_plus).subs(u, u_from_plus)) == 0

# Filling both colors gives a nondegenerate five-variable quadratic cone.
separated_relation = sp.expand(
    xm * ym
    - beta / alpha * xp * yp
    - (1 - beta / alpha) * a**2
)
assert sp.expand((xm * ym - R_minus).subs(u, u_from_plus) - separated_relation) == 0
separated_variables = (a, xp, yp, xm, ym)
separated_hessian = sp.hessian(separated_relation, separated_variables)
assert sp.simplify(separated_hessian.det()) != 0
assert all(
    sp.diff(separated_relation, variable).subs(
        {a: 0, xp: 0, yp: 0, xm: 0, ym: 0}
    )
    == 0
    for variable in separated_variables
)
print("PASS: the separated two-color fill is a nondegenerate quadric cone")

# Adding the unramified color C in separated form has dependent conormal
# rows at the common origin.
gamma_values = (sp.Integer(4), alpha, beta)
conormal_at_origin = sp.Matrix(
    [[0, -gamma, 0, 0, 0, 0, 0, 0] for gamma in gamma_values]
)
assert conormal_at_origin.rank() == 1
print("PASS: the naive three-color Cox fill is singular at the common vertex")


def verify_affine_linear_mask_ledger(auxiliary_count: int) -> None:
    """Verify det[h'+M'z,M] = p.(h'+M'z) and pM=0."""

    row_count = auxiliary_count + 1
    matrix_symbols = sp.symbols(f"m0:{row_count * auxiliary_count}")
    derivative_symbols = sp.symbols(f"d0:{row_count * auxiliary_count}")
    hprime_symbols = sp.symbols(f"h0:{row_count}")
    z_symbols = sp.symbols(f"z0:{auxiliary_count}")

    matrix = sp.Matrix(row_count, auxiliary_count, matrix_symbols)
    matrix_derivative = sp.Matrix(
        row_count, auxiliary_count, derivative_symbols
    )
    hprime = sp.Matrix(hprime_symbols)
    z_vector = sp.Matrix(z_symbols)
    first_column = hprime + matrix_derivative * z_vector

    signed_minors = sp.Matrix(
        1,
        row_count,
        [
            (-1) ** row
            * matrix[
                [index for index in range(row_count) if index != row],
                :,
            ].det()
            for row in range(row_count)
        ],
    )

    assert all(sp.expand(entry) == 0 for entry in signed_minors * matrix)
    direct = sp.Matrix.hstack(first_column, matrix).det()
    ledger = (signed_minors * first_column)[0]
    assert sp.expand(direct - ledger) == 0

    # Constancy in the auxiliary variables is exactly p M'=0.
    for index, variable in enumerate(z_symbols):
        assert sp.expand(sp.diff(direct, variable) - (signed_minors * matrix_derivative[:, index])[0]) == 0


verify_affine_linear_mask_ledger(2)
verify_affine_linear_mask_ledger(3)
print("PASS: affine-linear mask ledgers hold for two and three auxiliaries")

# Any polynomial thickening that restricts to (u,P,0,...,0) on the
# auxiliary zero section has block-triangular derivative there.
for auxiliary_count in (2, 3):
    top_right = sp.Matrix(
        2,
        auxiliary_count,
        sp.symbols(f"b{auxiliary_count}_0:{2 * auxiliary_count}"),
    )
    lower_right = sp.Matrix(
        auxiliary_count,
        auxiliary_count,
        sp.symbols(
            f"c{auxiliary_count}_0:{auxiliary_count * auxiliary_count}"
        ),
    )
    zero_section_derivative = sp.Matrix.vstack(
        sp.Matrix.hstack(
            sp.Matrix([u, P]).jacobian((a, u)),
            top_right,
        ),
        sp.Matrix.hstack(
            sp.zeros(auxiliary_count, 2),
            lower_right,
        ),
    )
    assert sp.expand(
        zero_section_derivative.det() + J * lower_right.det()
    ) == 0
print("PASS: zero-section thickenings retain the J5 Jacobian divisor")
print("PASS absolute D5 affine-modification frontier")
