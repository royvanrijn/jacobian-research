#!/usr/bin/env python3
"""Exact certificate for the cotangent and Weyl lifts of the foundational map."""

import sympy as sp

x, y, z = sp.symbols("x y z")
variables = (x, y, z)
u = 1 + x*y
F = sp.Matrix([
    u**3*z + y**2*u*(4 + 3*x*y),
    y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y),
    2*x - 3*x**2*y - x**3*z,
])

J = F.jacobian(variables)
B = -sp.Rational(1, 2) * J.adjugate()

assert sp.factor(J.det()) == -2
assert (J*B).applyfunc(sp.expand) == sp.eye(3)
assert (B*J).applyfunc(sp.expand) == sp.eye(3)


def derivation(vector_field: sp.Matrix, polynomial: sp.Expr) -> sp.Expr:
    """Apply a polynomial vector field to a polynomial."""
    return sp.expand(sum(
        vector_field[index] * sp.diff(polynomial, variable)
        for index, variable in enumerate(variables)
    ))


# The columns of B are the inverse-Jacobian derivations delta_i.
for i in range(3):
    for j in range(3):
        expected = 1 if i == j else 0
        assert sp.expand(derivation(B[:, i], F[j]) - expected) == 0

# Their pairwise Lie brackets vanish, hence the Weyl D-images commute.
for i in range(3):
    for j in range(3):
        bracket = sp.Matrix([
            sp.expand(
                derivation(B[:, i], B[row, j])
                - derivation(B[:, j], B[row, i])
            )
            for row in range(3)
        ])
        assert bracket == sp.zeros(3, 1)

# Exact preservation of the canonical one-form: J^T B^T p = p.
p = sp.Matrix(sp.symbols("p1 p2 p3"))
target_momentum = B.T*p
assert (J.T*target_momentum - p).applyfunc(sp.expand) == sp.zeros(3, 1)

# The block-triangular Jacobian of the six-variable lift has determinant one.
assert sp.expand(J.det() * B.T.det() - 1) == 0

# The rational three-point collision persists over arbitrary target momentum.
q = sp.Matrix(sp.symbols("q1 q2 q3"))
points = [
    (0, 0, -sp.Rational(1, 4)),
    (1, -sp.Rational(3, 2), sp.Rational(13, 2)),
    (-1, sp.Rational(3, 2), sp.Rational(13, 2)),
]
target = sp.Matrix([-sp.Rational(1, 4), 0, 0])
for point in points:
    substitution = dict(zip(variables, point))
    assert F.subs(substitution) == target

    source_momentum = J.subs(substitution).T*q
    lifted_momentum = B.subs(substitution).T*source_momentum
    assert (lifted_momentum - q).applyfunc(sp.expand) == sp.zeros(3, 1)

print("PASS: det(JF) = -2 and B = J^-1 is polynomial")
print("PASS: inverse-Jacobian derivations satisfy the Weyl relations")
print("PASS: the cotangent lift preserves theta and has determinant one")
print("PASS: the three-point collision persists for symbolic target momentum")