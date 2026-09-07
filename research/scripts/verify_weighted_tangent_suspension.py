#!/usr/bin/env python3
"""Exact certificate for the weighted tangent-suspension plane core."""

import sympy as sp


u, v, s, t, X, W, gamma = sp.symbols("u v s t X W gamma")

A20, A11, A30, A21, A40, A31 = sp.symbols(
    "A20 A11 A30 A21 A40 A31"
)
B01, B20, B11, B30, B21 = sp.symbols("B01 B20 B11 B30 B21")

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
C = 2 - 3 * u - v

# The invariant-plane coordinate change has the required determinant -1.
coordinate_matrix = sp.Matrix([[-3, -1], [-1, 0]])
assert coordinate_matrix.det() == -1

# Foundational invariant polynomials and their compact Poisson-square identity.
foundational = {
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
}
plane_substitution = {u: -t, v: 3 * t - s}
foundational_A = sp.factor(A.subs(foundational).subs(plane_substitution))
foundational_B = sp.factor(B.subs(foundational).subs(plane_substitution))
expected_A = (t - 1) * (s * (t - 1) ** 2 + t * (2 * t - 3))
expected_B = -3 * s * (t - 1) ** 2 - 2 * t * (3 * t - 4)
assert sp.expand(foundational_A - expected_A) == 0
assert sp.expand(foundational_B - expected_B) == 0

plane_C = s + 2
P = sp.expand(plane_C**2 * foundational_A)
Q = sp.expand(plane_C * foundational_B)
oriented_bracket = sp.diff(P, t) * sp.diff(Q, s) - sp.diff(P, s) * sp.diff(Q, t)
assert sp.factor(oriented_bracket) == -2 * (s + 2) ** 2

# The pair is the tangent map for H(W)=W^2(1-W), up to target scaling by 4.
H = W**2 * (1 - W)
tangent_substitution = {s: 2 * gamma - 2, t: 1 - W / gamma}
assert sp.factor(Q.subs(tangent_substitution) / 4 - (sp.diff(H, W) + gamma)) == 0
assert sp.factor(
    P.subs(tangent_substitution) / 4
    - (W * (sp.diff(H, W) + gamma) - H)
) == 0
assert sp.factor(4 * H - Q.subs(tangent_substitution) * W + P.subs(tangent_substitution)) == 0

# Rewrite the complete normalized ansatz in four univariate X-bands.
band_substitution = {u: -t, v: 3 * t - X + 2}
band_A = sp.expand(A.subs(band_substitution))
band_B = sp.expand(B.subs(band_substitution))
widehat_P = sp.expand(X**2 * band_A / 2)
band_Q = sp.expand(X * band_B)
p3 = sp.expand(widehat_P.coeff(X, 3))
p2 = sp.expand(widehat_P.coeff(X, 2))
q2 = sp.expand(band_Q.coeff(X, 2))
q1 = sp.expand(band_Q.coeff(X, 1))
assert sp.expand(widehat_P - p3 * X**3 - p2 * X**2) == 0
assert sp.expand(band_Q - q2 * X**2 - q1 * X) == 0

def wronskian(i, j, f, g):
    return sp.expand(i * f * sp.diff(g, t) - j * sp.diff(f, t) * g)


K4 = wronskian(3, 2, p3, q2)
K3 = sp.expand(wronskian(3, 1, p3, q1) + wronskian(2, 2, p2, q2))
K2 = wronskian(2, 1, p2, q1)
standard_bracket = sp.expand(
    sp.diff(widehat_P, X) * sp.diff(band_Q, t)
    - sp.diff(widehat_P, t) * sp.diff(band_Q, X)
)
assert sp.expand(standard_bracket - (K4 * X**4 + K3 * X**3 + K2 * X**2)) == 0

invariant_determinant = sp.Matrix(
    [
        (-2 * A, sp.diff(A, u), sp.diff(A, v)),
        (-B, sp.diff(B, u), sp.diff(B, v)),
        (C, sp.diff(C, u), sp.diff(C, v)),
    ]
).det()
transformed_keller_equation = sp.expand(
    (invariant_determinant + 2).subs(band_substitution)
)
assert sp.expand(
    X**2 * transformed_keller_equation
    + 2 * (standard_bracket - X**2)
) == 0

# The Wronskian cascade defines the same dual-number ideal as the earlier
# invariant-determinant computation.
variables = (
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
cascade_equations = []
for polynomial in (K4, K3, K2 - 1):
    cascade_equations.extend(sp.Poly(polynomial, t).all_coeffs())
cascade = sp.groebner(cascade_equations, *variables, order="grevlex", domain=sp.QQ)

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
    triangular_relations, *variables, order="grevlex", domain=sp.QQ
)
assert all(triangular.reduce(poly.as_expr())[1] == 0 for poly in cascade.polys)
assert all(cascade.reduce(poly.as_expr())[1] == 0 for poly in triangular.polys)

# Along the unique affine lift of the dual-number direction, every band defect
# is visibly quadratic.  The leading homogeneous layer already detects
# epsilon^2.
epsilon = sp.symbols("epsilon")
linear_lift = {
    A20: (7 * B01 + 27) / 12,
    A11: (B01 + 3) / 2,
    A30: (9 * B01 + 1) / 4,
    A21: B01,
    A40: (3 * B01 - 3) / 2,
    A31: (B01 - 1) / 2,
    B20: (11 * B01 - 9) / 2,
    B11: 3 * B01 - 3,
    B30: 6 * B01 - 9,
    B21: 2 * B01 - 3,
}
epsilon_substitution = {B01: 3 + epsilon}
assert sp.factor(K4.subs(linear_lift).subs(epsilon_substitution)) == (
    -epsilon**2 * (t - 2) * (t - 1) ** 2 / 4
)
assert sp.factor(K3.subs(linear_lift).subs(epsilon_substitution)) == (
    -epsilon**2 * (3 * t**4 - 12 * t**3 + 33 * t**2 - 64 * t + 36) / 24
)
assert sp.factor((K2 - 1).subs(linear_lift).subs(epsilon_substitution)) == (
    -epsilon**2 * (t - 2) * (t + 2) * (3 * t**2 - 8 * t + 12) / 48
)

assert sp.factor(p3.subs(foundational)) == (t - 1) ** 3 / 2
assert sp.factor(p2.subs(foundational)) == (t - 2) * (t - 1) / 2
assert sp.factor(q2.subs(foundational)) == -3 * (t - 1) ** 2
assert sp.expand(q1.subs(foundational) + 2 * (2 * t - 3)) == 0
assert K4.subs(foundational) == 0
assert K3.subs(foundational) == 0
assert K2.subs(foundational) == 1

print("PASS: oriented plane-core determinant equals the Poisson square")
print("PASS: foundational P,Q reduce to the tangent pencil H(W)-sW+t")
print("PASS: three Wronskian layers recover the normalized dual-number scheme")
print("PASS: the leading Wronskian layer already detects the quadratic obstruction")
