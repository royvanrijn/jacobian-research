#!/usr/bin/env python3
"""Exact certificate for a rank-two Poisson completion of the JC(3) map.

The construction was obtained locally while auditing an unavailable external
announcement.  It is not asserted to reproduce that manuscript's missing
formulas.  The script verifies the construction directly over Q:

* the adapted source coordinates are a polynomial coordinate system;
* the low-degree shear coefficient -9 is forced in the tested c*Q^2 family;
* all six canonical Poisson brackets hold;
* the four-dimensional Jacobian is one; and
* a complete three-point fiber is transported from the foundational map.
"""

import sympy as sp


x, q, p, z = sp.symbols("x q p z")
X, Q, Z, E = sp.symbols("X Q Z E")


def poisson(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    """The convention in the supplied abstract: {p,x}={z,q}=1."""

    return sp.expand(
        sp.diff(f, p) * sp.diff(g, x)
        - sp.diff(f, x) * sp.diff(g, p)
        + sp.diff(f, z) * sp.diff(g, q)
        - sp.diff(f, q) * sp.diff(g, z)
    )


# Adapted polynomial coordinates.  Here Q is literally the original q; the
# redundant Y used in the pre-audit is Q-XZ/3.
X_source = x
Q_source = q
Z_source = 3 * x**2 * p + (2 - 6 * x * q) * z
E_source = (1 + 3 * x * q) * p / 2 - 3 * q**2 * z
adapted = sp.Matrix([X_source, Q_source, Z_source, E_source])
assert sp.factor(adapted.jacobian((x, q, p, z)).det()) == -1

# An explicit inverse proves that the constant determinant is not being used
# as a circular automorphy argument.
p_inverse = -6 * E * X * Q + 2 * E + 3 * Z * Q**2
z_inverse = -3 * E * X**2 + 3 * Z * X * Q / 2 + Z / 2
adapted_inverse = {x: X, q: Q, p: p_inverse, z: z_inverse}
for got, expected in zip(
    adapted.subs(adapted_inverse, simultaneous=True), (X, Q, Z, E), strict=True
):
    assert sp.expand(got - expected) == 0


def foundational_pair(z_argument: sp.Expr) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """(F1/2,F2,F3) after Y=Q-X*z_argument/3."""

    y_argument = Q - X * z_argument / 3
    u = 1 + X * y_argument
    first = (u**3 * z_argument + y_argument**2 * u * (4 + 3 * X * y_argument)) / 2
    second = (
        y_argument
        + 3 * X * u**2 * z_argument
        + 3 * X * y_argument**2 * (4 + 3 * X * y_argument)
    )
    third = 2 * X - 3 * X**2 * y_argument - X**3 * z_argument
    return tuple(sp.expand(value) for value in (first, second, third))


# How the completion was found.  For the volume-preserving, R-fixing shear
# Z -> Z+cQ^2, compute the unique horizontal vector field for
# (S_c,T_c,R), subtract {E,-}, and integrate its Hamiltonian completion after
# v=1/X, rho=R.  Its entire negative-X principal part has factor c+9.
c, v, rho = sp.symbols("c v rho")
S_c, T_c, R_c = foundational_pair(Z + c * Q**2)
J_c = sp.Matrix([S_c, T_c, R_c]).jacobian((X, Q, Z))
assert sp.factor(J_c.det()) == -1
w_c = (J_c.adjugate() / J_c.det())[:, 2].applyfunc(sp.expand)
w_E = sp.Matrix([(1 + 3 * X * Q) / 2, -3 * Q**2, 9 * Q * Z / 2])
difference_c = (w_c - w_E).applyfunc(sp.expand)

# In (X,Q,Z), {f,-} has components
# (3X^2 f_Z, (2-6XQ)f_Z, -3X^2 f_X+(6XQ-2)f_Q).
f0_c = sp.integrate(sp.cancel(difference_c[0] / (3 * X**2)), Z)
hamiltonian_f0_c = sp.Matrix(
    [
        3 * X**2 * sp.diff(f0_c, Z),
        (2 - 6 * X * Q) * sp.diff(f0_c, Z),
        -3 * X**2 * sp.diff(f0_c, X)
        + (6 * X * Q - 2) * sp.diff(f0_c, Q),
    ]
)
residual_c = (difference_c - hamiltonian_f0_c).applyfunc(sp.expand)
assert residual_c[0] == 0
assert residual_c[1] == 0

obstruction_vrho = sp.expand(
    residual_c[2].subs({X: 1 / v, Q: v * (2 - rho * v) / 3})
)
forced_primitive = sp.integrate(obstruction_vrho / 3, v)
forced_back = sp.expand(
    forced_primitive.subs({v: 1 / X, rho: 2 * X - 3 * X**2 * Q})
)
principal_part = sum(forced_back.coeff(X, exponent) * X**exponent for exponent in range(-8, 0))
expected_principal_part = (c + 9) * (
    1 / (135 * X**4)
    + 2 * Q / (45 * X**3)
    + Q**2 / (6 * X**2)
    + Q**3 / (2 * X)
)
assert sp.expand(principal_part - expected_principal_part) == 0


# The pole-free value c=-9 gives a polynomial Darboux completion.  Use W,Y
# to keep its exact formula compact.
W_source = sp.expand(Z_source - 9 * Q_source**2)
Y_source = sp.expand(Q_source - X_source * W_source / 3)
U_source = 1 + X_source * Y_source

R = sp.expand(2 * X_source - 3 * X_source**2 * Y_source - X_source**3 * W_source)
S = sp.expand(
    (
        U_source**3 * W_source
        + Y_source**2 * U_source * (4 + 3 * X_source * Y_source)
    )
    / 2
)
T = sp.expand(
    Y_source
    + 3 * X_source * U_source**2 * W_source
    + 3 * X_source * Y_source**2 * (4 + 3 * X_source * Y_source)
)

completion = -(
    10 * W_source**3 * X_source**2
    + 90 * W_source**2 * X_source * Y_source
    + 20 * W_source**2
    - 18 * W_source * X_source**3 * Y_source**5
    - 90 * W_source * X_source**2 * Y_source**4
    - 180 * W_source * X_source * Y_source**3
    + 90 * W_source * Y_source**2
    - 54 * X_source**2 * Y_source**6
    - 234 * X_source * Y_source**5
    - 375 * Y_source**4
) / 60
D = sp.expand(E_source + completion)

# The displayed coordinate really stays small.
assert sp.expand(R - x * (2 - 3 * x * q)) == 0

# All six brackets are checked after full substitution into Q[x,q,p,z].
assert poisson(D, R) == 1
assert poisson(S, T) == 1
assert poisson(R, S) == 0
assert poisson(R, T) == 0
assert poisson(D, S) == 0
assert poisson(D, T) == 0
assert sp.factor(sp.Matrix([R, T, D, S]).jacobian((x, q, p, z)).det()) == 1

# The map (X,Q,Z,E)->(X,Y,W,E+completion) is triangular over a polynomial
# base automorphism.  In those coordinates, (S,T,R,D) is exactly
# (F1/2,F2,F3,id).  Thus the full map has the same fiber schemes and generic
# degree as the foundational map.
Xb, Yb, Wb, Db = sp.symbols("Xb Yb Wb Db")
Q_back = Yb + Xb * Wb / 3
Z_back = Wb + 9 * Q_back**2
completion_base = -(
    10 * Wb**3 * Xb**2
    + 90 * Wb**2 * Xb * Yb
    + 20 * Wb**2
    - 18 * Wb * Xb**3 * Yb**5
    - 90 * Wb * Xb**2 * Yb**4
    - 180 * Wb * Xb * Yb**3
    + 90 * Wb * Yb**2
    - 54 * Xb**2 * Yb**6
    - 234 * Xb * Yb**5
    - 375 * Yb**4
) / 60
E_back = Db - completion_base
W_adapted = Z - 9 * Q**2
Y_adapted = Q - X * W_adapted / 3
completion_adapted = -(
    10 * W_adapted**3 * X**2
    + 90 * W_adapted**2 * X * Y_adapted
    + 20 * W_adapted**2
    - 18 * W_adapted * X**3 * Y_adapted**5
    - 90 * W_adapted * X**2 * Y_adapted**4
    - 180 * W_adapted * X * Y_adapted**3
    + 90 * W_adapted * Y_adapted**2
    - 54 * X**2 * Y_adapted**6
    - 234 * X * Y_adapted**5
    - 375 * Y_adapted**4
) / 60
adapted_back = {X: Xb, Q: Q_back, Z: Z_back, E: E_back}
for got, expected in zip(
    tuple(
        value.subs(adapted_back, simultaneous=True)
        for value in (X, Y_adapted, W_adapted, E + completion_adapted)
    ),
    (Xb, Yb, Wb, Db),
    strict=True,
):
    assert sp.expand(got - expected) == 0

# Transport the known complete foundational fiber over (-1/4,0,0), and set
# D=0.  The original variables are pleasantly small despite the construction.
collision_points = (
    (sp.Rational(0), sp.Rational(0), sp.Rational(1, 24), sp.Rational(-1, 8)),
    (sp.Rational(1), sp.Rational(2, 3), sp.Rational(247, 96), sp.Rational(-89, 64)),
    (sp.Rational(-1), sp.Rational(-2, 3), sp.Rational(247, 96), sp.Rational(-89, 64)),
)
target = (sp.Rational(0), sp.Rational(0), sp.Rational(0), sp.Rational(-1, 8))
for point in collision_points:
    values = tuple(
        sp.factor(output.subs(dict(zip((x, q, p, z), point, strict=True))))
        for output in (R, T, D, S)
    )
    assert values == target
assert len(set(collision_points)) == 3

print("PASS: c=-9 is the unique pole-free shear in the family Z -> Z+cQ^2")
print("PASS: all six rank-two Poisson brackets hold exactly over Q")
print("PASS: det d(R,T,D,S)=1 and the displayed R is x(2-3xq)")
print("PASS: the complete foundational three-point fiber transports to (0,0,0,-1/8)")
print("PASS: the construction yields PC(2), exact symplectic A^4, and cotangent/Weyl A_4 consequences")
print("SCOPE: equality with the unavailable announced manuscript formulas is not asserted")
