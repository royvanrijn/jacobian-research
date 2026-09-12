#!/usr/bin/env python3
"""Exact pre-audit for the supplied rank-two Poisson abstract.

This script does not reconstruct the missing manuscript polynomials.  It
certifies the exact relation between the one displayed output R and the
foundational three-variable Keller map, and proves that the most immediate
candidate for the other two invariant outputs has no polynomial Darboux
completion.  Thus it records useful evidence without attributing formulas
that have not been obtained from the manuscript.
"""

import sympy as sp


x, q, p, z = sp.symbols("x q p z")
X, Y, Z, E = sp.symbols("X Y Z E")


def poisson(f: sp.Expr, g: sp.Expr) -> sp.Expr:
    """Bracket convention from the abstract: {p,x}={z,q}=1."""

    return sp.expand(
        sp.diff(f, p) * sp.diff(g, x)
        - sp.diff(f, x) * sp.diff(g, p)
        + sp.diff(f, z) * sp.diff(g, q)
        - sp.diff(f, q) * sp.diff(g, z)
    )


# The only output printed in the supplied abstract.
R = x * (2 - 3 * x * q)

# A polynomial source coordinate system adapted to its Hamiltonian flow.
Z_source = 3 * x**2 * p + (2 - 6 * x * q) * z
Y_source = q - x * Z_source / 3
E_source = (1 + 3 * x * q) * p / 2 - 3 * q**2 * z
X_source = x

source_coordinates = sp.Matrix([X_source, Y_source, Z_source, E_source])
assert sp.factor(source_coordinates.jacobian((x, q, p, z)).det()) == -1

# Explicit polynomial inverse.  This proves automorphy, rather than relying
# only on the constant Jacobian determinant.
q_inverse = Y + X * Z / 3
p_inverse = -6 * E * q_inverse * X + 2 * E + 3 * Z * q_inverse**2
z_inverse = -3 * E * X**2 + 3 * Z * q_inverse * X / 2 + Z / 2
inverse = {x: X, q: q_inverse, p: p_inverse, z: z_inverse}
for got, expected in zip(
    source_coordinates.subs(inverse, simultaneous=True), (X, Y, Z, E), strict=True
):
    assert sp.expand(got - expected) == 0

inverse_coordinates = sp.Matrix([X, q_inverse, p_inverse, z_inverse])
for got, expected in zip(
    inverse_coordinates.subs(
        {X: X_source, Y: Y_source, Z: Z_source, E: E_source},
        simultaneous=True,
    ),
    (x, q, p, z),
    strict=True,
):
    assert sp.expand(got - expected) == 0

# The repository's foundational determinant -2 map, in the new invariant
# variables.  Its third component is exactly the displayed R.
u = 1 + X * Y
F1 = sp.expand(u**3 * Z + Y**2 * u * (4 + 3 * X * Y))
F2 = sp.expand(Y + 3 * X * u**2 * Z + 3 * X * Y**2 * (4 + 3 * X * Y))
F3 = sp.expand(2 * X - 3 * X**2 * Y - X**3 * Z)
assert sp.expand(F3.subs({X: X_source, Y: Y_source, Z: Z_source}) - R) == 0

# The invariant quotient has its Nambu bracket.  Consequently the obvious
# scaled foundational outputs form one canonical pair and commute with R.
S0 = sp.expand((F1 / 2).subs({X: X_source, Y: Y_source, Z: Z_source}))
T0 = sp.expand(F2.subs({X: X_source, Y: Y_source, Z: Z_source}))
assert poisson(S0, T0) == 1
assert poisson(R, S0) == 0
assert poisson(R, T0) == 0
assert poisson(E_source, R) == 1

# They do not already give the four outputs claimed in the abstract.
assert poisson(E_source, S0) != 0
assert poisson(E_source, T0) != 0


# Prove that no correction D=E_source+f(X,Y,Z) can make this particular
# pair (S0,T0) into a full polynomial Darboux system.  Since the displayed
# coordinate change is an automorphism and {X,R}={Y,R}={Z,R}=0,
# every polynomial D with {D,R}=1 is necessarily of this form.
variables = (X, Y, Z)
S = F1 / 2
T = F2
r = F3
J = sp.Matrix([S, T, r]).jacobian(variables)
assert sp.factor(J.det()) == -1

# w is the unique vector field on the invariant quotient satisfying
# w(S)=w(T)=0 and w(r)=1.  w0 is the vector field {E_source,-} written in
# (X,Y,Z).  A correction f would have grad(f) x grad(r) = w-w0.
w = (J.adjugate() / J.det())[:, 2].applyfunc(sp.expand)
w0 = sp.Matrix(
    [
        (X**2 * Z + 3 * X * Y + 1) / 2,
        -(6 * X**2 * Z**2 + 24 * X * Y * Z + 18 * Y**2 + Z) / 6,
        3 * Z * (X * Z + 3 * Y) / 2,
    ]
)
difference = (w - w0).applyfunc(sp.expand)
grad_r = sp.Matrix([sp.diff(r, variable) for variable in variables])
assert sp.expand(difference.dot(grad_r)) == 0
assert sp.expand(
    sum(sp.diff(difference[index], variable) for index, variable in enumerate(variables))
) == 0

# The first cross-product equation is
#   (3 d/dZ - X d/dY)f = difference_X/X^2.
# Put Q=Y+XZ/3.  Integrating at fixed (X,Q) gives a particular f0; every
# possible correction is f0+h(X,Q).
Q = sp.symbols("Q")
first_rhs = sp.cancel(difference[0] / X**2)
first_rhs_at_q = sp.expand(first_rhs.subs(Y, Q - X * Z / 3))
f0_at_q = sp.integrate(first_rhs_at_q / 3, Z)
f0 = sp.expand(f0_at_q.subs(Q, Y + X * Z / 3))
grad_f0 = sp.Matrix([sp.diff(f0, variable) for variable in variables])
residual = (difference - grad_f0.cross(grad_r)).applyfunc(sp.expand)
assert residual[0] == 0

# The two remaining equations are equivalent to L(h)=obstruction, where
# L=-3X^2 d/dX+(6XQ-2)d/dQ.  The obstruction is recorded in compact form.
obstruction = sp.factor(residual[2].subs(Y, Q - X * Z / 3))
expected_obstruction = (
    Q**3 * (54 * Q**3 * X**3 + 189 * Q**2 * X**2 + 222 * Q * X + 89) / 2
)
assert sp.expand(obstruction - expected_obstruction) == 0
assert sp.expand(
    residual[1].subs(Y, Q - X * Z / 3) + X * obstruction / 3
) == 0

# In k[X,X^-1,Q], put v=1/X and rho=2X-3X^2Q.  Then L(v)=3 and
# L(rho)=0.  Hence every localized solution is H(v,rho)+C(rho), with
# C polynomial, and the following H is the forced antiderivative.
v, rho = sp.symbols("v rho")
L = lambda expression: sp.expand(
    -3 * X**2 * sp.diff(expression, X)
    + (6 * X * Q - 2) * sp.diff(expression, Q)
)
assert L(1 / X) == 3
assert L(2 * X - 3 * X**2 * Q) == 0

obstruction_vrho = sp.factor(
    expected_obstruction.subs({X: 1 / v, Q: v * (2 - rho * v) / 3})
)
H = sp.expand(sp.integrate(obstruction_vrho / 3, v))
assert sp.expand(3 * sp.diff(H, v) - obstruction_vrho) == 0

H_back = sp.factor(H.subs({v: 1 / X, rho: 2 * X - 3 * X**2 * Q}))
expected_H_back = (
    54 * Q**6 * X**6
    + 234 * Q**5 * X**5
    + 375 * Q**4 * X**4
    + 270 * Q**3 * X**3
    + 90 * Q**2 * X**2
    + 24 * Q * X
    + 4
) / (60 * X**4)
assert sp.expand(H_back - expected_H_back) == 0

# C(rho) has no negative X powers after rho=2X-3X^2Q, whereas H_back has
# the nonzero term 1/(15X^4).  It therefore cannot be made polynomial.
assert sp.limit(X**4 * H_back, X, 0) == sp.Rational(1, 15)

print("PASS: the displayed R is exactly foundational F3 after a polynomial source automorphism")
print("PASS: the obvious invariant outputs satisfy {S0,T0}=1 and commute with R")
print("PASS: the localized obstruction proves that this naive pair has no polynomial D-completion")
print("SCOPE: T,D,S from the announced manuscript remain unavailable and are not reconstructed here")
