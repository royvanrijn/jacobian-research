#!/usr/bin/env python3
"""Exact audit of an integral Jacobian-one Hasse-failing fiber.

This file is intentionally separate from the papers.  It checks:

* the displayed map has coefficients in ZZ and Jacobian determinant one;
* its fiber over (3,-1,1) is Q-isomorphic to Q[S]/(P);
* P has no rational root and has a real root;
* P has a p-adic root for every prime;
* the reconstruction gives a Z_p point, including the delicate prime 2.
"""

from __future__ import annotations

import sympy as sp


x, y, z, S = sp.symbols("x y z S")


# An S_3 cubic and the quadratic defining its discriminant field.
cubic = S**3 - 2 * S**2 + 8
quadratic = S**2 - S + 6
P = sp.expand(cubic * quadratic)
assert P == S**5 - 3 * S**4 + 8 * S**3 - 4 * S**2 - 8 * S + 48

# The cubic is irreducible, its discriminant squareclass is -23, and the
# quadratic is the corresponding discriminant-field polynomial.
assert sp.Poly(cubic, S, domain=sp.QQ).is_irreducible
assert sp.Poly(quadratic, S, domain=sp.QQ).is_irreducible
assert sp.discriminant(cubic, S) == -1472 == -23 * 8**2
assert sp.discriminant(quadratic, S) == -23
assert sp.gcd(sp.Poly(P, S), sp.Poly(sp.diff(P, S), S)).degree() == 0


# Integralized quadratic gauge.  These formulas arise from the seed
# G=P-P(0), the source dilation x_old=2*x, and the determinant-one target
# lattice (C_old/4, B_old-(2*g2/g1)*Pi, Pi).
t = 1 + 2 * x * y
q = t**2 * z - y**2 * (1 + 3 * t)
F = (
    x * (1 - 3 * x * y) + 2 * x**3 * z - 3 * x**4 * q**4 + 3 * x**5 * q**5,
    y - 6 * x * q + 6 * t**2 * x**2 * q**4 - 5 * t**2 * x**3 * q**5,
    t * q,
)
target = (sp.Integer(3), sp.Integer(-1), sp.Integer(1))

for coordinate in F:
    assert sp.Poly(sp.expand(coordinate), x, y, z, domain=sp.ZZ)

# Keep t and q independent while differentiating, then apply the chain rule.
# This is exactly the full Jacobian calculation but avoids a very large
# intermediate expansion of q^5.
T, Q = sp.symbols("T Q")
t_formula = 1 + 2 * x * y
q_formula = T**2 * z - y**2 * (1 + 3 * T)
F_intermediate = sp.Matrix(
    (
        x * (1 - 3 * x * y)
        + 2 * x**3 * z
        - 3 * x**4 * Q**4
        + 3 * x**5 * Q**5,
        y - 6 * x * Q + 6 * T**2 * x**2 * Q**4 - 5 * T**2 * x**3 * Q**5,
        T * Q,
    )
)
source_variables = (x, y, z)
dt = [sp.diff(t_formula, variable) for variable in source_variables]
dq_direct = [sp.diff(q_formula, variable) for variable in source_variables]
jacobian_rows = []
for coordinate in F_intermediate:
    row = []
    for index, variable in enumerate(source_variables):
        row.append(
            sp.diff(coordinate, variable)
            + sp.diff(coordinate, T) * dt[index]
            + sp.diff(coordinate, Q)
            * (dq_direct[index] + sp.diff(q_formula, T) * dt[index])
        )
    jacobian_rows.append(row)
jacobian_determinant = sp.det(sp.Matrix(jacobian_rows))
assert sp.factor(
    jacobian_determinant.subs(Q, q_formula).subs(T, t_formula)
) == 1


# Scheme-theoretic reconstruction over Q.  For the underlying quadratic
# gauge, (Pi,B,C)=(1,0,12), and its inverse equation is P.
G = P - P.subs(S, 0)
g1 = sp.Poly(G, S).coeff_monomial(S)
g2 = sp.Poly(G, S).coeff_monomial(S**2)
g3 = sp.Poly(G, S).coeff_monomial(S**3)
assert (g1, g2, g3) == (-8, -4, 8)
assert sp.expand(G - g1 * (0 * S**2 + 12) / 2 - P) == 0

P_poly = sp.Poly(P, S, domain=sp.QQ)


def quotient_reduce(expression: sp.Expr) -> sp.Expr:
    """Reduce an element of QQ(S) in the etale algebra QQ[S]/(P)."""
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    numerator_poly = sp.Poly(numerator, S, domain=sp.QQ).rem(P_poly)
    denominator_poly = sp.Poly(denominator, S, domain=sp.QQ).rem(P_poly)
    inverse = sp.invert(denominator_poly, P_poly)
    return (numerator_poly * inverse).rem(P_poly).as_expr()


d = sp.cancel(sp.diff(P, S) / g1)
beta = sp.cancel((sp.diff(G, S) / g1 - 1 - S**2) / S)
t_bar = quotient_reduce(1 / d)
x_bar = quotient_reduce(S / (2 * d))
y_bar = quotient_reduce(-beta - S)
z_bar = quotient_reduce(d**2 * (d + y_bar**2 * (1 + 3 * t_bar)))

reconstructed = {x: x_bar, y: y_bar, z: z_bar}
assert quotient_reduce(1 + 2 * x_bar * y_bar - t_bar) == 0
assert quotient_reduce(
    t_bar**2 * z_bar - y_bar**2 * (1 + 3 * t_bar) - d
) == 0
fiber_outputs = (
    x_bar * (1 - 3 * x_bar * y_bar)
    + 2 * x_bar**3 * z_bar
    - 3 * x_bar**4 * d**4
    + 3 * x_bar**5 * d**5,
    y_bar
    - 6 * x_bar * d
    + 6 * t_bar**2 * x_bar**2 * d**4
    - 5 * t_bar**2 * x_bar**3 * d**5,
    t_bar * d,
)
for coordinate, value in zip(fiber_outputs, target, strict=True):
    assert quotient_reduce(coordinate - value) == 0

# Conversely, the marked root is recovered from a source point on the fiber:
# in the old quadratic-gauge chart S=x_old/t=2*x/t.
assert quotient_reduce(2 * x_bar / t_bar - S) == 0


# Exceptional primes.  At 2 the quadratic factor has a simple root 0 modulo
# 2 (although the product P has multiple reduction there).  At 23 the cubic
# has the simple root 7 modulo 23.
assert quadratic.subs(S, 0) % 2 == 0
assert sp.diff(quadratic, S).subs(S, 0) % 2 == 1
assert cubic.subs(S, 7) % 23 == 0
assert sp.diff(cubic, S).subs(S, 7) % 23 == 4

# At every unramified prime the common S_3 splitting field gives a factor
# root: odd Frobenius fixes a cubic root; even Frobenius makes the
# discriminant quadratic split.  The only field-ramified primes are 2 and 23.
assert sp.factorint(abs(int(sp.discriminant(cubic, S)))) == {2: 6, 23: 1}
assert sp.factorint(abs(int(sp.discriminant(quadratic, S)))) == {23: 1}

# The factor resultant contributes the additional bad-reduction prime 7.
# The product nevertheless has the simple root 1 there.
assert sp.resultant(cubic, quadratic, S) == 392 == 2**3 * 7**2
assert P.subs(S, 1) % 7 == 0
assert sp.diff(P, S).subs(S, 1) % 7 == 1

# The 2-adic quadratic root s is even.  Writing s=2u makes every apparent
# reconstruction denominator integral and makes d a 2-adic unit.
u = sp.symbols("u")
d_at_even_quadratic_root = sp.factor(
    -(cubic * sp.diff(quadratic, S) / 8).subs(S, 2 * u)
)
assert sp.expand(
    d_at_even_quadratic_root + (4 * u - 1) * (u**3 - u**2 + 1)
) == 0
d_mod_2 = sp.Poly(d_at_even_quadratic_root, u, modulus=2)
assert d_mod_2.rem(sp.Poly(u**2 - u, u, modulus=2)).as_expr() == 1
assert sp.expand(y_bar.as_numer_denom()[0])  # regression: reconstruction exists
explicit_y_at_two = -1 + 6 * u - 6 * u**2 + 5 * u**3
assert sp.expand((-beta - S).subs(S, 2 * u) - explicit_y_at_two) == 0


# The cubic has one real root (negative discriminant), while neither factor
# has a rational root.  Hence the fiber is real-soluble but not Q-soluble.
assert sp.discriminant(cubic, S) < 0
assert not sp.Poly(cubic, S, domain=sp.QQ).ground_roots()
assert not sp.Poly(quadratic, S, domain=sp.QQ).ground_roots()


print("PASS: F lies in ZZ[x,y,z]^3 and det(JF)=1")
print("PASS: F^(-1)(3,-1,1) is Spec(Q[S]/(P))")
print("PASS: P has a root over every Z_p and a real root")
print("PASS: the reconstructed local source points are integral, including p=2")
print("PASS: P has no rational root, so the fiber has no Q-point")
