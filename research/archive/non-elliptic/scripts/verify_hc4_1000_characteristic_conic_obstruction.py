#!/usr/bin/env python3
"""Verify the transverse characteristic-conic obstruction in chart 1000.

Use

    u = a-rho*b-sigma,  v=b,
    p = u^2+q(v)

for a smooth geometrically irreducible conic with p_aa != 0.  If the conic
is characteristic for the Hamiltonian field of f_a, then

    f_a-c = p*ell

with ell affine-linear.  The normal-Hessian equation on p=0 is

    ell*D^2 = 2*kappa*p_u,
    D = p_u*(4*v^3-2*g_b)+p_b*(2*g_a-1).

The divisor of p_u on the projective conic consists of the two distinct
points cut out by u=0.  Hence p_u/ell can be a square only if ell=x*u.
This geometric parity step is recorded in the audit; the checker verifies
the exact polynomial consequences.

Caustic divisibility then gives the quartic normal form

    f = c*u + x*p^2/4 + 5*v^3/6 + alpha*v + beta.

On p=0, the remaining normal-Hessian equation becomes x*D^2=4*kappa, so D
must be constant.  For every cubic g, however, the normal remainder of D
modulo p has coefficient 8 on u*v^3.  Thus no transverse quartic
characteristic-conic component survives.
"""

from __future__ import annotations

from itertools import product

import sympy as sp


u, v, rho = sp.symbols("u v rho")
q0, q1, q2 = sp.symbols("q0 q1 q2")
x, c, alpha, beta = sp.symbols("x c alpha beta", nonzero=True)
q = q2 * v**2 + q1 * v + q0
p = u**2 + q


def derivative_a(expression: sp.Expr) -> sp.Expr:
    return sp.diff(expression, u)


def derivative_b(expression: sp.Expr) -> sp.Expr:
    return sp.diff(expression, v) - rho * sp.diff(expression, u)


def remainder_mod_p(expression: sp.Expr) -> sp.Expr:
    return sp.rem(
        sp.Poly(sp.expand(expression), u, domain=sp.EX),
        sp.Poly(p, u, domain=sp.EX),
    ).as_expr()


f = c * u + x * p**2 / 4 + sp.Rational(5, 6) * v**3 + alpha * v + beta
f_a = derivative_a(f)
f_aa = derivative_a(f_a)
f_ab = derivative_b(f_a)
f_bb = derivative_b(derivative_b(f))
L = 5 * v * f_aa - (f_aa * f_bb - f_ab**2)

assert sp.factor(f_a - c - x * u * p) == 0
assert sp.factor(remainder_mod_p(L)) == 0
assert sp.factor(remainder_mod_p(f_aa) + 2 * x * q) == 0

g_exponents = [
    powers
    for powers in product(range(4), repeat=2)
    if sum(powers) <= 3
]
g_coefficients = sp.symbols(f"g0:{len(g_exponents)}")
g = sum(
    coefficient * u**powers[0] * v**powers[1]
    for coefficient, powers in zip(
        g_coefficients, g_exponents, strict=True
    )
)
g_a = derivative_a(g)
g_b = derivative_b(g)
p_a = derivative_a(p)
p_b = derivative_b(p)
D = p_a * (4 * v**3 - 2 * g_b) + p_b * (2 * g_a - 1)
D_remainder = sp.expand(remainder_mod_p(D))

C = 4 * v**3 * f_aa - 2 * f_aa * g_b + 2 * f_ab * g_a - f_ab
assert sp.factor(remainder_mod_p(C - x * u * D)) == 0

D_polynomial = sp.Poly(D_remainder, u, v)
assert D_polynomial.coeff_monomial(u * v**3) == 8

# All g-dependent terms in D have total degree at most three before
# reduction by the monic quadratic p.  The degree-four term is exactly the
# displayed 8*u*v^3 and therefore cannot be cancelled.
degree_four_part = sum(
    coefficient * u**powers[0] * v**powers[1]
    for powers, coefficient in D_polynomial.terms()
    if sum(powers) == 4
)
assert sp.factor(degree_four_part - 8 * u * v**3) == 0


def main() -> None:
    print("PASS: f_a-c=x*u*p and the displayed quartic f has p | L")
    print("PASS: C=x*u*D modulo p")
    print("PASS: the degree-four part of D mod p is exactly 8*u*v^3")
    print("RESULT: D cannot be constant, contradicting x*D^2=4*kappa")
    print("        no transverse quartic characteristic conic survives")


if __name__ == "__main__":
    main()
