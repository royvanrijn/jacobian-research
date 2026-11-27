#!/usr/bin/env python3
"""Replay the odd-square-root truncation bracket in the tangent-conic (9,6) case."""

from sympy import Poly, expand, symbols


Z, a, b, c, d, e = symbols("Z a b c d e")
Da, Db, Dc, Dd, De = symbols("Da Db Dc Dd De")
S = Z**3 + a*Z + b
U = c*Z**2 + d*Z + e
P = S**3 + 3*S*U/2 + 3*(c**2*Z + 2*c*d)/8
Q = S**2 + U


def derivation(polynomial):
    return sum(
        polynomial.diff(coefficient)*derivative
        for coefficient, derivative in zip((a, b, c, d, e), (Da, Db, Dc, Dd, De))
    )


bracket = Poly(expand(derivation(P)*Q.diff(Z) - P.diff(Z)*derivation(Q)), Z)
I = -a*c**2 + 2*c*e + d**2
J = -2*a*c*d - b*c**2 + 2*d*e
assert bracket.degree() == 4
assert expand(bracket.nth(4) + 9*derivation(I)/4) == 0
assert expand(bracket.nth(3) + 9*derivation(J)/4) == 0

u, v, Du, Dv = symbols("u v Du Dv")
normalize = {
    a: 2*v + u**2,
    b: -2*u*(v + u**2),
    d: c*u,
    e: c*v,
    Da: 2*Dv + 2*u*Du,
    Db: -2*u*Dv - 2*(v + 3*u**2)*Du,
    Dd: Dc*u + c*Du,
    De: Dc*v + c*Dv,
}
assert expand(bracket.nth(1).subs(normalize, simultaneous=True)
              - 9*c**2*(c*Du + 3*u*Dc)/8) == 0
on_d_zero = {a: 2*e/c, b: 0, d: 0, Db: 0, Dd: 0,
             Da: 2*De/c - 2*e*Dc/c**2}
assert expand(bracket.nth(2).subs(on_d_zero, simultaneous=True)
              - 9*c**2*Dc/8 + 9*e*De/2) == 0
assert expand(bracket.nth(0).subs(on_d_zero, simultaneous=True)
              + 3*(c**2 + 8*e**2/c)*De/8) == 0

print("Tangent-conic (9,6): exact Laurent bracket identities verified.")
