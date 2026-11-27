#!/usr/bin/env python3
"""Replay the centered quintic-cubic bracket identities."""

from sympy import Poly, expand, symbols


s, a, b, c, d, e, f, g = symbols("s a b c d e f g")
Da, Db, Dc, Dd, De, Df, Dg = symbols("Da Db Dc Dd De Df Dg")
b0, e0, f0 = symbols("b0 e0 f0")
P = s**5 + a*s**4 + b*s**3 + e*s**2 + f*s + g
Q = s**3 + c*s + d
DP = Da*s**4 + Db*s**3 + De*s**2 + Df*s + Dg
DQ = Dc*s + Dd
bracket = Poly(DP*Q.diff(s) - P.diff(s)*DQ, s)


def derivation(polynomial):
    return polynomial.diff(c)*Dc + polynomial.diff(d)*Dd


values = {
    b: 5*c/3 + b0,
    e: 4*a*c/3 + 5*d/3 + e0,
    f: 5*c**2/9 + b0*c + 4*a*d/3 + f0,
    g: 2*a*c**2/9 + 10*c*d/9 + 2*e0*c/3 + b0*d,
}
substitution = dict(values)
substitution[Da] = 0
for derivative, coefficient in zip((Db, De, Df, Dg), (b, e, f, g), strict=True):
    substitution[derivative] = derivation(values[coefficient])

remainder = bracket.as_expr().subs(substitution, simultaneous=True)
L = 5*c**3 - 36*a*c*d - 27*f0*c - 45*d**2 - 54*e0*d
M = 4*a*c**3 + 15*c**2*d + 9*e0*c**2 - 18*a*d**2 - 27*f0*d
assert expand(remainder - s*derivation(L)/27 - derivation(M)/27) == 0

print("Sextic triple-line cubic mate: exact quintic-cubic bracket verified.")
