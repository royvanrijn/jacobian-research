#!/usr/bin/env python3
"""Replay centered quartic-quintic bracket coefficients and leading pole factors."""

from sympy import Poly, Rational, expand, symbols


Z, a, b, e, f, g, c, d, h = symbols("Z a b e f g c d h")
b0, e0, f0, g0, Dc, Dd, Dh = symbols("b0 e0 f0 g0 Dc Dd Dh")
b = 5*c/4 + b0
e = a*c + 5*d/4 + e0
f = 5*c**2/32 + 3*b0*c/4 + a*d + 5*h/4 + f0
g = 5*c*d/16 + e0*c/2 + 3*b0*d/4 + a*h + g0
P = Z**5 + a*Z**4 + b*Z**3 + e*Z**2 + f*Z + g
Q = Z**4 + c*Z**2 + d*Z + h


def derivative(value):
    return sum(value.diff(var)*dv for var, dv in ((c, Dc), (d, Dd), (h, Dh)))


bracket = Poly(derivative(P)*Q.diff(Z) - P.diff(Z)*derivative(Q), Z)
L2 = 3*b0*c**2/8 - 3*b0*h + 5*c**3/32 - c*(4*f0 + 5*h)/4 - 5*d**2/8 - 2*d*e0
L1 = 3*b0*c*d/4 + c**2*(15*d + 16*e0)/32 - d*(4*f0 + 5*h)/4 - 2*e0*h
C0 = ((5*d**2/16 + d*e0/2)*Dc + (3*b0*d/4 + 5*c*d/16)*Dd
      - (3*b0*c/4 + 5*c**2/32 + f0 + 5*h/4)*Dh)
assert expand(bracket.as_expr() - derivative(L2)*Z**2 - derivative(L1)*Z - C0) == 0

u = symbols("u")
h_pole = 3*c**2/8
dominant = (5*d**2*u*c/16 + 5*c*d*(3*u*d/2)/16
            - (5*c**2/32 + 5*h_pole/4)*(2*u*h_pole))
assert expand(dominant.subs(d**2, -c**3/2)) == -55*u*c**4/64
h_pole = c**2/8
assert expand(-(5*c**2/32 + 5*h_pole/4)*(2*u*h_pole)) == -5*u*c**4/64

print("Sextic triple line: quartic-quintic bracket and pole factors verified.")
