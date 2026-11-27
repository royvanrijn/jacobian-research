#!/usr/bin/env python3
"""Replay the centered sextic-quartic identities and leading pole comparisons."""

from sympy import Poly, expand, symbols


Z, a, c, d, h, k = symbols("Z a c d h k")
b0, e0, f0, g0, i0 = symbols("b0 e0 f0 g0 i0")
Dc, Dd, Dh, Dk = symbols("Dc Dd Dh Dk")
b = 3*c/2 + b0
e = 5*a*c/4 + 3*d/2 + e0
f = 3*c**2/8 + b0*c + 5*a*d/4 + 3*h/2 + f0
g = 5*a*c**2/32 + 3*c*d/4 + 3*e0*c/4 + b0*d + 5*a*h/4 + g0
i = 5*a*c*d/16 + b0*h - c**3/16 + c*f0/2 + 3*c*h/4 + 3*d**2/8 + 3*d*e0/4 + i0
P = Z**6 + a*Z**5 + b*Z**4 + e*Z**3 + f*Z**2 + g*Z + i
Q = Z**4 + c*Z**2 + d*Z + h


def derivative(value):
    return value.diff(c)*Dc + value.diff(d)*Dd + value.diff(k)*Dk


bracket = (
    (P.diff(c)*Dc + P.diff(d)*Dd + P.diff(h)*Dh)*Q.diff(Z)
    - P.diff(Z)*(Q.diff(c)*Dc + Q.diff(d)*Dd + Q.diff(h)*Dh)
)
shifted = expand(bracket.subs({h: k+c**2/4, Dh: Dk+c*Dc/2}, simultaneous=True))
L = (-5*a*c**3/32 - 5*a*c*k/4 - 5*a*d**2/8 - 3*e0*c**2/8
     - c*g0 - 2*d*f0 - 3*d*k - 3*e0*k)
M = (5*a*c**2*d/32 - 5*a*d*k/4 + 3*c*d**2/4 + 3*c*d*e0/4
     - d*g0 - 2*f0*k - 3*k**2/2)
N = ((-15*a*c**3/64 - 5*a*c*k/8 + 5*a*d**2/16 - 3*c**2*e0/8
      - c*g0/2 + d*f0/2 + 3*d*k/4)*Dc
     + (5*a*c*d/16 + 3*d**2/4 + 3*d*e0/4)*Dd
     + (-15*a*c**2/32 - 5*a*k/4 - 3*c*e0/4 - g0)*Dk)
assert expand(shifted - derivative(L)*Z**2 - derivative(M)*Z - N) == 0

u = symbols("u")
assert expand((-15*a*c**3/64 + 3*d*k/4).subs(d*k, -5*a*c**3/96)*u*c) == -35*a*u*c**4/128
assert expand((-3*c**2*e0/8 + 3*d*k/4).subs(d*k, -e0*c**2/8)*u*c) == -15*e0*u*c**3/32
assert expand((-g0*c/2 + 3*d*k/4).subs(d*k, -g0*c/3)*u*c) == -3*g0*u*c**2/4

T, DT, A0 = symbols("T DT A0")
last = N.subs({a: 0, e0: 0, g0: 0, k: T-2*f0/3, Dk: DT}, simultaneous=True)
assert expand(expand(last).subs(d*T, A0) - (3*A0*Dc/4 + 3*d**2*Dd/4)) == 0

S = Z**2 + c/2
Q0 = Q.subs(h, c**2/4 - 2*f0/3)
P0 = P.subs({a: 0, e0: 0, g0: 0, h: c**2/4 - 2*f0/3}, simultaneous=True)
assert expand(Q0 - (S**2 + d*Z - 2*f0/3)) == 0
assert expand(P0 - b0*Q0 - (S**3 + 3*d*Z*S/2 + 3*d**2/8 + i0)) == 0
print("Sextic triple line: sextic-quartic bracket and pole factors verified.")
