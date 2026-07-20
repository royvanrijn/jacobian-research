#!/usr/bin/env python3
"""Chart-complete converse: the quartic map is proper off V(C*Q4)."""

import sympy as sp


A, B, C, W, S = sp.symbols("A B C W S")
E = W**2-W**4-2*B*C*W+A*C**2
dE = sp.diff(E, W)
Q4 = (
    A-B**2
    + C**2*(27*B**4-36*A*B**2+8*A**2)
    + 16*A**3*C**4
)

# On C!=0, the discriminant is nonzero exactly off Q4=0.
assert sp.factor(sp.discriminant(E, W) + 16*C**2*Q4) == 0
assert sp.factor(sp.resultant(E, dE, W) + sp.discriminant(E, W)) == 0

# The projective inverse quartic has no root at infinity: its leading
# coefficient is the constant -1, independent of the target.
E_homogeneous = W**2*S**2-W**4-2*B*C*W*S**3+A*C**2*S**4
assert E_homogeneous.subs({W: 1, S: 0}) == -1

# Every finite root off the discriminant has dE!=0. The complete source
# reconstruction has no denominator beyond C and dE.
gamma = -dE/2
x = C/gamma
u = W/gamma
v = (u-1)/3
y = v/x
s_source = 1-4*v-gamma
z = s_source/x**2
denominators = tuple(
    sp.factor(sp.together(expression).as_numer_denom()[1])
    for expression in (x, y, z)
)
assert sp.factor(denominators[0]+dE/2) == 0
assert denominators[1] == 3*C
assert denominators[2] == 3*C**2

# Algebraically verify the reconstruction target identities modulo E.
x0, y0, z0 = sp.symbols("x0 y0 z0")
uu = 1+3*x0*y0
gg = 1-4*x0*y0-x0**2*z0
G = (
    sp.cancel((2*uu+uu**2-3*uu**4*gg**2)/x0**2),
    sp.cancel((1+uu-2*uu**3*gg**2)/x0),
    sp.expand(x0*gg),
)
for got, want in zip(G, (A, B, C)):
    numerator = sp.together(got.subs({x0: x, y0: y, z0: z})-want).as_numer_denom()[0]
    assert sp.factor(sp.rem(sp.Poly(numerator, W), sp.Poly(E, W)).as_expr()) == 0

print("PASS: off C*Q4=0 the inverse quartic has four finite simple roots")
print("PASS: the projective inverse has no root-at-infinity chart")
print("PASS: reconstruction denominators are units off C*Q4=0")
print("PASS: the quartic map is proper outside V(C) union V(Q4)")
