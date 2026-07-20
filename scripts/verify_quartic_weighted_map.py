#!/usr/bin/env python3
"""Exact quartic-sheet weighted map and its C!=0 root reconstruction."""

import sympy as sp


x, y, z = sp.symbols("x y z")
A, B, C, W = sp.symbols("A B C W")

u = 1 + 3*x*y
gamma = 1 - 4*x*y - x**2*z
w_source = u*gamma

# Integer normalization of the seed p(w)=w-2w^3, H=(w^2-w^4)/2.
G = (
    sp.cancel((2*u + u**2 - 3*u**4*gamma**2)/x**2),
    sp.cancel((1 + u - 2*u**3*gamma**2)/x),
    sp.expand(x*gamma),
)
assert all(sp.denom(component) == 1 for component in G)

determinant = sp.factor(sp.Matrix(G).jacobian((x, y, z)).det())
degrees = tuple(sp.Poly(component, x, y, z).total_degree() for component in G)
assert determinant == -6
assert degrees == (12, 11, 4)

points = ((1, 0, 0), (-1, 0, 2))
target = (0, 0, 1)
for point in points:
    image = tuple(sp.expand(component.subs(dict(zip((x, y, z), point)))) for component in G)
    assert image == target
assert points[0] != points[1]


# The inverse polynomial has degree four:
#
#   E(W)=W^2-W^4-2BCW+AC^2.
#
# On the source, E(w)=0 and E'(w)=-2 gamma=-2C/x.
E = W**2 - W**4 - 2*B*C*W + A*C**2
dE = sp.diff(E, W)
assert sp.factor(E.subs({A: G[0], B: G[1], C: G[2], W: w_source})) == 0
assert sp.factor(dE.subs({A: G[0], B: G[1], C: G[2], W: w_source}) + 2*gamma) == 0


# Conversely, for C!=0 every simple root reconstructs exactly one source point.
gamma_root = B*C - W + 2*W**3
assert sp.factor(gamma_root + dE/2) == 0
x_root = C/gamma_root
u_root = W/gamma_root
v_root = (u_root - 1)/3
y_root = v_root/x_root
s_root = 1 - 4*v_root - gamma_root
z_root = s_root/x_root**2

reconstruction = {x: x_root, y: y_root, z: z_root}
for got, want in zip(G, (A, B, C)):
    numerator = sp.factor(sp.together(got.subs(reconstruction) - want).as_numer_denom()[0])
    remainder = sp.rem(sp.Poly(numerator, W), sp.Poly(E, W)).as_expr()
    assert sp.factor(remainder) == 0

# These identities recover the construction coordinates, proving uniqueness.
assert sp.factor((1 + 3*x_root*y_root) - u_root) == 0
assert sp.factor((1 - 4*x_root*y_root - x_root**2*z_root) - gamma_root) == 0
assert sp.factor(u_root*gamma_root - W) == 0

# No hidden reconstruction denominator occurs on C!=0: gamma_root=-E'/2, so
# the only excluded roots are exactly the repeated roots E'=0.
assert sp.factor(x_root + 2*C/dE) == 0


# A concrete target has four distinct roots, certifying generic degree four.
# At (A,B,C)=(1,0,1), E=-(W^4-W^2-1), which is squarefree.
quartic_control = sp.factor(E.subs({A: 1, B: 0, C: 1}))
assert quartic_control == -W**4 + W**2 + 1
assert sp.gcd(sp.Poly(quartic_control, W), sp.Poly(sp.diff(quartic_control, W), W)).degree() == 0

print("PASS: quartic weighted map is polynomial with degrees", degrees)
print("PASS: det(DG) =", determinant, "and the stored collision is exact")
print("PASS: E(W)=W^2-W^4-2BCW+AC^2 is the inverse quartic")
print("PASS: every simple root for C!=0 reconstructs uniquely")
print("PASS: repeated roots are exactly the C!=0 reconstruction poles")
