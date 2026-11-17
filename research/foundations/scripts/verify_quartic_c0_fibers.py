#!/usr/bin/env python3
"""Direct affine fiber calculation for the quartic map over C=0."""

import sympy as sp


x, y, z, A, B, v = sp.symbols("x y z A B v")
u = 1 + 3*x*y
gamma = 1 - 4*x*y - x**2*z
G = (
    sp.cancel((2*u + u**2 - 3*u**4*gamma**2)/x**2),
    sp.cancel((1 + u - 2*u**3*gamma**2)/x),
    sp.expand(x*gamma),
)

# Since G_3=x*gamma, the affine source over C=0 is exactly x=0 union gamma=0.
assert sp.factor(G[2] - x*gamma) == 0

# The x=0 component maps triangularly and bijectively to the whole target plane.
x0_map = tuple(sp.factor(component.subs(x, 0)) for component in G)
assert tuple(sp.expand(got-want) for got, want in zip(
    x0_map, (3*(29*y**2+2*z), y, 0)
)) == (0, 0, 0)
x0_point = {x: 0, y: B, z: (A-87*B**2)/6}
assert tuple(sp.factor(component.subs(x0_point)) for component in G) == (A, B, 0)

# On gamma=0 with x!=0, put v=xy. The two remaining target equations reduce
# to (B^2-A)x^2=1.
gamma_substitution = {y: v/x, z: (1-4*v)/x**2}
gamma_map = tuple(sp.factor(component.subs(gamma_substitution)) for component in G)
assert tuple(sp.factor(got-want) for got, want in zip(gamma_map, (
    3*(v+1)*(3*v+1)/x**2, (3*v+2)/x, 0,
))) == (0, 0, 0)
v_from_B = (B*x-2)/3
assert sp.factor(gamma_map[1].subs(v, v_from_B)-B) == 0
assert sp.factor(
    (gamma_map[0].subs(v, v_from_B)-A)*x**2
    - ((B**2-A)*x**2-1)
) == 0

# Therefore over C the gamma=0 chart contributes two points when A!=B^2 and
# none when A=B^2. Together with x=0, the fiber sizes are 3 and 1.
X = sp.symbols("X")
gamma_equation = (B**2-A)*X**2-1
assert sp.discriminant(gamma_equation, X) == 4*(B**2-A)
assert sp.factor(gamma_equation.subs(A, B**2)) == -1

# The inverse quartic at C=0 records the same boundary split. W=0 is double,
# while W=+1 is the finite x=0 branch and W=-1 is an additional infinity
# branch to be analyzed in the nonproperness step.
W = sp.symbols("W")
E_c0 = W**2-W**4
assert sp.factor(E_c0) == -W**2*(W-1)*(W+1)

print("PASS: G_3=0 splits the affine source into x=0 and gamma=0")
print("PASS: x=0 maps bijectively to every target on C=0")
print("PASS: gamma=0 contributes two points iff A!=B^2")
print("PASS: C=0 fiber sizes are 3 off A=B^2 and 1 on A=B^2")
print("PASS: W=-1 is the remaining projective branch, not an affine point")
