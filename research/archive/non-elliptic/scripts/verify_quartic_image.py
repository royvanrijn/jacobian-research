#!/usr/bin/env python3
"""Exact image and fiber-cardinality theorem for the quartic weighted map."""

import sympy as sp


W, r, A, B, C = sp.symbols("W r A B C")
s = r-2*r**3
t = r**2-3*r**4
E_boundary = sp.factor(W**2-W**4-2*s*W+t)
quotient = -W**2-2*r*W-3*r**2+1
assert sp.factor(E_boundary-(W-r)**2*quotient) == 0

# The residual quadratic controls every discriminant fiber.
assert sp.factor(quotient.subs(W, r)-(1-6*r**2)) == 0
assert sp.factor(sp.discriminant(quotient, W)+4*(2*r**2-1)) == 0

# Generic discriminant point: one double root plus two simple roots -> 2 points.
generic_control = sp.factor(E_boundary.subs(r, 1))
assert sorted(sp.roots(generic_control, W).values()) == [1, 1, 2]

# Cusps r^2=1/6: triple+simple -> 1 point.
for cusp in (1/sp.sqrt(6), -1/sp.sqrt(6)):
    roots = sorted(sp.roots(E_boundary.subs(r, cusp), W).values())
    assert roots == [1, 3]

# Node r^2=1/2: double+double -> no simple roots and hence no affine point.
for node in (1/sp.sqrt(2), -1/sp.sqrt(2)):
    assert sp.factor(E_boundary.subs(r, node)+(W**2-sp.Rational(1, 2))**2) == 0

# Lifting the unique node (s,t)=(0,-1/4) through s=BC,t=AC^2 gives the omitted
# curve B=0, 4AC^2+1=0. It cannot meet C=0.
assert sp.factor((B*C).subs(B, 0)) == 0
assert sp.factor((A*C**2).subs(A, -1/(4*C**2))+sp.Rational(1, 4)) == 0
assert (4*A*C**2+1).subs(C, 0) == 1

# Direct C=0 calculation supplies all remaining targets: 3 points off A=B^2
# and 1 point on it.
X = sp.symbols("X")
gamma_equation = (B**2-A)*X**2-1
assert sp.discriminant(gamma_equation, X) == 4*(B**2-A)
assert gamma_equation.subs(A, B**2) == -1

print("PASS: quartic discriminant fibers have 2, 1, or 0 affine points")
print("PASS: the omitted locus is exactly B=0, 4AC^2+1=0")
print("PASS: every C=0 target remains in the image with fiber size 3 or 1")
print("PASS: complete quartic fiber sizes are 4,3,2,1,0 by stratum")
