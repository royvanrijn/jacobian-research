#!/usr/bin/env python3
"""Exact regression for the elementary transport skeleton shared by
JC(3), Long's GMC(3) witness, and the homogeneous GVC(3) witness.

This is a finite symbolic replay.  The all-m jet statement is proved in the
accompanying analysis by the endpoint-contact argument; here we check m<=6.
"""
from __future__ import annotations

import sympy as sp

p, q, x, y, z, t, u, w = sp.symbols("p q x y z t u w")

# Universal elementary transport frame.
M = sp.Matrix([[1 + p*q, q], [p, 1]])
E12 = sp.Matrix([[1, q], [0, 1]])
E21 = sp.Matrix([[1, 0], [p, 1]])
assert M.det() == 1
assert M == E12 * E21

Ap = sp.simplify(M.inv() * M.diff(p))
Aq = sp.simplify(M.inv() * M.diff(q))
assert Ap == sp.Matrix([[0, 0], [1, 0]])
assert Aq == sp.Matrix([[p, 1], [-p**2, -p]])
assert Ap**2 == sp.zeros(2)
assert Aq**2 == sp.zeros(2)
# Maurer-Cartan coefficient of dp wedge dq.
assert sp.simplify(Aq.diff(p) - Ap.diff(q) + Ap*Aq - Aq*Ap) == sp.zeros(2)

# Foundational JC(3) marked-root identity.
uxy = 1 + x*y
a = uxy**3*z + y**2*uxy*(4 + 3*x*y)
b = y + 3*x*uxy**2*z + 3*x*y**2*(4 + 3*x*y)
c = 2*x - 3*x**2*y - x**3*z
Q = c*uxy**3 - 2*uxy**2*x + b*uxy*x**2 - 2*a*x**3
assert sp.expand(Q) == 0
assert sp.det(sp.Matrix([[uxy, y], [x, 1]])) == 1

# Long/Hopf endpoint divisibility on t^2 + 2xy = 1, with a=1.
U_long = 1 + x
D_long = 1 - t**2*U_long**2
D_long_reduced = sp.expand(D_long.subs(t**2, 1 - 2*x*y))
assert sp.rem(sp.Poly(D_long_reduced, x), sp.Poly(x, x)) == 0

# Homogeneous GVC(3) polynomial lift identity.
rho = t**2 + x*y
A = rho + x**2
C = y*rho**2 - 2*x*t**2*rho - x**3*t**2
assert sp.expand(x*C - (rho**3 - t**2*A**2)) == 0
# On rho=1: P = (U/V)*(1-t^2 U^2)^2 with U=1+x^2, V=x^2.

# Same adjacent jet line for contact exponents s=1 (Long) and s=2 (GVC).
B = 1 + u
for s in (1, 2):
    for m in range(1, 7):
        J = sp.integrate((1 - w**2)**(s*m), (w, 0, B))
        K = sp.Poly(sp.expand(B**(m - 1) * J), u)
        Csm = sp.integrate((1 - w**2)**(s*m), (w, 0, 1))
        assert sp.simplify(K.coeff_monomial(u**m)) == 0
        assert sp.simplify(K.coeff_monomial(u**(m - 1)) - Csm) == 0

# Cohn's determinant-one frame: a flat polynomial frame, but known not to be
# elementary over k[p,q].  This demonstrates that flatness alone is too weak.
Cohn = sp.Matrix([[p**2, p*q - 1], [p*q + 1, q**2]])
assert sp.expand(Cohn.det()) == 1
Cp = sp.simplify(Cohn.inv() * Cohn.diff(p))
Cq = sp.simplify(Cohn.inv() * Cohn.diff(q))
assert sp.simplify(Cq.diff(p) - Cp.diff(q) + Cp*Cq - Cq*Cp) == sp.zeros(2)

print("PASS universal frame: M=E12(q)E21(p), det(M)=1")
print("PASS universal Maurer-Cartan connection is flat; both coordinate components are nilpotent")
print("PASS JC(3): marked root [1+xy:x] is the universal first column")
print("PASS Long/GMC(3): endpoint numerator is divisible by the pole coordinate")
print("PASS GVC(3): homogeneous polynomialization identity holds")
print("PASS Long s=1 and GVC s=2 share the same adjacent jet line for m=1..6")
print("PASS Cohn frame is flat and determinant one, showing flatness alone cannot be the invariant")
