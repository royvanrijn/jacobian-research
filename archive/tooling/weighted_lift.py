#!/usr/bin/env python3
"""Construct/check the quadratic weighted lift behind the announced example."""
import sympy as sp

x, y, z, w = sp.symbols("x y z w")
p = 2*w - 3*w**2
c = sp.Integer(1)
assert p.subs(w, 0) == 0 and p.subs(w, 1) == -c
assert sp.integrate(p, (w, 0, 1)) == 0

q = sp.integrate(w*sp.diff(p, w)/c, w)  # integration constant is zero
kappa = sp.diff(p, w).subs(w, 1)/c
a0 = -(1 + kappa)/(2 + kappa)
b0 = -sp.Rational(1, 2)
v, s = x*y, x**2*z
u = 1 + v
gamma = 1 + a0*v + b0*s
W = u*gamma
alpha = sp.cancel(u + q.subs(w, W)/gamma**2)
beta = sp.cancel(c + p.subs(w, W)/gamma)
coords = [sp.cancel(alpha/x**2), sp.cancel(beta/x), sp.cancel(x*gamma)]

for f in coords:
    _, den = sp.fraction(f)
    assert sp.expand(den) == 1, f"uncancelled denominator: {den}"
det = sp.factor(sp.Matrix(coords).jacobian((x, y, z)).det())
assert det == b0*c
print("PASS: endpoint and integral constraints")
print("PASS: all apparent poles cancel")
print("PASS: weighted lift determinant =", det)
print("component degrees =", tuple(sp.Poly(f, x, y, z).total_degree() for f in coords))

