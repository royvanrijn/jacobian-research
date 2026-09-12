#!/usr/bin/env python3
"""Verify the primitive cubic, rational reconstruction, and discriminant."""
import sympy as sp

x, y, z, A, B, C, T = sp.symbols("x y z A B C T")
u = 1 + x*y
a = u**3*z + y**2*u*(4 + 3*x*y)
b = y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y)
c = 2*x - 3*x**2*y - x**3*z
t = y + 1/x

assert sp.factor(b - (4*t + 2/x - 3*c*t**2)) == 0
assert sp.factor(2*a - (c*t**3 - 2*t**2 + b*t)) == 0
r = sp.factor(3*c*t**2 - 4*t + b)
assert sp.factor(r - 2/x) == 0

x_rec = 2/r
y_rec = t - r/2
z_rec = sp.Rational(5, 4)*r**2 - sp.Rational(3, 2)*t*r - c*r**3/8
assert sp.factor(x_rec - x) == 0
assert sp.factor(y_rec - y) == 0
assert sp.factor(z_rec - z) == 0

P = C*T**3 - 2*T**2 + B*T - 2*A
Q = 27*A**2*C**2 - 18*A*B*C + 16*A + B**3*C - B**2
assert sp.factor(sp.discriminant(P, T) + 4*Q) == 0
assert sp.factor(sp.discriminant(Q, A) + 4*(3*B*C - 4)**3) == 0
print("PASS: cubic identities and reconstruction")
print("PASS: Disc_T(P) = -4 Q and Disc_A(Q) = -4(3BC-4)^3")

