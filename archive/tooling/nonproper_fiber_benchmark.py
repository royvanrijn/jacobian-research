#!/usr/bin/env python3
"""Exact certificate for the one-parameter nonproper-fiber benchmark."""
import sympy as sp

x, y, z, s, t, w = sp.symbols("x y z s t w")
u = 1 + x*y
F = (
    u**3*z + y**2*u*(4 + 3*x*y),
    y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y),
    2*x - 3*x**2*y - x**3*z,
)
target = (-sp.Rational(1, 4) + s, 0, 0)
J = sp.Matrix(F).jacobian((x, y, z))

# The bounded path, valid for every s.
bounded = {x: 0, y: 0, z: -sp.Rational(1, 4) + s}
assert tuple(sp.expand(f.subs(bounded)) for f in F) == target

# Put s=1/4-t^2.  The signs t>0 and t<0 give the two unbounded paths.
unbounded = {
    x: -1/(2*t),
    y: 3*t,
    z: 26*t**2,
    s: sp.Rational(1, 4) - t**2,
}
assert tuple(sp.factor(f.subs(unbounded)) for f in F) == (-t**2, 0, 0)

# No affine singularity explains the loss of two solutions at the endpoint.
assert sp.factor(J.det()) == -2
J_path = J.subs(unbounded).applyfunc(sp.factor)
assert J_path == sp.Matrix([
    [sp.Rational(9, 2)*t**3, sp.Rational(3, 4)*t, -sp.Rational(1, 8)],
    [sp.Rational(3, 2)*t**2, sp.Rational(25, 4), -sp.Rational(3, 8)/t],
    [-sp.Rational(17, 2), -sp.Rational(3, 4)/t**2, sp.Rational(1, 8)/t**3],
])
assert J_path.inv().applyfunc(sp.factor) == sp.Matrix([
    [-sp.Rational(1, 4)/t**3, 0, -sp.Rational(1, 4)],
    [-sp.Rational(3, 2)/t, sp.Rational(1, 4), -sp.Rational(3, 4)*t**2],
    [-26, sp.Rational(3, 2)*t, -sp.Rational(27, 2)*t**3],
])

# Clearing the pole gives projective coordinates [X:Y:Z:W].  Both signs tend
# to the same point [-1:0:0:0], while their local parameter t has opposite sign.
projective = (-1, 6*t**2, 52*t**3, 2*t)
assert tuple(sp.limit(q, t, 0) for q in projective) == (-1, 0, 0, 0)
assert sp.simplify(projective[0]/projective[3] - unbounded[x]) == 0
assert sp.simplify(projective[1]/projective[3] - unbounded[y]) == 0
assert sp.simplify(projective[2]/projective[3] - unbounded[z]) == 0

print("PASS: one bounded branch and two branches t=+/-sqrt(1/4-s)")
print("PASS: unbounded asymptotics x=-1/(2t), y=3t, z=26t^2")
print("PASS: both unbounded branches limit to [-1:0:0:0]")
print("PASS: det(DF)=-2 on every affine path")
print("PASS: Jacobian conditioning can grow on the order of delta^-3")
