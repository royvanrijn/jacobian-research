#!/usr/bin/env python3
"""Exact finite certificate for the announced C^3 Keller counterexample."""
import sympy as sp

x, y, z = sp.symbols("x y z")
u = 1 + x*y
F = sp.Matrix([
    u**3*z + y**2*u*(4 + 3*x*y),
    y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y),
    2*x - 3*x**2*y - x**3*z,
])

determinant = sp.factor(F.jacobian((x, y, z)).det())
assert determinant == -2

points = [
    (0, 0, -sp.Rational(1, 4)),
    (1, -sp.Rational(3, 2), sp.Rational(13, 2)),
    (-1, sp.Rational(3, 2), sp.Rational(13, 2)),
]
target = sp.Matrix([-sp.Rational(1, 4), 0, 0])
for point in points:
    image = sp.simplify(F.subs(dict(zip((x, y, z), point))))
    assert image == target, (point, image)

degrees = tuple(sp.Poly(f, x, y, z).total_degree() for f in F)
assert degrees == (7, 6, 4)

# Exact Palais--Smale curve for L=16*||F-target||^2.
R = sp.symbols("R", nonzero=True, real=True)
escape_point = {x: R, y: -1/R, z: 5/R**2}
L = 16 * ((F - target).dot(F - target))
gradient_L = sp.Matrix([sp.diff(L, variable) for variable in (x, y, z)])
assert sp.simplify(F.subs(escape_point)) == sp.Matrix([0, 2/R, 0])
assert sp.simplify(L.subs(escape_point) - (1 + 64/R**2)) == 0
assert sp.simplify(gradient_L.subs(escape_point)) == sp.Matrix([
    -392/R**3,
    264/R,
    0,
])

print("PASS: det(JF) =", determinant)
print("PASS: all three rational points map to", tuple(target))
print("PASS: component total degrees =", degrees)
print("PASS: explicit Palais--Smale curve tends to level 1 at infinity")
