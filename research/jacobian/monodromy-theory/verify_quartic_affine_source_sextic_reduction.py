#!/usr/bin/env python3
"""Check the septic witness and finite support bounds in the triple-line case."""

from sympy import Matrix, Poly, Rational, factor, symbols


def quadratic_options(degree):
    options = set()
    for r_degree in (1, 2):
        for r_order in range(1, r_degree + 1):
            for j in range(-r_order, degree + 1):
                if r_order + 2*j >= 0 and 4*r_degree + 2*j <= degree:
                    options.add((r_degree, r_order, j))
    return options


for degree in (6, 7):
    assert quadratic_options(degree) == {(1, 1, 0), (1, 1, 1), (2, 2, -1)}


x, y, z, v, r = symbols("x y z v r")
u = 1 + x*y
P = u**3*z + y**2*u*(4 + 3*x*y)
Q = y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y)
R = 2*x - 3*x**2*y - x**3*z
F = Matrix([P, Q, R])
assert [Poly(f, x, y, z).total_degree() for f in F] == [7, 6, 4]
assert factor(F.jacobian((x, y, z)).det()) == -2
for point in [(0, 0, -Rational(1, 4)), (1, -Rational(3, 2), Rational(13, 2)),
              (-1, Rational(3, 2), Rational(13, 2))]:
    assert F.subs(dict(zip((x, y, z), point))) == Matrix([-Rational(1, 4), 0, 0])

shear = {y: -2*(v - x*y)/3, z: -2*y}
g = Poly(P.subs(shear, simultaneous=True), y)
h = Poly(Q.subs(shear, simultaneous=True), y)
assert (g.degree(), h.degree()) == (3, 2)
assert factor(g.LC().subs(v, (r - x)/x**2)) == -16*r*x**3/27
assert factor(h.LC().subs(v, (r - x)/x**2)) == -8*r*x**2/3

print("Quartic affine source: quadratic-mate bounds and septic witness verified.")
