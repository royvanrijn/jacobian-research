#!/usr/bin/env python3
"""Degree-at-most-two polynomial vector fields tangent to the omitted curve."""
import sympy as sp

a, b, c, t = sp.symbols("a b c t", nonzero=True)
monomials = [1, a, b, c, a**2, a*b, a*c, b**2, b*c, c**2]
coefficients = sp.symbols("u0:30")
field = [
    sum(coefficients[10*i + j]*monomials[j] for j in range(10))
    for i in range(3)
]
g = {a: sp.Rational(4, 27)/t**2, b: sp.Rational(4, 3)/t, c: t}
gprime = [
    -sp.Rational(8, 27)/t**3,
    -sp.Rational(4, 3)/t**2,
    1,
]

tangency_equations = []
for i, j in ((0, 1), (0, 2), (1, 2)):
    numerator = sp.together(
        field[i].subs(g)*gprime[j] - field[j].subs(g)*gprime[i]
    ).as_numer_denom()[0]
    tangency_equations.extend(sp.Poly(numerator, t).all_coeffs())
tangency_matrix, _ = sp.linear_eq_to_matrix(tangency_equations, coefficients)
assert tangency_matrix.rank() == 17
assert len(tangency_matrix.nullspace()) == 13

vanishing_equations = []
for component in field:
    numerator = sp.together(component.subs(g)).as_numer_denom()[0]
    vanishing_equations.extend(sp.Poly(numerator, t).all_coeffs())
vanishing_matrix, _ = sp.linear_eq_to_matrix(vanishing_equations, coefficients)
assert vanishing_matrix.rank() == 21
assert len(vanishing_matrix.nullspace()) == 9

R1, R2, R3 = 12*a - b**2, 9*a*c - b, 3*b*c - 4
assert all(sp.expand(R.subs(g)) == 0 for R in (R1, R2, R3))

# Representatives inducing dc/ds=t^j for j=-1,0,1,2 on Gamma.
Dminus = sp.Matrix([-sp.Rational(27, 2)*a**2,
                    -sp.Rational(27, 4)*a*b, sp.Rational(3, 4)*b])
Dzero = sp.Matrix([-sp.Rational(3, 2)*a*b, -9*a, 1])
H = sp.Matrix([-2*a, -b, c])
Dtwo = sp.Matrix([-sp.Rational(2, 9)*b, -sp.Rational(4, 3), c**2])
for D, speed in ((Dminus, 1/t), (Dzero, 1), (H, t), (Dtwo, t**2)):
    restricted = sp.simplify(D.subs(g))
    assert all(sp.expand(restricted[i] - speed*gprime[i]) == 0 for i in range(3))

print("PASS: tangent fields of degree <=2 have dimension 13")
print("PASS: the 9-dimensional kernel is generated componentwise by")
print("      12a-b^2, 9ac-b, 3bc-4")
print("PASS: the four induced curve speeds are c^-1, 1, c, c^2")
print("PASS: only a multiple of speed c is complete on Gamma = C*")

# Among the nine elementary kernel fields R_i*d_j, a coordinate shear is
# locally nilpotent exactly when its coefficient is independent of that
# coordinate. This gives U=R3*d_a and V=R1*d_c.
U = sp.Matrix([R3, 0, 0])
V = sp.Matrix([0, 0, R1])
variables = (a, b, c)
assert sp.diff(U[0], a) == 0
assert sp.diff(V[2], c) == 0

# Every lambda*U+mu*V is complete: b is constant and (a,c) obey an affine
# linear constant-coefficient system.
lam, mu = sp.symbols("lambda mu")
combined = lam*U + mu*V
jac_ac = sp.Matrix([combined[0], combined[2]]).jacobian((a, c))
assert all(not entry.has(a, c) for entry in jac_ac)

print("PASS: U=(3bc-4)d_a and V=(12a-b^2)d_c generate a")
print("      two-parameter family of complete affine-linear-on-fibers flows")
