#!/usr/bin/env python3
"""Exact certificates for the commuting inverse-Jacobian frame."""
import sympy as sp

x, y, z = sp.symbols("x y z")
alpha, beta, gamma = sp.symbols("alpha beta gamma")
u = 1 + x*y
F = sp.Matrix([
    u**3*z + y**2*u*(4 + 3*x*y),
    y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y),
    2*x - 3*x**2*y - x**3*z,
])
J = F.jacobian((x, y, z))
X = sp.simplify(J.inv())

assert sp.factor(J.det()) == -2
assert sp.simplify(J*X) == sp.eye(3)

variables = (x, y, z)
for i in range(3):
    for j in range(i + 1, 3):
        bracket = X[:, j].jacobian(variables)*X[:, i] - X[:, i].jacobian(variables)*X[:, j]
        assert all(sp.expand(entry) == 0 for entry in bracket)

r = sp.symbols("r", nonzero=True)
curve_r = {x: 1/r, y: -sp.Rational(3, 2)*r, z: sp.Rational(13, 2)*r**2}
assert sp.simplify(F.subs(curve_r)) == sp.Matrix([-r**2/4, 0, 0])

c = sp.symbols("c", nonzero=True)
gp = sp.Matrix([-sp.Rational(8, 27)/c**3, -sp.Rational(4, 3)/c**2, 1])
v = sp.Matrix([alpha, beta, gamma])
minors = [sp.together(v[i]*gp[j] - v[j]*gp[i]) for i in range(3) for j in range(i + 1, 3)]
numerators = [sp.Poly(m.as_numer_denom()[0], c) for m in minors]
coefficients = [q for p in numerators for q in p.all_coeffs()]
assert sp.solve(coefficients, (alpha, beta, gamma), dict=True) == [
    {alpha: 0, beta: 0, gamma: 0}
]

print("PASS: the inverse-Jacobian columns are polynomial and commute")
print("PASS: F(1/r,-3r/2,13r^2/2)=(-r^2/4,0,0)")
print("PASS: Gamma has no nonzero constant translation direction")

# Exact local cubic at a point g(c) of Gamma. With T=t0+U and target
# g(c)+h*v, it has the form c*U^3+h*(A0+A1*U+A2*U^2+gamma*U^3).
h, U = sp.symbols("h U")
t0 = sp.Rational(2, 3)/c
a0, b0 = sp.Rational(4, 27)/c**2, sp.Rational(4, 3)/c
P_local = sp.expand(
    (c + h*gamma)*(t0 + U)**3
    - 2*(t0 + U)**2
    + (b0 + h*beta)*(t0 + U)
    - 2*(a0 + h*alpha)
)
A0 = sp.factor(gamma*t0**3 + beta*t0 - 2*alpha)
A1 = sp.factor(3*gamma*t0**2 + beta)
A2 = sp.factor(3*gamma*t0)
assert sp.expand(P_local - (c*U**3 + h*(A0 + A1*U + A2*U**2 + gamma*U**3))) == 0

# A0=A1 is precisely the tangent line to Gamma at g(c).
tangent_solution = sp.solve([A0, A1], (alpha, beta), dict=True)
assert tangent_solution == [{
    alpha: -sp.Rational(8, 27)*gamma/c**3,
    beta: -sp.Rational(4, 3)*gamma/c**2,
}]

print("PASS: exact Gamma strata are controlled by A0, A1, A2")
print("      A0!=0: exponent 2/3; A0=0,A1!=0: exponent 1;")
print("      A0=A1=0,v!=0: exponent 2")

# A nonconstant complete combination: lift the weighted Euler field
# H=(-2a,-b,c) on the target.
H_of_F = sp.Matrix([-2*F[0], -F[1], F[2]])
Y = sp.Matrix([x, -y, -2*z])
assert all(sp.expand(q) == 0 for q in J*Y - H_of_F)

a, b, C = sp.symbols("a b C")
Q = 27*a**2*C**2 - 18*a*b*C + 16*a + b**3*C - b**2
HQ = sp.expand(-2*a*sp.diff(Q, a) - b*sp.diff(Q, b) + C*sp.diff(Q, C))
assert HQ == -2*Q

# Classify all affine target vector fields tangent to Gamma. The solution is
# the one-dimensional span of H.
M_symbols = sp.symbols("m00:03 m10:13 m20:23")
M = sp.Matrix(3, 3, M_symbols)
d_symbols = sp.symbols("d0:3")
d = sp.Matrix(d_symbols)
g = sp.Matrix([sp.Rational(4, 27)/c**2, sp.Rational(4, 3)/c, c])
affine_field = M*g + d
tangency_equations = []
for i, j in ((0, 1), (0, 2), (1, 2)):
    numerator = sp.together(
        affine_field[i]*gp[j] - affine_field[j]*gp[i]
    ).as_numer_denom()[0]
    tangency_equations.extend(sp.Poly(numerator, c).all_coeffs())
unknowns = list(M) + list(d)
affine_solutions = sp.linsolve(tangency_equations, unknowns)
free_parameter = M[2, 2]
expected = (-2*free_parameter, 0, 0, 0, -free_parameter, 0,
            0, 0, free_parameter, 0, 0, 0)
assert len(affine_solutions.args) == 1
assert tuple(affine_solutions.args[0]) == expected

print("PASS: (-2F1)X1-F2X2+F3X3 = x*d_x-y*d_y-2z*d_z is complete")
print("PASS: its target field spans all affine vector fields preserving Gamma")

# The other elementary complete target shear V=(12a-b^2)*d_c also has an
# incomplete lift. This exact branch maps to (0,1,1-r^2).
escape_v = {
    x: -(r + 1)/r,
    y: 1,
    z: -r*(r - 3),
}
assert sp.simplify(F.subs(escape_v)) == sp.Matrix([0, 1, 1 - r**2])
# Here 12a-b^2=-1, so the lifted target field is -X_3. The branch starts at
# r=1 and escapes at r=0, corresponding to finite flow time.

print("PASS: the lift of (12a-b^2)*d_c is incomplete via")
print("      F(-(r+1)/r,1,-r(r-3))=(0,1,1-r^2)")
