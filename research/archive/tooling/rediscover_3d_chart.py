#!/usr/bin/env python3
"""Rediscover the announced 3D map from a sparse rational-chart ansatz."""
import sympy as sp

x, y, z, t, r, c = sp.symbols("x y z t r c")
a20, a11, a30, b10, b01, b21 = sp.symbols("a20 a11 a30 b10 b01 b21")
R = 2*x - 3*x**2*y - x**3*z
Phi = (y+1/x, 2/x, R)
A = a20*t**2 + a11*r*t + a30*c*t**3
B = b10*t + b01*r + b21*c*t**2
Psi = sp.Matrix([A, B, c])

# Normalize the r coefficient of B. Demand the reciprocal chart Jacobian r/2;
# pole cancellation of both pullbacks supplies the remaining equations.
equations = [b01-1]
j_residual = sp.Poly(sp.expand(Psi.jacobian((t, r, c)).det()-r/2), t, r, c)
equations += j_residual.coeffs()
pullbacks = [sp.cancel(f.subs({t: Phi[0], r: Phi[1], c: Phi[2]}, simultaneous=True)) for f in (A, B)]
for f in pullbacks:
    num, den = sp.fraction(f)
    if den != 1:
        poly = sp.Poly(sp.expand(num), x, y, z)
        order = sp.Poly(den, x).degree()
        for monomial, coefficient in poly.terms():
            if monomial[0] < order:
                equations.append(coefficient)

# The sparse system retains scale/triangular symmetries. Fix B_t=4, the
# normalization selected by the smallest integral collision-bearing member.
equations.append(b10-4)
solutions = sp.solve(equations, (a20,a11,a30,b10,b01,b21), dict=True)
assert solutions
target = next(sol for sol in solutions if sol.get(b21) == -3)
A0, B0 = sp.expand(A.subs(target)), sp.expand(B.subs(target))
F = [sp.cancel(f.subs({t: Phi[0], r: Phi[1], c: Phi[2]}, simultaneous=True)) for f in (A0, B0, c)]
assert all(sp.denom(f) == 1 for f in F)
assert sp.factor(sp.Matrix(F).jacobian((x,y,z)).det()) == -2
print("chart Phi =", Phi)
print("J(Phi) =", sp.factor(sp.Matrix(Phi).jacobian((x,y,z)).det()))
print("rediscovered Psi =", (A0, B0, c))
print("J(Psi) =", sp.factor(sp.Matrix([A0,B0,c]).jacobian((t,r,c)).det()))
print("pullback degrees =", tuple(sp.Poly(f,x,y,z).total_degree() for f in F))
print("PASS: poles cancel and composite determinant is -2")

