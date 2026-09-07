#!/usr/bin/env python3
"""Generated SymPy input for the common-arithmetic-fibers explicit example."""
import sympy as sp

S, x, y, z = sp.symbols("S x y z")
factor_coefficients_ascending = [['-19', '0', '0', '1'], ['1', '1', '1']]
factors = [
    sum(sp.Rational(c) * S**i for i, c in enumerate(coefficients))
    for coefficients in factor_coefficients_ascending
]
P5 = sp.expand(sp.prod(factors))
G5 = sp.expand(P5 - P5.subs(S, 0))
t = 1 + x*y
q = t**2*z - 19*y**2*(3*t + 1)
F = (
    q*t,
    -5*q**5*t**2*x**3 - 4*q**4*t**2*x**2 + 38*q*t - 3*q*x + 19*y,
    3*q**5*x**5 + 2*q**4*x**4 + x**3*z + 19*x*(5 - 3*t),
)
normalized_target = tuple(map(sp.Rational, ['1', '0', '-2']))
integral_scaling = tuple(map(sp.Rational, ['1', '19', '19']))
integral_target = tuple(a*b for a, b in zip(integral_scaling, normalized_target))
inverse_polynomial = sp.expand(
    G5 - sp.Poly(G5, S).coeff_monomial(S)
    * (normalized_target[1]*S**2 + normalized_target[2]) / 2
)
assert P5 == S**5 + S**4 + S**3 - 19*S**2 - 19*S - 19
assert inverse_polynomial == P5
assert integral_target == tuple(map(sp.Rational, ['1', '0', '-38']))
assert sp.factor(sp.Matrix(F).jacobian((x, y, z)).det()) == sp.Rational("-722")
print("PASS: generated SymPy paper example")
