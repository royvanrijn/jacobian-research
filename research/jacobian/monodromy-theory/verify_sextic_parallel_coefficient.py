#!/usr/bin/env python3
"""Replay the centered bracket and pole coefficients for the sextic exclusion."""

from math import gcd

from sympy import Poly, cancel, symbols


def admissible_degree_pairs(total_degree):
    pairs = {}
    for m in range(3, total_degree):
        for n in range(2, m):
            divisor = gcd(m, n)
            degree_bound = min(
                (total_degree - m)//(m//divisor),
                (total_degree - n)//(n//divisor),
            )
            if degree_bound >= 1:
                pairs[m, n] = degree_bound
    return pairs


assert admissible_degree_pairs(6) == {(3, 2): 1, (4, 2): 1}
assert admissible_degree_pairs(7) == {(3, 2): 1, (4, 2): 1}

Y, q, Dq, K, DK, f, Df, h, Dh = symbols("Y q Dq K DK f Df h Dh")
bracket = Poly((Df*Y + Dh*Y**3)*(2*q*Y) - (f + 3*h*Y**2)*(Dq*Y**2 + DK), Y)
assert bracket.nth(0) == -f*DK
assert bracket.nth(2) == 2*q*Df - f*Dq - 3*h*DK

a, x, q0, f0, B0, kappa = symbols("a x q0 f0 B0 kappa", nonzero=True)

for multiplicity, exponent in ((2, 1), (3, 1), (4, 1), (4, 2)):
    A = a*x**multiplicity
    q = q0*x**exponent
    f = f0*x**(1 + exponent - multiplicity)
    h = cancel(f*A*(2*q*f.diff(x) - f*q.diff(x))/(3*kappa))
    B_square = 4*q0*kappa/(exponent*a*f0)
    if exponent == 1 and multiplicity > 2:
        coefficient = f + 3*h*B_square/(4*q**2)
        assert cancel(coefficient - (4 - 2*multiplicity)*f0*x**(2 - multiplicity)) == 0
    else:
        coefficient = f*B0/(2*q) + h*B0*B_square/(8*q**3)
        divisor = 3 if exponent == 1 else 6
        expected = f0*B0*x**(1 - multiplicity)/(divisor*q0)
        assert cancel(coefficient - expected) == 0

print("Sextic parallel coefficient: degree pairs and exact pole coefficients verified.")
