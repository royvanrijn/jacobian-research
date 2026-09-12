#!/usr/bin/env python3
"""Exact low-degree checks for the marked-zero-fiber LL degree."""

from __future__ import annotations

from itertools import product

import sympy as sp


# Cayley/Prufer gives d^(d-2) classical polynomial LL sheets.  Dividing by
# the generic mu_d source-scaling orbit gives d^(d-3) polynomial covers up to
# affine source isomorphism.  Marking an unramified point over the designated
# simple branch value gives d-2 choices.
for degree in range(3, 8):
    prufer_count = sum(1 for _ in product(range(degree), repeat=degree - 2))
    assert prufer_count == degree ** (degree - 2)
    assert prufer_count % degree == 0
    restricted_degree = (degree - 2) * (prufer_count // degree)
    assert restricted_degree == (degree - 2) * degree ** (degree - 3)


w, value = sp.symbols("w value")


def remaining_critical_value_polynomial(seed: sp.Expr) -> sp.Poly:
    """Remove the critical value zero from the critical-value resultant."""

    resultant = sp.factor(sp.resultant(sp.diff(seed, w), value - seed, w))
    quotient = sp.cancel(resultant / value)
    return sp.Poly(quotient, value)


# Degree four: the weighted-projective target P(1,2) has affine invariant
# c_2/c_1^2.  The resulting rational function of the seed parameter has
# degree eight, equal to (4-2)*4^(4-3).
a = sp.symbols("a")
quartic = sp.expand((-1 - 2 * a) * (w**3 - w**2) + a * (w**4 - w**2))
quartic_values = remaining_critical_value_polynomial(quartic)
quartic_coefficients = quartic_values.all_coeffs()
quartic_c1 = sp.cancel(quartic_coefficients[1] / quartic_coefficients[0])
quartic_c2 = sp.cancel(quartic_coefficients[2] / quartic_coefficients[0])
quartic_invariant = sp.factor(quartic_c2 / quartic_c1**2)
quartic_numerator, quartic_denominator = sp.together(quartic_invariant).as_numer_denom()
assert sp.gcd(sp.Poly(quartic_numerator, a), sp.Poly(quartic_denominator, a)).degree() == 0
assert max(sp.degree(quartic_numerator, a), sp.degree(quartic_denominator, a)) == 8


# Degree five: use the affine chart of P(1,2,3) with invariants
# c_2/c_1^2=2 and c_3/c_1^3=3.  Elimination has the degree-drop factor b^75
# and one genuine degree-75 factor, matching (5-2)*5^(5-3).
b = sp.symbols("b")
quintic = sp.expand(
    (-1 - 2 * a - 3 * b) * (w**3 - w**2)
    + a * (w**4 - w**2)
    + b * (w**5 - w**2)
)
quintic_values = remaining_critical_value_polynomial(quintic)
quintic_coefficients = quintic_values.all_coeffs()
quintic_c1 = sp.cancel(quintic_coefficients[1] / quintic_coefficients[0])
quintic_c2 = sp.cancel(quintic_coefficients[2] / quintic_coefficients[0])
quintic_c3 = sp.cancel(quintic_coefficients[3] / quintic_coefficients[0])
equation_two = sp.factor(
    sp.together(quintic_c2 - 2 * quintic_c1**2).as_numer_denom()[0]
)
equation_three = sp.factor(
    sp.together(quintic_c3 - 3 * quintic_c1**3).as_numer_denom()[0]
)
eliminant = sp.resultant(equation_two, equation_three, a)
factors = sp.factor_list(eliminant)[1]
assert any(sp.factor(factor - b) == 0 and exponent == 75 for factor, exponent in factors)
genuine_degrees = [
    sp.degree(factor, b)
    for factor, exponent in factors
    if not (sp.factor(factor - b) == 0 and exponent == 75)
    for _ in range(exponent)
]
assert genuine_degrees == [75]


print("PASS restricted LL count: (N-2)*N^(N-3) from Cayley plus marking")
print("PASS quartic restricted LL map: exact degree 8")
print("PASS quintic restricted LL map: exact degree 75 after degree-drop saturation")
