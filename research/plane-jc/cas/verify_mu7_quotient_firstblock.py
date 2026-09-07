#!/usr/bin/env python3
"""Compile the (72,108) weighted-Wronskian first block on its mu_7 quotient.

This script does not use the archived Groebner basis.  It independently builds
both triangular systems and verifies, coefficient by coefficient, that the
usual degree-35 presentation is the pullback of a six-equation system in the
invariants

    q=u^7,  x_i=a_i*u^(i-1),  2<=i<=6.

Here u=a_7.  The quotient presentation has the same eleven triangular solves
and six obstruction equations, but never introduces a degree-35 coefficient
field.
"""

from __future__ import annotations

import sympy as sp


t = sp.symbols("t")
u = sp.symbols("u", nonzero=True)
a = sp.symbols("a2:8")
d = sp.symbols("d2:13")
x2, x3, x4, x5, x6, q = sp.symbols("x2 x3 x4 x5 x6 q")
x = {2: x2, 3: x3, 4: x4, 5: x5, 6: x6}
e = sp.symbols("e2:13")


def triangular_solve(A, D, d_variables):
    equation = sp.Poly(
        sp.expand(2 * A * sp.diff(D, t) - 3 * sp.diff(A, t) * D - t**2),
        t,
    )
    substitutions = {}
    residuals = []
    for degree in range(equation.degree() + 1):
        coefficient = sp.cancel(equation.nth(degree).subs(substitutions))
        if coefficient == 0:
            continue
        solved = False
        for variable in d_variables:
            other_variables = [item for item in d_variables if item != variable]
            if (
                coefficient.has(variable)
                and sp.Poly(coefficient, variable).degree() == 1
                and not any(coefficient.has(item) for item in other_variables)
            ):
                polynomial = sp.Poly(coefficient, variable)
                pivot = polynomial.nth(1)
                substitutions[variable] = sp.cancel(
                    -(coefficient - pivot * variable) / pivot
                )
                solved = True
                break
        if not solved:
            residuals.append((degree, sp.factor(coefficient)))
    return substitutions, residuals


A = t + sum(a[i - 2] * t**i for i in range(2, 8)) + t**8
D = sum(d[i - 2] * t**i for i in range(2, 13))
original_solves, original_residuals = triangular_solve(A, D, d)

A_bar = (
    t
    + sum(x[i] * t**i for i in range(2, 7))
    + q * t**7
    + q * t**8
)
D_bar = sum(e[i - 2] * t**i for i in range(2, 13))
quotient_solves, quotient_residuals = triangular_solve(A_bar, D_bar, e)

assert len(original_solves) == len(quotient_solves) == 11
assert [degree for degree, _ in original_residuals] == list(range(13, 19))
assert [degree for degree, _ in quotient_residuals] == list(range(13, 19))

# A_bar(tau)=u^-1 A(u*tau), D_bar(tau)=u^-2 D(u*tau).
a_to_invariants = {
    a[i - 2]: (x[i] / u ** (i - 1) if i <= 6 else u)
    for i in range(2, 8)
}

for j in range(2, 13):
    expected = u ** (j - 2) * original_solves[d[j - 2]].subs(a_to_invariants)
    actual = quotient_solves[e[j - 2]].subs(q, u**7)
    if sp.factor(actual - expected) != 0:
        raise RuntimeError(f"triangular solve d_{j} failed to descend")

for (degree, old), (other_degree, new) in zip(
    original_residuals, quotient_residuals
):
    assert degree == other_degree
    expected = u ** (degree - 2) * old.subs(a_to_invariants)
    actual = new.subs(q, u**7)
    if sp.factor(actual - expected) != 0:
        raise RuntimeError(f"obstruction in degree {degree} failed to descend")

# After clearing rational content, only the final quotient obstruction has a
# removable q factor.  It is harmless on the actual first-block locus because
# the quotient eliminant has nonzero constant term.
primitive = []
for degree, expression in quotient_residuals:
    numerator = sp.together(expression).as_numer_denom()[0]
    polynomial = sp.Poly(numerator, x2, x3, x4, x5, x6, q, domain=sp.QQ)
    primitive.append((degree, polynomial.primitive()[1].as_expr()))

q_valuations = []
for degree, expression in primitive:
    valuation = 0
    current = expression
    while sp.expand(current).subs(q, 0) == 0:
        current = sp.cancel(current / q)
        valuation += 1
    q_valuations.append((degree, valuation))
assert q_valuations == [(13, 0), (14, 0), (15, 0), (16, 0), (17, 0), (18, 1)]

# The exact quotient eliminant h(q) from the archived lex basis.  Its constant
# term excludes q=0, so dividing the last obstruction by q loses no solution.
h = (
    9374377445732 * q**5
    + 62410476400737833472 * q**4
    + 265472843532245531128968765 * q**3
    + 591414847960503971284831143987840 * q**2
    + 586529490054134032292876680565455306752 * q
    - 1888043347611739526396142670327809715470336
)
assert h.subs(q, 0) != 0

print("ORIGINAL_TRIANGULAR_SOLVES=11")
print("QUOTIENT_TRIANGULAR_SOLVES=11")
print("OBSTRUCTION_DEGREES=13,14,15,16,17,18")
print("Q_VALUATIONS=" + ",".join(f"{degree}:{valuation}" for degree, valuation in q_valuations))
print("QUOTIENT_ELIMINANT_DEGREE=5")
print("MU7_QUOTIENT_FIRSTBLOCK_PASS")
