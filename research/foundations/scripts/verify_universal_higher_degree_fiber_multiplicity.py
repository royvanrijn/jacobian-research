#!/usr/bin/env python3
"""Exact top-jet invariant for universal multiplicity in degree N>=6."""

from __future__ import annotations

import sympy as sp


N = sp.symbols("N", integer=True, positive=True)
s, second_coefficient = sp.symbols("s second_coefficient")

# For j>=4 the quadratic-gauge stable weights are
# w_j=(1-j,-j).  The top three satisfy one universal second-difference
# relation.
def weight(index: sp.Expr) -> sp.Matrix:
    return sp.Matrix((1 - index, -index))


assert weight(N - 2) + weight(N) - 2 * weight(N - 1) == sp.zeros(2, 1)

# A centered monic degree-N polynomial has top terms
# T^N+c_(N-2)T^(N-2)+..., with c_(N-2)=-Tr(eta^2)/2.
# Translation T=s+S gives these top three seed coefficients.
g_top = sp.Integer(1)
g_next = N * s
g_next_two = sp.binomial(N, 2) * s**2 + second_coefficient

stable_invariant = sp.factor(
    g_next_two * g_top / g_next**2
)
expected = sp.factor(
    sp.Rational(1, 2) * (N - 1) / N
    + second_coefficient / (N**2 * s**2)
)
assert sp.simplify(stable_invariant - expected) == 0

# Nonzero second trace moment means second_coefficient is nonzero, so the
# invariant has an exact pole of order two at s=0.
pole_numerator = sp.factor(
    s**2 * stable_invariant
)
assert pole_numerator.subs(s, 0) == second_coefficient / N**2

# Check the concrete weight relation and invariant in the first degrees.
for degree in range(6, 13):
    assert (
        weight(degree - 2)
        + weight(degree)
        - 2 * weight(degree - 1)
    ) == sp.zeros(2, 1)
    concrete = sp.factor(stable_invariant.subs(N, degree))
    assert sp.factor((s**2 * concrete).subs(s, 0)) == (
        second_coefficient / degree**2
    )

# Degree six is the first uniform case.  Expand every lower coefficient to
# verify that none enters the separating top-jet invariant and that all clean
# coefficient failures are only a finite translation set.
T = sp.symbols("T")
c4, c3, c2, c1, c0 = sp.symbols("c4 c3 c2 c1 c0")
sextic = T**6 + c4 * T**4 + c3 * T**3 + c2 * T**2 + c1 * T + c0
translated_sextic = sp.Poly(
    sp.expand(sextic.subs(T, s + T) - sextic.subs(T, s)),
    T,
)
sextic_jets = {
    index: translated_sextic.coeff_monomial(T**index)
    for index in range(1, 7)
}
assert sextic_jets == {
    1: 6 * s**5 + 4 * c4 * s**3 + 3 * c3 * s**2 + 2 * c2 * s + c1,
    2: 15 * s**4 + 6 * c4 * s**2 + 3 * c3 * s + c2,
    3: 20 * s**3 + 4 * c4 * s + c3,
    4: 15 * s**2 + c4,
    5: 6 * s,
    6: 1,
}
sextic_invariant = sp.factor(
    sextic_jets[4] * sextic_jets[6] / sextic_jets[5] ** 2
)
assert sp.simplify(
    sextic_invariant - (sp.Rational(5, 12) + c4 / (36 * s**2))
) == 0
assert sp.simplify(
    sp.diff(sextic_invariant, s) + c4 / (18 * s**3)
) == 0
assert all(
    sp.Poly(sextic_jets[index], s).as_expr() != 0
    for index in (1, 3, 4, 5, 6)
)

print("PASS: top-three stable weights satisfy (1,-2,1) in every N>=6")
print("PASS: J_N=(N-1)/(2N)+c_(N-2)/(N^2*s^2)")
print("PASS: nonzero second trace moment forces a nonconstant stable invariant")
print("PASS: every degree-six lower-coefficient and clean-locus path is explicit")
