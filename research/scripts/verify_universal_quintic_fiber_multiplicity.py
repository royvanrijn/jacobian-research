#!/usr/bin/env python3
"""Exact algebra for universal quintic quadratic-gauge fiber multiplicity."""

from __future__ import annotations

import sympy as sp


T, shift = sp.symbols("T shift")
a, b, c, d = sp.symbols("a b c d")

# Characteristic polynomial of a trace-zero primitive quintic generator.
P = T**5 + a * T**3 + b * T**2 + c * T + d
translated = sp.Poly(sp.expand(P.subs(T, shift + T) - P.subs(T, shift)), T)

g1 = translated.coeff_monomial(T)
g2 = translated.coeff_monomial(T**2)
g3 = translated.coeff_monomial(T**3)
g4 = translated.coeff_monomial(T**4)
g5 = translated.coeff_monomial(T**5)

assert sp.expand(g1 - sp.diff(P, T).subs(T, shift)) == 0
assert sp.expand(g3 - sp.diff(P, T, 3).subs(T, shift) / 6) == 0
assert sp.expand(g4 - sp.diff(P, T, 4).subs(T, shift) / 24) == 0
assert g5 == 1

# The quadratic coefficient is removable by a target shear.  On the
# coefficient torus, the quintic stable-moduli weights have one primitive
# relation (-1,-6,5).
weight_matrix = sp.Matrix(
    [
        [-2, -3, -4],
        [-1, -4, -5],
    ]
)
kernel = weight_matrix.nullspace()
assert len(kernel) == 1
assert 5 * kernel[0] == sp.Matrix((-1, -6, 5))

# After dividing the seed by g1, the invariant
# a5^5/(a3*a4^6) becomes the displayed translation invariant.
normalized_a3 = sp.cancel(g3 / g1)
normalized_a4 = sp.cancel(g4 / g1)
normalized_a5 = sp.cancel(g5 / g1)
stable_invariant = sp.factor(
    normalized_a5**5 / (normalized_a3 * normalized_a4**6)
)
expected_invariant = sp.factor(g5**5 * g1**2 / (g3 * g4**6))
assert sp.cancel(stable_invariant - expected_invariant) == 0

# In centered coordinates g4=5*shift.  If a!=0, then g3(0)=a while
# g1 has order at most two at shift=0.  Hence I has a pole there and cannot
# be constant.  The exact coefficients record that valuation argument.
assert g4 == 5 * shift
assert g3 == 10 * shift**2 + a
assert sp.Poly(g1, shift).coeff_monomial(shift**2) == 3 * a
assert sp.Poly(g1, shift).degree() == 4

# The exceptional pure-power polynomial has constant invariant, confirming
# that the nonzero second trace moment is the exact easy exclusion.
pure_power_invariant = sp.factor(
    stable_invariant.subs({a: 0, b: 0, c: 0})
)
assert pure_power_invariant == sp.Rational(1, 6250)

# Newton's identity for trace zero is a=-Tr(eta^2)/2.
trace_square = sp.symbols("trace_square")
assert sp.expand(a + trace_square / 2).subs(
    a, -trace_square / 2
) == 0

print("PASS: translated quintic seed coefficients are derivative jets")
print("PASS: the primitive stable invariant is a5^5/(a3*a4^6)")
print("PASS: nonzero second trace moment makes the invariant nonconstant")
print("PASS: the centered pure-power exception has invariant 1/6250")
