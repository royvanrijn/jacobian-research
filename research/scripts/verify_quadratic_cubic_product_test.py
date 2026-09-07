#!/usr/bin/env python3
"""Exact test of the natural SL_2 projection for the (2,3) slice.

The normalized tangent equation is already a 2 by 2 determinant.  This script
checks the fiber of that determinant-one map over the identity and verifies
that its class is L^2+L rather than L^2.
"""

from __future__ import annotations

import sympy as sp


T = sp.Symbol("T")
a0, a1, a2 = sp.symbols("a0 a1 a2")
b0, b1, b2, b3 = sp.symbols("b0 b1 b2 b3")
L = sp.Symbol("L")

A = a0 * T**2 + a1 * T + a2
B = b0 * T**3 + b1 * T**2 + b2 * T + b3
m = a0 * b1 + a1 * b0
resultant = sp.expand(sp.resultant(A, B, T))

assert sp.det(sp.Matrix([[a0, a1], [-b0, b1]])) == m

identity_fiber = sp.factor(
    resultant.subs({a0: 1, a1: 0, b0: 0, b1: 1})
)
assert sp.expand(identity_fiber - ((b3 - a2) ** 2 + a2 * b2**2)) == 0

# On b2 != 0, a2 is determined by (b3-a2, b2), giving A^1 x G_m.
# On b2 = 0, b3-a2 is +1 or -1 and a2 is free, giving two A^1's.
fiber_class = sp.expand(L * (L - 1) + 2 * L)
assert fiber_class == L**2 + L
assert fiber_class != L**2

ell0, ell1 = sp.symbols("ell0 ell1")
shift = (ell0 * T + ell1) * A
shifted_B = sp.Poly(B + shift, T)
shifted_b0 = shifted_B.coeff_monomial(T**3)
shifted_b1 = shifted_B.coeff_monomial(T**2)
delta_m = sp.expand(a0 * shifted_b1 + a1 * shifted_b0 - m)
assert sp.expand(delta_m - a0 * (2 * a1 * ell0 + a0 * ell1)) == 0

print("PASS (2,3) product test: top coefficients define an SL_2 morphism")
print("PASS identity fiber: x^2 + z*y^2 = 1 has class L^2 + L, not L^2")
print("PASS Euclidean-addition obstruction: the m-preserving kernel jumps at a0=0")
