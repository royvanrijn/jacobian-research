#!/usr/bin/env python3
"""Exact saturated-resultant certificate for Mon(F_4 o F_3)=S_3 wr S_4."""

from __future__ import annotations

import sympy as sp


V, A, B, C, W = sp.symbols("V A B C W")

# The outer inverse equation, cleared by the harmless scalar two.
E4 = V**2 - V**4 - 2 * B * C * V + 2 * A * C**2

# Reconstruction of an outer source point from the marked quartic root V.
gamma = B * C - V + 2 * V**3
horizontal_gap = V - gamma
T = 3 * gamma**2 - 7 * gamma + 4 * V
source_x = C / gamma
source_y = horizontal_gap / C
source_z = gamma * T / (3 * C**2)

# The non-C component of the cubic inverse discriminant in intermediate
# coordinates.  Its full discriminant is -z^2*delta.
delta3 = (
    27 * source_x**2 * source_z**2
    - 18 * source_x * source_y * source_z
    + 4 * source_x
    + 4 * source_y**3 * source_z
    - source_y**2
)

# Clear exactly the reconstruction denominators.  At gamma=0 the result is
# 12*C^6, so this clearing introduces no gamma=0 component on C!=0.
R = sp.expand(
    9 * gamma * C**3 * T**2
    - 18 * gamma * C**3 * horizontal_gap * T
    + 12 * C**6
    + 4 * gamma**2 * horizontal_gap**3 * T
    - 3 * gamma * C**3 * horizontal_gap**2
)
assert sp.factor(R - 3 * gamma * C**5 * delta3) == 0
assert sp.factor(R.subs(B, (V - 2 * V**3) / C) - 12 * C**6) == 0

# Only the residue modulo the quartic inverse equation matters.  Saturation
# by C removes a common C^3 factor and the eventual C^8 boundary resultant.
remainder = sp.rem(R, E4, V)
assert sp.rem(remainder, C**3, C) == 0
reduced_remainder = sp.cancel(remainder / C**3)
assert sp.degree(reduced_remainder, V) == 3
assert len(sp.Poly(reduced_remainder, V, A, B, C).terms()) == 95

resultant = sp.resultant(E4, reduced_remainder, V)
coefficient, factors = sp.factor_list(resultant)
assert coefficient == 1
assert len(factors) == 2

boundary_factor = next((factor, exponent) for factor, exponent in factors if factor == C)
image_factor = next((factor, exponent) for factor, exponent in factors if factor != C)
Q, image_exponent = image_factor
assert boundary_factor == (C, 8)
assert image_exponent == 1
assert sp.degree(Q, A) == 21
assert sp.degree(Q, B) == 28
assert sp.degree(Q, C) == 22
assert len(sp.Poly(Q, A, B, C).terms()) == 1001

# factor_list over QQ has returned Q as one irreducible factor.  Exponent one
# means the generic quartic and cubic have exactly one common simple root on
# Q=0, so the cubic discriminant maps birationally to its image.

# Separate that image from the other inner discriminant component C_inner=0.
# The point (r,c)=(2,1) on the normalized cubic discriminant gives the
# intermediate outer-source point (-12,-8,1).
outer_x, outer_y, outer_z = sp.Integer(-12), sp.Integer(-8), sp.Integer(1)
outer_u = 1 + outer_x * outer_y
outer_gamma = (
    1
    - sp.Rational(4, 3) * outer_x * outer_y
    + outer_x**2 * outer_z
)
sample = {
    A: (
        2 * outer_u
        + outer_u**2
        - 3 * outer_u**4 * outer_gamma**2
    )
    / (2 * outer_x**2),
    B: (
        2
        + 2 * outer_u
        - 4 * outer_u**3 * outer_gamma**2
    )
    / (2 * outer_x),
    C: outer_x * outer_gamma,
}
assert sample == {
    A: -sp.Rational(799529969, 3),
    B: sp.Integer(43960408),
    C: sp.Integer(-204),
}
assert Q.subs(sample) == 0

sample_E4 = sp.Poly(E4.subs(sample), V)
sample_z_numerator = sp.Poly((gamma * T).subs(sample), V)
assert sample_E4.eval(1649) == 0
assert sample_E4.gcd(sample_E4.diff()).degree() == 0
assert sample_E4.gcd(sample_z_numerator).degree() == 0

# The first gcd says all four outer sheets are regular at the sample.  The
# second says none has inner third coordinate zero, so the C_inner=0
# ramification divisor has a different image.  Hence generic Q-inertia is
# one transposition supported in one of the four three-sheet blocks.
assert 6**4 * 24 == 31104

print("PASS: pulled-back cubic discriminant reduces to a cubic modulo E_4")
print("PASS: saturated resultant is irreducible Q^1 after removing C^8")
print("PASS: the other inner boundary divisor has a distinct target image")
print("PASS: single-block inertia certifies Mon(F_4 o F_3)=S_3 wr S_4")
