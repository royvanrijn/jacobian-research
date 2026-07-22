#!/usr/bin/env python3
"""Exact audit of the proposed distinguished root-one normalization mark.

The root-one sheet is intrinsic in the regular-reconstruction open, but is
disjoint from the ramification divisor.  A rerooted quartic demonstrates that
bare marked-root pencils identify their normalization coordinates by a
scaling which sends the new affine root one to an old extra boundary root.
This regression supports the gap statement in
extended-geometry/DECORATED_NORMALIZATION_INVARIANT.md; it does not disprove a
future generator-rigidity theorem for the full cover with its affine open.
"""

from __future__ import annotations

import sympy as sp


w, s, t, A, B, C = sp.symbols("w s t A B C")

# A normalized boundary-clean split quartic with roots 0 (double), 1, and 3.
H = sp.expand(sp.Rational(1, 2) * w**2 * (w - 1) * (w - 3))
assert H.subs(w, 0) == sp.diff(H, w).subs(w, 0) == 0
assert H.subs(w, 1) == 0
assert sp.diff(H, w).subs(w, 1) == -1
assert sp.gcd(sp.Poly(H, w), sp.Poly(sp.diff(H, w), w)).monic().as_expr() == w

# The root-one affine sheet cannot meet the ramification divisor over C=0.
inverse = H - B * C * w + A * C**2
assert inverse.subs({C: 0, w: 1}) == 0
assert sp.diff(inverse, w).subs({C: 0, w: 1}) == -1

# Its competitor at zero is a degree-two affine component over k(A,B).
h2 = sp.expand(H).coeff(w, 2)
zero_chart = h2 * w**2 - B * w + A
zero_discriminant = sp.discriminant(zero_chart, w)
assert sp.factor(zero_discriminant - (B**2 - 4 * h2 * A)) == 0
assert sp.Poly(zero_discriminant, A).degree() == 1  # hence it is not a square in k(A,B)

# Reroot at the old extra root a=3.  The new seed again has root one and
# derivative -1 there, but that root corresponds to the old boundary root 3.
a = sp.Integer(3)
Hprime_a = sp.diff(H, w).subs(w, a)
kappa = sp.cancel(-1 / (a * Hprime_a))
G = sp.expand(kappa * H.subs(w, a * w))
assert G.subs(w, 1) == 0
assert sp.diff(G, w).subs(w, 1) == -1

E_G = sp.expand(G - s * w + t)
s_H = sp.cancel(s / (kappa * a))
t_H = sp.cancel(t / kappa)
E_H_transported = sp.expand(kappa * (H.subs(w, a * w) - s_H * a * w + t_H))
assert sp.expand(E_G - E_H_transported) == 0

# The transported G root one is the H root a, and it is an extra-root
# boundary branch for H by the exact reconstruction criterion.
assert H.subs(w, a) == 0
assert sp.diff(H, w).subs(w, a) != 0
assert sp.expand(a + sp.diff(H, w).subs(w, a)) != 0
assert a != 1

print("PASS affine root-one audit: unique simple affine stratum is unramified")
print("PASS affine root-one audit: rerooting scales the bare discriminant normalization")
print("SCOPE: cross-stratum generator rigidity for the distinguished affine open remains open")
