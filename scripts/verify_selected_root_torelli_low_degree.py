#!/usr/bin/env python3
"""Exact degree-four/five tests for the selected-root Torelli audit."""

import sympy as sp


W, r = sp.symbols("W r")


# Degree four: every quadratic Hessian divisor has a reflection.  Its center
# agrees with the intrinsic boundary center only at the centered mu_2 seed.
H4 = -W**2 * (W - 1) * (W - r) / (1 - r)
K4 = sp.factor(sp.diff(H4, W, 2))
K4_poly = sp.Poly(K4, W)
midpoint = sp.factor(
    -K4_poly.all_coeffs()[1] / (2 * K4_poly.all_coeffs()[0])
)
assert sp.simplify(midpoint - (r + 1) / 4) == 0
reflection4 = sp.factor(2 * midpoint - W)
assert sp.rem(
    sp.Poly(K4.subs(W, reflection4), W),
    sp.Poly(K4, W),
) == 0
assert sp.solve(sp.Eq(2 * midpoint, 0), r) == [-1]

H4_centered = sp.factor(H4.subs(r, -1))
assert sp.expand(H4_centered - W**2 * (1 - W**2) / 2) == 0
assert sp.expand(sp.diff(H4_centered, W, 2) - (1 - 6 * W**2)) == 0
assert sp.factor(H4_centered.subs(W, -W) - H4_centered) == 0
assert {-1, 1} == set(sp.solve(H4_centered / W**2, W))


# Degree five, reflection locus: the Hessian divisor is invariant under
# W -> 4-W, but the intrinsic boundary center 0 is sent to 4.
H5_reflection = -W**2 * (W - 1) * (5 * W**2 - 45 * W + 106) / 66
K5_reflection = sp.factor(sp.diff(H5_reflection, W, 2))
assert sp.expand(
    K5_reflection + (W - 2) * (50 * W**2 - 200 * W + 53) / 33
) == 0
assert sp.factor(K5_reflection.subs(W, 4 - W) + K5_reflection) == 0
assert (4 - W).subs(W, 0) == 4
assert sp.factor(H5_reflection.subs(W, 1)) == 0
assert sp.factor(sp.diff(H5_reflection, W).subs(W, 1)) == -1
assert sp.discriminant(H5_reflection / W**2, W) != 0
assert sp.discriminant(K5_reflection, W) != 0


# Degree five, centered mu_3 locus: the boundary marks survive, but the group
# cycles the three primitive roots, so a chosen affine sheet kills it.
zeta = sp.symbols("zeta")
H5_centered = W**2 * (1 - W**3) / 3
K5_centered = sp.factor(sp.diff(H5_centered, W, 2))
assert sp.expand(
    K5_centered + sp.Rational(2, 3) * (10 * W**3 - 1)
) == 0
mu3_difference = sp.Poly(
    sp.expand(K5_centered.subs(W, zeta * W) - K5_centered),
    zeta,
)
assert mu3_difference.rem(sp.Poly(zeta**3 - 1, zeta)).is_zero
assert sp.expand(
    H5_centered / W**2 + (W - 1) * (W**2 + W + 1) / 3
) == 0


# Twice-primitive reconstruction and normalized rerooting are exact in both
# degrees.  This is the algebraic core of the all-degree argument.
for H, rho in (
    (H4_centered, sp.Integer(-1)),
    (H5_reflection, sp.Integer(1)),
    (H5_centered, sp.Integer(1)),
):
    K = sp.diff(H, W, 2)
    J = sp.integrate(sp.integrate(K, W), W)
    J = sp.expand(J - J.subs(W, 0) - W * sp.diff(J, W).subs(W, 0))
    assert sp.expand(J - H) == 0
    rerooted = sp.factor(-J.subs(W, rho * W) / (rho * sp.diff(J, W).subs(W, rho)))
    assert rerooted.subs(W, 0) == 0
    assert sp.diff(rerooted, W).subs(W, 0) == 0
    assert sp.factor(rerooted.subs(W, 1)) == 0
    assert sp.factor(sp.diff(rerooted, W).subs(W, 1)) == -1


print("PASS degree four: boundary center and affine sheet kill Hessian reflection")
print("PASS degree five: boundary center kills mu_2 and affine sheet kills mu_3")
print("PASS selected-root algebra: twice-primitive rerooting is normalized")
