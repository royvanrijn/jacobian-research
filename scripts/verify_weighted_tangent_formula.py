#!/usr/bin/env python3
"""Exact audit of the weighted tangent formula and its outer quotients.

This certificate deliberately starts from the two intrinsic identities

    B*C       = H'(W) + c*gamma,
    c*A*C**2  = W*(H'(W) + c*gamma) - H(W),

rather than attempting to parse the nested fractions in a rendered PDF.  It
proves the all-degree Jacobian factorization without expanding H, verifies the
coefficientwise primitive identity behind L_H, and performs a complete N=4
regression for H=W^4+W^3-2W^2.
"""

from __future__ import annotations

import sympy as sp


# ---------------------------------------------------------------------------
# Universal structural calculation: H is never expanded.
# ---------------------------------------------------------------------------

x, y, z, v, source_s = sp.symbols("x y z v source_s")
u, gamma, W, C = sp.symbols("u gamma W C", nonzero=True)
a, b, c = sp.symbols("a b c", nonzero=True)
A, B = sp.symbols("A B")
h, p, dp = sp.symbols("h p dp")  # H(W), H'(W), H''(W)

# First pass through invariant source coordinates.
J_source_invariants = sp.Matrix((x, x * y, x**2 * z)).jacobian((x, y, z)).det()
assert sp.factor(J_source_invariants - x**3) == 0

# Weighted source chart (x,v,S) -> (W,gamma,C).
gamma_source = 1 + a * v + b * source_s
W_source = (1 + v) * gamma_source
C_source = x * gamma_source
J_weighted_source = sp.Matrix((W_source, gamma_source, C_source)).jacobian(
    (x, v, source_s)
).det()
assert sp.factor(J_weighted_source - b * gamma_source**2) == 0

# Plane tangent map (W,gamma,C) -> (s,t,C), where
# s=H'(W)+c*gamma and t=W*s-H(W).
J_plane_tangent = sp.Matrix(
    (
        (dp, c, 0),
        (c * gamma + W * dp, c * W, 0),
        (0, 0, 1),
    )
).det()
assert sp.factor(J_plane_tangent + c**2 * gamma) == 0

# Target quotient map (A,B,C) -> (s,t,C)=(BC,cAC^2,C).
J_target_quotients = sp.Matrix((B * C, c * A * C**2, C)).jacobian((A, B, C)).det()
assert sp.factor(J_target_quotients + c * C**3) == 0

# Chain-rule ledger.  Since C=x*gamma, every nonconstant factor cancels.
J_total = sp.factor(
    (
        J_source_invariants
        * (b * gamma**2)
        * J_plane_tangent
        / J_target_quotients
    ).subs(C, x * gamma)
)
assert J_total == b * c


# ---------------------------------------------------------------------------
# The displayed nested fractions are exactly the solved intrinsic identities.
# ---------------------------------------------------------------------------

s_intrinsic = p + c * gamma
t_intrinsic = W * s_intrinsic - h
A_intrinsic = t_intrinsic / (c * C**2)
B_intrinsic = s_intrinsic / C

# c*L_H(W)=W*H'(W)-H(W), so after W=u*gamma and C=x*gamma:
L_at_W = (u * gamma * p - h) / c
A_displayed = (u + L_at_W / gamma**2) / x**2
B_displayed = (c + p / gamma) / x
assert sp.factor(
    A_intrinsic.subs({W: u * gamma, C: x * gamma}) - A_displayed
) == 0
assert sp.factor(B_intrinsic.subs(C, x * gamma) - B_displayed) == 0
assert sp.factor(h - B_intrinsic * C * W + c * A_intrinsic * C**2) == 0

# The primitive identity is coefficientwise and therefore valid in every
# degree: for a monomial h_k W^k, k>=2, integration gives
# L_k=(k-1)h_k W^k/c.
k = sp.symbols("k", integer=True, positive=True)
h_k = sp.symbols("h_k")
H_k = h_k * W**k
L_k = (k - 1) * h_k * W**k / c
assert sp.simplify(c * L_k - (W * sp.diff(H_k, W) - H_k)) == 0

# Polynomiality at the distinguished source boundary.  H(1)=0 and H'(1)=-c
# give L_H(1)=-1.  Also L_H'(1)=kappa=H''(1)/c.  The chosen value of a kills
# the remaining first-order v term, so the first numerator lies in (v^2,S),
# while the second lies in (v,S).  Under v=xy and S=x^2z these are divisible by
# x^2 and x respectively.
epsilon, kappa, quadratic_remainder = sp.symbols(
    "epsilon kappa quadratic_remainder"
)
a_endpoint = -(1 + kappa) / (2 + kappa)
u_endpoint = 1 + epsilon
gamma_endpoint = 1 + a_endpoint * epsilon
W_endpoint = sp.expand(u_endpoint * gamma_endpoint)
L_endpoint = (
    -1
    + kappa * (W_endpoint - 1)
    + quadratic_remainder * (W_endpoint - 1) ** 2
)
A_numerator_to_first_order = sp.series(
    u_endpoint + L_endpoint / gamma_endpoint**2,
    epsilon,
    0,
    2,
).removeO()
assert sp.factor(A_numerator_to_first_order) == 0
assert sp.factor(1 + kappa + a_endpoint * (2 + kappa)) == 0


# ---------------------------------------------------------------------------
# Complete quartic regression matching the reported N=4 test.
# ---------------------------------------------------------------------------

W4 = sp.symbols("W4")
H4 = W4**4 + W4**3 - 2 * W4**2
p4 = sp.diff(H4, W4)
c4 = sp.Integer(-3)
kappa4 = sp.cancel(sp.diff(H4, W4, 2).subs(W4, 1) / c4)
a4 = sp.cancel(-(1 + kappa4) / (2 + kappa4))
assert kappa4 == sp.Rational(-14, 3)
assert a4 == sp.Rational(-11, 8)

u4 = 1 + x * y
gamma4 = 1 + a4 * x * y + x**2 * z
marked_W4 = sp.expand(u4 * gamma4)
C4 = sp.expand(x * gamma4)
L4 = sp.integrate(
    W4 * sp.diff(H4, W4, 2) / c4,
    (W4, 0, W4),
)

A4 = sp.cancel((u4 + L4.subs(W4, marked_W4) / gamma4**2) / x**2)
B4 = sp.cancel((c4 + p4.subs(W4, marked_W4) / gamma4) / x)
assert all(sp.denom(component) == 1 for component in (A4, B4, C4))

J4 = sp.factor(sp.Matrix((A4, B4, C4)).jacobian((x, y, z)).det())
assert J4 == c4
assert sp.factor(B4 * C4 - (p4.subs(W4, marked_W4) + c4 * gamma4)) == 0
assert sp.factor(
    c4 * A4 * C4**2
    - (
        marked_W4 * (p4.subs(W4, marked_W4) + c4 * gamma4)
        - H4.subs(W4, marked_W4)
    )
) == 0
assert sp.factor(
    H4.subs(W4, marked_W4)
    - B4 * C4 * marked_W4
    + c4 * A4 * C4**2
) == 0

print("PASS: the intrinsic identities solve to the displayed outer quotients")
print("PASS: the all-degree chain-rule ledger gives det=b*c")
print("PASS: endpoint conditions clear x^2 and x coefficientwise")
print("PASS: N=4 expands polynomially with exact Jacobian -3")
print("PASS: H(W)-BCW+cAC^2 vanishes identically on the source")
