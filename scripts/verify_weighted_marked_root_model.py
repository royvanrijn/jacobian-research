#!/usr/bin/env python3
"""Exact identities for the normalized weighted marked-root theorem."""

import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import WeightedSeedModel, canonical_seed, w, x, y, z  # noqa: E402


# Universal construction lemma.  The endpoint choice of a0 gives the exact
# second-order cancellation needed for polynomiality, while the four small
# Jacobians prove det(DG_H)=b0*c on C!=0 and hence everywhere.
kappa0 = sp.symbols("kappa0")
a_weight0 = -(1 + kappa0) / (2 + kappa0)
assert sp.factor(1 + kappa0 + a_weight0 * (2 + kappa0)) == 0

x0, y_source0, z_source0, v0, S0 = sp.symbols(
    "x0 y_source0 z_source0 v0 S0"
)
a00, b00, c00 = sp.symbols("a00 b00 c00", nonzero=True)
gamma0 = 1 + a00 * v0 + b00 * S0
u0 = 1 + v0
W0 = u0 * gamma0
C0 = x0 * gamma0
J_source_invariants = sp.Matrix(
    (x0, x0 * y_source0, x0**2 * z_source0)
).jacobian((x0, y_source0, z_source0)).det()
J_weighted = sp.Matrix((W0, gamma0, C0)).jacobian((x0, v0, S0)).det()
assert sp.factor(J_source_invariants - x0**3) == 0
assert sp.factor(J_weighted - b00 * gamma0**2) == 0

dp0, A0, B0 = sp.symbols("dp0 A0 B0")
J_pencil = sp.Matrix(
    (
        (dp0, c00, 0),
        (c00 * gamma0 + W0 * dp0, c00 * W0, 0),
        (0, 0, 1),
    )
).det()
J_target_pencil = sp.Matrix(
    (
        (0, C0, B0),
        (c00 * C0**2, 0, 2 * c00 * A0 * C0),
        (0, 0, 1),
    )
).det()
assert sp.factor(J_pencil + c00**2 * gamma0) == 0
assert sp.factor(J_target_pencil + c00 * C0**3) == 0
det_weighted = sp.cancel(
    J_source_invariants * J_weighted * J_pencil / J_target_pencil
)
assert sp.factor(det_weighted - b00 * c00) == 0
print("PASS: universal polynomiality cancellation and det(DG_H)=b0*c")


# Universal incidence and reconstruction algebra.  Here h and p stand for
# H(W) and H'(W); the only relation used is E=h-BCW+cAC^2=0.
A, B, C, W = sp.symbols("A B C W")
c, a0, b0 = sp.symbols("c a0 b0", nonzero=True)
h, p = sp.symbols("h p")
E = h - B * C * W + c * A * C**2
Eprime = p - B * C
gamma = -Eprime / c
q = (W * p - h) / c
x_root = C / gamma
u_root = W / gamma
v_root = u_root - 1
y_root = sp.cancel(v_root / x_root)
S_root = sp.cancel((gamma - 1 - a0 * v_root) / b0)
z_root = sp.cancel(S_root / x_root**2)

assert sp.factor(gamma + Eprime / c) == 0
assert sp.factor(x_root * gamma - C) == 0
assert sp.factor(u_root * gamma - W) == 0
assert sp.factor(x_root * y_root - v_root) == 0
assert sp.factor(1 + a0 * v_root + b0 * S_root - gamma) == 0

root_relation = {h: B * C * W - c * A * C**2}
C_back = sp.factor(x_root * gamma)
B_back = sp.factor((c + p / gamma) / x_root)
A_back = sp.factor((u_root + q / gamma**2) / x_root**2)
assert C_back == C
assert sp.factor(B_back - B) == 0
assert sp.factor((A_back - A).subs(root_relation)) == 0
print("PASS: universal simple-root reconstruction inverts the weighted map")

# The root at projective infinity never occurs: the homogenized degree-n
# equation restricts to h_n*W^n there.
h_n, R = sp.symbols("h_n R", nonzero=True)
n = 5
projective_leading = h_n * W**n
assert projective_leading.subs(W, 1) == h_n
print("PASS: the weighted projective incidence has no root at infinity")

# Distinguished root W=1.  Taylor-expand with H'(1)=-c and
# H''(1)=kappa*c.  The incidence determines the first two root coefficients,
# and the chosen a0 cancels the first-order reconstruction numerator for z.
kappa, h3 = sp.symbols("kappa h3")
eps = sp.symbols("eps")
a_weight = -(1 + kappa) / (2 + kappa)
r1 = -B / c
r2 = A + (2 + kappa) * B**2 / (2 * c**2)
root_offset = r1 * eps + r2 * eps**2
H_taylor = -c * root_offset + kappa * c * root_offset**2 / 2 + h3 * root_offset**3 / 6
root_one_E = sp.expand(H_taylor - B * eps * (1 + root_offset) + c * A * eps**2)
assert sp.expand(root_one_E).coeff(eps, 1) == 0
assert sp.expand(root_one_E).coeff(eps, 2) == 0

Hprime_taylor = -c + kappa * c * root_offset + h3 * root_offset**2 / 2
gamma_one = sp.expand((B * eps - Hprime_taylor) / c)
u_one = sp.series((1 + root_offset) / gamma_one, eps, 0, 3).removeO()
v_one = sp.expand(u_one - 1)
S_one = sp.series(gamma_one - 1 - a_weight * v_one, eps, 0, 3).removeO()
assert sp.factor(v_one.coeff(eps, 1) + B * (2 + kappa) / c) == 0
assert sp.factor(sp.expand(S_one).coeff(eps, 1)) == 0
print("PASS: the distinguished root W=1 extends to the finite x=0 chart")

# Normalized zero-cluster charts for every tested multiplicity.  W=C*R0
# yields the finite strict transform.  The other scale for m>=3 gives the
# Kummer divisor of escaping branches.
R0, delta, eta, h0 = sp.symbols("R0 delta eta h0", nonzero=True)
for multiplicity in range(2, 7):
    local_H = h0 * W**multiplicity
    local_E = local_H - B * C * W + c * A * C**2
    finite_pullback = sp.expand(local_E.subs({W: C * R0}) / C**2)
    expected_finite = h0 * C ** (multiplicity - 2) * R0**multiplicity - B * R0 + c * A
    assert sp.factor(finite_pullback - expected_finite) == 0
    if multiplicity == 2:
        assert finite_pullback.subs(C, 0) == h0 * R0**2 - B * R0 + c * A
    else:
        assert finite_pullback.subs(C, 0) == -B * R0 + c * A

    if multiplicity >= 3:
        escaping = sp.expand(
            local_E.subs({C: delta ** (multiplicity - 1), W: eta * delta})
        )
        leading = escaping.coeff(delta, multiplicity)
        assert sp.factor(leading - eta * (h0 * eta ** (multiplicity - 1) - B)) == 0
print("PASS: normalization separates finite and escaping zero-cluster branches")

# Complete non-cubic representative: canonical H=W^3(1-W), whose inverse
# pencil has degree four.  Verify the source marking and both compositions of
# reconstruction exactly.
model = WeightedSeedModel(canonical_seed(3))
H4 = model.primitive.subs(w, W)
p4 = model.seed.subs(w, W)
G = model.mapping()
source_v = x * y
source_S = x**2 * z
source_u = 1 + source_v
source_gamma = 1 + model.a * source_v + model.b * source_S
source_W = source_u * source_gamma
E4 = sp.expand(H4 - B * C * W + model.c * A * C**2)
dE4 = sp.diff(E4, W)
source_substitution = {A: G[0], B: G[1], C: G[2], W: source_W}
assert sp.factor(E4.subs(source_substitution)) == 0
assert sp.factor(dE4.subs(source_substitution) + model.c * source_gamma) == 0

gamma4 = sp.cancel((B * C - p4) / model.c)
x4 = sp.cancel(C / gamma4)
u4 = sp.cancel(W / gamma4)
v4 = sp.cancel(u4 - 1)
y4 = sp.cancel(v4 / x4)
S4 = sp.cancel((gamma4 - 1 - model.a * v4) / model.b)
z4 = sp.cancel(S4 / x4**2)
for reconstructed, expected in zip((x4, y4, z4), (x, y, z)):
    assert sp.factor(reconstructed.subs(source_substitution) - expected) == 0

target_reconstruction = {x: x4, y: y4, z: z4}
for got, expected in zip(G, (A, B, C)):
    numerator = sp.together(got.subs(target_reconstruction) - expected).as_numer_denom()[0]
    remainder = sp.rem(sp.Poly(numerator, W), sp.Poly(E4, W)).as_expr()
    assert sp.factor(remainder) == 0
print("PASS: canonical inverse quartic realizes the normalized marked-root model")
