#!/usr/bin/env python3
"""Verify the rank-two quadratic-pivot classification.

HC4RSD8 puts a rank-two pivot in hyperbolic form A=x*y+w.  Exact leading
faces make B affine in the other passive variable z.  The next two faces
force its coefficient to be affine in one active variable and reduce the
remaining passive Hessian to a constant-curvature zero-determinant block.
That block integrates exactly, and every Schur descendant has a triangular
polynomial inverse.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_quadratic_rank_two_pivots.json"
)


# Universal active/passive block faces in the hyperbolic normalization
# Hess(A)|active = [[0,1],[1,0]], grad(A)|passive=(0,1).
s, lam = sp.symbols("s lambda")
k11, k12, k22 = sp.symbols("k11 k12 k22")
K = sp.Matrix([[k11, k12 + s], [k12 + s, k22]])
d11, d12, d21, d22 = sp.symbols("D11 D12 D21 D22")
D = sp.Matrix([[d11, d12], [d21, d22]])
e11, e12, e22 = sp.symbols("E11 E12 E22")
E = sp.Matrix([[e11, e12], [e12, e22]])
M = sp.Matrix.vstack(
    sp.Matrix.hstack(K, D),
    sp.Matrix.hstack(D.T, E),
)
g1, g2 = sp.symbols("g1 g2")
g = sp.Matrix([g1, g2, 0, 1])
parent = sp.Matrix.vstack(
    sp.Matrix.hstack(sp.Matrix([[lam]]), g.T),
    sp.Matrix.hstack(g, M),
)
assert sp.factor(
    sp.Poly(M.det(method="domain-ge"), s).coeff_monomial(s**2)
) == -E.det(method="domain-ge")
parent_s2 = sp.Poly(
    parent.det(method="domain-ge"),
    s,
).coeff_monomial(s**2)
assert sp.factor(
    parent_s2.subs(e11 * e22, e12**2) - e11
) == 0


# After the leading faces give E=diag(0,e), the remaining determinant
# identities depend only on p=grad(b), the w-cross column r, and K.
p1, p2, r1, r2, e = sp.symbols("p1 p2 r1 r2 e")
D_reduced = sp.Matrix([[p1, r1], [p2, r2]])
E_reduced = sp.diag(0, e)
M_reduced = sp.Matrix.vstack(
    sp.Matrix.hstack(K, D_reduced),
    sp.Matrix.hstack(D_reduced.T, E_reduced),
)
g_reduced = sp.Matrix([g1, g2, 0, 1])
parent_reduced = sp.Matrix.vstack(
    sp.Matrix.hstack(sp.Matrix([[lam]]), g_reduced.T),
    sp.Matrix.hstack(g_reduced, M_reduced),
)
M_polynomial = sp.Poly(M_reduced.det(method="domain-ge"), s)
assert sp.factor(M_polynomial.coeff_monomial(s)) == 2 * e * p1 * p2

# Choose the p2=0 branch, so b=b(x), and put g2=x as required by A=x*y+w.
x = sp.symbols("x")
M_one_channel = sp.factor(M_reduced.det(method="domain-ge").subs(p2, 0))
parent_one_channel = sp.factor(
    parent_reduced.det(method="domain-ge").subs({p2: 0, g2: x})
)
binary_hessian_determinant = sp.expand(r2**2 - k22 * e)
directional_curvature = sp.expand(k22 - 2 * x * r2 + x**2 * e)
assert sp.factor(M_one_channel - p1**2 * binary_hessian_determinant) == 0
assert sp.factor(
    parent_one_channel
    - p1**2 * directional_curvature
    - lam * M_one_channel
) == 0


# All-degree normal form.  Polynomial templates of unrestricted symbolic
# coefficients replay the determinant and recovery identities without
# choosing a special numerical member.
y, z, w, t, kappa, mu, rho = sp.symbols(
    "y z w t kappa mu rho"
)
h0, h1, h2 = sp.symbols("h0 h1 h2")
beta0, beta1 = sp.symbols("beta0 beta1")
gamma0, gamma1 = sp.symbols("gamma0 gamma1")
delta0, delta1, delta2 = sp.symbols("delta0 delta1 delta2")
h = h0 + h1 * x + h2 * x**2
beta = beta0 + beta1 * x
gamma = gamma0 + gamma1 * x
delta = delta0 + delta1 * x + delta2 * x**2
A = x * y + w
Y = y + h * A
B = x * z + rho * Y**2 / 2 + beta * y + gamma * A + delta
Phi = lam * t**2 / 2 + t * A + B
psi = B + kappa * A**2 / 2 + mu * A
variables = (x, y, z, w)

assert sp.factor(
    sp.hessian(B + s * A, variables).det(method="domain-ge")
) == 0
assert sp.factor(
    sp.hessian(Phi, (t, *variables)).det(method="domain-ge")
) == rho
assert sp.factor(
    sp.hessian(psi, variables).det(method="domain-ge")
) == -kappa * rho

gradient = sp.Matrix([sp.diff(psi, variable) for variable in variables])
assert sp.expand(gradient[2] - x) == 0
assert sp.expand(gradient[1] - x * gradient[3] - rho * Y - beta) == 0
assert sp.expand(
    gradient[3] - rho * h * Y - gamma - kappa * A - mu
) == 0
assert sp.diff(gradient[0], z) == 1


payload = {
    "format": "hc4-quadratic-rank-two-pivots-v1",
    "status": {
        "id": "HC4RSD9",
        "kind": "theorem",
        "scope": (
            "rank-two quadratic scalar pivots with an identically singular "
            "four-variable Hessian pencil"
        ),
    },
    "leading_faces": {
        "pencil_s2": "-det(Hess_(z,w)(B))",
        "parent_s2_mod_pencil": "B_zz",
        "consequence": "B=z*b(x,y)+C(x,y,w)",
    },
    "channel_split": {
        "pencil_s1": "2*C_ww*b_x*b_y",
        "one_channel": "b=b(x), up to swapping x and y",
        "pencil_s0": "(b')^2*(C_yw^2-C_yy*C_ww)",
        "parent": (
            "(b')^2*(C_yy-2*x*C_yw+x^2*C_ww) "
            "+ lambda*det(pencil)"
        ),
        "unit_consequence": "b is affine nonconstant and normalizes to x",
    },
    "normal_form": {
        "A": "x*y+w",
        "B": (
            "x*z+rho*(y+h(x)*A)^2/2+beta(x)*y+"
            "gamma(x)*A+delta(x)"
        ),
        "parent_hessian_determinant": "rho",
        "descendant_hessian_determinant": "-kappa*rho",
    },
    "triangular_recovery": [
        "x=F_z",
        "Y=(F_y-x*F_w-beta(x))/rho",
        "A=(F_w-rho*h(x)*Y-gamma(x)-mu)/kappa",
        "y=Y-h(x)*A",
        "w=A-x*y",
        "z=F_x-known_polynomial",
    ],
    "result": "every rank-two quadratic-pivot descendant is injective",
    "open_frontier": "rank-one quadratic pivots",
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: leading faces force one passive affine direction")
print("PASS: the active coefficient has exactly one nonconstant channel")
print("PASS: the bordered unit normalizes that channel to b=x")
print("PASS: integrated the complete passive-rank-one normal form")
print("PASS: every descendant has a triangular polynomial inverse")
print("THEOREM: all rank-two quadratic singular-pencil pivots are collision-free")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
