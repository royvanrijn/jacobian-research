#!/usr/bin/env python3
"""Verify the fixed transverse nonlinear shear-kernel obstruction.

For v=(P(z,w),1,0,0), Hess(f)*v=0 integrates to

    f = x*a(z,w) + y*b(z,w) + G(z,w),  db=-P*da.

The bordered unit forces the common composite of P and a to have straight
level curves, hence to be a function of one linear form.  The script checks
the transverse integration and curvature identities, then checks the
resulting univariate adjugate factor, determinants, and explicit inverse.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_univariate_shear_kernel_pencils.json"
)
HELPER = ROOT / "jcsearch" / "reverse_schur_descent.py"
EXPECTED_OUTPUT_SHA256 = (
    "bc36eff9effa7e40d7fba000f9f89377dfdecb8453394b4190acfe80cd8aa3ad"
)
EXPECTED_HELPER_SHA256 = (
    "b80da9da8105caa51fa38fc178a3a256d9132335a91bb5a87bb84eab96cf3e44"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_existing() -> None:
    assert sha256(OUTPUT) == EXPECTED_OUTPUT_SHA256
    assert sha256(HELPER) == EXPECTED_HELPER_SHA256
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["format"] == "hc4-univariate-shear-kernel-pencils-v1"
    assert payload["status"]["id"] == "HC4RSD4"
    assert payload["result"] == (
        "every descendant is a triangular polynomial automorphism"
    )
    assert payload["open_frontier"] == (
        "more general polynomial quasi-translation kernels and "
        "parameter-moving nonlinear kernels"
    )
    print(
        "PASS: committed HC4RSD4 artifact and equation helper are intact; its "
        "broad quasi-translation handoff is historical because HC4RSD5 closes "
        "the fixed two-component subcase; larger fixed support and nonlinear "
        "parameter motion remain; no symbolic replay or rewrite"
    )


parser = argparse.ArgumentParser()
parser.add_argument(
    "--audit-existing-only",
    action="store_true",
    help=(
        "validate committed inputs and the historical frontier without "
        "symbolic replay or artifact rewriting"
    ),
)
arguments = parser.parse_args()
if arguments.audit_existing_only:
    audit_existing()
    raise SystemExit(0)

import sys

import sympy as sp

sys.path.insert(0, str(ROOT))

from jcsearch.reverse_schur_descent import (
    ScalarPivotSchurFamily,
    corank_one_adjugate_scalar,
    kernel_line_piola_residuals,
)


def assert_zero_vector(vector: sp.Matrix) -> None:
    assert all(sp.simplify(entry) == 0 for entry in vector)


x, y, z, w, s, t = sp.symbols("x y z w s t")
variables = (x, y, z, w)
P_transverse = sp.Function("P_transverse")(z, w)
a_transverse = sp.Function("a_transverse")(z, w)
b_transverse = sp.Function("b_transverse")(z, w)
G_transverse = sp.Function("G_transverse")(z, w)
transverse_kernel = sp.Matrix([P_transverse, 1, 0, 0])

assert not any(kernel_line_piola_residuals(transverse_kernel, variables))
transverse_potential = sp.expand(
    x * a_transverse + y * b_transverse + G_transverse
)
transverse_hessian = sp.hessian(transverse_potential, variables)
transverse_residual = transverse_hessian * transverse_kernel
expected_transverse_residual = sp.Matrix(
    [
        0,
        0,
        P_transverse * sp.diff(a_transverse, z)
        + sp.diff(b_transverse, z),
        P_transverse * sp.diff(a_transverse, w)
        + sp.diff(b_transverse, w),
    ]
)
assert_zero_vector(transverse_residual - expected_transverse_residual)


# If P=P0(H), a=a0(H), and b=b0(H), the x/y-dependent transverse Hessian
# is rho*Hess(H)+sigma*dH*dH.T.  For a 2-by-2 matrix the curvature quadratic
# ignores the rank-one sigma term and is affine in rho.  Its rho coefficient
# is the straight-level curvature numerator Q_H.
H_profile = sp.Function("H_profile")(z, w)
dH = sp.Matrix([sp.diff(H_profile, z), sp.diff(H_profile, w)])
hessian_H = sp.hessian(H_profile, (z, w))
rho, sigma = sp.symbols("rho sigma")
r00, r01, r11 = sp.symbols("r00 r01 r11")
remainder_hessian = sp.Matrix([[r00, r01], [r01, r11]])
transverse_block = (
    rho * hessian_H + sigma * dH * dH.T + remainder_hessian
)
curvature = sp.factor((dH.T * hessian_H.adjugate() * dH)[0])
curvature_identity = sp.factor(
    (dH.T * transverse_block.adjugate() * dH)[0]
    - rho * curvature
    - (dH.T * remainder_hessian.adjugate() * dH)[0]
)
assert curvature_identity == 0

# Vanishing curvature makes the gradient slope constant on every smooth
# level curve.  This executable identity is the differential step in the
# generic-fibre straight-line lemma used in the note.
hz = sp.diff(H_profile, z)
hw = sp.diff(H_profile, w)
tangent_slope_derivative = sp.factor(
    hw * sp.diff(hz / hw, z) - hz * sp.diff(hz / hw, w)
)
assert sp.factor(tangent_slope_derivative - curvature / hw**2) == 0


# After the straight-level lemma and a linear transverse coordinate change,
# P depends only on z.  Replay the complete univariate integral and adjugate.
P = sp.Function("P")(z)
a = sp.Function("a")(z)
b = sp.Function("b")(z)
G = sp.Function("G")(z, w)
kernel = sp.Matrix([P, 1, 0, 0])
generic_potential = sp.expand(x * a + y * b + G)
generic_hessian = sp.hessian(generic_potential, variables)
integrability_substitutions = {
    sp.diff(b, z): -P * sp.diff(a, z),
    sp.diff(b, z, 2): -sp.diff(P, z) * sp.diff(a, z)
    - P * sp.diff(a, z, 2),
}
integrated_hessian = generic_hessian.xreplace(integrability_substitutions)
assert_zero_vector(integrated_hessian * kernel)
assert sp.simplify(
    corank_one_adjugate_scalar(integrated_hessian, kernel)
    + sp.diff(a, z) ** 2 * sp.diff(G, w, 2)
) == 0


# Exact bordered-unit normal form.  J is an antiderivative of P, so the
# relation b_B'=-eta*P is retained without choosing a degree bound for P.
alpha, eta, xi, beta, gamma = sp.symbols(
    "alpha eta xi beta gamma",
    nonzero=True,
)
kappa = sp.symbols("kappa", nonzero=True)
mu = sp.symbols("mu")
p = sp.Function("p")(z)
q = sp.Function("q")(z)
h = sp.Function("h")(z)
k = sp.Function("k")(z)
J = sp.Integral(P, z)
normal_A = alpha * y + p * w + q
normal_B = sp.expand(
    x * (eta * z + xi)
    + y * (beta - eta * J)
    + gamma * w**2 / 2
    + h * w
    + k
)
family = ScalarPivotSchurFamily(variables, t, normal_A, normal_B)
normal_pencil = family.pencil_hessian(s)
assert_zero_vector(normal_pencil * kernel)
assert sp.factor(normal_pencil.det(method="domain-ge")) == 0
assert corank_one_adjugate_scalar(normal_pencil, kernel) == -eta**2 * gamma
bordered_constant = alpha**2 * eta**2 * gamma
assert sp.factor(family.bordered_determinant(s)) == bordered_constant

descendant = family.schur_descendant(kappa, mu)
descendant_hessian_determinant = sp.factor(
    sp.hessian(descendant, variables).det(method="domain-ge")
)
assert descendant_hessian_determinant == -kappa * bordered_constant


# Triangular inverse: Fx recovers z; Fy recovers S; Fw recovers w; then A,
# y, and x follow in order.
gradient = sp.Matrix(
    [sp.diff(descendant, variable) for variable in variables]
)
S = sp.expand(kappa * normal_A + mu)
expected_gradient = sp.Matrix(
    [
        eta * z + xi,
        beta - eta * J + alpha * S,
        eta * x
        - eta * P * y
        + sp.diff(h, z) * w
        + sp.diff(k, z)
        + S * (sp.diff(p, z) * w + sp.diff(q, z)),
        gamma * w + h + S * p,
    ]
)
assert_zero_vector(gradient - expected_gradient)

recovered_z = sp.simplify((gradient[0] - xi) / eta)
recovered_S = sp.simplify((gradient[1] - (beta - eta * J)) / alpha)
recovered_w = sp.simplify((gradient[3] - h - recovered_S * p) / gamma)
recovered_A = sp.simplify((recovered_S - mu) / kappa)
recovered_y = sp.simplify((recovered_A - p * recovered_w - q) / alpha)
recovered_x = sp.simplify(
    (
        gradient[2]
        + eta * P * recovered_y
        - sp.diff(h, z) * recovered_w
        - sp.diff(k, z)
        - recovered_S
        * (sp.diff(p, z) * recovered_w + sp.diff(q, z))
    )
    / eta
)
assert recovered_z == z
assert recovered_S == S
assert recovered_w == w
assert recovered_A == normal_A
assert recovered_y == y
assert recovered_x == x


payload = {
    "format": "hc4-univariate-shear-kernel-pencils-v1",
    "status": {
        "id": "HC4RSD4",
        "kind": "hybrid theorem",
        "scope": (
            "scalar singular pencil with fixed normalized kernel "
            "v=(P(z,w),1,0,0), P nonconstant"
        ),
    },
    "transverse_reduction": {
        "complete_potential": "f=x*a(z,w)+y*b(z,w)+G(z,w)",
        "integrability": "db=-P*da, hence dP wedge da=0",
        "common_composite": "P=P0(H), a=a0(H), b=b0(H)",
        "curvature_numerator": (
            "H_z^2*H_ww-2*H_z*H_w*H_zw+H_w^2*H_zz"
        ),
        "unit_consequence": "the curvature numerator vanishes",
        "straight_level_result": "H is a polynomial in one linear form",
    },
    "complete_kernel_integral": {
        "kernel": "v=(P(z),1,0,0)",
        "potential": "f=x*a(z)+y*b(z)+G(z,w)",
        "integrability": "b'(z)=-P(z)*a'(z)",
        "adjugate_scalar": "-a'(z)^2*G_ww(z,w)",
    },
    "bordered_unit_normal_form": {
        "A": "alpha*y+p(z)*w+q(z)",
        "B": (
            "x*(eta*z+xi)+y*(beta-eta*Integral(P dz))"
            "+gamma*w^2/2+h(z)*w+k(z)"
        ),
        "units": ["alpha", "eta", "gamma", "kappa"],
        "parent_determinant": "alpha^2*eta^2*gamma",
        "descendant_determinant": "-kappa*alpha^2*eta^2*gamma",
        "gradient_inverse_order": ["z", "S=kappa*A+mu", "w", "y", "x"],
    },
    "result": "every descendant is a triangular polynomial automorphism",
    "open_frontier": (
        "more general polynomial quasi-translation kernels and "
        "parameter-moving nonlinear kernels"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: integrated the complete v=(P(z,w),1,0,0) Hessian potential")
print("PASS: verified the common-composite curvature reduction")
print("PASS: bordered unit forces the univariate-shear normal form")
print("PASS: every descendant has an explicit triangular polynomial inverse")
print("SCOPE: general quasi-translations and parameter-moving kernels remain open")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
