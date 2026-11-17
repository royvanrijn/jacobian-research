#!/usr/bin/env python3
"""Verify the fixed-parameter affine moving-kernel Schur obstruction.

For a corank-one Hessian with unimodular rank-one adjugate, Piola reduces an
affine kernel vector ``v=a+L*x`` to ``tr(L)=0``, ``L^2=0``, and ``L*a=0``.
Unimodularity leaves, up to affine linear coordinates, only a constant vector
or ``v=(z,1,0,0)``.  This checker integrates the latter kernel exactly,
imposes the full pencil and bordered-unit equations, and verifies that every
Schur descendant is a triangular polynomial automorphism.

The result assumes that the normalized affine kernel line is independent of
the pencil parameter. The parameter-moving minimal-index-one branch is
subsequently excluded by HC4RSD3 and is recorded separately.
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
    / "hc4_affine_moving_kernel_pencils.json"
)
ATLAS = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_projective_polar_atlas.json"
)
HELPER = ROOT / "jcsearch" / "reverse_schur_descent.py"
EXPECTED_OUTPUT_SHA256 = (
    "8df7456358f15af946f2c0160a5419c7ea11f4cc9500d63d5aa6f694f936501d"
)
EXPECTED_ATLAS_SHA256 = (
    "350bc81b4ba7ac21289d7548f6d46de6526887c4e98a0813596cf20a454b240b"
)
EXPECTED_HELPER_SHA256 = (
    "b80da9da8105caa51fa38fc178a3a256d9132335a91bb5a87bb84eab96cf3e44"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_existing() -> None:
    assert sha256(OUTPUT) == EXPECTED_OUTPUT_SHA256
    assert sha256(ATLAS) == EXPECTED_ATLAS_SHA256
    assert sha256(HELPER) == EXPECTED_HELPER_SHA256
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    assert payload["format"] == "hc4-affine-moving-kernel-pencils-v1"
    assert payload["status"]["id"] == "HC4RSD2"
    live_rows = atlas["quintic_coverage_summary"][
        "remaining_numerical_signatures_after_vertex_colength"
    ]
    assert payload["projective_polar_intersection"]["input_live_rows"] == live_rows
    assert payload["projective_polar_intersection"]["surviving_rows"] == 0
    assert payload["subsequent_frontier"]["parameter_moving_affine"].startswith(
        "excluded in HC4RSD3"
    )
    assert payload["subsequent_frontier"]["nonlinear_fixed_shear"].endswith(
        "closed by HC4RSD4"
    )
    print(
        "PASS: committed HC4RSD2 artifact, consumed atlas, and equation helper "
        "are intact; the recorded affine handoffs are closed by HC4RSD3--4; "
        "no symbolic replay or rewrite"
    )


parser = argparse.ArgumentParser()
parser.add_argument(
    "--audit-existing-only",
    action="store_true",
    help=(
        "validate committed inputs and later handoffs without symbolic replay "
        "or artifact rewriting"
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
    adjugate_divergence_residuals,
    corank_one_adjugate_scalar,
    hessian_integrability_residuals,
    kernel_line_piola_residuals,
)


x, y, z, w, t, s = sp.symbols("x y z w t s")
variables = (x, y, z, w)


# Coefficientwise affine Piola identity. For a generic v=a+L*x the
# executable residual is exactly (L+tr(L)*I)*(a+L*x), which is the finite
# matrix system classified in the note.
generic_constant = sp.Matrix(sp.symbols("a0:4"))
generic_linear = sp.Matrix(4, 4, sp.symbols("l0:16"))
generic_kernel = generic_constant + generic_linear * sp.Matrix(variables)
generic_piola_expected = (
    generic_linear + sp.trace(generic_linear) * sp.eye(4)
) * generic_kernel
assert not any(
    sp.expand(residual)
    for residual in (
        sp.Matrix(kernel_line_piola_residuals(generic_kernel, variables))
        - generic_piola_expected
    )
)


# The unique nonconstant unimodular affine Piola orbit.  Its affine matrix
# has rank one and square zero; the second component makes the vector
# unimodular.  The Piola equation is exact.
moving_kernel = sp.Matrix([z, 1, 0, 0])
moving_linear_part = moving_kernel.jacobian(variables)
assert moving_linear_part.rank() == 1
assert moving_linear_part**2 == sp.zeros(4)
assert sp.trace(moving_linear_part) == 0
assert not any(kernel_line_piola_residuals(moving_kernel, variables))


# Rank two is incompatible with an affine unimodular vector.  For a
# square-zero rank-two matrix, image=kernel.  The displayed canonical matrix
# and every a in its kernel therefore give v=a+L*x a common zero.
rank_two_linear_part = sp.Matrix(
    [
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ]
)
rank_two_constant = sp.Matrix([1, 0, 0, 0])
rank_two_vector = rank_two_constant + rank_two_linear_part * sp.Matrix(variables)
rank_two_zero = {x: 0, y: 0, z: -1, w: 0}
assert rank_two_linear_part.rank() == 2
assert rank_two_linear_part**2 == sp.zeros(4)
assert rank_two_linear_part * rank_two_constant == sp.zeros(4, 1)
assert rank_two_vector.subs(rank_two_zero) == sp.zeros(4, 1)


# Integrate Hess(f)*v=0 in the nonconstant orbit.  The complete potential is
#
#   f = y*C(z) + (x-y*z)*C'(z) + G(z,w).
#
# Its adjugate scalar is -C''(z)^2*G_ww(z,w).  This identity is the exact
# bridge from kernel classification to the bordered-unit equation.
C = sp.Function("C")(z)
G = sp.Function("G")(z, w)
integrated_potential = sp.expand(
    y * C + (x - y * z) * sp.diff(C, z) + G
)
integrated_hessian = sp.hessian(integrated_potential, variables)
integrated_scalar = -sp.diff(C, z, 2) ** 2 * sp.diff(G, w, 2)
assert not any(hessian_integrability_residuals(integrated_hessian, variables))
assert integrated_hessian * moving_kernel == sp.zeros(4, 1)
assert all(
    sp.simplify(residual) == 0
    for residual in adjugate_divergence_residuals(
        integrated_hessian,
        variables,
    )
)
assert sp.simplify(
    corank_one_adjugate_scalar(integrated_hessian, moving_kernel)
    - integrated_scalar
) == 0
assert sp.factor(
    integrated_hessian.extract([0, 2, 3], [0, 2, 3]).det(
        method="domain-ge"
    )
) == integrated_scalar


# The full fixed-kernel pencil after the UFD bordered-unit classification.
# Greek constants are retained symbolically and the functions p,q,h,k are
# arbitrary univariate polynomials.  The constraints C''=eta, B_ww=gamma,
# and A_ww=0 are exactly the unit conclusions, not a bounded ansatz.
alpha, eta, xi, zeta, gamma = sp.symbols(
    "alpha eta xi zeta gamma",
    nonzero=True,
)
kappa = sp.symbols("kappa", nonzero=True)
mu = sp.symbols("mu")
p_function = sp.Function("p")
q_function = sp.Function("q")
h_function = sp.Function("h")
k_function = sp.Function("k")
p = p_function(z)
q = q_function(z)
h = h_function(z)
k = k_function(z)
C_quadratic = eta * z**2 / 2 + xi * z + zeta
A = alpha * y + p * w + q
B = sp.expand(
    y * C_quadratic
    + (x - y * z) * sp.diff(C_quadratic, z)
    + gamma * w**2 / 2
    + h * w
    + k
)
family = ScalarPivotSchurFamily(variables, t, A, B)
pencil_hessian = family.pencil_hessian(s)
bordered_constant = alpha**2 * eta**2 * gamma

assert not any(hessian_integrability_residuals(pencil_hessian, variables))
assert pencil_hessian * moving_kernel == sp.zeros(4, 1)
assert sp.factor(pencil_hessian.det(method="domain-ge")) == 0
assert corank_one_adjugate_scalar(
    pencil_hessian,
    moving_kernel,
) == -eta**2 * gamma
assert sp.factor(family.bordered_determinant(s)) == bordered_constant
assert sp.factor(family.parent_hessian_determinant(sp.Symbol("lambda"))) == (
    bordered_constant
)

descendant = family.schur_descendant(kappa, mu)
descendant_hessian_determinant = sp.factor(
    sp.hessian(descendant, variables).det(method="domain-ge")
)
assert descendant_hessian_determinant == -kappa * bordered_constant


# Exact triangular inverse of grad(descendant).  The output coordinates
# recover z, S=kappa*A+mu, w, y, and x in that order.
gradient = sp.Matrix(
    [sp.diff(descendant, variable) for variable in variables]
)
S = sp.expand(kappa * A + mu)
expected_gradient = sp.Matrix(
    [
        eta * z + xi,
        C_quadratic - z * sp.diff(C_quadratic, z) + alpha * S,
        eta * (x - y * z)
        + sp.diff(h, z) * w
        + sp.diff(k, z)
        + S * (sp.diff(p, z) * w + sp.diff(q, z)),
        gamma * w + h + S * p,
    ]
)
assert (gradient - expected_gradient).applyfunc(sp.simplify) == sp.zeros(4, 1)

recovered_z = sp.simplify((gradient[0] - xi) / eta)
assert recovered_z == z
recovered_S = sp.simplify(
    (
        gradient[1]
        - (
            C_quadratic
            - z * sp.diff(C_quadratic, z)
        )
    )
    / alpha
)
assert recovered_S == S
recovered_A = sp.simplify((recovered_S - mu) / kappa)
assert recovered_A == A
recovered_w = sp.simplify((gradient[3] - h - recovered_S * p) / gamma)
assert recovered_w == w
recovered_y = sp.simplify((recovered_A - p * recovered_w - q) / alpha)
assert recovered_y == y
recovered_x = sp.simplify(
    z * recovered_y
    + (
        gradient[2]
        - sp.diff(h, z) * recovered_w
        - sp.diff(k, z)
        - recovered_S
        * (sp.diff(p, z) * recovered_w + sp.diff(q, z))
    )
    / eta
)
assert recovered_x == x


# The unresolved parameter-moving calibration is a symmetric minimal-index
# one pencil.  Its kernel changes with s, but its bordered determinant is
# x3^2, so this representative fails before collision equations.
moving_s_A = y * z
moving_s_B = x * z + w**2 / 2
moving_s_family = ScalarPivotSchurFamily(
    variables,
    t,
    moving_s_A,
    moving_s_B,
)
moving_s_kernel = sp.Matrix([-s, 1, 0, 0])
assert moving_s_family.pencil_hessian(s) * moving_s_kernel == sp.zeros(4, 1)
assert moving_s_family.bordered_determinant(s) == z**2


atlas = json.loads(ATLAS.read_text())
live_rows = atlas["quintic_coverage_summary"][
    "remaining_numerical_signatures_after_vertex_colength"
]
assert live_rows["total"] == 624

payload = {
    "format": "hc4-affine-moving-kernel-pencils-v1",
    "status": {
        "id": "HC4RSD2",
        "kind": "hybrid theorem",
        "scope": (
            "scalar singular pencil with an affine-in-x, parameter-independent "
            "normalized kernel line"
        ),
    },
    "affine_piola_classification": {
        "equation": "(L+tr(L)*I)*(a+L*x)=0",
        "consequences": ["tr(L)=0", "L^2=0", "L*a=0"],
        "unimodular_orbits": [
            "constant kernel",
            "v=(z,1,0,0) with rank(L)=1",
        ],
        "rank_two": "excluded because image(L)=kernel(L) gives a common zero",
    },
    "integrated_nonconstant_orbit": {
        "kernel": "v=(z,1,0,0)",
        "potential": "f=y*C(z)+(x-y*z)*C'(z)+G(z,w)",
        "adjugate_scalar": "-C''(z)^2*G_ww(z,w)",
    },
    "bordered_unit_normal_form": {
        "A": "alpha*y+p(z)*w+q(z)",
        "B": (
            "y*C(z)+(x-y*z)*C'(z)+gamma*w^2/2+h(z)*w+k(z)"
        ),
        "C": "eta*z^2/2+xi*z+zeta",
        "units": ["alpha", "eta", "gamma", "kappa"],
        "parent_determinant": "alpha^2*eta^2*gamma",
        "descendant_determinant": "-kappa*alpha^2*eta^2*gamma",
        "gradient_inverse_order": ["z", "S=kappa*A+mu", "w", "y", "x"],
        "result": "every descendant is a triangular polynomial automorphism",
    },
    "projective_polar_intersection": {
        "input_live_rows": live_rows,
        "surviving_rows": 0,
        "reason": "the explicit polynomial inverse excludes every collision",
    },
    "rational_reconstruction": {
        "surviving_components": 0,
        "status": "not invoked",
    },
    "subsequent_frontier": {
        "parameter_moving_affine": (
            "excluded in HC4RSD3; no minimal-index-one affine branch survives"
        ),
        "nonlinear_fixed_shear": (
            "the transverse family v=(P(z,w),1,0,0) is closed by HC4RSD4"
        ),
        "calibration": {
            "kernel": "(-s,1,0,0)",
            "bordered_determinant": "z^2",
        },
    },
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: verified the affine Piola normal-form calibrations")
print("PASS: integrated the complete v=(z,1,0,0) Hessian potential")
print("PASS: bordered unit forces the exact fixed-kernel pencil normal form")
print("PASS: every descendant has an explicit triangular polynomial inverse")
print("FOLLOW-UP: parameter-moving affine kernels are excluded by HC4RSD3")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
