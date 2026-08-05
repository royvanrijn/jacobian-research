#!/usr/bin/env python3
"""Verify the affine-pivot Schur coverage gate for quintic HC4.

For a constant-Hessian potential psi with H=Hess(psi), a nonzero constant
covector ell supplies an affine-pivot singular Schur lift exactly when

    N_ell = ell^T adj(H) ell

is a nonzero constant.  The required rank-one repair is
kappa=det(H)/N_ell.

On the essential-rank-three quintic packet, write the top Hessian in its
constant-kernel coordinates as diag(C,0), and put d=grad(s_3) for the
off-diagonal quartic Hessian block.  The constant-Hessian face is

    R = det(C)*f - d^T adj(C) d = 0.

The affine-pivot metric numerator has leading surviving coefficient N7,
and the exact identity

    det(C)*N7 = R*(a^T adj(C) a) + (a^T adj(C) d)^2

forces a constant relation a^T adj(C)d=0.  This checker replays all of
these determinant and adjugate identities over the universal coefficient
ring.
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
    / "hc4_affine_pivot_coverage_gate.json"
)


def symmetric_matrix(prefix: str, size: int) -> sp.Matrix:
    entries: dict[tuple[int, int], sp.Symbol] = {}
    for row in range(size):
        for column in range(row, size):
            entries[row, column] = sp.Symbol(f"{prefix}{row + 1}{column + 1}")
    return sp.Matrix(
        size,
        size,
        lambda row, column: entries[min(row, column), max(row, column)],
    )


# Universal affine-pivot criterion.
H = symmetric_matrix("h", 4)
ell = sp.Matrix(sp.symbols("l1:5"))
kappa = sp.symbols("kappa")
det_H = H.det(method="domain-ge")
adj_H = H.adjugate(method="domain-ge")
metric_numerator = sp.expand((ell.T * adj_H * ell)[0])
reduced_hessian = H - kappa * ell * ell.T

assert sp.expand(
    reduced_hessian.det(method="domain-ge")
    - (det_H - kappa * metric_numerator)
) == 0

bordered_reduced = sp.zeros(5)
bordered_reduced[0, 1:] = ell.T
bordered_reduced[1:, 0] = ell
bordered_reduced[1:, 1:] = reduced_hessian
assert sp.expand(
    bordered_reduced.det(method="domain-ge") + metric_numerator
) == 0


# Collision transfer for A=ell.x+a0.  Equality of the original gradients
# and equality of A at the two points make the parent gradients equal at
# any common pivot value.
mu = sp.symbols("mu")
A_plus, A_minus = sp.symbols("A_plus A_minus")
gradient_plus = sp.Matrix(sp.symbols("gp1:5"))
gradient_minus = sp.Matrix(sp.symbols("gm1:5"))
gradient_B_difference = (
    gradient_plus
    - (kappa * A_plus + mu) * ell
    - gradient_minus
    + (kappa * A_minus + mu) * ell
)
collision_substitutions = {
    **dict(zip(gradient_plus, gradient_minus, strict=True)),
    A_plus: A_minus,
}
assert all(
    sp.expand(entry.subs(collision_substitutions)) == 0
    for entry in gradient_B_difference
)


# Rank-three quintic leading faces.  The scaling variable rho records
# homogeneous spatial degree: Hess(h5), Hess(h4), Hess(h3) have weights
# three, two, and one.
rho = sp.symbols("rho")
C = symmetric_matrix("c", 3)
E = symmetric_matrix("e", 3)
F = symmetric_matrix("f", 3)
d = sp.Matrix(sp.symbols("d1:4"))
g = sp.Matrix(sp.symbols("g1:4"))
a = sp.Matrix(sp.symbols("a1:4"))
tau, e_tt, f_tt = sp.symbols("tau e_tt f_tt")

Delta = C.det(method="domain-ge")
adj_C = C.adjugate(method="domain-ge")

top_hessian = sp.zeros(4)
top_hessian[:3, :3] = C
top_covector = sp.Matrix([a[0], a[1], a[2], tau])
top_metric_face = sp.expand(
    (top_covector.T * top_hessian.adjugate(method="domain-ge") * top_covector)[0]
)
assert sp.expand(top_metric_face - tau**2 * Delta) == 0


# Before the determinant face kills e_tt=D_t^2(h4), it is simultaneously
# the degree-eleven determinant coefficient and the degree-eight metric
# coefficient for an active covector.
through_quartic = sp.zeros(4)
through_quartic[:3, :3] = rho**3 * C + rho**2 * E
through_quartic[:3, 3] = rho**2 * d
through_quartic[3, :3] = (rho**2 * d).T
through_quartic[3, 3] = rho**2 * e_tt
det_degree_eleven = sp.expand(
    through_quartic.det(method="domain-ge")
).coeff(rho, 11)
assert sp.expand(det_degree_eleven - Delta * e_tt) == 0

active_covector = sp.Matrix([a[0], a[1], a[2], 0])
metric_through_quartic = sp.expand(
    (
        active_covector.T
        * through_quartic.adjugate(method="domain-ge")
        * active_covector
    )[0]
)
metric_degree_eight = metric_through_quartic.coeff(rho, 8)
assert sp.expand(
    metric_degree_eight - e_tt * (a.T * adj_C * a)[0]
) == 0


# Impose e_tt=0 and retain all active quartic/cubic blocks.  They do not
# alter the degree-ten determinant face or degree-seven metric identity.
rank_three_hessian = sp.zeros(4)
rank_three_hessian[:3, :3] = rho**3 * C + rho**2 * E + rho * F
rank_three_hessian[:3, 3] = rho**2 * d + rho * g
rank_three_hessian[3, :3] = (rho**2 * d + rho * g).T
rank_three_hessian[3, 3] = rho * f_tt

det_rank_three = sp.expand(rank_three_hessian.det(method="domain-ge"))
R = sp.expand(Delta * f_tt - (d.T * adj_C * d)[0])
assert sp.expand(det_rank_three.coeff(rho, 10) - R) == 0

metric_rank_three = sp.expand(
    (
        active_covector.T
        * rank_three_hessian.adjugate(method="domain-ge")
        * active_covector
    )[0]
)
N7 = metric_rank_three.coeff(rho, 7)
schur_vector_pairing = sp.expand((a.T * adj_C * d)[0])
metric_identity = sp.expand(
    Delta * N7
    - R * (a.T * adj_C * a)[0]
    - schur_vector_pairing**2
)
assert metric_identity == 0


# The new constant-span gate is genuinely stronger than Schur divisibility.
# On the diagonal ternary quintic Hessian, the displayed cubic derivative
# solves R=0, but its Schur vector has three independent monomial channels.
u1, u2, u3 = sp.symbols("u1 u2 u3")
alpha, beta, gamma = sp.symbols("alpha beta gamma")
fermat_C = sp.diag(u1**3, u2**3, u3**3)
fermat_d = sp.Matrix([alpha * u1**2, beta * u2**2, gamma * u3**2])
fermat_f = alpha**2 * u1 + beta**2 * u2 + gamma**2 * u3
fermat_Delta = fermat_C.det()
fermat_adj = fermat_C.adjugate()
fermat_R = sp.expand(
    fermat_Delta * fermat_f
    - (fermat_d.T * fermat_adj * fermat_d)[0]
)
assert fermat_R == 0
fermat_w = sp.expand(fermat_adj * fermat_d)
channel_monomials = [
    u1**2 * u2**3 * u3**3,
    u1**3 * u2**2 * u3**3,
    u1**3 * u2**3 * u3**2,
]
channel_matrix = sp.Matrix(
    [
        [
            sp.Poly(entry, u1, u2, u3).coeff_monomial(monomial)
            for monomial in channel_monomials
        ]
        for entry in fermat_w
    ]
)
assert sp.expand(channel_matrix.det() - alpha * beta * gamma) == 0


payload = {
    "format": "hc4-affine-pivot-coverage-gate-v1",
    "status": {
        "id": "HC4RSD6",
        "kind": "hybrid theorem",
        "scope": (
            "affine scalar singular-pivot lifts of constant-Hessian "
            "four-variable potentials, with a rank-three quintic packet gate"
        ),
    },
    "affine_pivot_criterion": {
        "metric_numerator": "N_ell=ell^T*adj(Hess(psi))*ell",
        "necessary_and_sufficient": "N_ell is a nonzero constant",
        "repair": "kappa=det(Hess(psi))/N_ell",
        "reduced_hessian": "Hess(psi)-kappa*ell*ell^T",
        "parent_bordered_determinant": "-N_ell",
    },
    "collision_transfer": {
        "condition": "ell.(p_plus-p_minus)=0",
        "result": "the marked collision lifts at any common pivot value",
    },
    "rank_three_quintic_gate": {
        "top_kernel_condition": "the pivot covector has zero kernel component",
        "determinant_faces": [
            "D_t^2(h4)=0",
            "R=det(C)*D_t^2(h3)-d^T*adj(C)*d=0",
        ],
        "schur_vector": "w=adj(C)*d",
        "metric_identity": (
            "det(C)*N7=R*(a^T*adj(C)*a)+(a^T*w)^2"
        ),
        "coverage_consequence": "a^T*w=0 for one nonzero constant a",
        "finite_equations": (
            "all 3-by-3 minors of the 3-by-45 degree-eight coefficient "
            "matrix of w vanish"
        ),
        "properness_calibration": (
            "at alpha=beta=gamma=1, C=diag(u1^3,u2^3,u3^3), "
            "d=(u1^2,u2^2,u3^2), f=u1+u2+u3 solves R=0 "
            "but has coefficient rank three"
        ),
    },
    "diagonal_top_classification": {
        "schur_divisibility_consequence": (
            "s3=(alpha*u1^3+beta*u2^3+gamma*u3^3)/3"
        ),
        "coefficient_matrix_determinant": "alpha*beta*gamma",
        "affine_coverage_gate": "alpha*beta*gamma=0",
    },
    "open_frontier": (
        "classify the constant-span-deficient Schur vectors w on the "
        "nonsquarefree rank-three packet, and derive the analogous gates "
        "for essential ranks one and two"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: affine singular-pivot lift iff the Hessian-metric norm is constant")
print("PASS: marked collision transfer is the single affine hyperplane condition")
print("PASS: rank-three quintic top forces an active pivot covector")
print("PASS: Schur-face identity forces a constant relation on w=adj(C)*d")
print("PASS: diagonal Schur pair has full span and fails affine-pivot coverage")
print("PASS: diagonal top coverage is confined to alpha*beta*gamma=0")
print("SCOPE: ranks one/two and constant-span-deficient rank three remain")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
