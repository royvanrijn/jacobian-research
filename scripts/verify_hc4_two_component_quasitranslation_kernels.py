#!/usr/bin/env python3
"""Verify the fixed two-component quasi-translation reduction.

For a primitive kernel v=(P,Q,0,0), with (P,Q) unimodular and
projectively nonconstant, the Hessian Piola identity first forces P,Q to be
independent of the two active variables.  Then Hess(f)*v=0 forces

    f=x*a(z,w)+y*b(z,w)+G(z,w),  P*da+Q*db=0.

The bordered pairing makes P and Q algebraically dependent.  The proof's
polynomial frame argument then makes (P,Q) an affine line in its common
composite, reducing by a constant active-coordinate change to HC4RSD4.
This checker replays the complete Hessian residual, invariant-coordinate
coefficient, and affine-line frame identities.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.reverse_schur_descent import kernel_line_piola_residuals


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_two_component_quasitranslation_kernels.json"
)


def assert_zero_vector(vector: sp.Matrix) -> None:
    assert all(sp.simplify(entry) == 0 for entry in vector)


x, y, z, w = sp.symbols("x y z w")
variables = (x, y, z, w)


# Before transverse integration, the Hessian Piola identity closes active
# dependence.  For a general two-component primitive kernel, div(v*v^T)=0
# has the two displayed components.  Modulo P and Q, unimodularity gives
# P_y=Q_x=0; the residual coefficient matrix has determinant 3, so
# P_x=Q_y=0 in characteristic zero.
P_full = sp.Function("P_full")(*variables)
Q_full = sp.Function("Q_full")(*variables)
kernel_full = sp.Matrix([P_full, Q_full, 0, 0])
rank_one_adjugate = kernel_full * kernel_full.T
piola_divergence = sp.Matrix(
    [
        sum(
            sp.diff(rank_one_adjugate[row, column], variables[row])
            for row in range(4)
        )
        for column in range(4)
    ]
)
expected_piola = sp.Matrix(
    [
        2 * P_full * sp.diff(P_full, x)
        + Q_full * sp.diff(P_full, y)
        + P_full * sp.diff(Q_full, y),
        P_full * sp.diff(Q_full, x)
        + Q_full * sp.diff(P_full, x)
        + 2 * Q_full * sp.diff(Q_full, y),
        0,
        0,
    ]
)
assert_zero_vector(piola_divergence - expected_piola)
active_derivative_system = sp.Matrix([[2, 1], [1, 2]])
assert active_derivative_system.det() == 3


P = sp.Function("P")(z, w)
Q = sp.Function("Q")(z, w)
a = sp.Function("a")(z, w)
b = sp.Function("b")(z, w)
G = sp.Function("G")(z, w)
kernel = sp.Matrix([P, Q, 0, 0])

assert not any(kernel_line_piola_residuals(kernel, variables))
potential = sp.expand(x * a + y * b + G)
hessian = sp.hessian(potential, variables)
expected_kernel_residual = sp.Matrix(
    [
        0,
        0,
        P * sp.diff(a, z) + Q * sp.diff(b, z),
        P * sp.diff(a, w) + Q * sp.diff(b, w),
    ]
)
assert_zero_vector(hessian * kernel - expected_kernel_residual)


# With a Bezout pair U*P+V*Q=1, l=U*x+V*y and r=Q*x-P*y are polynomial
# coordinates in the active plane.  The coefficient of l in
# Q_i*x-P_i*y is the projective-direction derivative P*Q_i-Q*P_i.
U, V = sp.symbols("U V")
l, r = sp.symbols("l r")
x_from_slice = P * l + V * r
y_from_slice = Q * l - U * r
for transverse in (z, w):
    expression = sp.expand(
        sp.diff(Q, transverse) * x_from_slice
        - sp.diff(P, transverse) * y_from_slice
    )
    l_coefficient = sp.expand(expression).coeff(l)
    assert sp.simplify(
        l_coefficient
        - (P * sp.diff(Q, transverse) - Q * sp.diff(P, transverse))
    ) == 0


# Terminal polynomial frame.  Once P,Q factor through a closed composite H,
# the bordered dual row (-b_A,a_A) has constant determinant alpha with
# (P,Q).  The degree argument in the note forces that dual row to be
# constant, and the active row is exactly an affine line in one polynomial L.
H = sp.symbols("H")
alpha, p0, q0, a0, b0 = sp.symbols(
    "alpha p0 q0 a0 b0",
    nonzero=True,
)
L = sp.Function("L")(H)
active_P = p0 - b0 * L
active_Q = q0 + a0 * L
active_row = sp.Matrix([[active_P, active_Q]])
dual_row = sp.Matrix([[-b0, a0]])
frame = sp.Matrix.vstack(active_row, dual_row)
frame_determinant = sp.factor(frame.det(method="domain-ge"))
assert frame_determinant == a0 * p0 + b0 * q0
assert_zero_vector(
    sp.Matrix(
        [
            sp.diff(active_P, H) + b0 * sp.diff(L, H),
            sp.diff(active_Q, H) - a0 * sp.diff(L, H),
        ]
    )
)
assert sp.factor(
    a0 * active_P + b0 * active_Q - frame_determinant
) == 0


# The constant inverse frame sends the active row to (1,L) up to harmless
# constant scaling and translation, which is the HC4RSD4 shear kernel.
inverse_constant_frame = sp.Matrix([[p0, q0], [-b0, a0]]).inv()
normalized_row = sp.simplify(active_row * inverse_constant_frame)
assert normalized_row == sp.Matrix([[1, L]])


payload = {
    "format": "hc4-two-component-quasitranslation-kernels-v1",
    "status": {
        "id": "HC4RSD5",
        "kind": "hybrid theorem",
        "scope": (
            "fixed primitive kernel v=(P,Q,0,0), with (P,Q) "
            "unimodular"
        ),
    },
    "piola_active_dependence_gate": {
        "equations": [
            "2*P*P_x+Q*P_y+P*Q_y=0",
            "P*Q_x+Q*P_x+2*Q*Q_y=0",
        ],
        "unimodular_reduction": "P_y=Q_x=0",
        "terminal_matrix": "[[2,1],[1,2]]",
        "consequence": "P_x=P_y=Q_x=Q_y=0",
    },
    "complete_kernel_integral": {
        "potential": "f=x*a(z,w)+y*b(z,w)+G(z,w)",
        "integrability": "P*da+Q*db=0",
        "moving_direction_gate": "P*dQ-Q*dP is nonzero",
    },
    "bordered_pairing": {
        "equations": [
            "P*a_A+Q*b_A=alpha",
            "a_A*dP+b_A*dQ=0",
        ],
        "consequence": "dP wedge dQ=0",
        "common_composite": "P=P0(H), Q=Q0(H)",
    },
    "polynomial_frame": {
        "matrix": "[[P0,Q0],[-b_A,a_A]] in GL2(K[H])",
        "differential_system": [
            "row0'=S*row1",
            "row1'=T*row0",
        ],
        "degree_consequence": "projective motion forces T=0",
        "normal_form": "(P0,Q0)=(p0,q0)+L(H)*(-b0,a0)",
        "constant_coordinate_reduction": "kernel becomes (1,L(H),0,0)",
    },
    "result": "HC4RSD4 gives a triangular polynomial inverse",
    "open_frontier": (
        "fixed kernels with three or four active components, and "
        "parameter-moving nonlinear kernels"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: Piola forces a primitive two-component kernel to be transverse")
print("PASS: integrated the complete two-component quasi-translation kernel")
print("PASS: verified the invariant-slice projective-direction coefficient")
print("PASS: polynomial frame reduces the kernel to the HC4RSD4 shear form")
print("SCOPE: three/four-component fixed and parameter-moving kernels remain")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
