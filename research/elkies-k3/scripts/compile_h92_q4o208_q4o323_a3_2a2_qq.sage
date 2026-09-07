#!/usr/bin/env sage
"""Compile the corrected q4/o323 edge from its exact halved section.

This is a resolved-RR/chord-quartic calculation.  It tests the three old I4
supports for the single connected-component trace after the smooth P.O=1
saturation, then selects the unique trace whose binary-quartic Jacobian has
the prescribed A3+2A2 root type.  No elimination or Groebner basis is used.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
if not (ROOT / "MATH_STATUS.json").exists():
    ROOT = Path.cwd()
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
COMPACT = LOCAL / "q4o208-compact-weierstrass-qq.json"
HORIZONTAL = LOCAL / "q4o208-q4o323-horizontal-by-halving-qq.json"
MOD131 = LOCAL / "q4o208-physical-q4o323-horizontal-mod131.json"
MARKING = LOCAL / "q24-2a5-physical-q4o208-equation-marking-qq.json"
EDGE = GENERATED / "elkies-k3-h3-q4o208-physical-q4o323-corrected-a3-2a2-certificate.json"
SUFFIX = GENERATED / "elkies-k3-h3-q4o208-canonical-suffix-physical-nef-audit.json"
OUTPUT = LOCAL / "q4o208-q4o323-a3-2a2-rr-qq.json"
INPUTS = (COMPACT, HORIZONTAL, MOD131, MARKING, EDGE, SUFFIX)
started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficients(poly):
    return [str(value) for value in poly.list()]


compact = json.loads(COMPACT.read_text())
horizontal = json.loads(HORIZONTAL.read_text())
mod131 = json.loads(MOD131.read_text())
marking = json.loads(MARKING.read_text())
edge = json.loads(EDGE.read_text())
suffix = json.loads(SUFFIX.read_text())
assert compact["status"] == "PASS_EXACT_QQ_Q4O208_COMPACT_WEIERSTRASS_NORMALIZATION"
assert horizontal["status"] == "PASS_EXACT_QQ_Q4O323_HORIZONTAL_BY_MW_HALVING"
assert mod131["status"] == "PASS_EXACT_Q4O323_POLYNOMIAL_SECTION_SUBGROUP_OBSTRUCTION"
assert marking["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING"
assert edge["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
assert edge["candidate_id"]["label"] == "q4o323-physical-wall-corrected"
assert edge["child"]["root_data"] == [7, 24, 36]
assert suffix["status"] == "PASS_EXACT_Q4O208_CANONICAL_SUFFIX_PHYSICAL_WALL_CORRECTION"

R = PolynomialRing(QQ, "t")
t = R.gen()
A = R([QQ(value) for value in compact["compact_model"]["A_coefficients_low_to_high"]])
B = R([QQ(value) for value in compact["compact_model"]["B_coefficients_low_to_high"]])
x_record = horizontal["exact_QQ_horizontal"]["x"]
y_record = horizontal["exact_QQ_horizontal"]["y"]
Nx = R([QQ(value) for value in x_record["numerator_coefficients_low_to_high"]])
Dx = R([QQ(value) for value in x_record["denominator_coefficients_low_to_high"]])
Ny = R([QQ(value) for value in y_record["numerator_coefficients_low_to_high"]])
Dy = R([QQ(value) for value in y_record["denominator_coefficients_low_to_high"]])
assert Dx.is_square()
Z = Dx.sqrt()
dy_factors = list(Dy.factor())
assert len(dy_factors) == 1 and dy_factors[0][1] == 3
assert dy_factors[0][0] == Z
X, Y = Nx, Ny
assert (Z.degree(), X.degree(), Y.degree()) == (1, 6, 8)
assert Y**2 == X**3 + A*X*Z**4 + B*Z**6

target_class = vector(ZZ, mod131["target"]["NS_coordinates"])
physical_fibre = vector(ZZ, suffix["wall_correction"]["physical_fibre"])
old_zero = vector(ZZ, [-1, 1] + [0]*17)
vertical_residual = physical_fibre-old_zero-target_class
assert vertical_residual == vector(ZZ, [1, 0, 0, 0, 1, 1] + [0]*13)
second_cycle = [
    vector(ZZ, value) for value in
    marking["physical_fibres"]["second_old_I6_I4"]["components_in_cycle_order"]
]
assert vertical_residual == second_cycle[2]+second_cycle[3]
assert {
    record["label"]: record["support"]
    for record in compact["compact_model"]["reducible_fibres"]
}["second_old_I6_I4"] == "0"

# Ambient f=(a0+a1*t+a2*t^2+a3*t^3+b*Z*m)/Z^2.  Smoothness at
# the unique O-intersection Z=0 is a*X-b*Y == 0 mod Z^2.
smooth_rows = matrix(QQ, [
    [(t**degree*X).mod(Z**2)[order] for degree in range(4)]
    + [-Y.mod(Z**2)[order]]
    for order in range(2)
])
assert smooth_rows.rank() == 2
smooth_kernel = smooth_rows.right_kernel_matrix()
assert smooth_kernel.nrows() == 3

P = PolynomialRing(QQ, ("p", "q", "r"))
p, q, r = P.gens()
parameters = (p, q, r)
L = PolynomialRing(P, "l")
l = L.gen()
ambient = [
    sum(parameters[row]*smooth_kernel[row, column] for row in range(3))
    for column in range(5)
]
a_general = sum(ambient[degree]*l**degree for degree in range(4))
b_general = ambient[4]


def lift_poly(poly):
    return L([P(value) for value in poly.list()])


A_l, X_l, Y_l, Z_l = map(lift_poly, (A, X, Y, Z))
raw_general = (
    a_general**4 - 6*X_l*a_general**2*b_general**2
    + 8*Y_l*a_general*b_general**3 - 3*X_l**2*b_general**4
    - 4*A_l*b_general**4*Z_l**4
)
general_after_pole, remainder = raw_general.quo_rem(Z_l**4)
assert not remainder and general_after_pole.degree() == 8


def squared_linear_traces(value):
    return [
        factor for factor, exponent in P(value).factor()
        if factor.degree() == 1 and exponent == 2
    ]


trace_candidates = []
for support, value in (
    ("t=0", general_after_pole[0]),
    ("t=1", general_after_pole(L(1))),
    ("t=infinity", general_after_pole[8]),
):
    for trace in squared_linear_traces(value):
        trace_candidates.append((support, trace))


def compile_trace(support, trace):
    trace_row = vector(QQ, [trace.monomial_coefficient(value) for value in parameters])
    trace_kernel = matrix(QQ, [list(trace_row)]).right_kernel_matrix()
    h0 = trace_kernel*smooth_kernel
    assert h0.nrows() == 2 and h0.rank() == 2
    S = PolynomialRing(QQ, "u")
    u = S.gen()
    T = PolynomialRing(S, "t")
    tt = T.gen()
    a_u = sum(S(h0[0, degree]+u*h0[1, degree])*tt**degree for degree in range(4))
    b_u = S(h0[0, 4]+u*h0[1, 4])
    A_u, X_u, Y_u, Z_u = map(T, (A, X, Y, Z))
    raw = (
        a_u**4 - 6*X_u*a_u**2*b_u**2 + 8*Y_u*a_u*b_u**3
        - 3*X_u**2*b_u**4 - 4*A_u*b_u**4*Z_u**4
    )
    after_pole, pole_remainder = raw.quo_rem(Z_u**4)
    assert not pole_remainder
    factors = list(after_pole.factor())
    odd = [factor for factor, exponent in factors if exponent % 2]
    quartics = [factor for factor in odd if factor.degree() == 4]
    if len(quartics) != 1:
        return {"support": support, "trace": trace, "h0": h0, "factors": factors}
    quartic = quartics[0]
    values = list(quartic.list())
    if len(values) != 5:
        return {"support": support, "trace": trace, "h0": h0, "factors": factors}
    e, d, c, b, a = values
    I = 12*a*e-3*b*d+c**2
    J = 72*a*c*e+9*b*c*d-27*a*d**2-27*b**2*e-2*c**3
    A_child = -27*I
    B_child = -27*J
    Delta_child = -16*(4*A_child**3+27*B_child**2)
    delta_factors = list(Delta_child.factor())
    finite_repeated = sorted(
        (int(factor.degree()), int(exponent))
        for factor, exponent in delta_factors if exponent > 1
    )
    infinity_exponent = 24-int(Delta_child.degree())
    root_rank = sum(degree*(exponent-1) for degree, exponent in finite_repeated)
    if infinity_exponent > 1:
        root_rank += infinity_exponent-1
    return {
        "support": support, "trace": trace, "h0": h0, "factors": factors,
        "quartic": quartic, "A": A_child, "B": B_child, "Delta": Delta_child,
        "delta_factors": delta_factors, "finite_repeated": finite_repeated,
        "infinity_exponent": infinity_exponent, "root_rank": root_rank,
    }


compiled = [compile_trace(support, trace) for support, trace in trace_candidates]
hits = [record for record in compiled if record.get("root_rank") == 7]
selected_hits = [record for record in hits if record["support"] == "t=0"]
if len(selected_hits) != 1:
    summary = [
        (record["support"], str(record["trace"]), record.get("finite_repeated"),
         record.get("infinity_exponent"), record.get("root_rank"))
        for record in compiled
    ]
    raise ArithmeticError(f"expected one second-I4 A3+2A2 RR plane, found {len(selected_hits)}: {summary}")
selected = selected_hits[0]
assert sorted(
    (int(factor.degree()), int(exponent))
    for factor, exponent in selected["factors"]
) == [(1, 2), (1, 2), (4, 1)]
h0 = selected["h0"]
quartic = selected["quartic"]
A_child = selected["A"]
B_child = selected["B"]
Delta_child = selected["Delta"]
assert (A_child.degree(), B_child.degree(), Delta_child.degree()) == (8, 12, 20)
assert selected["finite_repeated"] == [(1, 3), (1, 3)]
assert selected["infinity_exponent"] == 4
nodal_factors = [
    factor for factor, exponent in selected["delta_factors"] if exponent == 1
]
assert len(nodal_factors) == 1 and nodal_factors[0].degree() == 14
assert nodal_factors[0].is_squarefree()
assert all(
    A_child.gcd(factor) == B_child.gcd(factor) == 1
    for factor, unused in selected["delta_factors"]
)
assert sum(factor.degree()*exponent for factor, exponent in selected["delta_factors"]) + 4 == 24

payload = {
    "schema": "elkies-k3.h3-q4o208-q4o323-a3-2a2-rr-qq.v1",
    "status": "PASS_EXACT_QQ_Q4O323_A3_2A2_RR_AND_JACOBIAN",
    "inputs": [
        {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)} for path in INPUTS
    ],
    "divisor": {
        "physical_fibre": list(map(int, physical_fibre)),
        "horizontal_NS_class": list(map(int, target_class)),
        "vertical_residual": list(map(int, vertical_residual)),
        "vertical_components": [
            "second_old_I6_I4_component_2",
            "second_old_I6_I4_component_3",
        ],
        "identity": "D = old_zero + horizontal_section + vertical_residual",
    },
    "resolved_RR": {
        "ambient_dimension": int(5),
        "smooth_saturation_rank": int(2),
        "resolved_component_trace_rank": int(1),
        "h0": int(2),
        "selected_support": selected["support"],
        "selected_trace_coefficients_in_p_q_r": [
            str(selected["trace"].monomial_coefficient(value)) for value in parameters
        ],
        "kernel_basis": [[str(value) for value in row] for row in h0.rows()],
    },
    "quartic": {
        "coefficients_in_t_low_to_high": [coefficients(value) for value in quartic.list()],
        "coefficient_degrees_in_u": [int(value.degree()) for value in quartic.list()],
    },
    "child": {
        "minimal_A_coefficients_low_to_high": coefficients(A_child),
        "minimal_B_coefficients_low_to_high": coefficients(B_child),
        "degrees_A_B_Delta": [int(A_child.degree()), int(B_child.degree()), int(Delta_child.degree())],
        "finite_fibres": [
            {"factor_coefficients_low_to_high": coefficients(factor),
             "degree": int(factor.degree()), "kodaira": f"I{int(exponent)}"}
            for factor, exponent in selected["delta_factors"]
        ],
        "infinity_kodaira": f"I{int(selected['infinity_exponent'])}",
        "ADE": "A3+2A2",
        "root_rank": int(selected["root_rank"]),
        "MW_rank_if_rho19": int(10),
        "lattice_root_data": edge["child"]["root_data"],
        "euler_number": int(24),
    },
    "method": {
        "trace_candidates": len(trace_candidates),
        "large_Groebner_required": False,
        "runtime_seconds": time.monotonic()-started,
    },
    "proof_boundary": (
        "Exact marked q4/o323 horizontal, resolved 5->3->2 RR plane, quartic, "
        "binary-quartic Jacobian and A3+2A2 fibre profile. Child zero/components and "
        "the outgoing q4/o207 equation marking remain separate gates."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q4O323RRQQ|ambient=5|smooth=3|h0=2|support={}|traces={}|"
    "ADE=A3+2A2|status={}|runtime={:.3f}|output={}".format(
        selected["support"], len(trace_candidates), payload["status"],
        payload["method"]["runtime_seconds"], OUTPUT,
    ),
    flush=True,
)
