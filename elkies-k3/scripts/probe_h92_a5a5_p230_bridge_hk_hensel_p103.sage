#!/usr/bin/env sage -python
"""Lift the reduced q4/orbit230 bridge ansatz at p=103.

status: EXPERIMENT, promoted to exact only after literal QQ identities
claim: 23-variable joint Hensel lift of Q0 and R=Q2+Q3 using recovered poles

Let Q0=-P1229 be the selected integral section branch.  For R=Q2+Q3 write

    x(R)=XR/L^2,  y(R)=YR/L^3,
    XR=x(Q0)L^2+H Z,  YR=y(Q0)L^3+H K,

where L is the reconstructed linear pole of R and Z is the reconstructed
quadratic pole of P230=Q0+R.  We jointly lift the integral section Q0 and the
unknowns deg(H)<=4, deg(K)<=5.  The resulting 28 coefficient equations have a
23 by 23 unit Jacobian minor at p=103.  No Groebner basis is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, Qp, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
SURFACE = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"
POOL = LOCAL / "q24-2a5-zero-pole-sections-p103.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--precision", type=int, default=2000)
parser.add_argument("--checkpoint-digits")
parser.add_argument("--p230-source", required=True)
parser.add_argument("--q2-plus-q3-source", required=True)
parser.add_argument(
    "--output",
    default="artifacts/local/elkies-k3/q24-2a5-p230-bridge-hk-hensel-p103.json",
)
args = parser.parse_args()
assert args.precision >= 40


def resolved_path(value):
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_monic_pole(path, expected_degrees):
    payload = json.loads(path.read_text())
    assert tuple(payload["projective_chart"]["degrees_X_Y_Z"]) == expected_degrees
    checkpoints = payload["hensel"]["checkpoints"]
    assert len(checkpoints) >= 2
    previous = checkpoints[-2]["rational_reconstructions"]
    current = checkpoints[-1]["rational_reconstructions"]
    x_count = expected_degrees[0] + 1
    y_count = expected_degrees[1] + 1
    z_count = expected_degrees[2]
    start = x_count + y_count
    assert len(previous) == len(current) == start + z_count
    coefficients = current[start:start + z_count]
    assert all(value is not None for value in coefficients)
    assert previous[start:start + z_count] == coefficients
    return [QQ(value) for value in coefficients] + [QQ(1)]


output = resolved_path(args.output)
p230_source = resolved_path(args.p230_source)
q2_plus_q3_source = resolved_path(args.q2_plus_q3_source)
precision = args.precision
checkpoints_requested = (
    tuple(int(value) for value in args.checkpoint_digits.split(","))
    if args.checkpoint_digits else (
        max(20, precision // 8), precision // 4,
        precision // 2, precision - 10,
    )
)
assert tuple(sorted(set(checkpoints_requested))) == checkpoints_requested
assert all(0 < value < precision for value in checkpoints_requested)

started = time.monotonic()
surface = json.loads(SURFACE.read_text())
pool = json.loads(POOL.read_text())
assert surface["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR"
assert pool["status"] == "PASS_BOUNDED_MOD103_ZERO_POLE_SECTION_ENUMERATION"

L_QQ_coefficients = stable_monic_pole(q2_plus_q3_source, (6, 9, 1))
Z_QQ_coefficients = stable_monic_pole(p230_source, (8, 12, 2))

p = ZZ(103)
A_QQ = [QQ(value) for value in surface["child"]["minimal_A_coefficients_low_to_high"]]
B_QQ = [QQ(value) for value in surface["child"]["minimal_B_coefficients_low_to_high"]]

Kp = Qp(p, prec=precision, type="capped-rel")
RT = PolynomialRing(Kp, "T")
T = RT.gen()
A = RT([Kp(value) for value in A_QQ])
B = RT([Kp(value) for value in B_QQ])
L = RT([Kp(value) for value in L_QQ_coefficients])
Z = RT([Kp(value) for value in Z_QQ_coefficients])

sections = {entry["index"]: entry for entry in pool["sections"]}
q0_seed = sections[114]
initial = (
    q0_seed["X_coefficients_low_to_high"]
    + q0_seed["Y_coefficients_low_to_high"]
    + [101, 0, 61, 83, 47]
    + [58, 17, 4, 85, 66, 98]
)
assert len(initial) == 23


def bridge_polynomials(values):
    x0 = RT(list(values[:5]))
    y0 = RT(list(values[5:12]))
    H = RT(list(values[12:17]))
    K = RT(list(values[17:23]))
    XR = x0 * L**2 + H * Z
    YR = y0 * L**3 + H * K
    return x0, y0, H, K, XR, YR


def residual(values):
    x0, y0, _, _, XR, YR = bridge_polynomials(values)
    q0_identity = y0**2 - x0**3 - A * x0 - B
    r_identity = YR**2 - XR**3 - A * XR * L**4 - B * L**6
    return vector(
        Kp,
        [q0_identity[index] for index in range(13)]
        + [r_identity[index] for index in range(15)],
    )


def jacobian(values):
    x0, y0, H, K, XR, YR = bridge_polynomials(values)
    q0_x_derivative = -3 * x0**2 - A
    q0_y_derivative = 2 * y0
    r_x_derivative = -3 * XR**2 - A * L**4
    r_y_derivative = 2 * YR
    zero = RT.zero()
    columns = [
        (q0_x_derivative * T**shift, r_x_derivative * L**2 * T**shift)
        for shift in range(5)
    ] + [
        (q0_y_derivative * T**shift, r_y_derivative * L**3 * T**shift)
        for shift in range(7)
    ] + [
        (zero, (r_x_derivative * Z + r_y_derivative * K) * T**shift)
        for shift in range(5)
    ] + [
        (zero, r_y_derivative * H * T**shift) for shift in range(6)
    ]
    return matrix(
        Kp, 28, 23,
        lambda row, column: (
            columns[column][0][row]
            if row < 13 else columns[column][1][row - 13]
        ),
    )


# Independent modular construction audit.
F = GF(p)
RTF = PolynomialRing(F, "T")
TF = RTF.gen()
to_F = lambda values: RTF([F(value.numerator()) / F(value.denominator()) for value in values])
A_F, B_F = to_F(A_QQ), to_F(B_QQ)
L_F, Z_F = to_F(L_QQ_coefficients), to_F(Z_QQ_coefficients)
assert L_F == RTF([41, 1]) and Z_F == RTF([6, 3, 1])
x0_F, y0_F = RTF(initial[:5]), RTF(initial[5:12])
H_F, K_F = RTF(initial[12:17]), RTF(initial[17:23])
assert y0_F**2 == x0_F**3 + A_F * x0_F + B_F
XR_F = x0_F * L_F**2 + H_F * Z_F
YR_F = y0_F * L_F**3 + H_F * K_F
assert XR_F == RTF([4, 99, 58, 28, 72, 29, 56])
assert YR_F == RTF([88, 43, 4, 38, 42, 71, 36, 32, 23, 9])
assert YR_F**2 == XR_F**3 + A_F * XR_F * L_F**4 + B_F * L_F**6
q0_xd_F, q0_yd_F = -3 * x0_F**2 - A_F, 2 * y0_F
r_xd_F, r_yd_F = -3 * XR_F**2 - A_F * L_F**4, 2 * YR_F
zero_F = RTF.zero()
columns_F = [
    (q0_xd_F * TF**shift, r_xd_F * L_F**2 * TF**shift)
    for shift in range(5)
] + [
    (q0_yd_F * TF**shift, r_yd_F * L_F**3 * TF**shift)
    for shift in range(7)
] + [
    (zero_F, (r_xd_F * Z_F + r_yd_F * K_F) * TF**shift)
    for shift in range(5)
] + [
    (zero_F, r_yd_F * H_F * TF**shift) for shift in range(6)
]
J_F = matrix(
    F, 28, 23,
    lambda row, column: (
        columns_F[column][0][row]
        if row < 13 else columns_F[column][1][row - 13]
    ),
)
assert J_F.rank() == 23
pivot_rows = tuple(map(int, J_F.transpose().pivots()))
assert pivot_rows == tuple(range(12)) + tuple(range(13, 24))
J_square_F = matrix(F, [J_F.row(index) for index in pivot_rows])
determinant = int(J_square_F.det())
assert determinant == 31


def minimum_valuation(entries):
    valuations = [value.valuation() for value in entries if value]
    return min(valuations) if valuations else precision


values = vector(Kp, [Kp(value).add_bigoh(1) for value in initial])
iterations = []
known_precision = 1
while known_precision < precision:
    working_precision = min(2 * known_precision, precision)
    values = vector(Kp, [
        Kp(value.lift()).add_bigoh(working_precision) for value in values
    ])
    all_residual = residual(values)
    square_residual = vector(Kp, [all_residual[index] for index in pivot_rows])
    full_jacobian = jacobian(values)
    square_jacobian = matrix(Kp, [
        full_jacobian.row(index) for index in pivot_rows
    ])
    correction = square_jacobian.solve_right(-square_residual)
    values += correction
    iterations.append({
        "working_precision_p_adic_digits": working_precision,
        "minimum_square_residual_valuation_before": int(
            minimum_valuation(square_residual)
        ),
        "minimum_correction_valuation": int(minimum_valuation(correction)),
        "minimum_full_residual_valuation_after": int(
            minimum_valuation(residual(values))
        ),
    })
    known_precision = working_precision
final_residual_valuations = [
    int(value.valuation()) if value else precision for value in residual(values)
]
assert min(final_residual_valuations) >= precision - 8


def reconstruct_at_digits(value, digits):
    modulus = p**digits
    residue = ZZ(value.lift()) % modulus
    try:
        reconstruction = ZZ(residue).rational_reconstruction(modulus)
    except (ArithmeticError, ValueError):
        reconstruction = None
    return residue, reconstruction


ring_QQ = PolynomialRing(QQ, "T")
A_exact, B_exact = ring_QQ(A_QQ), ring_QQ(B_QQ)
L_exact, Z_exact = ring_QQ(L_QQ_coefficients), ring_QQ(Z_QQ_coefficients)
checkpoint_records = []
previous = None
exact_record = None
for digits in checkpoints_requested:
    modulus = p**digits
    residues = []
    reconstructions = []
    for value in values:
        residue, reconstruction = reconstruct_at_digits(value, digits)
        residues.append(str(residue))
        reconstructions.append(reconstruction)
    stable = 0 if previous is None else sum(
        current is not None and current == old
        for current, old in zip(reconstructions, previous)
    )
    exact_identity = False
    if all(value is not None for value in reconstructions):
        x0_exact = ring_QQ(reconstructions[:5])
        y0_exact = ring_QQ(reconstructions[5:12])
        H_exact = ring_QQ(reconstructions[12:17])
        K_exact = ring_QQ(reconstructions[17:23])
        XR_exact = x0_exact * L_exact**2 + H_exact * Z_exact
        YR_exact = y0_exact * L_exact**3 + H_exact * K_exact
        exact_identity = (
            y0_exact**2 == x0_exact**3 + A_exact * x0_exact + B_exact
            and YR_exact**2
            == XR_exact**3 + A_exact * XR_exact * L_exact**4 + B_exact * L_exact**6
        )
        if exact_identity:
            raw_XP = (
                K_exact**2 - 2 * x0_exact * L_exact**2 * Z_exact**2
                - H_exact * Z_exact**3
            )
            XP_exact, remainder_XP = raw_XP.quo_rem(L_exact**2)
            assert not remainder_XP
            raw_YP = (
                -y0_exact * L_exact * Z_exact**3
                + K_exact * (x0_exact * Z_exact**2 - XP_exact)
            )
            YP_exact, remainder_YP = raw_YP.quo_rem(L_exact)
            assert not remainder_YP
            assert (
                YP_exact**2
                == XP_exact**3 + A_exact * XP_exact * Z_exact**4
                + B_exact * Z_exact**6
            )
            exact_record = {
                "Q0": {
                    "X_coefficients_low_to_high": [str(value) for value in x0_exact.list()],
                    "Y_coefficients_low_to_high": [str(value) for value in y0_exact.list()],
                    "Z_coefficients_low_to_high": ["1"],
                    "exact_Weierstrass_identity": True,
                },
                "H_coefficients_low_to_high": [str(value) for value in H_exact.list()],
                "K_coefficients_low_to_high": [str(value) for value in K_exact.list()],
                "Q2_plus_Q3": {
                    "X_coefficients_low_to_high": [str(value) for value in XR_exact.list()],
                    "Y_coefficients_low_to_high": [str(value) for value in YR_exact.list()],
                    "Z_coefficients_low_to_high": [str(value) for value in L_exact.list()],
                    "exact_Weierstrass_identity": True,
                },
                "P230": {
                    "X_coefficients_low_to_high": [str(value) for value in XP_exact.list()],
                    "Y_coefficients_low_to_high": [str(value) for value in YP_exact.list()],
                    "Z_coefficients_low_to_high": [str(value) for value in Z_exact.list()],
                    "exact_Weierstrass_identity": True,
                    "exact_group_law_construction": "Q0 + (Q2+Q3)",
                },
            }
    checkpoint_records.append({
        "p_adic_digits": digits,
        "modulus_bits": int(modulus.nbits()),
        "non_null_rational_reconstructions": sum(
            value is not None for value in reconstructions
        ),
        "stable_from_previous_checkpoint": stable,
        "exact_bridge_identity": exact_identity,
        "rational_reconstructions_Q0X_Q0Y_H_K": [
            None if value is None else str(value) for value in reconstructions
        ],
        "residues_0_to_p_power_minus_1": residues,
    })
    previous = reconstructions

status = (
    "PASS_EXACT_QQ_P230_BRIDGE_HK_SECTION"
    if exact_record is not None
    else "PASS_REGULAR_P103_P230_BRIDGE_HK_HENSEL_DIAGNOSTIC"
)
input_paths = (SURFACE, POOL, p230_source, q2_plus_q3_source)
payload = {
    "schema": "elkies-k3.q24-2a5-p230-bridge-hk-hensel-p103.v1",
    "status": status,
    "software": "SageMath 10.9 (conda-forge pinned repository environment)",
    "prime": int(p),
    "ansatz": {
        "Q0": "-P1229",
        "R": "Q2+Q3",
        "P230": "Q0+R",
        "identities": [
            "XR=x(Q0)*L^2+H*Z",
            "YR=y(Q0)*L^3+H*K",
        ],
        "degree_bounds_H_K": [4, 5],
        "variables": 23,
        "coefficient_equations": 28,
        "mod103_jacobian_rank": int(J_F.rank()),
        "selected_independent_equation_degrees": list(pivot_rows),
        "selected_jacobian_determinant_mod103": determinant,
        "fixed_L_coefficients_low_to_high": [str(value) for value in L_QQ_coefficients],
        "fixed_Z_coefficients_low_to_high": [str(value) for value in Z_QQ_coefficients],
        "seed_Q0X_Q0Y_H_K_mod103": initial,
    },
    "hensel": {
        "working_precision_p_adic_digits": precision,
        "iterations": iterations,
        "final_residual_valuations_by_degree": final_residual_valuations,
        "checkpoints": checkpoint_records,
    },
    "exact_QQ": exact_record,
    "large_Groebner_required": False,
    "proof_boundary": (
        "The stable L and Z reconstructions are construction inputs. Promotion "
        "requires reconstructed Q0,H,K and literal QQ Weierstrass and group-law "
        "identities, all recorded only when exact_QQ is non-null."
    ),
    "inputs": {
        "paths": [str(path) for path in input_paths],
        "sha256": {str(path): sha256(path) for path in input_paths},
    },
    "elapsed_seconds": round(time.monotonic() - started, 6),
}
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A5A5P230BRIDGEHK|rank=23|det=31|precision={}|iterations={}|exact={}|"
    "status={}|output={}".format(
        precision, len(iterations), exact_record is not None, status, output
    ),
    flush=True,
)
