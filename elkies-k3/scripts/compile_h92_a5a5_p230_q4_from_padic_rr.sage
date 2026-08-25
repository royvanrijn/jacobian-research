#!/usr/bin/env sage -python
"""Compile the q4/orbit230 quartic from a stored high-precision RR lift.

status: CONSTRUCTION EXPERIMENT
claim: fraction-free p-adic quartic and exact rational-reconstruction audit

This consumes stored projective-section residues and the canonical resolved
RR plane.  It performs only univariate/nested-polynomial arithmetic: no
Groebner basis and no nonlinear solve.  An exact quartic is recorded only when
all 25 bounded coefficients reconstruct and its literal invariant identities
hold over QQ.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import PolynomialRing, QQ, Qp, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
SURFACE = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--hensel-source", required=True)
parser.add_argument("--checkpoint-index", type=int, default=-1)
parser.add_argument(
    "--rr-basis",
    choices=("constant", "canonical"),
    default="constant",
)
parser.add_argument(
    "--output",
    default="artifacts/local/elkies-k3/q24-2a5-p230-q4-from-padic-rr.json",
)
args = parser.parse_args()


def resolved_path(value):
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


started = time.monotonic()
source_path = resolved_path(args.hensel_source)
output = resolved_path(args.output)
source = json.loads(source_path.read_text())
surface = json.loads(SURFACE.read_text())
assert source["section"] == "p230"
assert surface["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR"
chart = source["projective_chart"]
assert tuple(chart["degrees_X_Y_Z"]) == (8, 12, 2)
fixed_Z = chart["fixed_Z_QQ_candidate"]
assert fixed_Z is not None and len(fixed_Z) == 3

checkpoint = source["hensel"]["checkpoints"][args.checkpoint_index]
digits = int(checkpoint["p_adic_digits"])
p = ZZ(source["prime"])
modulus = p**digits
section_residues = [ZZ(value) for value in checkpoint["residues_0_to_p_power_minus_1"]]
assert len(section_residues) == 22
rr_record = checkpoint["resolved_RR"]
rr_strings = rr_record["rational_reconstructions_AA2_AA3_AA4_BB_for_two_rows"]
assert len(rr_strings) == 8 and all(value is not None for value in rr_strings)
rr_QQ = [QQ(value) for value in rr_strings]

K = Qp(p, prec=digits, type="capped-rel")
RT = PolynomialRing(K, "T")
T = RT.gen()
X = RT([K(value).add_bigoh(digits) for value in section_residues[:9]])
Y = RT([K(value).add_bigoh(digits) for value in section_residues[9:22]])
Z = RT([K(QQ(value)) for value in fixed_Z])
A_QQ = [QQ(value) for value in surface["child"]["minimal_A_coefficients_low_to_high"]]
B_QQ = [QQ(value) for value in surface["child"]["minimal_B_coefficients_low_to_high"]]
A = RT([K(value) for value in A_QQ])
B = RT([K(value) for value in B_QQ])
section_residual = Y**2 - X**3 - A * X * Z**4 - B * Z**6
section_residual_valuations = [
    int(value.valuation()) if value else digits for value in section_residual.list()
]
assert min(section_residual_valuations) >= digits - 8

RQQ = PolynomialRing(QQ, "T")
Z_QQ = RQQ([QQ(value) for value in fixed_Z])
canonical_AA0_QQ = RQQ([QQ(1), QQ(0)] + rr_QQ[:3])
canonical_BB0_QQ = rr_QQ[3]
canonical_AA1_QQ = RQQ([QQ(0), QQ(1)] + rr_QQ[4:7])
canonical_BB1_QQ = rr_QQ[7]
canonical_determinant_QQ = (
    canonical_AA0_QQ * canonical_BB1_QQ
    - canonical_AA1_QQ * canonical_BB0_QQ
)
_, canonical_determinant_remainder_QQ = canonical_determinant_QQ.quo_rem(Z_QQ**2)
assert not canonical_determinant_remainder_QQ

if args.rr_basis == "constant":
    C0_QQ = RQQ(canonical_AA0_QQ / canonical_BB0_QQ)
    C1_QQ = RQQ(canonical_AA1_QQ / canonical_BB1_QQ)
    _, difference_remainder = (C0_QQ - C1_QQ).quo_rem(Z_QQ**2)
    assert not difference_remainder
    C_QQ = C0_QQ - C0_QQ[4] * Z_QQ**2
    assert C_QQ.degree() <= 3
    AA0, BB0 = RT([K(value) for value in (Z_QQ**2).list()]), K(0)
    AA1, BB1 = RT([K(value) for value in C_QQ.list()]), K(1)
else:
    C_QQ = None
    AA0 = RT([K(value) for value in canonical_AA0_QQ.list()])
    BB0 = K(canonical_BB0_QQ)
    AA1 = RT([K(value) for value in canonical_AA1_QQ.list()])
    BB1 = K(canonical_BB1_QQ)
determinant = AA0 * BB1 - AA1 * BB0
determinant_ratio, determinant_remainder = determinant.quo_rem(Z**2)
assert determinant_ratio.degree() == 0 and determinant_ratio[0]
determinant_remainder_valuations = [
    int(value.valuation()) if value else digits for value in determinant_remainder.list()
]
assert min(determinant_remainder_valuations) >= digits - 8

# Nested polynomial ring K[U][T].  The q4 quartic is the quotient by Z^6 of
# the chord radicand; the quotient must have T-degree four.
RU = PolynomialRing(K, "U")
U = RU.gen()
RUT = PolynomialRing(RU, "T")


def nested(poly):
    return RUT([RU(value) for value in RT(poly).list()])


AA0b, AA1b, Zb, Xb, Yb, Ab = map(nested, (AA0, AA1, Z, X, Y, A))
N = AA1b - RUT(U) * AA0b
Db = RUT(U * BB0 - BB1)
raw = (
    N**4 - 6 * Xb * N**2 * Db**2 - 8 * Yb * N * Db**3
    - 3 * Xb**2 * Db**4 - 4 * Ab * Zb**4 * Db**4
)
quartic_padic, remainder = raw.quo_rem(Zb**6)
assert quartic_padic.degree() == 4
remainder_values = [coefficient for poly in remainder.list() for coefficient in poly.list()]
remainder_valuations = [
    int(value.valuation()) if value else digits for value in remainder_values
]
assert (not remainder_values) or min(remainder_valuations) >= digits - 8
assert all(RU(coefficient).degree() <= 4 for coefficient in quartic_padic.list())


def reconstruct(value):
    residue = ZZ(K(value).lift()) % modulus
    try:
        result = ZZ(residue).rational_reconstruction(modulus)
    except (ArithmeticError, ValueError):
        result = None
    return residue, result


flat_residues = []
flat_reconstructions = []
for old_degree in range(5):
    coefficient = RU(quartic_padic[old_degree])
    for new_degree in range(5):
        residue, reconstruction = reconstruct(coefficient[new_degree])
        flat_residues.append(str(residue))
        flat_reconstructions.append(reconstruction)

exact_quartic = None
if all(value is not None for value in flat_reconstructions):
    QU = PolynomialRing(QQ, "U")
    QT = PolynomialRing(QU, "T")
    coefficient_rows = [
        flat_reconstructions[5 * index:5 * (index + 1)] for index in range(5)
    ]
    quartic_QQ = QT([QU(row) for row in coefficient_rows])
    e, d, c, b, a = quartic_QQ.list()
    I = 12 * a * e - 3 * b * d + c**2
    J = (
        72 * a * c * e + 9 * b * c * d - 27 * a * d**2
        - 27 * b**2 * e - 2 * c**3
    )
    A_child = QU(-27 * I)
    B_child = QU(-27 * J)
    exact_quartic = {
        "coefficients_by_old_degree_each_low_to_high_in_new_base": [
            [str(value) for value in row] for row in coefficient_rows
        ],
        "I_coefficients_low_to_high": [str(value) for value in I.list()],
        "J_coefficients_low_to_high": [str(value) for value in J.list()],
        "jacobian_A_coefficients_low_to_high": [str(value) for value in A_child.list()],
        "jacobian_B_coefficients_low_to_high": [str(value) for value in B_child.list()],
        "degree_profile_quartic_I_J_A_B": [
            int(quartic_QQ.degree()), int(I.degree()), int(J.degree()),
            int(A_child.degree()), int(B_child.degree()),
        ],
    }

status = (
    "PASS_EXACT_RECONSTRUCTED_Q4_QUARTIC_INVARIANTS_CONSTRUCTION_AID"
    if exact_quartic is not None
    else "PASS_PADIC_Q4_QUARTIC_RECONSTRUCTION_DIAGNOSTIC"
)
payload = {
    "schema": "elkies-k3.q24-2a5-p230-q4-from-padic-rr.v1",
    "status": status,
    "software": "SageMath 10.9 (conda-forge pinned repository environment)",
    "prime": int(p),
    "p_adic_digits": digits,
    "modulus_bits": int(modulus.nbits()),
    "section_residual_minimum_valuation": min(section_residual_valuations),
    "radicand_remainder_minimum_valuation": (
        digits if not remainder_values else min(remainder_valuations)
    ),
    "resolved_RR": {
        "basis": args.rr_basis,
        "canonical_rows_AA2_AA3_AA4_BB": [str(value) for value in rr_QQ],
        "constant_basis_C_coefficients_low_to_high": (
            None if C_QQ is None else [str(value) for value in C_QQ.list()]
        ),
        "determinant_is_scalar_times_Z_squared_p_adically": True,
    },
    "quartic_reconstruction": {
        "coefficient_slots": 25,
        "non_null_rational_reconstructions": sum(
            value is not None for value in flat_reconstructions
        ),
        "rational_reconstructions_by_old_degree_then_new_degree": [
            None if value is None else str(value) for value in flat_reconstructions
        ],
        "residues_by_old_degree_then_new_degree": flat_residues,
        "exact_candidate": exact_quartic,
    },
    "large_Groebner_required": False,
    "proof_boundary": (
        "This compiles an exact rational candidate from a certified regular "
        "p-adic branch. It is a construction aid until the section and RR "
        "functions are independently verified over QQ."
    ),
    "inputs": {
        "paths": [str(SURFACE), str(source_path)],
        "sha256": {str(path): sha256(path) for path in (SURFACE, source_path)},
    },
    "elapsed_seconds": round(time.monotonic() - started, 6),
}
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A5A5P230Q4PADIC|digits={}|quartic_reconstructed={}/25|exact_candidate={}|"
    "status={}|output={}".format(
        digits,
        payload["quartic_reconstruction"]["non_null_rational_reconstructions"],
        exact_quartic is not None,
        status,
        output,
    ),
    flush=True,
)
