#!/usr/bin/env sage -python
"""Certify exact q4/orbit230 section and H0 from scaled X reconstruction.

status: EXACT
claim: exact QQ P230, exact two-dimensional H0(O+P230), no Groebner basis

The child model has denominator profile D_A=43*u^4, D_B=43^2*u^6.
After x'=u^2*x and y'=u^3*y, reconstruct only X' from the stored regular
p-adic branch.  The exact Weierstrass right side is a square, which recovers
Y' and avoids reconstructing its substantially larger coefficients.  The
resolved RR plane is then checked by literal collision congruences over QQ.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, inverse_mod, lcm, matrix, power_mod


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
SURFACE = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"
WORD = LOCAL / "q24-2a5-q4o230-horizontal-word.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--hensel-source",
    default=(
        "artifacts/local/elkies-k3/"
        "q24-2a5-p230-projective-hensel-p103-prec300000.json"
    ),
)
parser.add_argument(
    "--output",
    default="artifacts/local/elkies-k3/q24-2a5-p230-scaled-x-qq.json",
)
args = parser.parse_args()


def resolved_path(value):
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


started = time.monotonic()
hensel_path = resolved_path(args.hensel_source)
output = resolved_path(args.output)
surface = json.loads(SURFACE.read_text())
word = json.loads(WORD.read_text())
hensel = json.loads(hensel_path.read_text())
assert surface["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR"
assert word["status"] == "PASS_EXACT_Q24_2A5_Q4O230_LOW_POLE_HORIZONTAL_WORD"
assert not any(word["q4_divisor_correction"]["effective_vertical_NS_coordinates"])
assert hensel["section"] == "p230"
assert tuple(hensel["projective_chart"]["degrees_X_Y_Z"]) == (8, 12, 2)

R = PolynomialRing(QQ, "T")
T = R.gen()
A = R([QQ(value) for value in surface["child"]["minimal_A_coefficients_low_to_high"]])
B = R([QQ(value) for value in surface["child"]["minimal_B_coefficients_low_to_high"]])
denominator_A = lcm([value.denominator() for value in A])
denominator_B = lcm([value.denominator() for value in B])
assert denominator_B**2 == ZZ(43) * denominator_A**3
assert denominator_A % 43 == 0 and denominator_B % (43**2) == 0
u_A, exact_u_A = (denominator_A // 43).nth_root(4, truncate_mode=True)
u_B, exact_u_B = (denominator_B // (43**2)).nth_root(6, truncate_mode=True)
assert exact_u_A and exact_u_B and u_A == u_B
u = ZZ(u_A)
A_scaled = R([u**4 * value for value in A])
B_scaled = R([u**6 * value for value in B])
assert all(value.denominator().divides(43) for value in A_scaled)
assert all(value.denominator().divides(43**2) for value in B_scaled)

checkpoint = hensel["hensel"]["checkpoints"][-1]
p = ZZ(hensel["prime"])
digits = int(checkpoint["p_adic_digits"])
modulus = p**digits
residues = [ZZ(value) for value in checkpoint["residues_0_to_p_power_minus_1"]]
assert len(residues) == 22
Z = R([QQ(value) for value in hensel["projective_chart"]["fixed_Z_QQ_candidate"]])
assert Z.degree() == 2 and Z.is_monic()

X_scaled_coefficients = []
for residue in residues[:9]:
    scaled_residue = residue * power_mod(u, 2, modulus) % modulus
    X_scaled_coefficients.append(
        ZZ(scaled_residue).rational_reconstruction(modulus)
    )
X_scaled = R(X_scaled_coefficients)
square_rhs = X_scaled**3 + A_scaled * X_scaled * Z**4 + B_scaled * Z**6
assert square_rhs.degree() == 24 and square_rhs.is_square()
Y_scaled = R(square_rhs.sqrt())

# Select the unique sign matching the pinned modular branch.
F = GF(p)
RF = PolynomialRing(F, "T")
seed_Y = RF(hensel["projective_chart"]["seed_X_Y_Z"]["Y_coefficients_low_to_high"])
target_Y_scaled = seed_Y * F(u)**3


def reduce_polynomial(poly):
    return RF([
        F(value.numerator()) / F(value.denominator()) for value in R(poly).list()
    ])


if reduce_polynomial(Y_scaled) != target_Y_scaled:
    Y_scaled = -Y_scaled
assert reduce_polynomial(Y_scaled) == target_Y_scaled

X = R([value / u**2 for value in X_scaled])
Y = R([value / u**3 for value in Y_scaled])
assert (X.degree(), Y.degree(), Z.degree()) == (8, 12, 2)
assert Y**2 == X**3 + A * X * Z**4 + B * Z**6

# A substantial p-adic branch agreement audit, beyond the modulo-p sign.  The
# exact QQ identities carry the proof; there is no need to repeat inversions at
# the full two-million-bit reconstruction modulus.
agreement_digits = min(digits, 1000)
agreement_modulus = p**agreement_digits


def residue_mod(value):
    value = QQ(value)
    assert value.denominator() % p
    return (
        ZZ(value.numerator())
        * inverse_mod(ZZ(value.denominator()), agreement_modulus)
        % agreement_modulus
    )


assert [residue_mod(value) for value in X.list()] == [
    value % agreement_modulus for value in residues[:9]
]
assert [residue_mod(value) for value in Y.list()] == [
    value % agreement_modulus for value in residues[9:22]
]

# Exact canonical RR plane and literal H0 collision calculation.
rr_strings = checkpoint["resolved_RR"][
    "rational_reconstructions_AA2_AA3_AA4_BB_for_two_rows"
]
assert len(rr_strings) == 8 and all(value is not None for value in rr_strings)
rr_values = [QQ(value) for value in rr_strings]
AA0 = R([1, 0] + rr_values[:3])
BB0 = rr_values[3]
AA1 = R([0, 1] + rr_values[4:7])
BB1 = rr_values[7]
collision_modulus = Z**2
assert not (AA0 * X - BB0 * Y) % collision_modulus
assert not (AA1 * X - BB1 * Y) % collision_modulus

ambient_pairs = (
    [(T**degree, R.zero()) for degree in range(5)]
    + [(R.zero(), R.one())]
)
collision_remainders = [
    R((AA * X - BB * Y) % collision_modulus) for AA, BB in ambient_pairs
]
collision_matrix = matrix(QQ, [
    [remainder[degree] for remainder in collision_remainders]
    for degree in range(4)
])
assert collision_matrix.rank() == 4
canonical_rows = matrix(QQ, [
    [AA0[index] for index in range(5)] + [BB0],
    [AA1[index] for index in range(5)] + [BB1],
])
assert canonical_rows.rank() == 2
assert collision_matrix * canonical_rows.transpose() == 0

determinant = AA0 * BB1 - AA1 * BB0
determinant_ratio, determinant_remainder = determinant.quo_rem(Z**2)
assert not determinant_remainder and determinant_ratio.degree() == 0
constant_C = R(AA0 / BB0)
constant_C -= constant_C[4] * Z**2
assert constant_C.degree() <= 3
constant_rows = matrix(QQ, [
    [(Z**2)[index] for index in range(5)] + [0],
    [constant_C[index] for index in range(5)] + [1],
])
assert constant_rows.row_space() == canonical_rows.row_space()


def height_profile(polynomials):
    values = [value for poly in polynomials for value in R(poly).list()]
    return {
        "maximum_numerator_bits": int(max(abs(value.numerator()).nbits() for value in values)),
        "maximum_denominator_bits": int(max(value.denominator().nbits() for value in values)),
        "maximum_rational_bits": int(max(
            max(abs(value.numerator()).nbits(), value.denominator().nbits())
            for value in values
        )),
    }


input_paths = (SURFACE, WORD, hensel_path)
payload = {
    "schema": "elkies-k3.q24-2a5-p230-scaled-x-qq.v1",
    "status": "PASS_EXACT_QQ_P230_SECTION_AND_H0",
    "software": "SageMath 10.9 (conda-forge pinned repository environment)",
    "method": {
        "surface_scaling": "x'=u^2*x, y'=u^3*y, A'=u^4*A, B'=u^6*B",
        "denominator_identities": ["D_A=43*u^4", "D_B=43^2*u^6"],
        "u": str(u),
        "u_bits": int(u.nbits()),
        "reconstructed_coordinates": "X' only",
        "Y_recovery": "exact polynomial square root with sign selected mod 103",
        "large_Groebner_required": False,
    },
    "P230": {
        "X_coefficients_low_to_high": [str(value) for value in X.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y.list()],
        "Z_coefficients_low_to_high": [str(value) for value in Z.list()],
        "degrees_X_Y_Z": [8, 12, 2],
        "P_dot_O": 2,
        "exact_Weierstrass_identity": True,
        "p_adic_branch_agreement_digits": agreement_digits,
        "height_profile": height_profile((X, Y, Z)),
        "scaled_height_profile": height_profile((X_scaled, Y_scaled, Z)),
    },
    "H0": {
        "divisor": "O + P230",
        "vertical_correction": "zero",
        "ambient_pairs": "deg(AA)<=4, BB constant",
        "ambient_dimension": 6,
        "collision_modulus": "Z^2",
        "collision_rank": int(collision_matrix.rank()),
        "dimension": 2,
        "canonical_basis": [
            {
                "AA_coefficients_low_to_high": [str(AA0[index]) for index in range(5)],
                "BB": str(BB0),
            },
            {
                "AA_coefficients_low_to_high": [str(AA1[index]) for index in range(5)],
                "BB": str(BB1),
            },
        ],
        "constant_function_basis": [
            {
                "AA_coefficients_low_to_high": [str((Z**2)[index]) for index in range(5)],
                "BB": "0",
            },
            {
                "AA_coefficients_low_to_high": [str(constant_C[index]) for index in range(5)],
                "BB": "1",
            },
        ],
        "determinant_identity": "AA0*BB1-AA1*BB0 = scalar*Z^2",
        "determinant_scalar": str(determinant_ratio[0]),
        "exact_collision_congruences": True,
    },
    "proof_boundary": (
        "This certifies the exact P230 equation and the complete two-dimensional "
        "H0(O+P230). Quartic/Jacobian compilation and child fibre verification "
        "remain separate downstream gates."
    ),
    "inputs": {
        "paths": [str(path) for path in input_paths],
        "sha256": {str(path): sha256(path) for path in input_paths},
    },
    "elapsed_seconds": round(time.monotonic() - started, 6),
}
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A5A5P230QQ|Xbits={}|Ybits={}|H0=2|collision_rank=4|status={}|output={}".format(
        payload["P230"]["height_profile"]["maximum_rational_bits"],
        payload["P230"]["scaled_height_profile"]["maximum_rational_bits"],
        payload["status"], output,
    ),
    flush=True,
)
