#!/usr/bin/env sage -python
"""Recover exact polynomial P1229 from its regular p-adic branch.

Use the already-certified global scaling x'=u^2*x, y'=u^3*y.  Reconstruct
only the five coefficients of X', then recover Y' as the exact square root of
X'^3+A'X'+B'.  This is a low-degree Hensel/reconstruction certificate and
uses no Groebner basis.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, inverse_mod, power_mod


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
SURFACE = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"
WORD = LOCAL / "q24-2a5-q6o1307-horizontal-word.json"
SCALE = LOCAL / "q24-2a5-p230-scaled-x-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--hensel-source",
    default="artifacts/local/elkies-k3/q24-2a5-p1229-hensel-p103-prec300000.json",
)
parser.add_argument(
    "--output",
    default="artifacts/local/elkies-k3/q24-2a5-p1229-scaled-x-qq.json",
)
args = parser.parse_args()


def resolved(value):
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


hensel_path = resolved(args.hensel_source)
output = resolved(args.output)
surface = json.loads(SURFACE.read_text())
word = json.loads(WORD.read_text())
scale = json.loads(SCALE.read_text())
hensel = json.loads(hensel_path.read_text())
assert surface["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR"
assert word["status"] == "PASS_EXACT_Q24_2A5_Q6O1307_LOW_POLE_HORIZONTAL_WORD"
assert scale["status"] == "PASS_EXACT_QQ_P230_SECTION_AND_H0"
assert hensel["status"] == "PASS_TARGETED_REGULAR_MOD103_HENSEL_LIFTS"
assert hensel["hensel"]["P1229_full_lift_indices"] == [115]

R = PolynomialRing(QQ, "T")
A = R([QQ(value) for value in surface["child"]["minimal_A_coefficients_low_to_high"]])
B = R([QQ(value) for value in surface["child"]["minimal_B_coefficients_low_to_high"]])
u = ZZ(scale["method"]["u"])
A_scaled = R([u**4 * value for value in A])
B_scaled = R([u**6 * value for value in B])

probe = hensel["hensel"]["unique_selected_probe"]
assert probe["section_index"] == 115 and probe["full_thirteen_equation_hensel_lift"]
checkpoint = probe["checkpoints"][-1]
p = ZZ(hensel["prime"])
digits = int(checkpoint["p_adic_digits"])
modulus = p**digits
residues = [ZZ(value) for value in checkpoint["residues_0_to_p_power_minus_1"]]
assert len(residues) == 12

X_scaled = R([
    ZZ(residue * power_mod(u, 2, modulus) % modulus).rational_reconstruction(modulus)
    for residue in residues[:5]
])
square_rhs = X_scaled**3 + A_scaled * X_scaled + B_scaled
assert square_rhs.degree() == 12 and square_rhs.is_square()
Y_scaled = R(square_rhs.sqrt())

F = GF(p)
RF = PolynomialRing(F, "T")
seed_Y = RF(probe["Y_coefficients_low_to_high_mod103"])
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
assert (X.degree(), Y.degree()) == (4, 6)
assert Y**2 == X**3 + A * X + B

# Audit a substantial p-adic prefix in addition to the exact identity.
agreement_digits = min(digits, 1000)
agreement_modulus = p**agreement_digits


def residue_mod(value):
    value = QQ(value)
    return (
        ZZ(value.numerator())
        * inverse_mod(ZZ(value.denominator()), agreement_modulus)
        % agreement_modulus
    )


assert [residue_mod(value) for value in X.list()] == [
    value % agreement_modulus for value in residues[:5]
]
assert [residue_mod(value) for value in Y.list()] == [
    value % agreement_modulus for value in residues[5:]
]


def height_profile(polynomials):
    values = [value for polynomial in polynomials for value in polynomial.list()]
    return {
        "maximum_numerator_bits": int(max(abs(value.numerator()).nbits() for value in values)),
        "maximum_denominator_bits": int(max(value.denominator().nbits() for value in values)),
        "maximum_rational_bits": int(max(
            max(abs(value.numerator()).nbits(), value.denominator().nbits()) for value in values
        )),
    }


inputs = (SURFACE, WORD, SCALE, hensel_path)
payload = {
    "schema": "elkies-k3.q24-2a5-p1229-scaled-x-qq.v1",
    "status": "PASS_EXACT_QQ_P1229_POLYNOMIAL_SECTION",
    "method": {
        "surface_scaling": "x'=u^2*x, y'=u^3*y, A'=u^4*A, B'=u^6*B",
        "u_bits": int(u.nbits()),
        "reconstructed_coordinates": "X' only",
        "Y_recovery": "exact polynomial square root with sign selected mod 103",
        "large_Groebner_required": False,
    },
    "P1229": {
        "X_coefficients_low_to_high": [str(value) for value in X.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y.list()],
        "degrees_X_Y_Z": [4, 6, 0],
        "P_dot_O": 0,
        "exact_Weierstrass_identity": True,
        "p_adic_branch_agreement_digits": agreement_digits,
        "NS_coordinates": word["sections"]["q6_orbit1229"]["effective_section"],
        "height_profile": height_profile((X, Y)),
        "scaled_height_profile": height_profile((X_scaled, Y_scaled)),
    },
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in inputs
        },
    },
}
output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "P1229QQ|X=4|Y=6|Z=0|bits={}|status={}|output={}".format(
        payload["P1229"]["height_profile"]["maximum_rational_bits"],
        payload["status"], output,
    ),
    flush=True,
)
