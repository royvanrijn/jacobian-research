#!/usr/bin/env sage -python
"""Recover the physical q10 horizontal section by a short exact MW word.

status: ACTIVE_PROOF
claim: exact QQ horizontal section for the physical-nef 2A5-to-3A3 q10 edge
outputs: artifacts/local/elkies-k3/q24-2a5-direct-q10-horizontal-qq.json

Only the regular zero-pole branch 112 is newly reconstructed.  The target is

  2 P_aff - 2 P1 + 3 P1229 + 2 P230 + 2 P14.

All remaining points are already exact.  Elliptic group law replaces a direct
singular projective lift; no Groebner basis is used.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
SURFACE = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"
Q10 = LOCAL / "q24-2a5-direct-physical-q10-certificate.json"
SHELL = LOCAL / "q24-2a5-zero-pole-shell-match-p103.json"
POOL = LOCAL / "q24-2a5-zero-pole-sections-p103.json"
Q4_WORD = LOCAL / "q24-2a5-q4o230-horizontal-word.json"
P230 = LOCAL / "q24-2a5-p230-scaled-x-qq.json"
P1229 = LOCAL / "q24-2a5-p1229-scaled-x-qq.json"
P146 = LOCAL / "q24-2a5-p146-p1307-scaled-x-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--hensel-source",
    default="artifacts/local/elkies-k3/q24-2a5-p14-hensel-p103-prec300000.json",
)
parser.add_argument(
    "--output",
    default="artifacts/local/elkies-k3/q24-2a5-direct-q10-horizontal-qq.json",
)
args = parser.parse_args()


def resolved(value):
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


hensel_path = resolved(args.hensel_source)
output = resolved(args.output)
surface = json.loads(SURFACE.read_text())
q10 = json.loads(Q10.read_text())
shell = json.loads(SHELL.read_text())
pool = json.loads(POOL.read_text())
q4_word = json.loads(Q4_WORD.read_text())
p230 = json.loads(P230.read_text())
p1229 = json.loads(P1229.read_text())
p146 = json.loads(P146.read_text())
hensel = json.loads(hensel_path.read_text())
assert surface["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR"
assert q10["status"] == "PASS_EXACT_PHYSICAL_NEF_Q10_CURRENT_3A3_PRESENTATION"
assert shell["status"] == "PASS_EXHAUSTIVE_MOD103_ZERO_POLE_SHELL_EMBEDDINGS_CANONICAL_MARKING"
assert pool["status"] == "PASS_BOUNDED_MOD103_ZERO_POLE_SECTION_ENUMERATION"
assert p230["status"] == "PASS_EXACT_QQ_P230_SECTION_AND_H0"
assert p1229["status"] == "PASS_EXACT_QQ_P1229_POLYNOMIAL_SECTION"
assert p146["status"] == "PASS_EXACT_QQ_P146_AND_P1307_SHORT_MW_WORDS"
assert hensel["status"] == "PASS_TARGETED_REGULAR_MOD103_HENSEL_LIFTS"
assert hensel["hensel"]["full_lift_indices"] == [112]

R = PolynomialRing(QQ, "T")
K = R.fraction_field()
T = R.gen()
A = R([QQ(value) for value in surface["child"]["minimal_A_coefficients_low_to_high"]])
B = R([QQ(value) for value in surface["child"]["minimal_B_coefficients_low_to_high"]])
E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])
u = ZZ(p230["method"]["u"])
A_scaled = u**4 * A
B_scaled = u**6 * B
p = ZZ(hensel["prime"])
F = GF(p)
RF = PolynomialRing(F, "T")


def reduce_polynomial(polynomial):
    return RF([
        F(value.numerator()) / F(value.denominator())
        for value in R(polynomial).list()
    ])


probe = hensel["hensel"]["candidate_probes"][0]
assert probe["section_index"] == 112
assert probe["mod103_jacobian_rank"] == 12
checkpoint = probe["checkpoints"][-1]
digits = int(checkpoint["p_adic_digits"])
modulus = p**digits
residues = [ZZ(value) for value in checkpoint["residues_0_to_p_power_minus_1"]]
assert len(residues) == 12
X14_scaled = R([
    ZZ(residue * power_mod(u, 2, modulus) % modulus).rational_reconstruction(modulus)
    for residue in residues[:5]
])
rhs14_scaled = X14_scaled**3 + A_scaled * X14_scaled + B_scaled
assert rhs14_scaled.degree() == 12 and rhs14_scaled.is_square()
Y14_scaled = R(rhs14_scaled.sqrt())
seed14 = pool["sections"][112]
target_Y14_scaled = RF(seed14["Y_coefficients_low_to_high"]) * F(u)**3
if reduce_polynomial(Y14_scaled) != target_Y14_scaled:
    Y14_scaled = -Y14_scaled
assert reduce_polynomial(Y14_scaled) == target_Y14_scaled
X14 = R([value / u**2 for value in X14_scaled])
Y14 = R([value / u**3 for value in Y14_scaled])
assert (X14.degree(), Y14.degree()) == (4, 6)
assert Y14**2 == X14**3 + A * X14 + B
P14 = E(K(X14), K(Y14))


def point_from_payload(record):
    X = R([QQ(value) for value in record["X_coefficients_low_to_high"]])
    Y = R([QQ(value) for value in record["Y_coefficients_low_to_high"]])
    Z = R([QQ(value) for value in record.get("Z_coefficients_low_to_high", ["1"])])
    assert Y**2 == X**3 + A * X * Z**4 + B * Z**6
    return E(K(X / Z**2), K(Y / Z**3))


Paff = point_from_payload(p146["polynomial_inputs"]["old_A11_affine"])
P1 = point_from_payload(p146["polynomial_inputs"]["P1"])
P1229_point = point_from_payload(p1229["P1229"])
P230_point = point_from_payload(p230["P230"])
Pdirect = 2 * Paff - 2 * P1 + 3 * P1229_point + 2 * P230_point + 2 * P14


def projectivize(point):
    x, y = point.xy()
    denominator_x = R(x.denominator()).monic()
    assert denominator_x.is_square()
    Z = R(denominator_x.sqrt()).monic()
    assert R(y.denominator()).monic() == Z**3
    X = R(x * Z**2)
    Y = R(y * Z**3)
    assert (X.degree(), Y.degree(), Z.degree()) == (14, 21, 5)
    assert Y**2 == X**3 + A * X * Z**4 + B * Z**6
    return X, Y, Z


Xdirect, Ydirect, Zdirect = projectivize(Pdirect)

# Full NS quotient identity, including the exact trivial-lattice correction.
classes = shell["exact_shell"]["indexed_component_nef_classes"]
target_NS = vector(ZZ, q10["RR_profile"]["horizontal_section"])
P230_NS = vector(ZZ, q4_word["q4_orbit230_horizontal"]["effective_section"])
word_NS = (
    2 * vector(ZZ, classes[45]["NS_coordinates"])
    - 2 * vector(ZZ, classes[1]["NS_coordinates"])
    + 3 * vector(ZZ, classes[24]["NS_coordinates"])
    + 2 * P230_NS
    + 2 * vector(ZZ, classes[14]["NS_coordinates"])
)
correction = target_NS - word_NS
assert entries(correction) == [
    -5, -6, 1, 0, 1, 2, 1, -4, -2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0,
]

# Independent compatible mod-103 group word.  Mapping 35 is the unique one
# that simultaneously gives P1229=115 and the regular P230 triple.
mapping = shell["all_complete_mappings_exact_index_to_modular_index"][35]
assert [mapping[str(index)] for index in (45, 1, 24, 14)] == [49, 30, 115, 112]
assert [mapping[str(index)] for index in (0, 2, 3)] == [114, 62, 36]
EF = EllipticCurve(RF.fraction_field(), [
    0, 0, 0,
    RF(pool["surface_mod_103"]["A_coefficients_low_to_high"]),
    RF(pool["surface_mod_103"]["B_coefficients_low_to_high"]),
])


def modular_point(index):
    record = pool["sections"][index]
    return EF(
        RF(record["X_coefficients_low_to_high"]),
        RF(record["Y_coefficients_low_to_high"]),
    )


P230_mod = modular_point(114) + modular_point(62) + modular_point(36)
Pdirect_mod = (
    2 * modular_point(49) - 2 * modular_point(30) + 3 * modular_point(115)
    + 2 * P230_mod + 2 * modular_point(112)
)
x_mod, y_mod = Pdirect_mod.xy()
Z_mod = RF(x_mod.denominator()).monic().sqrt().monic()
assert RF(y_mod.denominator()).monic() == Z_mod**3
assert reduce_polynomial(Xdirect) == RF(x_mod * Z_mod**2)
assert reduce_polynomial(Ydirect) == RF(y_mod * Z_mod**3)
assert reduce_polynomial(Zdirect) == Z_mod


def height_profile(polynomials):
    values = [value for polynomial in polynomials for value in R(polynomial).list()]
    return {
        "maximum_numerator_bits": int(max(abs(value.numerator()).nbits() for value in values)),
        "maximum_denominator_bits": int(max(value.denominator().nbits() for value in values)),
        "maximum_rational_bits": int(max(
            max(abs(value.numerator()).nbits(), value.denominator().nbits())
            for value in values
        )),
    }


inputs = (SURFACE, Q10, SHELL, POOL, Q4_WORD, P230, P1229, P146, hensel_path)
payload = {
    "schema": "elkies-k3.q24-2a5-direct-q10-horizontal-qq.v1",
    "status": "PASS_EXACT_QQ_PHYSICAL_Q10_HORIZONTAL_SECTION",
    "method": {
        "MW_word": "2*P_aff-2*P1+3*P1229+2*P230+2*P14",
        "new_regular_mod103_branch": 112,
        "P14_exact_shell_index": 14,
        "compatible_shell_mapping_index": 35,
        "target_minus_word_NS_coordinates": entries(correction),
        "direct_projective_Hensel_lift_avoided": True,
        "large_Groebner_required": False,
    },
    "P14": {
        "X_coefficients_low_to_high": [str(value) for value in X14.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y14.list()],
        "degrees_X_Y_Z": [4, 6, 0],
        "height_profile": height_profile((X14, Y14)),
        "p_adic_branch_agreement_digits": min(digits, 1000),
    },
    "physical_q10_horizontal": {
        "X_coefficients_low_to_high": [str(value) for value in Xdirect.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Ydirect.list()],
        "Z_coefficients_low_to_high": [str(value) for value in Zdirect.list()],
        "degrees_X_Y_Z": [14, 21, 5],
        "P_dot_O": 5,
        "NS_coordinates": entries(target_NS),
        "exact_Weierstrass_identity": True,
        "independent_mod103_group_word_match": True,
        "height_profile": height_profile((Xdirect, Ydirect, Zdirect)),
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
    "A5DIRECTQ10QQ|P14_bits={}|horizontal_bits={}|degrees=14,21,5|"
    "status={}|output={}".format(
        payload["P14"]["height_profile"]["maximum_rational_bits"],
        payload["physical_q10_horizontal"]["height_profile"]["maximum_rational_bits"],
        payload["status"], output,
    ),
    flush=True,
)
