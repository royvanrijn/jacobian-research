#!/usr/bin/env sage -python
"""Recover exact P146 and P1307 by short Mordell--Weil words.

The missing simple-pole direction is not regularly liftable from its direct
mod-103 projective seed.  Instead use the exact old-A11 affine section and
the certified lattice identity

    P146 = P_affine + P1 + P32  (modulo the trivial lattice).

P1 and P32 are regular polynomial branches.  Reconstruct only their scaled
X coordinates, recover Y by an exact square root, and use exact elliptic
addition.  Finally form P1307=P146+P1229.  No Groebner basis is used.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    EllipticCurve, GF, PolynomialRing, QQ, ZZ, inverse_mod, power_mod, vector,
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
SURFACE = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"
MARKING = LOCAL / "q24-a11-to-2a5-q8-equation-marking-qq.json"
SHELL = LOCAL / "q24-2a5-zero-pole-shell-match-p103.json"
WORD = LOCAL / "q24-2a5-q6o1307-horizontal-word.json"
SCALE = LOCAL / "q24-2a5-p230-scaled-x-qq.json"
P1229 = LOCAL / "q24-2a5-p1229-scaled-x-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--hensel-source",
    default="artifacts/local/elkies-k3/q24-2a5-p1-p32-hensel-p103-prec300000.json",
)
parser.add_argument(
    "--output",
    default="artifacts/local/elkies-k3/q24-2a5-p146-p1307-scaled-x-qq.json",
)
args = parser.parse_args()


def resolved(value):
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


hensel_path = resolved(args.hensel_source)
output = resolved(args.output)
surface = json.loads(SURFACE.read_text())
marking = json.loads(MARKING.read_text())
shell = json.loads(SHELL.read_text())
word = json.loads(WORD.read_text())
scale = json.loads(SCALE.read_text())
p1229 = json.loads(P1229.read_text())
hensel = json.loads(hensel_path.read_text())
assert surface["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_RESOLVED_RR"
assert marking["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_EQUATION_MARKING"
assert shell["status"] == "PASS_EXHAUSTIVE_MOD103_ZERO_POLE_SHELL_EMBEDDINGS_CANONICAL_MARKING"
assert word["status"] == "PASS_EXACT_Q24_2A5_Q6O1307_LOW_POLE_HORIZONTAL_WORD"
assert scale["status"] == "PASS_EXACT_QQ_P230_SECTION_AND_H0"
assert p1229["status"] == "PASS_EXACT_QQ_P1229_POLYNOMIAL_SECTION"
assert hensel["status"] == "PASS_TARGETED_REGULAR_MOD103_HENSEL_LIFTS"
assert hensel["hensel"]["full_lift_indices"] == [30, 71]

R = PolynomialRing(QQ, "T")
K = R.fraction_field()
T = R.gen()
A = R([QQ(value) for value in surface["child"]["minimal_A_coefficients_low_to_high"]])
B = R([QQ(value) for value in surface["child"]["minimal_B_coefficients_low_to_high"]])
E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])
u = ZZ(scale["method"]["u"])
A_scaled = R([u**4 * value for value in A])
B_scaled = R([u**6 * value for value in B])
p = ZZ(hensel["prime"])
F = GF(p)
RF = PolynomialRing(F, "T")


def reduce_polynomial(polynomial):
    return RF([
        F(value.numerator()) / F(value.denominator())
        for value in R(polynomial).list()
    ])


def recover_polynomial_section(section_index):
    probe = next(
        item for item in hensel["hensel"]["candidate_probes"]
        if item["section_index"] == section_index
    )
    assert probe["full_thirteen_equation_hensel_lift"]
    checkpoint = probe["checkpoints"][-1]
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
    seed_Y = RF(probe["Y_coefficients_low_to_high_mod103"])
    target_Y_scaled = seed_Y * F(u)**3
    if reduce_polynomial(Y_scaled) != target_Y_scaled:
        Y_scaled = -Y_scaled
    assert reduce_polynomial(Y_scaled) == target_Y_scaled
    X = R([value / u**2 for value in X_scaled])
    Y = R([value / u**3 for value in Y_scaled])
    assert (X.degree(), Y.degree()) == (4, 6)
    assert Y**2 == X**3 + A * X + B
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
    return X, Y, X_scaled, Y_scaled, agreement_digits


P1_data = recover_polynomial_section(30)
P32_data = recover_polynomial_section(71)
P1 = E(K(P1_data[0]), K(P1_data[1]))
P32 = E(K(P32_data[0]), K(P32_data[1]))

affine_record = marking["old_A11_affine_section_on_component9_pointed_child"]
assert affine_record["exact_child_identity"]
assert affine_record["x"]["denominator_coefficients_low_to_high"] == ["1"]
assert affine_record["y"]["denominator_coefficients_low_to_high"] == ["1"]
X_affine = R([QQ(value) for value in affine_record["x"]["numerator_coefficients_low_to_high"]])
Y_affine = R([QQ(value) for value in affine_record["y"]["numerator_coefficients_low_to_high"]])
assert Y_affine**2 == X_affine**3 + A * X_affine + B
P_affine = E(K(X_affine), K(Y_affine))

classes = shell["exact_shell"]["indexed_component_nef_classes"]
assert classes[45]["NS_coordinates"] == affine_record["NS_coordinates_in_selected_child_basis"]
target_P146_NS = vector(ZZ, word["sections"]["q4_orbit146"]["effective_section"])
word_P146_NS = sum(
    (vector(ZZ, classes[index]["NS_coordinates"]) for index in (1, 32, 45)),
    vector(ZZ, [0] * 19),
)
P146_trivial_correction = target_P146_NS - word_P146_NS
assert list(P146_trivial_correction) == [
    -1, -2, 1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
]


def projectivize(point, expected_degrees):
    x, y = point.xy()
    denominator_x = R(x.denominator()).monic()
    assert denominator_x.is_square()
    Z = R(denominator_x.sqrt()).monic()
    assert R(y.denominator()).monic() == Z**3
    X = R(x * Z**2)
    Y = R(y * Z**3)
    assert (X.degree(), Y.degree(), Z.degree()) == expected_degrees
    assert Y**2 == X**3 + A * X * Z**4 + B * Z**6
    return X, Y, Z


P146 = P_affine + P1 + P32
X146, Y146, Z146 = projectivize(P146, (6, 9, 1))

X1229 = R([QQ(value) for value in p1229["P1229"]["X_coefficients_low_to_high"]])
Y1229 = R([QQ(value) for value in p1229["P1229"]["Y_coefficients_low_to_high"]])
P1229_point = E(K(X1229), K(Y1229))
P1307 = P146 + P1229_point
X1307, Y1307, Z1307 = projectivize(P1307, (10, 15, 3))

# The direct P146 seed is singular as a lifting chart, but it remains a useful
# independent reduction check on the exact word constructed above.
assert reduce_polynomial(X146) == RF([51, 81, 65, 33, 68, 52, 27])
assert reduce_polynomial(Y146) == RF([39, 13, 22, 37, 30, 94, 45, 102, 63, 22])
assert reduce_polynomial(Z146) == RF([23, 1])


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


def point_payload(X, Y, Z, ns_coordinates, P_dot_O):
    return {
        "X_coefficients_low_to_high": [str(value) for value in X.list()],
        "Y_coefficients_low_to_high": [str(value) for value in Y.list()],
        "Z_coefficients_low_to_high": [str(value) for value in Z.list()],
        "degrees_X_Y_Z": [int(X.degree()), int(Y.degree()), int(Z.degree())],
        "P_dot_O": P_dot_O,
        "NS_coordinates": ns_coordinates,
        "exact_Weierstrass_identity": True,
        "height_profile": height_profile((X, Y, Z)),
    }


inputs = (SURFACE, MARKING, SHELL, WORD, SCALE, P1229, hensel_path)
payload = {
    "schema": "elkies-k3.q24-2a5-p146-p1307-scaled-x-qq.v1",
    "status": "PASS_EXACT_QQ_P146_AND_P1307_SHORT_MW_WORDS",
    "method": {
        "P146_word": "P_affine + P1 + P32",
        "P146_word_exact_indices": [45, 1, 32],
        "P146_target_minus_word_NS_coordinates": [
            int(value) for value in P146_trivial_correction
        ],
        "P1307_word": "P146 + P1229",
        "regular_mod103_polynomial_indices_lifted": [30, 71],
        "direct_singular_P146_projective_lift_avoided": True,
        "large_Groebner_required": False,
    },
    "polynomial_inputs": {
        "P1": point_payload(P1_data[0], P1_data[1], R.one(), classes[1]["NS_coordinates"], 0),
        "P32": point_payload(P32_data[0], P32_data[1], R.one(), classes[32]["NS_coordinates"], 0),
        "old_A11_affine": point_payload(
            X_affine, Y_affine, R.one(), classes[45]["NS_coordinates"], 0
        ),
    },
    "P146": point_payload(
        X146, Y146, Z146,
        word["sections"]["q4_orbit146"]["effective_section"], 1,
    ),
    "P1307": point_payload(
        X1307, Y1307, Z1307,
        word["sections"]["q6_orbit1307"]["effective_section"], 3,
    ),
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
    "P146P1307QQ|P146_bits={}|P1307_bits={}|status={}|output={}".format(
        payload["P146"]["height_profile"]["maximum_rational_bits"],
        payload["P1307"]["height_profile"]["maximum_rational_bits"],
        payload["status"], output,
    ),
    flush=True,
)
