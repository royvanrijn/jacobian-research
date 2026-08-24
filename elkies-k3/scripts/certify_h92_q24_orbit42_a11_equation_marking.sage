#!/usr/bin/env sage -python
"""Bind the exact orbit42 A11 equation to its equation-side lattice frame.

The exact resolved-RR certificate constructs the A11 child but leaves the
identity-shell marking pinned through good reduction.  This script replays
that boundary at p=100003.  It restricts the exact q6 pencil to the eighteen
exact identity-class sections, computes their new-fibre degrees, and compares
the ordered degree fingerprint with every A11 neighbour of the R3-zero D12
frame and every pointed shell isometry.

The selected C10 physical orientation is orbit64 with shell mapping 7.  The
only other match is orbit65 with mapping 6, obtained by exchanging the two
spinor arms.  The output is a marking certificate, not a new equation lift.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, block_diagonal_matrix, identity_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, default=100003)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-a11-equation-marking-orbit64-mod100003.json",
)
args = parser.parse_args()

RR_PATH = LOCAL / "q24-d12-to-a11-orbit42-resolved-rr-qq.json"
PARENT_PATH = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
CANDIDATES_PATH = LOCAL / "q24-orbit42-exact-section-candidates-qq.json"
ZERO_PATH = LOCAL / "q24-orbit42-rational-zero-pole-sections-qq.json"
MODEL_PATH = LOCAL / "q24-orbit42-zero-pole-seeds-mod-43.json"
IDENTITY_PATH = LOCAL / "q24-orbit42-identity-halving-audit.json"
MATCHING_PATH = LOCAL / "q24-orbit42-identity-halving-qq.json"
NEIGHBOURS_PATH = LOCAL / "q24-downstream-lift/d12-c10a-zero-q6-all.json"
FRAME_PATH = LOCAL / "q24-downstream-lift/d12-c10a-zero-frame.txt"
CHILD_FRAME_PATH = (
    LOCAL
    / "q24-downstream-lift/d12-c10a-zero-q6-frames/"
    "q6-o0064-r11-n132-d12-ad4a027cb197.txt"
)

paths = (
    RR_PATH,
    PARENT_PATH,
    CANDIDATES_PATH,
    ZERO_PATH,
    MODEL_PATH,
    IDENTITY_PATH,
    MATCHING_PATH,
    NEIGHBOURS_PATH,
    FRAME_PATH,
    CHILD_FRAME_PATH,
)
for path in paths:
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

rr = json.loads(RR_PATH.read_text())
parent = json.loads(PARENT_PATH.read_text())
candidates = json.loads(CANDIDATES_PATH.read_text())
zero = json.loads(ZERO_PATH.read_text())
model = json.loads(MODEL_PATH.read_text())["exact_model"]
identity = json.loads(IDENTITY_PATH.read_text())
matching = json.loads(MATCHING_PATH.read_text())
neighbours = json.loads(NEIGHBOURS_PATH.read_text())

assert rr["status"] == "PASS_EXACT_Q24_D12_Q6_A11_COMPONENT_VALUATION_RR"
assert candidates["status"] == "PASS_EXACT_Q42_ORBIT42_SECTION_CANDIDATES_QQ"
assert zero["status"] == "PASS_EXACT_Q42_RATIONAL_ZERO_POLE_SECTIONS_QQ"
assert identity["status"] == "PASS_Q42_ORBIT42_IDENTITY_HALVING_LATTICE_GATE"
assert matching["status"] == "Q42_IDENTITY_HALVING_HAS_NO_A11_CHORD"

p = ZZ(args.prime)
F = GF(p)
RU = PolynomialRing(F, "u")
u = RU.gen()
KU = RU.fraction_field()
RV = PolynomialRing(F, "V")
V = RV.gen()
KV = RV.fraction_field()


def reduce_q(value):
    value = QQ(value)
    if value.denominator() % p == 0:
        raise ArithmeticError(f"bad denominator at p={p}")
    return F(value.numerator()) / F(value.denominator())


def reduce_poly_string(value, source_variable, target_ring):
    source_ring = PolynomialRing(QQ, source_variable)
    poly = source_ring(value)
    return target_ring([reduce_q(coefficient) for coefficient in poly.list()])


def reduce_rational_string(value):
    VQ = PolynomialRing(QQ, "V")
    KQ = VQ.fraction_field()
    item = KQ(value)
    return KV(
        RV([reduce_q(coefficient) for coefficient in VQ(item.numerator()).list()])
    ) / KV(
        RV([reduce_q(coefficient) for coefficient in VQ(item.denominator()).list()])
    )


u_of_V = reduce_rational_string(rr["coordinate_change"]["u_of_V"])
x_scale = reduce_rational_string(rr["coordinate_change"]["x_scale"])
y_scale = reduce_rational_string(rr["coordinate_change"]["y_scale"])
A_parent = RV([
    reduce_q(value) for value in parent["child"]["minimal_A_coefficients_low_to_high"]
])
B_parent = RV([
    reduce_q(value) for value in parent["child"]["minimal_B_coefficients_low_to_high"]
])


def evaluate_u(value, argument):
    value = KU(value)
    return KV(RV(value.numerator())(argument)) / KV(RV(value.denominator())(argument))


selected = candidates["candidates"][rr["selected_section"]["candidate_index"]]
X_shell = RU([reduce_q(value) for value in selected["X_coefficients_low_to_high"]])
Y_shell = RU([reduce_q(value) for value in selected["Y_coefficients_low_to_high"]])
Z_shell = RU([reduce_q(value) for value in selected["Z_coefficients_low_to_high"]])
x_parent = x_scale * evaluate_u(KU(X_shell) / KU(Z_shell**2), u_of_V)
y_parent = y_scale * evaluate_u(KU(Y_shell) / KU(Z_shell**3), u_of_V)
E_parent = EllipticCurve(KV, [0, 0, 0, KV(A_parent), KV(B_parent)])
P_parent = E_parent(x_parent, y_parent)

Z_parent_denominator = RV(x_parent.denominator()).monic()
Z_parent = RV.one()
for factor, multiplicity in Z_parent_denominator.factor():
    if int(multiplicity) % 2:
        raise ArithmeticError("selected x denominator is not a square")
    Z_parent *= factor.monic() ** (int(multiplicity) // 2)
X_parent = RV(x_parent * Z_parent**2)
Y_parent = RV(y_parent * Z_parent**3)
assert E_parent(KV(X_parent) / KV(Z_parent**2), KV(Y_parent) / KV(Z_parent**3)) == P_parent

alpha = reduce_q(model["I8star_root"])
collision_modulus = Z_parent**2
X_inverse = X_parent.inverse_mod(collision_modulus)
rr_pairs = []
for BB in (RV.one(), V):
    AA = RV((BB * Y_parent * X_inverse) % collision_modulus)
    AA -= AA(alpha) / Z_parent(alpha) ** 2 * Z_parent**2
    assert AA(alpha) == 0
    assert (AA * X_parent - BB * Y_parent) % collision_modulus == 0
    rr_pairs.append((AA, BB))

(AA0, BB0), (AA1, BB1) = rr_pairs
a0, b0 = KV(AA0) / KV(Z_parent**2), KV(BB0) / KV(Z_parent)
a1, b1 = KV(AA1) / KV(Z_parent**2), KV(BB1) / KV(Z_parent)

equation_degrees = []
for row in zero["sections"]:
    x_shell = RU([reduce_q(value) for value in row["x_coefficients_low_to_high"]])
    y_shell = RU([reduce_q(value) for value in row["y_coefficients_low_to_high"]])
    x_curve = x_scale * evaluate_u(x_shell, u_of_V)
    y_curve = y_scale * evaluate_u(y_shell, u_of_V)
    assert y_curve**2 == x_curve**3 + KV(A_parent) * x_curve + KV(B_parent)
    if x_curve == x_parent:
        raise ArithmeticError("identity-shell curve coincides with the marked orbit42 section")
    chord = (y_curve + y_parent) / (x_curve - x_parent)
    new_base = (a1 + b1 * chord) / (a0 + b0 * chord)
    numerator = RV(new_base.numerator())
    denominator = RV(new_base.denominator())
    common = numerator.gcd(denominator)
    numerator //= common
    denominator //= common
    equation_degrees.append(int(max(numerator.degree(), denominator.degree())))

expected_equation_degrees = [6, 4, 2, 8, 7, 3, 9, 1, 2, 8, 6, 4, 7, 3, 3, 7, 9, 1]
assert equation_degrees == expected_equation_degrees

G0 = matrix(ZZ, [
    [ZZ(value) for value in line.split()]
    for line in FRAME_PATH.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
assert G0.dimensions() == (17, 17)
root = G0[:12, :12]
coupling = G0[:12, 12:]
old_section_classes = []
for values in identity["exact_model_R3_zero"]["identity_vectors"]:
    z = vector(ZZ, values)
    rr_coefficients = -(z * coupling.transpose()) * root.inverse()
    assert all(value in ZZ for value in rr_coefficients)
    old_section_classes.append(
        vector(ZZ, [1, 1] + list(map(ZZ, rr_coefficients)) + list(z))
    )

shell_mappings = matching["matching"]["mappings_abstract_to_equation"]
matches = []
profiles = []
for record in neighbours["neighbors"]:
    if record.get("child_ade") != "A11":
        continue
    adapted_basis = matrix(ZZ, record["child_root_adapted_basis"])
    neighbour_basis = matrix(ZZ, record["neighbor_basis"])
    transition = block_diagonal_matrix(identity_matrix(ZZ, 2), adapted_basis) * neighbour_basis
    assert abs(transition.det()) == 1
    inverse = transition.inverse()
    abstract_degrees = [int((section * inverse)[1]) for section in old_section_classes]
    matching_indices = []
    for mapping_index, mapping in enumerate(shell_mappings):
        ordered = [None] * len(mapping)
        for abstract_index, equation_index in enumerate(mapping):
            ordered[equation_index] = abstract_degrees[abstract_index]
        if ordered == equation_degrees:
            matching_indices.append(mapping_index)
            matches.append({
                "orbit_index": int(record["orbit_index"]),
                "mapping_index": int(mapping_index),
            })
    profiles.append({
        "orbit_index": int(record["orbit_index"]),
        "abstract_identity_shell_degrees": abstract_degrees,
        "matching_indices": matching_indices,
    })

assert matches == [
    {"orbit_index": 64, "mapping_index": 7},
    {"orbit_index": 65, "mapping_index": 6},
]

selected_frame = matrix(ZZ, [
    [ZZ(value) for value in line.split()]
    for line in CHILD_FRAME_PATH.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
orbit64_record = next(
    record for record in neighbours["neighbors"] if record["orbit_index"] == 64
)
assert selected_frame == matrix(ZZ, orbit64_record["child_root_adapted_frame"])

selected_match = matches[0]
payload = {
    "schema": "elkies-k3.h3-q24-a11-equation-marking-orbit64.v1",
    "status": "PASS_Q42_A11_EQUATION_MARKING_ORBIT64_MOD100003",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in paths],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths
        },
    },
    "good_reduction_prime": int(p),
    "equation_identity_shell_new_fibre_degrees": equation_degrees,
    "a11_profiles": profiles,
    "matches": matches,
    "selected": {
        **selected_match,
        "source_frame": str(FRAME_PATH.relative_to(ROOT)),
        "child_frame": str(CHILD_FRAME_PATH.relative_to(ROOT)),
        "physical_orientation": "C10",
    },
    "spinor_conjugate": {
        "orbit_index": 65,
        "mapping_index": 6,
        "physical_orientation": "C11",
    },
    "proof_boundary": (
        "The exact QQ pencil and exact identity-shell sections are reduced at the pinned "
        "good prime 100003. Coprime numerator/denominator degrees at good reduction certify "
        "the characteristic-zero new-fibre degree fingerprint. It selects orbit64/mapping7 "
        "in the C10 orientation, with orbit65/mapping6 as the spinor conjugate. This binds "
        "the exact A11 equation to a lattice frame; it does not construct the next q8 pencil."
    ),
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q42A11MARK|prime={}|degrees={}|matches=64:7,65:6|selected=64:7|status={}".format(
        p,
        ",".join(map(str, equation_degrees)),
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output}", flush=True)
