#!/usr/bin/env sage -python
"""Compose the physical q10 2A5 exit with the pinned R17 suffix."""

import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
Q10 = LOCAL / "q24-2a5-direct-physical-q10-certificate.json"
EFFECTIVE_ZERO = LOCAL / "q24-2a5-direct-physical-q10-effective-c5-zero-certificate.json"
Q8_MARKING = LOCAL / "q24-a11-to-2a5-q8-equation-marking-qq.json"
ZERO_FRAME = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
CURRENT_3A3 = GENERATED / "elkies-k3-h3-current_3A3-marked-frame.json"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = GENERATED / "elkies-k3-h3-a5a5-direct-physical-q10-promoted-route-certificate.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(item) for item in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


q10 = json.loads(Q10.read_text())
effective_zero = json.loads(EFFECTIVE_ZERO.read_text())
q8_marking = json.loads(Q8_MARKING.read_text())
zero_frame = json.loads(ZERO_FRAME.read_text())
current_3a3 = json.loads(CURRENT_3A3.read_text())
manifest = json.loads(MANIFEST.read_text())
assert q10["status"] == "PASS_EXACT_PHYSICAL_NEF_Q10_CURRENT_3A3_PRESENTATION"
assert effective_zero["status"] == "PASS_EXACT_PHYSICAL_Q10_EFFECTIVE_C5_ZERO_TO_PINNED_R17"
assert q8_marking["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_EQUATION_MARKING"
assert zero_frame["status"] == "PASS_EXACT_A11_Q8_ORBIT12_EXPLICIT_ZERO_FRAMES"
assert current_3a3["status"] == "PASS_EXACT_CURRENT_SUFFIX_STAGE_MARKING"
assert manifest["status"] == "PASS_H3_R17_BACKWARD_EXACT_LIFT_MANIFEST"

parent_frame = matrix(ZZ, zero_frame["selected"]["frame"])
current_frame_path = ROOT / current_3a3["frame_output"]
current_frame = load_matrix(current_frame_path)
pinned_frame = load_matrix(PINNED_FRAME)
g_parent = block_diagonal_matrix(U2, -parent_frame)
g_current = block_diagonal_matrix(U2, -current_frame)
g_pinned = block_diagonal_matrix(U2, -pinned_frame)

parent_to_child = matrix(ZZ, effective_zero["transport"]["parent_to_effective_zero_child_basis"])
child_to_parent = matrix(ZZ, effective_zero["transport"]["effective_zero_child_to_parent_basis"])
child_frame = matrix(ZZ, effective_zero["selection"]["child_frame"])
g_child = block_diagonal_matrix(U2, -child_frame)
assert parent_to_child * child_to_parent == matrix.identity(ZZ, 19)
assert child_to_parent * parent_to_child == matrix.identity(ZZ, 19)
assert parent_to_child * g_parent * parent_to_child.transpose() == g_child
assert abs(parent_to_child.det()) == 1

current_in_child = matrix(
    ZZ, effective_zero["current_3A3_identification"]["current_3A3_basis_in_effective_zero_child"]
)
child_in_current = matrix(
    ZZ, effective_zero["current_3A3_identification"]["effective_zero_child_basis_in_current_3A3"]
)
assert current_in_child * child_in_current == matrix.identity(ZZ, 19)
assert current_in_child * g_child * current_in_child.transpose() == g_current

pinned_in_current = matrix(ZZ, current_3a3["pinned_R17_basis_in_source"])
assert abs(pinned_in_current.det()) == 1
assert pinned_in_current * g_current * pinned_in_current.transpose() == g_pinned
pinned_in_child = matrix(
    ZZ, effective_zero["endpoint"]["canonical_pinned_basis_in_effective_zero_child"]
)
assert pinned_in_child * g_child * pinned_in_child.transpose() == g_pinned
pinned_in_parent = pinned_in_child * parent_to_child
assert abs(pinned_in_parent.det()) == 1
assert pinned_in_parent * g_parent * pinned_in_parent.transpose() == g_pinned

# Compose through the exact component-9-zero realization of the already-lifted
# A11 q8/orbit12 edge.  Rows of this matrix are the explicit 2A5 basis in the
# equation-A11 coordinates.
equation_to_parent = matrix(
    ZZ, zero_frame["selected"]["equation_A11_to_explicit_zero_basis"]
)
parent_to_equation = equation_to_parent.inverse().change_ring(ZZ)
g_equation = parent_to_equation * g_parent * parent_to_equation.transpose()
current_in_equation = current_in_child * parent_to_child * equation_to_parent
equation_in_current = current_in_equation.inverse().change_ring(ZZ)
assert abs(current_in_equation.det()) == 1
assert current_in_equation * g_equation * current_in_equation.transpose() == g_current
pinned_in_equation = pinned_in_child * parent_to_child * equation_to_parent
equation_in_pinned = pinned_in_equation.inverse().change_ring(ZZ)
assert abs(pinned_in_equation.det()) == 1
assert pinned_in_equation * g_equation * pinned_in_equation.transpose() == g_pinned

profile = q10["RR_profile"]
pairings = q10["physical_weyl_repair"]["known_effective_curve_pairings"]
degree_zero = sum(value == 0 for value in pairings.values())
degree_one = sum(value == 1 for value in pairings.values())
terms = {
    "P_dot_O": 900 * int(profile["P_dot_O"]),
    "horizontal_degree": 250 * int(q10["physical_weyl_repair"]["old_fibre_degree"]),
    "RR_ambient": 120 * int(profile["expected_RR_ambient"]),
    "vertical_layers": 60 * int(profile["vertical_connected_layers"]),
    "vertical_support": 25 * sum(
        value != 0 for value in profile["vertical_coefficients_in_physical_root_adapted_frame"]
    ),
    "child_root_count": int(q10["landing"]["root_data"][1]),
    "coordinate_growth": max(abs(value) for value in q10["physical_weyl_repair"]["repaired_fibre"]),
    "no_explicit_degree_one_curve": 4000 if degree_one == 0 else 0,
    "explicit_degree_one_credit": -500 * min(degree_one, 6),
    "explicit_degree_zero_credit": -100 * min(degree_zero, 12),
}
score = int(sum(terms.values()))
assert degree_zero == 8 and degree_one == 4 and score == 4471

# Retain the older named-curve convention as a conservative comparison.  It
# omitted the two current-I6 affine components even though both are exact
# physical curves in the equation marking.
conservative_pairings = {
    name: value for name, value in pairings.items()
    if name not in {"first_I6_affine_component", "second_I6_affine_component"}
}
conservative_zero = sum(value == 0 for value in conservative_pairings.values())
conservative_one = sum(value == 1 for value in conservative_pairings.values())
conservative_terms = dict(terms)
conservative_terms["no_explicit_degree_one_curve"] = 4000 if conservative_one == 0 else 0
conservative_terms["explicit_degree_one_credit"] = -500 * min(conservative_one, 6)
conservative_terms["explicit_degree_zero_credit"] = -100 * min(conservative_zero, 12)
conservative_score = int(sum(conservative_terms.values()))
assert conservative_score == 5071

withdrawn_direct_score = 13518
improvement = withdrawn_direct_score - score
improvement_percent = QQ(100 * improvement) / withdrawn_direct_score
suffix = [
    {"parent": step["parent"], "child": step["child"], "q": int(step["q"]), "orbit": int(step["orbit"])}
    for step in manifest["forward_steps"][4:]
]

inputs = (Q10, EFFECTIVE_ZERO, Q8_MARKING, ZERO_FRAME, CURRENT_3A3, current_frame_path, MANIFEST, PINNED_FRAME)
payload = {
    "schema": "elkies-k3.h3-a5a5-direct-physical-q10-promoted-route.v1",
    "status": "PASS_EXACT_PROMOTED_PHYSICAL_Q10_ROUTE_TO_PINNED_R17",
    "promotion": {
        "promote_as_lifting_target": True,
        "switch_after": "the exact A11 q8/orbit12 equation reaches component-9-zero 2A5/MW7",
        "next_equation_lift": "physical q10 degree-two pencil with effective old_A11_component_5 zero, P.O=5, and expected RR ambient 15",
        "resume_after_lift": "canonical current_3A3, then the unchanged certified suffix to pinned R17",
        "withdrawn_targets": [
            "q4/orbit230--q6/orbit1315 pseudo-zero route (4199)",
            "q6/orbit1307 component10 route (10334)",
            "stored nonphysical q104 comparator (13518)",
        ],
    },
    "splice": {
        "parent": "equation-explicit component-9-zero 2A5/MW7",
        "child": "canonical current 3A3/MW8",
        "q": 10,
        "old_fibre_degree": 2,
        "physical_reflection_count": int(q10["physical_weyl_repair"]["reflection_count"]),
        "equation_effective_zero": effective_zero["selection"]["zero"],
        "canonical_chamber_zero_rejected": True,
        "P_dot_O": int(profile["P_dot_O"]),
        "expected_RR_ambient": int(profile["expected_RR_ambient"]),
        "vertical_connected_layers": int(profile["vertical_connected_layers"]),
        "explicit_degree_zero_curve_count": degree_zero,
        "explicit_degree_one_curve_count": degree_one,
        "operational_equation_cost_score": score,
        "operational_equation_cost_terms": terms,
        "conservative_old_named_curve_score": conservative_score,
        "conservative_old_named_curve_terms": conservative_terms,
        "withdrawn_q104_score": withdrawn_direct_score,
        "absolute_improvement": improvement,
        "relative_improvement_percent_exact": str(improvement_percent),
        "strict_improvement": True,
        "lattice_certificate": str(Q10.relative_to(ROOT)),
        "effective_zero_certificate": str(EFFECTIVE_ZERO.relative_to(ROOT)),
    },
    "landing_current_3A3_identification": {
        "effective_zero_child_basis_in_parent": rows(parent_to_child),
        "parent_basis_in_effective_zero_child": rows(child_to_parent),
        "current_3A3_basis_in_effective_zero_child": rows(current_in_child),
        "effective_zero_child_basis_in_current_3A3": rows(child_in_current),
        "forward_determinant": int(current_in_child.det()),
        "inverse_determinant": int(child_in_current.det()),
        "gram_exactly_aligned": True,
        "root_data": q10["landing"]["root_data"],
        "MW_rank": q10["landing"]["MW_rank"],
    },
    "existing_suffix_after_landing": suffix,
    "full_route_q_sequence_from_A11": [8, 10] + [item["q"] for item in suffix],
    "endpoint": {
        "name": "pinned_R17",
        "root_data": [0, 0, 1],
        "MW_rank": 17,
        "canonical_pinned_basis_in_equation_A11": rows(pinned_in_equation),
        "equation_A11_basis_in_canonical_pinned": rows(equation_in_pinned),
        "forward_determinant": int(pinned_in_equation.det()),
        "inverse_determinant": int(equation_in_pinned.det()),
        "gram_identification": "U plus negative pinned rank17_gram.txt exactly",
    },
    "proof_boundary": (
        "The q10 edge has an exact primitive nef isotropic fibre, an equation-effective "
        "component-5 zero, known physical-component, "
        "all-section, and complete finite horizontal-wall gates, an exact marked U, exact "
        "3A3 root/MW data, and bidirectional determinant-one NS transport. The full composed "
        "basis identifies the endpoint with pinned R17, not merely by ADE/MW. The score is a "
        "deterministic compiler estimate; the q10 characteristic-zero RR pencil, quartic, "
        "Jacobian, fibres, and equation marking remain to be compiled."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in inputs
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "PROMOTEDQ10|score={}|conservative={}|old={}|improvement={}|landing_det={}|"
    "endpoint_det={}|status={}|output={}".format(
        score, conservative_score, withdrawn_direct_score, improvement,
        current_in_child.det(), pinned_in_equation.det(), payload["status"], OUTPUT,
    ),
    flush=True,
)
