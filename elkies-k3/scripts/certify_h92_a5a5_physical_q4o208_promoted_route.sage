#!/usr/bin/env sage -python
"""Promote the fully marked physical q4/orbit208 route to pinned R17."""

import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
CERT = GENERATED / "elkies-k3-h3-a5a5-physical-q4o208-to-pinned-r17-certificate.json"
CURRENT = GENERATED / "elkies-k3-h3-current_3A3-marked-frame.json"
MANIFEST = ROOT / "artifacts/local/elkies-k3/h3-r17-backward-exact-lift-manifest.json"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
DIRECT_FRONTIER = GENERATED / "elkies-k3-h3-a5a5-physical-component-chamber-d2-q4q6q8q10-3a3-only-frontier.json"
DIRECT_COST = GENERATED / "elkies-k3-h3-a5a5-physical-component-chamber-d2-q4q6q8q10-3a3-only-equation-cost.json"
RR_QQ = ROOT / "artifacts/local/elkies-k3/q24-2a5-physical-q4o208-rr-qq.json"
MARKING_QQ = ROOT / "artifacts/local/elkies-k3/q24-2a5-physical-q4o208-equation-marking-qq.json"
SUFFIX_NEF = GENERATED / "elkies-k3-h3-q4o208-canonical-suffix-physical-nef-audit.json"
OUTPUT = GENERATED / "elkies-k3-h3-a5a5-physical-q4o208-promoted-route-certificate.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


cert = json.loads(CERT.read_text())
current = json.loads(CURRENT.read_text())
manifest = json.loads(MANIFEST.read_text())
direct_frontier = json.loads(DIRECT_FRONTIER.read_text())
direct_cost = json.loads(DIRECT_COST.read_text())
rr_qq = json.loads(RR_QQ.read_text())
marking_qq = json.loads(MARKING_QQ.read_text())
suffix_nef = json.loads(SUFFIX_NEF.read_text())
assert cert["status"] == "PASS_EXACT_PHYSICAL_Q4O208_3A3_TO_PINNED_R17"
assert current["status"] == "PASS_EXACT_CURRENT_SUFFIX_STAGE_MARKING"
assert manifest["status"] == "PASS_H3_R17_BACKWARD_EXACT_LIFT_MANIFEST"
assert direct_frontier["status"] == "PASS_EXACT_MARKED_ROOT_ADAPTED_FRONTIER_RANKING"
assert direct_cost["status"] == "PASS_EXACT_MARKED_FRONTIER_EQUATION_COST_SCORING"
assert rr_qq["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_3A3_RR_AND_JACOBIAN"
assert marking_qq["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING"
assert marking_qq["selected_zero"] == "old_A11_component_5"
assert marking_qq["root_lattice"]["type"] == "3A3"
assert marking_qq["root_lattice"]["determinant"] == 64
assert len(marking_qq["physical_fibres"]) == 3
assert all(
    row["identity_component_index"] in range(4)
    for row in marking_qq["physical_fibres"].values()
)
assert marking_qq["transport"]["inverse_exact"]
assert marking_qq["transport"]["Gram_transport_exact"]
assert suffix_nef["status"] == "PASS_EXACT_Q4O208_CANONICAL_SUFFIX_PHYSICAL_WALL_CORRECTION"
assert not suffix_nef["canonical_pullback_to_C5_equation_frame"]["nef_in_physical_equation_chamber"]
assert suffix_nef["wall_correction"]["number_of_reflections"] == 1
assert rr_qq["divisor"] == {
    "P_dot_O": 0,
    "horizontal": "P1229",
    "identity": "F_q4=O+P1229+C10+C8",
    "old_fibre_degree": 2,
    "q": 4,
}
assert rr_qq["resolved_RR"]["kernel_dimension"] == 2
assert rr_qq["resolved_RR"]["condition_rank"] == 2
assert len(rr_qq["resolved_RR"]["ambient_basis"]) == 4
assert rr_qq["quartic"]["degree"] == 4
assert rr_qq["child"]["degrees_A_B_Delta"] == [8, 12, 24]
assert rr_qq["child"]["finite_fibres"] == [
    {"count": 3, "kodaira": "I4"},
    {"count": 12, "kodaira": "I1"},
]
assert rr_qq["child"]["infinity"]["kodaira"] == "smooth"
assert rr_qq["child"]["euler_number"] == 24
assert rr_qq["child"]["ADE"] == "3A3"
assert direct_frontier["input_candidate_count"] == 181
assert direct_frontier["full_nef_candidate_count"] == 56
assert direct_cost["best_candidate"]["candidate_id"] == {
    "q": 4, "old_fibre_degree": 2, "orbit_index": 208,
}
assert direct_cost["best_candidate"]["equation_cost_score"] == -1412
assert cert["compiler_profile"]["operational_equation_cost_score"] == -1412
assert cert["compiler_profile"]["expected_RR_ambient"] == 4
assert cert["compiler_profile"]["P_dot_O"] == 0
assert cert["selection"]["zero"] == "old_A11_component_5"
assert cert["selection"]["child_root_data"] == [9, 36, 64]

source_frame = load_matrix(ROOT / "artifacts/generated-results/elkies-k3-h3-a5a5-physical-component-chamber-frame.txt")
source_gram = block_diagonal_matrix(U2, -source_frame)
pinned_frame = load_matrix(PINNED_FRAME)
pinned_gram = block_diagonal_matrix(U2, -pinned_frame)
pinned_in_source = matrix(ZZ, cert["endpoint"]["pinned_R17_basis_in_parent"])
source_in_pinned = matrix(ZZ, cert["endpoint"]["parent_basis_in_pinned_R17"])
assert pinned_in_source * source_in_pinned == matrix.identity(ZZ, 19)
assert pinned_in_source * source_gram * pinned_in_source.transpose() == pinned_gram
assert abs(pinned_in_source.det()) == 1

suffix = [
    {"parent": step["parent"], "child": step["child"], "q": int(step["q"]), "orbit": int(step["orbit"])}
    for step in manifest["forward_steps"][4:]
]
old_score = 4471
new_score = int(cert["compiler_profile"]["operational_equation_cost_score"])
gross_score = sum(
    value for name, value in cert["compiler_profile"]["operational_equation_cost_terms"].items()
    if not name.endswith("_credit") and not name.endswith("_penalty")
)
assert gross_score == 1388

inputs = (
    CERT, CURRENT, MANIFEST, PINNED_FRAME, DIRECT_FRONTIER, DIRECT_COST,
    RR_QQ, MARKING_QQ, SUFFIX_NEF,
)
payload = {
    "schema": "elkies-k3.h3-a5a5-physical-q4o208-promoted-route.v1",
    "status": "PASS_EXACT_PROMOTED_PHYSICAL_Q4O208_ROUTE_TO_PINNED_R17",
    "promotion": {
        "promote_as_lifting_target": True,
        "switch_after": "the exact A11 q8/orbit12 equation reaches component-9-zero 2A5/MW7",
        "next_equation_lift": (
            "physical q4/orbit208 degree-two pencil from the literal I4 divisor "
            "old_zero + P1229 + C10 + C8; use effective C5 as the child zero"
        ),
        "resume_after_lift": (
            "canonical current 3A3 lattice landing; do not compile the raw inherited suffix "
            "until its physical chamber correction and effective zero pass"
        ),
        "supersedes_lifting_target": "physical q10/RR15 route with score 4471",
        "direct_3A3_q4_q6_q8_q10_closure": {
            "exact_candidates": 181, "full_nef_candidates": 56,
            "unique_cost_leader": "q4/orbit208", "leader_score": -1412,
        },
        "equation_lift": {
            "status": rr_qq["status"],
            "resolved_RR_dimensions": [4, 2, 2],
            "quartic_degree": 4,
            "jacobian_degrees_A_B_Delta": [8, 12, 24],
            "fibres": "3I4 + 12I1; infinity smooth",
            "equation_marking_status": marking_qq["status"],
            "selected_zero": marking_qq["selected_zero"],
            "opposite_exact_section": "first_I6_affine_component",
            "physical_I4_supports_mod_103": sorted(
                row["support_mod_103"] for row in marking_qq["physical_fibres"].values()
            ),
            "next_gate": (
                "certify the one-wall-corrected 3A3-to-A3+2A2 class against all sections and "
                "select an equation-effective zero, or pivot to a cheaper physical exit"
            ),
        },
    },
    "splice": {
        "parent": cert["parent"], "child": "canonical current 3A3/MW8",
        "q": 4, "orbit": 208, "old_fibre_degree": 2,
        "equation_effective_zero": cert["selection"]["zero"],
        "P_dot_O": 0, "expected_RR_ambient": 4, "vertical_connected_layers": 2,
        "literal_special_I4": cert["compiler_profile"]["literal_special_I4"],
        "compiled_resolved_RR": {
            "ambient_dimension": 4,
            "condition_rank": 2,
            "h0": 2,
            "maximum_kernel_rational_bits": rr_qq["resolved_RR"]["maximum_kernel_rational_bits"],
        },
        "compiled_jacobian": {
            "degrees_A_B_Delta": rr_qq["child"]["degrees_A_B_Delta"],
            "finite_fibres": rr_qq["child"]["finite_fibres"],
            "infinity": rr_qq["child"]["infinity"],
            "euler_number": rr_qq["child"]["euler_number"],
        },
        "explicit_degree_zero_curve_count": len(cert["compiler_profile"]["explicit_degree_zero_curves"]),
        "explicit_degree_one_curve_count": len(cert["compiler_profile"]["explicit_degree_one_curves"]),
        "operational_equation_cost_score": new_score,
        "gross_positive_compiler_burden": gross_score,
        "superseded_q10_score": old_score,
        "absolute_score_reduction": old_score - new_score,
        "strictly_below_withdrawn_4199_threshold": new_score < 4199,
        "lattice_certificate": str(CERT.relative_to(ROOT)),
    },
    "existing_suffix_after_landing": suffix,
    "physical_suffix_gate": {
        "raw_canonical_first_edge_nef": False,
        "negative_explicit_wall": "old_zero=old_A11_component_9",
        "known_wall_correction_reflections": 1,
        "known_curve_nonnegative_after_correction": True,
        "all_section_nef_and_effective_zero_still_required": True,
        "audit": str(SUFFIX_NEF.relative_to(ROOT)),
    },
    "full_route_q_sequence_from_A11": [8, 4] + [item["q"] for item in suffix],
    "endpoint": {
        "name": "pinned_R17", "root_data": [0, 0, 1], "MW_rank": 17,
        "forward_determinant": int(pinned_in_source.det()),
        "inverse_determinant": int(source_in_pinned.det()),
        "gram_identification": "U plus negative pinned rank17_gram.txt exactly",
    },
    "proof_boundary": (
        "This promotion inherits the exact physical nef and horizontal-wall gates, literal "
        "four-curve I4 divisor, equation-effective C5 marked U, complete root data, bidirectional "
        "unimodular transports, and pinned endpoint from the dedicated certificate. The score is "
        "a compiler-planning estimate. The exact characteristic-zero H0 plane, quartic, "
        "Jacobian, fibre classification, equation-effective C5 pointing, opposite affine "
        "section, all three physical I4 cycles, and the full old-curve-to-child equation "
        "marking now pass. The raw canonical first suffix fibre crosses one explicit effective "
        "wall in this physical equation chamber; its one-reflection correction is nonnegative "
        "on all inherited explicit curves, but all-section nefness and an effective zero remain "
        "separate gates before any suffix equation compilation."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "PROMOTEDQ4O208|score={}|gross={}|old_q10={}|reduction={}|RR=4|PO=0|"
    "endpoint_det={}|status={}|output={}".format(
        new_score, gross_score, old_score, old_score - new_score,
        pinned_in_source.det(), payload["status"], OUTPUT,
    ), flush=True,
)
