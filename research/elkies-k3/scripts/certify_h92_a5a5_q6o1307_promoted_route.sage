#!/usr/bin/env sage -python
"""Compose the strict-cost q6/o1307 zero loop with the pinned suffix.

status: ACTIVE_PROOF
claim: full marked-NS q6/o1307 splice and exact pinned-R17 endpoint
inputs: certified splice edges, current-3A3 marking, pinned suffix and Gram
outputs: artifacts/generated-results/elkies-k3-h3-a5a5-q6o1307-promoted-route-certificate.json
"""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
A11_EDGE = GENERATED / "elkies-k3-h3-a11-q8-orbit12-lattice-certificate.json"
FIRST = GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q6-orbit1307-lattice-certificate.json"
ZERO = GENERATED / "elkies-k3-h3-a5a5-q6-o1307-explicit-zero-frames.json"
START = GENERATED / "elkies-k3-h3-a5a5-q6o1307-suffix-marking.json"
RETURN = GENERATED / "elkies-k3-h3-a5a5-q6o1307-q4-return-a5a5-certificate.json"
EXIT = GENERATED / "elkies-k3-h3-a5a5-q6o1307-loop-current-3a3-certificate.json"
LOOPS = GENERATED / "elkies-k3-h3-a5a5-zero-changing-loop-search.json"
DIRECT = GENERATED / "elkies-k3-h3-a5a5-current-route-equation-cost-audit.json"
SUFFIX = GENERATED / "elkies-k3-h3-pinned-r17-current-suffix-marking.json"
CURRENT_3A3 = GENERATED / "elkies-k3-h3-current_3A3-marked-frame.json"
FINGERPRINT = LOCAL / "q24-a11-q8-construction-fingerprint.json"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
OUTPUT = GENERATED / "elkies-k3-h3-a5a5-q6o1307-promoted-route-certificate.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


data = {path: json.loads(path.read_text()) for path in (
    A11_EDGE, FIRST, ZERO, START, RETURN, EXIT, LOOPS, DIRECT, SUFFIX, CURRENT_3A3, FINGERPRINT, MANIFEST
)}
a11 = data[A11_EDGE]
first = data[FIRST]
zero = data[ZERO]
start = data[START]
return_edge = data[RETURN]
exit_edge = data[EXIT]
loops = data[LOOPS]
direct = data[DIRECT]
suffix = data[SUFFIX]
current_3a3 = data[CURRENT_3A3]
fingerprint = data[FINGERPRINT]
manifest = data[MANIFEST]

assert a11["status"] == "PASS_EXACT_A11_Q8_EQUATION_COST_LATTICE_CERTIFICATE"
assert first["status"] == "PASS_EXACT_A5A5_EXPLICIT_ZERO_CANDIDATE_LATTICE_CERTIFICATE"
assert first["edge"]["exact_horizontal_nef_gate"] and not first["edge"]["exact_negative_horizontal_walls"]
assert zero["status"] == "PASS_EXACT_CANDIDATE_A11_COMPONENT_EXPLICIT_ZERO_FRAMES"
assert zero["selected_component_index"] == 10
assert start["status"] == "PASS_EXACT_A5A5_CANDIDATE_SUFFIX_MARKING"
for edge in (return_edge, exit_edge):
    assert edge["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
    assert edge["first_edge_nef_audit"]["nef_in_selected_component_chamber"]
    assert edge["first_edge_exact_horizontal_nef_gate"]
    assert not edge["first_edge_exact_negative_horizontal_walls"]
assert suffix["status"] == "PASS_EXACT_PINNED_R17_CURRENT_SUFFIX_MARKING"
assert current_3a3["status"] == "PASS_EXACT_CURRENT_SUFFIX_STAGE_MARKING"

start_frame_path = ROOT / start["frame_output"]
return_frame_path = ROOT / return_edge["frame_output"]
landing_frame_path = ROOT / exit_edge["frame_output"]
start_frame = load_matrix(start_frame_path)
return_frame = load_matrix(return_frame_path)
landing_frame = load_matrix(landing_frame_path)
g_start = block_diagonal_matrix(U2, -start_frame)
g_return = block_diagonal_matrix(U2, -return_frame)
g_landing = block_diagonal_matrix(U2, -landing_frame)
t_return = matrix(ZZ, return_edge["source_to_child_basis"])
t_exit = matrix(ZZ, exit_edge["source_to_child_basis"])
assert t_return * g_start * t_return.transpose() == g_return
assert t_exit * g_return * t_exit.transpose() == g_landing
cumulative = t_exit * t_return
inverse = cumulative.inverse().change_ring(ZZ)
assert abs(cumulative.det()) == 1

# Full landing identification with the certified current 3A3 basis.
current_3a3_basis_start = matrix(
    ZZ, start["current_suffix_stage_bases_in_root_adapted_hub"]["current_3A3"]
)
current_3a3_basis_landing = current_3a3_basis_start * inverse
canonical_3a3_frame = load_matrix(ROOT / current_3a3["frame_output"])
g_canonical_3a3 = block_diagonal_matrix(U2, -canonical_3a3_frame)
assert abs(current_3a3_basis_landing.det()) == 1
assert current_3a3_basis_landing * g_landing * current_3a3_basis_landing.transpose() == g_canonical_3a3
assert vector(ZZ, current_3a3_basis_landing.row(0)) == vector(ZZ, [1, 0] + [0] * 17)

# Independent full endpoint identification with the repository-pinned R17 basis.
historical_in_equation = block_diagonal_matrix(
    identity_matrix(ZZ, 2),
    matrix(ZZ, fingerprint["selected"]["frame_isometry_historical_basis_in_equation_coordinates"]),
)
pinned_in_historical_a11 = matrix(
    ZZ, suffix["current_suffix_stages"]["current_A11"]["pinned_R17_basis_in_stage"]
)
pinned_basis_equation = pinned_in_historical_a11 * historical_in_equation
equation_to_start = matrix(ZZ, start["equation_A11_to_root_adapted_hub_basis"])
pinned_basis_start = pinned_basis_equation * equation_to_start.inverse().change_ring(ZZ)
pinned_basis_landing = pinned_basis_start * inverse
pinned = load_matrix(PINNED_FRAME)
g_pinned = block_diagonal_matrix(U2, -pinned)
assert abs(pinned_basis_landing.det()) == 1
assert pinned_basis_landing * g_landing * pinned_basis_landing.transpose() == g_pinned

best = loops["ranked_loops"][0]
assert best["first_edge_candidate_id"] == {"q": 6, "old_fibre_degree": 2, "orbit_index": 1307}
assert best["explicit_zero_curve"] == "old_A11_component_10"
assert best["loop_q_sequence"] == [6, 4, 6]
assert best["loop_old_fibre_degrees"] == [2, 2, 2]
loop_score = int(best["loop_equation_cost_score"])
direct_score = int(loops["direct_q104_comparator_score"])
assert loop_score == 10334 and direct_score == 13518 and loop_score < direct_score
improvement = direct_score - loop_score
improvement_percent = QQ(100 * improvement) / direct_score
negative_direct = [
    name for name, value in direct["direct_equation_cost_profile"]["named_explicit_curve_degrees"].items()
    if value < 0
]
assert len(negative_direct) == 2

current_suffix_after_3a3 = [
    {"parent": step["parent"], "child": step["child"], "q": int(step["q"]), "orbit": int(step["orbit"])}
    for index, step in enumerate(manifest["forward_steps"])
    if index >= 4
]

inputs = (
    A11_EDGE, FIRST, ZERO, START, RETURN, EXIT, LOOPS, DIRECT, SUFFIX,
    CURRENT_3A3, FINGERPRINT, MANIFEST, PINNED_FRAME, start_frame_path, return_frame_path, landing_frame_path,
)
payload = {
    "schema": "elkies-k3.h3-a5a5-q6o1307-promoted-route-certificate.v1",
    "status": "PASS_EXACT_PROMOTED_EQUATION_COST_ROUTE_TO_PINNED_R17",
    "promotion": {
        "promote_as_lifting_target": True,
        "switch_after": "A11 q8 orbit12 reaches the equation-explicit 2A5/MW7 zero",
        "replacement": "replace the direct q104 2A5-to-3A3 presentation by q6/o1307, q4 return, q6 exit",
        "resume_at": "the fully identified current_3A3 stage, then use the existing certified suffix",
    },
    "new_splice": {
        "nodes": [
            "current 2A5/MW7 explicit orbit12 zero",
            "A1+A3+A5/MW8 q6 orbit1307 with old_A11_component_10 zero",
            "current 2A5/MW7 with changed zero",
            "current 3A3/MW8",
        ],
        "q_sequence": [6, 4, 6],
        "old_fibre_degrees": [2, 2, 2],
        "edge_certificates": [
            str(FIRST.relative_to(ROOT)), str(RETURN.relative_to(ROOT)), str(EXIT.relative_to(ROOT)),
        ],
        "equation_cost_score": loop_score,
        "direct_q104_comparator_score": direct_score,
        "absolute_score_improvement": improvement,
        "relative_score_improvement_percent_exact": str(improvement_percent),
        "relative_score_improvement_percent_decimal": "{:.6f}".format(float(improvement_percent)),
        "strict_score_improvement": True,
        "direct_q104_negative_named_explicit_curves": negative_direct,
        "all_new_edges_exact_nef": True,
    },
    "landing_current_3A3_identification": {
        "canonical_current_3A3_basis_in_landing": rows(current_3a3_basis_landing),
        "landing_basis_in_canonical_current_3A3": rows(current_3a3_basis_landing.inverse().change_ring(ZZ)),
        "forward_determinant": int(current_3a3_basis_landing.det()),
        "inverse_determinant": int(current_3a3_basis_landing.inverse().det()),
        "fibre_exactly_aligned": True,
        "gram_exactly_aligned": True,
        "root_data": exit_edge["child"]["root_data"],
        "MW_rank": exit_edge["child"]["mw_rank"],
    },
    "existing_suffix_after_landing": current_suffix_after_3a3,
    "full_route_q_sequence_from_A11": [8, 6, 4, 6, 4, 4, 4, 4, 4, 4, 6],
    "endpoint": {
        "name": "pinned_R17",
        "root_data": [0, 0, 1],
        "MW_rank": 17,
        "canonical_pinned_basis_in_landing": rows(pinned_basis_landing),
        "landing_basis_in_canonical_pinned": rows(pinned_basis_landing.inverse().change_ring(ZZ)),
        "forward_determinant": int(pinned_basis_landing.det()),
        "inverse_determinant": int(pinned_basis_landing.inverse().det()),
        "gram_identification": "U plus negative pinned rank17_gram.txt exactly",
    },
    "proof_boundary": (
        "The new three-edge splice has exact primitive isotropic fibres, complete component, "
        "affine, all-section and finite horizontal-wall nef gates, marked U, exact root/MW "
        "data, and determinant-one transports. Its landing is identified by a full basis, "
        "not ADE/MW alone, with the existing current 3A3 stage. The unchanged exact suffix "
        "then identifies the endpoint with pinned R17. Compiler scores are deterministic "
        "planning estimates; the characteristic-zero equation lifts remain to be executed."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("PROMOTED1307|splice_q=6,4,6|score={}|direct={}|improvement={}|percent={}|landing_det={}|endpoint_det={}|status={}|output={}".format(
    loop_score, direct_score, improvement, improvement_percent,
    current_3a3_basis_landing.det(), pinned_basis_landing.det(), payload["status"], OUTPUT,
), flush=True)
