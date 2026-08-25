#!/usr/bin/env sage -python
"""Certify the D13 q4/o11 zero loop and compose it to pinned R17.

status: ACTIVE_PROOF
claim: exact D13-to-D12 cost-improving splice and full pinned endpoint
inputs: exhaustive loop ranking, canonical equation-D13/D12/pinned transport,
        measured direct q24 RR artifact, and promoted A11 suffix certificate
outputs: artifacts/generated-results/elkies-k3-h3-d13-q4o11-promoted-route-certificate.json
"""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
SEARCH = GENERATED / "elkies-k3-h3-d13-zero-changing-d12-presentations.json"
D13_PATH = LOCAL / "q24-equation-d13-to-pinned-r17.json"
DIRECT_RR = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
A11_PROMOTED = GENERATED / "elkies-k3-h3-a5a5-q6o1307-promoted-route-certificate.json"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = GENERATED / "elkies-k3-h3-d13-q4o11-promoted-route-certificate.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(ZZ, [[ZZ(x) for x in line.split()] for line in path.read_text().splitlines()
                       if line.strip() and not line.lstrip().startswith("#")])


def rows(value):
    return [[int(x) for x in row] for row in value.rows()]


search = json.loads(SEARCH.read_text())
d13_path = json.loads(D13_PATH.read_text())
direct_rr = json.loads(DIRECT_RR.read_text())
a11_promoted = json.loads(A11_PROMOTED.read_text())
assert search["status"] == "PASS_EXACT_D13_ZERO_CHANGING_D12_PRESENTATION_SEARCH"
assert d13_path["status"] == "PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH"
assert direct_rr["status"] == "PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR"
assert a11_promoted["status"] == "PASS_EXACT_PROMOTED_EQUATION_COST_ROUTE_TO_PINNED_R17"

best = search["ranked_presentations"][0]
assert best["first_edge_candidate_id"] == {"q": 4, "old_fibre_degree": 2, "orbit_index": 11}
assert best["explicit_zero_curve"] == "old_D13_component_5"
assert best["q_sequence"] == [4, 4, 24] and best["old_fibre_degrees"] == [2, 2, 2]
assert best["first_edge_child"]["root_data"] == [14, 184, 16]
assert best["explicit_child_root_data"] == [14, 184, 16]
assert best["exit_child_root_data"] == [12, 264, 4]
for prefix in ("first_edge", "return", "exit"):
    assert best[f"{prefix}_nef_audit"]["nef"]
    assert best[f"{prefix}_exact_negative_horizontal_walls"] == []

source = matrix(ZZ, d13_path["equation_d13_frame"])
explicit = matrix(ZZ, best["explicit_child_frame"])
returned = matrix(ZZ, best["returned_frame"])
exit_child = matrix(ZZ, best["exit_child_frame"])
g_source = block_diagonal_matrix(U2, -source)
g_explicit = block_diagonal_matrix(U2, -explicit)
g_returned = block_diagonal_matrix(U2, -returned)
g_exit = block_diagonal_matrix(U2, -exit_child)
s_to_e = matrix(ZZ, best["source_to_explicit_child_basis"])
e_to_r = matrix(ZZ, best["explicit_child_to_returned_D13_basis"])
r_to_d12 = matrix(ZZ, best["returned_D13_to_exit_D12_basis"])
assert abs(s_to_e.det()) == abs(e_to_r.det()) == abs(r_to_d12.det()) == 1
assert s_to_e * g_source * s_to_e.transpose() == g_explicit
assert e_to_r * g_explicit * e_to_r.transpose() == g_returned
assert r_to_d12 * g_returned * r_to_d12.transpose() == g_exit
cumulative_return = e_to_r * s_to_e
assert cumulative_return == matrix(ZZ, best["source_to_returned_D13_basis"])
assert abs(cumulative_return.det()) == 1

# Identify the landing with the exact stored current D12 basis, not just ADE.
canonical_d12_in_source = matrix(ZZ, d13_path["steps"][0]["transition"])
source_to_returned_inverse = cumulative_return.inverse().change_ring(ZZ)
exit_inverse = r_to_d12.inverse().change_ring(ZZ)
canonical_d12_in_exit = canonical_d12_in_source * source_to_returned_inverse * exit_inverse
canonical_d12_frame = matrix(ZZ, d13_path["q24"]["child_frame"])
g_canonical_d12 = block_diagonal_matrix(U2, -canonical_d12_frame)
assert abs(canonical_d12_in_exit.det()) == 1
assert canonical_d12_in_exit * g_exit * canonical_d12_in_exit.transpose() == g_canonical_d12
assert vector(ZZ, canonical_d12_in_exit.row(0)) == vector(ZZ, [1, 0] + [0] * 17)

# Independent endpoint identification through the full equation-D13 marking.
pinned_basis_source = matrix(ZZ, d13_path["equation_d13_to_pinned_r17_transition"])
pinned_basis_exit = pinned_basis_source * source_to_returned_inverse * exit_inverse
pinned = load_matrix(PINNED_FRAME)
g_pinned = block_diagonal_matrix(U2, -pinned)
assert abs(pinned_basis_exit.det()) == 1
assert pinned_basis_exit * g_exit * pinned_basis_exit.transpose() == g_pinned

raw_direct_score = int(search["direct_q24"]["score"])
raw_new_score = int(best["total_equation_cost_score"])
direct_estimated_rr = int(search["direct_q24"]["profile"]["expected_RR_ambient"])
direct_measured_rr = int(direct_rr["rr"]["ambient_dimension"])
assert (raw_direct_score, raw_new_score, direct_estimated_rr, direct_measured_rr) == (28485, 25323, 61, 56)
calibrated_direct_score = raw_direct_score - 120 * (direct_estimated_rr - direct_measured_rr)
calibrated_improvement = calibrated_direct_score - raw_new_score
calibrated_percent = QQ(100 * calibrated_improvement) / calibrated_direct_score
assert calibrated_direct_score == 27885 and calibrated_improvement == 2562

old_combined = calibrated_direct_score + int(a11_promoted["new_splice"]["direct_q104_comparator_score"])
new_combined = raw_new_score + int(a11_promoted["new_splice"]["equation_cost_score"])
combined_improvement = old_combined - new_combined
combined_percent = QQ(100 * combined_improvement) / old_combined
assert (old_combined, new_combined, combined_improvement) == (41403, 35657, 5746)

inputs = (SEARCH, D13_PATH, DIRECT_RR, A11_PROMOTED, PINNED_FRAME)
payload = {
    "schema": "elkies-k3.h3-d13-q4o11-promoted-route-certificate.v1",
    "status": "PASS_EXACT_PROMOTED_D13_EQUATION_COST_ROUTE_TO_PINNED_R17",
    "promotion": {
        "promote_as_lifting_target": True,
        "switch_after": "the existing characteristic-zero D13/MW4 equation frontier",
        "replacement": "replace the direct D13-to-D12 presentation by q4/o11, q4 return, q24 exit",
        "resume_at": "the fully identified current D12/MW5 stage; retain the q6/o1307 A11 suffix promotion",
    },
    "new_D13_splice": {
        "nodes": ["equation D13/MW4", "D5+D9/MW3 q4 orbit11 with old_D13_component_5 zero",
                  "equation D13/MW4 with changed zero", "current D12/MW5"],
        "q_sequence": [4, 4, 24],
        "old_fibre_degrees": [2, 2, 2],
        "profiles": {name: best[f"{name}_profile"] for name in ("first_edge", "return", "exit")},
        "raw_equation_cost_score": raw_new_score,
        "direct_q24_raw_comparator_score": raw_direct_score,
        "direct_q24_measured_RR_ambient": direct_measured_rr,
        "direct_q24_calibrated_comparator_score": calibrated_direct_score,
        "calibrated_score_improvement": calibrated_improvement,
        "calibrated_relative_improvement_percent_exact": str(calibrated_percent),
        "calibrated_relative_improvement_percent_decimal": "{:.6f}".format(float(calibrated_percent)),
        "all_new_edges_exact_nef": True,
    },
    "combined_bottleneck_comparison": {
        "old_direct_q24_plus_direct_q104": old_combined,
        "new_D13_splice_plus_A11_q6o1307_splice": new_combined,
        "absolute_improvement": combined_improvement,
        "relative_improvement_percent_exact": str(combined_percent),
        "relative_improvement_percent_decimal": "{:.6f}".format(float(combined_percent)),
    },
    "landing_current_D12_identification": {
        "canonical_current_D12_basis_in_landing": rows(canonical_d12_in_exit),
        "landing_basis_in_canonical_current_D12": rows(canonical_d12_in_exit.inverse().change_ring(ZZ)),
        "forward_determinant": int(canonical_d12_in_exit.det()),
        "inverse_determinant": int(canonical_d12_in_exit.inverse().det()),
        "fibre_exactly_aligned": True,
        "gram_exactly_aligned": True,
        "root_data": [12, 264, 4],
        "MW_rank": 5,
    },
    "full_route_q_sequence_from_D13": [4, 4, 24, 6] + a11_promoted["full_route_q_sequence_from_A11"],
    "full_route_q_sequence_from_H3": [6, 8, 4, 4, 24, 6] + a11_promoted["full_route_q_sequence_from_A11"],
    "endpoint": {
        "name": "pinned_R17", "root_data": [0, 0, 1], "MW_rank": 17,
        "canonical_pinned_basis_in_landing": rows(pinned_basis_exit),
        "landing_basis_in_canonical_pinned": rows(pinned_basis_exit.inverse().change_ring(ZZ)),
        "forward_determinant": int(pinned_basis_exit.det()),
        "inverse_determinant": int(pinned_basis_exit.inverse().det()),
        "gram_identification": "U plus negative pinned rank17_gram.txt exactly",
    },
    "proof_boundary": (
        "All three new fibres have exact marked U, primitive isotropic classes, complete component, affine, "
        "all-section and finite horizontal-wall gates, exact roots, and bidirectional unimodular transports. "
        "The D12 landing and pinned endpoint use full determinant-one bases. The direct comparator is calibrated "
        "to the measured resolved-RR ambient dimension 56. The new edge dimensions and total compiler scores "
        "remain planning estimates until their characteristic-zero equation lifts are executed."
    ),
    "inputs": {"paths": [str(path.relative_to(ROOT)) for path in inputs],
               "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                           for path in inputs}},
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("PROMOTEDD13|splice_q=4,4,24|new={}|direct_calibrated={}|saving={}|percent={}|"
      "landing_det={}|endpoint_det={}|status={}|output={}".format(
          raw_new_score, calibrated_direct_score, calibrated_improvement, calibrated_percent,
          canonical_d12_in_exit.det(), pinned_basis_exit.det(), payload["status"], OUTPUT
      ), flush=True)
