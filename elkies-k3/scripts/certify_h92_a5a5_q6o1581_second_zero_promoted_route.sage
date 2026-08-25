#!/usr/bin/env sage -python
"""Certify the q6/o1581 second zero change and compose it to pinned R17.

status: ACTIVE_PROOF
claim: strict second-zero improvement of the promoted A11 suffix
outputs: artifacts/generated-results/elkies-k3-h3-a5a5-q6o1581-second-zero-promoted-route-certificate.json
"""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
SEARCH = GENERATED / "elkies-k3-h3-a5a5-q6o1307-second-zero-changing-3a3-presentations.json"
SOURCE_MARKING = GENERATED / "elkies-k3-h3-a5a5-q6o1307-q4-return-a5a5-certificate.json"
DIRECT_EXIT = GENERATED / "elkies-k3-h3-a5a5-q6o1307-loop-current-3a3-certificate.json"
OLD_PROMOTED = GENERATED / "elkies-k3-h3-a5a5-q6o1307-promoted-route-certificate.json"
D13_PROMOTED = GENERATED / "elkies-k3-h3-d13-q4o11-promoted-route-certificate.json"
CURRENT_3A3 = GENERATED / "elkies-k3-h3-current_3A3-marked-frame.json"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = GENERATED / "elkies-k3-h3-a5a5-q6o1581-second-zero-promoted-route-certificate.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(ZZ, [[ZZ(x) for x in line.split()] for line in path.read_text().splitlines()
                       if line.strip() and not line.lstrip().startswith("#")])


def rows(value):
    return [[int(x) for x in row] for row in value.rows()]


search = json.loads(SEARCH.read_text())
source_marking = json.loads(SOURCE_MARKING.read_text())
direct_exit = json.loads(DIRECT_EXIT.read_text())
old_promoted = json.loads(OLD_PROMOTED.read_text())
d13_promoted = json.loads(D13_PROMOTED.read_text())
current_3a3 = json.loads(CURRENT_3A3.read_text())
assert search["status"] == "PASS_EXACT_A5A5_ZERO_CHANGING_3A3_PRESENTATION_SEARCH"
assert source_marking["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
assert direct_exit["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
assert old_promoted["status"] == "PASS_EXACT_PROMOTED_EQUATION_COST_ROUTE_TO_PINNED_R17"
assert d13_promoted["status"] == "PASS_EXACT_PROMOTED_D13_EQUATION_COST_ROUTE_TO_PINNED_R17"

best = search["ranked_presentations"][0]
assert best["first_edge_candidate_id"] == {"q": 6, "old_fibre_degree": 2, "orbit_index": 1581}
assert best["explicit_zero_curve"] == "old_A5A5_component_0"
assert best["q_sequence"] == [6, 4, 8] and best["old_fibre_degrees"] == [2, 2, 2]
assert best["explicit_child_root_data"] == [9, 32, 96]
assert best["exit_child_root_data"] == [9, 36, 64]
assert best["inherited_explicit_curve_gate"] == "PASS"
for prefix in ("first_edge", "return", "exit"):
    assert best[f"{prefix}_nef_audit"]["nef"]
    assert best[f"{prefix}_exact_negative_horizontal_walls"] == []

source_frame_path = ROOT / source_marking["frame_output"]
source = load_matrix(source_frame_path)
explicit = matrix(ZZ, best["explicit_child_frame"])
returned = matrix(ZZ, best["returned_frame"])
exit_child = matrix(ZZ, best["exit_child_frame"])
g_source = block_diagonal_matrix(U2, -source)
g_explicit = block_diagonal_matrix(U2, -explicit)
g_returned = block_diagonal_matrix(U2, -returned)
g_exit = block_diagonal_matrix(U2, -exit_child)
t_first = matrix(ZZ, best["source_to_explicit_child_basis"])
t_return = matrix(ZZ, best["explicit_child_to_returned_A5A5_basis"])
t_exit = matrix(ZZ, best["returned_A5A5_to_exit_3A3_basis"])
for transition in (t_first, t_return, t_exit):
    assert abs(transition.det()) == 1
assert t_first * g_source * t_first.transpose() == g_explicit
assert t_return * g_explicit * t_return.transpose() == g_returned
assert t_exit * g_returned * t_exit.transpose() == g_exit
t_new = t_exit * t_return * t_first
assert abs(t_new.det()) == 1

# Identify the new landing with the same exact current-3A3 basis as the old
# direct exit, then transport the already pinned full endpoint basis through
# that full landing bridge.
t_direct = matrix(ZZ, direct_exit["source_to_child_basis"])
direct_landing_frame = load_matrix(ROOT / direct_exit["frame_output"])
g_direct_landing = block_diagonal_matrix(U2, -direct_landing_frame)
assert t_direct * g_source * t_direct.transpose() == g_direct_landing
direct_to_new = t_direct * t_new.inverse().change_ring(ZZ)
assert direct_to_new * g_exit * direct_to_new.transpose() == g_direct_landing

canonical_old = matrix(
    ZZ, old_promoted["landing_current_3A3_identification"]["canonical_current_3A3_basis_in_landing"]
)
canonical_new = canonical_old * direct_to_new
canonical_frame = load_matrix(ROOT / current_3a3["frame_output"])
g_canonical = block_diagonal_matrix(U2, -canonical_frame)
assert abs(canonical_new.det()) == 1
assert canonical_new * g_exit * canonical_new.transpose() == g_canonical
assert vector(ZZ, canonical_new.row(0)) == vector(ZZ, [1, 0] + [0] * 17)

pinned_old = matrix(ZZ, old_promoted["endpoint"]["canonical_pinned_basis_in_landing"])
pinned_new = pinned_old * direct_to_new
pinned_frame = load_matrix(PINNED_FRAME)
g_pinned = block_diagonal_matrix(U2, -pinned_frame)
assert abs(pinned_new.det()) == 1
assert pinned_new * g_exit * pinned_new.transpose() == g_pinned

old_loop_score = int(old_promoted["new_splice"]["equation_cost_score"])
old_exit_raw = int(search["direct_q6_current"]["score"])
old_exit_inherited = int(search["direct_q6_current"]["inherited_explicit_operational_score"])
prefix_score = old_loop_score - old_exit_raw
new_second_raw = int(best["total_equation_cost_score"])
new_second_floor = int(best["inherited_explicit_equation_cost"]["operational_total_score"])
new_raw_score = prefix_score + new_second_raw
new_conservative_score = prefix_score + new_second_floor
old_conservative_score = prefix_score + old_exit_inherited
assert (old_loop_score, old_exit_raw, old_exit_inherited, prefix_score) == (10334, 3519, 3019, 6815)
assert (new_second_raw, new_second_floor, new_raw_score, new_conservative_score) == (1172, 1730, 7987, 8545)
assert new_raw_score < old_loop_score and new_conservative_score < old_conservative_score

direct_q104 = int(old_promoted["new_splice"]["direct_q104_comparator_score"])
raw_improvement = direct_q104 - new_raw_score
conservative_improvement = old_conservative_score - new_conservative_score
raw_percent = QQ(100 * raw_improvement) / direct_q104
conservative_percent = QQ(100 * conservative_improvement) / old_conservative_score

old_combined = int(d13_promoted["combined_bottleneck_comparison"]["old_direct_q24_plus_direct_q104"])
new_combined = int(d13_promoted["new_D13_splice"]["raw_equation_cost_score"]) + new_raw_score
combined_improvement = old_combined - new_combined
combined_percent = QQ(100 * combined_improvement) / old_combined
assert (old_combined, new_combined, combined_improvement) == (41403, 33310, 8093)

old_a11_q = old_promoted["full_route_q_sequence_from_A11"]
assert old_a11_q[:4] == [8, 6, 4, 6]
new_a11_q = old_a11_q[:3] + [6, 4, 8] + old_a11_q[4:]
new_h3_q = [6, 8, 4, 4, 24, 6] + new_a11_q

inputs = (
    SEARCH, SOURCE_MARKING, DIRECT_EXIT, OLD_PROMOTED, D13_PROMOTED,
    CURRENT_3A3, PINNED_FRAME, source_frame_path,
)
payload = {
    "schema": "elkies-k3.h3-a5a5-q6o1581-second-zero-promoted-route-certificate.v1",
    "status": "PASS_EXACT_PROMOTED_SECOND_ZERO_EQUATION_COST_ROUTE_TO_PINNED_R17",
    "promotion": {
        "promote_as_lifting_target": True,
        "switch_after": "the q6/o1307 q4 return reaches its equation-explicit changed-zero 2A5/MW7 model",
        "replacement": "replace the q6 current-3A3 exit by q6/o1581, q4 return, q8 exit",
        "resume_at": "the exactly identified current_3A3 stage and unchanged pinned suffix",
    },
    "second_zero_splice": {
        "nodes": [
            "q6/o1307 returned 2A5/MW7",
            "A1+A2+2A3/MW8 q6/o1581 with old_A5A5_component_0 zero",
            "2A5/MW7 with a second changed zero",
            "current 3A3/MW8",
        ],
        "q_sequence": [6, 4, 8],
        "old_fibre_degrees": [2, 2, 2],
        "profiles": {name: best[f"{name}_profile"] for name in ("first_edge", "return", "exit")},
        "inherited_explicit_curve_degrees": best["inherited_explicit_curve_degrees"],
        "raw_score": new_second_raw,
        "inherited_explicit_horizontal_floor_score": new_second_floor,
        "direct_exit_raw_score": old_exit_raw,
        "direct_exit_inherited_explicit_horizontal_floor_score": old_exit_inherited,
        "all_new_edges_exact_nef": True,
    },
    "a11_splice_comparison": {
        "direct_q104_comparator_score": direct_q104,
        "new_raw_score": new_raw_score,
        "raw_improvement": raw_improvement,
        "raw_relative_improvement_percent_exact": str(raw_percent),
        "raw_relative_improvement_percent_decimal": "{:.6f}".format(float(raw_percent)),
        "old_conservative_score": old_conservative_score,
        "new_conservative_score": new_conservative_score,
        "conservative_improvement": conservative_improvement,
        "conservative_relative_improvement_percent_exact": str(conservative_percent),
        "conservative_relative_improvement_percent_decimal": "{:.6f}".format(float(conservative_percent)),
    },
    "combined_bottleneck_comparison": {
        "old_direct_q24_plus_direct_q104": old_combined,
        "new_D13_splice_plus_second_zero_A11_splice": new_combined,
        "absolute_improvement": combined_improvement,
        "relative_improvement_percent_exact": str(combined_percent),
        "relative_improvement_percent_decimal": "{:.6f}".format(float(combined_percent)),
    },
    "landing_current_3A3_identification": {
        "canonical_current_3A3_basis_in_landing": rows(canonical_new),
        "landing_basis_in_canonical_current_3A3": rows(canonical_new.inverse().change_ring(ZZ)),
        "forward_determinant": int(canonical_new.det()),
        "inverse_determinant": int(canonical_new.inverse().det()),
        "fibre_exactly_aligned": True,
        "gram_exactly_aligned": True,
        "root_data": [9, 36, 64],
        "MW_rank": 8,
    },
    "full_route_q_sequence_from_A11": new_a11_q,
    "full_route_q_sequence_from_H3": new_h3_q,
    "endpoint": {
        "name": "pinned_R17", "root_data": [0, 0, 1], "MW_rank": 17,
        "canonical_pinned_basis_in_landing": rows(pinned_new),
        "landing_basis_in_canonical_pinned": rows(pinned_new.inverse().change_ring(ZZ)),
        "forward_determinant": int(pinned_new.det()),
        "inverse_determinant": int(pinned_new.inverse().det()),
        "gram_identification": "U plus negative pinned rank17_gram.txt exactly",
    },
    "proof_boundary": (
        "All three second-zero fibres have exact marked U, primitive isotropic classes, complete component, "
        "affine, all-section and finite horizontal-wall gates, exact root data, and bidirectional unimodular "
        "transports. The current-3A3 landing and pinned endpoint are full-basis identifications. Raw and "
        "inherited-explicit horizontal-floor scores are deterministic planning estimates, not measured runtimes."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("PROMOTED1581|second_q=6,4,8|raw={}|floor={}|old={}|combined={}|saving={}|"
      "landing_det={}|endpoint_det={}|status={}|output={}".format(
          new_second_raw, new_second_floor, old_exit_raw, new_combined, combined_improvement,
          canonical_new.det(), pinned_new.det(), payload["status"], OUTPUT,
      ), flush=True)
