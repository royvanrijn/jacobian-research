#!/usr/bin/env sage -python
"""Certify the q4/o11 changed-zero replacement of the first H3 q8 edge."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
SOURCE_MARKING = GENERATED / "elkies-k3-h3-first-q8-source-marking.json"
SEARCH = GENERATED / "elkies-k3-h3-first-q8-zero-changing-d13-presentations.json"
HISTORICAL_D13_FRAME = ROOT / "elkies-k3/data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt"
EQUATION_TO_PINNED = ROOT / "artifacts/local/elkies-k3/q24-equation-d13-to-pinned-r17.json"
D13_PROMOTED = GENERATED / "elkies-k3-h3-d13-q4o11-promoted-route-certificate.json"
A11_PROMOTED = GENERATED / "elkies-k3-h3-a5a5-q4o230-q6o1315-promoted-route-certificate.json"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = GENERATED / "elkies-k3-h3-first-q8-q4o11-promoted-route-certificate.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(ZZ, [[ZZ(x) for x in line.split()] for line in path.read_text().splitlines()
                       if line.strip() and not line.lstrip().startswith("#")])


def rows(value):
    return [[int(x) for x in row] for row in value.rows()]


paths = (SOURCE_MARKING, SEARCH, HISTORICAL_D13_FRAME, EQUATION_TO_PINNED,
         D13_PROMOTED, A11_PROMOTED, PINNED_FRAME)
source_marking, search, equation_to_pinned, d13, a11 = (
    json.loads(path.read_text()) for path in
    (SOURCE_MARKING, SEARCH, EQUATION_TO_PINNED, D13_PROMOTED, A11_PROMOTED)
)

assert source_marking["status"] == "PASS_EXACT_FIRST_Q8_SOURCE_MARKING"
assert source_marking["root_data"] == [14, 312, 3]
assert search["status"] == "PASS_EXACT_E8E6_ZERO_CHANGING_D13_PRESENTATION_SEARCH"
assert equation_to_pinned["status"] == "PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH"
assert d13["status"] == "PASS_EXACT_PROMOTED_D13_EQUATION_COST_ROUTE_TO_PINNED_R17"
assert a11["status"] == "PASS_EXACT_PROMOTED_DOUBLE_ZERO_EQUATION_COST_ROUTE_TO_PINNED_R17"

selected = next(
    row for row in search["ranked_presentations"]
    if row["first_edge_candidate_id"] == {"q": 4, "old_fibre_degree": 2, "orbit_index": 11}
    and row["explicit_zero_curve"] == "old_E8E6_component_1"
)
assert selected is search["ranked_presentations"][0]
assert selected["q_sequence"] == [4, 4, 4]
assert selected["old_fibre_degrees"] == [2, 2, 2]
assert selected["explicit_child_root_data"] == [14, 172, 24]
assert selected["exit_child_root_data"] == [13, 312, 4]
assert selected["inherited_explicit_curve_gate"] == "PASS"
for prefix in ("first_edge", "return", "exit"):
    assert selected[f"{prefix}_nef_audit"]["nef"]
    assert selected[f"{prefix}_exact_negative_horizontal_walls"] == []

# All matrices have rows equal to the new marked basis in coordinates of the old one.
source_frame = load_matrix(ROOT / source_marking["frame_output"])
g_source = block_diagonal_matrix(U2, -source_frame)
explicit_frame = matrix(ZZ, selected["explicit_child_frame"])
returned_frame = matrix(ZZ, selected["returned_frame"])
landing_frame = matrix(ZZ, selected["exit_child_frame"])
g_explicit = block_diagonal_matrix(U2, -explicit_frame)
g_returned = block_diagonal_matrix(U2, -returned_frame)
g_landing = block_diagonal_matrix(U2, -landing_frame)
t_first = matrix(ZZ, selected["source_to_explicit_child_basis"])
t_return = matrix(ZZ, selected["explicit_child_to_returned_E8E6_basis"])
t_exit = matrix(ZZ, selected["returned_E8E6_to_exit_D13_basis"])
assert t_first * g_source * t_first.transpose() == g_explicit
assert t_return * g_explicit * t_return.transpose() == g_returned
assert t_exit * g_returned * t_exit.transpose() == g_landing
transitions = (t_first, t_return, t_exit)
assert all(abs(item.det()) == 1 for item in transitions)
assert all(item.inverse() in MatrixSpace(ZZ, 19) for item in transitions)
t_source_to_landing = t_exit * t_return * t_first
assert abs(t_source_to_landing.det()) == 1
assert t_source_to_landing * g_source * t_source_to_landing.transpose() == g_landing

# Identify the landing with the historical selected D13 marking, and then with
# the current equation-D13 marking, using full bases rather than ADE/MW data.
historical_in_source = matrix(ZZ, source_marking["equation_D13_basis_in_root_adapted_hub"])
historical_frame = load_matrix(HISTORICAL_D13_FRAME)
g_historical = block_diagonal_matrix(U2, -historical_frame)
assert historical_in_source * g_source * historical_in_source.transpose() == g_historical
assert abs(historical_in_source.det()) == 1
historical_in_landing = historical_in_source * t_source_to_landing.inverse().change_ring(ZZ)
assert historical_in_landing * g_landing * historical_in_landing.transpose() == g_historical

historical_in_equation = matrix(ZZ, equation_to_pinned["historical_d13_basis_in_equation_d13"])
equation_frame = matrix(ZZ, equation_to_pinned["equation_d13_frame"])
g_equation = block_diagonal_matrix(U2, -equation_frame)
assert historical_in_equation * g_equation * historical_in_equation.transpose() == g_historical
assert abs(historical_in_equation.det()) == 1
equation_in_landing = historical_in_equation.inverse().change_ring(ZZ) * historical_in_landing
assert equation_in_landing * g_landing * equation_in_landing.transpose() == g_equation
assert abs(equation_in_landing.det()) == 1
assert vector(ZZ, equation_in_landing.row(0)) == vector(ZZ, [1, 0] + [0] * 17)

# Compose the independently exact canonical endpoint marking.  The cheaper
# continuation itself is the two promoted segment certificates asserted above.
pinned_in_equation = matrix(ZZ, equation_to_pinned["equation_d13_to_pinned_r17_transition"])
pinned_frame = load_matrix(PINNED_FRAME)
g_pinned = block_diagonal_matrix(U2, -pinned_frame)
assert pinned_in_equation * g_equation * pinned_in_equation.transpose() == g_pinned
assert abs(pinned_in_equation.det()) == 1
pinned_in_landing = pinned_in_equation * equation_in_landing
assert pinned_in_landing * g_landing * pinned_in_landing.transpose() == g_pinned
assert abs(pinned_in_landing.det()) == 1

direct_score = int(search["direct_q8"]["inherited_explicit_operational_score"])
new_raw = int(selected["total_equation_cost_score"])
new_operational = int(selected["inherited_explicit_equation_cost"]["operational_total_score"])
assert (direct_score, new_raw, new_operational) == (5802, 3394, 3961)
assert selected["inherited_explicit_equation_cost"]["operational_edge_scores"] == [1509, 500, 1952]
assert new_operational < direct_score

old_h3_q = a11["full_route_q_sequence_from_H3"]
assert old_h3_q[:2] == [6, 8]
new_h3_q = [6, 4, 4, 4] + old_h3_q[2:]
continuation_score = int(a11["combined_bottleneck_comparison"][
    "new_D13_splice_plus_double_zero_A11_splice_horizontal_floor"
])
old_extended = direct_score + continuation_score
new_extended = new_operational + continuation_score
assert (continuation_score, old_extended, new_extended) == (29522, 35324, 33483)

payload = {
    "schema": "elkies-k3.h3-first-q8-changed-zero-promoted-route-certificate.v1",
    "status": "PASS_EXACT_PROMOTED_FIRST_Q8_TRIPLE_Q4_ROUTE_TO_PINNED_R17",
    "promotion": {
        "promote_as_lifting_target": True,
        "switch_after": "the initial H3 q6 reaches the equation-explicit E8+E6/MW3 frame",
        "replacement": "q4/orbit11 with old_E8E6_component_1 as zero, q4 return, q4 exit to equation D13",
        "resume_at": "equation D13, then use the certified q4/o11 D13 and q4/o230 A11 optimized continuation",
        "supersedes_first_q8_only": True,
    },
    "first_q8_replacement": {
        "q_sequence": [4, 4, 4],
        "old_fibre_degrees": [2, 2, 2],
        "nodes": ["E8+E6/MW3", "A2+D5+E7/MW3",
                  "E8+E6/MW3 with changed zero", "equation D13/MW4"],
        "explicit_zero_curve": selected["explicit_zero_curve"],
        "root_data": [[14, 312, 3], [14, 172, 24], [14, 312, 3], [13, 312, 4]],
        "profiles": {prefix: selected[f"{prefix}_profile"]
                     for prefix in ("first_edge", "return", "exit")},
        "inherited_explicit_curve_degrees": selected["inherited_explicit_curve_degrees"],
        "raw_score": new_raw,
        "operational_score": new_operational,
        "direct_q8_operational_comparator_score": direct_score,
        "absolute_improvement": direct_score - new_operational,
        "relative_improvement_percent_decimal": "31.730438",
        "relative_improvement_percent_exact": "92050/2901",
        "all_edges_primitive_nef_isotropic": True,
        "all_transports_bidirectional_unimodular": True,
        "full_explicit_curve_gate": "PASS",
    },
    "equation_D13_identification": {
        "canonical_equation_D13_basis_in_landing": rows(equation_in_landing),
        "landing_basis_in_canonical_equation_D13": rows(equation_in_landing.inverse().change_ring(ZZ)),
        "forward_determinant": int(equation_in_landing.det()),
        "inverse_determinant": int(equation_in_landing.inverse().det()),
        "fibre_exactly_aligned": True,
        "gram_exactly_aligned": True,
        "root_data": [13, 312, 4],
        "MW_rank": 4,
    },
    "combined_equation_cost_comparison": {
        "unchanged_D13_plus_A11_continuation_score": continuation_score,
        "old_first_q8_plus_continuation": old_extended,
        "new_triple_q4_plus_continuation": new_extended,
        "absolute_improvement": old_extended - new_extended,
    },
    "full_route_q_sequence_from_H3": new_h3_q,
    "continuation_certificates": {
        "D13_splice": str(D13_PROMOTED.relative_to(ROOT)),
        "A11_splice_and_pinned_suffix": str(A11_PROMOTED.relative_to(ROOT)),
        "both_exact_and_unchanged": True,
    },
    "endpoint": {
        "name": "pinned_R17",
        "root_data": [0, 0, 1],
        "MW_rank": 17,
        "canonical_pinned_basis_in_landing_D13": rows(pinned_in_landing),
        "landing_D13_basis_in_canonical_pinned": rows(pinned_in_landing.inverse().change_ring(ZZ)),
        "forward_determinant": int(pinned_in_landing.det()),
        "inverse_determinant": int(pinned_in_landing.inverse().det()),
        "gram_identification": "U plus negative pinned rank17_gram.txt exactly",
    },
    "proof_boundary": (
        "All three replacement fibres are primitive nef isotropic marked-U classes. "
        "Their component, affine-component, section, and finite horizontal-wall gates pass exactly. "
        "Every 19-by-19 NS transport and inverse is integral unimodular. The landing is identified "
        "with current equation D13 and pinned rootless MW17/R17 by full exact bases. The D13 and "
        "A11 continuation certificates remain unchanged. Equation-cost scores are deterministic "
        "compiler-planning estimates, not measured characteristic-zero runtimes."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in paths],
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in paths},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("PROMOTED_FIRST_Q8|q=4,4,4|operational={}|direct={}|saving={}|combined={}|"
      "landing_det={}|endpoint_det={}|status={}|output={}".format(
          new_operational, direct_score, direct_score - new_operational, new_extended,
          equation_in_landing.det(), pinned_in_landing.det(), payload["status"], OUTPUT,
      ), flush=True)
