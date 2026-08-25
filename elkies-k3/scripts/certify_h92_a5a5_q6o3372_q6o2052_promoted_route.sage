#!/usr/bin/env sage -python
"""Certify the two-zero q3372/q2052 A11 splice through pinned R17."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--variant", choices=("q3372", "q230"), default="q3372")
args = parser.parse_args()
CONFIG = {
    "q3372": {
        "first_id": {"q": 6, "old_fibre_degree": 2, "orbit_index": 3372},
        "first_zero": "old_A11_component_5", "first_root": [10, 46, 96],
        "returned": "elkies-k3-h3-a5a5-q6o3372-c5-returned-marking.json",
        "second": "elkies-k3-h3-a5a5-q6o3372-c5-second-zero-changing-3a3-presentations.json",
        "second_id": {"q": 6, "old_fibre_degree": 2, "orbit_index": 2052},
        "second_zero": "old_A5A5_component_6", "second_root": [9, 32, 96],
        "scores": (2021, 2774, -506, 1730, 1515, 4504),
        "previous": "elkies-k3-h3-a5a5-q6o1581-second-zero-promoted-route-certificate.json",
        "previous_scores": (7987, 8545),
        "output": "elkies-k3-h3-a5a5-q6o3372-q6o2052-promoted-route-certificate.json",
        "status": "PASS_EXACT_PROMOTED_DOUBLE_ZERO_EQUATION_COST_ROUTE_TO_PINNED_R17",
    },
    "q230": {
        "first_id": {"q": 4, "old_fibre_degree": 2, "orbit_index": 230},
        "first_zero": "old_A11_component_10", "first_root": [10, 52, 60],
        "returned": "elkies-k3-h3-a5a5-q4o230-c10-returned-marking.json",
        "second": "elkies-k3-h3-a5a5-q4o230-c10-second-zero-changing-3a3-presentations.json",
        "second_id": {"q": 6, "old_fibre_degree": 2, "orbit_index": 1315},
        "second_zero": "old_A5A5_component_1", "second_root": [9, 30, 108],
        "scores": (2296, 2296, -338, 1903, 1958, 4199),
        "previous": "elkies-k3-h3-a5a5-q6o3372-q6o2052-promoted-route-certificate.json",
        "previous_scores": (1515, 4504),
        "output": "elkies-k3-h3-a5a5-q4o230-q6o1315-promoted-route-certificate.json",
        "status": "PASS_EXACT_PROMOTED_DOUBLE_ZERO_EQUATION_COST_ROUTE_TO_PINNED_R17",
    },
}
config = CONFIG[args.variant]
A11_EDGE = GENERATED / "elkies-k3-h3-a11-q8-orbit12-lattice-certificate.json"
SOURCE_ZERO = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
A11_MARKING = GENERATED / "elkies-k3-h3-current_A11-marked-frame.json"
FIRST_SEARCH = GENERATED / "elkies-k3-h3-a5a5-zero-changing-loop-search.json"
FIRST_COST = GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4q6-equation-cost.json"
RETURNED_MARKING = GENERATED / config["returned"]
SECOND_SEARCH = GENERATED / config["second"]
PREVIOUS = GENERATED / config["previous"]
OLD_PROMOTED = GENERATED / "elkies-k3-h3-a5a5-q6o1307-promoted-route-certificate.json"
OLD_START = GENERATED / "elkies-k3-h3-a5a5-q6o1307-suffix-marking.json"
OLD_RETURN = GENERATED / "elkies-k3-h3-a5a5-q6o1307-q4-return-a5a5-certificate.json"
OLD_EXIT = GENERATED / "elkies-k3-h3-a5a5-q6o1307-loop-current-3a3-certificate.json"
D13_PROMOTED = GENERATED / "elkies-k3-h3-d13-q4o11-promoted-route-certificate.json"
CURRENT_3A3 = GENERATED / "elkies-k3-h3-current_3A3-marked-frame.json"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = GENERATED / config["output"]
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(ZZ, [[ZZ(x) for x in line.split()] for line in path.read_text().splitlines()
                       if line.strip() and not line.lstrip().startswith("#")])


def rows(value):
    return [[int(x) for x in row] for row in value.rows()]


paths = (A11_EDGE, SOURCE_ZERO, A11_MARKING, FIRST_SEARCH, FIRST_COST,
         RETURNED_MARKING, SECOND_SEARCH, PREVIOUS, OLD_PROMOTED, OLD_START,
         OLD_RETURN, OLD_EXIT, D13_PROMOTED, CURRENT_3A3)
data = {path: json.loads(path.read_text()) for path in paths}
a11_edge, source_zero, a11_marking = data[A11_EDGE], data[SOURCE_ZERO], data[A11_MARKING]
first_search, first_cost = data[FIRST_SEARCH], data[FIRST_COST]
returned_marking, second_search = data[RETURNED_MARKING], data[SECOND_SEARCH]
previous, old_promoted = data[PREVIOUS], data[OLD_PROMOTED]
old_start, old_return, old_exit = data[OLD_START], data[OLD_RETURN], data[OLD_EXIT]
d13, current_3a3 = data[D13_PROMOTED], data[CURRENT_3A3]

assert a11_edge["status"] == "PASS_EXACT_A11_Q8_EQUATION_COST_LATTICE_CERTIFICATE"
assert a11_edge["edge"]["nef"] and a11_edge["edge"]["divisibility"] == 1
assert source_zero["status"] == "PASS_EXACT_A11_Q8_ORBIT12_EXPLICIT_ZERO_FRAMES"
assert a11_marking["status"] == "PASS_EXACT_CURRENT_SUFFIX_STAGE_MARKING"
assert first_search["status"] == "PASS_EXACT_ZERO_CHANGING_LOOP_SEARCH"
assert returned_marking["status"] == "PASS_EXACT_A5A5_ZERO_LOOP_RETURNED_MARKING"
assert second_search["status"] == "PASS_EXACT_A5A5_ZERO_CHANGING_3A3_PRESENTATION_SEARCH"
assert previous["status"] in {
    "PASS_EXACT_PROMOTED_SECOND_ZERO_EQUATION_COST_ROUTE_TO_PINNED_R17",
    "PASS_EXACT_PROMOTED_DOUBLE_ZERO_EQUATION_COST_ROUTE_TO_PINNED_R17",
}
assert old_promoted["status"] == "PASS_EXACT_PROMOTED_EQUATION_COST_ROUTE_TO_PINNED_R17"
assert old_start["status"] == "PASS_EXACT_A5A5_CANDIDATE_SUFFIX_MARKING"
assert old_return["status"] == old_exit["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
assert d13["status"] == "PASS_EXACT_PROMOTED_D13_EQUATION_COST_ROUTE_TO_PINNED_R17"

first = next(row for row in first_search["ranked_loops"]
             if row["first_edge_candidate_id"] == config["first_id"]
             and row["explicit_zero_curve"] == config["first_zero"])
first_gate = next(row for row in first_cost["ranked_candidates"]
                  if row["candidate_id"] == first["first_edge_candidate_id"])
assert first_gate["full_declared_nef_gate"] == "PASS"
assert first["return_nef_audit"]["nef_in_selected_component_chamber"]
assert first["return_exact_negative_horizontal_walls"] == []
assert first["explicit_child_root_data"] == config["first_root"]

second = second_search["ranked_presentations"][0]
assert second["first_edge_candidate_id"] == config["second_id"]
assert second["explicit_zero_curve"] == config["second_zero"]
assert second["q_sequence"] == [6, 4, 4]
assert second["old_fibre_degrees"] == [2, 2, 2]
assert second["explicit_child_root_data"] == config["second_root"]
assert second["exit_child_root_data"] == [9, 36, 64]
assert second["inherited_explicit_curve_gate"] == "PASS"
for prefix in ("first_edge", "return", "exit"):
    assert second[f"{prefix}_nef_audit"]["nef"]
    assert second[f"{prefix}_exact_negative_horizontal_walls"] == []

# The q8 source supplies the exact equation-A11 Gram and first marked U.
zero = source_zero["selected"]
source_frame = matrix(ZZ, zero["frame"])
g_source = block_diagonal_matrix(U2, -source_frame)
t_q8 = matrix(ZZ, zero["equation_A11_to_explicit_zero_basis"])
t_q8_inv = t_q8.inverse().change_ring(ZZ)
g_equation = t_q8_inv * g_source * t_q8_inv.transpose()
assert abs(t_q8.det()) == 1

# First changed-zero loop: q6/o3372 and q4 return.
t_equation_to_explicit1 = matrix(ZZ, first["equation_A11_to_explicit_child_basis"])
t_first1 = t_equation_to_explicit1 * t_q8_inv
t_return1 = matrix(ZZ, first["return_transition"])
explicit1_frame = -(t_first1 * g_source * t_first1.transpose())[2:, 2:]
g_explicit1 = block_diagonal_matrix(U2, -explicit1_frame)
returned1_frame = matrix(ZZ, first["returned_A5A5_frame"])
g_returned1 = block_diagonal_matrix(U2, -returned1_frame)
assert t_first1 * g_source * t_first1.transpose() == g_explicit1
assert t_return1 * g_explicit1 * t_return1.transpose() == g_returned1
t_equation_to_returned1 = t_return1 * t_equation_to_explicit1
assert t_equation_to_returned1 == matrix(
    ZZ, returned_marking["equation_A11_to_root_adapted_hub_basis"]
)

# Second changed-zero loop and its q4 current-3A3 exit.
explicit2_frame = matrix(ZZ, second["explicit_child_frame"])
returned2_frame = matrix(ZZ, second["returned_frame"])
exit_frame = matrix(ZZ, second["exit_child_frame"])
g_explicit2 = block_diagonal_matrix(U2, -explicit2_frame)
g_returned2 = block_diagonal_matrix(U2, -returned2_frame)
g_exit = block_diagonal_matrix(U2, -exit_frame)
t_first2 = matrix(ZZ, second["source_to_explicit_child_basis"])
t_return2 = matrix(ZZ, second["explicit_child_to_returned_A5A5_basis"])
t_exit = matrix(ZZ, second["returned_A5A5_to_exit_3A3_basis"])
assert t_first2 * g_returned1 * t_first2.transpose() == g_explicit2
assert t_return2 * g_explicit2 * t_return2.transpose() == g_returned2
assert t_exit * g_returned2 * t_exit.transpose() == g_exit

transitions = (t_q8, t_first1, t_return1, t_first2, t_return2, t_exit)
assert all(abs(item.det()) == 1 for item in transitions)
assert all(abs(item.inverse().change_ring(ZZ).det()) == 1 for item in transitions)
t_equation_to_exit = t_exit * t_return2 * t_first2 * t_equation_to_returned1
assert abs(t_equation_to_exit.det()) == 1

# Full-basis identification through the independently certified old 3A3 landing.
t_old_start = matrix(ZZ, old_start["equation_A11_to_root_adapted_hub_basis"])
t_old_return = matrix(ZZ, old_return["source_to_child_basis"])
t_old_exit = matrix(ZZ, old_exit["source_to_child_basis"])
old_start_frame = load_matrix(ROOT / old_start["frame_output"])
old_return_frame = load_matrix(ROOT / old_return["frame_output"])
old_exit_frame = load_matrix(ROOT / old_exit["frame_output"])
g_old_start = block_diagonal_matrix(U2, -old_start_frame)
g_old_return = block_diagonal_matrix(U2, -old_return_frame)
g_old_exit = block_diagonal_matrix(U2, -old_exit_frame)
assert t_old_start * g_equation * t_old_start.transpose() == g_old_start
assert t_old_return * g_old_start * t_old_return.transpose() == g_old_return
assert t_old_exit * g_old_return * t_old_exit.transpose() == g_old_exit
t_equation_to_old_exit = t_old_exit * t_old_return * t_old_start
old_to_new = t_equation_to_old_exit * t_equation_to_exit.inverse().change_ring(ZZ)
assert old_to_new * g_exit * old_to_new.transpose() == g_old_exit
canonical_3a3_landing = matrix(
    ZZ, old_promoted["landing_current_3A3_identification"]["canonical_current_3A3_basis_in_landing"]
) * old_to_new
canonical_3a3_frame = load_matrix(ROOT / current_3a3["frame_output"])
g_canonical_3a3 = block_diagonal_matrix(U2, -canonical_3a3_frame)
assert abs(canonical_3a3_landing.det()) == 1
assert canonical_3a3_landing * g_exit * canonical_3a3_landing.transpose() == g_canonical_3a3
assert vector(ZZ, canonical_3a3_landing.row(0)) == vector(ZZ, [1, 0] + [0] * 17)

pinned_landing = matrix(ZZ, old_promoted["endpoint"]["canonical_pinned_basis_in_landing"]) * old_to_new
pinned_frame = load_matrix(PINNED_FRAME)
g_pinned = block_diagonal_matrix(U2, -pinned_frame)
assert abs(pinned_landing.det()) == 1
assert pinned_landing * g_exit * pinned_landing.transpose() == g_pinned

first_raw = int(first["first_edge_equation_cost_score"] + first["return_equation_cost_score"])
first_floor = sum(max(500, int(value)) for value in
                  (first["first_edge_equation_cost_score"], first["return_equation_cost_score"]))
second_raw = int(second["total_equation_cost_score"])
second_floor = int(second["inherited_explicit_equation_cost"]["operational_total_score"])
new_raw = first_raw + second_raw
new_floor = first_floor + second_floor
assert (first_raw, first_floor, second_raw, second_floor, new_raw, new_floor) == config["scores"]
if "a11_splice_comparison" in previous:
    previous_raw = int(previous["a11_splice_comparison"]["new_raw_score"])
    previous_floor = int(previous["a11_splice_comparison"]["new_conservative_score"])
    direct = int(previous["a11_splice_comparison"]["direct_q104_comparator_score"])
else:
    previous_raw = int(previous["a11_splice"]["raw_score"])
    previous_floor = int(previous["a11_splice"]["inherited_explicit_horizontal_floor_score"])
    direct = int(previous["a11_splice"]["direct_q104_comparator_score"])
assert (previous_raw, previous_floor) == config["previous_scores"]
assert direct == 13518
assert new_floor < previous_floor

d13_score = int(d13["new_D13_splice"]["raw_equation_cost_score"])
old_combined = int(d13["combined_bottleneck_comparison"]["old_direct_q24_plus_direct_q104"])
new_combined_raw = d13_score + new_raw
new_combined_floor = d13_score + new_floor
assert (d13_score, old_combined) == (25323, 41403)
assert new_combined_raw == d13_score + new_raw
assert new_combined_floor == d13_score + new_floor

old_a11_q = previous["full_route_q_sequence_from_A11"]
suffix_q = old_a11_q[6:]
new_a11_q = [8, config["first_id"]["q"], 4, 6, 4, 4] + suffix_q
new_h3_q = [6, 8, 4, 4, 24, 6] + new_a11_q

inputs = paths + (PINNED_FRAME, ROOT / old_start["frame_output"], ROOT / old_return["frame_output"],
                  ROOT / old_exit["frame_output"], ROOT / current_3a3["frame_output"])
payload = {
    "schema": "elkies-k3.h3-a5a5-double-zero-promoted-route-certificate.v1",
    "status": config["status"],
    "promotion": {
        "promote_as_lifting_target": True,
        "switch_after": "A11 q8/orbit12 reaches its equation-explicit 2A5/MW7 zero",
        "replacement": "q{}/o{}, q4 return, q6/o{}, q4 return, q4 exit".format(
            config["first_id"]["q"], config["first_id"]["orbit_index"],
            config["second_id"]["orbit_index"],
        ),
        "resume_at": "the fully identified current_3A3 stage and unchanged pinned suffix",
    },
    "a11_splice": {
        "q_sequence": [config["first_id"]["q"], 4, 6, 4, 4],
        "old_fibre_degrees": [2, 2, 2, 2, 2],
        "nodes": ["2A5/MW7", first["first_edge_child"]["ade"] + "/MW" +
                  str(first["first_edge_child"]["mw_rank"]), "2A5/MW7 changed zero",
                  second["first_edge_child"]["ade"] + "/MW" +
                  str(second["first_edge_child"]["mw_rank"]),
                  "2A5/MW7 second changed zero", "current 3A3/MW8"],
        "first_zero": config["first_zero"],
        "second_zero": config["second_zero"],
        "root_data": [[10, 60, 36], config["first_root"], [10, 60, 36],
                      config["second_root"], [10, 60, 36], [9, 36, 64]],
        "profiles": {
            "first_return": first["return_profile"],
            "second_first": second["first_edge_profile"],
            "second_return": second["return_profile"],
            "exit": second["exit_profile"],
        },
        "inherited_explicit_curve_degrees": second["inherited_explicit_curve_degrees"],
        "raw_score": new_raw,
        "inherited_explicit_horizontal_floor_score": new_floor,
        "direct_q104_comparator_score": direct,
        "previous_promoted_raw_score": previous_raw,
        "previous_promoted_horizontal_floor_score": previous_floor,
        "strict_operational_score_improvement": True,
        "raw_credit_score_improvement": new_raw < previous_raw,
        "all_edges_exact_nef": True,
        "all_transports_bidirectional_unimodular": True,
    },
    "combined_bottleneck_comparison": {
        "old_direct_q24_plus_direct_q104": old_combined,
        "new_D13_splice_plus_double_zero_A11_splice_raw": new_combined_raw,
        "new_D13_splice_plus_double_zero_A11_splice_horizontal_floor": new_combined_floor,
        "raw_absolute_improvement": old_combined - new_combined_raw,
        "horizontal_floor_absolute_improvement": old_combined - new_combined_floor,
    },
    "landing_current_3A3_identification": {
        "canonical_current_3A3_basis_in_landing": rows(canonical_3a3_landing),
        "landing_basis_in_canonical_current_3A3": rows(canonical_3a3_landing.inverse().change_ring(ZZ)),
        "forward_determinant": int(canonical_3a3_landing.det()),
        "inverse_determinant": int(canonical_3a3_landing.inverse().det()),
        "fibre_exactly_aligned": True, "gram_exactly_aligned": True,
        "root_data": [9, 36, 64], "MW_rank": 8,
    },
    "full_route_q_sequence_from_A11": new_a11_q,
    "full_route_q_sequence_from_H3": new_h3_q,
    "endpoint": {
        "name": "pinned_R17", "root_data": [0, 0, 1], "MW_rank": 17,
        "canonical_pinned_basis_in_landing": rows(pinned_landing),
        "landing_basis_in_canonical_pinned": rows(pinned_landing.inverse().change_ring(ZZ)),
        "forward_determinant": int(pinned_landing.det()),
        "inverse_determinant": int(pinned_landing.inverse().det()),
        "gram_identification": "U plus negative pinned rank17_gram.txt exactly",
    },
    "proof_boundary": (
        "The initial q8 and all five replacement fibres are primitive nef isotropic marked-U classes. "
        "All replacement gates include components, affine components, all sections, and finite horizontal walls. "
        "Every NS transport and inverse is integral unimodular; current 3A3 and pinned R17 are identified by "
        "full exact bases. Cost scores are deterministic compiler-planning estimates, not measured runtimes."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                   for path in inputs},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("PROMOTED_DOUBLE_ZERO|variant={}|q={},4,6,4,4|raw={}|floor={}|combined_raw={}|saving={}|"
      "landing_det={}|endpoint_det={}|status={}|output={}".format(
          args.variant, config["first_id"]["q"], new_raw, new_floor,
          new_combined_raw, old_combined - new_combined_raw,
          canonical_3a3_landing.det(), pinned_landing.det(), payload["status"], OUTPUT,
      ), flush=True)
