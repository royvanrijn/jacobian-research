#!/usr/bin/env python3
"""Certify a low-pole MW word for the promoted q4/orbit230 edge.

status: ACTIVE_PROOF
claim: exact fixed-component reduction and MW word in the pointed 2A5 frame
output: artifacts/local/elkies-k3/q24-2a5-q4o230-horizontal-word.json

The shared physical-fibre reductions for P146 and P208 come from the existing
q6/o1307 word certificate.  This checker only reduces the new q4/o230 target
and verifies the resulting integral NS identity.  No Groebner basis or
finite-field search enters.
"""

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
MARKING = LOCAL / "q24-a11-to-2a5-q8-equation-marking-qq.json"
SCORES = GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4q6-equation-cost.json"
OLD_WORD = LOCAL / "q24-2a5-q6o1307-horizontal-word.json"
ROUTE = GENERATED / "elkies-k3-h3-a5a5-q4o230-q6o1315-promoted-route-certificate.json"
FRAME = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frame.txt"
OUTPUT = LOCAL / "q24-2a5-q4o230-horizontal-word.json"


def matrix_from_text(path):
    return [
        [int(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def add(*vectors):
    return [sum(values) for values in zip(*vectors)]


def scale(multiplier, vector):
    return [multiplier * value for value in vector]


def pairing(left, gram, right):
    return sum(
        left[row] * gram[row][column] * right[column]
        for row in range(len(left))
        for column in range(len(right))
    )


marking = json.loads(MARKING.read_text())
scores = json.loads(SCORES.read_text())
old_word = json.loads(OLD_WORD.read_text())
route = json.loads(ROUTE.read_text())
frame = matrix_from_text(FRAME)
assert marking["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_EQUATION_MARKING"
assert scores["status"] == "PASS_EXACT_A5A5_EXPLICIT_ZERO_EQUATION_COST_SCORING"
assert old_word["status"] == "PASS_EXACT_Q24_2A5_Q6O1307_LOW_POLE_HORIZONTAL_WORD"
assert route["status"] == "PASS_EXACT_PROMOTED_DOUBLE_ZERO_EQUATION_COST_ROUTE_TO_PINNED_R17"
assert route["a11_splice"]["q_sequence"] == [4, 4, 6, 4, 4]

gram = [[0] * 19 for _ in range(19)]
gram[0][1] = gram[1][0] = 1
for row in range(17):
    for column in range(17):
        gram[row + 2][column + 2] = -frame[row][column]

old_fibre = [1, 0] + [0] * 17
old_zero = [-1, 1] + [0] * 17
old_components = [
    marking["physical_2A5"]["child_coordinates"][f"old_A11_component_{index}"]
    for index in range(11)
]
chains = marking["physical_2A5"]["chains"]
root_indices = [index for chain in chains for index in chain]
root_components = [old_components[index] for index in root_indices]
affine_components = [
    add(old_fibre, scale(-1, add(*[old_components[index] for index in chain])))
    for chain in chains
]
fibre_components = root_components + affine_components
fibre_component_labels = [
    f"old_A11_component_{index}" for index in root_indices
] + ["first_I6_affine_component", "second_I6_affine_component"]

record = next(
    item for item in scores["ranked_candidates"]
    if item["candidate_id"] == {"q": 4, "old_fibre_degree": 2, "orbit_index": 230}
)
raw_target = record["horizontal"]["section"]
removed_indices = [1, 2, 3]
removed_sum = add(*[fibre_components[index] for index in removed_indices])
target = add(raw_target, scale(-1, removed_sum))
target_pairings = [pairing(target, gram, component) for component in fibre_components]
assert pairing(target, gram, target) == -2
assert pairing(target, gram, old_fibre) == 1
assert pairing(target, gram, old_zero) == 2
assert all(value >= 0 for value in target_pairings)

p146 = old_word["sections"]["q4_orbit146"]["effective_section"]
p208 = old_word["sections"]["q4_orbit208"]["effective_section"]
word_difference = add(target, scale(-1, p146), scale(-1, p208))
trivial_correction = add(
    scale(-1, old_fibre),
    scale(-1, old_zero),
    old_components[0],
    old_components[3],
    scale(-1, old_components[5]),
    scale(-1, old_components[10]),
)
assert word_difference == trivial_correction

raw_vertical = [0, 0] + record["horizontal"]["vertical"] + [0] * 7
effective_vertical = add(raw_vertical, removed_sum)
new_fibre = add(old_zero, raw_target, raw_vertical)
assert new_fibre == add(old_zero, target, effective_vertical)
assert pairing(new_fibre, gram, new_fibre) == 0
assert pairing(new_fibre, gram, old_fibre) == 2

payload = {
    "schema": "elkies-k3.q24-2a5-q4o230-horizontal-word.v1",
    "status": "PASS_EXACT_Q24_2A5_Q4O230_LOW_POLE_HORIZONTAL_WORD",
    "active_route": {
        "q_sequence_to_current_3A3": [4, 4, 6, 4, 4],
        "nodes": route["a11_splice"]["nodes"],
        "operational_horizontal_floor_score": route["a11_splice"][
            "inherited_explicit_horizontal_floor_score"
        ],
    },
    "pointed_parent": {
        "zero": "old_A11_component_9",
        "root_type": "2A5",
        "physical_chains_by_old_A11_component_index": chains,
    },
    "q4_orbit230_horizontal": {
        "raw_planning_representative": raw_target,
        "removed_fixed_fibre_components_in_order": [
            fibre_component_labels[index] for index in removed_indices
        ],
        "effective_section": target,
        "P_dot_O": 2,
        "weierstrass_projective_degrees_X_Y_Z": [8, 12, 2],
        "all_twelve_I6_component_pairings": dict(
            zip(fibre_component_labels, target_pairings)
        ),
    },
    "mordell_weil_word": {
        "identity_in_NS_modulo_trivial_lattice": "P230 = P146 + P208",
        "target_minus_word_sum_NS_coordinates": word_difference,
        "exact_trivial_lattice_correction": "-F - O + C0 + C3 - C5 - C10",
        "exact_trivial_lattice_correction_NS_coordinates": trivial_correction,
        "summand_pole_orders": [1, 1],
        "torsion_ambiguity": False,
    },
    "q4_divisor_correction": {
        "fibre_NS_coordinates": new_fibre,
        "decomposition": "new_fibre = old_zero + effective_P230 + effective_vertical",
        "old_planning_vertical_root_coordinates": record["horizontal"]["vertical"],
        "effective_vertical_NS_coordinates": effective_vertical,
    },
    "compiler_plan": {
        "exact_sections_to_recover": ["q4_orbit146", "q4_orbit208"],
        "pole_orders": [1, 1],
        "construction": "recover two simple-pole sections, then use fraction-free elliptic addition",
        "direct_q4_target_ansatz_avoided": True,
        "large_Groebner_required": False,
        "next_gate": "recover and equation-mark P146 and P208 over QQ",
    },
    "proof_boundary": (
        "This proves the exact NS fixed-component reduction, q4 divisor correction, "
        "and Mordell--Weil quotient word. It does not yet construct P146 or P208 as "
        "characteristic-zero equations or compile the q4 genus-one pencil."
    ),
    "inputs": {
        "paths": [
            str(path.relative_to(ROOT))
            for path in (MARKING, SCORES, OLD_WORD, ROUTE, FRAME)
        ],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (MARKING, SCORES, OLD_WORD, ROUTE, FRAME)
        },
    },
}
LOCAL.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O230WORD|word=P146+P208|poles=1,1|"
    f"removed={tuple(removed_indices)}|status={payload['status']}|output={OUTPUT}",
    flush=True,
)
