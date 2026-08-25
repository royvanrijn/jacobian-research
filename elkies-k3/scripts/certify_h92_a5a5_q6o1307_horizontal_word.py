#!/usr/bin/env python3
"""Certify the low-pole Mordell--Weil word for the q6/orbit1307 exit.

status: ACTIVE_PROOF
claim: exact fixed-component reductions and MW word in the pointed 2A5 frame
inputs: exact q8 equation marking and q4/q6 lattice/cost certificates
outputs: artifacts/local/elkies-k3/q24-2a5-q6o1307-horizontal-word.json

This is deliberately dependency-free.  It uses only exact integral NS
coordinates; no section ansatz, Groebner basis, or finite-field search enters.
"""

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
MARKING = LOCAL / "q24-a11-to-2a5-q8-equation-marking-qq.json"
SCORES = GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4q6-equation-cost.json"
EDGE = GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q6-orbit1307-lattice-certificate.json"
FRAME = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frame.txt"
OUTPUT = LOCAL / "q24-2a5-q6o1307-horizontal-word.json"


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


def fixed_component_reduce(section, components, gram):
    """Remove forced physical fibre components in a deterministic chamber."""
    current = list(section)
    removed = []
    while True:
        negative = sorted(
            (pairing(current, gram, component), index)
            for index, component in enumerate(components)
            if pairing(current, gram, component) < 0
        )
        if not negative:
            return current, removed
        intersection, index = negative[0]
        assert intersection == -1, (intersection, index)
        current = add(current, scale(-1, components[index]))
        removed.append(index)


marking = json.loads(MARKING.read_text())
scores = json.loads(SCORES.read_text())
edge = json.loads(EDGE.read_text())
frame = matrix_from_text(FRAME)
assert marking["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_EQUATION_MARKING"
assert scores["status"] == "PASS_EXACT_A5A5_EXPLICIT_ZERO_EQUATION_COST_SCORING"
assert edge["status"] == "PASS_EXACT_A5A5_EXPLICIT_ZERO_CANDIDATE_LATTICE_CERTIFICATE"
assert len(frame) == 17 and all(len(row) == 17 for row in frame)

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
assert old_components[9] == old_zero
chains = marking["physical_2A5"]["chains"]
root_indices = [index for chain in chains for index in chain]
assert sorted(root_indices + [9]) == list(range(11))

for index, component in enumerate(old_components):
    assert pairing(component, gram, component) == -2
    assert pairing(component, gram, old_fibre) == (1 if index == 9 else 0)

for chain in chains:
    for left_position, left_index in enumerate(chain):
        for right_position, right_index in enumerate(chain):
            expected = -2 if left_position == right_position else (
                1 if abs(left_position - right_position) == 1 else 0
            )
            assert pairing(old_components[left_index], gram, old_components[right_index]) == expected
assert all(
    pairing(old_components[left], gram, old_components[right]) == 0
    for left in chains[0]
    for right in chains[1]
)

# Each A5 chain omits the affine component of its I6 fibre.  It is essential
# to include those two components in the chamber reduction: some raw planner
# representatives have become nef against the ten roots while remaining
# negative against an omitted affine component.
root_components = [old_components[index] for chain in chains for index in chain]
affine_components = [
    add(old_fibre, scale(-1, add(*[old_components[index] for index in chain])))
    for chain in chains
]
fibre_components = root_components + affine_components
fibre_component_labels = [
    f"old_A11_component_{index}" for chain in chains for index in chain
] + ["first_I6_affine_component", "second_I6_affine_component"]
for component in affine_components:
    assert pairing(component, gram, component) == -2
    assert pairing(component, gram, old_fibre) == 0
assert all(
    sum(pairing(section, gram, component) for component in fibre_components[offset:offset + 5])
    + pairing(section, gram, affine_components[offset // 5])
    == pairing(section, gram, old_fibre)
    for section in (old_zero,)
    for offset in (0, 5)
)


def candidate(q, orbit):
    candidate_id = {"q": q, "old_fibre_degree": 2, "orbit_index": orbit}
    return next(
        item for item in scores["ranked_candidates"]
        if item["candidate_id"] == candidate_id
    )


records = {
    "q6_orbit1307": candidate(6, 1307),
    "q6_orbit1229": candidate(6, 1229),
    "q4_orbit32": candidate(4, 32),
    "q4_orbit478": candidate(4, 478),
    "q4_orbit146": candidate(4, 146),
    "q4_orbit208": candidate(4, 208),
}
expected_reductions = {
    "q6_orbit1307": [3, 2, 1, 0],
    "q6_orbit1229": [],
    "q4_orbit32": [3, 2, 1, 0, 10],
    "q4_orbit478": [3, 2, 1, 0],
    "q4_orbit146": [3, 2, 1, 0],
    "q4_orbit208": [2, 1, 0, 3, 2, 1],
}
effective = {}
section_payload = {}
for name, record in records.items():
    raw = record["horizontal"]["section"]
    reduced, removed = fixed_component_reduce(raw, fibre_components, gram)
    assert removed == expected_reductions[name]
    assert pairing(reduced, gram, reduced) == -2
    assert pairing(reduced, gram, old_fibre) == 1
    assert pairing(reduced, gram, old_zero) == record["horizontal"]["P_dot_O"]
    component_pairings = [pairing(reduced, gram, component) for component in fibre_components]
    assert all(value >= 0 for value in component_pairings)
    effective[name] = reduced
    pole_order = record["horizontal"]["P_dot_O"]
    section_payload[name] = {
        "raw_planning_representative": raw,
        "removed_fixed_fibre_component_indices_in_order": removed,
        "removed_fixed_fibre_components_in_order": [fibre_component_labels[index] for index in removed],
        "effective_section": reduced,
        "P_dot_O": pole_order,
        "weierstrass_projective_degrees_X_Y_Z": [4 + 2 * pole_order, 6 + 3 * pole_order, pole_order],
        "all_twelve_I6_component_pairings": dict(zip(fibre_component_labels, component_pairings)),
    }

target = effective["q6_orbit1307"]
primary_word_sum = add(effective["q4_orbit146"], effective["q6_orbit1229"])
primary_word_difference = add(target, scale(-1, primary_word_sum))
assert primary_word_difference == scale(-1, old_zero)

alternate_word_sum = add(
    effective["q4_orbit32"],
    effective["q4_orbit478"],
    effective["q4_orbit146"],
    scale(-1, effective["q4_orbit208"]),
)
alternate_word_difference = add(target, scale(-1, alternate_word_sum))
alternate_trivial_correction = add(
    scale(2, old_fibre), scale(-1, old_zero),
    scale(-1, old_components[3]), scale(-1, old_components[4]),
    scale(-1, old_components[5]), scale(-1, old_components[10]),
)
assert alternate_word_difference == alternate_trivial_correction

# The planner decomposed the q6 fibre using a reducible section representative.
# Transfer the removed components to its vertical part and check the full class.
raw_target = records["q6_orbit1307"]["horizontal"]["section"]
raw_vertical = [0, 0] + records["q6_orbit1307"]["horizontal"]["vertical"] + [0] * 7
removed_sum = add(*[fibre_components[index] for index in expected_reductions["q6_orbit1307"]])
effective_vertical = add(raw_vertical, removed_sum)
new_fibre = edge["marked_U"]["fibre_in_parent"]
assert new_fibre == add(old_zero, raw_target, raw_vertical)
assert new_fibre == add(old_zero, target, effective_vertical)
effective_vertical_root_coordinates = effective_vertical[2:12]
assert effective_vertical_root_coordinates == [0, 0, 0, 0, 0, 1, 2, 1, 2, 2]

payload = {
    "schema": "elkies-k3.q24-2a5-q6o1307-horizontal-word.v1",
    "status": "PASS_EXACT_Q24_2A5_Q6O1307_LOW_POLE_HORIZONTAL_WORD",
    "pointed_parent": {
        "zero": "old_A11_component_9",
        "root_type": "2A5",
        "physical_chains_by_old_A11_component_index": chains,
    },
    "sections": section_payload,
    "mordell_weil_word": {
        "primary_identity_in_NS_modulo_trivial_lattice": "P1307 = P146 + P1229",
        "primary_target_minus_word_sum_NS_coordinates": primary_word_difference,
        "primary_exact_trivial_lattice_correction": "-O",
        "primary_exact_trivial_lattice_correction_NS_coordinates": scale(-1, old_zero),
        "alternate_identity_in_NS_modulo_trivial_lattice": "P1307 = P32 + P478 + P146 - P208",
        "alternate_target_minus_word_sum_NS_coordinates": alternate_word_difference,
        "alternate_exact_trivial_lattice_correction": "2F - O - C3 - C4 - C5 - C10",
        "alternate_exact_trivial_lattice_correction_NS_coordinates": alternate_trivial_correction,
        "torsion_ambiguity": False,
    },
    "q6_divisor_correction": {
        "fibre_NS_coordinates": new_fibre,
        "decomposition": "new_fibre = old_zero + effective_P1307 + effective_vertical",
        "old_planning_vertical_root_coordinates": records["q6_orbit1307"]["horizontal"]["vertical"],
        "effective_vertical_root_coordinates": effective_vertical_root_coordinates,
        "effective_vertical_layers": max(effective_vertical_root_coordinates),
        "effective_vertical_support": sum(value != 0 for value in effective_vertical_root_coordinates),
    },
    "compiler_plan": {
        "exact_sections_to_recover": ["q4_orbit146", "q6_orbit1229"],
        "pole_orders": [1, 0],
        "construction": "recover one simple-pole and one polynomial section, then use fraction-free elliptic addition",
        "q6_orbit1229_neighbor_nef_status_is_irrelevant": (
            "its horizontal class is already an effective section; only that section, not its rejected neighbour fibre, is used"
        ),
        "direct_q6_target_ansatz_avoided": True,
        "large_Groebner_required": False,
        "next_gate": "recover and equation-mark the two primary low-pole sections over QQ",
    },
    "proof_boundary": (
        "This certificate proves the exact NS fixed-component reductions, the q6 divisor "
        "correction, and the Mordell--Weil quotient word. It does not yet construct the four "
        "low-pole sections as characteristic-zero equations or the q6 genus-one pencil."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (MARKING, SCORES, EDGE, FRAME)],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (MARKING, SCORES, EDGE, FRAME)
        },
    },
}
LOCAL.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q6O1307WORD|word=P146+P1229|poles=1,0|alternate=P32+P478+P146-P208|"
    f"vertical={effective_vertical_root_coordinates}|status={payload['status']}|output={OUTPUT}",
    flush=True,
)
