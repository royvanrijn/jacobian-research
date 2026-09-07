#!/usr/bin/env sage -python
"""Audit and repair q6/orbit1307 against the exact physical 2I6 marking.

The historical candidate passed an affine gate in an abstract A5+A5 root
chamber.  In the equation marking its fibre has negative intersection with
the first *physical* I6 affine component.  Replay the three required physical
Weyl reflections, compose the full NS transport, and list the genuinely
degree-one explicit old-A11 curves.  This is lattice/marking work only; it
does not compile the repaired RR pencil or a downstream continuation.
"""

import hashlib
import json
from pathlib import Path

from sage.all import (
    ZZ, block_diagonal_matrix, gcd, identity_matrix, matrix, vector,
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
MARKING = LOCAL / "q24-a11-to-2a5-q8-equation-marking-qq.json"
WORD = LOCAL / "q24-2a5-q6o1307-horizontal-word.json"
EDGE = GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q6-orbit1307-lattice-certificate.json"
SCORES = GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4q6-equation-cost.json"
FRAME = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frame.txt"
OUTPUT = LOCAL / "q24-2a5-q6o1307-physical-nef-audit.json"


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


marking = json.loads(MARKING.read_text())
word = json.loads(WORD.read_text())
edge = json.loads(EDGE.read_text())
scores = json.loads(SCORES.read_text())
assert marking["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_EQUATION_MARKING"
assert word["status"] == "PASS_EXACT_Q24_2A5_Q6O1307_LOW_POLE_HORIZONTAL_WORD"
assert edge["status"] == "PASS_EXACT_A5A5_EXPLICIT_ZERO_CANDIDATE_LATTICE_CERTIFICATE"
assert scores["status"] == "PASS_EXACT_A5A5_EXPLICIT_ZERO_EQUATION_COST_SCORING"

frame = matrix(ZZ, [
    [ZZ(value) for value in line.split()]
    for line in FRAME.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
gram = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -frame)
old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
components = {
    index: vector(ZZ, marking["physical_2A5"]["child_coordinates"][
        f"old_A11_component_{index}"
    ])
    for index in range(11)
}
chains = marking["physical_2A5"]["chains"]
affine = {
    index: old_fibre - sum(
        (components[item] for item in chain), vector(ZZ, [0] * 19)
    )
    for index, chain in enumerate(chains)
}
physical_curves = [components[index] for index in range(11)] + [affine[0], affine[1]]
physical_labels = [f"old_A11_component_{index}" for index in range(11)] + [
    "first_I6_affine_component", "second_I6_affine_component",
]
assert all(curve * gram * curve == -2 for curve in physical_curves)

stored_fibre = vector(ZZ, edge["marked_U"]["fibre_in_parent"])
assert list(stored_fibre) == word["q6_divisor_correction"]["fibre_NS_coordinates"]
score = next(
    item for item in scores["ranked_candidates"]
    if item["candidate_id"] == {"q": 6, "old_fibre_degree": 2, "orbit_index": 1307}
)
assert score["full_declared_nef_gate"] == "PASS"
assert score["parent_affine_component_pairings"] == [2, 0]

stored_physical_pairings = [int(stored_fibre * gram * curve) for curve in physical_curves]
assert stored_physical_pairings[-2:] == [-1, 0]


def reflection(curve):
    # Row-coordinate action v |-> v+(v.curve)curve for curve^2=-2.
    answer = identity_matrix(ZZ, 19) + (gram * curve.column()) * matrix(ZZ, [list(curve)])
    assert answer * gram * answer.transpose() == gram
    assert answer.det() == -1
    return answer


repaired_fibre = stored_fibre
weyl = identity_matrix(ZZ, 19)
reflection_labels = []
for unused in range(30):
    negative = [
        (int(repaired_fibre * gram * curve), index)
        for index, curve in enumerate(physical_curves)
        if repaired_fibre * gram * curve < 0
    ]
    if not negative:
        break
    value, index = min(negative)
    assert value == -1
    curve = physical_curves[index]
    step = reflection(curve)
    repaired_fibre = repaired_fibre * step
    weyl = weyl * step
    reflection_labels.append(physical_labels[index])
else:
    raise ArithmeticError("physical I6 Weyl reduction did not terminate")

assert reflection_labels == [
    "first_I6_affine_component", "old_A11_component_0", "old_A11_component_3",
]
assert repaired_fibre == stored_fibre * weyl
assert repaired_fibre * gram * repaired_fibre == 0
assert repaired_fibre * gram * old_fibre == 2
assert gcd([abs(int(value)) for value in gram * repaired_fibre]) == 1
repaired_pairings = [int(repaired_fibre * gram * curve) for curve in physical_curves]
assert all(value >= 0 for value in repaired_pairings)
assert list(repaired_fibre) == [
    3, 2, 0, 1, 1, -3, -3, 1, -4, -4, -2, -2, -1, 0, 0, 0, -1, -1, 1,
]

horizontal = vector(ZZ, word["sections"]["q6_orbit1307"]["effective_section"])
vertical = repaired_fibre - old_zero - horizontal
assert vertical == -(components[3] + components[4] + components[5])
degree_one = [
    physical_labels[index]
    for index, value in enumerate(repaired_pairings)
    if value == 1 and index < 11
]
assert degree_one == [
    "old_A11_component_3", "old_A11_component_5", "old_A11_component_9",
]
assert repaired_pairings[10] == 0  # old_A11_component_10 is not a section.

# Apply the same parent Weyl isometry to every row of the stored child basis.
# The abstract child frame/root data stay unchanged, while the parent marking
# and the list of genuinely available explicit zeros change.
stored_parent_to_child = matrix(ZZ, edge["transport"]["parent_to_child_basis"])
stored_equation_to_child = matrix(ZZ, edge["transport"]["equation_A11_to_child_basis"])
parent_to_equation = stored_parent_to_child.inverse().change_ring(ZZ) * stored_equation_to_child
repaired_parent_to_child = stored_parent_to_child * weyl
repaired_child_to_parent = repaired_parent_to_child.inverse().change_ring(ZZ)
repaired_equation_to_child = repaired_parent_to_child * parent_to_equation
repaired_child_to_equation = repaired_equation_to_child.inverse().change_ring(ZZ)
assert vector(ZZ, repaired_parent_to_child.row(0)) == repaired_fibre
assert abs(repaired_parent_to_child.det()) == 1
assert repaired_parent_to_child * gram * repaired_parent_to_child.transpose() == (
    stored_parent_to_child * gram * stored_parent_to_child.transpose()
)

inputs = (MARKING, WORD, EDGE, SCORES, FRAME)
payload = {
    "schema": "elkies-k3.q24-2a5-q6o1307-physical-nef-audit.v1",
    "status": "PASS_EXACT_Q6O1307_PHYSICAL_WEYL_REPAIR_REJECT_C10_ZERO",
    "selection": {
        "candidate_id": {"q": 6, "old_fibre_degree": 2, "orbit_index": 1307},
    },
    "stored_candidate_audit": {
        "stored_fibre": list(map(int, stored_fibre)),
        "stored_abstract_affine_pairings": score["parent_affine_component_pairings"],
        "actual_physical_I6_pairings": dict(zip(physical_labels, stored_physical_pairings)),
        "failure": "negative intersection with first_I6_affine_component",
    },
    "physical_weyl_repair": {
        "reflection_sequence": reflection_labels,
        "repaired_fibre": list(map(int, repaired_fibre)),
        "repaired_physical_I6_pairings": dict(zip(physical_labels, repaired_pairings)),
        "old_fibre_degree": 2,
        "primitive": True,
        "isotropic": True,
        "horizontal_section": "q6_orbit1307 / P1307",
        "vertical_correction": "-(old_A11_component_3+old_A11_component_4+old_A11_component_5)",
        "vertical_connected_layers": 1,
        "expected_RR_dimensions": [9, 3, 2],
        "genuinely_degree_one_old_A11_components": degree_one,
        "old_A11_component_10_degree": repaired_pairings[10],
    },
    "transport": {
        "parent_weyl_isometry": rows(weyl),
        "child_to_parent_basis": rows(repaired_child_to_parent),
        "parent_to_child_basis": rows(repaired_parent_to_child),
        "child_to_equation_A11_basis": rows(repaired_child_to_equation),
        "equation_A11_to_child_basis": rows(repaired_equation_to_child),
        "forward_determinant": int(repaired_parent_to_child.det()),
        "inverse_determinant": int(repaired_child_to_parent.det()),
        "Gram_transport_exact": True,
    },
    "child": edge["child"],
    "proof_boundary": (
        "This exact audit rejects the stored physical chamber and component-10 zero, "
        "then supplies the physical-I6 Weyl-repaired fibre and full unimodular transport. "
        "It proves nonnegativity only on the twelve displayed physical I6 components; a "
        "finite horizontal-wall audit, an RR lift, and a continuation from one of components "
        "3, 5, or 9 remain separate gates."
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
    "Q6O1307PHYSICALNEF|stored_affine=-1,0|reflections={}|degree1={}|"
    "C10degree={}|RR=9,3,2|status={}|output={}".format(
        ",".join(reflection_labels), ",".join(degree_one), repaired_pairings[10],
        payload["status"], OUTPUT,
    ),
    flush=True,
)
