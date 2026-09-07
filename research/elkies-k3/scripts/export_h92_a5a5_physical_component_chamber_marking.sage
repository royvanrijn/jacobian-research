#!/usr/bin/env sage -python
"""Export the actual physical-component chamber of the component-9-zero 2A5 hub.

The older ``physical-source`` artifact retained an abstract A5+A5 root basis.
Here the first ten frame vectors are the ten effective nonidentity components
of the two physical I6 fibres, including both affine components.  Dominant
Weyl representatives in this frame therefore satisfy the physical component
gate before any candidate scoring.
"""

import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, identity_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
ZERO = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
PHYSICAL = LOCAL / "q24-a11-to-2a5-q8-equation-marking-qq.json"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
FINGERPRINT = LOCAL / "q24-a11-q8-construction-fingerprint.json"
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
MISMATCH = GENERATED / "elkies-k3-h3-a11-quintic-bridge-zero-mismatch.json"
P230_WORD = LOCAL / "q24-2a5-q4o230-horizontal-word.json"
P1229 = LOCAL / "q24-2a5-p1229-scaled-x-qq.json"
P146 = LOCAL / "q24-2a5-p146-p1307-scaled-x-qq.json"
Q10 = LOCAL / "q24-2a5-direct-physical-q10-certificate.json"
Q10_EFFECTIVE = LOCAL / "q24-2a5-direct-physical-q10-effective-c5-zero-certificate.json"
FRAME_OUTPUT = GENERATED / "elkies-k3-h3-a5a5-physical-component-chamber-frame.txt"
OUTPUT = GENERATED / "elkies-k3-h3-a5a5-physical-component-chamber-marking.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


zero = json.loads(ZERO.read_text())
physical = json.loads(PHYSICAL.read_text())
manifest = json.loads(MANIFEST.read_text())
fingerprint = json.loads(FINGERPRINT.read_text())
crossovers = json.loads(CROSSOVERS.read_text())
mismatch = json.loads(MISMATCH.read_text())
p230_word = json.loads(P230_WORD.read_text())
p1229 = json.loads(P1229.read_text())
p146 = json.loads(P146.read_text())
q10 = json.loads(Q10.read_text())
q10_effective = json.loads(Q10_EFFECTIVE.read_text())
assert zero["status"] == "PASS_EXACT_A11_Q8_ORBIT12_EXPLICIT_ZERO_FRAMES"
assert zero["selected_zero_curve"] == "old_A11_component_9"
assert physical["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_EQUATION_MARKING"
assert manifest["status"] == "PASS_H3_R17_BACKWARD_EXACT_LIFT_MANIFEST"
assert crossovers["status"] == "PASS_EXACT_MARKED_TARGET_CROSSOVER_AUDIT"
assert mismatch["status"] == "REJECT_A11_QUINTIC_BRIDGE_ZERO_MISMATCH"
assert p230_word["status"] == "PASS_EXACT_Q24_2A5_Q4O230_LOW_POLE_HORIZONTAL_WORD"
assert p1229["status"] == "PASS_EXACT_QQ_P1229_POLYNOMIAL_SECTION"
assert p146["status"] == "PASS_EXACT_QQ_P146_AND_P1307_SHORT_MW_WORDS"
assert q10["status"] == "PASS_EXACT_PHYSICAL_NEF_Q10_CURRENT_3A3_PRESENTATION"
assert q10_effective["status"] == "PASS_EXACT_PHYSICAL_Q10_EFFECTIVE_C5_ZERO_TO_PINNED_R17"

selected = zero["selected"]
abstract_frame_path = ROOT / zero["selected_frame_output"]
abstract_frame = load_matrix(abstract_frame_path)
assert abstract_frame == matrix(ZZ, selected["frame"])
abstract_gram = block_diagonal_matrix(U2, -abstract_frame)
old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
components = {
    index: vector(ZZ, physical["physical_2A5"]["child_coordinates"][
        f"old_A11_component_{index}"
    ])
    for index in range(11)
}
chains = physical["physical_2A5"]["chains"]
affines = [
    old_fibre - sum((components[index] for index in chain), vector(ZZ, [0] * 19))
    for chain in chains
]
assert components[9] * abstract_gram * old_fibre == 1

# Select every physical component not met by the equation zero C9.  Each I6
# contributes four old A11 components and its affine component.
effective_simple = []
effective_simple_names = []
for chain_index, chain in enumerate(chains):
    for index in chain:
        curve = components[index]
        if old_zero * abstract_gram * curve == 0:
            effective_simple.append(curve)
            effective_simple_names.append(f"old_A11_component_{index}")
    effective_simple.append(affines[chain_index])
    effective_simple_names.append(f"I6_affine_{chain_index}")
assert effective_simple_names == [
    "old_A11_component_0", "old_A11_component_3", "old_A11_component_4",
    "old_A11_component_5", "I6_affine_0", "old_A11_component_1",
    "old_A11_component_2", "old_A11_component_6", "old_A11_component_7",
    "I6_affine_1",
]

root_rows = matrix(ZZ, [-curve[2:] for curve in effective_simple])
smith, _, smith_right = root_rows.smith_form()
assert tuple(abs(smith[index, index]) for index in range(10)) == (1,) * 10
completion = smith_right.inverse().change_ring(ZZ)
adaptation = root_rows.stack(completion[10:])
assert abs(adaptation.det()) == 1
physical_frame = adaptation * abstract_frame * adaptation.transpose()
change = block_diagonal_matrix(identity_matrix(ZZ, 2), adaptation)
physical_gram = block_diagonal_matrix(U2, -physical_frame)
assert physical_gram == change * abstract_gram * change.transpose()
physical_cartan = physical_frame[:10, :10]
assert physical_cartan.det() == 36
assert all(physical_cartan[index, index] == 2 for index in range(10))
assert all(
    physical_cartan[left, right] in (0, -1)
    for left in range(10) for right in range(10) if left != right
)

FRAME_OUTPUT.write_text(
    "# Physical-component A5+A5 frame; generated by "
    "export_h92_a5a5_physical_component_chamber_marking.sage\n"
    + "\n".join(" ".join(map(str, row)) for row in rows(physical_frame)) + "\n"
)

# Full equation-A11 and suffix transports.
abstract_in_equation = matrix(ZZ, selected["equation_A11_to_explicit_zero_basis"])
physical_in_equation = change * abstract_in_equation
equation_in_physical = physical_in_equation.inverse().change_ring(ZZ)
assert abs(physical_in_equation.det()) == 1
historical_in_equation = block_diagonal_matrix(
    identity_matrix(ZZ, 2),
    matrix(ZZ, fingerprint["selected"]["frame_isometry_historical_basis_in_equation_coordinates"]),
)
stage_key_by_label = {
    "2A5/MW7": "current_A5A5", "3A3/MW8": "current_3A3",
    "A3+2A2/MW10": "current_A3_2A2", "5A1/MW12": "current_5A1",
    "4A1/MW13": "current_4A1", "3A1/MW14": "current_3A1",
    "2A1/MW15": "current_2A1", "A1/MW16": "current_A1",
    "rootless/MW17": "pinned_R17",
}
cumulative = identity_matrix(ZZ, 19)
suffix_bases_equation = {"current_A11": historical_in_equation}
for index, step in enumerate(manifest["forward_steps"]):
    if index < 2:
        continue
    cumulative = matrix(ZZ, step["transition"]) * cumulative
    suffix_bases_equation[stage_key_by_label[step["child"]]] = cumulative * historical_in_equation
suffix_bases_physical = {
    name: basis * equation_in_physical for name, basis in suffix_bases_equation.items()
}
assert all(abs(basis.det()) == 1 for basis in suffix_bases_physical.values())

targets_equation = {
    record["target"]: vector(ZZ, record["target_fibre_in_state"])
    for record in crossovers["records"] if record["state"] == "equation_A11"
}
targets_physical = {
    name: value * equation_in_physical for name, value in targets_equation.items()
}
for name, basis in suffix_bases_physical.items():
    targets_physical[name] = vector(ZZ, basis.row(0))
assert targets_physical["current_A5A5"] == old_fibre
assert all(value * physical_gram * value == 0 for value in targets_physical.values())

# The historical current-3A3 target above is the raw q104 chamber class.  Add
# the physically Weyl-reduced q10 landing and compose the unchanged suffix in
# that exact marked basis.  Keeping both names prevents another accidental
# comparison against the nonphysical raw class.
q10_basis_abstract = matrix(ZZ, q10["landing"]["parent_to_current_3A3_basis"])
q10_basis_physical = q10_basis_abstract * change.inverse().change_ring(ZZ)
assert abs(q10_basis_physical.det()) == 1
assert vector(ZZ, q10_basis_physical.row(0)) == vector(
    ZZ, q10["physical_weyl_repair"]["repaired_fibre"]
) * change.inverse().change_ring(ZZ)
physical_suffix_bases = {"physical_q10_current_3A3": q10_basis_physical}
suffix_cumulative = identity_matrix(ZZ, 19)
for step in manifest["forward_steps"][4:]:
    suffix_cumulative = matrix(ZZ, step["transition"]) * suffix_cumulative
    key = "physical_q10_" + stage_key_by_label[step["child"]]
    physical_suffix_bases[key] = suffix_cumulative * q10_basis_physical
assert all(abs(basis.det()) == 1 for basis in physical_suffix_bases.values())
for name, basis in physical_suffix_bases.items():
    suffix_bases_physical[name] = basis
    targets_physical[name] = vector(ZZ, basis.row(0))
assert targets_physical["physical_q10_current_3A3"][0:2] == vector(ZZ, [5, 2])
assert all(value * physical_gram * value == 0 for value in targets_physical.values())

def to_physical(value):
    return vector(ZZ, value) * change.inverse().change_ring(ZZ)


# Only curves already explicit over QQ are included.  The twelve fibre
# components are equation-marked; C9 and O are exact sections; the remaining
# section classes have exact polynomial representatives in the named inputs.
explicit_abstract = {
    "old_zero": old_zero,
    **{f"old_A11_component_{index}": curve for index, curve in components.items()},
    "first_I6_affine_component": affines[0],
    "second_I6_affine_component": affines[1],
    "old_A11_affine": vector(ZZ, physical["old_A11_affine_section_on_component9_pointed_child"][
        "NS_coordinates_in_selected_child_basis"
    ]),
    "oldI9_A0": vector(
        ZZ, mismatch["correct_selected_R3_transport"]["oldI9_A0"]["child_coordinates"]
    ) * abstract_in_equation.inverse().change_ring(ZZ),
    "close_P24": vector(
        ZZ, mismatch["correct_selected_R3_transport"]["close_P24"]["child_coordinates"]
    ) * abstract_in_equation.inverse().change_ring(ZZ),
    "P230": vector(ZZ, p230_word["q4_orbit230_horizontal"]["effective_section"]),
    "P1229": vector(ZZ, p1229["P1229"]["NS_coordinates"]),
    "P146": vector(ZZ, p146["P146"]["NS_coordinates"]),
    "P1307": vector(ZZ, p146["P1307"]["NS_coordinates"]),
    "P1": vector(ZZ, p146["polynomial_inputs"]["P1"]["NS_coordinates"]),
    "P32": vector(ZZ, p146["polynomial_inputs"]["P32"]["NS_coordinates"]),
}
assert all(curve * abstract_gram * curve == -2 for curve in explicit_abstract.values())
explicit_physical = {name: to_physical(curve) for name, curve in explicit_abstract.items()}
assert all(curve * physical_gram * curve == -2 for curve in explicit_physical.values())
assert all(to_physical(curve) == -vector(ZZ, [0, 0] + [1 if j == i else 0 for j in range(17)])
           for i, curve in enumerate(effective_simple))

inputs = (
    ZERO, PHYSICAL, MANIFEST, FINGERPRINT, CROSSOVERS, MISMATCH,
    P230_WORD, P1229, P146, Q10, Q10_EFFECTIVE, abstract_frame_path,
)
payload = {
    "schema": "elkies-k3.h3-a5a5-physical-component-chamber-marking.v1",
    "status": "PASS_EXACT_A5A5_PHYSICAL_COMPONENT_CHAMBER_MARKING",
    "hub": "equation_effective_component9_zero_physical_component_chamber_2A5",
    "root_data": selected["root_data"],
    "physical_simple_roots": effective_simple_names,
    "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "frame_sha256": sha256(FRAME_OUTPUT),
    "physical_component_chamber_basis_in_equation_A11": rows(physical_in_equation),
    "equation_A11_basis_in_physical_component_chamber": rows(equation_in_physical),
    "equation_A11_to_root_adapted_hub_basis": rows(physical_in_equation),
    "root_adapted_hub_to_equation_A11_basis": rows(equation_in_physical),
    "target_fibres_in_root_adapted_hub": {
        name: entries(value) for name, value in targets_physical.items()
    },
    "current_suffix_stage_bases_in_root_adapted_hub": {
        name: rows(basis) for name, basis in suffix_bases_physical.items()
    },
    "equation_explicit_curves_in_child": {
        name: entries(value) for name, value in explicit_physical.items()
    },
    "known_exact_QQ_sections": [
        "old_zero", "old_A11_component_9", "old_A11_affine", "P230",
        "P1229", "P146", "P1307", "P1", "P32",
    ],
    "prefix_operational_score": 0,
    "proof_boundary": (
        "Exact primitive physical A5+A5 component basis, full equation-A11 transport, "
        "named reverse targets and suffix bases, and already-explicit QQ curves. Outgoing "
        "neighbours and equation-cost claims require separate exhaustive artifacts."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in inputs},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A5PHYSICALCHAMBER|root=A5+A5|explicit={}|sections={}|targets={}|"
    "current3A3_degree={}|pinned_degree={}|status={}|output={}".format(
        len(explicit_physical), len(payload["known_exact_QQ_sections"]), len(targets_physical),
        targets_physical["current_3A3"][1], targets_physical["pinned_R17"][1],
        payload["status"], OUTPUT,
    ), flush=True,
)
