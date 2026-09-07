#!/usr/bin/env sage
"""Transport the fixed suffix through the physical q323 wall reflection.

status: ACTIVE_PROOF
claim: the one-wall-corrected q323 edge carries a reflected copy of the fixed suffix
inputs: physical q4/o208 marking, q323 wall audit, component-2 child marking
outputs: artifacts/local/elkies-k3/q4o323-reflected-fixed-suffix-component2-marking.json

The canonical q323 fibre is not nef in the equation-effective q4/o208
chamber.  Its physical representative is s_C(F), where C is the old-zero/C9
wall.  Applying s_C only to F destroys the next marked degree.  Applying the
same integral isometry to every later fixed-suffix fibre preserves every
intersection and is the correct marking transport to test.  This script does
only exact lattice arithmetic.
"""

import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, identity_matrix, matrix, vector


HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
if not (ROOT / "MATH_STATUS.json").is_file():
    ROOT = Path.cwd().resolve()
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
SOURCE = GENERATED / "elkies-k3-h3-q4o208-physical-3a3-marking.json"
AUDIT = GENERATED / "elkies-k3-h3-q4o208-canonical-suffix-physical-nef-audit.json"
CHILD = GENERATED / (
    "elkies-k3-h3-q4o208-corrected-a3-2a2-"
    "old_a11_component_2-marking.json"
)
OUTPUT = LOCAL / "q4o323-reflected-fixed-suffix-component2-marking.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


source = json.loads(SOURCE.read_text())
audit = json.loads(AUDIT.read_text())
child = json.loads(CHILD.read_text())
assert source["status"] == "PASS_EXACT_Q4O208_PHYSICAL_3A3_MARKING"
assert audit["status"] == "PASS_EXACT_Q4O208_CANONICAL_SUFFIX_PHYSICAL_WALL_CORRECTION"
assert child["status"] == "PASS_EXACT_CORRECTED_A3_2A2_PHYSICAL_EFFECTIVE_ZERO_MARKING"
assert child["zero"] == "old_A11_component_2"

source_frame = load_matrix(ROOT / source["frame_output"])
g_source = block_diagonal_matrix(U2, -source_frame)
# The audit records the wall in the pre-adaptation C5 child coordinates.  In
# this script's root-adapted physical 3A3 coordinates the same curve is the
# explicitly transported old_zero/C9 class.
wall = vector(ZZ, source["equation_explicit_curves_in_child"]["old_zero"])
assert wall == vector(ZZ, source["equation_explicit_curves_in_child"]["old_A11_component_9"])
assert wall * g_source * wall == -2
reflection = identity_matrix(ZZ, 19) + (g_source * wall.column()) * matrix(ZZ, [list(wall)])
assert reflection * g_source * reflection.transpose() == g_source
assert reflection.det() == -1 and reflection * reflection == identity_matrix(ZZ, 19)

source_targets = {
    name: vector(ZZ, value)
    for name, value in source["target_fibres_in_root_adapted_hub"].items()
}
reflected_source_targets = {
    name: value * reflection for name, value in source_targets.items()
}

suffix = (
    "current_A3_2A2",
    "current_5A1",
    "current_4A1",
    "current_3A1",
    "current_2A1",
    "current_A1",
    "current_rootless",
)
child_in_source = matrix(ZZ, child["basis_in_source"])
canonical_q323 = source_targets["current_A3_2A2"]
physical_q323 = child_in_source.row(0)
assert source_targets["current_A3_2A2"] == canonical_q323
assert reflected_source_targets["current_A3_2A2"] == physical_q323

adjacent = []
for name in suffix:
    value = reflected_source_targets[name]
    assert value * g_source * value == 0
for parent, target in zip(suffix, suffix[1:]):
    degree = ZZ(reflected_source_targets[parent] * g_source * reflected_source_targets[target])
    assert degree == 2
    adjacent.append({"parent": parent, "child": target, "old_fibre_degree": int(degree)})

source_in_child = matrix(ZZ, child["source_in_basis"])
assert source_in_child == child_in_source.inverse().change_ring(ZZ)
g_child = child_in_source * g_source * child_in_source.transpose()
child_frame = load_matrix(ROOT / child["frame_output"])
assert g_child == block_diagonal_matrix(U2, -child_frame)
reflected_child_targets = {
    name: value * source_in_child
    for name, value in reflected_source_targets.items()
}
assert reflected_child_targets["current_A3_2A2"] == vector(ZZ, [1] + [0] * 18)
assert reflected_child_targets["current_5A1"][1] == 2

# The reflected q207 divisor is in the correct isometry class but still has
# fixed vertical components in the component-2 equation chamber.  Replay its
# exact affine-Weyl reduction and apply the same right action to every marked
# target.  The three root components are A3 on (0,2,1), then A2 on (3,4) and
# A2 on (5,6); the displayed highest roots are checked against the Cartan
# block before use.
cartan = child_frame[:7, :7]
assert cartan == matrix(ZZ, (
    (2, 0, -1, 0, 0, 0, 0),
    (0, 2, -1, 0, 0, 0, 0),
    (-1, -1, 2, 0, 0, 0, 0),
    (0, 0, 0, 2, -1, 0, 0),
    (0, 0, 0, -1, 2, 0, 0),
    (0, 0, 0, 0, 0, 2, -1),
    (0, 0, 0, 0, 0, -1, 2),
))
component_roots = {
    "simple_{}".format(index): vector(
        ZZ, [0, 0] + [-ZZ(index == other) for other in range(17)]
    )
    for index in range(7)
}
highest = (
    (1, 1, 1, 0, 0, 0, 0),
    (0, 0, 0, 1, 1, 0, 0),
    (0, 0, 0, 0, 0, 1, 1),
)
for index, top in enumerate(highest):
    assert vector(ZZ, top) * cartan * vector(ZZ, top) == 2
    component_roots["affine_{}".format(index)] = vector(
        ZZ, [1, 0] + list(top) + [0] * 10
    )
reduction_word = (
    "simple_0", "simple_2", "simple_4", "simple_3", "affine_0",
    "simple_0", "simple_1", "simple_2", "affine_1",
)
q207_action = identity_matrix(ZZ, 19)
q207_reflections = []
for wall_name in reduction_word:
    root = component_roots[wall_name]
    value = reflected_child_targets["current_5A1"]
    pairing = ZZ(value * g_child * root)
    assert pairing < 0
    step = identity_matrix(ZZ, 19) + (g_child * root.column()) * matrix(ZZ, [list(root)])
    assert step * g_child * step.transpose() == g_child
    reflected_child_targets = {
        name: target * step for name, target in reflected_child_targets.items()
    }
    q207_action *= step
    q207_reflections.append({"wall": wall_name, "pairing": int(pairing)})
assert q207_action * g_child * q207_action.transpose() == g_child
assert reflected_child_targets["current_A3_2A2"] == vector(ZZ, [1] + [0] * 18)
q207_fibre = reflected_child_targets["current_5A1"]
assert q207_fibre[1] == 2 and q207_fibre * g_child * q207_fibre == 0
assert all(q207_fibre * g_child * root >= 0 for root in component_roots.values())

# Among the shortest lifts of the q207 MW coset, this representative has no
# vertical root residual at all.  It gives the especially simple physical
# divisor D=O+P-4F.  Record its construction degree on the q4/o208 parent;
# this is the input size for a possible fibrewise Abel recovery.
q207_section = vector(ZZ, [
    11, 1, -2, 3, 3, 4, 3, 9, 0, 3, 0, 4, 2, -4, 2, -1, 1, 1, -1,
])
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
old_fibre = vector(ZZ, [1, 0] + [0] * 17)
assert q207_section * g_child * q207_section == -2
assert q207_section * g_child * old_fibre == 1
assert q207_section * g_child * old_zero == 10
assert q207_fibre == old_zero + q207_section - 4 * old_fibre
section_in_q4o208 = q207_section * child_in_source
assert section_in_q4o208[1] == 16

root = child_frame[:7, :7]
coupling = child_frame[:7, 7:]
tail = child_frame[7:, 7:]
height = tail - coupling.transpose() * root.inverse() * coupling
section_mw = vector(ZZ, q207_section[-10:])
assert section_mw * height * section_mw == ZZ(65) / 3

# The whole remaining suffix was acted on, so all adjacent intersections stay
# two after the second physical correction as well.
for parent, target in zip(suffix, suffix[1:]):
    assert (
        reflected_child_targets[parent] * g_child * reflected_child_targets[target]
    ) == 2

payload = {
    "schema": "elkies-k3.h3-q4o323-reflected-fixed-suffix-marking.v1",
    "status": "PASS_EXACT_Q4O323_REFLECTED_FIXED_SUFFIX_MARKING",
    "hub": "q4o323_component2_zero_reflected_fixed_suffix",
    "zero": child["zero"],
    "root_data": child["root_data"],
    "ade": child["ade"],
    "frame_output": child["frame_output"],
    "frame_sha256": child["frame_sha256"],
    "basis_in_source": child["basis_in_source"],
    "source_in_basis": child["source_in_basis"],
    "equation_A11_to_root_adapted_hub_basis": child[
        "equation_A11_to_root_adapted_hub_basis"
    ],
    "physical_simple_components_in_source": child[
        "physical_simple_components_in_source"
    ],
    "equation_explicit_curves_in_child": child["equation_explicit_curves_in_child"],
    "target_fibres_in_root_adapted_hub": {
        name: entries(value) for name, value in reflected_child_targets.items()
    },
    "fixed_suffix_transport": {
        "wall": audit["wall_correction"]["wall"],
        "wall_class_in_q4o208_source": entries(wall),
        "reflection_right_action_in_q4o208_source": rows(reflection),
        "reflection_determinant": int(reflection.det()),
        "suffix": list(suffix),
        "adjacent_intersections": adjacent,
        "q323_fibre_equals_reflected_current_A3_2A2": True,
        "reflected_current_5A1_old_fibre_degree": int(
            reflected_child_targets["current_5A1"][1]
        ),
        "q207_component_reduction": {
            "reflections": q207_reflections,
            "right_action_in_q323_child": rows(q207_action),
            "action_determinant": int(q207_action.det()),
            "physical_q": int(q207_fibre[0] * q207_fibre[1]),
            "physical_fibre": entries(q207_fibre),
            "component_nef": True,
            "equation_preflight": {
                "horizontal_section": entries(q207_section),
                "horizontal_height": "65/3",
                "P_dot_O": int(10),
                "vertical_residual": [int(0)] * 7,
                "fibre_twist": int(-4),
                "identity": "D=O+P-4F",
                "smooth_chord_ambient_estimate": int(22),
                "section_in_q4o208_source": entries(section_in_q4o208),
                "q4o208_parent_degree": int(section_in_q4o208[1]),
                "q4o208_parent_a_minus_b": int(
                    section_in_q4o208[0] - section_in_q4o208[1]
                ),
            },
        },
    },
    "prefix_operational_score": child.get("prefix_operational_score"),
    "proof_boundary": (
        "Exact integral marking transport plus the nine-reflection vertical fixed-component "
        "replay for the next divisor. It proves that carrying both physical corrections "
        "through the whole fixed suffix preserves every adjacent degree-two pairing. The "
        "separate marked-candidate certificate must still prove the finite horizontal-wall "
        "nef gate and construct the 5A1 child frame."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (SOURCE, AUDIT, CHILD)],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in (SOURCE, AUDIT, CHILD)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O323REFLECTEDSUFFIX|zero={}|edges={}|next_degree={}|det={}|status={}|output={}".format(
        child["zero"], len(adjacent), reflected_child_targets["current_5A1"][1],
        reflection.det(), payload["status"], OUTPUT,
    ),
    flush=True,
)
