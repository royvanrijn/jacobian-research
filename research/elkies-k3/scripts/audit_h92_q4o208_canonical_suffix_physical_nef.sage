#!/usr/bin/env sage -python
"""Audit the inherited canonical 3A3 suffix in the physical q4/o208 chamber.

The abstract q4/o208 landing is identified with the stored canonical 3A3
frame by an integral isometry.  Such an isometry need not preserve the nef
chamber.  Pull the first canonical suffix fibre back to the equation-effective
C5 frame and test it against every inherited physical curve.  It crosses the
wall of old_zero=old_A11_component_9 once.  Reflecting in that effective
(-2)-curve gives the physical-chamber representative.
"""

import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
ROUTE = GENERATED / "elkies-k3-h3-a5a5-physical-q4o208-to-pinned-r17-certificate.json"
MARKING = GENERATED / "elkies-k3-h3-a5a5-physical-component-chamber-marking.json"
EQUATION_MARKING = LOCAL / "q24-2a5-physical-q4o208-equation-marking-qq.json"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
CURRENT = GENERATED / "elkies-k3-h3-current_3A3-marked-frame.json"
OUTPUT = GENERATED / "elkies-k3-h3-q4o208-canonical-suffix-physical-nef-audit.json"
INPUTS = (ROUTE, MARKING, EQUATION_MARKING, MANIFEST, CURRENT)

route = json.loads(ROUTE.read_text())
marking = json.loads(MARKING.read_text())
equation_marking = json.loads(EQUATION_MARKING.read_text())
manifest = json.loads(MANIFEST.read_text())
current = json.loads(CURRENT.read_text())
assert route["status"] == "PASS_EXACT_PHYSICAL_Q4O208_3A3_TO_PINNED_R17"
assert marking["status"] == "PASS_EXACT_A5A5_PHYSICAL_COMPONENT_CHAMBER_MARKING"
assert equation_marking["status"] == "PASS_EXACT_QQ_PHYSICAL_Q4O208_C5_EQUATION_MARKING"
assert manifest["status"] == "PASS_H3_R17_BACKWARD_EXACT_LIFT_MANIFEST"
assert current["status"] == "PASS_EXACT_CURRENT_SUFFIX_STAGE_MARKING"

U2 = matrix(ZZ, ((0, 1), (1, 0)))
g_child = block_diagonal_matrix(U2, -matrix(ZZ, route["selection"]["child_frame"]))
canonical_in_child = matrix(
    ZZ, route["transport"]["current_3A3_basis_in_effective_zero_child"]
)
parent_to_child = matrix(
    ZZ, route["transport"]["effective_zero_child_to_parent_basis"]
)
current_frame = matrix(ZZ, [
    [ZZ(value) for value in line.split()]
    for line in (ROOT / current["frame_output"]).read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
assert canonical_in_child * g_child * canonical_in_child.transpose() == (
    block_diagonal_matrix(U2, -current_frame)
)

step = manifest["forward_steps"][4]
assert step["parent"] == "3A3/MW8" and step["child"] == "A3+2A2/MW10"
assert (int(step["q"]), int(step["orbit"])) == (4, 323)
canonical_fibre = vector(ZZ, step["new_fibre_in_parent"])
canonical_horizontal = vector(ZZ, step["horizontal"]["section_class"])
canonical_zero = vector(ZZ, step["new_zero_in_parent"])
fibre = vector(ZZ, canonical_fibre * canonical_in_child)
horizontal = vector(ZZ, canonical_horizontal * canonical_in_child)
new_zero = vector(ZZ, canonical_zero * canonical_in_child)
assert fibre * g_child * fibre == 0
assert horizontal * g_child * horizontal == -2
assert horizontal * g_child * fibre == 0

curves = {
    name: vector(ZZ, value) * parent_to_child
    for name, value in marking["equation_explicit_curves_in_child"].items()
}
for fibre_name, record in equation_marking["physical_fibres"].items():
    for index, value in enumerate(record["components_in_cycle_order"]):
        curves[f"{fibre_name}_component_{index}"] = vector(ZZ, value)

# Remove aliases such as old_zero and old_A11_component_9 while retaining a
# deterministic preferred name for every physical class.
unique_curves = {}
for name, curve in sorted(curves.items()):
    unique_curves.setdefault(tuple(curve), (name, curve))
curves = {name: curve for name, curve in unique_curves.values()}


def degrees(divisor):
    return {name: int(divisor * g_child * curve) for name, curve in curves.items()}


canonical_degrees = degrees(fibre)
negative = {name: value for name, value in canonical_degrees.items() if value < 0}
assert negative == {"old_A11_component_9": -1}
wall = curves["old_A11_component_9"]
assert wall * g_child * wall == -2


def reflect(value, root):
    return vector(ZZ, value + (value * g_child * root) * root)


physical_fibre = reflect(fibre, wall)
physical_horizontal = reflect(horizontal, wall)
physical_zero = reflect(new_zero, wall)
physical_degrees = degrees(physical_fibre)
assert all(value >= 0 for value in physical_degrees.values())
assert physical_fibre * g_child * physical_fibre == 0
assert physical_horizontal * g_child * physical_fibre == 0
assert physical_horizontal[1] == 1
assert physical_fibre[1] == 2

degree_zero = sorted(name for name, value in physical_degrees.items() if value == 0)
degree_one = sorted(name for name, value in physical_degrees.items() if value == 1)
assert len(degree_zero) == 7
expected_degree_one = [
    "first_old_I6_I4_component_1",
    "first_old_I6_I4_component_2",
    "old_A11_component_1",
    "old_A11_component_10",
    "old_A11_component_2",
    "old_A11_component_9",
]
if degree_one != expected_degree_one:
    raise ArithmeticError(f"unexpected explicit degree-one curves: {degree_one}")

payload = {
    "schema": "elkies-k3.h3-q4o208-canonical-suffix-physical-nef-audit.v1",
    "status": "PASS_EXACT_Q4O208_CANONICAL_SUFFIX_PHYSICAL_WALL_CORRECTION",
    "abstract_suffix_edge": {
        "parent": step["parent"], "child": step["child"],
        "q": int(step["q"]), "orbit": int(step["orbit"]),
    },
    "canonical_pullback_to_C5_equation_frame": {
        "fibre": [int(value) for value in fibre],
        "horizontal": [int(value) for value in horizontal],
        "new_zero": [int(value) for value in new_zero],
        "negative_explicit_curve_intersections": negative,
        "nef_in_physical_equation_chamber": False,
    },
    "wall_correction": {
        "wall": "old_zero=old_A11_component_9",
        "wall_class": [int(value) for value in wall],
        "reflection_formula": "s_C(x)=x+(x.C)C for C^2=-2",
        "number_of_reflections": 1,
        "physical_fibre": [int(value) for value in physical_fibre],
        "physical_horizontal": [int(value) for value in physical_horizontal],
        "physical_new_zero": [int(value) for value in physical_zero],
        "physical_fibre_square": 0,
        "old_fibre_degree": int(physical_fibre[1]),
        "nonnegative_on_all_inherited_explicit_curves": True,
        "explicit_curve_degrees": physical_degrees,
        "explicit_degree_zero_curves": degree_zero,
        "explicit_degree_one_curves": degree_one,
    },
    "large_Groebner_required": False,
    "proof_boundary": (
        "The abstract canonical q4/orbit323 suffix fibre is not nef in the actual C5 equation "
        "chamber. One exact reflection in the inherited effective old-zero/C9 curve produces "
        "a class nonnegative on every currently explicit physical curve and preserves the "
        "abstract lattice edge. Full nefness against all effective curves, an effective child "
        "zero frame, and the characteristic-zero equation lift remain separate gates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O208SUFFIXNEF|canonical_negative=old_zero:-1|reflections=1|"
    "physical_known_nonnegative=1|degree0={}|degree1={}|status={}|output={}".format(
        len(degree_zero), len(degree_one), payload["status"], OUTPUT,
    ),
    flush=True,
)
