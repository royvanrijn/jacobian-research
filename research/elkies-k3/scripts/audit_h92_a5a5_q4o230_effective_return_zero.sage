#!/usr/bin/env sage-python
"""Audit the effective q4/o230 return zero against the stored route chamber."""

import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
RETURN = LOCAL / "q24-a1a4a5-to-2a5-q4-return-resolved-rr-qq.json"
MARKING = LOCAL / "q24-2a5-to-a1a4a5-q4o230-equation-marking-qq.json"
PRESENTATIONS = GENERATED / (
    "elkies-k3-h3-a5a5-q4o230-c10-second-zero-changing-3a3-presentations.json"
)
PARENT_FRAME = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frame.txt"
OUTPUT = GENERATED / "elkies-k3-h3-a5a5-q4o230-effective-return-zero-audit.json"
INPUTS = (RETURN, MARKING, PRESENTATIONS, PARENT_FRAME)


def load_frame(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


returned = json.loads(RETURN.read_text())
marking = json.loads(MARKING.read_text())
presentations = json.loads(PRESENTATIONS.read_text())
assert returned["status"] == "PASS_EXACT_Q24_A1A4A5_Q4_RETURN_2A5_CHANGED_ZERO"
assert marking["status"] == "PASS_EXACT_Q24_2A5_Q4O230_COMPONENT10_EQUATION_MARKING"

selected = next(
    row for row in presentations["ranked_presentations"]
    if row["first_edge_candidate_id"] == {
        "q": 6, "old_fibre_degree": 2, "orbit_index": 1315,
    }
    and row["explicit_zero_curve"] == "old_A5A5_component_1"
)

u2 = matrix(ZZ, ((0, 1), (1, 0)))
parent_frame = load_frame(PARENT_FRAME)
g_parent = block_diagonal_matrix(u2, -parent_frame)
stored_to_equation = matrix(
    ZZ, returned["transport"]["stored_promoted_return_to_equation_adapted_return_basis"]
)
equation_to_child = matrix(
    ZZ, returned["transport"]["equation_adapted_return_to_component10_child_basis"]
)
child_to_parent = matrix(
    ZZ, marking["transport"]["component10_zero_child_to_component9_zero_parent_basis"]
)
stored_to_parent = stored_to_equation * equation_to_child * child_to_parent

fibre = vector(ZZ, [1, 0] + [0] * 17)
zero = vector(ZZ, [-1, 1] + [0] * 17)
pseudo_zero_parent = zero * stored_to_parent
effective_p230_parent = zero * equation_to_child * child_to_parent
q6_fibre_stored = vector(ZZ, selected["source_to_explicit_child_basis"][0])
q6_fibre_parent = q6_fibre_stored * stored_to_parent
vertical_difference = pseudo_zero_parent - effective_p230_parent

assert entries(pseudo_zero_parent) == [
    27, 1, -2, 3, 3, -8, -9, 1, -8, -8, -6, -4, -2, 0, 1, -1, -1, -1, -1,
]
assert pseudo_zero_parent * g_parent * zero == 26
assert effective_p230_parent * g_parent * zero == 2
assert effective_p230_parent * g_parent * fibre == 1
# The parent root rank is ten: fibre plus the first ten complement coordinates
# generate the vertical subgroup used here.
assert vertical_difference[1] == 0
assert all(value == 0 for value in vertical_difference[12:])
assert q6_fibre_parent * g_parent * q6_fibre_parent == 0
assert q6_fibre_parent * g_parent * fibre == 2
assert q6_fibre_parent * g_parent * effective_p230_parent == 54
assert q6_fibre_parent * g_parent * zero == 58

payload = {
    "schema": "elkies-k3.h3-a5a5-q4o230-effective-return-zero-audit.v1",
    "status": "PASS_EXACT_Q4O230_EFFECTIVE_ZERO_INVALIDATES_Q6O1315_COST",
    "q4_return": {
        "equation_exact": True,
        "effective_changed_zero": "forward I2 nonidentity component (P230 branch)",
        "effective_P230_in_parent_2A5": entries(effective_p230_parent),
        "effective_P230_dot_original_zero": 2,
    },
    "stored_route_chamber": {
        "zero_in_parent_2A5": entries(pseudo_zero_parent),
        "zero_dot_original_zero": 26,
        "difference_from_effective_P230": entries(vertical_difference),
        "difference_is_vertical_in_old_fibre_and_2A5_roots": True,
        "same_MW_quotient_class_as_P230": True,
        "interpretation": "Weyl/chamber pseudo-zero, not the effective equation section",
    },
    "q6_orbit1315": {
        "fibre_in_stored_returned_2A5": entries(q6_fibre_stored),
        "fibre_in_parent_2A5": entries(q6_fibre_parent),
        "primitive_isotropic_and_old_fibre_degree_two": True,
        "intersection_with_effective_P230": 54,
        "intersection_with_original_zero": 58,
    },
    "cost_consequence": {
        "withdraw_operational_score": 4199,
        "withdraw_q6_orbit1315_as_lifting_target": True,
        "lattice_path_remains_exact": True,
        "q4_return_equation_remains_exact": True,
        "required_next_step": "rerank returned frontier using equation-effective curves only",
    },
    "proof_boundary": (
        "Exact integral transport and intersection audit. This invalidates the former "
        "equation-cost interpretation, not the primitive-nef lattice path. It does not "
        "select a replacement suffix."
    ),
    "reproduce_command": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/audit_h92_a5a5_q4o230_effective_return_zero.sage"
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
    "Q4O230ZERO|pseudo_PO=26|effective_PO=2|q6_P230=54|q6_O=58|"
    "withdraw=4199|status={}|output={}".format(payload["status"], OUTPUT),
    flush=True,
)
