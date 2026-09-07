#!/usr/bin/env sage -python
"""Reframe the physical q10 pencil with an equation-effective section zero."""

import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
Q10 = LOCAL / "q24-2a5-direct-physical-q10-certificate.json"
MARKING = LOCAL / "q24-a11-to-2a5-q8-equation-marking-qq.json"
ZERO_FRAME = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
CURRENT_3A3 = GENERATED / "elkies-k3-h3-current_3A3-marked-frame.json"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = LOCAL / "q24-2a5-direct-physical-q10-effective-c5-zero-certificate.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(item) for item in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


q10 = json.loads(Q10.read_text())
marking = json.loads(MARKING.read_text())
zero_frame = json.loads(ZERO_FRAME.read_text())
current = json.loads(CURRENT_3A3.read_text())
assert q10["status"] == "PASS_EXACT_PHYSICAL_NEF_Q10_CURRENT_3A3_PRESENTATION"
assert marking["status"] == "PASS_EXACT_Q24_A11_Q8_2A5_EQUATION_MARKING"
assert zero_frame["status"] == "PASS_EXACT_A11_Q8_ORBIT12_EXPLICIT_ZERO_FRAMES"
assert current["status"] == "PASS_EXACT_CURRENT_SUFFIX_STAGE_MARKING"

frame = matrix(ZZ, zero_frame["selected"]["frame"])
gram = block_diagonal_matrix(U2, -frame)
fibre = vector(ZZ, q10["physical_weyl_repair"]["repaired_fibre"])
old_fibre = vector(ZZ, [1, 0] + [0] * 17)
components = {
    index: vector(ZZ, marking["physical_2A5"]["child_coordinates"][
        f"old_A11_component_{index}"
    ])
    for index in range(11)
}
affines = [
    old_fibre - sum(
        (components[index] for index in chain), vector(ZZ, [0] * 19)
    )
    for chain in marking["physical_2A5"]["chains"]
]
physical_curves = {
    **{f"old_A11_component_{index}": curve for index, curve in components.items()},
    "first_I6_affine_component": affines[0],
    "second_I6_affine_component": affines[1],
}

# The historical/canonical 3A3 basis still carries a chamber pseudo-zero.
canonical_parent_to_child = matrix(
    ZZ, q10["landing"]["parent_to_current_3A3_basis"]
)
canonical_zero = vector(
    ZZ, canonical_parent_to_child.row(1) - canonical_parent_to_child.row(0)
)
canonical_zero_pairings = {
    name: int(canonical_zero * gram * curve)
    for name, curve in physical_curves.items()
}
assert canonical_zero * gram * canonical_zero == -2
assert canonical_zero * gram * fibre == 1
assert [name for name, value in canonical_zero_pairings.items() if value < 0] == [
    "old_A11_component_0", "old_A11_component_2", "old_A11_component_4",
    "old_A11_component_7", "second_I6_affine_component",
]

degree_one = {
    name: curve for name, curve in physical_curves.items()
    if curve * gram * fibre == 1
}
assert list(degree_one) == [
    "old_A11_component_1", "old_A11_component_5",
    "old_A11_component_6", "first_I6_affine_component",
]


def child_basis(section):
    mate = fibre + section
    kernel = matrix(ZZ, [list(fibre * gram), list(mate * gram)]).right_kernel_matrix()
    basis = matrix(ZZ, [list(fibre), list(mate)] + [list(row) for row in kernel.rows()])
    assert abs(basis.det()) == 1
    child_frame = -(kernel * gram * kernel.transpose())
    assert basis * gram * basis.transpose() == block_diagonal_matrix(U2, -child_frame)
    return basis, child_frame


candidates = []
for name, section in degree_one.items():
    basis, child_frame = child_basis(section)
    root_result = pari(child_frame).qfminim(2)
    half_roots = matrix(ZZ, root_result[2]).transpose()
    root_module = half_roots.row_module()
    root_basis = root_module.basis_matrix()
    root_gram = root_basis * child_frame * root_basis.transpose()
    root_data = [int(root_module.rank()), int(root_result[0]), int(root_gram.det())]
    assert root_data == [9, 36, 64]
    candidates.append({
        "zero": name,
        "section": section,
        "basis": basis,
        "child_frame": child_frame,
        "root_basis": root_basis,
        "root_gram": root_gram,
        "root_data": root_data,
        "frame_coordinate_growth": int(max(abs(item) for item in child_frame.list())),
    })

candidates.sort(key=lambda item: (item["frame_coordinate_growth"], item["zero"]))
selected = candidates[0]
assert selected["zero"] == "old_A11_component_5"
assert [(item["zero"], item["frame_coordinate_growth"]) for item in candidates] == [
    ("old_A11_component_5", 58),
    ("old_A11_component_6", 62),
    ("old_A11_component_1", 146),
    ("first_I6_affine_component", 156),
]

parent_to_child = selected["basis"]
child_to_parent = parent_to_child.inverse().change_ring(ZZ)
child_gram = block_diagonal_matrix(U2, -selected["child_frame"])
assert vector(ZZ, parent_to_child.row(1) - parent_to_child.row(0)) == degree_one[selected["zero"]]

# Identify the effective-zero child with both canonical current 3A3 and pinned
# R17 by complete bases, independently of ADE/MW labels.
current_frame = load_matrix(ROOT / current["frame_output"])
current_gram = block_diagonal_matrix(U2, -current_frame)
pinned_frame = load_matrix(PINNED_FRAME)
pinned_gram = block_diagonal_matrix(U2, -pinned_frame)
pinned_in_current = matrix(ZZ, current["pinned_R17_basis_in_source"])
pinned_in_parent = pinned_in_current * canonical_parent_to_child
pinned_in_child = pinned_in_parent * child_to_parent
assert abs(pinned_in_child.det()) == 1
assert pinned_in_child * child_gram * pinned_in_child.transpose() == pinned_gram
current_in_child = pinned_in_current.inverse().change_ring(ZZ) * pinned_in_child
child_in_current = current_in_child.inverse().change_ring(ZZ)
assert abs(current_in_child.det()) == 1
assert current_in_child * child_gram * current_in_child.transpose() == current_gram

inputs = (Q10, MARKING, ZERO_FRAME, CURRENT_3A3, ROOT / current["frame_output"], PINNED_FRAME)
payload = {
    "schema": "elkies-k3.q24-2a5-direct-physical-q10-effective-zero.v1",
    "status": "PASS_EXACT_PHYSICAL_Q10_EFFECTIVE_C5_ZERO_TO_PINNED_R17",
    "canonical_zero_audit": {
        "zero": entries(canonical_zero),
        "physical_curve_pairings": canonical_zero_pairings,
        "negative_physical_curves": [
            name for name, value in canonical_zero_pairings.items() if value < 0
        ],
        "equation_effective": False,
    },
    "effective_degree_one_candidates": [
        {
            "zero": item["zero"],
            "section": entries(item["section"]),
            "frame_coordinate_growth": item["frame_coordinate_growth"],
            "root_data": item["root_data"],
        }
        for item in candidates
    ],
    "selection": {
        "zero": selected["zero"],
        "reason": "minimum exact child-frame coefficient growth among the four physical degree-one curves",
        "section": entries(selected["section"]),
        "marked_U": {
            "fibre_in_parent": entries(fibre),
            "mate_in_parent": entries(parent_to_child.row(1)),
            "zero_in_parent": entries(parent_to_child.row(1) - parent_to_child.row(0)),
        },
        "child_root_data": selected["root_data"],
        "child_root_basis": rows(selected["root_basis"]),
        "child_root_gram": rows(selected["root_gram"]),
        "child_frame": rows(selected["child_frame"]),
    },
    "transport": {
        "parent_to_effective_zero_child_basis": rows(parent_to_child),
        "effective_zero_child_to_parent_basis": rows(child_to_parent),
        "forward_determinant": int(parent_to_child.det()),
        "inverse_determinant": int(child_to_parent.det()),
        "Gram_transport_exact": True,
    },
    "current_3A3_identification": {
        "current_3A3_basis_in_effective_zero_child": rows(current_in_child),
        "effective_zero_child_basis_in_current_3A3": rows(child_in_current),
        "forward_determinant": int(current_in_child.det()),
        "inverse_determinant": int(child_in_current.det()),
        "Gram_transport_exact": True,
    },
    "endpoint": {
        "name": "pinned_R17",
        "canonical_pinned_basis_in_effective_zero_child": rows(pinned_in_child),
        "effective_zero_child_basis_in_canonical_pinned": rows(
            pinned_in_child.inverse().change_ring(ZZ)
        ),
        "forward_determinant": int(pinned_in_child.det()),
        "inverse_determinant": int(pinned_in_child.inverse().det()),
        "Gram_transport_exact": True,
    },
    "proof_boundary": (
        "Exact rejection of the canonical chamber pseudo-zero, exhaustive audit of the four "
        "known physical degree-one curves, exact marked U with effective component C5 as zero, "
        "3A3 root data, bidirectional unimodular NS transports, and full pinned-R17 endpoint "
        "identification. The q10 characteristic-zero pencil and equation marking remain to be compiled."
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
    "Q10EFFECTIVEZERO|canonical_effective=0|candidates=4|selected=C5|frame_max=58|"
    "landing_det={}|endpoint_det={}|status={}|output={}".format(
        current_in_child.det(), pinned_in_child.det(), payload["status"], OUTPUT,
    ),
    flush=True,
)
