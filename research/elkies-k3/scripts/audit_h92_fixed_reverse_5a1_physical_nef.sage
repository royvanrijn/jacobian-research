#!/usr/bin/env sage
"""Audit the physical-nef fixed reverse q4/orbit52 fibre on 4A1.

The input 4A1 marking is the exact four-reflection q114 marking, not the raw
abstract corridor frame.  Reduce q52 against all four physical I2 affine
cycles, then apply the all-section closest-vector and complete finite
horizontal-wall gates of Proposition C2.  No equation algebra is performed.
"""

import hashlib
import json
from pathlib import Path

from sage.all import (
    IntegralLattice,
    QQ,
    ZZ,
    block_diagonal_matrix,
    block_matrix,
    identity_matrix,
    matrix,
    pari,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
CURRENT = LOCAL / "fixed-reverse-4a1-physical-nef-audit.json"
POINTING = LOCAL / "fixed-reverse-4a1-pointing-qq.json"
OUTPUT = LOCAL / "fixed-reverse-5a1-physical-nef-audit.json"


def read_json(path):
    return json.loads(path.read_text())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def values(value):
    return [int(entry) for entry in value]


def pairing(left, gram, right):
    return ZZ(left * gram * right)


manifest = read_json(MANIFEST)
current = read_json(CURRENT)
pointing = read_json(POINTING)
assert current["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_PHYSICAL_NEF"
assert pointing["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_POINTING"

U2 = matrix(ZZ, [[0, 1], [1, 0]])
G4 = matrix(ZZ, current["full_marked_transport"]["physical_4A1_positive_frame"])
Q4 = block_diagonal_matrix(U2, -G4)
assert G4[:4, :4] == 2 * identity_matrix(ZZ, 4)
E19 = identity_matrix(ZZ, 19)
old_fibre = vector(ZZ, E19.row(0))

# The pointing gate contracts the reflected equation curves.  In the exact
# physical 4A1 basis their effective signs are e2,-e3,-e4,e5.
simple_components = [
    vector(ZZ, E19.row(2)),
    -vector(ZZ, E19.row(3)),
    -vector(ZZ, E19.row(4)),
    vector(ZZ, E19.row(5)),
]
affine_components = [old_fibre - component for component in simple_components]
component_walls = simple_components + affine_components
component_names = [
    "old_I2_1_nonidentity",
    "old_I2_2_nonidentity",
    "old_I2_3_nonidentity",
    "old_I2_4_nonidentity",
    "old_I2_1_affine",
    "old_I2_2_affine",
    "old_I2_3_affine",
    "old_I2_4_affine",
]
assert all(pairing(curve, Q4, curve) == -2 for curve in component_walls)

step = manifest["forward_steps"][-5]
assert (step["parent"], step["child"], step["q"], step["orbit"]) == (
    "5A1/MW12", "4A1/MW13", 4, 52,
)
T52 = matrix(ZZ, step["transition"])
raw_parent_basis = T52.inverse().change_ring(ZZ)
raw_fibre = vector(ZZ, raw_parent_basis.row(0))
assert pairing(raw_fibre, Q4, raw_fibre) == 0 and raw_fibre[1] == 2


def reflection_matrix(root):
    action = identity_matrix(ZZ, 19) + (Q4 * root.column()) * matrix(ZZ, [list(root)])
    assert action * Q4 * action.transpose() == Q4
    assert action.det() == -1
    return action


def reduce_components(fibre):
    value = vector(ZZ, fibre)
    action = identity_matrix(ZZ, 19)
    record = []
    for unused in range(10000):
        negative = next((
            (name, curve, pairing(value, Q4, curve))
            for name, curve in zip(component_names, component_walls)
            if pairing(value, Q4, curve) < 0
        ), None)
        if negative is None:
            return value, action, record
        name, curve, degree = negative
        step_action = reflection_matrix(curve)
        value *= step_action
        action *= step_action
        record.append({"wall": name, "pairing": int(degree)})
    raise RuntimeError("q52 physical component reduction did not terminate")


def all_section_minimum(fibre):
    degree = ZZ(fibre[1])
    center = vector(QQ, fibre[2:]) / degree
    closest = vector(ZZ, next(IntegralLattice(G4).enumerate_close_vectors(center)))
    distance = (closest - center) * G4 * (closest - center)
    return distance, degree * (distance - 2) / 2


def negative_horizontal_walls(fibre):
    degree = ZZ(fibre[1])
    w = vector(ZZ, fibre[2:])
    walls = []
    for old_degree in range(1, int(degree) + 1):
        m = ZZ(old_degree)
        cross = -degree * m * G4 * w.column()
        augmented = block_matrix(ZZ, [
            [degree ** 2 * G4, cross],
            [cross.transpose(), matrix(ZZ, [[m ** 2 * (w * G4 * w) + 1]])],
        ])
        result = pari(augmented).qfminim(2 * degree ** 2 - 1)
        normalized = set()
        for raw in matrix(ZZ, result[2]).transpose().rows():
            if abs(raw[-1]) != 1:
                continue
            value = raw if raw[-1] == 1 else -raw
            normalized.add(tuple(value))
        for value in normalized:
            x = vector(ZZ, value[:-1])
            x_norm = ZZ(x * G4 * x)
            if (x_norm - 2) % (2 * m):
                continue
            k = ZZ((x_norm - 2) // (2 * m))
            intersection = ZZ(
                (w * G4 * w // (2 * degree)) * m
                + degree * k - w * G4 * x
            )
            if intersection < 0:
                walls.append({
                    "old_fibre_degree": int(m),
                    "curve": [int(k), int(m)] + values(x),
                    "intersection": int(intersection),
                })
    return sorted(walls, key=lambda item: (item["intersection"], item["curve"]))


physical_fibre, component_action, reflection_record = reduce_components(raw_fibre)
component_degrees = {
    name: int(pairing(physical_fibre, Q4, curve))
    for name, curve in zip(component_names, component_walls)
}
distance, section_minimum = all_section_minimum(physical_fibre)
horizontal_walls = negative_horizontal_walls(physical_fibre)

physical_parent_basis = raw_parent_basis * component_action
current_to_physical_parent = physical_parent_basis.inverse().change_ring(ZZ)
parent_gram = physical_parent_basis * Q4 * physical_parent_basis.transpose()
assert parent_gram[:2, :2] == U2 and not parent_gram[:2, 2:]
physical_zero = vector(
    ZZ, physical_parent_basis.row(1) - physical_parent_basis.row(0)
)
physical_roots = [
    vector(ZZ, physical_parent_basis.row(index))
    for index in range(2, 7)
]
assert all(pairing(root, Q4, physical_fibre) == 0 for root in physical_roots)
assert all(pairing(root, Q4, physical_zero) == 0 for root in physical_roots)
assert abs(physical_parent_basis.det()) == 1

status = (
    "PASS_EXACT_QQ_FIXED_REVERSE_5A1_PHYSICAL_NEF"
    if min(component_degrees.values()) >= 0
    and section_minimum >= 0
    and not horizontal_walls
    else "FAIL_PHYSICAL_NEF_AUDIT"
)
payload = {
    "schema": "elkies-k3.fixed-reverse-5a1-physical-nef-audit.v1",
    "status": status,
    "edge": "5A1/MW12 --q4 orbit52--> 4A1/MW13, audited in reverse",
    "raw_fibre_in_4A1_coordinates": values(raw_fibre),
    "raw_component_degrees": {
        name: int(pairing(raw_fibre, Q4, curve))
        for name, curve in zip(component_names, component_walls)
    },
    "physical_component_reduction": {
        "reflection_sequence": reflection_record,
        "action": rows(component_action),
        "determinant": int(component_action.det()),
        "reduced_fibre_in_4A1_coordinates": values(physical_fibre),
        "old_fibre_degree": int(physical_fibre[1]),
        "q": int(physical_fibre[0] * physical_fibre[1]),
        "isotropic": pairing(physical_fibre, Q4, physical_fibre) == 0,
        "component_degrees": component_degrees,
    },
    "all_section_gate": {
        "closest_vector_distance": str(distance),
        "minimum_section_intersection": str(section_minimum),
        "pass": bool(section_minimum >= 0),
    },
    "finite_horizontal_wall_gate": {
        "proposition": "RANK_MUTATION_AND_LIFT_THEOREMS.md Proposition C2",
        "negative_walls": horizontal_walls,
        "pass": not horizontal_walls,
    },
    "full_marked_transport": {
        "physical_5A1_basis_in_4A1_coordinates": rows(physical_parent_basis),
        "4A1_to_physical_5A1_basis": rows(current_to_physical_parent),
        "basis_determinant": int(physical_parent_basis.det()),
        "mutually_inverse": True,
        "physical_zero_in_4A1_coordinates": values(physical_zero),
        "physical_5A1_root_classes_in_4A1_coordinates": [
            values(root) for root in physical_roots
        ],
        "physical_5A1_positive_frame": rows(-parent_gram[2:, 2:]),
    },
    "proof_boundary": (
        "Exact physical component reduction, all-section closest-vector gate, "
        "complete finite horizontal-wall enumeration, and full marked NS "
        "transport. The horizontal section, equation and pointing are separate."
    ),
    "inputs": {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (MANIFEST, CURRENT, POINTING)
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDREVERSE5A1NEF|status={}|reflections={}|q={}|section_min={}|walls={}|output={}".format(
        status, len(reflection_record),
        payload["physical_component_reduction"]["q"], section_minimum,
        len(horizontal_walls), OUTPUT,
    )
)
