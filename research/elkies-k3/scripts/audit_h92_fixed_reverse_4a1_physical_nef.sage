#!/usr/bin/env sage
"""Audit the physical-nef representative of fixed reverse q4/orbit114.

This is lattice-only.  It reconstructs the exact 3A1 Neron--Severi frame,
enumerates the 77 already-constructed degree-one sections, reduces the stored
orbit114 fibre against the three physical I2 affine cycles, and applies the
complete finite old-horizontal wall test of Proposition C2.
"""

import hashlib
import itertools
import json
from pathlib import Path

from sage.all import (
    IntegralLattice,
    QQ,
    ZZ,
    block_diagonal_matrix,
    block_matrix,
    ceil,
    floor,
    identity_matrix,
    matrix,
    pari,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
BRIDGE = LOCAL / "q24-equation-d13-to-pinned-r17.json"
HORIZONTAL = LOCAL / "fixed-reverse-4a1-horizontal-from-3a1-qq.json"
OUTPUT = LOCAL / "fixed-reverse-4a1-physical-nef-audit.json"


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
bridge = read_json(BRIDGE)
horizontal = read_json(HORIZONTAL)

G1 = matrix(ZZ, bridge["final_a1_frame"])
step_2a1_to_a1 = manifest["forward_steps"][-2]
step_3a1_to_2a1 = manifest["forward_steps"][-3]
step_4a1_to_3a1 = manifest["forward_steps"][-4]
assert [step["orbit"] for step in (
    step_2a1_to_a1,
    step_3a1_to_2a1,
    step_4a1_to_3a1,
)] == [981, 498, 114]
T981 = matrix(ZZ, step_2a1_to_a1["transition"])
T498 = matrix(ZZ, step_3a1_to_2a1["transition"])
T114 = matrix(ZZ, step_4a1_to_3a1["transition"])
U2 = matrix(ZZ, [[0, 1], [1, 0]])
Q1 = block_diagonal_matrix(U2, -G1)
Q2 = T981.inverse().change_ring(ZZ) * Q1 * T981.inverse().change_ring(ZZ).transpose()
Q3 = T498.inverse().change_ring(ZZ) * Q2 * T498.inverse().change_ring(ZZ).transpose()
assert Q3[:2, :2] == U2 and not Q3[:2, 2:]
G3 = (-Q3[2:, 2:]).change_ring(ZZ)
assert G3[:3, :3] == 2 * identity_matrix(ZZ, 3)

# Reproduce the exact 77 degree-one section classes used by the horizontal
# constructor, without evaluating any elliptic group law.
G2 = (-Q2[2:, 2:]).change_ring(ZZ)
R2 = G2[:2, :2]
coupling = G2[:2, 2:]
H2 = G2[2:, 2:] - coupling.transpose() * R2.inverse() * coupling
H2_twice = (2 * H2).change_ring(ZZ)
sections = []
for shell in IntegralLattice(H2_twice).short_vectors(13):
    for tail_value in shell:
        tail = vector(ZZ, tail_value)
        raw_tail = vector(ZZ, [0, 0] + list(tail))
        dual = vector(QQ, (raw_tail * G2)[:2]) * R2.inverse()
        root_choices = [
            sorted({ZZ(floor(-entry)), ZZ(ceil(-entry))})
            for entry in dual
        ]
        for root_coordinates in itertools.product(*root_choices):
            frame = vector(ZZ, list(root_coordinates) + list(tail))
            norm = ZZ(frame * G2 * frame)
            if norm < 4 or (norm - 4) % 2:
                continue
            section_2a1 = vector(ZZ, [(norm - 4) // 2 + 1, 1] + list(frame))
            section_3a1 = vector(ZZ, section_2a1 * T498)
            if section_3a1[1] == 1:
                sections.append(section_3a1)
assert len(sections) == 77
assert all(pairing(section, Q3, section) == -2 for section in sections)

E19 = identity_matrix(ZZ, 19)
old_fibre = vector(ZZ, E19.row(0))
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
simple_components = [
    -vector(ZZ, E19.row(index))
    for index in (2, 3, 4)
]
affine_components = [old_fibre - component for component in simple_components]
component_walls = simple_components + affine_components
component_names = [
    "old_I2_1_nonidentity",
    "old_I2_2_nonidentity",
    "old_I2_3_nonidentity",
    "old_I2_1_affine",
    "old_I2_2_affine",
    "old_I2_3_affine",
]
assert all(pairing(curve, Q3, curve) == -2 for curve in component_walls)

raw_fibre = vector(ZZ, E19.row(0) * T114.inverse().change_ring(ZZ))
assert values(raw_fibre) == horizontal["fixed_edge"]["reverse_fibre_in_3A1_coordinates"]
assert pairing(raw_fibre, Q3, raw_fibre) == 0


def reflection_matrix(root):
    action = identity_matrix(ZZ, 19) + (Q3 * root.column()) * matrix(ZZ, [list(root)])
    assert action * Q3 * action.transpose() == Q3
    assert action.det() == -1
    return action


def reduce_components(fibre):
    value = vector(ZZ, fibre)
    action = identity_matrix(ZZ, 19)
    record = []
    for unused in range(10000):
        negative = next((
            (name, curve, pairing(value, Q3, curve))
            for name, curve in zip(component_names, component_walls)
            if pairing(value, Q3, curve) < 0
        ), None)
        if negative is None:
            return value, action, record
        name, curve, degree = negative
        step = reflection_matrix(curve)
        value *= step
        action *= step
        record.append({"wall": name, "pairing": int(degree)})
    raise RuntimeError("physical component reduction did not terminate")


def minimum_section_intersection(fibre):
    degree = ZZ(fibre[1])
    center = vector(QQ, fibre[2:]) / degree
    closest = vector(ZZ, next(IntegralLattice(G3).enumerate_close_vectors(center)))
    distance = (closest - center) * G3 * (closest - center)
    return distance, degree * (distance - 2) / 2


def negative_horizontal_walls(fibre):
    degree = ZZ(fibre[1])
    w = vector(ZZ, fibre[2:])
    walls = []
    for old_degree in range(1, int(degree) + 1):
        m = ZZ(old_degree)
        cross = -degree * m * G3 * w.column()
        augmented = block_matrix(ZZ, [
            [degree ** 2 * G3, cross],
            [cross.transpose(), matrix(ZZ, [[m ** 2 * (w * G3 * w) + 1]])],
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
            x_norm = ZZ(x * G3 * x)
            if (x_norm - 2) % (2 * m):
                continue
            k = ZZ((x_norm - 2) // (2 * m))
            intersection = ZZ(
                (w * G3 * w // (2 * degree)) * m
                + degree * k - w * G3 * x
            )
            if intersection < 0:
                walls.append({
                    "old_fibre_degree": int(m),
                    "curve": [int(k), int(m)] + values(x),
                    "intersection": int(intersection),
                })
    return sorted(walls, key=lambda item: (item["intersection"], item["curve"]))


physical_fibre, component_action, component_reflections = reduce_components(raw_fibre)
raw_parent_basis_in_3a1 = T114.inverse().change_ring(ZZ)
physical_parent_basis_in_3a1 = raw_parent_basis_in_3a1 * component_action
three_a1_to_physical_parent = physical_parent_basis_in_3a1.inverse().change_ring(ZZ)
physical_parent_gram = physical_parent_basis_in_3a1 * Q3 * physical_parent_basis_in_3a1.transpose()
assert physical_parent_gram[:2, :2] == U2 and not physical_parent_gram[:2, 2:]
physical_zero = vector(
    ZZ, physical_parent_basis_in_3a1.row(1) - physical_parent_basis_in_3a1.row(0)
)
assert physical_zero == simple_components[2]
physical_root_classes = [
    vector(ZZ, physical_parent_basis_in_3a1.row(index))
    for index in range(2, 6)
]
assert all(pairing(root, Q3, physical_fibre) == 0 for root in physical_root_classes)
assert all(pairing(root, Q3, physical_zero) == 0 for root in physical_root_classes)
assert abs(physical_parent_basis_in_3a1.det()) == 1
assert (
    physical_parent_basis_in_3a1 * three_a1_to_physical_parent
    == identity_matrix(ZZ, 19)
)
component_degrees = {
    name: int(pairing(physical_fibre, Q3, curve))
    for name, curve in zip(component_names, component_walls)
}
known_section_degrees = [int(pairing(physical_fibre, Q3, section)) for section in sections]
distance, section_minimum = minimum_section_intersection(physical_fibre)
horizontal_walls = negative_horizontal_walls(physical_fibre)

known_extra_curves = [old_zero]
known_extra_names = ["old_zero"]
raw_stored_root_degrees = {}
for index, record in enumerate(horizontal["effective_4A1_horizontal_roots_on_3A1_source"]):
    raw_curve = vector(ZZ, record["class_in_3A1_coordinates"])
    raw_stored_root_degrees[
        "stored_new_4A1_horizontal_root_{}".format(index + 1)
    ] = int(pairing(physical_fibre, Q3, raw_curve))
    known_extra_names.append("physical_new_4A1_horizontal_root_{}".format(index + 1))
    known_extra_curves.append(raw_curve * component_action)
for index, raw_values in enumerate(
    horizontal["remaining_vertical_4A1_roots"]["classes_in_3A1_coordinates"]
):
    raw_curve = vector(ZZ, raw_values)
    raw_stored_root_degrees[
        "stored_new_4A1_vertical_root_{}".format(index + 1)
    ] = int(pairing(physical_fibre, Q3, raw_curve))
    known_extra_names.append("physical_new_4A1_vertical_root_{}".format(index + 1))
    known_extra_curves.append(raw_curve * component_action)
extra_degrees = {
    name: int(pairing(physical_fibre, Q3, curve))
    for name, curve in zip(known_extra_names, known_extra_curves)
}
assert all(extra_degrees[name] == 0 for name in known_extra_names if "4A1" in name)
expected_effective_roots = [
    physical_root_classes[0],
    -physical_root_classes[1],
    -physical_root_classes[2],
    physical_root_classes[3],
]
assert {
    tuple(curve) for curve in known_extra_curves[1:]
} == {
    tuple(curve) for curve in expected_effective_roots
}

status = (
    "PASS_EXACT_QQ_FIXED_REVERSE_4A1_PHYSICAL_NEF"
    if min(component_degrees.values()) >= 0
    and min(known_section_degrees) >= 0
    and section_minimum >= 0
    and not horizontal_walls
    else "FAIL_PHYSICAL_NEF_AUDIT"
)
payload = {
    "schema": "elkies-k3.fixed-reverse-4a1-physical-nef-audit.v1",
    "status": status,
    "edge": "4A1/MW13 --q4 orbit114--> 3A1/MW14, audited in reverse",
    "raw_fibre_in_3A1_coordinates": values(raw_fibre),
    "raw_component_degrees": {
        name: int(pairing(raw_fibre, Q3, curve))
        for name, curve in zip(component_names, component_walls)
    },
    "physical_component_reduction": {
        "reflection_sequence": component_reflections,
        "action": rows(component_action),
        "determinant": int(component_action.det()),
        "reduced_fibre_in_3A1_coordinates": values(physical_fibre),
        "old_fibre_degree": int(physical_fibre[1]),
        "q": int(physical_fibre[0] * physical_fibre[1]),
        "isotropic": pairing(physical_fibre, Q3, physical_fibre) == 0,
        "component_degrees": component_degrees,
    },
    "full_marked_transport": {
        "physical_4A1_basis_in_3A1_coordinates": rows(physical_parent_basis_in_3a1),
        "3A1_to_physical_4A1_basis": rows(three_a1_to_physical_parent),
        "basis_determinant": int(physical_parent_basis_in_3a1.det()),
        "mutually_inverse": True,
        "physical_zero_in_3A1_coordinates": values(physical_zero),
        "physical_zero_is_third_old_I2_nonidentity_component": True,
        "physical_4A1_root_classes_in_3A1_coordinates": [
            values(root) for root in physical_root_classes
        ],
        "physical_4A1_positive_frame": rows(-physical_parent_gram[2:, 2:]),
    },
    "known_effective_curve_gate": {
        "enumerated_degree_one_section_count": len(sections),
        "enumerated_section_minimum": min(known_section_degrees),
        "enumerated_section_degree_histogram": {
            str(value): known_section_degrees.count(value)
            for value in sorted(set(known_section_degrees))
        },
        "extra_curve_degrees": extra_degrees,
        "raw_stored_root_degrees_before_reflection": raw_stored_root_degrees,
        "root_class_correction": "apply the recorded physical component action to every abstract q114 root class",
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
    "proof_boundary": (
        "Exact lattice reconstruction, physical affine-I2 component reduction, "
        "all-section closest-vector gate, and complete finite old-horizontal "
        "wall enumeration. Equation construction and the corrected horizontal "
        "section identity are separate."
    ),
    "inputs": {
        str(path.relative_to(ROOT)): sha256(path)
        for path in (MANIFEST, BRIDGE, HORIZONTAL)
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDREVERSE4A1NEF|status={}|reflections={}|q={}|section_min={}|walls={}|output={}".format(
        status,
        len(component_reflections),
        payload["physical_component_reduction"]["q"],
        section_minimum,
        len(horizontal_walls),
        OUTPUT,
    )
)
