#!/usr/bin/env sage -python
"""Certify the q6/orbit84 equation-cost lateral move from orbit1991 MW4."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
U2 = matrix(ZZ, ((0, 1), (1, 0)))
SCORES = GENERATED / "elkies-k3-h3-o1991-explicit-zero-equation-cost-neighbors.json"
NEIGHBORS = GENERATED / "elkies-k3-h3-a11-o1991-explicit-zero-degree2-neighbors.json"
ZERO_FRAMES = GENERATED / "elkies-k3-h3-a11-q8-orbit1991-explicit-zero-frames.json"
OUTPUT = GENERATED / "elkies-k3-h3-o1991-q6-orbit84-lattice-certificate.json"
FRAME_OUTPUT = GENERATED / "elkies-k3-h3-o1991-q6-orbit84-frame.txt"
INPUTS = (SCORES, NEIGHBORS, ZERO_FRAMES)


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


def connected_components(cartan):
    unseen = set(range(cartan.nrows()))
    result = []
    while unseen:
        todo = [unseen.pop()]
        component = []
        while todo:
            node = todo.pop()
            component.append(node)
            for other in tuple(unseen):
                if cartan[node, other]:
                    unseen.remove(other)
                    todo.append(other)
        result.append(tuple(sorted(component)))
    return tuple(result)


def highest_roots(cartan):
    half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
    roots = tuple(half) + tuple(-item for item in half)
    result = []
    for component in connected_components(cartan):
        candidates = [
            item
            for item in roots
            if all(value >= 0 for value in item)
            and all(index in component or item[index] == 0 for index in range(cartan.nrows()))
        ]
        result.append(max(candidates, key=lambda item: sum(item)))
    return tuple(result)


scores = json.loads(SCORES.read_text())
neighbors = json.loads(NEIGHBORS.read_text())
zero_frames = json.loads(ZERO_FRAMES.read_text())
assert scores["status"] == "PASS_EXACT_O1991_EXPLICIT_ZERO_EQUATION_COST_SCORING"
best = scores["best_candidate"]
assert best["candidate_id"] == {"q": 6, "old_fibre_degree": 2, "orbit_index": 84}
raw = next(
    item for item in neighbors["neighbors"]
    if int(item["q"]) == 6 and int(item["orbit_index"]) == 84
)
parent_path = ROOT / neighbors["frame"]
parent = load_matrix(parent_path)
g_parent = block_diagonal_matrix(U2, -parent)
root_rank = 13
root = parent[:root_rank, :root_rank]
witness = vector(ZZ, raw["witness"])
fibre = vector(ZZ, raw["fiber"])
old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
assert fibre == vector(ZZ, [3, 2] + list(witness))
assert fibre * g_parent * fibre == 0 and gcd(tuple(g_parent * fibre)) == 1
assert fibre * g_parent * old_fibre == 2 and fibre * g_parent * old_zero == 1

component_pairings = tuple(map(ZZ, raw["dominant_labels"]))
affine_pairings = tuple(
    ZZ(2 - highest * vector(ZZ, component_pairings))
    for highest in highest_roots(root)
)
assert all(value >= 0 for value in component_pairings + affine_pairings)

center = vector(QQ, witness) / 2
closest = vector(ZZ, next(IntegralLattice(parent).enumerate_close_vectors(center)))
closest_distance = (closest - center) * parent * (closest - center)
assert closest_distance >= 2
assert witness * parent * witness % 4 == 0

neighbor_basis = matrix(ZZ, raw["neighbor_basis"])
adapted_basis = matrix(ZZ, raw["child_root_adapted_basis"])
child = matrix(ZZ, raw["child_root_adapted_frame"])
transition = block_diagonal_matrix(identity_matrix(ZZ, 2), adapted_basis) * neighbor_basis
inverse = transition.inverse().change_ring(ZZ)
g_child = block_diagonal_matrix(U2, -child)
assert abs(transition.det()) == 1
assert transition * g_parent * transition.transpose() == g_child
assert inverse * g_child * inverse.transpose() == g_parent
assert tuple(map(int, raw["child_root_data"])) == (13, 188, 24)
assert raw["child_ade"] == "A1+A2+D10" and int(raw["child_mw_rank"]) == 4

equation_a11_to_parent = matrix(
    ZZ, zero_frames["selected"]["equation_A11_to_explicit_zero_basis"]
)
equation_a11_to_child = transition * equation_a11_to_parent
child_to_equation_a11 = equation_a11_to_child.inverse().change_ring(ZZ)
assert abs(equation_a11_to_child.det()) == 1

FRAME_OUTPUT.write_text(
    "# q6 orbit84 D10+A1+A2/MW4 equation-cost lateral child\n"
    + "\n".join(" ".join(map(str, item)) for item in child.rows())
    + "\n"
)
payload = {
    "schema": "elkies-k3.h3-o1991-q6-orbit84-lattice-certificate.v1",
    "status": "PASS_EXACT_O1991_Q6_ORBIT84_LATTICE_CERTIFICATE",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS] + [str(parent_path.relative_to(ROOT))],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS + (parent_path,)
        },
    },
    "selection": best,
    "edge": {
        "q": 6,
        "factorization": [3, 2],
        "orbit_index": 84,
        "primitive_nef_isotropic_fibre": entries(fibre),
        "component_pairings": entries(vector(ZZ, component_pairings)),
        "affine_pairings": entries(vector(ZZ, affine_pairings)),
        "closest_section_vector": entries(closest),
        "closest_section_distance": str(closest_distance),
        "minimum_section_intersection": str(closest_distance - 2),
        "bisection_parity_exclusion": True,
        "nef": True,
    },
    "marked_U": {
        "fibre_in_parent": entries(transition.row(0)),
        "isotropic_mate_in_parent": entries(transition.row(1)),
        "zero_in_parent": entries(transition.row(1) - transition.row(0)),
        "gram": [[0, 1], [1, 0]],
    },
    "child": {
        "ade": raw["child_ade"],
        "mw_rank": int(raw["child_mw_rank"]),
        "root_data": raw["child_root_data"],
        "frame": rows(child),
        "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
        "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    },
    "transport": {
        "parent_to_child_basis": rows(transition),
        "child_to_parent_basis": rows(inverse),
        "equation_A11_to_child_basis": rows(equation_a11_to_child),
        "child_to_equation_A11_basis": rows(child_to_equation_a11),
        "forward_determinant": int(transition.det()),
        "inverse_determinant": int(inverse.det()),
    },
    "route_status": (
        "Certified cheap lateral edge only. No child-to-pinned-R17 continuation "
        "has yet been certified, so the active lifting target is unchanged."
    ),
    "proof_boundary": (
        "Exact full-nef lattice/marked-U/root/transport certificate. The displayed "
        "RR ambient remains an equation-cost estimate; no equation lift is claimed."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "O1991Q6O84|PO={}|RR={}|closest={}|child=A1+A2+D10/MW4|"
    "root=13,188,24|det_fwd={}|det_inv={}|nef=1|status={}".format(
        best["horizontal"]["P_dot_O"], best["expected_RR_ambient"],
        closest_distance, transition.det(), inverse.det(), payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{OUTPUT.resolve()}", flush=True)
