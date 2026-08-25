#!/usr/bin/env sage -python
"""Certify the semistable MW2 -> pinned R17 reverse suffix chamber-wise."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "elkies-k3/data/fibrations"
RANK17 = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
TRANSPORT = ROOT / "artifacts/generated-results/elkies-k3-h3-semistable-mw2-pinned-transport.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h3-semistable-mw2-reverse-suffix-nef.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


def components(cartan):
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


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    half = [vector(ZZ, column) for column in matrix(ZZ, result[2]).columns()]
    roots = tuple(half + [-root for root in half])
    root_basis = matrix(ZZ, [list(root) for root in roots]).row_module().basis_matrix()
    root_gram = root_basis * gram * root_basis.transpose()
    return roots, root_basis, (root_basis.rank(), count, abs(ZZ(root_gram.det())))


def deterministic_simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    positive = [root for root in roots if next(value for value in root if value != 0) > 0]
    positive_set = {tuple(root) for root in positive}
    simple = [root for root in positive if not any(tuple(root - left) in positive_set for left in positive)]
    simple = matrix(ZZ, [list(root) for root in simple])
    assert simple.nrows() == simple.rank() == data[0]
    return simple, simple * gram * simple.transpose()


def root_adaptation(gram):
    unused, root_basis, invariants = roots_and_data(gram)
    root_rank = invariants[0]
    smith, smith_left, smith_right = root_basis.smith_form()
    assert smith == smith_left * root_basis * smith_right
    assert tuple(abs(smith[index, index]) for index in range(root_rank)) == (1,) * root_rank
    simple, cartan = deterministic_simple_roots(gram)
    completion = smith_right.inverse()
    adapted_basis = simple.stack(completion[root_rank:])
    adapted = adapted_basis * gram * adapted_basis.transpose()
    coupling = adapted[:root_rank, root_rank:]
    height = adapted[root_rank:, root_rank:] - coupling.transpose() * cartan.inverse() * coupling
    scale = lcm(entry.denominator() for entry in height.list())
    quotient_change = block_diagonal_matrix(
        identity_matrix(ZZ, root_rank), matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram()).transpose()
    )
    adapted_basis = quotient_change * adapted_basis
    adapted = adapted_basis * gram * adapted_basis.transpose()
    return adapted, adapted_basis, invariants


def highest_roots(cartan):
    half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
    roots = tuple(half) + tuple(-item for item in half)
    result = []
    for component in components(cartan):
        candidates = [item for item in roots if all(value >= 0 for value in item)
                      and all(index in component or item[index] == 0 for index in range(cartan.nrows()))]
        result.append(max(candidates, key=lambda item: sum(item)))
    return tuple(result)


def nef_audit(source, fibre):
    g_source = block_diagonal_matrix(U2, -source)
    assert fibre * g_source * fibre == 0
    assert gcd(tuple(g_source * fibre)) == 1
    adapted, adapted_basis, root_data = root_adaptation(source)
    change = block_diagonal_matrix(identity_matrix(ZZ, 2), adapted_basis)
    fibre_adapted = fibre * change.inverse().change_ring(ZZ)
    assert fibre_adapted * block_diagonal_matrix(U2, -adapted) * fibre_adapted == 0
    degree = ZZ(fibre_adapted[1])
    assert degree > 0
    root_rank = root_data[0]
    labels = vector(ZZ, fibre_adapted[2:]) * adapted[:, :root_rank]
    affine = [degree - top * labels for top in highest_roots(adapted[:root_rank, :root_rank])]
    center = vector(QQ, fibre_adapted[2:]) / degree
    closest = vector(ZZ, next(IntegralLattice(adapted).enumerate_close_vectors(center)))
    distance = (closest - center) * adapted * (closest - center)
    minimum_section = (QQ(degree) / 2) * distance - degree
    return {
        "primitive_nef_isotropic_fibre": entries(fibre),
        "old_fibre_degree": int(degree),
        "component_pairings": entries(labels),
        "affine_pairings": entries(affine),
        "closest_section_vector": entries(closest),
        "closest_section_distance": str(distance),
        "minimum_section_intersection": str(minimum_section),
        "root_data": list(map(int, root_data)),
        "nef": bool(min(labels) >= 0 and min(affine) >= 0 and minimum_section >= 0),
    }


def reduce_component_chamber(source, fibre):
    adapted, adapted_basis, root_data = root_adaptation(source)
    root_rank = root_data[0]
    change = block_diagonal_matrix(identity_matrix(ZZ, 2), adapted_basis)
    value = vector(ZZ, fibre * change.inverse().change_ring(ZZ))
    g_adapted = block_diagonal_matrix(U2, -adapted)
    cartan = adapted[:root_rank, :root_rank]
    wall_roots = []
    wall_names = []
    for index in range(root_rank):
        wall_roots.append(vector(ZZ, [0, 0] + [-ZZ(index == other) for other in range(17)]))
        wall_names.append(f"simple_{index}")
    for component_index, top in enumerate(highest_roots(cartan)):
        wall_roots.append(vector(ZZ, [1, 0] + list(top) + [0] * (17 - root_rank)))
        wall_names.append(f"affine_{component_index}")
    assert all(root * g_adapted * root == -2 for root in wall_roots)
    reflections = []
    for unused in range(100000):
        pairings = [ZZ(value * g_adapted * root) for root in wall_roots]
        negative = next((index for index, pairing in enumerate(pairings) if pairing < 0), None)
        if negative is None:
            break
        pairing = pairings[negative]
        value += pairing * wall_roots[negative]
        reflections.append({"wall": wall_names[negative], "pairing": int(pairing)})
    else:
        raise RuntimeError("component Weyl reduction did not terminate")
    reduced = vector(ZZ, value * change)
    return reduced, reflections


def construct_child(source, fibre):
    ns = block_diagonal_matrix(U2, -source)
    pairings = list(ns * fibre)
    current = ZZ(0)
    mate = [ZZ(0)] * 19
    for index, pairing in enumerate(pairings):
        if pairing == 0:
            continue
        divisor, left, right = xgcd(current, ZZ(pairing))
        mate = [left * entry for entry in mate]
        mate[index] += right
        current = divisor
    assert abs(current) == 1
    if current == -1:
        mate = [-entry for entry in mate]
    mate = vector(ZZ, mate)
    mate -= ZZ(mate * ns * mate // 2) * fibre
    kernel = matrix(ZZ, [list(fibre * ns), list(mate * ns)]).right_kernel_matrix()
    child = -(kernel * ns * kernel.transpose())
    basis = matrix(ZZ, [list(fibre), list(mate)] + [list(row) for row in kernel.rows()])
    assert abs(basis.det()) == 1
    assert basis * ns * basis.transpose() == block_diagonal_matrix(U2, -child)
    return child, basis


def qfisom_row(source, target):
    raw = pari(source).qfisom(pari(target))
    assert raw != 0
    value = matrix(ZZ, raw)
    for candidate in (value, value.transpose(), value.inverse(), value.inverse().transpose()):
        candidate = matrix(ZZ, candidate)
        if abs(candidate.det()) == 1 and candidate * source * candidate.transpose() == target:
            return candidate
    raise ArithmeticError("qfisom result had no verified row convention")


transport = json.loads(TRANSPORT.read_text())
frame_paths = (
    DATA / "mw2_a5_a4_a3a3_frame.txt",
    DATA / "mw3_d6_a5_a3_frame.txt",
    DATA / "q25_mw4_frame.txt",
    DATA / "q25_mw7_frame.txt",
    RANK17,
)
frames = list(map(load_matrix, frame_paths))
forward = [
    matrix(ZZ, transport["steps"][3]["semistable_to_D6_transport"]),
    matrix(ZZ, transport["steps"][2]["transport"]).inverse().change_ring(ZZ),
    matrix(ZZ, transport["steps"][1]["transport"]).inverse().change_ring(ZZ),
    matrix(ZZ, transport["steps"][0]["transport"]).inverse().change_ring(ZZ),
]
labels = ("q14_to_D6", "inverse_q4_to_MW4", "inverse_q4_to_MW7", "inverse_q25_to_pinned_R17")
records = []
cumulative = identity_matrix(ZZ, 19)
for index, (label, transition) in enumerate(zip(labels, forward)):
    g_source = block_diagonal_matrix(U2, -frames[index])
    g_target = block_diagonal_matrix(U2, -frames[index + 1])
    assert abs(transition.det()) == 1
    assert transition * g_source * transition.transpose() == g_target
    original_fibre = vector(ZZ, transition.row(0))
    fibre, reflections = reduce_component_chamber(frames[index], original_fibre)
    raw_child, raw_transition = construct_child(frames[index], fibre)
    endpoint_isometry = qfisom_row(raw_child, frames[index + 1])
    transition = block_diagonal_matrix(identity_matrix(ZZ, 2), endpoint_isometry) * raw_transition
    assert transition * g_source * transition.transpose() == g_target
    audit = nef_audit(frames[index], fibre)
    records.append({
        "label": label,
        "source_frame": str(frame_paths[index].relative_to(ROOT)),
        "target_frame": str(frame_paths[index + 1].relative_to(ROOT)),
        "transport": rows(transition),
        "inverse_transport": rows(transition.inverse().change_ring(ZZ)),
        "forward_determinant": int(transition.det()),
        "inverse_determinant": int(transition.inverse().det()),
        "original_inverse_transport_fibre": entries(original_fibre),
        "component_weyl_reflections": reflections,
        "raw_child_to_pinned_target_isometry": rows(endpoint_isometry),
        "nef_audit": audit,
    })
    cumulative = transition * cumulative

all_nef = all(record["nef_audit"]["nef"] for record in records)
payload = {
    "schema": "elkies-k3.h3-semistable-mw2-reverse-suffix-nef.v1",
    "status": "PASS_EXACT_SEMISTABLE_MW2_TO_PINNED_R17_NEF_SUFFIX" if all_nef else "REJECT_NONNEF_REVERSE_SUFFIX",
    "route": ["A5+A4+2A3/MW2", "D6+A5+A3/MW3", "D4+A3+2A2+2A1/MW4", "A3+7A1/MW7", "pinned rootless/MW17"],
    "steps": records,
    "composite_semistable_to_pinned_R17": rows(cumulative),
    "composite_pinned_R17_to_semistable": rows(cumulative.inverse().change_ring(ZZ)),
    "composite_forward_determinant": int(cumulative.det()),
    "composite_inverse_determinant": int(cumulative.inverse().det()),
    "proof_boundary": "Exact marked U at every edge, complete root-component and all-section CVP nef checks, full integral transports both directions, and exact pinned rank17 endpoint. No equation model is asserted.",
    "inputs": {
        "paths": [str(TRANSPORT.relative_to(ROOT))] + [str(path.relative_to(ROOT)) for path in frame_paths],
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in (TRANSPORT,) + frame_paths},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
for record in records:
    audit = record["nef_audit"]
    print(f"SEMREVNEF|step={record['label']}|degree={audit['old_fibre_degree']}|min_section={audit['minimum_section_intersection']}|nef={int(audit['nef'])}", flush=True)
print(f"SEMREVNEF|status={payload['status']}|output={OUTPUT}", flush=True)
