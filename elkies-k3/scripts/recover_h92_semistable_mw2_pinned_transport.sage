#!/usr/bin/env sage -python
"""Recover and certify the pinned q25,4,4,14 path to A5+A4+2A3/MW2."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "elkies-k3/data/fibrations"
RANK17 = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h3-semistable-mw2-pinned-transport.json"
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


def transition(parent, a, b, witness):
    ns = block_diagonal_matrix(U2, -parent)
    fibre = vector(ZZ, [a, b] + list(witness))
    assert fibre * ns * fibre == 0
    pairings = list(ns * fibre)
    current = ZZ(0)
    mate = [ZZ(0)] * 19
    for index, value in enumerate(pairings):
        if value == 0:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        mate = [left * entry for entry in mate]
        mate[index] += right
        current = divisor
    if abs(current) != 1:
        return None
    if current == -1:
        mate = [-entry for entry in mate]
    mate = vector(ZZ, mate)
    square = ZZ(mate * ns * mate)
    assert square % 2 == 0
    mate -= (square // 2) * fibre
    kernel = matrix(ZZ, [list(fibre * ns), list(mate * ns)]).right_kernel_matrix()
    child = -(kernel * ns * kernel.transpose())
    basis = matrix(ZZ, [list(fibre), list(mate)] + [list(row) for row in kernel.rows()])
    assert abs(basis.det()) == 1
    assert basis * ns * basis.transpose() == block_diagonal_matrix(U2, -child)
    return child, basis


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


def dominant_weights(cartan, component, bound):
    inverse = cartan.matrix_from_rows_and_columns(component, component).inverse()
    result = []

    def recurse(prefix, norm):
        index = len(prefix)
        if index == len(component):
            result.append((tuple(prefix), norm))
            return
        value = 0
        while True:
            added = inverse[index, index] * value**2
            added += 2 * value * sum(inverse[index, previous] * prefix[previous] for previous in range(index))
            if norm + added > bound:
                break
            recurse(prefix + [value], norm + added)
            value += 1

    recurse([], QQ(0))
    return tuple(result)


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    if count == 0:
        return (), matrix(ZZ, 0, gram.nrows()), (0, 0, 1)
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
    assert abs(adapted_basis.det()) == 1
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


def qfisom_row(source, target):
    """Return Q with Q*source*Q^t=target, or None."""
    raw = pari(source).qfisom(pari(target))
    if raw == 0:
        return None
    value = matrix(ZZ, raw)
    candidates = (value, value.transpose(), value.inverse(), value.inverse().transpose())
    for candidate in candidates:
        candidate = matrix(ZZ, candidate)
        if abs(candidate.det()) == 1 and candidate * source * candidate.transpose() == target:
            return candidate
    return None


paths = (
    RANK17,
    DATA / "q25_mw7_frame.txt",
    DATA / "q25_mw4_frame.txt",
    DATA / "mw3_d6_a5_a3_frame.txt",
    DATA / "mw2_a5_a4_a3a3_frame.txt",
)
frames = list(map(load_matrix, paths))
known = (
    (5, 5, (-1, 0, -4, 3, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, -3, 0, 0)),
    (2, 2, (-1, -2, 1, 0, 1, 1, 2, -3, 0, -2, 0, 1, 0, 0, -1, 0, 0)),
    (2, 2, (-1, -1, 2, 0, -2, 0, -1, 0, -1, 0, 0, -1, -1, 1, 1, 0, 0)),
)
transitions = []
step_records = []
for index, (a, b, values) in enumerate(known):
    witness = vector(ZZ, values)
    result = transition(frames[index], a, b, witness)
    assert result is not None
    child, basis = result
    assert child == frames[index + 1]
    transitions.append(basis)
    step_records.append({"q": a * b, "factorization": [a, b], "witness": entries(witness), "transport": rows(basis)})

# The discovery beam was not retained.  Recover the useful reverse q14 orbit
# from the stored semistable endpoint to the pinned D6 reverse hub.  This is
# the direction needed by a future A11 -> semistable -> pinned route.
parent = frames[4]
target = frames[3]
adapted_parent, adapted_basis, parent_root_data = root_adaptation(parent)
assert parent_root_data == (15, 74, 480)
parent_change = block_diagonal_matrix(identity_matrix(ZZ, 2), adapted_basis)
root_rank = 15
cartan = adapted_parent[:root_rank, :root_rank]
coupling = adapted_parent[:root_rank, root_rank:]
height = adapted_parent[root_rank:, root_rank:] - coupling.transpose() * cartan.inverse() * coupling
target_norm = ZZ(28)
height_scale = lcm(entry.denominator() for entry in height.list())
mw_result = pari((height_scale * height).change_ring(ZZ)).qfminim(height_scale * target_norm)
mw_map = {tuple([0] * (17 - root_rank)): vector(ZZ, [0] * (17 - root_rank))}
for column in matrix(ZZ, mw_result[2]).columns():
    for sign in (1, -1):
        value = sign * vector(ZZ, column)
        if value * height * value <= target_norm:
            canonical = min(tuple(value), tuple(-value))
            mw_map[canonical] = vector(ZZ, canonical)
mw_vectors = tuple(sorted(mw_map.values(), key=lambda value: (value * height * value, tuple(value))))
component_list = components(cartan)
weight_lists = tuple(dominant_weights(cartan, component, QQ(target_norm)) for component in component_list)
combined = {}


def combine(index, choices, norm):
    if index == len(weight_lists):
        combined.setdefault(norm, []).append(tuple(choices))
        return
    for values, weight_norm in weight_lists[index]:
        if norm + weight_norm <= target_norm:
            combine(index + 1, choices + [(values, weight_norm)], norm + weight_norm)


combine(0, [], QQ(0))
cartan_inverse = cartan.inverse()
highest = highest_roots(cartan)
matches = []
orbit_count = component_nef_count = section_nef_count = root_data_hit_count = 0
seen = set()
for mw in mw_vectors:
    for choices in combined.get(target_norm - mw * height * mw, ()):
        labels = vector(ZZ, [0] * root_rank)
        for component, (values, unused_norm) in zip(component_list, choices):
            for index, value in zip(component, values):
                labels[index] = value
        root_coordinates = cartan_inverse * (labels - coupling * mw)
        if not all(value in ZZ for value in root_coordinates):
            continue
        witness_adapted = vector(ZZ, list(root_coordinates) + list(mw))
        if tuple(witness_adapted) in seen:
            continue
        seen.add(tuple(witness_adapted))
        orbit_count += 1
        for a, b in ((2, 7), (7, 2)):
            affine = [ZZ(b - top * labels) for top in highest]
            if min(affine) < 0:
                continue
            component_nef_count += 1
            center = vector(QQ, witness_adapted) / b
            closest = vector(ZZ, next(IntegralLattice(adapted_parent).enumerate_close_vectors(center)))
            distance = (closest - center) * adapted_parent * (closest - center)
            minimum_section_intersection = (QQ(b) / 2) * distance - b
            if minimum_section_intersection < 0:
                continue
            section_nef_count += 1
            fibre_parent = vector(ZZ, [a, b] + list(witness_adapted)) * parent_change
            assert tuple(fibre_parent[:2]) == (a, b)
            witness = vector(ZZ, fibre_parent[2:])
            result = transition(parent, a, b, witness)
            if result is None:
                continue
            child, raw_basis = result
            unused_roots, unused_basis, child_root_data = roots_and_data(child)
            if child_root_data != (14, 102, 96):
                continue
            root_data_hit_count += 1
            isometry = qfisom_row(child, target)
            if isometry is not None:
                target_basis = block_diagonal_matrix(identity_matrix(ZZ, 2), isometry) * raw_basis
                matches.append((tuple(witness_adapted), a, b, witness, target_basis, entries(labels), list(map(int, affine)), str(distance), str(minimum_section_intersection), rows(isometry)))
assert matches
matches.sort(key=lambda item: item[0])
unused_adapted_witness, reverse_a, reverse_b, witness, semistable_to_d6, labels, affine, closest_distance, minimum_section_intersection, endpoint_isometry = matches[0]
d6_to_semistable = semistable_to_d6.inverse().change_ring(ZZ)
transitions.append(d6_to_semistable)
step_records.append({
    "q": 14,
    "certified_nef_direction": "semistable_MW2_to_D6_reverse_hub",
    "reverse_factorization": [reverse_a, reverse_b],
    "reverse_witness_in_semistable_frame": entries(witness),
    "reverse_component_pairings": labels,
    "reverse_affine_pairings": affine,
    "reverse_closest_section_distance": closest_distance,
    "reverse_minimum_section_intersection": minimum_section_intersection,
    "raw_reverse_child_to_pinned_D6_frame": endpoint_isometry,
    "semistable_to_D6_transport": rows(semistable_to_d6),
    "D6_to_semistable_transport": rows(d6_to_semistable),
})

cumulative = identity_matrix(ZZ, 19)
for index, basis in enumerate(transitions):
    cumulative = basis * cumulative
    assert cumulative * block_diagonal_matrix(U2, -frames[0]) * cumulative.transpose() == block_diagonal_matrix(U2, -frames[index + 1])
inverse = cumulative.inverse().change_ring(ZZ)
assert abs(cumulative.det()) == abs(inverse.det()) == 1
assert inverse * block_diagonal_matrix(U2, -frames[-1]) * inverse.transpose() == block_diagonal_matrix(U2, -frames[0])

# Exact component/root data for endpoint identification.
endpoint_frame = frames[-1]
root_result = pari(endpoint_frame).qfminim(2)
root_half = matrix(ZZ, root_result[2]).columns()
roots = tuple(vector(ZZ, root) for root in root_half) + tuple(-vector(ZZ, root) for root in root_half)
root_basis = matrix(ZZ, [list(root) for root in roots]).row_module().basis_matrix()
root_gram = root_basis * endpoint_frame * root_basis.transpose()
assert (root_basis.rank(), ZZ(root_result[0]), abs(ZZ(root_gram.det()))) == (15, 74, 480)

payload = {
    "schema": "elkies-k3.h3-semistable-mw2-pinned-transport.v1",
    "status": "PASS_EXACT_PINNED_R17_SEMISTABLE_MW2_TRANSPORT_WITH_NEF_REVERSE_Q14",
    "route": [
        "pinned rootless/MW17",
        "q25 A3+7A1/MW7",
        "q4 D4+A3+2A2+2A1/MW4",
        "q4 D6+A5+A3/MW3",
        "q14 A5+A4+2A3/MW2",
    ],
    "q14_recovery": {
        "mw_projection_representatives": len(mw_vectors),
        "dominant_weyl_orbits": orbit_count,
        "component_nef_orbits": component_nef_count,
        "all_section_nef_orbits": section_nef_count,
        "root_data_hits": root_data_hit_count,
        "exact_endpoint_isometry_matches": len(matches),
        "selected_reverse_factorization": [reverse_a, reverse_b],
        "selected_reverse_witness": entries(witness),
    },
    "steps": step_records,
    "endpoint": {
        "ade": "A5+A4+2A3",
        "mw_rank": 2,
        "root_data": [15, 74, 480],
        "root_basis_in_endpoint_frame": rows(root_basis),
        "root_gram": rows(root_gram),
    },
    "transport": {
        "pinned_R17_to_semistable_basis": rows(cumulative),
        "semistable_to_pinned_R17_basis": rows(inverse),
        "forward_determinant": int(cumulative.det()),
        "inverse_determinant": int(inverse.det()),
    },
    "pinned_R17_marked_U_in_semistable": {
        "fibre": entries(inverse.row(0)),
        "isotropic_mate": entries(inverse.row(1)),
        "zero": entries(inverse.row(1) - inverse.row(0)),
        "gram": [[0, 1], [1, 0]],
    },
    "proof_boundary": "Exact primitive isotropic classes, raw child equality/isometry at every edge, endpoint roots, and integral transports both ways. The semistable-to-D6 q14 direction has a complete component/all-section nef certificate. Nefness of the remaining inverse q4,q4,q25 suffix is not re-certified here, and no equation model is asserted for the semistable endpoint.",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in paths],
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "SEMISTABLEPIN|q14_witness={}|shell={}|matches={}|det_fwd={}|det_inv={}|root=15,74,480|status={}".format(
        ",".join(map(str, witness)), orbit_count, len(matches), cumulative.det(), inverse.det(), payload["status"]
    ),
    flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)
