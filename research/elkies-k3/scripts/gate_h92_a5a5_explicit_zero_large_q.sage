#!/usr/bin/env sage -python
"""Exhaustively gate larger-q exits from the explicit-zero orbit12 child.

This is deliberately a pre-certificate search.  It enumerates dominant Weyl
orbits exactly, transports every fibre to the physical equation-A11 marking,
and rejects any fibre meeting an already explicit (-2)-curve negatively.
Survivors still require child-root adaptation, the all-section CVP nef test,
and full marked-U certification before they can enter a route.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
U2 = matrix(ZZ, ((0, 1), (1, 0)))

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--q", type=int, action="append", required=True)
parser.add_argument(
    "--output",
    type=Path,
    default=GENERATED / "elkies-k3-h3-a5a5-explicit-zero-large-q-gate.json",
)
args = parser.parse_args()

FRAME = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frame.txt"
O12_ZERO = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
ZERO_MISMATCH = GENERATED / "elkies-k3-h3-a11-quintic-bridge-zero-mismatch.json"
D12_FRAME = LOCAL / "q24-downstream-lift/d12-c10a-zero-frame.txt"
Q6_SHELL = LOCAL / "q24-downstream-lift/d12-c10a-zero-q6-all.json"
IDENTITY = LOCAL / "q24-orbit42-identity-halving-audit.json"
MATCHING = LOCAL / "q24-orbit42-identity-halving-qq.json"
INPUTS = (FRAME, O12_ZERO, CROSSOVERS, ZERO_MISMATCH, D12_FRAME, Q6_SHELL, IDENTITY, MATCHING)
OUTPUT = args.output.resolve()


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


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
    weights = []

    def recurse(prefix, norm):
        index = len(prefix)
        if index == len(component):
            weights.append((tuple(prefix), norm))
            return
        value = 0
        while True:
            added = inverse[index, index] * value**2
            added += 2 * value * sum(
                inverse[index, previous] * prefix[previous]
                for previous in range(index)
            )
            if norm + added > bound:
                break
            recurse(prefix + [value], norm + added)
            value += 1

    recurse([], QQ(0))
    return tuple(weights)


def highest_roots(cartan):
    half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
    roots = tuple(half) + tuple(-item for item in half)
    result = []
    for component in components(cartan):
        candidates = [
            item for item in roots
            if all(value >= 0 for value in item)
            and all(index in component or item[index] == 0 for index in range(cartan.nrows()))
        ]
        result.append(max(candidates, key=lambda item: sum(item)))
    return tuple(result)


frame = load_matrix(FRAME)
root_rank = 10
cartan = frame[:root_rank, :root_rank]
coupling = frame[:root_rank, root_rank:]
tail = frame[root_rank:, root_rank:]
height = tail - coupling.transpose() * cartan.inverse() * coupling
g_parent = block_diagonal_matrix(U2, -frame)

o12_zero = json.loads(O12_ZERO.read_text())
a11_to_parent = matrix(ZZ, o12_zero["selected"]["equation_A11_to_explicit_zero_basis"])
g_a11 = a11_to_parent.inverse() * g_parent * a11_to_parent.inverse().transpose()
assert a11_to_parent * g_a11 * a11_to_parent.transpose() == g_parent

# Physical equation-A11 curves, including all 18 exact identity-shell sections.
zero_mismatch = json.loads(ZERO_MISMATCH.read_text())
old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
old_simple = [
    vector(ZZ, [0, 0] + [-ZZ(index == node) for index in range(17)])
    for node in range(11)
]
old_affine = old_fibre + vector(ZZ, [0, 0] + [1] * 11 + [0] * 6)
a0 = vector(ZZ, zero_mismatch["correct_selected_R3_transport"]["oldI9_A0"]["child_coordinates"])
p24 = vector(ZZ, zero_mismatch["correct_selected_R3_transport"]["close_P24"]["child_coordinates"])

d12 = load_matrix(D12_FRAME)
q6_shell = json.loads(Q6_SHELL.read_text())
identity = json.loads(IDENTITY.read_text())
matching = json.loads(MATCHING.read_text())
selected_q6 = next(item for item in q6_shell["neighbors"] if int(item["orbit_index"]) == 64)
d12_to_a11 = block_diagonal_matrix(
    identity_matrix(ZZ, 2), matrix(ZZ, selected_q6["child_root_adapted_basis"])
) * matrix(ZZ, selected_q6["neighbor_basis"])
d12_inverse = d12_to_a11.inverse().change_ring(ZZ)
d12_root = d12[:12, :12]
d12_coupling = d12[:12, 12:]
shell = []
for values in identity["exact_model_R3_zero"]["identity_vectors"]:
    z = vector(ZZ, values)
    root_coefficients = -(z * d12_coupling.transpose()) * d12_root.inverse()
    section = vector(ZZ, [1, 1] + list(map(ZZ, root_coefficients)) + list(z))
    shell.append(section * d12_inverse)
mapping = matching["matching"]["mappings_abstract_to_equation"][7]
reordered = [None] * 18
for abstract_index, equation_index in enumerate(mapping):
    reordered[equation_index] = shell[abstract_index]
shell = reordered

curves = shell + [old_zero] + old_simple + [old_affine, a0, p24]
names = (
    [f"shell_S{i}" for i in range(18)]
    + ["old_A11_zero"]
    + [f"old_A11_component_{i}" for i in range(11)]
    + ["old_A11_affine", "oldI9_A0", "close_P24"]
)
assert all(curve * g_a11 * curve == -2 for curve in curves)

crossovers = json.loads(CROSSOVERS.read_text())
targets = {
    item["target"]: vector(ZZ, item["target_fibre_in_state"])
    for item in crossovers["records"] if item["state"] == "equation_A11"
}

component_weight_lists = None
component_list = components(cartan)
cartan_inverse = cartan.inverse()
highest = highest_roots(cartan)
height_scale = lcm(entry.denominator() for entry in height.list())
scaled_height = (height_scale * height).change_ring(ZZ)
records = []
summaries = []

for q in sorted(set(args.q)):
    if q <= 0 or q % 2:
        raise ValueError("q must be a positive even integer for old-fibre degree 2")
    target_norm = ZZ(2 * q)
    mw_result = pari(scaled_height).qfminim(height_scale * target_norm)
    mw_map = {}
    for column in matrix(ZZ, mw_result[2]).columns():
        for sign in (1, -1):
            value = sign * vector(ZZ, column)
            if value == 0 or value * height * value > target_norm:
                continue
            canonical = min(tuple(value), tuple(-value))
            mw_map[canonical] = vector(ZZ, canonical)
    mw_vectors = tuple(sorted(mw_map.values(), key=lambda value: (value * height * value, tuple(value))))

    component_weight_lists = tuple(
        dominant_weights(cartan, component, QQ(target_norm))
        for component in component_list
    )
    combined = {}

    def combine(index, choices, norm):
        if index == len(component_weight_lists):
            combined.setdefault(norm, []).append(tuple(choices))
            return
        for values, weight_norm in component_weight_lists[index]:
            if norm + weight_norm <= target_norm:
                combine(index + 1, choices + [(values, weight_norm)], norm + weight_norm)

    combine(0, [], QQ(0))
    seen = set()
    orbit_count = 0
    survivor_count = 0
    for mw in mw_vectors:
        mw_norm = mw * height * mw
        for choices in combined.get(target_norm - mw_norm, ()):
            labels = vector(ZZ, [0] * root_rank)
            for component, (values, _) in zip(component_list, choices):
                for index, value in zip(component, values):
                    labels[index] = value
            root_coordinates = cartan_inverse * (labels - coupling * mw)
            if not all(value in ZZ for value in root_coordinates):
                continue
            witness = vector(ZZ, list(root_coordinates) + list(mw))
            key = tuple(witness)
            if key in seen:
                continue
            seen.add(key)
            orbit_count += 1
            affine_pairings = [int(2 - top * labels) for top in highest]
            if min(affine_pairings) < 0:
                continue
            fibre_parent = vector(ZZ, [q // 2, 2] + list(witness))
            assert fibre_parent * g_parent * fibre_parent == 0
            fibre_a11 = fibre_parent * a11_to_parent
            degrees = [int(curve * g_a11 * fibre_a11) for curve in curves]
            if min(degrees) < 0:
                continue
            survivor_count += 1
            marked = {
                target: int(fibre_a11 * g_a11 * target_fibre)
                for target, target_fibre in targets.items()
            }
            record = {
                "candidate_id": {"q": int(q), "old_fibre_degree": 2, "orbit_index": orbit_count},
                "fibre_in_parent": list(map(int, fibre_parent)),
                "fibre_in_equation_A11": list(map(int, fibre_a11)),
                "mw_projection": list(map(int, mw)),
                "dominant_labels": list(map(int, labels)),
                "parent_affine_component_pairings": affine_pairings,
                "explicit_curve_degrees": dict(zip(names, degrees)),
                "explicit_degree_zero_curves": [names[i] for i, value in enumerate(degrees) if value == 0],
                "explicit_degree_one_curves": [names[i] for i, value in enumerate(degrees) if value == 1],
                "marked_target_degrees": marked,
                "coordinate_growth_max": int(max(abs(value) for value in fibre_parent)),
            }
            records.append(record)
    summaries.append({
        "q": int(q),
        "mw_projection_representatives": len(mw_vectors),
        "pari_vector_count": int(mw_result[0]),
        "dominant_orbits": orbit_count,
        "explicit_curve_nef_survivors": survivor_count,
    })
    print(
        f"A5A5LARGEGATE|q={q}|mw={len(mw_vectors)}|orbits={orbit_count}|survivors={survivor_count}|status=PASS",
        flush=True,
    )

records.sort(key=lambda item: (
    item["marked_target_degrees"]["pinned_R17"],
    -len(item["explicit_degree_one_curves"]),
    item["coordinate_growth_max"],
))
payload = {
    "schema": "elkies-k3.h3-a5a5-explicit-zero-large-q-gate.v1",
    "status": "PASS_EXACT_A5A5_EXPLICIT_ZERO_LARGE_Q_GATE",
    "q_values": sorted(set(args.q)),
    "summaries": summaries,
    "survivors": records,
    "proof_boundary": (
        "Exhaustive dominant-Weyl enumeration and exact nonnegative intersections with the 18 identity-shell "
        "sections, physical A11 zero/components/affine component, A0, and P24. Child root/component data, "
        "all-section nefness, marked U, and full transports are not certified here."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in INPUTS},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"A5A5LARGEGATE|total_survivors={len(records)}|status={payload['status']}|output={OUTPUT}", flush=True)
