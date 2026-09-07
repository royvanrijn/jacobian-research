#!/usr/bin/env sage -python
"""Score the complete equation-marked A5+A5 q4 shell for compiler cost."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
NEIGHBORS = GENERATED / "elkies-k3-h3-a11-o12-q4-degree2-all.json"
O12_CERT = GENERATED / "elkies-k3-h3-a11-q8-orbit12-lattice-certificate.json"
TARGET_RANKING = GENERATED / "elkies-k3-h3-a5a5-marked-target-neighbor-ranking.json"
PARENT_D12 = LOCAL / "q24-downstream-lift/d12-c10a-zero-frame.txt"
Q6 = LOCAL / "q24-downstream-lift/d12-c10a-zero-q6-all.json"
IDENTITY = LOCAL / "q24-orbit42-identity-halving-audit.json"
MATCHING = LOCAL / "q24-orbit42-identity-halving-qq.json"
ZERO_MISMATCH = GENERATED / "elkies-k3-h3-a11-quintic-bridge-zero-mismatch.json"
OUTPUT = GENERATED / "elkies-k3-h3-a5a5-equation-cost-neighbors.json"
INPUTS = (NEIGHBORS, O12_CERT, TARGET_RANKING, PARENT_D12, Q6, IDENTITY, MATCHING, ZERO_MISMATCH)
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


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


neighbors = json.loads(NEIGHBORS.read_text())
o12 = json.loads(O12_CERT.read_text())
target_ranking = json.loads(TARGET_RANKING.read_text())
q6 = json.loads(Q6.read_text())
identity = json.loads(IDENTITY.read_text())
matching = json.loads(MATCHING.read_text())
zero_mismatch = json.loads(ZERO_MISMATCH.read_text())
assert neighbors["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert target_ranking["status"] == "PASS_EXACT_A5A5_MARKED_TARGET_NEIGHBOR_RANKING"

parent = load_matrix(ROOT / neighbors["frame"])
g_parent = block_diagonal_matrix(U2, -parent)
root_rank = 10
root = parent[:root_rank, :root_rank]
coupling = parent[:root_rank, root_rank:]
height = parent[root_rank:, root_rank:] - coupling.transpose() * root.inverse() * coupling
root_lattice = IntegralLattice(root)
root_edges = [
    (left, right) for left in range(root_rank) for right in range(left + 1, root_rank)
    if root[left, right] == -1
]
equation_to_parent = matrix(ZZ, o12["transport"]["parent_to_child_basis"])
a11_shell = json.loads((LOCAL / "q24-a11-orbit64-q8-all.json").read_text())
g_a11 = block_diagonal_matrix(U2, -load_matrix(ROOT / a11_shell["frame"]))
assert equation_to_parent * g_a11 * equation_to_parent.transpose() == g_parent

# Reconstruct the 18 exact identity-shell sections in equation A11 coordinates.
d12 = load_matrix(PARENT_D12)
selected_q6 = next(item for item in q6["neighbors"] if int(item["orbit_index"]) == 64)
d12_to_a11 = block_diagonal_matrix(
    identity_matrix(ZZ, 2), matrix(ZZ, selected_q6["child_root_adapted_basis"])
) * matrix(ZZ, selected_q6["neighbor_basis"])
d12_root_rank = 12
d12_root = d12[:d12_root_rank, :d12_root_rank]
d12_coupling = d12[:d12_root_rank, d12_root_rank:]
shell = []
for values in identity["exact_model_R3_zero"]["identity_vectors"]:
    z = vector(ZZ, values)
    root_coefficients = -(z * d12_coupling.transpose()) * d12_root.inverse()
    section = vector(ZZ, [1, 1] + list(map(ZZ, root_coefficients)) + list(z))
    shell.append(section * d12_to_a11.inverse().change_ring(ZZ))
mapping = matching["matching"]["mappings_abstract_to_equation"][7]
reordered = [None] * 18
for abstract_index, equation_index in enumerate(mapping):
    reordered[equation_index] = shell[abstract_index]
shell = reordered
a0 = vector(ZZ, zero_mismatch["correct_selected_R3_transport"]["oldI9_A0"]["child_coordinates"])
p24 = vector(ZZ, zero_mismatch["correct_selected_R3_transport"]["close_P24"]["child_coordinates"])
old_a11_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_a11_zero = vector(ZZ, [-1, 1] + [0] * 17)
a11_simple = [
    vector(ZZ, [0, 0] + [-ZZ(index == node) for index in range(17)])
    for node in range(11)
]
a11_affine = old_a11_fibre + vector(
    ZZ, [0, 0] + list(highest_roots((-g_a11[2:, 2:])[:11, :11])[0]) + [0] * 6
)
explicit_curves = shell + [a0, p24, old_a11_zero] + a11_simple + [a11_affine]
explicit_names = (
    [f"shell_S{index}" for index in range(18)]
    + ["oldI9_A0", "close_P24", "old_A11_zero"]
    + [f"old_A11_component_{index}" for index in range(11)]
    + ["old_A11_affine"]
)
assert len(explicit_curves) == len(explicit_names) == 33
assert all(curve * g_a11 * curve == -2 for curve in explicit_curves)


def connected_layers(coefficients):
    magnitudes = [abs(ZZ(value)) for value in coefficients]
    total = 0
    previous = 0
    for level in sorted(set(value for value in magnitudes if value)):
        active = set(index for index, value in enumerate(magnitudes) if value >= level)
        count = 0
        while active:
            count += 1
            todo = [active.pop()]
            while todo:
                node = todo.pop()
                for left, right in root_edges:
                    other = right if left == node else left if right == node else None
                    if other in active:
                        active.remove(other)
                        todo.append(other)
        total += (level - previous) * count
        previous = level
    return int(total)


profile_cache = {}


def best_horizontal(z):
    key = tuple(z)
    if key in profile_cache:
        return profile_cache[key]
    z = vector(ZZ, z)
    base = vector(ZZ, [0] * root_rank + list(z))
    dual = vector(QQ, base * parent[:, :root_rank]) * root.inverse()
    iterator = root_lattice.enumerate_close_vectors(-dual)
    minimum = None
    profiles = []
    for unused in range(100000):
        shift = vector(ZZ, next(iterator))
        lifted = base + vector(ZZ, list(shift) + [0] * 7)
        norm = QQ(lifted * parent * lifted)
        if minimum is None:
            minimum = norm
        elif norm > minimum:
            break
        pole = (norm - 4) / 2
        if pole in ZZ and pole >= 0:
            section = vector(ZZ, [ZZ(pole) + 1, 1] + list(lifted))
            profiles.append((section, ZZ(pole)))
    assert profiles
    profiles.sort(key=lambda item: (item[1], tuple(item[0])))
    profile_cache[key] = profiles[0]
    return profiles[0]


target_degree_by_orbit = {
    int(item["candidate_id"]["orbit_index"]): int(item["marked_target_degrees"]["pinned_R17"])
    for item in target_ranking["rankings_top_100"]["pinned_R17"]
}
parent_highest = highest_roots(root)
records = []
for raw in neighbors["neighbors"]:
    fibre = vector(ZZ, raw["fiber"])
    fibre_equation = fibre * equation_to_parent
    section, pole = best_horizontal(raw["mw_projection"])
    residual = fibre - vector(ZZ, [-1, 1] + [0] * 17) - section
    vertical = vector(ZZ, residual[2:2 + root_rank])
    layers = connected_layers(vertical)
    explicit_degrees = [int(curve * g_a11 * fibre_equation) for curve in explicit_curves]
    negative = [explicit_names[index] for index, value in enumerate(explicit_degrees) if value < 0]
    degree_zero = [explicit_names[index] for index, value in enumerate(explicit_degrees) if value == 0]
    degree_one = [explicit_names[index] for index, value in enumerate(explicit_degrees) if value == 1]
    labels = vector(ZZ, raw["dominant_labels"])
    affine = [int(2 - highest * labels) for highest in parent_highest]
    negative_affine = [index for index, value in enumerate(affine) if value < 0]
    rr = 2 + 2 * int(pole) + layers
    terms = {
        "explicit_non_nef_penalty": 1000000 if negative or negative_affine else 0,
        "P_dot_O": 900 * int(pole),
        "horizontal_degree": 500,
        "RR_ambient": 120 * rr,
        "vertical_layers": 60 * layers,
        "vertical_support": 25 * sum(value != 0 for value in vertical),
        "child_root_count": int(raw["child_root_data"][1]),
        "coordinate_growth": max(abs(int(value)) for value in fibre),
        "no_explicit_degree_one_curve": 4000 if not degree_one else 0,
        "explicit_degree_one_credit": -500 * min(len(degree_one), 6),
        "explicit_degree_zero_credit": -100 * min(len(degree_zero), 12),
    }
    record = {
        "candidate_id": {"q": 4, "old_fibre_degree": 2, "orbit_index": int(raw["orbit_index"])},
        "child": {"ade": raw["child_ade"], "mw_rank": int(raw["child_mw_rank"]), "root_data": raw["child_root_data"]},
        "declared_explicit_curve_nef_gate": "PASS" if not negative and not negative_affine else "REJECT",
        "horizontal": {"P_dot_O": int(pole), "section": entries(section), "vertical": entries(vertical), "vertical_layers": layers},
        "expected_RR_ambient": rr,
        "explicit_curve_degrees": dict(zip(explicit_names, explicit_degrees)),
        "explicit_degree_zero_curves": degree_zero,
        "explicit_degree_one_curves": degree_one,
        "negative_explicit_curves": negative,
        "parent_affine_component_pairings": affine,
        "negative_parent_affine_components": negative_affine,
        "pinned_R17_degree": target_degree_by_orbit.get(int(raw["orbit_index"])),
        "equation_cost_terms": terms,
        "equation_cost_score": int(sum(terms.values())),
    }
    records.append(record)

records.sort(key=lambda item: (item["equation_cost_score"], item["pinned_R17_degree"] or 10**100, item["candidate_id"]["orbit_index"]))
payload = {
    "schema": "elkies-k3.h3-a5a5-equation-cost-neighbors.v1",
    "status": "PASS_EXACT_A5A5_EQUATION_COST_SCORING",
    "inputs": {"paths": [str(path.relative_to(ROOT)) for path in INPUTS], "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in INPUTS}},
    "candidate_count": len(records),
    "ranked_candidates": records,
    "proof_boundary": "Exact curve intersections and lattice decompositions; RR and weighted totals are planning estimates. Full section-wall nefness is separate.",
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
best = records[0]
print("A5A5COST|best=o{}|score={}|PO={}|RR={}|deg0={}|deg1={}|pinned={}|child={}/MW{}|status={}".format(
    best["candidate_id"]["orbit_index"], best["equation_cost_score"], best["horizontal"]["P_dot_O"],
    best["expected_RR_ambient"], len(best["explicit_degree_zero_curves"]), len(best["explicit_degree_one_curves"]),
    best["pinned_R17_degree"], best["child"]["ade"], best["child"]["mw_rank"], payload["status"]), flush=True)
print(f"OUTPUT|{OUTPUT.resolve()}", flush=True)
