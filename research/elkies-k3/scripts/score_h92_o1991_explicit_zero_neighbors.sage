#!/usr/bin/env sage -python
"""Score q4..q10 neighbours of the explicit-zero orbit1991 MW4 frame.

All intersections are exact.  The RR size and weighted total are planning
estimates.  Candidate nefness is screened against every transported exact
curve and every physical affine root component; full section-wall nefness is
left to the retained-candidate certifier.
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
parser.add_argument(
    "--output",
    type=Path,
    default=GENERATED / "elkies-k3-h3-o1991-explicit-zero-equation-cost-neighbors.json",
)
parser.add_argument("--retain", type=int, default=300)
args = parser.parse_args()

ZERO_FRAMES = GENERATED / "elkies-k3-h3-a11-q8-orbit1991-explicit-zero-frames.json"
NEIGHBORS = GENERATED / "elkies-k3-h3-a11-o1991-explicit-zero-degree2-neighbors.json"
D12_FRAME = LOCAL / "q24-downstream-lift/d12-c10a-zero-frame.txt"
Q6 = LOCAL / "q24-downstream-lift/d12-c10a-zero-q6-all.json"
IDENTITY = LOCAL / "q24-orbit42-identity-halving-audit.json"
MATCHING = LOCAL / "q24-orbit42-identity-halving-qq.json"
ZERO_MISMATCH = GENERATED / "elkies-k3-h3-a11-quintic-bridge-zero-mismatch.json"
INPUTS = (ZERO_FRAMES, NEIGHBORS, D12_FRAME, Q6, IDENTITY, MATCHING, ZERO_MISMATCH)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")


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


def graph_edges(cartan):
    return [
        (left, right)
        for left in range(cartan.nrows())
        for right in range(left + 1, cartan.nrows())
        if cartan[left, right] == -1
    ]


def connected_count(edges, active):
    active = set(active)
    count = 0
    while active:
        count += 1
        todo = [active.pop()]
        while todo:
            node = todo.pop()
            for left, right in edges:
                other = right if left == node else left if right == node else None
                if other in active:
                    active.remove(other)
                    todo.append(other)
    return count


def vertical_layers(coefficients, edges):
    magnitudes = [abs(ZZ(item)) for item in coefficients]
    total = previous = ZZ(0)
    for level in sorted(set(item for item in magnitudes if item)):
        active = [index for index, item in enumerate(magnitudes) if item >= level]
        total += (level - previous) * connected_count(edges, active)
        previous = level
    return int(total)


zero_frames = json.loads(ZERO_FRAMES.read_text())
neighbors = json.loads(NEIGHBORS.read_text())
q6 = json.loads(Q6.read_text())
identity = json.loads(IDENTITY.read_text())
matching = json.loads(MATCHING.read_text())
zero_mismatch = json.loads(ZERO_MISMATCH.read_text())
assert zero_frames["status"] == "PASS_EXACT_A11_Q8_ORBIT1991_EXPLICIT_ZERO_FRAMES"
assert neighbors["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert zero_mismatch["status"] == "REJECT_A11_QUINTIC_BRIDGE_ZERO_MISMATCH"

frame = load_matrix(ROOT / neighbors["frame"])
root_rank = 13
root = frame[:root_rank, :root_rank]
coupling = frame[:root_rank, root_rank:]
tail = frame[root_rank:, root_rank:]
height = tail - coupling.transpose() * root.inverse() * coupling
root_lattice = IntegralLattice(root)
edges = graph_edges(root)
highest = highest_roots(root)
g_parent = block_diagonal_matrix(U2, -frame)

# Reconstruct exact equation-A11 shell sections.
d12 = load_matrix(D12_FRAME)
selected_q6 = next(item for item in q6["neighbors"] if int(item["orbit_index"]) == 64)
d12_to_a11 = block_diagonal_matrix(
    identity_matrix(ZZ, 2), matrix(ZZ, selected_q6["child_root_adapted_basis"])
) * matrix(ZZ, selected_q6["neighbor_basis"])
d12_inverse = d12_to_a11.inverse().change_ring(ZZ)
d12_root = d12[:12, :12]
d12_coupling = d12[:12, 12:]
abstract_sections = []
for values in identity["exact_model_R3_zero"]["identity_vectors"]:
    z = vector(ZZ, values)
    root_coefficients = -(z * d12_coupling.transpose()) * d12_root.inverse()
    section = vector(ZZ, [1, 1] + list(map(ZZ, root_coefficients)) + list(z))
    abstract_sections.append(section * d12_inverse)
shell_mapping = matching["matching"]["mappings_abstract_to_equation"][7]
a11_shell = [None] * 18
for abstract_index, equation_index in enumerate(shell_mapping):
    a11_shell[equation_index] = abstract_sections[abstract_index]
assert all(item is not None for item in a11_shell)

# Add the full old A11 fibre, old zero, and the correctly transported A0/P24
# curves.  These are all actual characteristic-zero curves, not abstract roots.
g_a11 = block_diagonal_matrix(U2, -matrix(ZZ, selected_q6["child_root_adapted_frame"]))
old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
a11_simple = [
    vector(ZZ, [0, 0] + [-ZZ(index == node) for index in range(17)])
    for node in range(11)
]
a11_affine = old_fibre + vector(ZZ, [0, 0] + list(highest_roots(matrix(ZZ, selected_q6["child_root_adapted_frame"])[:11, :11])[0]) + [0] * 6)
correct = zero_mismatch["correct_selected_R3_transport"]
a0 = vector(ZZ, correct["oldI9_A0"]["child_coordinates"])
p24 = vector(ZZ, correct["close_P24"]["child_coordinates"])
named_a11_curves = (
    [(f"shell_S{index}", item) for index, item in enumerate(a11_shell)]
    + [("old_A11_zero", old_zero), ("old_A11_affine", a11_affine)]
    + [(f"old_A11_component_{index}", item) for index, item in enumerate(a11_simple)]
    + [("oldI9_A0", a0), ("close_P24", p24)]
)
assert all(item * g_a11 * item == -2 for unused, item in named_a11_curves)

equation_to_parent = matrix(ZZ, zero_frames["selected"]["equation_A11_to_explicit_zero_basis"])
parent_inverse = equation_to_parent.inverse().change_ring(ZZ)
exact_curves = [(name, item * parent_inverse) for name, item in named_a11_curves]
assert all(item * g_parent * item == -2 for unused, item in exact_curves)
exact_sections = [(name, item) for name, item in exact_curves if item[1] == 1]
zero_hits = [name for name, item in exact_curves if item == vector(ZZ, [-1, 1] + [0] * 17)]
assert len(zero_hits) == 1
section_mw = [vector(ZZ, item[-4:]) for unused, item in exact_sections if any(item[-4:])]
known_section_lattice = (
    matrix(ZZ, section_mw).row_module() if section_mw else matrix(ZZ, 0, 4).row_module()
)

profile_cache = {}


def section_profiles(z):
    z = vector(ZZ, z)
    key = tuple(z)
    if key in profile_cache:
        return profile_cache[key]
    horizontal_height = QQ(z * height * z)
    base = vector(ZZ, [0] * root_rank + list(z))
    dual = vector(QQ, base * frame[:, :root_rank]) * root.inverse()
    iterator = root_lattice.enumerate_close_vectors(-dual)
    minimum = None
    result = []
    for unused in range(100000):
        shift = vector(ZZ, next(iterator))
        lifted = base + vector(ZZ, list(shift) + [0] * 4)
        norm = QQ(lifted * frame * lifted)
        if minimum is None:
            minimum = norm
        elif norm > minimum:
            break
        pole = (norm - 4) / 2
        if pole in ZZ and pole >= 0:
            result.append((lifted, horizontal_height, norm - horizontal_height, ZZ(pole)))
    assert result
    profile_cache[key] = result
    return result


records = []
for raw in neighbors["neighbors"]:
    fibre = vector(ZZ, raw["fiber"])
    z = vector(ZZ, raw["mw_projection"])
    degree = int(raw["old_fiber_degree"])
    profiles = []
    for lifted, horizontal_height, correction, pole in section_profiles(z):
        section = vector(ZZ, [pole + 1, 1] + list(lifted))
        residual = fibre - (degree - 1) * vector(ZZ, [-1, 1] + [0] * 17) - section
        assert residual[1] == 0 and not any(residual[2 + root_rank:])
        vertical = vector(ZZ, residual[2:2 + root_rank])
        profiles.append(
            {
                "height": str(horizontal_height),
                "local_correction": str(correction),
                "P_dot_O": int(pole),
                "section": entries(section),
                "vertical": entries(vertical),
                "fibre_twist": int(residual[0]),
                "vertical_support": sum(item != 0 for item in vertical),
                "vertical_L1": int(sum(abs(item) for item in vertical)),
                "vertical_layers": vertical_layers(vertical, edges),
            }
        )
    profiles.sort(key=lambda item: (item["P_dot_O"], item["vertical_layers"], item["vertical_support"], item["vertical_L1"], tuple(item["section"])))
    horizontal = profiles[0]

    curve_degrees = {name: int(item * g_parent * fibre) for name, item in exact_curves}
    negative_curves = sorted(name for name, value in curve_degrees.items() if value < 0)
    explicit_degree_zero = sorted(name for name, value in curve_degrees.items() if value == 0)
    explicit_degree_one = sorted(name for name, value in curve_degrees.items() if value == 1)
    labels = vector(ZZ, raw["dominant_labels"])
    affine_pairings = [int(degree - highest_root * labels) for highest_root in highest]
    negative_affine = [index for index, value in enumerate(affine_pairings) if value < 0]
    declared_nef = not negative_curves and not negative_affine

    in_known = bool(known_section_lattice.rank() and z in known_section_lattice)
    rank_gap = 0 if in_known else 4 - int(known_section_lattice.rank())
    rr = 2 + 2 * horizontal["P_dot_O"] + horizontal["vertical_layers"]
    child_root_count = int(raw["child_root_data"][1])
    coordinate_growth = max(abs(int(item)) for item in fibre)
    terms = {
        "declared_non_nef_penalty": 1000000 if not declared_nef else 0,
        "unspanned_horizontal_rank_gap": 5000 * rank_gap,
        "P_dot_O": 900 * horizontal["P_dot_O"],
        "horizontal_degree": 250 * degree,
        "RR_ambient": 120 * rr,
        "vertical_layers": 60 * horizontal["vertical_layers"],
        "vertical_support": 25 * horizontal["vertical_support"],
        "child_root_count": child_root_count,
        "coordinate_growth": coordinate_growth,
        "no_explicit_degree_one_curve": 3000 if not explicit_degree_one else 0,
        "explicit_degree_one_credit": -500 * min(len(explicit_degree_one), 4),
        "explicit_degree_zero_credit": -100 * min(len(explicit_degree_zero), 8),
    }
    records.append(
        {
            "candidate_id": {
                "q": int(raw["q"]),
                "old_fibre_degree": degree,
                "orbit_index": int(raw["orbit_index"]),
            },
            "declared_curve_and_affine_nef_gate": "PASS" if declared_nef else "REJECT",
            "child": {
                "ade": raw["child_ade"],
                "mw_rank": int(raw["child_mw_rank"]),
                "root_data": raw["child_root_data"],
            },
            "fibre": entries(fibre),
            "mw_projection": entries(z),
            "horizontal": horizontal,
            "expected_RR_ambient": rr,
            "exact_curve_degrees": curve_degrees,
            "negative_exact_curves": negative_curves,
            "explicit_degree_zero_curves": explicit_degree_zero,
            "explicit_degree_one_curves": explicit_degree_one,
            "affine_component_pairings": affine_pairings,
            "negative_affine_components": negative_affine,
            "known_section_subgroup": {
                "rank": int(known_section_lattice.rank()),
                "target_in_subgroup": in_known,
                "unspanned_rank_gap": rank_gap,
            },
            "equation_cost_terms": terms,
            "equation_cost_score": sum(terms.values()),
        }
    )

records.sort(
    key=lambda item: (
        item["equation_cost_score"],
        item["horizontal"]["P_dot_O"],
        -len(item["explicit_degree_one_curves"]),
        item["expected_RR_ambient"],
        item["candidate_id"]["q"],
        item["candidate_id"]["orbit_index"],
    )
)
retained = records[:args.retain]
payload = {
    "schema": "elkies-k3.h3-o1991-explicit-zero-equation-cost-neighbors.v1",
    "status": "PASS_EXACT_O1991_EXPLICIT_ZERO_EQUATION_COST_SCORING",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in INPUTS},
    },
    "parent": {
        "ade": "D10+A1+A2",
        "mw_rank": 4,
        "exact_curve_count": len(exact_curves),
        "exact_section_names": [name for name, unused in exact_sections],
        "exact_nonzero_section_subgroup_rank": int(known_section_lattice.rank()),
    },
    "search_summaries": neighbors["summaries"],
    "candidate_count": len(records),
    "retained_count": len(retained),
    "best_candidate": retained[0],
    "retained_candidates": retained,
    "proof_boundary": (
        "Exact stored-shell scoring, horizontal closest-root lifts, transported QQ "
        "curve intersections, and physical affine-component rejects. The RR total "
        "is estimated, and PASS candidates still require the full section-wall "
        "nef test plus independent transport/endpoint certification."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
best = retained[0]
print(
    "O1991EQCOST|candidates={}|best=q{}o{}|child={}/MW{}|PO={}|RR={}|"
    "deg0={}|deg1={}|known_rank={}|score={}|status={}".format(
        len(records), best["candidate_id"]["q"], best["candidate_id"]["orbit_index"],
        best["child"]["ade"], best["child"]["mw_rank"], best["horizontal"]["P_dot_O"],
        best["expected_RR_ambient"], len(best["explicit_degree_zero_curves"]),
        len(best["explicit_degree_one_curves"]), best["known_section_subgroup"]["rank"],
        best["equation_cost_score"], payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
