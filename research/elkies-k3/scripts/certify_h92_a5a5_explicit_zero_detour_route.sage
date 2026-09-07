#!/usr/bin/env sage -python
"""Certify and cost the explicit-zero A5A5 loop route to pinned R17."""

import hashlib
import json
from pathlib import Path

from sage.all import *
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
FINGERPRINT = ROOT / "artifacts/local/elkies-k3/q24-a11-q8-construction-fingerprint.json"
PINNED_SUFFIX = GENERATED / "elkies-k3-h3-pinned-r17-current-suffix-marking.json"
A11_EDGE = GENERATED / "elkies-k3-h3-a11-q8-orbit12-lattice-certificate.json"
A5_EDGE = GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q6-orbit3372-lattice-certificate.json"
START = GENERATED / "elkies-k3-h3-a5a5-q6o3372-suffix-marking.json"
FIRST_COST = GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4q6-equation-cost.json"
DIRECT_COST = GENERATED / "elkies-k3-h3-a5a5-current-route-equation-cost-audit.json"
OUTPUT = GENERATED / "elkies-k3-h3-a5a5-explicit-zero-detour-route-certificate.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))

EDGES = (
    ("return_current_A5A5", GENERATED / "elkies-k3-h3-a5a5-q6o3372-q6-pinned-best-certificate.json"),
    ("current_A5A5_to_current_3A3", GENERATED / "elkies-k3-h3-a5a5-q6o3372-q6-current-3a3-certificate.json"),
    ("current_3A3_to_current_A3_2A2", GENERATED / "elkies-k3-h3-detour-current-a3-2a2-certificate.json"),
    ("current_A3_2A2_to_current_5A1", GENERATED / "elkies-k3-h3-detour-current-5a1-certificate.json"),
    ("current_5A1_to_current_4A1", GENERATED / "elkies-k3-h3-detour-current-4a1-certificate.json"),
    ("current_4A1_to_current_3A1", GENERATED / "elkies-k3-h3-detour-current-3a1-certificate.json"),
    ("current_3A1_to_current_2A1", GENERATED / "elkies-k3-h3-detour-current-2a1-certificate.json"),
    ("current_2A1_to_current_A1", GENERATED / "elkies-k3-h3-detour-current-a1-certificate.json"),
    ("current_A1_to_pinned_R17", GENERATED / "elkies-k3-h3-detour-pinned-r17-certificate.json"),
)


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


def root_components(cartan):
    unseen = set(range(cartan.nrows()))
    answer = []
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
        answer.append(tuple(sorted(component)))
    return tuple(answer)


def ade_label(cartan):
    labels = []
    for component in root_components(cartan):
        sub = cartan.matrix_from_rows_and_columns(component, component)
        rank = len(component)
        degrees = [sum(sub[i, j] == -1 for j in range(rank)) for i in range(rank)]
        determinant = abs(sub.det())
        if max(degrees, default=0) <= 2 and determinant == rank + 1:
            labels.append("A{}".format(rank))
        elif determinant == 4 and rank >= 4:
            labels.append("D{}".format(rank))
        elif determinant == 3 and rank == 6:
            labels.append("E6")
        elif determinant == 2 and rank == 7:
            labels.append("E7")
        elif determinant == 1 and rank == 8:
            labels.append("E8")
        else:
            labels.append("rank{}_det{}".format(rank, determinant))
    return "+".join(sorted(labels)) if labels else "rootless"


def vertical_layers(coefficients, cartan):
    edges = [
        (i, j) for i in range(cartan.nrows()) for j in range(i + 1, cartan.nrows())
        if cartan[i, j] == -1
    ]
    magnitudes = [abs(ZZ(value)) for value in coefficients]
    previous = total = 0
    for level in sorted(set(value for value in magnitudes if value)):
        active = set(index for index, value in enumerate(magnitudes) if value >= level)
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
        total += (level - previous) * count
        previous = level
    return int(total)


def rr_profile(fibre, frame, root_rank):
    root = frame[:root_rank, :root_rank]
    coupling = frame[:root_rank, root_rank:]
    base = vector(ZZ, [0] * root_rank + list(fibre[2 + root_rank:]))
    if root_rank:
        dual = vector(QQ, base * frame[:, :root_rank]) * root.inverse()
        iterator = IntegralLattice(root).enumerate_close_vectors(-dual)
    else:
        iterator = iter((vector(ZZ, []),))
    choices = []
    minimum = None
    for shift in iterator:
        lifted = base + vector(ZZ, list(shift) + [0] * (17 - root_rank))
        norm = QQ(lifted * frame * lifted)
        if minimum is None:
            minimum = norm
        elif norm > minimum:
            break
        pole = (norm - 4) / 2
        if pole in ZZ and pole >= 0:
            choices.append((ZZ(pole), lifted))
    assert choices
    pole, lifted = min(choices, key=lambda item: (item[0], tuple(item[1])))
    section = vector(ZZ, [pole + 1, 1] + list(lifted))
    old_zero = vector(ZZ, [-1, 1] + [0] * 17)
    vertical = vector(ZZ, (fibre - old_zero - section)[2:2 + root_rank])
    layers = vertical_layers(vertical, root)
    return {
        "P_dot_O": int(pole),
        "section": entries(section),
        "vertical": entries(vertical),
        "vertical_layers": layers,
        "vertical_support": sum(value != 0 for value in vertical),
        "expected_RR_ambient": 2 + 2 * int(pole) + layers,
    }


start = json.loads(START.read_text())
a11_edge = json.loads(A11_EDGE.read_text())
a5_edge = json.loads(A5_EDGE.read_text())
assert start["status"] == "PASS_EXACT_A5A5_CANDIDATE_SUFFIX_MARKING"
assert a11_edge["status"] == "PASS_EXACT_A11_Q8_EQUATION_COST_LATTICE_CERTIFICATE"
assert a5_edge["status"] == "PASS_EXACT_A5A5_EXPLICIT_ZERO_CANDIDATE_LATTICE_CERTIFICATE"

start_frame_path = ROOT / start["frame_output"]
start_frame = load_matrix(start_frame_path)
g_start = block_diagonal_matrix(U2, -start_frame)
equation_to_start = matrix(ZZ, start["equation_A11_to_root_adapted_hub_basis"])
assert abs(equation_to_start.det()) == 1

cumulative = identity_matrix(ZZ, 19)
source_frame = start_frame
source_root_rank = int(start["root_data"][0])
records = []
inputs = [A11_EDGE, A5_EDGE, START, FIRST_COST, DIRECT_COST, FINGERPRINT, PINNED_SUFFIX, start_frame_path, PINNED_FRAME]
for label, path in EDGES:
    certificate = json.loads(path.read_text())
    assert certificate["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
    assert certificate["first_edge_nef_audit"]["nef_in_selected_component_chamber"]
    transition = matrix(ZZ, certificate["source_to_child_basis"])
    inverse = matrix(ZZ, certificate["child_to_source_basis"])
    assert transition * inverse == identity_matrix(ZZ, 19)
    g_source = block_diagonal_matrix(U2, -source_frame)
    fibre = vector(ZZ, transition.row(0))
    assert fibre * g_source * fibre == 0 and gcd(tuple(g_source * fibre)) == 1
    child_frame_path = ROOT / certificate["frame_output"]
    child_frame = load_matrix(child_frame_path)
    g_child = block_diagonal_matrix(U2, -child_frame)
    assert transition * g_source * transition.transpose() == g_child
    cumulative = transition * cumulative
    assert matrix(ZZ, certificate["equation_A11_to_child_basis"]) == cumulative * equation_to_start

    nef = certificate["first_edge_nef_audit"]
    component_degrees = nef["component_pairings"] + nef["affine_pairings"]
    profile = rr_profile(fibre, source_frame, source_root_rank)
    child_root_rank = int(certificate["child"]["root_data"][0])
    child_cartan = child_frame[:child_root_rank, :child_root_rank]
    terms = {
        "P_dot_O": 900 * profile["P_dot_O"],
        "horizontal_degree": 250 * int(fibre[1]),
        "RR_ambient": 120 * profile["expected_RR_ambient"],
        "vertical_layers": 60 * profile["vertical_layers"],
        "vertical_support": 25 * profile["vertical_support"],
        "child_root_count": int(certificate["child"]["root_data"][1]),
        "coordinate_growth": max(abs(int(value)) for value in fibre),
        "no_explicit_degree_one_component": 4000 if 1 not in component_degrees else 0,
        "explicit_degree_one_component_credit": -500 * min(component_degrees.count(1), 6),
        "explicit_degree_zero_component_credit": -100 * min(component_degrees.count(0), 12),
    }
    records.append({
        "label": label,
        "certificate": str(path.relative_to(ROOT)),
        "q": int(fibre[0] * fibre[1]),
        "old_fibre_degree": int(fibre[1]),
        "fibre_in_source": entries(fibre),
        "nef_audit": nef,
        "source_root_rank": source_root_rank,
        "child_root_data": certificate["child"]["root_data"],
        "child_ADE": ade_label(child_cartan),
        "child_MW_rank": 17 - child_root_rank,
        "horizontal_RR_profile": profile,
        "explicit_component_degree_zero_count": component_degrees.count(0),
        "explicit_component_degree_one_count": component_degrees.count(1),
        "equation_cost_terms": terms,
        "equation_cost_score": int(sum(terms.values())),
        "transport_determinant": int(transition.det()),
    })
    inputs.extend((path, child_frame_path))
    source_frame = child_frame
    source_root_rank = child_root_rank

# Exact endpoint identification against the repository's pinned full basis.
g_final = block_diagonal_matrix(U2, -source_frame)
pinned_suffix = json.loads(PINNED_SUFFIX.read_text())
fingerprint = json.loads(FINGERPRINT.read_text())
assert pinned_suffix["status"] == "PASS_EXACT_PINNED_R17_CURRENT_SUFFIX_MARKING"
historical_in_equation = block_diagonal_matrix(
    identity_matrix(ZZ, 2),
    matrix(ZZ, fingerprint["selected"]["frame_isometry_historical_basis_in_equation_coordinates"]),
)
pinned_in_historical = matrix(
    ZZ, pinned_suffix["current_suffix_stages"]["current_A11"]["pinned_R17_basis_in_stage"]
)
pinned_basis_equation = pinned_in_historical * historical_in_equation
pinned_basis_start = pinned_basis_equation * equation_to_start.inverse().change_ring(ZZ)
pinned_basis_final = pinned_basis_start * cumulative.inverse().change_ring(ZZ)
pinned = load_matrix(PINNED_FRAME)
g_pinned = block_diagonal_matrix(U2, -pinned)
assert abs(pinned_basis_final.det()) == 1
assert pinned_basis_final * g_final * pinned_basis_final.transpose() == g_pinned
assert vector(ZZ, pinned_basis_final.row(0)) == vector(ZZ, [1, 0] + [0] * 17)
assert records[-1]["child_root_data"] == [0, 0, 1]

first_cost = json.loads(FIRST_COST.read_text())
first = next(
    record for record in first_cost["ranked_candidates"]
    if record["candidate_id"] == {"q": 6, "old_fibre_degree": 2, "orbit_index": 3372}
)
direct = json.loads(DIRECT_COST.read_text())["direct_equation_cost_profile"]
direct_zero = sum(value == 0 for value in direct["named_explicit_curve_degrees"].values())
direct_one = sum(value == 1 for value in direct["named_explicit_curve_degrees"].values())
direct_terms = {
    "P_dot_O": 900 * direct["P_dot_O"],
    "horizontal_degree": 250 * 2,
    "RR_ambient": 120 * direct["expected_RR_ambient"],
    "vertical_layers": 60 * direct["vertical_layers"],
    "vertical_support": 25 * sum(value != 0 for value in direct["vertical"]),
    "child_root_count": 36,
    "coordinate_growth": max(abs(value) for value in json.loads(DIRECT_COST.read_text())["explicit_zero_edge"]["dominant_fibre"]),
    "no_explicit_degree_one_curve": 4000 if not direct_one else 0,
    "explicit_degree_one_credit": -500 * min(direct_one, 6),
    "explicit_degree_zero_credit": -100 * min(direct_zero, 12),
}

loop_records = [first] + records[:2]
payload = {
    "schema": "elkies-k3.h3-a5a5-explicit-zero-detour-route-certificate.v1",
    "status": "PASS_EXACT_LATTICE_CERTIFIED_EXPLICIT_ZERO_DETOUR_TO_PINNED_R17",
    "promotion": {
        "promote_as_lifting_target": False,
        "reason": (
            "The route is fully lattice-certified and avoids the direct q104 "
            "explicit-curve gate failure, but its current deterministic raw cost "
            "score is not yet strictly below the q104 comparator."
        ),
    },
    "prefix": {
        "A11_to_2A5": {"q": 8, "old_fibre_degree": 2, "orbit": 12, "certificate": str(A11_EDGE.relative_to(ROOT))},
        "2A5_to_A1_A1_A3_A5": {"q": 6, "old_fibre_degree": 2, "orbit": 3372, "certificate": str(A5_EDGE.relative_to(ROOT))},
    },
    "detour_and_suffix_edges": records,
    "zero_changing_loop": {
        "edges": ["2A5_to_A1_A1_A3_A5", records[0]["label"], records[1]["label"]],
        "q_sequence": [6, 6, 12],
        "old_fibre_degrees": [2, 2, 2],
        "equation_cost_score": int(first["equation_cost_score"] + sum(item["equation_cost_score"] for item in records[:2])),
        "direct_q104_comparator_score": int(sum(direct_terms.values())),
        "direct_q104_comparator_terms": direct_terms,
        "direct_q104_negative_named_explicit_curves": [
            name for name, value in direct["named_explicit_curve_degrees"].items() if value < 0
        ],
        "direct_q104_explicit_curve_gate": "REJECT",
        "detour_exact_nef_gate": "PASS",
        "strict_score_improvement": first["equation_cost_score"] + sum(item["equation_cost_score"] for item in records[:2]) < sum(direct_terms.values()),
    },
    "full_route_q_sequence_from_A11": [8, 6] + [item["q"] for item in records],
    "full_route_old_fibre_degrees_from_A11": [2] * (2 + len(records)),
    "full_route_equation_cost_score_after_A11": int(first["equation_cost_score"] + sum(item["equation_cost_score"] for item in records)),
    "endpoint": {
        "name": "pinned_R17",
        "root_data": records[-1]["child_root_data"],
        "ADE": records[-1]["child_ADE"],
        "MW_rank": records[-1]["child_MW_rank"],
        "canonical_pinned_basis_in_final_child": rows(pinned_basis_final),
        "final_child_basis_in_canonical_pinned": rows(pinned_basis_final.inverse().change_ring(ZZ)),
        "forward_determinant": int(pinned_basis_final.det()),
        "inverse_determinant": int(pinned_basis_final.inverse().det()),
        "canonical_fibre_equals_final_marked_fibre": True,
        "gram_identification": "U plus negative pinned rank17_gram.txt exactly",
    },
    "proof_boundary": (
        "Every edge is an exact primitive nef isotropic class with exact component and "
        "all-section gates, marked U, full determinant-one NS transport in both directions, "
        "and exact root/MW data. The endpoint full basis is identified with pinned R17. "
        "Equation-cost scores are deterministic compiler estimates, not completed equation lifts."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in dict.fromkeys(inputs)],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in dict.fromkeys(inputs)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A5DETOUR|edges={}|q={}|loop_score={}|direct_q104_score={}|strict_improvement={}|"
    "endpoint={}/MW{}|det={}|status={}|output={}".format(
        len(payload["full_route_q_sequence_from_A11"]),
        ",".join(map(str, payload["full_route_q_sequence_from_A11"])),
        payload["zero_changing_loop"]["equation_cost_score"],
        payload["zero_changing_loop"]["direct_q104_comparator_score"],
        int(payload["zero_changing_loop"]["strict_score_improvement"]),
        payload["endpoint"]["ADE"], payload["endpoint"]["MW_rank"],
        payload["endpoint"]["forward_determinant"], payload["status"], OUTPUT,
    ),
    flush=True,
)
