#!/usr/bin/env sage -python
"""Certify an explicit degree-two neighbor of an exact marked H92 frame.

status: ACTIVE_PROOF
claim: marked U, exact nef audit, roots, and bidirectional unimodular transport
inputs: a marked source frame, candidate fibre, marking and expected target
outputs: caller-selected exact JSON certificate and child frame
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def rational_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    if count == 0:
        return (), matrix(ZZ, 0, gram.nrows()), (0, 0, 1)
    half = [vector(ZZ, column) for column in matrix(ZZ, result[2]).columns()]
    roots = tuple(half + [-root for root in half])
    basis = matrix(ZZ, [list(root) for root in roots]).row_module().basis_matrix()
    root_gram = basis * gram * basis.transpose()
    return roots, basis, (basis.rank(), count, abs(ZZ(root_gram.det())))


def deterministic_simple_roots(gram):
    roots, _, data = roots_and_data(gram)
    positive = [root for root in roots if next(value for value in root if value) > 0]
    positive_set = {tuple(root) for root in positive}
    simple = [
        root for root in positive
        if not any(tuple(root - left) in positive_set for left in positive)
    ]
    result = matrix(ZZ, [list(root) for root in simple])
    assert result.nrows() == result.rank() == data[0]
    return result, result * gram * result.transpose()


def root_adaptation(child):
    _, root_basis, data = roots_and_data(child)
    root_rank = data[0]
    if root_rank == 0:
        scale = ZZ(1)
        lll = matrix(ZZ, pari((scale * child).change_ring(ZZ)).qflllgram())
        basis = lll.transpose()
        adapted = basis * child * basis.transpose()
        return adapted, basis, adapted.change_ring(QQ), data
    smith, _, smith_right = root_basis.smith_form()
    assert tuple(abs(smith[index, index]) for index in range(root_rank)) == (1,) * root_rank
    simple, cartan = deterministic_simple_roots(child)
    completion = smith_right.inverse()
    basis = simple.stack(completion[root_rank:])
    assert abs(basis.det()) == 1
    adapted = basis * child * basis.transpose()
    coupling = adapted[:root_rank, root_rank:]
    tail = adapted[root_rank:, root_rank:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling
    scale = lcm(entry.denominator() for entry in height.list())
    lll = matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram())
    quotient = block_diagonal_matrix(identity_matrix(ZZ, root_rank), lll.transpose())
    basis = quotient * basis
    adapted = basis * child * basis.transpose()
    coupling = adapted[:root_rank, root_rank:]
    tail = adapted[root_rank:, root_rank:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling
    return adapted, basis, height, data


def bezout_vector_for_pairing(ns, fibre):
    current = ZZ(0)
    result = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(ns * fibre):
        if value == 0:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        result = [left * entry for entry in result]
        result[index] += right
        current = divisor
    assert abs(current) == 1
    if current == -1:
        result = [-entry for entry in result]
    return vector(ZZ, result)


def child_frame(ns, fibre):
    mate = bezout_vector_for_pairing(ns, fibre)
    mate -= ZZ(mate * ns * mate // 2) * fibre
    kernel = matrix(ZZ, [list(fibre * ns), list(mate * ns)]).right_kernel_matrix()
    child = -(kernel * ns * kernel.transpose())
    transition = matrix(ZZ, [list(fibre), list(mate)] + list(kernel.rows()))
    assert abs(transition.det()) == 1
    assert transition * ns * transition.transpose() == block_diagonal_matrix(U2, -child)
    return child, transition


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
    if cartan.nrows() == 0:
        return ()
    half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
    roots = tuple(half) + tuple(-item for item in half)
    answer = []
    for component in components(cartan):
        candidates = [
            root for root in roots
            if all(value >= 0 for value in root)
            and all(index in component or root[index] == 0 for index in range(cartan.nrows()))
        ]
        answer.append(max(candidates, key=lambda root: sum(root)))
    return tuple(answer)


def nef_profile(fibre, frame, root_rank):
    degree = ZZ(fibre[1])
    labels = vector(ZZ, fibre[2:]) * frame[:, :root_rank]
    cartan = frame[:root_rank, :root_rank]
    affine = [degree - top * labels for top in highest_roots(cartan)]
    if degree == 0:
        assert vector(ZZ, fibre[2:]).is_zero() and fibre[0] > 0
        return {
            "component_pairings": list(map(int, labels)),
            "affine_pairings": list(map(int, affine)),
            "minimum_section_distance": None,
            "minimum_section_intersection": str(fibre[0]),
            "nef_in_selected_component_chamber": True,
            "same_fibre_ray": True,
        }
    center = vector(QQ, fibre[2:]) / degree
    closest = vector(ZZ, next(IntegralLattice(frame).enumerate_close_vectors(center)))
    distance = (closest - center) * frame * (closest - center)
    minimum_section = degree * (distance - 2) / 2
    return {
        "component_pairings": list(map(int, labels)),
        "affine_pairings": list(map(int, affine)),
        "minimum_section_distance": str(distance),
        "minimum_section_intersection": str(minimum_section),
        "nef_in_selected_component_chamber": bool(
            min(tuple(labels) + tuple(affine) + (ZZ(0),)) >= 0
            and minimum_section >= 0
        ),
    }


def negative_horizontal_walls(fibre, frame):
    """Enumerate every negative old-horizontal (-2)-curve wall exactly."""
    degree = ZZ(fibre[1])
    w = vector(ZZ, fibre[2:])
    walls = []
    for old_degree in range(1, int(degree) + 1):
        m = ZZ(old_degree)
        cross = -degree * m * frame * w.column()
        augmented = block_matrix(ZZ, [
            [degree**2 * frame, cross],
            [cross.transpose(), matrix(ZZ, [[m**2 * (w * frame * w) + 1]])],
        ])
        result = pari(augmented).qfminim(2 * degree**2 - 1)
        normalized = set()
        for raw in matrix(ZZ, result[2]).transpose().rows():
            if abs(raw[-1]) != 1:
                continue
            value = raw if raw[-1] == 1 else -raw
            normalized.add(tuple(value))
        for value in normalized:
            x = vector(ZZ, value[:-1])
            x_norm = ZZ(x * frame * x)
            if (x_norm - 2) % (2 * m):
                continue
            k = ZZ((x_norm - 2) // (2 * m))
            intersection = ZZ(
                (w * frame * w // (2 * degree)) * m
                + degree * k - w * frame * x
            )
            if intersection < 0:
                walls.append({
                    "old_fibre_degree": int(m),
                    "curve": [int(k), int(m)] + list(map(int, x)),
                    "intersection": int(intersection),
                })
    return sorted(walls, key=lambda item: (item["intersection"], item["curve"]))


def reduce_component_chamber(fibre, frame, root_rank):
    if root_rank == 0:
        return vector(ZZ, fibre), identity_matrix(ZZ, 19), []
    g = block_diagonal_matrix(U2, -frame)
    cartan = frame[:root_rank, :root_rank]
    wall_roots = [
        vector(ZZ, [0, 0] + [-ZZ(index == other) for other in range(17)])
        for index in range(root_rank)
    ]
    wall_names = [f"simple_{index}" for index in range(root_rank)]
    for component_index, top in enumerate(highest_roots(cartan)):
        wall_roots.append(vector(ZZ, [1, 0] + list(top) + [0] * (17 - root_rank)))
        wall_names.append(f"affine_{component_index}")
    value = vector(ZZ, fibre)
    action = identity_matrix(ZZ, 19)
    reflections = []
    for unused in range(100000):
        pairings = [ZZ(value * g * root) for root in wall_roots]
        negative = next((index for index, pairing in enumerate(pairings) if pairing < 0), None)
        if negative is None:
            break
        root = wall_roots[negative]
        pairing = pairings[negative]
        step = identity_matrix(ZZ, 19) + (g * root.column()) * matrix(ZZ, [list(root)])
        value *= step
        action *= step
        reflections.append({"wall": wall_names[negative], "pairing": int(pairing)})
    else:
        raise RuntimeError("affine-Weyl chamber reduction did not terminate")
    assert action * g * action.transpose() == g
    return value, action, reflections


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source-marking", type=Path, required=True)
parser.add_argument("--source-frame", type=Path, required=True)
parser.add_argument("--fibre", required=True,
                    help="comma-separated 19-coordinate fibre in the source frame")
parser.add_argument("--candidate-label", required=True)
parser.add_argument("--target", default="orbit12")
parser.add_argument("--frame-output", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
SOURCE_MARKING = args.source_marking.resolve()
SOURCE_FRAME = args.source_frame.resolve()
FRAME_OUTPUT = args.frame_output.resolve()
OUTPUT = args.output.resolve()

marking = json.loads(SOURCE_MARKING.read_text())
assert marking["status"] in {
    "PASS_EXACT_REVERSE_HUB_EQUATION_MARKING",
    "PASS_EXACT_PINNED_R17_TARGETED_CANDIDATE_CERTIFICATE",
    "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE",
    "PASS_EXACT_MARKED_FRONTIER_CANDIDATE_CHECKPOINT",
    "PASS_EXACT_A5A5_CANDIDATE_SUFFIX_MARKING",
    "PASS_EXACT_A5A5_PHYSICAL_SOURCE_MARKING",
    "PASS_EXACT_A5A5_PHYSICAL_COMPONENT_CHAMBER_MARKING",
    "PASS_EXACT_Q4O208_PHYSICAL_3A3_MARKING",
    "PASS_EXACT_Q4O1584_PHYSICAL_EFFECTIVE_ZERO_MARKING",
    "PASS_EXACT_Q4O164_PHYSICAL_EFFECTIVE_ZERO_MARKING",
    "PASS_EXACT_CORRECTED_A3_2A2_PHYSICAL_EFFECTIVE_ZERO_MARKING",
    "PASS_EXACT_PHYSICAL_Q8_5A1_EFFECTIVE_ZERO_MARKING",
    "PASS_EXACT_PHYSICAL_AN_EFFECTIVE_ZERO_MARKING",
    "PASS_EXACT_PHYSICAL_A2_EFFECTIVE_ZERO_MARKING",
}
source = load_matrix(SOURCE_FRAME)
g_source = block_diagonal_matrix(U2, -source)
fibre = vector(ZZ, [ZZ(value) for value in args.fibre.split(",")])
assert len(fibre) == 19 and fibre[1] > 0 and fibre * g_source * fibre == 0
assert gcd(tuple(g_source * fibre)) == 1
source_root_rank = int(
    marking["root_data"][0]
    if "root_data" in marking
    else marking.get("child", {}).get("root_data", [0])[0]
)
first_edge_nef = nef_profile(fibre, source, source_root_rank)
assert first_edge_nef["nef_in_selected_component_chamber"]
first_edge_horizontal_walls = negative_horizontal_walls(fibre, source)
assert not first_edge_horizontal_walls

target_key = (
    "target_fibres_in_root_adapted_hub"
    if "target_fibres_in_root_adapted_hub" in marking
    else "target_fibres_in_child"
)
targets_source = {
    name: vector(ZZ, value) for name, value in marking[target_key].items()
}
assert args.target in targets_source and "pinned_R17" in targets_source
degrees = {name: int(fibre * g_source * value) for name, value in targets_source.items()}

raw_child, raw_transition = child_frame(g_source, fibre)
adapted, adaptation, height, root_data = root_adaptation(raw_child)
initial_transition = block_diagonal_matrix(identity_matrix(ZZ, 2), adaptation) * raw_transition
g_child = block_diagonal_matrix(U2, -adapted)
assert initial_transition * g_source * initial_transition.transpose() == g_child
initial_inverse = initial_transition.inverse().change_ring(ZZ)
anchor_initial = targets_source["pinned_R17"] * initial_inverse
anchor, chamber_action, chamber_reflections = reduce_component_chamber(
    anchor_initial, adapted, root_data[0]
)
transition = chamber_action.inverse().change_ring(ZZ) * initial_transition
inverse = transition.inverse().change_ring(ZZ)
assert abs(transition.det()) == 1
assert transition * g_source * transition.transpose() == g_child
targets_child = {name: value * inverse for name, value in targets_source.items()}
assert targets_child["pinned_R17"] == anchor

source_basis_key = (
    "equation_A11_to_root_adapted_hub_basis"
    if "equation_A11_to_root_adapted_hub_basis" in marking
    else "equation_A11_to_child_basis"
)
equation_to_source = matrix(ZZ, marking[source_basis_key])
equation_to_child = transition * equation_to_source
assert abs(equation_to_child.det()) == 1
profiles = {
    name: {
        "q": int(value[0] * value[1]),
        "old_fibre_degree": int(value[1]),
        "P_dot_O": int(value[0] - value[1]),
        "fibre_in_child": list(map(int, value)),
        "nef_audit": nef_profile(value, adapted, root_data[0]),
    }
    for name, value in targets_child.items()
}

certificate_status = (
    "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
    if fibre[1] == 2 else "PASS_EXACT_MARKED_ISOTROPIC_CANDIDATE_CERTIFICATE"
)
FRAME_OUTPUT.write_text(
    f"# exact marked degree-{fibre[1]} candidate child\n"
    + "\n".join(" ".join(map(str, row)) for row in adapted.rows()) + "\n"
)
inputs = (SOURCE_MARKING, SOURCE_FRAME)
payload = {
    "schema": (
        "elkies-k3.h3-marked-degree-two-candidate-certificate.v1"
        if fibre[1] == 2 else "elkies-k3.h3-marked-isotropic-candidate-certificate.v1"
    ),
    "status": certificate_status,
    "source_hub": marking.get("source_hub", marking.get("hub", "marked_child")),
    "candidate_id": {
        "label": args.candidate_label,
        "q": int(fibre[0] * fibre[1]),
        "old_fibre_degree": int(fibre[1]),
    },
    "selected_ranking_target": args.target,
    "marked_target_degrees": degrees,
    "first_edge_nef_audit": first_edge_nef,
    "first_edge_exact_negative_horizontal_walls": first_edge_horizontal_walls,
    "first_edge_exact_horizontal_nef_gate": True,
    "child": {
        "root_data": list(map(int, root_data)),
        "mw_rank": 17 - int(root_data[0]),
        "mw_height": rational_rows(height),
    },
    "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "source_to_child_basis": rows(transition),
    "child_to_source_basis": rows(inverse),
    "equation_A11_to_child_basis": rows(equation_to_child),
    "child_to_equation_A11_basis": rows(equation_to_child.inverse().change_ring(ZZ)),
    "target_fibres_in_child": {
        name: list(map(int, value)) for name, value in targets_child.items()
    },
    "target_profiles": profiles,
    "child_component_chamber": {
        "anchor": "pinned_R17",
        "anchor_nef_audit": profiles["pinned_R17"]["nef_audit"],
        "affine_weyl_reflections": chamber_reflections,
        "right_coordinate_action": rows(chamber_action),
    },
    "proof_boundary": (
        f"Exact primitive-nef degree-{fibre[1]} edge, including the complete finite "
        "old-horizontal wall audit, full child roots and MW rank, exact "
        "marked U, and determinant-one NS transports in both directions. Equation-A11 "
        "and pinned-R17 markings are composed exactly. A marked target profile is a "
        "second nef edge only when its recorded component/all-section audit passes."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in inputs
        },
    },
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
profile = profiles[args.target]
print(
    "MARKEDD2CERT|label={}|q={}|child_root={}|MW={}|target_degree={}|"
    "target_nef={}|det={}|status={}|output={}".format(
        args.candidate_label, fibre[0] * fibre[1], ",".join(map(str, root_data)),
        17 - root_data[0], degrees[args.target],
        int(profile["nef_audit"]["nef_in_selected_component_chamber"]),
        transition.det(), payload["status"], OUTPUT,
    ),
    flush=True,
)
