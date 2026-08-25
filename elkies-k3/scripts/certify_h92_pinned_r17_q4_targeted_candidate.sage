#!/usr/bin/env sage -python
"""Certify a targeted pinned-R17 neighbor and its marked hub transports."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
RANKING = GENERATED / "elkies-k3-h3-pinned-r17-q4-degree2-targeted-ranking.json"
MARKING = GENERATED / "elkies-k3-h3-pinned-r17-equation-marking.json"
CROSSOVERS = GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json"
SEMISTABLE_NEF = GENERATED / "elkies-k3-h3-semistable-mw2-reverse-suffix-nef.json"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
U2 = matrix(ZZ, ((0, 1), (1, 0)))

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--target", default="q25_mw7")
parser.add_argument("--rank", type=int, default=0)
parser.add_argument(
    "--fibre",
    help="optional comma-separated 19-coordinate pinned-R17 fibre",
)
parser.add_argument("--candidate-label", default="explicit")
parser.add_argument("--frame-output", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
FRAME_OUTPUT = args.frame_output.resolve()
OUTPUT = args.output.resolve()
INPUTS = (RANKING, MARKING, CROSSOVERS, SEMISTABLE_NEF, PINNED_FRAME)


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
    kernel = matrix(
        ZZ, [list(fibre * ns), list(mate * ns)]
    ).right_kernel_matrix()
    child = -(kernel * ns * kernel.transpose())
    transition = matrix(ZZ, [list(fibre), list(mate)] + list(kernel.rows()))
    assert abs(transition.det()) == 1
    assert transition * ns * transition.transpose() == block_diagonal_matrix(U2, -child)
    return child, transition


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
    simple = matrix(ZZ, [list(root) for root in simple])
    assert simple.nrows() == simple.rank() == data[0]
    return simple, simple * gram * simple.transpose()


def root_adaptation(child):
    _, root_basis, data = roots_and_data(child)
    root_rank = data[0]
    smith, _, smith_right = root_basis.smith_form()
    assert tuple(abs(smith[index, index]) for index in range(root_rank)) == (1,) * root_rank
    simple, cartan = deterministic_simple_roots(child)
    completion = smith_right.inverse()
    adapted_basis = simple.stack(completion[root_rank:])
    assert abs(adapted_basis.det()) == 1
    adapted = adapted_basis * child * adapted_basis.transpose()
    coupling = adapted[:root_rank, root_rank:]
    tail = adapted[root_rank:, root_rank:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling
    scale = lcm(entry.denominator() for entry in height.list())
    lll = matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram())
    quotient = block_diagonal_matrix(identity_matrix(ZZ, root_rank), lll.transpose())
    adapted_basis = quotient * adapted_basis
    adapted = adapted_basis * child * adapted_basis.transpose()
    coupling = adapted[:root_rank, root_rank:]
    tail = adapted[root_rank:, root_rank:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling
    return adapted, adapted_basis, height, data


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
    answer = []
    for component in components(cartan):
        candidates = [
            root for root in roots
            if all(value >= 0 for value in root)
            and all(index in component or root[index] == 0 for index in range(cartan.nrows()))
        ]
        answer.append(max(candidates, key=lambda root: sum(root)))
    return answer


def nef_profile(fibre, frame, root_rank):
    degree = ZZ(fibre[1])
    labels = vector(ZZ, fibre[2:]) * frame[:, :root_rank]
    cartan = frame[:root_rank, :root_rank]
    affine = [degree - top * labels for top in highest_roots(cartan)]
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
            min(tuple(labels) + tuple(affine)) >= 0 and minimum_section >= 0
        ),
    }


def reduce_component_chamber(fibre, frame, root_rank):
    """Return the affine-Weyl-reduced fibre and its right coordinate action."""
    g = block_diagonal_matrix(U2, -frame)
    cartan = frame[:root_rank, :root_rank]
    wall_roots = [
        vector(ZZ, [0, 0] + [-ZZ(index == other) for other in range(17)])
        for index in range(root_rank)
    ]
    wall_names = [f"simple_{index}" for index in range(root_rank)]
    for component_index, top in enumerate(highest_roots(cartan)):
        wall_roots.append(
            vector(ZZ, [1, 0] + list(top) + [0] * (17 - root_rank))
        )
        wall_names.append(f"affine_{component_index}")
    assert all(root * g * root == -2 for root in wall_roots)
    value = vector(ZZ, fibre)
    action = identity_matrix(ZZ, 19)
    reflections = []
    for unused in range(100000):
        pairings = [ZZ(value * g * root) for root in wall_roots]
        negative = next(
            (index for index, pairing in enumerate(pairings) if pairing < 0), None
        )
        if negative is None:
            break
        root = wall_roots[negative]
        pairing = pairings[negative]
        step = identity_matrix(ZZ, 19) + (g * root.column()) * matrix(ZZ, [list(root)])
        assert step * g * step.transpose() == g and abs(step.det()) == 1
        assert value * step == value + pairing * root
        value *= step
        action *= step
        reflections.append({"wall": wall_names[negative], "pairing": int(pairing)})
    else:
        raise RuntimeError("affine-Weyl chamber reduction did not terminate")
    assert action * g * action.transpose() == g
    return value, action, reflections


def refine_component_chamber(anchor, secondary, frame, root_rank):
    """Orient walls invisible to a nef anchor using a second known nef fibre."""
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
    anchor = vector(ZZ, anchor)
    value = vector(ZZ, secondary)
    action = identity_matrix(ZZ, 19)
    reflections = []
    for unused in range(100000):
        pairings = [ZZ(value * g * root) for root in wall_roots]
        negative = next(
            (
                index for index, pairing in enumerate(pairings)
                if pairing < 0 and anchor * g * wall_roots[index] == 0
            ),
            None,
        )
        if negative is None:
            break
        root = wall_roots[negative]
        pairing = pairings[negative]
        step = identity_matrix(ZZ, 19) + (g * root.column()) * matrix(ZZ, [list(root)])
        value *= step
        action *= step
        reflections.append({"wall": wall_names[negative], "pairing": int(pairing)})
    else:
        raise RuntimeError("secondary chamber refinement did not terminate")
    assert anchor * action == anchor
    return value, action, reflections


ranking = json.loads(RANKING.read_text())
marking = json.loads(MARKING.read_text())
crossovers = json.loads(CROSSOVERS.read_text())
semistable_nef = json.loads(SEMISTABLE_NEF.read_text())
assert ranking["status"] == "PASS_EXACT_PINNED_R17_Q4_TARGETED_RANKING"
assert marking["status"] == "PASS_EXACT_REVERSE_HUB_EQUATION_MARKING"
assert crossovers["status"] == "PASS_EXACT_MARKED_TARGET_CROSSOVER_AUDIT"
assert semistable_nef["status"] == "PASS_EXACT_SEMISTABLE_MW2_TO_PINNED_R17_NEF_SUFFIX"
pinned = load_matrix(PINNED_FRAME)
g_pinned = block_diagonal_matrix(U2, -pinned)
targets_pinned = {
    name: vector(ZZ, value)
    for name, value in marking["target_fibres_in_root_adapted_hub"].items()
}


def negative_horizontal_walls(w, degree):
    """Exactly enumerate all negative old-horizontal (-2)-curve walls."""
    walls = []
    w_norm = ZZ(w * pinned * w)
    for old_degree in range(1, int(degree) + 1):
        m = ZZ(old_degree)
        cross = -degree * m * pinned * w.column()
        augmented = block_matrix(ZZ, [
            [degree**2 * pinned, cross],
            [cross.transpose(), matrix(ZZ, [[m**2 * w_norm + 1]])],
        ])
        result = pari(augmented).qfminim(2 * degree**2 - 1)
        half = matrix(ZZ, result[2]).transpose().rows()
        normalized = set()
        for candidate in half:
            if abs(candidate[-1]) != 1:
                continue
            value = candidate if candidate[-1] == 1 else -candidate
            normalized.add(tuple(value))
        for value in normalized:
            x = vector(ZZ, value[:-1])
            x_norm = ZZ(x * pinned * x)
            if (x_norm - 2) % (2 * m):
                continue
            k = ZZ((x_norm - 2) // (2 * m))
            intersection = ZZ((w_norm // (2 * degree)) * m + degree * k - w * pinned * x)
            y = degree * x - m * w
            assert y * pinned * y == 2 * degree * m * intersection + 2 * degree**2
            if intersection < 0:
                walls.append({
                    "old_fibre_degree": int(m),
                    "curve": [int(k), int(m)] + list(map(int, x)),
                    "intersection": int(intersection),
                    "congruence_vector_norm": int(y * pinned * y),
                })
    return sorted(walls, key=lambda item: (item["intersection"], item["curve"]))


if args.fibre:
    fibre = vector(ZZ, [ZZ(value) for value in args.fibre.split(",")])
    assert len(fibre) == 19 and fibre[1] > 0 and fibre * g_pinned * fibre == 0
    q_value = ZZ(fibre[0] * fibre[1])
    marked = {
        name: int(fibre * g_pinned * value) for name, value in targets_pinned.items()
    }
    center = vector(QQ, fibre[2:]) / fibre[1]
    closest = vector(ZZ, next(IntegralLattice(pinned).enumerate_close_vectors(center)))
    distance = (closest - center) * pinned * (closest - center)
    minimum_section = fibre[1] * (distance - 2) / 2
    assert minimum_section >= 0
    horizontal_walls = negative_horizontal_walls(vector(ZZ, fibre[2:]), ZZ(fibre[1]))
    assert not horizontal_walls
    selected = {
        "candidate_id": {
            "q": int(q_value),
            "old_fibre_degree": int(fibre[1]),
            "label": args.candidate_label,
        },
        "fibre_in_pinned_R17": list(map(int, fibre)),
        "minimum_section_intersection": str(minimum_section),
        "minimum_section_distance": str(distance),
        "exact_negative_horizontal_walls": horizontal_walls,
        "exact_horizontal_nef_gate": True,
        "marked_target_degrees": marked,
    }
else:
    selected = ranking["rankings_top_200"][args.target][args.rank]
fibre = vector(ZZ, selected["fibre_in_pinned_R17"])
child, neighbor = child_frame(g_pinned, fibre)
adapted, adaptation, height, root_data = root_adaptation(child)
initial_transition = block_diagonal_matrix(identity_matrix(ZZ, 2), adaptation) * neighbor
initial_inverse = initial_transition.inverse().change_ring(ZZ)
g_child = block_diagonal_matrix(U2, -adapted)
assert abs(initial_transition.det()) == 1
assert initial_transition * g_pinned * initial_transition.transpose() == g_child

anchor_initial = targets_pinned["pinned_R17"] * initial_inverse
anchor_reduced, chamber_action, chamber_reflections = reduce_component_chamber(
    anchor_initial, adapted, root_data[0]
)
orbit12_after_anchor = (
    targets_pinned["orbit12"] * initial_inverse * chamber_action
)
unused_orbit12, refinement_action, refinement_reflections = refine_component_chamber(
    anchor_reduced, orbit12_after_anchor, adapted, root_data[0]
)
chamber_action *= refinement_action
chamber_reflections.extend(
    {**item, "anchor": "orbit12"} for item in refinement_reflections
)
# If x_new=x_old*S, then the new basis is S^{-1} times the old basis.
transition = chamber_action.inverse().change_ring(ZZ) * initial_transition
inverse = transition.inverse().change_ring(ZZ)
assert transition * g_pinned * transition.transpose() == g_child
targets_child = {name: value * inverse for name, value in targets_pinned.items()}
assert targets_child["pinned_R17"] == anchor_reduced
assert {
    name: int(fibre * g_pinned * value) for name, value in targets_pinned.items()
} == selected["marked_target_degrees"]

target_bases = {
    name: matrix(ZZ, value)
    for name, value in crossovers["reverse_target_bases_in_pinned_R17"].items()
}
# The crossover artifact fixes target fibres exactly but its auxiliary q25
# zero/components precede the later chamber-wise nef closeout.  Use the full
# physical q25 basis from that closeout; its first row is the same marked fibre.
target_bases["q25_mw7"] = matrix(
    ZZ, semistable_nef["steps"][3]["inverse_transport"]
)
assert target_bases["q25_mw7"].row(0) == vector(
    ZZ, marking["target_fibres_in_root_adapted_hub"]["q25_mw7"]
)
target_transports = {}
for name, basis in target_bases.items():
    basis_in_child = basis * inverse
    assert abs(basis_in_child.det()) == 1
    target_g = basis_in_child * g_child * basis_in_child.transpose()
    raw_target_frame = -target_g[2:, 2:]
    target_root_data = roots_and_data(raw_target_frame)[2]
    target_root_adaptation = identity_matrix(ZZ, 17)
    if target_root_data[0]:
        target_frame, target_root_adaptation, unused_height, target_root_data = (
            root_adaptation(raw_target_frame)
        )
        basis_in_child = block_diagonal_matrix(
            identity_matrix(ZZ, 2), target_root_adaptation
        ) * basis_in_child
    else:
        target_frame = raw_target_frame
    assert (
        basis_in_child * g_child * basis_in_child.transpose()
        == block_diagonal_matrix(U2, -target_frame)
    )
    target_chamber_reflections = []
    target_anchor_audit = None
    if target_root_data[0]:
        pinned_anchor_in_target = (
            targets_child["pinned_R17"] * basis_in_child.inverse()
        )
        pinned_anchor_reduced, target_action, target_chamber_reflections = (
            reduce_component_chamber(
                pinned_anchor_in_target, target_frame, target_root_data[0]
            )
        )
        basis_in_child = target_action.inverse().change_ring(ZZ) * basis_in_child
        assert (
            targets_child["pinned_R17"] * basis_in_child.inverse()
            == pinned_anchor_reduced
        )
        target_anchor_audit = nef_profile(
            pinned_anchor_reduced, target_frame, target_root_data[0]
        )
    child_fibre_in_target = vector(ZZ, [1, 0] + [0] * 17) * basis_in_child.inverse()
    target_degree = ZZ(child_fibre_in_target[1])
    target_q = ZZ(child_fibre_in_target[0] * child_fibre_in_target[1])
    assert vector(ZZ, child_fibre_in_target[2:]) * target_frame * vector(
        ZZ, child_fibre_in_target[2:]
    ) == 2 * target_q
    target_transports[name] = {
        "target_basis_in_child": rows(basis_in_child),
        "child_basis_in_target": rows(basis_in_child.inverse().change_ring(ZZ)),
        "forward_determinant": int(basis_in_child.det()),
        "target_root_data": list(map(int, target_root_data)),
        "raw_target_to_root_adapted_frame_basis": rows(target_root_adaptation),
        "component_chamber_anchor": "pinned_R17",
        "component_chamber_anchor_nef_audit": target_anchor_audit,
        "component_chamber_anchor_fibre": (
            None
            if target_root_data[0] == 0
            else list(map(int, pinned_anchor_reduced))
        ),
        "component_chamber_affine_weyl_reflections": target_chamber_reflections,
        "reverse_edge_profile": {
            "q": int(target_q),
            "old_fibre_degree": int(target_degree),
            "P_dot_O": int(child_fibre_in_target[0] - target_degree),
            "child_fibre_in_target": list(map(int, child_fibre_in_target)),
            "nef_audit": (
                nef_profile(child_fibre_in_target, target_frame, target_root_data[0])
                if target_root_data[0]
                else {
                    "minimum_section_intersection": selected["minimum_section_intersection"],
                    "rootless_component_gate": "vacuous",
                    "nef_in_selected_component_chamber": True,
                }
            ),
        },
    }

profiles = {}
for name, target in targets_child.items():
    q_value = ZZ(target[0] * target[1])
    assert vector(ZZ, target[2:]) * adapted * vector(ZZ, target[2:]) == 2 * q_value
    profiles[name] = {
        "q": int(q_value),
        "old_fibre_degree": int(target[1]),
        "P_dot_O": int(target[0] - target[1]),
        "fibre_in_child": list(map(int, target)),
        "nef_audit": nef_profile(target, adapted, root_data[0]),
    }

selected_target_nef = profiles[args.target]["nef_audit"][
    "nef_in_selected_component_chamber"
]
selected_reverse_edge_nef = (
    target_transports[args.target]["reverse_edge_profile"]["nef_audit"][
        "nef_in_selected_component_chamber"
    ]
    if args.target in target_transports
    else None
)

FRAME_OUTPUT.write_text(
    "# pinned-R17 targeted neighbor child\n"
    + "\n".join(" ".join(map(str, row)) for row in adapted.rows()) + "\n"
)
payload = {
    "schema": "elkies-k3.h3-pinned-r17-targeted-candidate-certificate.v1",
    "status": "PASS_EXACT_PINNED_R17_TARGETED_CANDIDATE_CERTIFICATE",
    "source_hub": "pinned_R17",
    "candidate_id": selected["candidate_id"],
    "selected_ranking_target": args.target,
    "selected_ranking_index": args.rank,
    "selected_target_nef_in_child_chamber": selected_target_nef,
    "selected_reverse_edge_nef_in_target_chamber": selected_reverse_edge_nef,
    "marked_target_degrees": selected["marked_target_degrees"],
    "child": {
        "root_data": list(map(int, root_data)),
        "mw_rank": 17 - int(root_data[0]),
        "mw_height": rational_rows(height),
    },
    "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "source_to_child_basis": rows(transition),
    "child_to_source_basis": rows(inverse),
    "equation_A11_to_child_basis": rows(
        transition * matrix(ZZ, marking["equation_A11_to_root_adapted_hub_basis"])
    ),
    "target_fibres_in_child": {
        name: list(map(int, value)) for name, value in targets_child.items()
    },
    "target_transports": target_transports,
    "target_profiles": profiles,
    "first_edge_nef_audit": {
        "minimum_section_distance": selected.get("minimum_section_distance"),
        "minimum_section_intersection": selected["minimum_section_intersection"],
        "rootless_component_gate": "vacuous",
        "negative_degree_two_curve_parity_exclusion": bool(fibre[1] == 2),
        "exact_negative_horizontal_walls": selected.get("exact_negative_horizontal_walls", []),
        "exact_horizontal_nef_gate": selected.get("exact_horizontal_nef_gate", fibre[1] == 2),
        "nef": True,
    },
    "child_component_chamber": {
        "anchor": "pinned_R17",
        "anchor_nef_audit": nef_profile(
            targets_child["pinned_R17"], adapted, root_data[0]
        ),
        "affine_weyl_reflections": chamber_reflections,
        "right_coordinate_action": rows(chamber_action),
    },
    "proof_boundary": (
        "Exact first-edge primitive-nef certificate, marked U, full child root data, "
        "and determinant-one transports both ways. Full bases of the known reverse "
        "hubs are transported into the child and checked integrally. A target profile "
        "is a certified second edge only when its component and all-section nef audit passes."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
profile = profiles[args.target]
print(
    "R17TARGETCERT|target={}|candidate={}|child_root={}|MW={}|target_q={}|"
    "target_degree={}|target_PO={}|target_nef={}|det={}|status={}|output={}".format(
        args.target, selected["candidate_id"], ",".join(map(str, root_data)),
        17 - root_data[0], profile["q"], profile["old_fibre_degree"],
        profile["P_dot_O"], int(profile["nef_audit"]["nef_in_selected_component_chamber"]),
        transition.det(), payload["status"], OUTPUT,
    ),
    flush=True,
)
