#!/usr/bin/env sage -python
"""Choose the cheapest equation-explicit zero on a certified marked child.

The candidate zeros are source simple/affine fibre components having degree
one over the new fibre.  For each zero, the child fibre's physical components
are re-based by deleting the component met by that zero in every reducible
fibre, preserving the declared effective chamber instead of choosing an
arbitrary abstract positive-root system.
"""

import argparse
import hashlib
import json
from collections import deque
from pathlib import Path

from sage.all import *
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]
U2 = matrix(ZZ, ((0, 1), (1, 0)))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--certificate", type=Path, required=True)
parser.add_argument("--source-frame", type=Path, required=True)
parser.add_argument("--next-target", required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--frame-output", type=Path, required=True)
args = parser.parse_args()
CERTIFICATE = args.certificate.resolve()
SOURCE_FRAME = args.source_frame.resolve()
OUTPUT = args.output.resolve()
FRAME_OUTPUT = args.frame_output.resolve()


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


def components(cartan):
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


def highest_roots(cartan):
    if cartan.nrows() == 0:
        return ()
    half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
    roots = tuple(half) + tuple(-root for root in half)
    answer = []
    for component in components(cartan):
        candidates = [
            root for root in roots
            if all(value >= 0 for value in root)
            and all(index in component or root[index] == 0 for index in range(cartan.nrows()))
        ]
        answer.append(max(candidates, key=lambda root: sum(root)))
    return tuple(answer)


def component_curves(frame, root_rank):
    cartan = frame[:root_rank, :root_rank]
    answer = [
        ("simple_{}".format(index), vector(ZZ, [0, 0] + [-ZZ(index == other) for other in range(17)]))
        for index in range(root_rank)
    ]
    answer.extend(
        ("affine_{}".format(index), vector(ZZ, [1, 0] + list(top) + [0] * (17 - root_rank)))
        for index, top in enumerate(highest_roots(cartan))
    )
    return answer


def reduce_component_chamber(fibre, frame, root_rank):
    if root_rank == 0:
        return vector(ZZ, fibre), identity_matrix(ZZ, 19), []
    g = block_diagonal_matrix(U2, -frame)
    walls = [curve for _, curve in component_curves(frame, root_rank)]
    names = [name for name, _ in component_curves(frame, root_rank)]
    value = vector(ZZ, fibre)
    action = identity_matrix(ZZ, 19)
    reflections = []
    for _ in range(100000):
        pairings = [ZZ(value * g * root) for root in walls]
        negative = next((index for index, pairing in enumerate(pairings) if pairing < 0), None)
        if negative is None:
            return value, action, reflections
        root = walls[negative]
        pairing = pairings[negative]
        step = identity_matrix(ZZ, 19) + (g * root.column()) * matrix(ZZ, [list(root)])
        value *= step
        action *= step
        reflections.append({"wall": names[negative], "pairing": int(pairing)})
    raise RuntimeError("affine-Weyl chamber reduction did not terminate")


def common_component_chamber(section, target, frame, root_rank):
    """Find a declared component chamber containing both effective classes."""
    g = block_diagonal_matrix(U2, -frame)
    walls = [curve for _, curve in component_curves(frame, root_rank)]
    for weight in tuple(ZZ(2) ** exponent for exponent in range(0, 31)):
        _, action, reflections = reduce_component_chamber(
            weight * vector(ZZ, target) + vector(ZZ, section), frame, root_rank
        )
        section_value = vector(ZZ, section) * action
        target_value = vector(ZZ, target) * action
        section_pairings = [ZZ(section_value * g * wall) for wall in walls]
        target_pairings = [ZZ(target_value * g * wall) for wall in walls]
        if (
            min(section_pairings + [ZZ(0)]) >= 0
            and min(target_pairings + [ZZ(0)]) >= 0
            and all(value in (0, 1) for value in section_pairings)
        ):
            return section_value, target_value, action, reflections, int(weight)
    # A weighted anchor can miss a chamber lying across a target-positive wall.
    # Search the simultaneous affine-Weyl orbit, branching only at a wall that
    # is currently negative for one of the two required effective classes.
    start_section = vector(ZZ, section)
    start_target = vector(ZZ, target)
    queue = deque([(start_section, start_target, identity_matrix(ZZ, 19), [])])
    visited = {(tuple(start_section), tuple(start_target))}
    for _ in range(200000):
        if not queue:
            break
        section_value, target_value, action, reflections = queue.popleft()
        section_pairings = [ZZ(section_value * g * wall) for wall in walls]
        target_pairings = [ZZ(target_value * g * wall) for wall in walls]
        if min(section_pairings + [ZZ(0)]) >= 0 and min(target_pairings + [ZZ(0)]) >= 0:
            assert all(value in (0, 1) for value in section_pairings)
            return section_value, target_value, action, reflections, 0
        negative_walls = sorted(set(
            [index for index, value in enumerate(section_pairings) if value < 0]
            + [index for index, value in enumerate(target_pairings) if value < 0]
        ))
        for index in negative_walls:
            root = walls[index]
            step = identity_matrix(ZZ, 19) + (g * root.column()) * matrix(ZZ, [list(root)])
            next_section = section_value * step
            next_target = target_value * step
            key = (tuple(next_section), tuple(next_target))
            if key in visited:
                continue
            visited.add(key)
            queue.append((
                next_section, next_target, action * step,
                reflections + [{"wall": component_curves(frame, root_rank)[index][0],
                                "section_pairing": int(section_pairings[index]),
                                "target_pairing": int(target_pairings[index])}],
            ))
    return None


def nef_profile(fibre, frame, root_rank):
    degree = ZZ(fibre[1])
    labels = vector(ZZ, fibre[2:]) * frame[:, :root_rank]
    affine = [degree - top * labels for top in highest_roots(frame[:root_rank, :root_rank])]
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
            min(tuple(labels) + tuple(affine) + (ZZ(0),)) >= 0 and minimum_section >= 0
        ),
    }


def negative_horizontal_walls(fibre, frame):
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
                (w * frame * w // (2 * degree)) * m + degree * k - w * frame * x
            )
            if intersection < 0:
                walls.append({"old_fibre_degree": int(m), "intersection": int(intersection)})
    return walls


def vertical_layers(coefficients, cartan):
    edges = [(i, j) for i in range(cartan.nrows()) for j in range(i + 1, cartan.nrows()) if cartan[i, j] == -1]
    magnitudes = [abs(ZZ(value)) for value in coefficients]
    previous = total = 0
    for level in sorted(set(value for value in magnitudes if value)):
        active = set(i for i, value in enumerate(magnitudes) if value >= level)
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
    base = vector(ZZ, [0] * root_rank + list(fibre[2 + root_rank:]))
    dual = vector(QQ, base * frame[:, :root_rank]) * root.inverse()
    iterator = IntegralLattice(root).enumerate_close_vectors(-dual)
    minimum = None
    choices = []
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
        "vertical_layers": layers,
        "vertical_support": sum(value != 0 for value in vertical),
        "expected_RR_ambient": 2 + 2 * int(pole) + layers,
    }


certificate = json.loads(CERTIFICATE.read_text())
assert certificate["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
source_frame = load_matrix(SOURCE_FRAME)
source_root_rank = len(certificate["first_edge_nef_audit"]["component_pairings"])
source_g = block_diagonal_matrix(U2, -source_frame)
transition = matrix(ZZ, certificate["source_to_child_basis"])
inverse = matrix(ZZ, certificate["child_to_source_basis"])
child_frame_path = ROOT / certificate["frame_output"]
child_frame = load_matrix(child_frame_path)
child_g = block_diagonal_matrix(U2, -child_frame)
child_root_rank = int(certificate["child"]["root_data"][0])
assert transition * source_g * transition.transpose() == child_g
targets_child = {name: vector(ZZ, value) for name, value in certificate["target_fibres_in_child"].items()}
assert args.next_target in targets_child and "pinned_R17" in targets_child

source_components = component_curves(source_frame, source_root_rank)
degree_one_source = [
    (name, curve) for name, curve in source_components
    if curve * source_g * vector(ZZ, transition.row(0)) == 1
]
child_components = component_curves(child_frame, child_root_rank)
fibre_child = vector(ZZ, [1, 0] + [0] * 17)
candidates = []
for zero_name, zero_source in degree_one_source:
    zero_child = zero_source * inverse
    assert zero_child * child_g * zero_child == -2 and zero_child * child_g * fibre_child == 1
    common = common_component_chamber(
        zero_child, targets_child[args.next_target], child_frame, child_root_rank
    )
    if common is None:
        continue
    zero_child, _, chamber_action, chamber_reflections, chamber_weight = common
    chamber_basis = chamber_action.inverse().change_ring(ZZ)
    aligned_targets_child = {
        name: value * chamber_action for name, value in targets_child.items()
    }
    equation_to_aligned_child = chamber_basis * matrix(
        ZZ, certificate["equation_A11_to_child_basis"]
    )
    mate = zero_child + fibre_child
    complement = matrix(ZZ, [list(fibre_child * child_g), list(mate * child_g)]).right_kernel_matrix()
    split = matrix(ZZ, [list(fibre_child), list(mate)] + list(complement.rows()))
    split_inverse = split.inverse().change_ring(ZZ)
    raw_frame = -(complement * child_g * complement.transpose())

    # Keep precisely the effective components not met by the selected zero.
    physical_roots = []
    physical_names = []
    for component_name, component in child_components:
        pairing = ZZ(zero_child * child_g * component)
        assert pairing in (0, 1)
        if pairing == 0:
            component_new = component * split_inverse
            assert component_new[:2] == vector(ZZ, [0, 0])
            physical_roots.append(vector(ZZ, component_new[2:]))
            physical_names.append(component_name)
    assert len(physical_roots) == child_root_rank
    physical_root_basis = matrix(ZZ, [list(root) for root in physical_roots])
    assert physical_root_basis.rank() == child_root_rank
    smith, _, right = physical_root_basis.smith_form()
    assert tuple(abs(smith[index, index]) for index in range(child_root_rank)) == (1,) * child_root_rank
    completion = right.inverse()
    adaptation = physical_root_basis.stack(completion[child_root_rank:])
    adapted = adaptation * raw_frame * adaptation.transpose()
    cartan = adapted[:child_root_rank, :child_root_rank]
    coupling = adapted[:child_root_rank, child_root_rank:]
    height = adapted[child_root_rank:, child_root_rank:] - coupling.transpose() * cartan.inverse() * coupling
    scale = lcm(value.denominator() for value in height.list())
    lll = matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram())
    quotient = block_diagonal_matrix(identity_matrix(ZZ, child_root_rank), lll.transpose())
    adaptation = quotient * adaptation
    adapted = adaptation * raw_frame * adaptation.transpose()
    reframing = block_diagonal_matrix(identity_matrix(ZZ, 2), adaptation) * split
    reframing_inverse = reframing.inverse().change_ring(ZZ)
    adapted_g = block_diagonal_matrix(U2, -adapted)
    assert abs(reframing.det()) == 1
    assert reframing * child_g * reframing.transpose() == adapted_g

    targets = {
        name: value * reframing_inverse
        for name, value in aligned_targets_child.items()
    }
    next_fibre = targets[args.next_target]
    if next_fibre[1] <= 0:
        continue
    nef = nef_profile(next_fibre, adapted, child_root_rank)
    walls = negative_horizontal_walls(next_fibre, adapted)
    exact_nef = nef["nef_in_selected_component_chamber"] and not walls
    profile = rr_profile(next_fibre, adapted, child_root_rank)
    degrees = nef["component_pairings"] + nef["affine_pairings"]
    terms = {
        "P_dot_O": 900 * profile["P_dot_O"],
        "horizontal_degree": 250 * int(next_fibre[1]),
        "RR_ambient": 120 * profile["expected_RR_ambient"],
        "vertical_layers": 60 * profile["vertical_layers"],
        "vertical_support": 25 * profile["vertical_support"],
        "coordinate_growth": max(abs(int(value)) for value in next_fibre),
        "no_explicit_degree_one_component": 4000 if 1 not in degrees else 0,
        "explicit_degree_one_component_credit": -500 * min(degrees.count(1), 6),
        "explicit_degree_zero_component_credit": -100 * min(degrees.count(0), 12),
    }
    candidates.append({
        "zero_source_component": zero_name,
        "zero_in_certified_child": entries(zero_child),
        "child_component_chamber_reflections": chamber_reflections,
        "child_component_chamber_anchor_target_weight": chamber_weight,
        "child_component_chamber_basis": rows(chamber_basis),
        "physical_nonidentity_components": physical_names,
        "frame": rows(adapted),
        "reframing_child_to_explicit_zero_basis": rows(reframing),
        "explicit_zero_to_certified_child_basis": rows(reframing_inverse),
        "equation_A11_to_explicit_zero_basis": rows(reframing * equation_to_aligned_child),
        "target_fibres_in_root_adapted_hub": {name: entries(value) for name, value in targets.items()},
        "next_target": args.next_target,
        "next_fibre": entries(next_fibre),
        "next_q": int(next_fibre[0] * next_fibre[1]),
        "next_old_fibre_degree": int(next_fibre[1]),
        "next_nef_audit": nef,
        "next_exact_negative_horizontal_walls": walls,
        "next_exact_nef_gate": exact_nef,
        "next_horizontal_RR_profile": profile,
        "next_equation_cost_terms": terms,
        "next_equation_cost_score_without_child_root_count": int(sum(terms.values())),
        "frame_max_abs": max(abs(int(value)) for value in adapted.list()),
        "frame_L1": sum(abs(int(value)) for value in adapted.list()),
    })

candidates.sort(key=lambda item: (
    not item["next_exact_nef_gate"], item["next_equation_cost_score_without_child_root_count"],
    item["next_q"], item["frame_max_abs"], item["frame_L1"], item["zero_source_component"],
))
assert candidates
selected = candidates[0]
selected_frame = matrix(ZZ, selected["frame"])
FRAME_OUTPUT.write_text(
    "# marked child explicit zero {}\n".format(selected["zero_source_component"])
    + "\n".join(" ".join(map(str, row)) for row in selected_frame.rows()) + "\n"
)
inputs = (CERTIFICATE, SOURCE_FRAME, child_frame_path)
payload = {
    "schema": "elkies-k3.h3-marked-child-explicit-component-zero.v1",
    "status": "PASS_EXACT_MARKED_CHILD_EXPLICIT_ZERO_SELECTION",
    "source_certificate": str(CERTIFICATE.relative_to(ROOT)),
    "source_hub": certificate.get("source_hub"),
    "child": certificate["child"],
    "root_data": certificate["child"]["root_data"],
    "next_target": args.next_target,
    "candidate_count": len(candidates),
    "exact_nef_candidate_count": sum(item["next_exact_nef_gate"] for item in candidates),
    "selection_rule": "exact-nef first, then equation cost, q, frame growth, and component name",
    "selected": selected,
    "candidates": candidates,
    "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "equation_A11_to_root_adapted_hub_basis": selected["equation_A11_to_explicit_zero_basis"],
    "root_adapted_hub_to_equation_A11_basis": rows(matrix(ZZ, selected["equation_A11_to_explicit_zero_basis"]).inverse().change_ring(ZZ)),
    "target_fibres_in_root_adapted_hub": selected["target_fibres_in_root_adapted_hub"],
    "proof_boundary": (
        "Exact reframing by an already-explicit source fibre component. The effective "
        "child component chamber is preserved by deleting precisely the components met "
        "by the new zero. The selected next target has exact component, affine, all-section, "
        "and finite horizontal-wall nef audits. Its child transport is certified separately."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("CHILDZERO|zeros={}|exact_nef={}|selected={}|next={}|q={}|degree={}|cost={}|status={}|output={}".format(
    len(candidates), payload["exact_nef_candidate_count"], selected["zero_source_component"],
    args.next_target, selected["next_q"], selected["next_old_fibre_degree"],
    selected["next_equation_cost_score_without_child_root_count"], payload["status"], OUTPUT,
), flush=True)
