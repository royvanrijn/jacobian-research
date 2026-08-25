#!/usr/bin/env sage -python
"""Search every compiler-gated q4/q6 A5A5 exit for a cheaper zero loop.

status: ACTIVE_SEARCH
claim: exhaustive equation-cost ranking over the retained q4/q6 first-edge pool
inputs: explicit-zero gate/cost artifacts, orbit12 zero frame, pinned manifest
outputs: artifacts/generated-results/elkies-k3-h3-a5a5-zero-changing-loop-search.json

For each already-explicit degree-one old-A11 component on a gate-passing
child, reframe that child with the component as zero, test the inverse return
to the marked current 2A5 fibre, rebuild the returned 2A5 frame, and test the
marked current 3A3 fibre there.  Retained edges pass component, affine,
all-section, and exact finite horizontal-wall nef gates.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
GATE = GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4q6-explicit-curve-gate.json"
COST = GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4q6-equation-cost.json"
SOURCE_ZERO = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
FINGERPRINT = LOCAL / "q24-a11-q8-construction-fingerprint.json"
DEFAULT_OUTPUT = GENERATED / "elkies-k3-h3-a5a5-zero-changing-loop-search.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--limit", type=int, default=0, help="cost-ranked first-edge limit; zero means all passes")
parser.add_argument("--retain", type=int, default=100)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
OUTPUT = args.output.resolve()


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


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
    return result


def root_adaptation(child):
    _, root_basis, data = roots_and_data(child)
    rank = data[0]
    if rank == 0:
        lll = matrix(ZZ, pari(child.change_ring(ZZ)).qflllgram())
        basis = lll.transpose()
        return basis * child * basis.transpose(), basis, data
    smith, _, right = root_basis.smith_form()
    assert tuple(abs(smith[i, i]) for i in range(rank)) == (1,) * rank
    simple = deterministic_simple_roots(child)
    completion = right.inverse()
    basis = simple.stack(completion[rank:])
    adapted = basis * child * basis.transpose()
    cartan = adapted[:rank, :rank]
    coupling = adapted[:rank, rank:]
    height = adapted[rank:, rank:] - coupling.transpose() * cartan.inverse() * coupling
    scale = lcm(value.denominator() for value in height.list())
    lll = matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram())
    basis = block_diagonal_matrix(identity_matrix(ZZ, rank), lll.transpose()) * basis
    adapted = basis * child * basis.transpose()
    assert abs(basis.det()) == 1
    return adapted, basis, data


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
    return child, transition


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
        for raw in matrix(ZZ, result[2]).transpose().rows():
            if abs(raw[-1]) != 1:
                continue
            value = raw if raw[-1] == 1 else -raw
            x = vector(ZZ, value[:-1])
            x_norm = ZZ(x * frame * x)
            if (x_norm - 2) % (2 * m):
                continue
            k = ZZ((x_norm - 2) // (2 * m))
            intersection = ZZ((w * frame * w // (2 * degree)) * m + degree * k - w * frame * x)
            if intersection < 0:
                walls.append({"old_fibre_degree": int(m), "intersection": int(intersection)})
    return walls


def reduce_component_chamber(fibre, frame, root_rank):
    if root_rank == 0:
        return vector(ZZ, fibre), identity_matrix(ZZ, 19)
    g = block_diagonal_matrix(U2, -frame)
    cartan = frame[:root_rank, :root_rank]
    walls = [
        vector(ZZ, [0, 0] + [-ZZ(index == other) for other in range(17)])
        for index in range(root_rank)
    ]
    walls.extend(
        vector(ZZ, [1, 0] + list(top) + [0] * (17 - root_rank))
        for top in highest_roots(cartan)
    )
    value = vector(ZZ, fibre)
    action = identity_matrix(ZZ, 19)
    for _ in range(100000):
        pairings = [ZZ(value * g * root) for root in walls]
        negative = next((i for i, pairing in enumerate(pairings) if pairing < 0), None)
        if negative is None:
            return value, action
        root = walls[negative]
        step = identity_matrix(ZZ, 19) + (g * root.column()) * matrix(ZZ, [list(root)])
        value *= step
        action *= step
    raise RuntimeError("affine-Weyl reduction did not terminate")


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


def cost_score(fibre, frame, root_rank, child_root_data, nef):
    profile = rr_profile(fibre, frame, root_rank)
    degrees = nef["component_pairings"] + nef["affine_pairings"]
    terms = {
        "P_dot_O": 900 * profile["P_dot_O"],
        "horizontal_degree": 250 * int(fibre[1]),
        "RR_ambient": 120 * profile["expected_RR_ambient"],
        "vertical_layers": 60 * profile["vertical_layers"],
        "vertical_support": 25 * profile["vertical_support"],
        "child_root_count": int(child_root_data[1]),
        "coordinate_growth": max(abs(int(value)) for value in fibre),
        "no_explicit_degree_one_component": 4000 if 1 not in degrees else 0,
        "explicit_degree_one_component_credit": -500 * min(degrees.count(1), 6),
        "explicit_degree_zero_component_credit": -100 * min(degrees.count(0), 12),
    }
    return profile, terms, int(sum(terms.values()))


gate = json.loads(GATE.read_text())
cost = json.loads(COST.read_text())
source_zero = json.loads(SOURCE_ZERO.read_text())
manifest = json.loads(MANIFEST.read_text())
fingerprint = json.loads(FINGERPRINT.read_text())
assert gate["status"] == "PASS_EXACT_A5A5_Q6Q8_EXPLICIT_CURVE_GATE"
assert cost["status"] == "PASS_EXACT_A5A5_EXPLICIT_ZERO_EQUATION_COST_SCORING"
source = source_zero["selected"]
source_frame = matrix(ZZ, source["frame"])
g_source = block_diagonal_matrix(U2, -source_frame)
equation_to_source = matrix(ZZ, source["equation_A11_to_explicit_zero_basis"])
source_to_equation = equation_to_source.inverse().change_ring(ZZ)

historical_in_equation = block_diagonal_matrix(
    identity_matrix(ZZ, 2),
    matrix(ZZ, fingerprint["selected"]["frame_isometry_historical_basis_in_equation_coordinates"]),
)
cumulative = identity_matrix(ZZ, 19)
targets_equation = {"current_A5A5": vector(ZZ, historical_in_equation.row(0))}
for index, step in enumerate(manifest["forward_steps"]):
    if index < 2:
        continue
    cumulative = matrix(ZZ, step["transition"]) * cumulative
    if step["child"] == "2A5/MW7":
        targets_equation["current_A5A5"] = vector(ZZ, (cumulative * historical_in_equation).row(0))
    elif step["child"] == "3A3/MW8":
        targets_equation["current_3A3"] = vector(ZZ, (cumulative * historical_in_equation).row(0))
        break

# Use the pinned ray only to choose deterministic component chambers.
crossovers = json.loads((GENERATED / "elkies-k3-h3-a11-candidate-target-crossovers.json").read_text())
targets_equation["pinned_R17"] = vector(ZZ, next(
    item["target_fibre_in_state"] for item in crossovers["records"]
    if item["state"] == "equation_A11" and item["target"] == "pinned_R17"
))

gate_by_id = {tuple(sorted(item["candidate_id"].items())): item for item in gate["survivors"]}
first_edges = [item for item in cost["ranked_candidates"] if item["full_declared_nef_gate"] == "PASS"]
if args.limit:
    first_edges = first_edges[:args.limit]

records = []
counts = {
    "first_edges": len(first_edges), "explicit_zeros": 0, "return_component_nef": 0,
    "return_exact_nef": 0, "exit_component_nef": 0, "exit_exact_nef": 0,
}
for first_index, scored in enumerate(first_edges):
    raw_gate = gate_by_id[tuple(sorted(scored["candidate_id"].items()))]
    raw = raw_gate["source_neighbor_record"]
    transition_raw = block_diagonal_matrix(
        identity_matrix(ZZ, 2), matrix(ZZ, raw["child_root_adapted_basis"])
    ) * matrix(ZZ, raw["neighbor_basis"])
    raw_frame = matrix(ZZ, raw["child_root_adapted_frame"])
    g_raw = block_diagonal_matrix(U2, -raw_frame)
    assert transition_raw * g_source * transition_raw.transpose() == g_raw
    equation_to_raw = transition_raw * equation_to_source
    raw_to_equation = equation_to_raw.inverse().change_ring(ZZ)
    for zero_name in raw_gate["explicit_degree_one_curves"]:
        if not zero_name.startswith("old_A11_component_"):
            continue
        counts["explicit_zeros"] += 1
        node = int(zero_name.rsplit("_", 1)[1])
        curve_equation = vector(ZZ, [0, 0] + [-ZZ(index == node) for index in range(17)])
        section = curve_equation * raw_to_equation
        assert section * g_raw * section == -2 and section * g_raw * vector(ZZ, [1, 0] + [0] * 17) == 1
        fibre_raw = vector(ZZ, [1, 0] + [0] * 17)
        mate = section + fibre_raw
        complement = matrix(ZZ, [list(fibre_raw * g_raw), list(mate * g_raw)]).right_kernel_matrix()
        split = matrix(ZZ, [list(fibre_raw), list(mate)] + list(complement.rows()))
        raw_positive = -(complement * g_raw * complement.transpose())
        explicit_frame, adaptation, explicit_root_data = root_adaptation(raw_positive)
        reframing = block_diagonal_matrix(identity_matrix(ZZ, 2), adaptation) * split
        equation_to_explicit = reframing * equation_to_raw
        explicit_to_equation = equation_to_explicit.inverse().change_ring(ZZ)
        g_explicit = block_diagonal_matrix(U2, -explicit_frame)
        assert equation_to_explicit * (source_to_equation * g_source * source_to_equation.transpose()) * equation_to_explicit.transpose() == g_explicit

        target_explicit = {
            name: value * explicit_to_equation for name, value in targets_equation.items()
        }
        return_fibre = target_explicit["current_A5A5"]
        return_nef = nef_profile(return_fibre, explicit_frame, int(explicit_root_data[0]))
        if not return_nef["nef_in_selected_component_chamber"]:
            continue
        counts["return_component_nef"] += 1
        return_walls = negative_horizontal_walls(return_fibre, explicit_frame)
        if return_walls:
            continue
        counts["return_exact_nef"] += 1

        returned_raw, returned_transition_raw = child_frame(g_explicit, return_fibre)
        returned_frame, returned_adaptation, returned_root_data = root_adaptation(returned_raw)
        if tuple(returned_root_data) != (10, 60, 36):
            continue
        initial_return_transition = block_diagonal_matrix(
            identity_matrix(ZZ, 2), returned_adaptation
        ) * returned_transition_raw
        initial_return_inverse = initial_return_transition.inverse().change_ring(ZZ)
        anchor = target_explicit["pinned_R17"] * initial_return_inverse
        _, chamber_action = reduce_component_chamber(anchor, returned_frame, 10)
        return_transition = chamber_action.inverse().change_ring(ZZ) * initial_return_transition
        return_inverse = return_transition.inverse().change_ring(ZZ)
        g_returned = block_diagonal_matrix(U2, -returned_frame)
        assert return_transition * g_explicit * return_transition.transpose() == g_returned
        exit_fibre = target_explicit["current_3A3"] * return_inverse
        assert exit_fibre * g_returned * exit_fibre == 0
        if exit_fibre[1] <= 0 or gcd(tuple(g_returned * exit_fibre)) != 1:
            continue
        exit_nef = nef_profile(exit_fibre, returned_frame, 10)
        if not exit_nef["nef_in_selected_component_chamber"]:
            continue
        counts["exit_component_nef"] += 1
        exit_walls = negative_horizontal_walls(exit_fibre, returned_frame)
        if exit_walls:
            continue
        counts["exit_exact_nef"] += 1

        exit_child_raw, _ = child_frame(g_returned, exit_fibre)
        _, _, exit_child_root_data = root_adaptation(exit_child_raw)
        if tuple(exit_child_root_data) != (9, 36, 64):
            continue
        return_profile, return_terms, return_score = cost_score(
            return_fibre, explicit_frame, int(explicit_root_data[0]), returned_root_data, return_nef
        )
        exit_profile, exit_terms, exit_score = cost_score(
            exit_fibre, returned_frame, 10, exit_child_root_data, exit_nef
        )
        records.append({
            "first_edge_candidate_id": scored["candidate_id"],
            "first_edge_equation_cost_score": scored["equation_cost_score"],
            "first_edge_child": scored["child"],
            "explicit_zero_curve": zero_name,
            "explicit_child_root_data": list(map(int, explicit_root_data)),
            "equation_A11_to_explicit_child_basis": rows(equation_to_explicit),
            "return_fibre_in_explicit_child": entries(return_fibre),
            "return_nef_audit": return_nef,
            "return_exact_negative_horizontal_walls": return_walls,
            "return_profile": return_profile,
            "return_equation_cost_terms": return_terms,
            "return_equation_cost_score": return_score,
            "return_transition": rows(return_transition),
            "returned_A5A5_frame": rows(returned_frame),
            "exit_3A3_fibre_in_returned_A5A5": entries(exit_fibre),
            "exit_q": int(exit_fibre[0] * exit_fibre[1]),
            "exit_old_fibre_degree": int(exit_fibre[1]),
            "exit_nef_audit": exit_nef,
            "exit_exact_negative_horizontal_walls": exit_walls,
            "exit_profile": exit_profile,
            "exit_equation_cost_terms": exit_terms,
            "exit_equation_cost_score": exit_score,
            "loop_equation_cost_score": int(scored["equation_cost_score"] + return_score + exit_score),
            "loop_q_sequence": [scored["candidate_id"]["q"], int(return_fibre[0] * return_fibre[1]), int(exit_fibre[0] * exit_fibre[1])],
            "loop_old_fibre_degrees": [2, int(return_fibre[1]), int(exit_fibre[1])],
            "coordinate_growth_max": max(
                max(abs(int(value)) for value in return_fibre),
                max(abs(int(value)) for value in exit_fibre),
            ),
        })
    if (first_index + 1) % 25 == 0:
        print("ZEROLOOPPROGRESS|first_edges={}|zeros={}|exact_returns={}|exact_exits={}".format(
            first_index + 1, counts["explicit_zeros"], counts["return_exact_nef"], counts["exit_exact_nef"]
        ), flush=True)

records.sort(key=lambda item: (
    item["loop_equation_cost_score"], item["coordinate_growth_max"],
    item["first_edge_candidate_id"]["q"], item["first_edge_candidate_id"]["orbit_index"],
    item["explicit_zero_curve"],
))
inputs = (GATE, COST, SOURCE_ZERO, MANIFEST, FINGERPRINT)
payload = {
    "schema": "elkies-k3.h3-a5a5-zero-changing-loop-search.v1",
    "status": "PASS_EXACT_ZERO_CHANGING_LOOP_SEARCH",
    "search_parameters": {"first_edge_limit": args.limit, "retained": args.retain},
    "counts": counts,
    "direct_q104_comparator_score": 13518,
    "strict_cost_winner_count": sum(item["loop_equation_cost_score"] < 13518 for item in records),
    "ranked_loops": records[:args.retain],
    "proof_boundary": (
        "Every retained loop uses an already-explicit old-A11 component as child zero. "
        "Both return and 3A3 exit pass exact component, affine, all-section, and finite "
        "horizontal-wall nef gates; full composed route certification is separate."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
best = records[0] if records else None
print("ZEROLOOP|records={}|strict={}|best_score={}|best={}|status={}|output={}".format(
    len(records), payload["strict_cost_winner_count"], None if best is None else best["loop_equation_cost_score"],
    None if best is None else "q{}o{}:{}".format(best["first_edge_candidate_id"]["q"], best["first_edge_candidate_id"]["orbit_index"], best["explicit_zero_curve"]),
    payload["status"], OUTPUT,
), flush=True)
