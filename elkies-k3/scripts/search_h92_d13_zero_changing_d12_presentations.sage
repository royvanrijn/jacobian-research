#!/usr/bin/env sage -python
"""Search cheap zero loops that re-present a fixed next-route fibre.

status: ACTIVE_SEARCH
claim: exhaustive zero-changing loop search over one supplied compact frontier
inputs: exact equation marking/frame and a marked neighbour frontier
outputs: a mode-specific JSON artifact under artifacts/generated-results/

Each first child is re-zeroed by every source-fibre component of degree one.
The original source fibre is then taken back, and the exact marked exit fibre
is priced in the returned source marking.  Every retained edge passes
component, affine, all-section, and Proposition-C2 finite horizontal-wall
gates.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
D13_MARKING = GENERATED / "elkies-k3-h3-equation-d13-marking.json"
D13_FRAME = GENERATED / "elkies-k3-h3-equation-d13-root-adapted-frame.txt"
D13_FRONTIER = GENERATED / "elkies-k3-h3-equation-d13-q4q6q8-marked-frontier-adapted-compact.json"
D13_OUTPUT = GENERATED / "elkies-k3-h3-d13-zero-changing-d12-presentations.json"
D12_MARKING = GENERATED / "elkies-k3-h3-current_D12-marked-frame.json"
D12_FRONTIER = GENERATED / "elkies-k3-h3-current-d12-q4q6q8-current-a11-marked-frontier.json"
D12_OUTPUT = GENERATED / "elkies-k3-h3-d12-zero-changing-a11-presentations.json"
A11_MARKING = GENERATED / "elkies-k3-h3-equation-a11-marking.json"
A11_FRONTIER = GENERATED / "elkies-k3-h3-a11-equation-cost-neighbors-all.json"
A11_OUTPUT = GENERATED / "elkies-k3-h3-a11-zero-changing-q8-presentations.json"
A5_MARKING = GENERATED / "elkies-k3-h3-a5a5-q6o1307-q4-return-a5a5-certificate.json"
A5_FRONTIER = GENERATED / "elkies-k3-h3-a5a5-q6o1307-returned-q4q6q8-current-3a3-frontier.json"
A5_OUTPUT = GENERATED / "elkies-k3-h3-a5a5-q6o1307-second-zero-changing-3a3-presentations.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--mode", choices=("e8e6", "d13", "d12", "a11", "a5"), default="d13")
parser.add_argument("--marking", type=Path,
                    help="override the mode's exact source marking")
parser.add_argument("--frame", type=Path,
                    help="override the mode's exact source frame")
parser.add_argument("--frontier", type=Path)
parser.add_argument("--limit", type=int, default=0)
parser.add_argument("--retain", type=int, default=100)
parser.add_argument("--output", type=Path)
parser.add_argument(
    "--max-operational-total", type=int,
    help=("prune a presentation once its inherited-explicit first edge plus the "
          "two unavoidable degree-two horizontal floors cannot beat this total"),
)
args = parser.parse_args()
if args.max_operational_total is not None and args.max_operational_total <= 1000:
    parser.error("--max-operational-total must exceed the two remaining 500-point floors")
if args.mode == "e8e6":
    if not args.marking:
        parser.error("--mode e8e6 requires --marking")
    MARKING = args.marking.resolve()
    FRAME = ROOT / json.loads(MARKING.read_text())["frame_output"]
    if not args.frontier or not args.output:
        parser.error("--mode e8e6 requires --frontier and --output")
    FRONTIER = args.frontier.resolve()
    OUTPUT = args.output.resolve()
    TARGET_KEY, SOURCE_LABEL, EXIT_LABEL, DIRECT_LABEL = "equation_D13", "E8E6", "D13", "q8"
elif args.mode == "d13":
    MARKING, FRAME = D13_MARKING, D13_FRAME
    FRONTIER = (args.frontier or D13_FRONTIER).resolve()
    OUTPUT = (args.output or D13_OUTPUT).resolve()
    TARGET_KEY, SOURCE_LABEL, EXIT_LABEL, DIRECT_LABEL = "current_0_D12", "D13", "D12", "q24"
elif args.mode == "d12":
    MARKING = D12_MARKING
    FRAME = ROOT / json.loads(D12_MARKING.read_text())["frame_output"]
    FRONTIER = (args.frontier or D12_FRONTIER).resolve()
    OUTPUT = (args.output or D12_OUTPUT).resolve()
    TARGET_KEY, SOURCE_LABEL, EXIT_LABEL, DIRECT_LABEL = "current_A11", "D12", "A11", "q6_orbit42"
elif args.mode == "a11":
    MARKING = A11_MARKING
    FRAME = ROOT / json.loads(A11_MARKING.read_text())["frame_output"]
    FRONTIER = (args.frontier or A11_FRONTIER).resolve()
    OUTPUT = (args.output or A11_OUTPUT).resolve()
    TARGET_KEY, SOURCE_LABEL, EXIT_LABEL, DIRECT_LABEL = "orbit12", "A11", "2A5", "q8_orbit12"
else:
    MARKING = A5_MARKING
    FRAME = ROOT / json.loads(A5_MARKING.read_text())["frame_output"]
    FRONTIER = (args.frontier or A5_FRONTIER).resolve()
    OUTPUT = (args.output or A5_OUTPUT).resolve()
    TARGET_KEY, SOURCE_LABEL, EXIT_LABEL, DIRECT_LABEL = "current_3A3", "A5A5", "3A3", "q6_current"

if args.marking:
    MARKING = args.marking.resolve()
if args.frame:
    FRAME = args.frame.resolve()
elif args.marking:
    FRAME = ROOT / json.loads(MARKING.read_text())["frame_output"]


def load_matrix(path):
    return matrix(ZZ, [[ZZ(x) for x in line.split()] for line in path.read_text().splitlines()
                       if line.strip() and not line.lstrip().startswith("#")])


def entries(value):
    return [int(x) for x in vector(ZZ, value)]


def rows(value):
    return [[int(x) for x in row] for row in value.rows()]


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    if count == 0:
        return (), matrix(ZZ, 0, gram.nrows()), (0, 0, 1)
    half = [vector(ZZ, column) for column in matrix(ZZ, result[2]).columns()]
    roots = tuple(half + [-root for root in half])
    basis = matrix(ZZ, [list(root) for root in roots]).row_module().basis_matrix()
    return roots, basis, (basis.rank(), count, abs(ZZ((basis * gram * basis.transpose()).det())))


def deterministic_simple_roots(gram):
    roots, _, data = roots_and_data(gram)
    positive = [root for root in roots if next(x for x in root if x) > 0]
    positive_set = {tuple(root) for root in positive}
    simple = matrix(ZZ, [list(root) for root in positive
                         if not any(tuple(root - left) in positive_set for left in positive)])
    assert simple.nrows() == simple.rank() == data[0]
    return simple


def root_adaptation(child):
    _, root_basis, data = roots_and_data(child)
    rank = data[0]
    if rank == 0:
        basis = matrix(ZZ, pari(child.change_ring(ZZ)).qflllgram()).transpose()
        return basis * child * basis.transpose(), basis, data, matrix(ZZ, 0, child.nrows())
    smith, _, right = root_basis.smith_form()
    smith_diagonal = tuple(abs(smith[i, i]) for i in range(rank))
    simple = deterministic_simple_roots(child)
    if smith_diagonal == (1,) * rank:
        # Preserve the historical primitive-root frame exactly.
        basis = simple.stack(right.inverse()[rank:])
    else:
        # In a torsion fibration the root lattice is nonprimitive.  Start the
        # full frame with its saturation and retain the actual simple roots as
        # embedded vectors in that unimodular frame.
        saturated = matrix(ZZ, root_basis.row_module(ZZ).saturation().basis_matrix())
        sat_smith, _, sat_right = saturated.smith_form()
        assert tuple(abs(sat_smith[i, i]) for i in range(rank)) == (1,) * rank
        basis = saturated.stack(sat_right.inverse()[rank:])
    adapted = basis * child * basis.transpose()
    simple_adapted = simple * basis.inverse().change_ring(ZZ)
    cartan = simple_adapted * adapted * simple_adapted.transpose()
    coupling = adapted[:rank, rank:]
    saturated_block = adapted[:rank, :rank]
    height = adapted[rank:, rank:] - coupling.transpose() * saturated_block.inverse() * coupling
    scale = lcm(x.denominator() for x in height.list())
    lll = matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram())
    basis = block_diagonal_matrix(identity_matrix(ZZ, rank), lll.transpose()) * basis
    adapted = basis * child * basis.transpose()
    simple_adapted = simple * basis.inverse().change_ring(ZZ)
    assert abs(basis.det()) == 1
    assert simple_adapted * adapted * simple_adapted.transpose() == cartan
    return adapted, basis, data, simple_adapted


def bezout_vector_for_pairing(ns, fibre):
    current = ZZ(0)
    answer = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(ns * fibre):
        if value == 0:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        answer = [left * x for x in answer]
        answer[index] += right
        current = divisor
    assert abs(current) == 1
    return vector(ZZ, answer if current == 1 else [-x for x in answer])


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
        candidates = [root for root in roots if all(x >= 0 for x in root)
                      and all(i in component or root[i] == 0 for i in range(cartan.nrows()))]
        answer.append(max(candidates, key=lambda root: sum(root)))
    return tuple(answer)


def nef_profile(fibre, frame, simple_roots):
    degree = ZZ(fibre[1])
    cartan = simple_roots * frame * simple_roots.transpose()
    labels = vector(ZZ, fibre[2:]) * frame * simple_roots.transpose()
    affine = [degree - top * labels for top in highest_roots(cartan)]
    center = vector(QQ, fibre[2:]) / degree
    closest = vector(ZZ, next(IntegralLattice(frame).enumerate_close_vectors(center)))
    distance = (closest - center) * frame * (closest - center)
    minimum = degree * (distance - 2) / 2
    return {
        "component_pairings": list(map(int, labels)),
        "affine_pairings": list(map(int, affine)),
        "minimum_section_intersection": str(minimum),
        "nef": bool(min(tuple(labels) + tuple(affine) + (ZZ(0),)) >= 0 and minimum >= 0),
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
            if abs(raw[-1]) == 1:
                value = raw if raw[-1] == 1 else -raw
                normalized.add(tuple(value))
        for value in normalized:
            x = vector(ZZ, value[:-1])
            norm = ZZ(x * frame * x)
            if (norm - 2) % (2 * m):
                continue
            k = ZZ((norm - 2) // (2 * m))
            intersection = ZZ((w * frame * w // (2 * degree)) * m + degree * k - w * frame * x)
            if intersection < 0:
                walls.append({"old_fibre_degree": int(m), "intersection": int(intersection)})
    return walls


def reduce_component_chamber(fibre, frame, simple_roots):
    g = block_diagonal_matrix(U2, -frame)
    cartan = simple_roots * frame * simple_roots.transpose()
    walls = [vector(ZZ, [0, 0] + list(-root)) for root in simple_roots.rows()]
    walls += [vector(ZZ, [1, 0] + list(top * simple_roots))
              for top in highest_roots(cartan)]
    value = vector(ZZ, fibre)
    action = identity_matrix(ZZ, 19)
    for _ in range(100000):
        pairings = [ZZ(value * g * root) for root in walls]
        negative = next((i for i, x in enumerate(pairings) if x < 0), None)
        if negative is None:
            return value, action
        root = walls[negative]
        step = identity_matrix(ZZ, 19) + (g * root.column()) * matrix(ZZ, [list(root)])
        value *= step
        action *= step
    raise RuntimeError("affine-Weyl reduction did not terminate")


def vertical_layers(coefficients, cartan):
    edges = [(i, j) for i in range(cartan.nrows()) for j in range(i + 1, cartan.nrows())
             if cartan[i, j] == -1]
    magnitudes = [abs(ZZ(x)) for x in coefficients]
    previous = total = 0
    for level in sorted(set(x for x in magnitudes if x)):
        active = set(i for i, x in enumerate(magnitudes) if x >= level)
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


def rr_profile(fibre, frame, simple_roots):
    root = simple_roots * frame * simple_roots.transpose()
    coefficients = vector(ZZ, fibre[2:])
    dual = vector(QQ, coefficients * frame * simple_roots.transpose()) * root.inverse()
    iterator = IntegralLattice(root).enumerate_close_vectors(-dual)
    choices = []
    admissible_norm = None
    for shift in iterator:
        lifted = coefficients + vector(ZZ, shift) * simple_roots
        norm = QQ(lifted * frame * lifted)
        if admissible_norm is not None and norm > admissible_norm:
            break
        pole = (norm - 4) / 2
        if pole in ZZ and pole >= 0:
            if admissible_norm is None:
                admissible_norm = norm
            choices.append((ZZ(pole), lifted))
    assert choices
    pole, lifted = min(choices, key=lambda item: (item[0], tuple(item[1])))
    section = vector(ZZ, [pole + 1, 1] + list(lifted))
    vertical_class = vector(ZZ, (fibre - vector(ZZ, [-1, 1] + [0] * 17) - section)[2:])
    vertical = vector(
        QQ,
        vertical_class * simple_roots.transpose()
        * (simple_roots * simple_roots.transpose()).inverse(),
    )
    assert vertical * simple_roots == vertical_class and all(value in ZZ for value in vertical)
    vertical = vector(ZZ, vertical)
    layers = vertical_layers(vertical, root)
    return {"P_dot_O": int(pole), "vertical_layers": layers,
            "vertical_support": sum(x != 0 for x in vertical),
            "expected_RR_ambient": 2 + 2 * int(pole) + layers}


def cost_score(fibre, frame, simple_roots, child_root_data, nef):
    profile = rr_profile(fibre, frame, simple_roots)
    degrees = nef["component_pairings"] + nef["affine_pairings"]
    terms = {
        "P_dot_O": 900 * profile["P_dot_O"],
        "horizontal_degree": 250 * int(fibre[1]),
        "RR_ambient": 120 * profile["expected_RR_ambient"],
        "vertical_layers": 60 * profile["vertical_layers"],
        "vertical_support": 25 * profile["vertical_support"],
        "child_root_count": int(child_root_data[1]),
        "coordinate_growth": max(abs(int(x)) for x in fibre),
        "no_explicit_degree_one_component": 4000 if 1 not in degrees else 0,
        "explicit_degree_one_component_credit": -500 * min(degrees.count(1), 6),
        "explicit_degree_zero_component_credit": -100 * min(degrees.count(0), 12),
    }
    return profile, terms, int(sum(terms.values()))


def replace_availability_terms(terms, explicit_degrees):
    """Re-score using only a supplied, provenance-tracked explicit curve set."""
    answer = dict(terms)
    for key in ("no_explicit_degree_one_component", "explicit_degree_one_component_credit",
                "explicit_degree_zero_component_credit"):
        answer.pop(key)
    answer.update({
        "no_explicit_degree_one_component": 4000 if 1 not in explicit_degrees else 0,
        "explicit_degree_one_component_credit": -500 * min(explicit_degrees.count(1), 6),
        "explicit_degree_zero_component_credit": -100 * min(explicit_degrees.count(0), 12),
    })
    return answer, int(sum(answer.values()))


marking = json.loads(MARKING.read_text())
frontier = json.loads(FRONTIER.read_text())
source = load_matrix(FRAME)
g_source = block_diagonal_matrix(U2, -source)
target_container = (
    marking["target_fibres_in_root_adapted_hub"]
    if "target_fibres_in_root_adapted_hub" in marking else marking["target_fibres_in_child"]
)
targets = {name: vector(ZZ, value) for name, value in target_container.items()}
target = targets[TARGET_KEY]
pinned = targets["pinned_R17"]
source_fibre = vector(ZZ, [1, 0] + [0] * 17)
source_root_record = marking["root_data"] if "root_data" in marking else marking["child"]["root_data"]
source_root_data = tuple(map(ZZ, source_root_record))
source_simple_roots = identity_matrix(ZZ, 17)[:int(source_root_data[0])]
source_cartan = source[:int(source_root_data[0]), :int(source_root_data[0])]
source_component_curves = [
    (f"old_{SOURCE_LABEL}_component_{node}",
     vector(ZZ, [0, 0] + [-ZZ(i == node) for i in range(17)]))
    for node in range(int(source_root_data[0]))
]
source_component_curves += [
    (f"old_{SOURCE_LABEL}_affine_{index}",
     vector(ZZ, [1, 0] + list(top) + [0] * (17 - int(source_root_data[0]))))
    for index, top in enumerate(highest_roots(source_cartan))
]
source_known_curves = source_component_curves + [
    (f"old_{SOURCE_LABEL}_zero", vector(ZZ, [-1, 1] + [0] * 17))
]

# The stored suffix ray is in the exact D13 basis but its root coordinates use
# the suffix chamber.  Reorient only the D13 root basis before pricing it; this
# is a full integral isometry and leaves the absolute fibre class unchanged.
target_direct, target_direct_action = reduce_component_chamber(
    target, source, source_simple_roots
)
target_nef = nef_profile(target_direct, source, source_simple_roots)
assert target_nef["nef"] and not negative_horizontal_walls(target_direct, source)
direct_child_raw, _ = child_frame(g_source, target_direct)
_, _, direct_child_root_data, _ = root_adaptation(direct_child_raw)
direct_profile, direct_terms, direct_score = cost_score(
    target_direct, source, source_simple_roots, direct_child_root_data, target_nef
)
direct_known_degrees = [int(curve * g_source * target_direct) for _, curve in source_known_curves]
direct_known_terms, direct_known_score = replace_availability_terms(direct_terms, direct_known_degrees)

if args.mode in ("e8e6", "d13", "d12", "a5"):
    first_edges = frontier["ranked_candidates"]
else:
    first_edges = [row for row in frontier["retained_candidates"]
                   if row["declared_curve_nef_gate"] == "PASS"]
first_edges = first_edges[:args.limit or None]
records = []
first_edge_operational_minimum = None
counts = {"first_edges": len(first_edges), "first_exact_nef": 0, "explicit_zeros": 0,
          "generic_child_rebuilds": 0,
          "first_candidate_operational_lower_bound_pruned": 0,
          "first_edge_operational_bound_pruned": 0,
          "nonprimitive_explicit_root_spans": 0,
          "return_component_nef": 0, "return_exact_nef": 0,
          "exit_component_nef": 0, "exit_exact_nef": 0}
for first_index, candidate in enumerate(first_edges):
    raw = candidate.get("source_neighbor_record", candidate)
    fibre = vector(ZZ, raw.get("fiber", raw.get("fibre")))
    if args.mode in ("e8e6", "d13", "d12", "a5"):
        physical_degrees = candidate["component_pairings"] + candidate["affine_pairings"]
    else:
        physical_degrees = candidate["explicit_curve_degrees"]["physical_old_A11_fibre_components"]
    explicit_zeros = [item for item, degree in zip(source_component_curves, physical_degrees)
                      if degree == 1]
    if not explicit_zeros:
        continue

    # Price the inherited first edge once from exact data already stored in
    # the marked frontier.  The score is independent of which degree-one
    # component is selected as zero.  Expensive independent nef/C2 replay and
    # child rebuilding are deferred until this exact score plus the two
    # unavoidable horizontal floors can still meet the requested bound.
    frontier_nef = {
        "component_pairings": candidate["component_pairings"],
        "affine_pairings": candidate["affine_pairings"],
    }
    first_profile, first_terms, first_score = cost_score(
        fibre, source, source_simple_roots, candidate["child"]["root_data"], frontier_nef
    )
    inherited_first_degrees = [int(curve * g_source * fibre)
                               for _, curve in source_known_curves]
    inherited_first_terms, inherited_first_score = replace_availability_terms(
        first_terms, inherited_first_degrees
    )
    candidate_terms = candidate.get("equation_cost_terms", {})
    first_coset_penalty = int(candidate_terms.get("finite_index_five_coset", 0)) + int(
        candidate_terms.get("missing_rank_direction", 0)
    )
    inherited_first_terms["target_coset_mod_already_explicit_sections"] = first_coset_penalty
    inherited_first_score += first_coset_penalty
    inherited_first_operational = max(
        inherited_first_terms["horizontal_degree"], inherited_first_score
    )
    first_edge_summary = {
        "candidate_id": candidate["candidate_id"],
        "explicit_zero_choice_count": len(explicit_zeros),
        "profile": first_profile,
        "inherited_explicit_terms": inherited_first_terms,
        "inherited_explicit_score": inherited_first_score,
        "operational_score": inherited_first_operational,
    }
    if (first_edge_operational_minimum is None
            or inherited_first_operational < first_edge_operational_minimum["operational_score"]):
        first_edge_operational_minimum = first_edge_summary
    if (args.max_operational_total is not None
            and inherited_first_operational + 1000 >= args.max_operational_total):
        counts["explicit_zeros"] += len(explicit_zeros)
        counts["first_edge_operational_bound_pruned"] += len(explicit_zeros)
        counts["first_candidate_operational_lower_bound_pruned"] += 1
        continue

    first_nef = nef_profile(fibre, source, source_simple_roots)
    if not first_nef["nef"] or negative_horizontal_walls(fibre, source):
        continue
    counts["first_exact_nef"] += 1
    # Rebuild the marked child uniformly.  Compact frontier rows deliberately
    # omit bases for a few high-root states, while the fibre itself determines
    # the same primitive U and full child lattice exactly.
    raw_frame, raw_transition = child_frame(g_source, fibre)
    counts["generic_child_rebuilds"] += 1
    g_raw = block_diagonal_matrix(U2, -raw_frame)
    assert raw_transition * g_source * raw_transition.transpose() == g_raw
    raw_inverse = raw_transition.inverse().change_ring(ZZ)

    for zero_name, curve_source in explicit_zeros:
        counts["explicit_zeros"] += 1
        section = curve_source * raw_inverse
        assert section * g_raw * section == -2 and section * g_raw * source_fibre == 1
        mate = section + source_fibre
        complement = matrix(ZZ, [list(source_fibre * g_raw), list(mate * g_raw)]).right_kernel_matrix()
        split = matrix(ZZ, [list(source_fibre), list(mate)] + list(complement.rows()))
        explicit_raw = -(complement * g_raw * complement.transpose())
        explicit_frame, explicit_adaptation, explicit_root_data, explicit_simple_roots = root_adaptation(
            explicit_raw
        )
        assert tuple(explicit_root_data) == tuple(candidate["child"]["root_data"])
        if abs(matrix(ZZ, deterministic_simple_roots(explicit_raw)).row_module(ZZ).index_in(
                matrix(ZZ, deterministic_simple_roots(explicit_raw)).row_module(ZZ).saturation()
        )) != 1:
            counts["nonprimitive_explicit_root_spans"] += 1
        explicit_reframe = block_diagonal_matrix(identity_matrix(ZZ, 2), explicit_adaptation) * split
        source_to_explicit = explicit_reframe * raw_transition
        explicit_inverse = source_to_explicit.inverse().change_ring(ZZ)
        g_explicit = block_diagonal_matrix(U2, -explicit_frame)
        assert source_to_explicit * g_source * source_to_explicit.transpose() == g_explicit

        # Root adaptation chooses an abstract child chamber.  Reorient that
        # basis to the known nef original-source return ray before applying gates.
        return_initial = source_fibre * explicit_inverse
        return_fibre, explicit_chamber_action = reduce_component_chamber(
            return_initial, explicit_frame, explicit_simple_roots
        )
        source_to_explicit = explicit_chamber_action.inverse().change_ring(ZZ) * source_to_explicit
        explicit_inverse = source_to_explicit.inverse().change_ring(ZZ)
        assert return_fibre == source_fibre * explicit_inverse
        assert source_to_explicit * g_source * source_to_explicit.transpose() == g_explicit
        return_nef = nef_profile(return_fibre, explicit_frame, explicit_simple_roots)
        if not return_nef["nef"]:
            continue
        counts["return_component_nef"] += 1
        return_walls = negative_horizontal_walls(return_fibre, explicit_frame)
        if return_walls:
            continue
        counts["return_exact_nef"] += 1
        returned_raw, returned_transition_raw = child_frame(g_explicit, return_fibre)
        returned_frame, returned_adaptation, returned_root_data, returned_simple_roots = root_adaptation(
            returned_raw
        )
        if tuple(returned_root_data) != tuple(source_root_data):
            continue
        initial_return = block_diagonal_matrix(identity_matrix(ZZ, 2), returned_adaptation) * returned_transition_raw
        initial_inverse = initial_return.inverse().change_ring(ZZ)
        target_explicit = target * explicit_inverse
        anchor_initial = target_explicit * initial_inverse
        _, chamber_action = reduce_component_chamber(
            anchor_initial, returned_frame, returned_simple_roots
        )
        return_transition = chamber_action.inverse().change_ring(ZZ) * initial_return
        return_inverse = return_transition.inverse().change_ring(ZZ)
        g_returned = block_diagonal_matrix(U2, -returned_frame)
        assert return_transition * g_explicit * return_transition.transpose() == g_returned
        exit_fibre = target_explicit * return_inverse
        exit_nef = nef_profile(exit_fibre, returned_frame, returned_simple_roots)
        if not exit_nef["nef"]:
            continue
        counts["exit_component_nef"] += 1
        exit_walls = negative_horizontal_walls(exit_fibre, returned_frame)
        if exit_walls:
            continue
        counts["exit_exact_nef"] += 1
        exit_child_raw, exit_transition_raw = child_frame(g_returned, exit_fibre)
        exit_child_frame, exit_adaptation, exit_child_root_data, exit_child_simple_roots = root_adaptation(
            exit_child_raw
        )
        if tuple(exit_child_root_data) != tuple(direct_child_root_data):
            continue
        exit_transition = block_diagonal_matrix(identity_matrix(ZZ, 2), exit_adaptation) * exit_transition_raw
        g_exit_child = block_diagonal_matrix(U2, -exit_child_frame)
        assert exit_transition * g_returned * exit_transition.transpose() == g_exit_child

        return_profile, return_terms, return_score = cost_score(
            return_fibre, explicit_frame, explicit_simple_roots, returned_root_data, return_nef
        )
        exit_profile, exit_terms, exit_score = cost_score(
            exit_fibre, returned_frame, returned_simple_roots, exit_child_root_data, exit_nef
        )
        total_transition = return_transition * source_to_explicit
        total_inverse = total_transition.inverse().change_ring(ZZ)
        pinned_returned = pinned * total_inverse
        inherited_return_degrees = [int((curve * explicit_inverse) * g_explicit * return_fibre)
                                    for _, curve in source_known_curves]
        inherited_exit_degrees = [int((curve * total_inverse) * g_returned * exit_fibre)
                                  for _, curve in source_known_curves]
        inherited_curve_gate = all(
            degree >= 0
            for degrees in (inherited_first_degrees, inherited_return_degrees, inherited_exit_degrees)
            for degree in degrees
        )
        inherited_return_terms, inherited_return_score = replace_availability_terms(
            return_terms, inherited_return_degrees
        )
        inherited_exit_terms, inherited_exit_score = replace_availability_terms(
            exit_terms, inherited_exit_degrees
        )
        inherited_total_score = inherited_first_score + inherited_return_score + inherited_exit_score
        inherited_operational_scores = [
            inherited_first_operational,
            max(inherited_return_terms["horizontal_degree"], inherited_return_score),
            max(inherited_exit_terms["horizontal_degree"], inherited_exit_score),
        ]
        inherited_operational_total = sum(inherited_operational_scores)
        total_score = first_score + return_score + exit_score
        record = {
            "first_edge_candidate_id": candidate["candidate_id"],
            "first_edge_child": candidate["child"],
            "explicit_zero_curve": zero_name,
            "first_edge_nef_audit": first_nef,
            "first_edge_exact_negative_horizontal_walls": [],
            "first_edge_profile": first_profile,
            "first_edge_score": first_score,
            "first_edge_terms": first_terms,
            "return_fibre_in_explicit_child": entries(return_fibre),
            "return_nef_audit": return_nef,
            "return_exact_negative_horizontal_walls": return_walls,
            "return_profile": return_profile,
            "return_score": return_score,
            "return_terms": return_terms,
            f"exit_{EXIT_LABEL}_fibre_in_returned_{SOURCE_LABEL}": entries(exit_fibre),
            "exit_nef_audit": exit_nef,
            "exit_exact_negative_horizontal_walls": exit_walls,
            "exit_profile": exit_profile,
            "exit_score": exit_score,
            "exit_terms": exit_terms,
            "q_sequence": [int(fibre[0] * fibre[1]), int(return_fibre[0] * return_fibre[1]),
                           int(exit_fibre[0] * exit_fibre[1])],
            "old_fibre_degrees": [int(fibre[1]), int(return_fibre[1]), int(exit_fibre[1])],
            "total_equation_cost_score": total_score,
            "strict_improvement_over_direct": total_score < direct_score,
            "inherited_explicit_curve_gate": "PASS" if inherited_curve_gate else "REJECT",
            "inherited_explicit_curve_degrees": {
                "names": [name for name, _ in source_known_curves],
                "first_edge": inherited_first_degrees,
                "return": inherited_return_degrees,
                "exit": inherited_exit_degrees,
            },
            "inherited_explicit_equation_cost": {
                "first_edge_terms": inherited_first_terms,
                "first_edge_score": inherited_first_score,
                "return_terms": inherited_return_terms,
                "return_score": inherited_return_score,
                "exit_terms": inherited_exit_terms,
                "exit_score": inherited_exit_score,
                "total_score": inherited_total_score,
                "operational_edge_scores": inherited_operational_scores,
                "operational_total_score": inherited_operational_total,
                "strict_improvement_over_direct": (
                    inherited_curve_gate and inherited_operational_total
                    < max(direct_known_terms["horizontal_degree"], direct_known_score)
                ),
            },
            "explicit_child_root_data": list(map(int, explicit_root_data)),
            "explicit_child_frame": rows(explicit_frame),
            "source_to_explicit_child_basis": rows(source_to_explicit),
            "explicit_child_to_source_basis": rows(explicit_inverse),
            f"explicit_child_to_returned_{SOURCE_LABEL}_basis": rows(return_transition),
            f"returned_{SOURCE_LABEL}_to_explicit_child_basis": rows(return_inverse),
            f"source_to_returned_{SOURCE_LABEL}_basis": rows(total_transition),
            f"returned_{SOURCE_LABEL}_to_source_basis": rows(total_inverse),
            "transport_determinant": int(total_transition.det()),
            f"pinned_R17_fibre_in_returned_{SOURCE_LABEL}": entries(pinned_returned),
            "returned_frame": rows(returned_frame),
            "exit_child_root_data": list(map(int, exit_child_root_data)),
            "exit_child_frame": rows(exit_child_frame),
            f"returned_{SOURCE_LABEL}_to_exit_{EXIT_LABEL}_basis": rows(exit_transition),
            f"exit_{EXIT_LABEL}_to_returned_{SOURCE_LABEL}_basis": rows(
                exit_transition.inverse().change_ring(ZZ)
            ),
        }
        # Preserve the original D13 artifact schema consumed by the standalone
        # promoted-route certificate while also exposing the generic predicate.
        if args.mode == "d13":
            record["strict_improvement_over_direct_q24"] = record["strict_improvement_over_direct"]
        records.append(record)
    if (first_index + 1) % 20 == 0:
        print("{}ZEROLOOPPROGRESS|first={}|zeros={}|returns={}|exits={}".format(
            SOURCE_LABEL,
            first_index + 1, counts["explicit_zeros"], counts["return_exact_nef"], counts["exit_exact_nef"]
        ), flush=True)

exit_key = f"exit_{EXIT_LABEL}_fibre_in_returned_{SOURCE_LABEL}"
records.sort(key=lambda row: (row["total_equation_cost_score"], max(abs(x) for x in row[exit_key]),
                              row["first_edge_candidate_id"]["q"], row["first_edge_candidate_id"]["orbit_index"],
                              row["explicit_zero_curve"]))
if args.mode in ("e8e6", "a11", "a5"):
    records.sort(key=lambda row: (
        row["inherited_explicit_curve_gate"] != "PASS",
        row["inherited_explicit_equation_cost"]["operational_total_score"],
        row["inherited_explicit_equation_cost"]["total_score"],
        row["total_equation_cost_score"],
        row["first_edge_candidate_id"]["orbit_index"],
        row["explicit_zero_curve"],
    ))
inputs = (MARKING, FRAME, FRONTIER)
payload = {
    "schema": f"elkies-k3.h3-{SOURCE_LABEL.lower()}-zero-changing-{EXIT_LABEL.lower()}-presentations.v1",
    "status": f"PASS_EXACT_{SOURCE_LABEL}_ZERO_CHANGING_{EXIT_LABEL}_PRESENTATION_SEARCH",
    "search_parameters": {"frontier": str(FRONTIER.relative_to(ROOT)), "limit": args.limit,
                          "retained": args.retain,
                          "max_operational_total": args.max_operational_total},
    "counts": counts,
    "first_edge_operational_minimum": first_edge_operational_minimum,
    f"direct_{DIRECT_LABEL}": {"profile": direct_profile, "terms": direct_terms, "score": direct_score,
                   "stored_suffix_fibre": entries(target), "component_nef_coordinates": entries(target_direct),
                   "component_chamber_action": rows(target_direct_action),
                   "inherited_explicit_curve_degrees": direct_known_degrees,
                   "inherited_explicit_terms": direct_known_terms,
                   "inherited_explicit_score": direct_known_score,
                   "inherited_explicit_operational_score": max(
                       direct_known_terms["horizontal_degree"], direct_known_score
                   ),
                   "child_root_data": list(map(int, direct_child_root_data))},
    "raw_strict_winner_count": sum(row["strict_improvement_over_direct"] for row in records),
    "strict_winner_count": sum(
        row["inherited_explicit_equation_cost"]["strict_improvement_over_direct"]
        if args.mode in ("e8e6", "a11", "a5") else row["strict_improvement_over_direct"]
        for row in records
    ),
    "ranked_presentations": records[:args.retain],
    "proof_boundary": (
        "The search is exhaustive only over the supplied compact frontier. Every retained first, return, "
        "and fixed exit fibre passes exact component, affine, all-section, and finite horizontal-wall "
        "gates. Full standalone certificates and composition with the later pinned route remain separate. "
        "Nonprimitive full root spans are handled by a saturated unimodular frame with the actual "
        "simple-root lattice retained as an embedded sublattice. "
        "Equation-cost scores are deterministic planning estimates, not measured runtimes."
    ),
    "inputs": {"paths": [str(path.relative_to(ROOT)) for path in inputs],
               "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                           for path in inputs}},
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
best = records[0] if records else None
best_score = None if best is None else (
    best["inherited_explicit_equation_cost"]["operational_total_score"]
    if args.mode in ("e8e6", "a11", "a5")
    else best["total_equation_cost_score"]
)
print("{}ZEROLOOP|records={}|strict={}|direct={}|best={}|score={}|status={}|output={}".format(
    SOURCE_LABEL, len(records), payload["strict_winner_count"], direct_score,
    None if best is None else "q{}o{}:{}".format(best["first_edge_candidate_id"]["q"],
                                                 best["first_edge_candidate_id"]["orbit_index"],
                                                 best["explicit_zero_curve"]),
    best_score, payload["status"], OUTPUT
), flush=True)
