#!/usr/bin/env sage -python
"""Fully certify a compiler-cheap q4/q6 exit from explicit-zero A5+A5.

status: ACTIVE_PROOF
claim: exact marked U, finite-wall nefness, roots and NS transports for one exit
inputs: explicit-zero q4/q6 gate and equation-cost artifacts
outputs: caller-selected exact JSON certificate and child frame
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
ZERO_FRAME = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--q", type=int, default=4)
parser.add_argument("--orbit", type=int, default=32)
parser.add_argument("--gate", type=Path, default=GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4q6-explicit-curve-gate.json")
parser.add_argument("--scores", type=Path, default=GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4q6-equation-cost.json")
parser.add_argument("--output", type=Path, default=GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4-orbit32-lattice-certificate.json")
parser.add_argument("--frame-output", type=Path, default=GENERATED / "elkies-k3-h3-a5a5-explicit-zero-q4-orbit32-frame.txt")
args = parser.parse_args()
GATE = args.gate.resolve()
SCORES = args.scores.resolve()
OUTPUT = args.output.resolve()
FRAME_OUTPUT = args.frame_output.resolve()
INPUTS = (GATE, SCORES, ZERO_FRAME)
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(ZZ, [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
                       if line.strip() and not line.lstrip().startswith("#")])


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


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
        candidates = [item for item in roots if all(value >= 0 for value in item)
                      and all(index in component or item[index] == 0 for index in range(cartan.nrows()))]
        result.append(max(candidates, key=lambda item: sum(item)))
    return tuple(result)


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
                (w * frame * w // (2 * degree)) * m + degree * k - w * frame * x
            )
            if intersection < 0:
                walls.append({
                    "old_fibre_degree": int(m),
                    "curve": [int(k), int(m)] + list(map(int, x)),
                    "intersection": int(intersection),
                })
    return sorted(walls, key=lambda item: (item["intersection"], item["curve"]))


def matrix_rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def bezout_vector_for_pairing(ns, fiber):
    pairings = list(ns * fiber)
    current = ZZ(0)
    result = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(pairings):
        if value == 0:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        result = [left * entry for entry in result]
        result[index] += right
        current = divisor
    if abs(current) != 1:
        return None
    if current == -1:
        result = [-entry for entry in result]
    return vector(ZZ, result)


def child_frame(ns, fiber, determinant):
    mate = bezout_vector_for_pairing(ns, fiber)
    assert mate is not None
    mate_square = ZZ(mate * ns * mate)
    assert mate_square % 2 == 0
    mate -= (mate_square // 2) * fiber
    kernel = matrix(ZZ, [list(fiber * ns), list(mate * ns)]).right_kernel_matrix()
    child = -(kernel * ns * kernel.transpose())
    neighbor_basis = matrix(ZZ, [list(fiber), list(mate)] + [list(row) for row in kernel.rows()])
    assert child.is_positive_definite() and child.det() == determinant
    assert abs(neighbor_basis.det()) == 1
    assert neighbor_basis * ns * neighbor_basis.transpose() == block_diagonal_matrix(U2, -child)
    return child, neighbor_basis


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


def component_name(cartan, component):
    block = cartan.matrix_from_rows_and_columns(component, component)
    rank = block.nrows()
    determinant = abs(ZZ(block.det()))
    root_count = ZZ(pari(block).qfminim(2)[0])
    if determinant == rank + 1 and root_count == rank * (rank + 1):
        return f"A{rank}"
    if rank >= 4 and determinant == 4 and root_count == 2 * rank * (rank - 1):
        return f"D{rank}"
    return {(6, 3, 72): "E6", (7, 2, 126): "E7", (8, 1, 240): "E8"}.get(
        (rank, determinant, root_count), f"R{rank}d{determinant}n{root_count}"
    )


def ade_name(cartan):
    return "+".join(component_name(cartan, component) for component in components(cartan))


def root_adaptation(child):
    unused, root_basis, invariants = roots_and_data(child)
    root_rank = invariants[0]
    if root_rank == 0:
        adapted_basis = matrix(ZZ, pari(child).qflllgram()).transpose()
        return adapted_basis * child * adapted_basis.transpose(), adapted_basis, invariants
    smith, smith_left, smith_right = root_basis.smith_form()
    assert smith == smith_left * root_basis * smith_right
    assert tuple(abs(smith[index, index]) for index in range(root_rank)) == (1,) * root_rank
    simple, cartan = deterministic_simple_roots(child)
    completion = smith_right.inverse()
    adapted_basis = simple.stack(completion[root_rank:])
    adapted = adapted_basis * child * adapted_basis.transpose()
    coupling = adapted[:root_rank, root_rank:]
    height = adapted[root_rank:, root_rank:] - coupling.transpose() * cartan.inverse() * coupling
    scale = lcm(entry.denominator() for entry in height.list())
    quotient_change = block_diagonal_matrix(
        identity_matrix(ZZ, root_rank), matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram()).transpose()
    )
    adapted_basis = quotient_change * adapted_basis
    adapted = adapted_basis * child * adapted_basis.transpose()
    assert abs(adapted_basis.det()) == 1
    return adapted, adapted_basis, invariants


gate = json.loads(GATE.read_text())
scores = json.loads(SCORES.read_text())
zero_frame = json.loads(ZERO_FRAME.read_text())
assert scores["status"] in (
    "PASS_EXACT_A5A5_EXPLICIT_ZERO_EQUATION_COST_SCORING",
    "PASS_EXACT_A5A5_EXPLICIT_ZERO_LARGE_Q_COST_SCORING",
)
selection = next(item for item in scores["ranked_candidates"]
                 if item["candidate_id"] == {"q": args.q, "old_fibre_degree": 2, "orbit_index": args.orbit})
assert selection.get("full_declared_nef_gate", "PASS") == "PASS"
gate_item = next(item for item in gate["survivors"]
                 if int(item["candidate_id"]["q"]) == args.q and int(item["candidate_id"]["orbit_index"]) == args.orbit)
raw = gate_item.get("source_neighbor_record", gate_item)

parent_path = ROOT / zero_frame["selected_frame_output"]
parent = load_matrix(parent_path)
g_parent = block_diagonal_matrix(U2, -parent)
fibre = vector(ZZ, raw["fiber"] if "fiber" in raw else raw["fibre_in_parent"])
witness = vector(ZZ, raw["witness"] if "witness" in raw else fibre[2:])
old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
assert args.q % 2 == 0
assert fibre == vector(ZZ, [args.q // 2, 2] + list(witness))
assert fibre * g_parent * fibre == 0 and gcd(tuple(g_parent * fibre)) == 1
assert fibre * g_parent * old_fibre == 2
assert fibre * g_parent * old_zero == args.q // 2 - 2

labels = tuple(map(ZZ, raw["dominant_labels"]))
affine = tuple(ZZ(2 - top * vector(ZZ, labels)) for top in highest_roots(parent[:10, :10]))
if "parent_affine_component_pairings" in selection:
    assert tuple(map(int, affine)) == tuple(selection["parent_affine_component_pairings"])
assert all(value >= 0 for value in labels + affine)

# Every section is indexed by a frame vector x. Its intersection with this
# degree-two fibre is (x-w/2)^2-2, so one exact CVP checks all section walls.
center = vector(QQ, witness) / 2
closest = vector(ZZ, next(IntegralLattice(parent).enumerate_close_vectors(center)))
closest_distance = (closest - center) * parent * (closest - center)
assert closest_distance >= 2
assert (witness * parent * witness - 2) % 4 != 0
horizontal_walls = negative_horizontal_walls(fibre, parent)
assert not horizontal_walls

if "child_root_adapted_basis" in raw:
    neighbor_basis = matrix(ZZ, raw["neighbor_basis"])
    adapted_basis = matrix(ZZ, raw["child_root_adapted_basis"])
    child = matrix(ZZ, raw["child_root_adapted_frame"])
    child_root_data = tuple(map(int, raw["child_root_data"]))
    child_ade = raw["child_ade"]
    child_mw_rank = int(raw["child_mw_rank"])
else:
    unadapted_child, neighbor_basis = child_frame(g_parent, fibre, abs(ZZ(parent.det())))
    child, adapted_basis, child_root_data = root_adaptation(unadapted_child)
    child_mw_rank = 17 - int(child_root_data[0])
    child_ade = ade_name(child[:child_root_data[0], :child_root_data[0]]) if child_root_data[0] else "rootless"
transition = block_diagonal_matrix(identity_matrix(ZZ, 2), adapted_basis) * neighbor_basis
inverse = transition.inverse().change_ring(ZZ)
g_child = block_diagonal_matrix(U2, -child)
assert abs(transition.det()) == 1
assert transition * g_parent * transition.transpose() == g_child
assert inverse * g_child * inverse.transpose() == g_parent
if "child" in selection:
    assert child_ade == selection["child"]["ade"]
    assert child_mw_rank == int(selection["child"]["mw_rank"])
    assert child_root_data == tuple(map(int, selection["child"]["root_data"]))

equation_a11_to_parent = matrix(ZZ, zero_frame["selected"]["equation_A11_to_explicit_zero_basis"])
equation_a11_to_child = transition * equation_a11_to_parent
child_to_equation_a11 = equation_a11_to_child.inverse().change_ring(ZZ)
assert abs(equation_a11_to_child.det()) == 1
FRAME_OUTPUT.write_text(f"# q{args.q} orbit{args.orbit} {child_ade}/MW{child_mw_rank} compiler-cost child\n" +
                        "\n".join(" ".join(map(str, row)) for row in child.rows()) + "\n")

payload = {
    "schema": "elkies-k3.h3-a5a5-explicit-zero-candidate-lattice-certificate.v1",
    "status": "PASS_EXACT_A5A5_EXPLICIT_ZERO_CANDIDATE_LATTICE_CERTIFICATE",
    "inputs": {"paths": [str(path.relative_to(ROOT)) for path in INPUTS] + [str(parent_path.relative_to(ROOT))], "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in INPUTS + (parent_path,)}},
    "selection": selection,
    "edge": {"q": args.q, "factorization": [args.q // 2, 2], "orbit_index": args.orbit, "primitive_nef_isotropic_fibre": entries(fibre),
             "component_pairings": list(map(int, labels)), "affine_pairings": list(map(int, affine)),
             "closest_section_vector": entries(closest), "closest_section_distance": str(closest_distance),
             "minimum_section_intersection": str(closest_distance - 2), "bisection_parity_exclusion": True,
             "exact_negative_horizontal_walls": horizontal_walls,
             "exact_horizontal_nef_gate": True, "nef": True},
    "marked_U": {"fibre_in_parent": entries(transition.row(0)), "isotropic_mate_in_parent": entries(transition.row(1)),
                 "zero_in_parent": entries(transition.row(1) - transition.row(0)), "gram": [[0, 1], [1, 0]]},
    "child": {"ade": child_ade, "mw_rank": child_mw_rank, "root_data": list(map(int, child_root_data)),
              "frame": rows(child), "frame_output": str(FRAME_OUTPUT.relative_to(ROOT)), "frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest()},
    "transport": {"parent_to_child_basis": rows(transition), "child_to_parent_basis": rows(inverse),
                  "equation_A11_to_child_basis": rows(equation_a11_to_child), "child_to_equation_A11_basis": rows(child_to_equation_a11),
                  "forward_determinant": int(transition.det()), "inverse_determinant": int(inverse.det())},
    "route_status": "Certified cheap edge only; no continuation to pinned R17 is yet certified, so the lifting target remains unchanged.",
    "proof_boundary": "Exact component, all-section, and complete finite horizontal-wall nef gates; marked U, roots, and unimodular transports. RR and equation cost remain planning estimates.",
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("A5A5CERT|q={}|orbit={}|PO={}|RR={}|closest={}|child={}/MW{}|root={}|det_fwd={}|det_inv={}|nef=1|status={}".format(
    args.q, args.orbit, selection["horizontal"]["P_dot_O"], selection["expected_RR_ambient"], closest_distance,
    child_ade, child_mw_rank, ",".join(map(str, child_root_data)),
    transition.det(), inverse.det(), payload["status"]), flush=True)
print(f"OUTPUT|{OUTPUT.resolve()}", flush=True)
