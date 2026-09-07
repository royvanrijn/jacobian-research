#!/usr/bin/env sage -python
"""Certify the best declared-curve-nef A11 equation-cost exit, q8 orbit849.

The candidate is selected by the separate exhaustive cost scorer.  This
checker independently proves primitivity, full nefness, the marked U split,
unimodular forward/inverse NS transport, complete child roots, and a
root-adapted child frame.  It does not claim a continuation to pinned R17 or
an equation lift.
"""

import argparse
import hashlib
import itertools
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
U2 = matrix(ZZ, ((0, 1), (1, 0)))

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--orbit", type=int, default=849)
parser.add_argument(
    "--output",
    type=Path,
    default=GENERATED / "elkies-k3-h3-a11-q8-orbit849-lattice-certificate.json",
)
parser.add_argument(
    "--frame-output",
    type=Path,
    default=GENERATED / "elkies-k3-h3-a11-q8-orbit849-frame.txt",
)
args = parser.parse_args()

SCORES = GENERATED / "elkies-k3-h3-a11-equation-cost-neighbors.json"
NEIGHBORS = LOCAL / "q24-a11-orbit64-q8-all.json"
INPUTS = (SCORES, NEIGHBORS)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")


def load_gram(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


def rational_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def connected_components(cartan):
    unseen = set(range(cartan.nrows()))
    result = []
    while unseen:
        todo = [min(unseen)]
        unseen.remove(todo[0])
        component = []
        while todo:
            index = todo.pop()
            component.append(index)
            for other in tuple(unseen):
                if cartan[index, other]:
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


def neighbor_frame(ns, fibre):
    mate = bezout_vector_for_pairing(ns, fibre)
    mate -= (mate * ns * mate // 2) * fibre
    complement = matrix(
        ZZ, [list(fibre * ns), list(mate * ns)]
    ).right_kernel_matrix()
    basis = matrix(ZZ, [list(fibre), list(mate)] + [list(item) for item in complement.rows()])
    assert fibre * ns * mate == 1 and mate * ns * mate == 0
    assert abs(basis.det()) == 1
    child = -(complement * ns * complement.transpose())
    assert basis * ns * basis.transpose() == block_diagonal_matrix(U2, -child)
    return child, basis, mate


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    half = [vector(ZZ, column) for column in matrix(ZZ, result[2]).columns()]
    roots = tuple(half + [-item for item in half])
    root_basis = matrix(ZZ, [list(item) for item in roots]).row_module().basis_matrix()
    root_gram = root_basis * gram * root_basis.transpose()
    return roots, root_basis, (
        root_basis.rank(),
        count,
        abs(ZZ(root_gram.det())),
    )


def deterministic_simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    positive = [item for item in roots if next(value for value in item if value != 0) > 0]
    positive_set = {tuple(item) for item in positive}
    simple = [
        item
        for item in positive
        if not any(tuple(item - left) in positive_set for left in positive)
    ]
    simple = matrix(ZZ, [list(item) for item in simple])
    assert simple.nrows() == simple.rank() == data[0]
    cartan = simple * gram * simple.transpose()
    return simple, cartan


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
    return "+".join(component_name(cartan, component) for component in connected_components(cartan))


def root_adaptation(child):
    unused, root_basis, invariants = roots_and_data(child)
    root_rank = invariants[0]
    smith, left, right = root_basis.smith_form()
    assert smith == left * root_basis * right
    assert tuple(abs(smith[index, index]) for index in range(root_rank)) == (1,) * root_rank
    simple, cartan = deterministic_simple_roots(child)
    completion = right.inverse()
    adapted_basis = simple.stack(completion[root_rank:])
    assert abs(adapted_basis.det()) == 1
    adapted = adapted_basis * child * adapted_basis.transpose()
    coupling = adapted[:root_rank, root_rank:]
    tail = adapted[root_rank:, root_rank:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling
    scale = lcm(value.denominator() for value in height.list())
    lll = matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram())
    quotient_change = block_diagonal_matrix(identity_matrix(ZZ, root_rank), lll.transpose())
    adapted_basis = quotient_change * adapted_basis
    adapted = adapted_basis * child * adapted_basis.transpose()
    coupling = adapted[:root_rank, root_rank:]
    tail = adapted[root_rank:, root_rank:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling
    assert abs(adapted_basis.det()) == 1
    return adapted, adapted_basis, height, cartan


def closest_root_squared(root, center):
    standard = CartanMatrix(["A", 11])
    ambient = matrix(ZZ, 11, 12)
    for index in range(11):
        ambient[index, index] = 1
        ambient[index, index + 1] = -1
    assert ambient * ambient.transpose() == standard
    isometry = matrix(ZZ, pari(standard).qfisom(pari(root)))
    assert abs(isometry.det()) == 1
    assert isometry.transpose() * root * isometry == standard
    target = ambient.transpose() * isometry.inverse() * center
    choices = []
    for value in target:
        lower, upper = floor(value), ceil(value)
        choices.append((lower,) if lower == upper else (lower, upper))
    minimum = None
    for candidate in itertools.product(*choices):
        if sum(candidate) != 0:
            continue
        difference = vector(QQ, candidate) - target
        value = difference * difference
        minimum = value if minimum is None or value < minimum else minimum
    assert minimum is not None
    return minimum


def section_distance_profile(frame, witness):
    root_rank = 11
    root = frame[:root_rank, :root_rank]
    root_mw = frame[:root_rank, root_rank:]
    height = frame[root_rank:, root_rank:] - frame[root_rank:, :root_rank] * root.inverse() * root_mw
    z = vector(ZZ, witness[root_rank:])
    denominator = lcm(value.denominator() for value in height.list())
    scaled_height = (denominator * height).change_ring(ZZ)
    short = pari(scaled_height).qfminim(8 * denominator - 1)
    half = matrix(ZZ, short[2]).transpose().rows()
    candidates = [vector(ZZ, [0] * len(z))] + list(half) + [-item for item in half]
    candidates = [
        item
        for item in candidates
        if all((item[index] - z[index]) % 2 == 0 for index in range(len(z)))
        and QQ(item * scaled_height * item) / (4 * denominator) < 2
    ]
    root_coordinate = vector(QQ, witness[:root_rank])
    distances = []
    for item in candidates:
        quotient_difference = vector(QQ, (item + z) / 2) - vector(QQ, z) / 2
        center = root_coordinate / 2 - root.inverse() * root_mw * quotient_difference
        distances.append(
            closest_root_squared(root, center)
            + QQ(item * scaled_height * item) / (4 * denominator)
        )
    return tuple(sorted(distances))


scores = json.loads(SCORES.read_text())
neighbors = json.loads(NEIGHBORS.read_text())
assert scores["status"] == "PASS_EXACT_A11_EQUATION_COST_SCORING"
best = scores["best_candidate"]
if int(best["candidate_id"]["orbit_index"]) != args.orbit:
    best = next(
        item for item in scores["retained_candidates"] + list(scores["named_construction_candidates"].values())
        if int(item["candidate_id"]["orbit_index"]) == args.orbit
    )
assert best["candidate_id"] == {"q": 8, "old_fibre_degree": 2, "orbit_index": args.orbit}
record = next(item for item in neighbors["neighbors"] if int(item["orbit_index"]) == args.orbit)
frame_path = ROOT / neighbors["frame"]
frame = load_gram(frame_path)
assert frame.det() == 948 and frame.nrows() == 17
assert frame[:11, :11] == matrix(ZZ, frame[:11, :11])

witness = vector(ZZ, record["witness"])
fibre = vector(ZZ, record["fiber"])
ns = block_diagonal_matrix(U2, -frame)
old_fibre = vector(ZZ, [1, 0] + [0] * 17)
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
assert witness * frame * witness == 16
assert fibre == vector(ZZ, [4, 2] + list(witness))
assert fibre * ns * fibre == 0 and fibre * ns * old_fibre == 2 and fibre * ns * old_zero == 2
assert gcd(tuple(ns * fibre)) == 1

effective_simple = tuple(
    vector(ZZ, [0, 0] + [-ZZ(index == node) for index in range(17)])
    for node in range(11)
)
component_pairings = tuple(fibre * ns * item for item in effective_simple)
affine_curves = tuple(
    old_fibre + vector(ZZ, [0, 0] + list(highest) + [0] * 6)
    for highest in highest_roots(frame[:11, :11])
)
affine_pairings = tuple(fibre * ns * item for item in affine_curves)
assert component_pairings == tuple(record["dominant_labels"])
assert affine_pairings == tuple(
    map(ZZ, best["explicit_curve_degrees"]["physical_old_A11_fibre_components"][11:])
)
assert all(value >= 0 for value in component_pairings + affine_pairings)

section_distances = section_distance_profile(frame, witness)
print(f"A11O{args.orbit}_SECTION_PROFILE|{','.join(map(str, section_distances)) or 'empty'}", flush=True)
assert not section_distances or min(section_distances) >= 2
# Universal degree-two parity exclusion for a negative bisection.
assert (witness * frame * witness - 2) % 4 != 0

child, neighbor_basis, mate = neighbor_frame(ns, fibre)
roots, root_basis, root_data = roots_and_data(child)
assert child.det() == 948
assert tuple(root_data[:2]) == tuple(map(int, record["child_root_data"][:2]))
adapted, adapted_basis, child_height, child_cartan = root_adaptation(child)
child_ade = ade_name(child_cartan)
child_mw_rank = 17 - int(root_data[0])
transition = block_diagonal_matrix(identity_matrix(ZZ, 2), adapted_basis) * neighbor_basis
inverse = transition.inverse()
assert abs(transition.det()) == 1 and inverse.change_ring(ZZ) == inverse
inverse = inverse.change_ring(ZZ)
assert transition * ns * transition.transpose() == block_diagonal_matrix(U2, -adapted)
assert inverse * block_diagonal_matrix(U2, -adapted) * inverse.transpose() == ns

args.frame_output.parent.mkdir(parents=True, exist_ok=True)
args.frame_output.write_text(
    f"# equation-cost A11 q8 orbit{args.orbit} child: {child_ade}/MW{child_mw_rank}\n"
    + "\n".join(" ".join(map(str, item)) for item in adapted.rows())
    + "\n"
)

payload = {
    "schema": "elkies-k3.h3-a11-q8-orbit849-lattice-certificate.v1",
    "status": "PASS_EXACT_A11_Q8_EQUATION_COST_LATTICE_CERTIFICATE",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS] + [str(frame_path.relative_to(ROOT))],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS + (frame_path,)
        },
    },
    "selection": best,
    "parent": {"ade": "A11", "mw_rank": 6, "frame": str(frame_path.relative_to(ROOT))},
    "edge": {
        "q": 8,
        "factorization": [4, 2],
        "old_fibre_degree": 2,
        "orbit_index": args.orbit,
        "primitive_nef_isotropic_fibre": entries(fibre),
        "fibre_square": 0,
        "divisibility": 1,
        "old_zero_intersection": 2,
        "component_pairings": entries(vector(ZZ, component_pairings)),
        "affine_component_pairings": entries(vector(ZZ, affine_pairings)),
        "section_distance_profile": [str(value) for value in section_distances],
        "minimum_section_distance": None if not section_distances else str(min(section_distances)),
        "bisection_parity_exclusion": True,
        "nef": True,
    },
    "marked_U": {
        "fibre_in_parent": entries(fibre),
        "isotropic_mate_in_parent": entries(mate),
        "zero_in_parent": entries(mate - fibre),
        "gram": [[0, 1], [1, 0]],
    },
    "child": {
        "ade": child_ade,
        "mw_rank": child_mw_rank,
        "root_data": list(map(int, root_data)),
        "frame_output": str(args.frame_output.resolve().relative_to(ROOT)),
        "frame_sha256": hashlib.sha256(args.frame_output.read_bytes()).hexdigest(),
        "frame": rows(adapted),
        "mw_height": rational_rows(child_height),
        "root_cartan": rows(child_cartan),
    },
    "transport": {
        "parent_to_child_basis": rows(transition),
        "child_to_parent_basis": rows(inverse),
        "forward_determinant": int(transition.det()),
        "inverse_determinant": int(inverse.det()),
        "forward_inverse_exact": transition * inverse == identity_matrix(ZZ, 19),
    },
    "route_status": (
        "Fully certified lattice edge and equation-cost lateral A11 node. "
        "No continuation from this marked child to pinned R17 is yet certified, "
        "so this is not a promoted lifting target."
    ),
    "proof_boundary": (
        "Exact primitive/nef/marked-U/unimodular-transport/root certificate. "
        "The equation-cost values remain planning estimates. This artifact does "
        "not construct the q8 pencil or identify a pinned-R17 continuation."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A11O{}|q=8|old_degree=2|PO={}|RR={}|section_min={}|"
    "child={}/MW{}|root={}|det_fwd={}|det_inv={}|nef=1|status={}".format(
        args.orbit,
        best["horizontal"]["P_dot_O"],
        best["expected_RR_ambient"],
        "empty" if not section_distances else min(section_distances),
        child_ade,
        child_mw_rank,
        ",".join(map(str, root_data)),
        transition.det(),
        inverse.det(),
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
