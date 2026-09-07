#!/usr/bin/env sage -python
"""Re-mark the orbit1991 child using its two explicit degree-one components.

The generic neighbour splitter supplies an abstract zero.  For equation work,
either of the two old A11 fibre components meeting the new fibre once is a
strictly better marked zero because its curve is already explicit.  This
checker constructs and root-adapts both full integral U splittings, then
selects one by a deterministic frame-coefficient convention.
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
    default=GENERATED / "elkies-k3-h3-a11-q8-orbit1991-explicit-zero-frames.json",
)
parser.add_argument(
    "--frame-output",
    type=Path,
    default=GENERATED / "elkies-k3-h3-a11-q8-orbit1991-explicit-zero-frame.txt",
)
args = parser.parse_args()

CERT = GENERATED / "elkies-k3-h3-a11-q8-orbit1991-lattice-certificate.json"
NEIGHBORS = LOCAL / "q24-a11-orbit64-q8-all.json"
INPUTS = (CERT, NEIGHBORS)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")


def entries(value):
    return [int(item) for item in vector(ZZ, value)]


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


def rational_rows(value):
    return [[str(item) for item in row] for row in value.rows()]


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    half = [vector(ZZ, column) for column in matrix(ZZ, result[2]).columns()]
    roots = tuple(half + [-item for item in half])
    root_basis = matrix(ZZ, [list(item) for item in roots]).row_module().basis_matrix()
    root_gram = root_basis * gram * root_basis.transpose()
    return roots, root_basis, (root_basis.rank(), count, abs(ZZ(root_gram.det())))


def deterministic_simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    positive = [item for item in roots if next(value for value in item if value != 0) > 0]
    positive_set = {tuple(item) for item in positive}
    simple = matrix(
        ZZ,
        [
            list(item)
            for item in positive
            if not any(tuple(item - left) in positive_set for left in positive)
        ],
    )
    assert simple.nrows() == simple.rank() == data[0]
    return simple, simple * gram * simple.transpose()


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


def component_name(cartan, component):
    block = cartan.matrix_from_rows_and_columns(component, component)
    rank = block.nrows()
    determinant = abs(ZZ(block.det()))
    count = ZZ(pari(block).qfminim(2)[0])
    if determinant == rank + 1 and count == rank * (rank + 1):
        return f"A{rank}"
    if rank >= 4 and determinant == 4 and count == 2 * rank * (rank - 1):
        return f"D{rank}"
    return {(6, 3, 72): "E6", (7, 2, 126): "E7", (8, 1, 240): "E8"}.get(
        (rank, determinant, count), f"R{rank}d{determinant}n{count}"
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
    return adapted, adapted_basis, height, cartan, invariants


certificate = json.loads(CERT.read_text())
neighbors = json.loads(NEIGHBORS.read_text())
assert certificate["status"] == "PASS_EXACT_A11_Q8_EQUATION_COST_LATTICE_CERTIFICATE"
record = next(item for item in neighbors["neighbors"] if int(item["orbit_index"]) == 1991)
degree_one_nodes = [index for index, value in enumerate(record["dominant_labels"]) if int(value) == 1]
assert degree_one_nodes == [8, 10]

parent_to_abstract_child = matrix(ZZ, certificate["transport"]["parent_to_child_basis"])
parent_frame_path = ROOT / neighbors["frame"]
parent_frame = matrix(
    ZZ,
    [
        [ZZ(value) for value in line.split()]
        for line in parent_frame_path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ],
)
g_parent = block_diagonal_matrix(U2, -parent_frame)
abstract_child = matrix(ZZ, certificate["child"]["frame"])
g_abstract = block_diagonal_matrix(U2, -abstract_child)
abstract_inverse = parent_to_abstract_child.inverse().change_ring(ZZ)
fibre = vector(ZZ, [1, 0] + [0] * 17)

candidates = []
for node in degree_one_nodes:
    parent_curve = vector(ZZ, [0, 0] + [-ZZ(index == node) for index in range(17)])
    section = parent_curve * abstract_inverse
    assert section * g_abstract * section == -2 and section * g_abstract * fibre == 1
    mate = section + fibre
    assert mate * g_abstract * mate == 0 and mate * g_abstract * fibre == 1
    complement = matrix(
        ZZ, [list(fibre * g_abstract), list(mate * g_abstract)]
    ).right_kernel_matrix()
    split = matrix(ZZ, [list(fibre), list(mate)] + [list(item) for item in complement.rows()])
    assert abs(split.det()) == 1
    raw = -(complement * g_abstract * complement.transpose())
    adapted, adapted_basis, height, cartan, root_data = root_adaptation(raw)
    explicit_transition = block_diagonal_matrix(identity_matrix(ZZ, 2), adapted_basis) * split
    total = explicit_transition * parent_to_abstract_child
    inverse = total.inverse().change_ring(ZZ)
    assert abs(total.det()) == 1
    # The two exact identities below are the operative transport checks.
    g_explicit = block_diagonal_matrix(U2, -adapted)
    assert explicit_transition * g_abstract * explicit_transition.transpose() == g_explicit
    assert total * g_parent * total.transpose() == g_explicit
    assert inverse * g_explicit * inverse.transpose() == g_parent
    assert explicit_transition.inverse().change_ring(ZZ) * g_explicit * explicit_transition.inverse().change_ring(ZZ).transpose() == g_abstract
    candidates.append(
        {
            "old_A11_component_index": node,
            "section_in_abstract_child": entries(section),
            "frame": adapted,
            "height": height,
            "cartan": cartan,
            "root_data": root_data,
            "abstract_child_to_explicit_zero_basis": explicit_transition,
            "equation_A11_to_explicit_zero_basis": total,
            "explicit_zero_to_equation_A11_basis": inverse,
            "frame_max_abs": max(abs(item) for item in adapted.list()),
            "frame_L1": sum(abs(item) for item in adapted.list()),
        }
    )

candidates.sort(
    key=lambda item: (
        item["frame_max_abs"],
        item["frame_L1"],
        item["old_A11_component_index"],
    )
)
selected = candidates[0]
print(f"A11O1991ZERO_ADE|{ade_name(selected['cartan'])}", flush=True)
assert ade_name(selected["cartan"]) == "D10+A1+A2"
assert tuple(selected["root_data"]) == (13, 188, 24)

args.frame_output.parent.mkdir(parents=True, exist_ok=True)
args.frame_output.write_text(
    "# orbit1991 A1+A2+D10/MW4 frame with explicit old-A11 component as zero\n"
    + "\n".join(" ".join(map(str, item)) for item in selected["frame"].rows())
    + "\n"
)

def serialized(item):
    return {
        "old_A11_component_index": item["old_A11_component_index"],
        "section_in_abstract_child": item["section_in_abstract_child"],
        "frame": rows(item["frame"]),
        "mw_height": rational_rows(item["height"]),
        "root_cartan": rows(item["cartan"]),
        "root_data": list(map(int, item["root_data"])),
        "abstract_child_to_explicit_zero_basis": rows(item["abstract_child_to_explicit_zero_basis"]),
        "equation_A11_to_explicit_zero_basis": rows(item["equation_A11_to_explicit_zero_basis"]),
        "explicit_zero_to_equation_A11_basis": rows(item["explicit_zero_to_equation_A11_basis"]),
        "forward_determinant": int(item["equation_A11_to_explicit_zero_basis"].det()),
        "inverse_determinant": int(item["explicit_zero_to_equation_A11_basis"].det()),
        "frame_max_abs": int(item["frame_max_abs"]),
        "frame_L1": int(item["frame_L1"]),
    }

payload = {
    "schema": "elkies-k3.h3-a11-q8-orbit1991-explicit-zero-frames.v1",
    "status": "PASS_EXACT_A11_Q8_ORBIT1991_EXPLICIT_ZERO_FRAMES",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "candidate_count": len(candidates),
    "selection_rule": "minimum frame max-absolute coefficient, then frame L1, then component index",
    "selected_component_index": selected["old_A11_component_index"],
    "selected_frame_output": str(args.frame_output.resolve().relative_to(ROOT)),
    "selected_frame_sha256": hashlib.sha256(args.frame_output.read_bytes()).hexdigest(),
    "selected": serialized(selected),
    "candidates": [serialized(item) for item in candidates],
    "proof_boundary": (
        "Exact full-rank NS changes of basis for both already-explicit degree-one "
        "old A11 components. This selects a compiler marking but does not execute "
        "the equation lift or certify any downstream route to pinned R17."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A11O1991ZERO|candidates=2|selected_component={}|child=A1+A2+D10/MW4|"
    "frame_max={}|det_fwd={}|det_inv={}|status={}".format(
        selected["old_A11_component_index"],
        selected["frame_max_abs"],
        selected["equation_A11_to_explicit_zero_basis"].det(),
        selected["explicit_zero_to_equation_A11_basis"].det(),
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
