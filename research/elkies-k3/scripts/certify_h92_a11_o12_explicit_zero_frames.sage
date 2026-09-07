#!/usr/bin/env sage -python
"""Reframe equation-side orbit12 using its two explicit degree-one curves."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
CERT = GENERATED / "elkies-k3-h3-a11-q8-orbit12-lattice-certificate.json"
NEIGHBORS = LOCAL / "q24-a11-orbit64-q8-all.json"
OUTPUT = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frames.json"
FRAME_OUTPUT = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frame.txt"
INPUTS = (CERT, NEIGHBORS)
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(ZZ, [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
                       if line.strip() and not line.lstrip().startswith("#")])


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


def rational_rows(value):
    return [[str(item) for item in row] for row in value.rows()]


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    half = [vector(ZZ, column) for column in matrix(ZZ, result[2]).columns()]
    roots = tuple(half + [-item for item in half])
    basis = matrix(ZZ, [list(item) for item in roots]).row_module().basis_matrix()
    root_gram = basis * gram * basis.transpose()
    return roots, basis, (basis.rank(), count, abs(ZZ(root_gram.det())))


def simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    positive = [item for item in roots if next(value for value in item if value) > 0]
    positive_set = {tuple(item) for item in positive}
    simple = matrix(ZZ, [list(item) for item in positive
                         if not any(tuple(item - left) in positive_set for left in positive)])
    assert simple.nrows() == simple.rank() == data[0]
    return simple, simple * gram * simple.transpose()


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


def ade(cartan):
    names = []
    for component in components(cartan):
        block = cartan.matrix_from_rows_and_columns(component, component)
        rank = block.nrows()
        determinant = abs(ZZ(block.det()))
        count = ZZ(pari(block).qfminim(2)[0])
        if determinant == rank + 1 and count == rank * (rank + 1):
            names.append(f"A{rank}")
        elif rank >= 4 and determinant == 4 and count == 2 * rank * (rank - 1):
            names.append(f"D{rank}")
        else:
            names.append({(6, 3, 72): "E6", (7, 2, 126): "E7", (8, 1, 240): "E8"}[(rank, determinant, count)])
    return "+".join(names)


def adapt(child):
    unused, root_basis, data = roots_and_data(child)
    rank = data[0]
    smith, left, right = root_basis.smith_form()
    assert smith == left * root_basis * right
    assert tuple(abs(smith[index, index]) for index in range(rank)) == (1,) * rank
    simple, cartan = simple_roots(child)
    basis = simple.stack(right.inverse()[rank:])
    adapted = basis * child * basis.transpose()
    height = adapted[rank:, rank:] - adapted[rank:, :rank] * cartan.inverse() * adapted[:rank, rank:]
    scale = lcm(value.denominator() for value in height.list())
    lll = matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram())
    quotient = block_diagonal_matrix(identity_matrix(ZZ, rank), lll.transpose())
    basis = quotient * basis
    adapted = basis * child * basis.transpose()
    height = adapted[rank:, rank:] - adapted[rank:, :rank] * cartan.inverse() * adapted[:rank, rank:]
    assert abs(basis.det()) == 1
    return adapted, basis, height, cartan, data


certificate = json.loads(CERT.read_text())
neighbors = json.loads(NEIGHBORS.read_text())
assert certificate["status"] == "PASS_EXACT_A11_Q8_EQUATION_COST_LATTICE_CERTIFICATE"
record = next(item for item in neighbors["neighbors"] if int(item["orbit_index"]) == 12)
assert [index for index, value in enumerate(record["dominant_labels"]) if int(value) == 1] == [9]

parent_frame = load_matrix(ROOT / neighbors["frame"])
g_parent = block_diagonal_matrix(U2, -parent_frame)
parent_to_abstract = matrix(ZZ, certificate["transport"]["parent_to_child_basis"])
abstract = matrix(ZZ, certificate["child"]["frame"])
g_abstract = block_diagonal_matrix(U2, -abstract)
abstract_inverse = parent_to_abstract.inverse().change_ring(ZZ)
fibre = vector(ZZ, [1, 0] + [0] * 17)
curves = {
    "old_A11_component_9": vector(ZZ, [0, 0] + [-ZZ(index == 9) for index in range(17)]),
    "old_A11_affine": vector(ZZ, [1, 0] + [1] * 11 + [0] * 6),
}
assert all(curve * g_parent * curve == -2 and curve * g_parent * parent_to_abstract.row(0) == 1
           for curve in curves.values())

candidates = []
for name, parent_curve in curves.items():
    section = parent_curve * abstract_inverse
    mate = section + fibre
    assert section * g_abstract * section == -2
    assert section * g_abstract * fibre == 1 and mate * g_abstract * mate == 0
    complement = matrix(ZZ, [list(fibre * g_abstract), list(mate * g_abstract)]).right_kernel_matrix()
    split = matrix(ZZ, [list(fibre), list(mate)] + [list(item) for item in complement.rows()])
    raw = -(complement * g_abstract * complement.transpose())
    child, adapted_basis, height, cartan, root_data = adapt(raw)
    explicit_transition = block_diagonal_matrix(identity_matrix(ZZ, 2), adapted_basis) * split
    total = explicit_transition * parent_to_abstract
    inverse = total.inverse().change_ring(ZZ)
    g_child = block_diagonal_matrix(U2, -child)
    assert abs(total.det()) == 1
    assert total * g_parent * total.transpose() == g_child
    assert inverse * g_child * inverse.transpose() == g_parent
    candidates.append({
        "explicit_zero_curve": name,
        "section_in_abstract_child": [int(value) for value in section],
        "frame": child,
        "mw_height": height,
        "root_cartan": cartan,
        "root_data": root_data,
        "equation_A11_to_explicit_zero_basis": total,
        "explicit_zero_to_equation_A11_basis": inverse,
        "frame_max_abs": max(abs(value) for value in child.list()),
        "frame_L1": sum(abs(value) for value in child.list()),
    })

candidates.sort(key=lambda item: (item["frame_max_abs"], item["frame_L1"], item["explicit_zero_curve"]))
selected = candidates[0]
assert ade(selected["root_cartan"]) == "A5+A5"
assert tuple(selected["root_data"]) == (10, 60, 36)
FRAME_OUTPUT.write_text(
    f"# orbit12 A5+A5/MW7 with explicit zero {selected['explicit_zero_curve']}\n"
    + "\n".join(" ".join(map(str, row)) for row in selected["frame"].rows()) + "\n"
)


def serialize(item):
    return {
        "explicit_zero_curve": item["explicit_zero_curve"],
        "section_in_abstract_child": item["section_in_abstract_child"],
        "frame": rows(item["frame"]),
        "mw_height": rational_rows(item["mw_height"]),
        "root_cartan": rows(item["root_cartan"]),
        "root_data": list(map(int, item["root_data"])),
        "equation_A11_to_explicit_zero_basis": rows(item["equation_A11_to_explicit_zero_basis"]),
        "explicit_zero_to_equation_A11_basis": rows(item["explicit_zero_to_equation_A11_basis"]),
        "forward_determinant": int(item["equation_A11_to_explicit_zero_basis"].det()),
        "inverse_determinant": int(item["explicit_zero_to_equation_A11_basis"].det()),
        "frame_max_abs": int(item["frame_max_abs"]),
        "frame_L1": int(item["frame_L1"]),
    }


payload = {
    "schema": "elkies-k3.h3-a11-q8-orbit12-explicit-zero-frames.v1",
    "status": "PASS_EXACT_A11_Q8_ORBIT12_EXPLICIT_ZERO_FRAMES",
    "inputs": {"paths": [str(path.relative_to(ROOT)) for path in INPUTS], "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in INPUTS}},
    "candidate_count": len(candidates),
    "selection_rule": "minimum frame max-absolute coefficient, then frame L1",
    "selected_zero_curve": selected["explicit_zero_curve"],
    "selected_frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "selected_frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "selected": serialize(selected),
    "candidates": [serialize(item) for item in candidates],
    "proof_boundary": "Exact full NS reframings using both already-explicit degree-one curves; no downstream route is promoted here.",
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("A11O12ZERO|selected={}|frame_max={}|frame_L1={}|det_fwd={}|det_inv={}|status={}".format(
    selected["explicit_zero_curve"], selected["frame_max_abs"], selected["frame_L1"],
    selected["equation_A11_to_explicit_zero_basis"].det(), selected["explicit_zero_to_equation_A11_basis"].det(), payload["status"]), flush=True)
print(f"OUTPUT|{OUTPUT.resolve()}", flush=True)
