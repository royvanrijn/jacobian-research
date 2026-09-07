#!/usr/bin/env sage -python
"""Reframe a certified child using named old-A11 components as explicit zeros."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
U2 = matrix(ZZ, ((0, 1), (1, 0)))
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--certificate", type=Path, required=True)
parser.add_argument("--a11-component", type=int, action="append", required=True)
parser.add_argument("--expected-ade", required=True)
parser.add_argument("--expected-root-data", required=True, help="rank,count,determinant")
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--frame-output", type=Path, required=True)
args = parser.parse_args()
CERT = args.certificate.resolve()
OUTPUT = args.output.resolve()
FRAME_OUTPUT = args.frame_output.resolve()


def rows(value):
    return [[int(item) for item in row] for row in value.rows()]


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
    result = []
    exceptional = {(6, 3, 72): "E6", (7, 2, 126): "E7", (8, 1, 240): "E8"}
    for component in components(cartan):
        block = cartan.matrix_from_rows_and_columns(component, component)
        rank = block.nrows()
        determinant = abs(ZZ(block.det()))
        count = ZZ(pari(block).qfminim(2)[0])
        if determinant == rank + 1 and count == rank * (rank + 1):
            result.append(f"A{rank}")
        elif rank >= 4 and determinant == 4 and count == 2 * rank * (rank - 1):
            result.append(f"D{rank}")
        else:
            result.append(exceptional[(rank, determinant, count)])
    return "+".join(result)


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
    basis = block_diagonal_matrix(identity_matrix(ZZ, rank), lll.transpose()) * basis
    adapted = basis * child * basis.transpose()
    height = adapted[rank:, rank:] - adapted[rank:, :rank] * cartan.inverse() * adapted[:rank, rank:]
    assert abs(basis.det()) == 1
    return adapted, basis, height, cartan, data


certificate = json.loads(CERT.read_text())
child = (
    matrix(ZZ, certificate["child"]["frame"])
    if "frame" in certificate["child"]
    else matrix(
        ZZ,
        [[ZZ(value) for value in line.split()]
         for line in (ROOT / certificate["frame_output"]).read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )
)
g_child = block_diagonal_matrix(U2, -child)
equation_to_child = matrix(
    ZZ,
    certificate.get("transport", certificate)["equation_A11_to_child_basis"],
)
child_to_equation = equation_to_child.inverse().change_ring(ZZ)
fibre = vector(ZZ, [1, 0] + [0] * 17)
expected_root_data = tuple(map(ZZ, args.expected_root_data.split(",")))

candidates = []
for node in args.a11_component:
    curve_equation = vector(ZZ, [0, 0] + [-ZZ(index == node) for index in range(17)])
    section = curve_equation * child_to_equation
    assert section * g_child * section == -2 and section * g_child * fibre == 1
    mate = section + fibre
    complement = matrix(ZZ, [list(fibre * g_child), list(mate * g_child)]).right_kernel_matrix()
    split = matrix(ZZ, [list(fibre), list(mate)] + [list(item) for item in complement.rows()])
    raw = -(complement * g_child * complement.transpose())
    reframed, adapted_basis, height, cartan, root_data = adapt(raw)
    reframing = block_diagonal_matrix(identity_matrix(ZZ, 2), adapted_basis) * split
    total = reframing * equation_to_child
    inverse = total.inverse().change_ring(ZZ)
    g_reframed = block_diagonal_matrix(U2, -reframed)
    g_equation = child_to_equation * g_child * child_to_equation.transpose()
    assert abs(total.det()) == 1
    # Direct identities in the certified child and equation markings.
    assert reframing * g_child * reframing.transpose() == g_reframed
    assert total * g_equation * total.transpose() == g_reframed
    assert inverse * g_reframed * inverse.transpose() == g_equation
    actual_ade = ade(cartan)
    assert sorted(actual_ade.split("+")) == sorted(args.expected_ade.split("+")) and tuple(root_data) == expected_root_data
    candidates.append({
        "a11_component_index": node,
        "section_in_child": [int(value) for value in section],
        "frame": reframed,
        "mw_height": height,
        "root_cartan": cartan,
        "root_data": root_data,
        "child_to_explicit_zero_basis": reframing,
        "equation_A11_to_explicit_zero_basis": total,
        "explicit_zero_to_equation_A11_basis": inverse,
        "frame_max_abs": max(abs(value) for value in reframed.list()),
        "frame_L1": sum(abs(value) for value in reframed.list()),
    })

candidates.sort(key=lambda item: (item["frame_max_abs"], item["frame_L1"], item["a11_component_index"]))
selected = candidates[0]
FRAME_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
FRAME_OUTPUT.write_text(f"# {args.expected_ade} explicit zero old_A11_component_{selected['a11_component_index']}\n" +
                        "\n".join(" ".join(map(str, row)) for row in selected["frame"].rows()) + "\n")


def serialize(item):
    return {
        "a11_component_index": item["a11_component_index"],
        "section_in_child": item["section_in_child"],
        "frame": rows(item["frame"]),
        "mw_height": [[str(value) for value in row] for row in item["mw_height"].rows()],
        "root_cartan": rows(item["root_cartan"]),
        "root_data": list(map(int, item["root_data"])),
        "child_to_explicit_zero_basis": rows(item["child_to_explicit_zero_basis"]),
        "equation_A11_to_explicit_zero_basis": rows(item["equation_A11_to_explicit_zero_basis"]),
        "explicit_zero_to_equation_A11_basis": rows(item["explicit_zero_to_equation_A11_basis"]),
        "forward_determinant": int(item["equation_A11_to_explicit_zero_basis"].det()),
        "inverse_determinant": int(item["explicit_zero_to_equation_A11_basis"].det()),
        "frame_max_abs": int(item["frame_max_abs"]),
        "frame_L1": int(item["frame_L1"]),
    }


payload = {
    "schema": "elkies-k3.h3-candidate-a11-component-explicit-zero-frames.v1",
    "status": "PASS_EXACT_CANDIDATE_A11_COMPONENT_EXPLICIT_ZERO_FRAMES",
    "input_certificate": str(CERT.relative_to(ROOT)),
    "input_sha256": hashlib.sha256(CERT.read_bytes()).hexdigest(),
    "candidate_count": len(candidates),
    "selection_rule": "minimum frame max-absolute coefficient, then frame L1, then A11 component index",
    "selected_component_index": selected["a11_component_index"],
    "selected_frame_output": str(FRAME_OUTPUT.relative_to(ROOT)),
    "selected_frame_sha256": hashlib.sha256(FRAME_OUTPUT.read_bytes()).hexdigest(),
    "selected": serialize(selected),
    "candidates": [serialize(item) for item in candidates],
    "proof_boundary": "Exact full NS reframings by the supplied already-explicit old-A11 component curves; downstream route certification is separate.",
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("CANDIDATEZERO|selected_component={}|frame_max={}|frame_L1={}|det_fwd={}|det_inv={}|status={}".format(
    selected["a11_component_index"], selected["frame_max_abs"], selected["frame_L1"],
    selected["equation_A11_to_explicit_zero_basis"].det(), selected["explicit_zero_to_equation_A11_basis"].det(), payload["status"]), flush=True)
print(f"OUTPUT|{OUTPUT}", flush=True)
