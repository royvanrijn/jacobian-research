#!/usr/bin/env sage -python
"""Export every equation-effective q4/o208 child zero as a marked 3A3 state."""

import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, identity_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
SOURCE_MARKING = GENERATED / "elkies-k3-h3-a5a5-physical-component-chamber-marking.json"
ROUTE = GENERATED / "elkies-k3-h3-a5a5-physical-q4o208-to-pinned-r17-certificate.json"
CURRENT = GENERATED / "elkies-k3-h3-current_3A3-marked-frame.json"
OUTPUT = GENERATED / "elkies-k3-h3-q4o208-child-effective-zero-markings.json"
OUTPUT_PREFIX = "elkies-k3-h3-q4o208-child-effective-zero-"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


source_marking = json.loads(SOURCE_MARKING.read_text())
route = json.loads(ROUTE.read_text())
current = json.loads(CURRENT.read_text())
assert source_marking["status"] == "PASS_EXACT_A5A5_PHYSICAL_COMPONENT_CHAMBER_MARKING"
assert route["status"] == "PASS_EXACT_PHYSICAL_Q4O208_3A3_TO_PINNED_R17"
assert current["status"] == "PASS_EXACT_CURRENT_SUFFIX_STAGE_MARKING"

source_frame_path = ROOT / source_marking["frame_output"]
current_frame_path = ROOT / current["frame_output"]
source_frame = load_matrix(source_frame_path)
current_frame = load_matrix(current_frame_path)
source_gram = block_diagonal_matrix(U2, -source_frame)
current_gram = block_diagonal_matrix(U2, -current_frame)
fibre = vector(ZZ, route["fibre"]["class_in_parent"])
curves = {
    name: vector(ZZ, value)
    for name, value in source_marking["equation_explicit_curves_in_child"].items()
}
degree_one = sorted(name for name, curve in curves.items() if curve * source_gram * fibre == 1)
assert degree_one == sorted(item["zero"] for item in route["effective_zero_candidates"])

target_vectors = {
    name: vector(ZZ, value)
    for name, value in source_marking["target_fibres_in_root_adapted_hub"].items()
}


def zero_presentation(name):
    section = curves[name]
    mate = fibre + section
    kernel = matrix(ZZ, [list(fibre * source_gram), list(mate * source_gram)]).right_kernel_matrix()
    raw_basis = matrix(ZZ, [list(fibre), list(mate)] + [list(row) for row in kernel.rows()])
    assert abs(raw_basis.det()) == 1
    raw_frame = -(kernel * source_gram * kernel.transpose())
    qiso = matrix(ZZ, pari(raw_frame).qfisom(pari(current_frame)))
    assert qiso and qiso.transpose() * current_frame * qiso == raw_frame
    current_tail_in_raw = qiso.transpose().inverse().change_ring(ZZ)
    current_in_raw = block_diagonal_matrix(identity_matrix(ZZ, 2), current_tail_in_raw)
    basis_in_source = current_in_raw * raw_basis
    source_in_basis = basis_in_source.inverse().change_ring(ZZ)
    assert basis_in_source * source_gram * basis_in_source.transpose() == current_gram
    assert basis_in_source.row(0) == fibre
    assert basis_in_source.row(1) == mate

    transformed_curves = {
        curve_name: entries(curve * source_in_basis)
        for curve_name, curve in curves.items()
    }
    transformed_targets = {
        target_name: entries(target * source_in_basis)
        for target_name, target in target_vectors.items()
    }
    assert all(
        vector(ZZ, value) * current_gram * vector(ZZ, value) == -2
        for value in transformed_curves.values()
    )
    assert vector(ZZ, transformed_curves[name]) == vector(ZZ, [-1, 1] + [0] * 17)
    for target_name, target in target_vectors.items():
        transformed = vector(ZZ, transformed_targets[target_name])
        assert transformed * current_gram * transformed == 0
        assert transformed * current_gram * vector(ZZ, [1, 0] + [0] * 17) == (
            target * source_gram * fibre
        )

    slug = name.lower()
    path = GENERATED / f"{OUTPUT_PREFIX}{slug}-marking.json"
    payload = {
        "schema": "elkies-k3.h3-q4o208-child-effective-zero-marking.v1",
        "status": "PASS_EXACT_Q4O208_CHILD_EFFECTIVE_ZERO_MARKING",
        "hub": f"q4o208_child_3A3_zero_{name}",
        "zero": name,
        "root_data": [9, 36, 64],
        "frame_output": str(current_frame_path.relative_to(ROOT)),
        "target_fibres_in_root_adapted_hub": transformed_targets,
        "equation_explicit_curves_in_child": transformed_curves,
        "child_basis_in_physical_2A5_parent": rows(basis_in_source),
        "physical_2A5_parent_basis_in_child": rows(source_in_basis),
        "prefix_operational_score": route["compiler_profile"]["operational_equation_cost_score"],
        "proof_boundary": (
            "Exact re-zeroing of the already certified physical q4/orbit208 fibre on an "
            "already-explicit degree-one curve.  The marked U and both full NS bases are "
            "integral unimodular; all target fibres and explicit curves are transported "
            "exactly.  This is a lattice/equation-availability search state, not a compiled "
            "new equation pointing for zeros other than C5."
        ),
        "inputs": {
            "paths": [str(path.relative_to(ROOT)) for path in (SOURCE_MARKING, ROUTE, CURRENT)],
            "sha256": {
                str(path.relative_to(ROOT)): sha256(path)
                for path in (SOURCE_MARKING, ROUTE, CURRENT)
            },
        },
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return {
        "zero": name,
        "marking": str(path.relative_to(ROOT)),
        "basis_determinant": int(basis_in_source.det()),
        "explicit_degree_zero_count": sum(
            vector(ZZ, value)[1] == 0 for value in transformed_curves.values()
        ),
        "explicit_section_count": sum(
            vector(ZZ, value)[1] == 1 for value in transformed_curves.values()
        ),
    }


presentations = [zero_presentation(name) for name in degree_one]
payload = {
    "schema": "elkies-k3.h3-q4o208-child-effective-zero-markings.v1",
    "status": "PASS_EXACT_Q4O208_ALL_EFFECTIVE_ZERO_MARKINGS",
    "presentations": presentations,
    "proof_boundary": (
        "Exhaustive over the four already-explicit physical curves of degree one over the "
        "certified q4/orbit208 fibre.  Each presentation uses the same root-adapted 3A3 "
        "positive frame but a distinct marked U and transported explicit-curve set."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"Q4O208ZEROS|count={len(presentations)}|status={payload['status']}|output={OUTPUT}")
