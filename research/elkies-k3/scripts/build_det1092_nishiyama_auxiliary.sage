#!/usr/bin/env sage
"""Construct the determinant-1092 Nishiyama auxiliary for the X1092 census.

The output is only the auxiliary/gluing gate.  It proves that the recovered
rootless rank-17 frame can be primitively complemented inside an even
unimodular rank-24 lattice.  It does not enumerate its Niemeier embeddings or
classify the rootless J2 frames of X1092.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

from sage.all import IntegralLattice, QQ, ZZ, block_diagonal_matrix, identity_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
PARENT = ROOT / "artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_parent_v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/det1092_nishiyama_auxiliary_v1.json"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def primitive_rows(value):
    diagonal = value.smith_form()[0]
    return all(abs(diagonal[index, index]) == 1 for index in range(value.nrows()))


def build():
    parent = json.loads(PARENT.read_text())
    assert parent["status"] == "EXPLICIT_MODEL_AND_BASIS_INPUTS"
    frame = matrix(QQ, parent["generic_height_gram"])
    assert frame.nrows() == 17 and frame.is_positive_definite() and frame.det() == 1092

    # D5 with two integral vectors.  This fixed presentation was obtained by
    # solving the exact Schur-complement equation det(K)=1092, then checking
    # discriminant anti-isometry.  There is no random or bounded selection in
    # the certificate itself.
    auxiliary = matrix(ZZ, [
        [2, -1, 0, 0, 0, 1, -1],
        [-1, 2, -1, 0, 0, 1, 0],
        [0, -1, 2, -1, -1, 0, 1],
        [0, 0, -1, 2, 0, 0, 1],
        [0, 0, -1, 0, 2, -1, 0],
        [1, 1, 0, 0, -1, 20, -2],
        [-1, 0, 1, 1, 0, -2, 22],
    ])
    assert auxiliary.is_positive_definite() and auxiliary.det() == 1092

    auxiliary_group = IntegralLattice(auxiliary).discriminant_group()
    frame_group = IntegralLattice(frame).discriminant_group()
    assert tuple(auxiliary_group.invariants()) == (1092,)
    assert tuple(frame_group.invariants()) == (1092,)
    q_auxiliary = auxiliary_group.gram_matrix_quadratic()[0, 0]
    q_frame = frame_group.gram_matrix_quadratic()[0, 0]
    anti_units = [
        unit for unit in range(1092)
        if math.gcd(unit, 1092) == 1
        and ((q_auxiliary * unit * unit + q_frame) / 2).denominator() == 1
    ]
    assert anti_units

    split = block_diagonal_matrix(auxiliary, frame)
    auxiliary_generator = vector(QQ, auxiliary_group.gen(0).lift())
    frame_generator = vector(QQ, frame_group.gen(0).lift())
    unit = anti_units[0]
    ambient = IntegralLattice(split).overlattice([
        vector(QQ, list(unit * auxiliary_generator) + list(frame_generator))
    ])
    ambient_gram = ambient.gram_matrix()
    assert ambient.rank() == 24 and ambient_gram.det() == 1
    assert all(entry.denominator() == 1 for entry in ambient_gram.list())
    ambient_gram = matrix(ZZ, ambient_gram)
    assert all(ambient_gram[index, index] % 2 == 0 for index in range(24))

    split_auxiliary = identity_matrix(ZZ, 24)[:7]
    auxiliary_coordinates = split_auxiliary * ambient.basis_matrix().inverse()
    assert all(entry.denominator() == 1 for entry in auxiliary_coordinates.list())
    auxiliary_coordinates = matrix(ZZ, auxiliary_coordinates)
    assert primitive_rows(auxiliary_coordinates)
    assert auxiliary_coordinates * ambient_gram * auxiliary_coordinates.transpose() == auxiliary
    complement_coordinates = (auxiliary_coordinates * ambient_gram).right_kernel_matrix()
    complement_gram = complement_coordinates * ambient_gram * complement_coordinates.transpose()
    assert complement_gram.nrows() == 17 and complement_gram.det() == 1092
    assert primitive_rows(complement_coordinates)
    isometry = pari(frame).qfisom(pari(complement_gram))
    assert isometry != 0
    assert int(pari(complement_gram).qfminim(2)[0]) == 0

    return {
        "schema": "elkies-k3.det1092-nishiyama-auxiliary.v1",
        "status": "PASS_EXACT_AUXILIARY_AND_UNIMODULAR_GLUE",
        "inputs": {str(PARENT.relative_to(ROOT)): sha256(PARENT)},
        "auxiliary": {
            "gram": rows(auxiliary),
            "determinant": int(auxiliary.det()),
            "discriminant_invariants": list(map(int, auxiliary_group.invariants())),
            "discriminant_generator_q": str(q_auxiliary),
        },
        "frame": {
            "determinant": int(frame.det()),
            "discriminant_invariants": list(map(int, frame_group.invariants())),
            "discriminant_generator_q": str(q_frame),
            "minimum_root_count": int(pari(frame).qfminim(2)[0]),
        },
        "anti_isometry_units_mod_1092": anti_units,
        "chosen_anti_isometry_unit": unit,
        "unimodular_glue": {
            "rank": 24,
            "determinant": int(ambient_gram.det()),
            "even": True,
            "primitive_auxiliary": True,
            "primitive_complement": True,
            "complement_determinant": int(complement_gram.det()),
            "complement_isometric_to_recovered_frame": True,
            "complement_root_count": int(pari(complement_gram).qfminim(2)[0]),
        },
        "next_gate": "Enumerate every primitive embedding of this auxiliary in all 23 Niemeier lattices and deduplicate precisely the rootless rank17 complements by integral isometry.",
        "claim_boundary": "This produces one valid Nishiyama auxiliary and gluing. It is not a complete Niemeier embedding enumeration, X1092 J2 classification, marked-U realization, equation construction, or specialization search.",
    }


def canonical(value):
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    payload = canonical(build())
    if arguments.build:
        assert not arguments.output.exists()
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_bytes(payload)
        print(f"WROTE {arguments.output.relative_to(ROOT)}")
    assert arguments.output.read_bytes() == payload
    print("PASS exact determinant1092 Nishiyama auxiliary and unimodular glue")
