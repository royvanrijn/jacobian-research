#!/usr/bin/env sage
"""Build the exact Leech/Co0 ambient backend from the Atlas representation.

The GAP AtlasRep package supplies the characteristic-zero 24-dimensional
matrix representation of 2.Co1 = Co0.  Its invariant symmetric form is
one-dimensional.  Solving that invariant-form system, primitively scaling the
positive generator, and checking determinant/minimum identifies the Leech
lattice in the representation basis.

This pins the separate Leech ambient and its Co0 generators.  It does not
enumerate primitive rank-seven sublattices or claim any determinant-band
completeness.

status: EXACT_BACKEND_FOUNDATION_NOT_EMBEDDING_ENUMERATION
output: artifacts/generated-results/elkies-k3-leech-co0-backend-v1.json
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

from sage.all import QQ, ZZ, gcd, lcm, libgap, matrix, pari, vector


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-leech-co0-backend-v1.json"


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def matrix_digest(value):
    encoded = "\n".join(" ".join(map(str, row)) for row in value.rows()) + "\n"
    return hashlib.sha256(encoded.encode()).hexdigest()


def symmetric_coordinates(rank):
    return [(row, column) for row in range(rank) for column in range(row, rank)]


def invariant_form(generators):
    rank = generators[0].nrows()
    coordinates = symmetric_coordinates(rank)
    equations = []
    # AtlasRep matrices act on column vectors.  For every upper-triangular
    # entry, expand A^t G A-G linearly in the symmetric coordinates of G.
    for generator in generators:
        for target_row, target_column in coordinates:
            equation = []
            for source_row, source_column in coordinates:
                coefficient = (
                    generator[source_row, target_row]
                    * generator[source_column, target_column]
                )
                if source_row != source_column:
                    coefficient += (
                        generator[source_column, target_row]
                        * generator[source_row, target_column]
                    )
                if (source_row, source_column) == (target_row, target_column):
                    coefficient -= 1
                equation.append(coefficient)
            equations.append(equation)
    kernel = matrix(QQ, equations).right_kernel_matrix()
    assert kernel.nrows() == 1
    invariant = kernel[0]
    denominator = lcm([entry.denominator() for entry in invariant])
    integral = vector(ZZ, [entry * denominator for entry in invariant])
    content = gcd([abs(entry) for entry in integral if entry])
    integral = integral / content
    gram = matrix(ZZ, rank, rank)
    for entry, (row, column) in zip(integral, coordinates):
        gram[row, column] = entry
        gram[column, row] = entry
    if not gram.is_positive_definite():
        gram = -gram
    assert gram.is_positive_definite()
    return gram, len(equations), len(coordinates)


def build():
    assert libgap.LoadPackage("atlasrep") is not libgap.fail
    group = libgap.AtlasGroup(
        "2.Co1", libgap.IsMatrixGroup, libgap.Dimension, 24, libgap.Characteristic, 0
    )
    generators = [
        matrix(ZZ, generator.sage())
        for generator in group.GeneratorsOfGroup()
    ]
    assert len(generators) == 2
    assert all(generator.nrows() == generator.ncols() == 24 for generator in generators)

    gram, equation_count, variable_count = invariant_form(generators)
    assert gram.det() == 1
    assert all(entry % 2 == 0 for entry in gram.diagonal())
    assert all(
        generator.transpose() * gram * generator == gram
        for generator in generators
    )

    shell = pari(gram).qfminim(4)
    signed_minimal_vectors = int(shell[0])
    representatives = matrix(ZZ, shell[2].sage()).ncols()
    assert signed_minimal_vectors == 196560
    assert 2 * representatives == signed_minimal_vectors
    group_order = int(group.Size())
    assert group_order == 8315553613086720000

    return {
        "schema": "elkies-k3.leech-co0-backend.v1",
        "status": "PASS_EXACT_LEECH_GRAM_AND_CO0_ATLAS_ACTION_BACKEND_FOUNDATION",
        "proof_scope": {
            "proved": (
                "The AtlasRep 24-dimensional characteristic-zero representation "
                "of 2.Co1 has a one-dimensional invariant symmetric-form space. "
                "Its primitive positive integral generator is even, unimodular, "
                "positive definite, has minimum 4 and 196560 minimal vectors; "
                "hence it is the Leech lattice. Both displayed generators preserve "
                "the Gram matrix and generate Co0 of the certified order."
            ),
            "not_proved": (
                "No primitive rank-seven sublattice, Co0 embedding orbit, ternary "
                "transcendental mate, or determinant-band completeness is enumerated."
            ),
        },
        "atlas_representation": {
            "group_name": "2.Co1",
            "identified_group": "Co0",
            "dimension": 24,
            "characteristic": 0,
            "generator_count": len(generators),
            "group_order": group_order,
            "generator_sha256": [matrix_digest(item) for item in generators],
            "generators": [rows(item) for item in generators],
        },
        "invariant_form_calculation": {
            "symmetric_variables": variable_count,
            "linear_equations": equation_count,
            "kernel_dimension": 1,
            "action_convention": "A^t G A = G on column vectors",
        },
        "leech_lattice": {
            "rank": 24,
            "gram": rows(gram),
            "gram_sha256": matrix_digest(gram),
            "determinant": int(gram.det()),
            "even": True,
            "minimum_squared_norm": 4,
            "signed_minimal_vectors": signed_minimal_vectors,
            "unoriented_minimal_pairs": representatives,
            "roots": 0,
        },
        "enumeration_contract": {
            "orbit_group": "Co0",
            "weyl_group": "trivial",
            "primitive_auxiliary_rank": 7,
            "complement_rank": 17,
            "complements_automatically_rootless": True,
            "determinant_bound": 5000,
            "state": "AMBIENT_READY_EMBEDDING_ORBITS_OPEN",
        },
        "reproduction": {
            "command": (
                "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
                "elkies-k3/scripts/build_leech_co0_backend.sage"
            ),
            "check_command": (
                "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
                "elkies-k3/scripts/build_leech_co0_backend.sage --check"
            ),
        },
    }


parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
payload = build()
encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"

if arguments.check:
    if not arguments.output.exists() or arguments.output.read_text() != encoded:
        raise SystemExit("Leech/Co0 backend artifact is stale")
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(encoded)

print(
    "LEECHCO0|det=1|min=4|minvecs=196560|group_order={}|"
    "invariant_kernel=1|status=PASS".format(
        payload["atlas_representation"]["group_order"]
    )
)
