#!/usr/bin/env sage-python
"""Certify the integral involution-glue calculus in the two E6 examples.

The rational Mordell--Weil height pairing is multiplied by 12 throughout, so
all displayed lattices are even and integral.  The script exhausts every
isotropic graph between the relevant 2-torsion eigendiscriminant subgroups,
constructs the corresponding overlattice, and classifies the results by
exact integral isometry.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

from sage.all import GF, QQ, ZZ, block_diagonal_matrix, identity_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-integral-character-glue-calculus-v1.json"
)


def digest(relative_path):
    return hashlib.sha256((ROOT / relative_path).read_bytes()).hexdigest()


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def lattice_intrinsics(gram):
    minimum_data = pari(gram).qfminim(gram.nrows() * 32)
    columns = matrix(ZZ, minimum_data[2].sage()).columns()
    norms = [int(vector(ZZ, column) * gram * vector(ZZ, column)) for column in columns]
    minimum = min(norms)
    minimum_pairs = sum(norm == minimum for norm in norms)
    smith = [
        abs(int(value))
        for value in gram.smith_form()[0].diagonal()
        if abs(int(value)) > 1
    ]
    return {
        "rank": gram.nrows(),
        "determinant_absolute": abs(int(gram.det())),
        "smith_invariants": smith,
        "minimum": minimum,
        "minimum_vectors_signed": 2 * minimum_pairs,
        "root_count_signed": 2 * sum(norm == 2 for norm in norms),
        "automorphism_group_order": int(pari(gram).qfauto()[0]),
    }


def qfisometric(left, right):
    return pari(left).qfisom(pari(right)) != 0


def nonzero_binary_vectors(rank):
    return [
        vector(ZZ, entries)
        for entries in itertools.product((0, 1), repeat=rank)
        if any(entries)
    ]


def gl2_f2():
    result = []
    for entries in itertools.product((0, 1), repeat=4):
        value = matrix(GF(2), 2, 2, entries)
        if value.is_invertible():
            result.append(matrix(ZZ, value))
    assert len(result) == 6
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    e6_21_path = (
        "artifacts/generated-results/"
        "elkies-k3-e6-ii-rank3-quadratic-base-change-v1.json"
    )
    e6_22_path = (
        "artifacts/generated-results/"
        "elkies-k3-e6-rank4-linear-chord-incidence-v1.json"
    )
    e6_21 = json.loads((ROOT / e6_21_path).read_text())
    e6_22 = json.loads((ROOT / e6_22_path).read_text())
    assert e6_21["status"] == "PASS_EXACT_E6_II_RANK_SUM_3_RHO19_ROOTLESS_IMPOSSIBLE"
    assert e6_22["status"] == "PASS_EXACT_E6_RANK4_INCIDENCE_DESCENT"

    plus = matrix(ZZ, [[16, -8], [-8, 16]])
    minus_21 = matrix(ZZ, [[8]])
    minus_22 = matrix(ZZ, [[88, 16], [16, 88]])
    pure_21 = block_diagonal_matrix(plus, minus_21)
    pure_22 = block_diagonal_matrix(plus, minus_22)

    # All nonzero half-classes used below are isotropic in the scaled even
    # discriminant forms.
    plus_classes = nonzero_binary_vectors(2)
    for value in plus_classes:
        assert (value * plus * value) % 8 == 0
    assert minus_21[0, 0] % 8 == 0
    for value in nonzero_binary_vectors(2):
        assert (value * minus_22 * value) % 8 == 0

    index_two_grams = []
    index_two_records = []
    for value in plus_classes:
        companion = next(
            candidate
            for candidate in plus_classes
            if abs(matrix(ZZ, [value, candidate]).det()) == 1
        )
        half_sum = vector(QQ, [value[0] / 2, value[1] / 2, QQ(1) / 2])
        change = matrix(
            QQ,
            [
                [companion[0], companion[1], 0],
                [0, 0, 1],
                list(half_sum),
            ],
        )
        assert abs(change.det()) == QQ(1) / 2
        gram = change * pure_21 * change.transpose()
        assert gram in matrix(ZZ, 3, 3).parent()
        gram = gram.change_ring(ZZ)
        assert all(gram[index, index] % 2 == 0 for index in range(3))
        index_two_grams.append(gram)
        index_two_records.append(
            {
                "plus_line_generator_mod_2": list(map(int, value)),
                "gram": rows(gram),
                "intrinsics": lattice_intrinsics(gram),
            }
        )
    assert all(qfisometric(index_two_grams[0], gram) for gram in index_two_grams)

    index_four_grams = []
    index_four_records = []
    for graph in gl2_f2():
        change = matrix(
            QQ,
            [
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [QQ(1) / 2, 0, QQ(graph[0, 0]) / 2, QQ(graph[0, 1]) / 2],
                [0, QQ(1) / 2, QQ(graph[1, 0]) / 2, QQ(graph[1, 1]) / 2],
            ],
        )
        assert abs(change.det()) == QQ(1) / 4
        gram = change * pure_22 * change.transpose()
        assert gram in matrix(ZZ, 4, 4).parent()
        gram = gram.change_ring(ZZ)
        assert all(gram[index, index] % 2 == 0 for index in range(4))
        index_four_grams.append(gram)
        index_four_records.append(
            {
                "graph_matrix_over_F2": rows(graph),
                "gram": rows(gram),
                "intrinsics": lattice_intrinsics(gram),
            }
        )
    assert all(qfisometric(index_four_grams[0], gram) for gram in index_four_grams)

    actual_scaled_22 = matrix(
        ZZ,
        [
            [16, -8, 8, -4],
            [-8, 16, -4, 8],
            [8, -4, 26, 2],
            [-4, 8, 2, 26],
        ],
    )
    assert qfisometric(index_four_grams[0], actual_scaled_22)
    assert abs(pure_21.det()) == 1536
    assert abs(index_two_grams[0].det()) == 384
    assert abs(pure_22.det()) == 1437696
    assert abs(index_four_grams[0].det()) == 89856

    payload = {
        "schema": "elkies-k3.integral-character-glue-calculus.v1",
        "status": "PASS_EXACT_INVOLUTION_GRAPH_GLUE_CLASSIFICATION",
        "theorem": {
            "statement": (
                "For an even lattice L with involution, L/(L+ plus L-) is "
                "killed by two and embeds in both eigendiscriminant groups as "
                "an isotropic anti-isometry graph. Conversely every such graph "
                "defines an even involution-stable overlattice."
            ),
            "root_warning": (
                "Graph glue only enlarges the pure eigensum, so it cannot remove "
                "a norm-two vector already present there."
            ),
        },
        "inputs": {e6_21_path: digest(e6_21_path), e6_22_path: digest(e6_22_path)},
        "E6_2_plus_1": {
            "pure_character_lattice": {
                "gram": rows(pure_21),
                "intrinsics": lattice_intrinsics(pure_21),
            },
            "actual_glue_index": 1,
            "alternative_nonzero_graph_count": 3,
            "alternative_integral_isometry_class_count": 1,
            "alternative_index_two_graphs": index_two_records,
            "conclusion": (
                "There are exactly two graph-glue isometry types: the actual "
                "index-one pure sum and one index-two type. The latter has "
                "smaller minimum and is a negative control, not a root killer."
            ),
        },
        "E6_2_plus_2": {
            "pure_character_lattice": {
                "gram": rows(pure_22),
                "intrinsics": lattice_intrinsics(pure_22),
            },
            "required_glue_index": 4,
            "full_graph_count": 6,
            "integral_isometry_class_count": 1,
            "full_graphs": index_four_records,
            "actual_scaled_saturated_gram": rows(actual_scaled_22),
            "conclusion": (
                "The six full graph identifications form one integral isometry "
                "class. Hence the observed two half-sums exhaust the index-four "
                "character-glue possibilities up to integral isometry."
            ),
        },
        "proof_boundary": {
            "proved": (
                "Complete finite graph enumeration and exact integral-isometry "
                "classification after scaling the rational height pairing by 12."
            ),
            "not_proved": (
                "The scaled lattice classification does not independently prove "
                "geometric section descent; that is supplied by the two input "
                "equation certificates."
            ),
        },
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "certify_integral_character_glue_calculus.sage --check"
        ),
    }

    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    if arguments.check:
        if not output.exists() or output.read_text() != serialized:
            raise SystemExit(f"stale or missing artifact: {output}")
        print("PASS integral character glue calculus")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(serialized)
    try:
        print(output.relative_to(ROOT))
    except ValueError:
        print(output)


if __name__ == "__main__":
    main()
