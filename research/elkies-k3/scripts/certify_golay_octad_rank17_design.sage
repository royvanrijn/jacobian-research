#!/usr/bin/env sage
"""Certify a rank-17 lattice designed from seven Golay octads.

The proposal layer is the support-intersection geometry of the extended
binary Golay code.  Every accepted assertion after the octads are loaded is
an exact integral-lattice computation: Construction A, primitive closure,
orthogonal complement, short shell, discriminant module, and the ternary
K3-realizability genus gate.

status: EXACT_COMPUTATION_NOT_STATUS_PROMOTED
claim: the pinned septuple defines a primitive determinant-720 rank-seven
  auxiliary in N(24A1), whose saturated rank-17 complement is rootless, has
  minimum 4 and 3064 norm-four vectors, and admits the displayed signature
  (2,1) discriminant-form mate.
inputs: elkies-k3/data/lattice-foundry/golay-octad-rank17-det720-v1.json
outputs: artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json
"""

import argparse
import json
from collections import Counter
from itertools import combinations, product
from pathlib import Path

from sage.all import GF, Genus, QQ, ZZ, QuadraticForm, identity_matrix, matrix, pari, vector
from sage.coding.golay_code import GolayCode
from sage.quadratic_forms.genera.genus import genera


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_INPUT = (
    ROOT
    / "elkies-k3/data/lattice-foundry/golay-octad-rank17-det720-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json"
)


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def rows_as_strings(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def discriminant_form_key(genus):
    normal = genus.discriminant_form().normal_form()
    return {
        "invariants": list(map(int, normal.invariants())),
        "quadratic_gram": rows_as_strings(normal.gram_matrix_quadratic()),
        "value_module": str(normal.value_module_qf()),
    }


def construction_a_model():
    code = GolayCode(GF(2), extended=True)
    generator = code.generator_matrix().echelon_form()
    assert code.length() == 24 and code.dimension() == 12
    assert code.minimum_distance() == 8
    assert generator[:, :12] == identity_matrix(GF(2), 12)

    # Rows are physical coordinates relative to orthogonal roots e_i^2=2.
    physical_basis = matrix(QQ, 24, 24)
    for index in range(12):
        physical_basis[index] = vector(
            QQ, [QQ(int(entry)) / 2 for entry in generator[index]]
        )
    for index in range(12):
        physical_basis[12 + index, 12 + index] = 1

    gram = physical_basis * (2 * identity_matrix(QQ, 24)) * physical_basis.transpose()
    assert gram in matrix(ZZ, 24, 24).parent()
    gram = matrix(ZZ, gram)
    assert gram.det() == 1 and gram.is_positive_definite()
    root_shell = pari(gram).qfminim(2)
    assert int(root_shell[0]) == 48
    return code, physical_basis, gram


def octad_mask(support):
    assert len(support) == 8 and len(set(support)) == 8
    assert all(1 <= entry <= 24 for entry in support)
    return sum(1 << (entry - 1) for entry in support)


def octad_ambient_coordinates(mask, inverse_physical_basis):
    physical = vector(
        QQ, [QQ(1) / 2 if (mask >> index) & 1 else 0 for index in range(24)]
    )
    coordinates = physical * inverse_physical_basis
    assert all(entry in ZZ for entry in coordinates)
    return vector(ZZ, coordinates)


def primitive_closure(basis):
    original = matrix(ZZ, basis)
    saturated = original.row_module(ZZ).saturation().basis_matrix()
    coordinates = original * saturated.pseudoinverse()
    assert all(entry.denominator() == 1 for entry in coordinates.list())
    index = abs(int(matrix(ZZ, coordinates).det()))
    return saturated, index


def short_shell(gram, bound):
    result = pari(gram).qfminim(bound)
    representatives = [
        vector(ZZ, column)
        for column in matrix(ZZ, result[2].sage()).columns()
    ]
    assert 2 * len(representatives) == int(result[0])
    by_norm = Counter(int(item * gram * item) for item in representatives)
    return {
        "minimum_squared_norm": min(by_norm),
        "signed_counts": {
            str(norm): 2 * count for norm, count in sorted(by_norm.items())
        },
    }


def combinatorial_norm_four_count(code, auxiliary_masks):
    # Type I: two coordinate roots with arbitrary signs.
    coordinate_pair_vectors = 0
    for left, right in combinations(range(24), 2):
        for left_sign, right_sign in product((-1, 1), repeat=2):
            if all(
                left_sign * ((mask >> left) & 1)
                + right_sign * ((mask >> right) & 1)
                == 0
                for mask in auxiliary_masks
            ):
                coordinate_pair_vectors += 1

    # Type II: every choice of signs on a Golay octad. These are exactly the
    # norm-four vectors in the nontrivial glue cosets of N(24A1).
    signed_octad_vectors = 0
    octad_count = 0
    for word in code:
        if word.hamming_weight() != 8:
            continue
        octad_count += 1
        support = [index for index, entry in enumerate(word) if entry]
        for signs in product((-1, 1), repeat=8):
            if all(
                sum(
                    sign
                    for coordinate, sign in zip(support, signs)
                    if (mask >> coordinate) & 1
                )
                == 0
                for mask in auxiliary_masks
            ):
                signed_octad_vectors += 1
    assert octad_count == 759
    return {
        "signed_coordinate_pair_vectors": coordinate_pair_vectors,
        "signed_octad_glue_vectors": signed_octad_vectors,
        "signed_total": coordinate_pair_vectors + signed_octad_vectors,
    }


def certify(payload):
    code, physical_basis, ambient_gram = construction_a_model()
    supports = payload["octads"]
    assert len(supports) == 7
    masks = [octad_mask(support) for support in supports]
    assert len(set(masks)) == 7

    for mask in masks:
        word = vector(GF(2), [(mask >> index) & 1 for index in range(24)])
        assert word in code and word.hamming_weight() == 8

    coverage = [
        sum((mask >> coordinate) & 1 for mask in masks)
        for coordinate in range(24)
    ]
    assert min(coverage) >= 1

    raw_intersection_gram = matrix(
        ZZ,
        7,
        7,
        lambda row, column: (
            4
            if row == column
            else (masks[row] & masks[column]).bit_count() // 2
        ),
    )
    assert raw_intersection_gram.det() == 720

    inverse_physical_basis = physical_basis.inverse()
    raw_auxiliary_basis = matrix(
        ZZ,
        [
            octad_ambient_coordinates(mask, inverse_physical_basis)
            for mask in masks
        ],
    )
    assert (
        raw_auxiliary_basis
        * ambient_gram
        * raw_auxiliary_basis.transpose()
        == raw_intersection_gram
    )

    auxiliary_basis, closure_index = primitive_closure(raw_auxiliary_basis)
    assert closure_index == 1
    auxiliary_gram = auxiliary_basis * ambient_gram * auxiliary_basis.transpose()
    assert auxiliary_gram.det() == 720

    complement_basis = (auxiliary_basis * ambient_gram).right_kernel_matrix()
    assert complement_basis.nrows() == 17
    frame_gram = complement_basis * ambient_gram * complement_basis.transpose()
    assert frame_gram.det() == auxiliary_gram.det() == 720
    assert frame_gram.is_positive_definite()

    shell = short_shell(frame_gram, 4)
    assert shell["minimum_squared_norm"] == 4
    assert shell["signed_counts"].get("2", 0) == 0
    assert shell["signed_counts"]["4"] == 3064
    combinatorial_shell = combinatorial_norm_four_count(code, masks)
    assert combinatorial_shell == {
        "signed_coordinate_pair_vectors": 24,
        "signed_octad_glue_vectors": 3040,
        "signed_total": 3064,
    }

    smith = frame_gram.smith_form()[0]
    smith_invariants = [
        abs(int(smith[index, index]))
        for index in range(frame_gram.nrows())
        if abs(int(smith[index, index])) > 1
    ]
    assert smith_invariants == [2, 6, 60]

    target_form = discriminant_form_key(Genus(frame_gram))
    ternary_genera = genera((2, 1), frame_gram.det(), even=True)
    matches = [
        genus
        for genus in ternary_genera
        if discriminant_form_key(genus) == target_form
    ]
    assert len(ternary_genera) == 20 and len(matches) == 1
    ternary_gram = matrix(ZZ, matches[0].representative())
    assert ternary_gram.det() == -720
    assert QuadraticForm(QQ, ternary_gram).signature() == 1

    return {
        "schema": "elkies-k3.golay-octad-rank17-design.v1",
        "status": "PASS_EXACT_GOLAY_OCTAD_RANK17_DET720_LATTICE_DESIGN",
        "proof_scope": {
            "proved": (
                "The seven pinned supports are Golay octads. Their glue "
                "vectors primitively span a determinant-720 rank-seven "
                "auxiliary in N(24A1). The saturated orthogonal complement "
                "has rank 17, minimum 4, and 3064 norm-four vectors. Its "
                "discriminant form has an exact even signature-(2,1) mate, "
                "so U plus the negative frame primitively embeds in the K3 "
                "lattice."
            ),
            "not_proved": (
                "The bounded proposal search was not exhaustive and no "
                "determinant or theta-shell optimality is claimed. No "
                "rational K3 family, arithmetic marking, Weierstrass "
                "equation, section basis over QQ(t), or specialization-rank "
                "consequence is constructed."
            ),
        },
        "input": payload,
        "golay_code": {
            "parameters": [24, 12, 8],
            "octad_count": 759,
            "construction_a_ambient": "N(24A1)",
            "ambient_rank": 24,
            "ambient_determinant": int(ambient_gram.det()),
            "ambient_signed_roots": 48,
        },
        "support_design": {
            "hex_masks": [format(mask, "06x") for mask in masks],
            "coordinate_coverage_multiplicities": coverage,
            "raw_octad_intersection_gram": rows(raw_intersection_gram),
            "raw_determinant": int(raw_intersection_gram.det()),
        },
        "auxiliary": {
            "rank": 7,
            "primitive_closure_index": closure_index,
            "determinant": int(auxiliary_gram.det()),
            "ambient_basis": rows(auxiliary_basis),
            "gram": rows(auxiliary_gram),
        },
        "frame": {
            "rank": 17,
            "determinant": int(frame_gram.det()),
            "minimum_squared_norm": shell["minimum_squared_norm"],
            "signed_short_shell_through_norm_4": shell["signed_counts"],
            "independent_combinatorial_norm_four_count": combinatorial_shell,
            "norm_four_vectors": shell["signed_counts"]["4"],
            "norm_four_unoriented_pairs": shell["signed_counts"]["4"] // 2,
            "smith_invariants_greater_than_one": smith_invariants,
            "ambient_basis": rows(complement_basis),
            "gram": rows(frame_gram),
            "discriminant_form_normal_key": target_form,
        },
        "k3_realizability": {
            "candidate_NS": "U + frame(-1)",
            "NS_signature": [1, 18],
            "NS_determinant": -720,
            "matching_even_ternary_genera": len(matches),
            "all_even_ternary_genera_at_determinant": len(ternary_genera),
            "transcendental_signature": [2, 1],
            "transcendental_gram": rows(ternary_gram),
            "transcendental_determinant": int(ternary_gram.det()),
            "discriminant_form_relation": "q_T = q_frame = -q_NS",
        },
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

input_payload = json.loads(arguments.input.read_text())
result = certify(input_payload)
rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"

if arguments.check:
    assert arguments.output.read_text() == rendered
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(rendered)

print(
    "GOLAYOCTAD17|det=720|min=4|norm4=3064|pairs=1532|"
    "smith=2,6,60|ternary_matches=1|status=PASS_EXACT"
)
