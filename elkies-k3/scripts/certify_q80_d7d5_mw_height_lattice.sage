#!/usr/bin/env sage
"""Certify the saturated MW5 height lattice at the generic Q80 D7+D5 node.

The existing chamber verifier constructs the exact second-q4 frame, its
primitive D7+D5 root lattice, and the pinned norm-eight third-q12 horizontal.
This wrapper reuses those objects and records the full saturated rank-five
orthogonal projection lattice.  It is an equation-search fingerprint: modular
section generators must recover this height lattice before a height-eight
shell can be called complete.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


WRAPPER = Path(__file__).resolve()
ROOT = WRAPPER.parents[2]
SOURCE = ROOT / "elkies-k3/scripts/analyze_q80_second_neighbor_chamber.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/q80-d7d5-mw5-height-lattice.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


# The source verifier derives the marked frame and target from the pinned Q80
# route.  Set its __file__ explicitly so its repository-relative inputs resolve
# identically whether this wrapper is run directly or through Sage load().
globals()["__file__"] = str(SOURCE)
load(str(SOURCE))
globals()["__file__"] = str(WRAPPER)

assert third_simple.nrows() == 12
assert third_root_gram.det() == 16
smith_diagonal = third_simple.smith_form()[0]
assert all(smith_diagonal[index, index] == 1 for index in range(12))

projection = (
    identity_matrix(QQ, 17)
    - second_frame
    * third_simple.transpose()
    * third_root_gram.inverse()
    * third_simple
)
projection_denominator = lcm(value.denominator() for value in projection.list())
scaled_projection = (projection_denominator * projection).change_ring(ZZ)
projected_integer_basis = scaled_projection.row_module().basis_matrix()
assert projected_integer_basis.nrows() == projected_integer_basis.rank() == 5

mw_basis = projected_integer_basis / projection_denominator
height_gram = mw_basis * second_frame * mw_basis.transpose()
assert height_gram.det() == QQ(237) / 4
height_scale = lcm(value.denominator() for value in height_gram.list())
lll_change = (height_scale * height_gram).change_ring(ZZ).LLL_gram().transpose()
assert abs(lll_change.det()) == 1
mw_basis = lll_change * mw_basis
height_gram = mw_basis * second_frame * mw_basis.transpose()
expected_height_gram = matrix(
    QQ,
    [
        [1, QQ(1) / 4, QQ(1) / 4, 0, QQ(1) / 4],
        [QQ(1) / 4, QQ(5) / 4, 0, 0, QQ(-1) / 2],
        [QQ(1) / 4, 0, QQ(7) / 4, 0, QQ(1) / 4],
        [0, 0, 0, 4, -2],
        [QQ(1) / 4, QQ(-1) / 2, QQ(1) / 4, -2, QQ(35) / 4],
    ],
)
assert height_gram == expected_height_gram
assert height_gram.det() == QQ(237) / 4

target_coordinates = (
    third_mw_projection
    * second_frame
    * mw_basis.transpose()
    * height_gram.inverse()
)
assert all(value in ZZ for value in target_coordinates)
target_coordinates = vector(ZZ, target_coordinates)
assert tuple(target_coordinates) == (-1, 1, -1, 1, 0)
target_height = (target_coordinates * height_gram * target_coordinates.column())[0]
assert target_height == 8
assert third_mw_projection == target_coordinates * mw_basis
assert third_section_pole == 2
assert all(value == 0 for _, value in third_section_component_pairings[:12])
assert tuple(value for _, value in third_section_component_pairings[12:]) == (1, 1)

# The fixed-u modular worker recovers the following rank-three polynomial
# subgroup.  Enumerate every integral embedding into MW5.  The largest target
# norm is twelve in the 4-scaled form, so PARI qfminim(12) is a complete finite
# vector list for all three basis images.
polynomial_height_gram = matrix(
    QQ,
    [
        [QQ(11) / 4, QQ(-1) / 2, QQ(-5) / 4],
        [QQ(-1) / 2, 3, QQ(-1) / 2],
        [QQ(-5) / 4, QQ(-1) / 2, 1],
    ],
)
scaled_height_gram = (4 * height_gram).change_ring(ZZ)
scaled_polynomial_gram = (4 * polynomial_height_gram).change_ring(ZZ)
short_half = matrix(ZZ, pari(scaled_height_gram).qfminim(12)[2]).transpose()
short_vectors = []
for row in short_half.rows():
    short_vectors.extend((vector(ZZ, row), -vector(ZZ, row)))
vectors_by_norm = {
    norm: [row for row in short_vectors if row * scaled_height_gram * row == norm]
    for norm in (4, 11, 12)
}
embeddings = []
for first in vectors_by_norm[11]:
    for second in vectors_by_norm[12]:
        for third in vectors_by_norm[4]:
            candidate = matrix(ZZ, [first, second, third])
            if candidate * scaled_height_gram * candidate.transpose() == scaled_polynomial_gram:
                embeddings.append(candidate)
assert len(embeddings) == 2 and embeddings[1] == -embeddings[0]
assert tuple(tuple(row) for row in embeddings[0].rows()) == (
    (1, 1, 0, 0, 0),
    (1, -1, -1, 0, 0),
    (-1, 0, 0, 0, 0),
)

polynomial_shell_coordinates = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (2, 1, 4),
    (1, 1, 2),
    (-1, -1, -1),
    (-1, 0, -1),
    (-1, 0, -2),
)
embedding_records = []
for embedding in embeddings:
    pairings = tuple(
        target_coordinates * height_gram * vector(ZZ, row)
        for row in embedding.rows()
    )
    assert pairings[:2] == (0, 0) and abs(pairings[2]) == 1
    orientation = int(pairings[2])
    ordered_po = [
        4 - orientation * coordinates[2]
        for coordinates in polynomial_shell_coordinates
    ]
    embedding_records.append(
        {
            "basis_images_in_MW5_height_coordinates": [
                list(map(int, row)) for row in embedding.rows()
            ],
            "H_pairings_with_polynomial_basis": list(map(str, pairings)),
            "ordered_P_dot_O_of_H_minus_eight_polynomial_points": ordered_po,
            "sorted_intersection_fingerprint": sorted(ordered_po),
        }
    )
assert {
    tuple(record["sorted_intersection_fingerprint"])
    for record in embedding_records
} == {
    (0, 2, 3, 4, 4, 5, 5, 6),
    (2, 3, 3, 4, 4, 5, 6, 8),
}

payload = {
    "schema": "elkies-k3.q80-d7d5-mw5-height-lattice.v1",
    "status": "PASS_EXACT_Q80_D7D5_SATURATED_MW5_AND_Q12_TARGET",
    "frame": {
        "roots": "D7+D5",
        "root_rank": 12,
        "root_determinant": 16,
        "root_lattice_primitive": True,
        "MW_rank": 5,
        "MW_torsion_order": 1,
        "NS_determinant": int(second_frame.det()),
        "projection_denominator": int(projection_denominator),
        "height_basis_in_frame_coordinates": rational_rows(mw_basis),
        "height_gram": rational_rows(height_gram),
        "height_determinant": str(height_gram.det()),
        "determinant_replay": "948/16=237/4",
    },
    "third_q12_horizontal": {
        "height_basis_coordinates": list(map(int, target_coordinates)),
        "frame_projection_coordinates": list(map(str, third_mw_projection)),
        "canonical_height": str(target_height),
        "P_dot_O": int(third_section_pole),
        "simple_root_component_pairings": [
            [name, int(value)] for name, value in third_section_component_pairings[:12]
        ],
        "affine_identity_component_pairings": [
            [name, int(value)] for name, value in third_section_component_pairings[12:]
        ],
        "component_profile": "identity at both I3* and I1* fibres",
    },
    "polynomial_rank3_subgroup": {
        "height_gram": rational_rows(polynomial_height_gram),
        "height_determinant": str(polynomial_height_gram.det()),
        "eight_points_in_basis_coordinates": [
            list(row) for row in polynomial_shell_coordinates
        ],
        "complete_embedding_count": len(embedding_records),
        "embeddings": embedding_records,
        "completeness_method": (
            "Enumerate every MW5 vector of 4-scaled norm at most 12 with "
            "PARI qfminim, then test the full 4-scaled rank-three Gram."
        ),
    },
    "inputs": {
        "source_verifier": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": sha256(SOURCE),
        },
        "wrapper": {
            "path": str(WRAPPER.relative_to(ROOT)),
            "sha256": sha256(WRAPPER),
        },
    },
    "claim_boundary": (
        "This is an exact lattice and marking certificate. It supplies the "
        "rank-five height-lattice and q12-horizontal fingerprints required by "
        "the modular equation search, but it does not construct a section or "
        "third-q12 equation over a finite field or over QQ."
    ),
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/certify_q80_d7d5_mw_height_lattice.sage"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
output = args.output.resolve()
if args.check:
    if not output.exists() or output.read_text() != serialized:
        raise SystemExit("stale Q80 D7+D5 MW5 height-lattice certificate")
else:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(serialized)

print(
    "Q80D7D5MW|rank=5|det=237/4|target=(-1,1,-1,1,0)|"
    "target_height=8|target_PO=2|status="
    "PASS_EXACT_Q80_D7D5_SATURATED_MW5_AND_Q12_TARGET",
    flush=True,
)
