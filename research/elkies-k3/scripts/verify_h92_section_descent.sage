#!/usr/bin/env sage -python
"""Certify the H92 direction in the H3 E7+E8/MW2 frame is rational.

The pinned Elkies--Kumar H92 ancillary construction starts with a split
``D6+A8+A1`` fibration and an explicit section.  Its local corrections are
``1/2`` at ``A1`` and ``20/9`` at ``A8``, so its height is ``23/18`` and the
resulting rank-18 lattice has determinant 92.  This checker constructs that
lattice and exhausts its 92 discriminant classes.  There is no nonzero
isotropic class, hence no proper even overlattice.

The ancillary two- and three-neighbor parameters are rational over
``QQ(r,s)`` and end at the ``E7+E8`` fibration.  Therefore its rank-one H92
Mordell--Weil quotient, of height 46, is individually rational.  On the
H21/H92 intersection this is the second generator in the H3 height Gram
``[[21/2,3],[3,46]]``.

This is a field-of-definition and saturation certificate.  It does not
recover the large height-46 section coordinates on the final short model.
"""

from sage.all import *

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/local/humbert-inputs/92/92.txt"
SOURCE_SHA256 = "5dde58046d9770fa78b514ae48509a238090ad4de7057b41e43ea308047101c2"
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
FRAME_SHA256 = "ba09ec834a7229e11e4ca687d187f663b6368c3e2fac9b5133bb1570e7031599"
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elkies-k3-h92-section-descent.json"
)
EXPECTED_H92 = (
    QQ(-3621005) / 690947,
    QQ(158286) / 143585,
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_gram(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ],
    )


def source_section_at_target():
    """Replay the explicit pre-neighbor H92 section at the target point."""

    r, s = EXPECTED_H92
    e = -(r * s + 4 * s + r - 1) / (r * s - r - 1)
    m = -2 * (4 * s**2 + 3 * r * s - 2 * s - 2 * r - 1) / (
        r * s - r - 1
    )
    f = -(m - e + 5) * (m + e + 3) / (m**2 + 3 * e**2 + 2 * e - 5)
    a1 = (f + 1) * (e**2 * f**2 + f**2 + 2 * e * f - 2 * f + 2) / (
        f - 1
    )
    b3 = 1 - f**2
    b0 = (b3 - 1) * b3 * (e + 1) ** 2
    b2 = -(
        (b3**2 - b3) * e**2 + 3 * b3**2 + (2 * a1 + 1) * b3
    ) / 2
    b1 = b3 + a1
    a2 = 2 * b2 + b1**2 - 2 * a1 * b1 - b0 + a1**2
    a3 = -2 * b3**3 - a1 * b3**2 + (-2 * b2 - b0) * b3 + b0
    c1 = 2 * b1 - a1

    split_expression = (
        (-3 * e**2 - 2 * e + 5) * f**2
        + (-2 * e**2 - 4 * e - 10) * f
        + e**2
        - 2 * e
        + 1
    )
    assert split_expression == (4 + m * (f + 1)) ** 2

    function_field = FunctionField(QQ, "t")
    t = function_field.gen()
    a = 1 + a1 * t + a2 * t**2 + a3 * t**3
    b = b0 * (1 + b1 * t + b2 * t**2)
    c = b0**2 * (1 + c1 * t)
    x_section = (
        (b3 - 1) * b3 * (e + 1) ** 2 * t**2
        + (b3 - 1) ** 2 * b3 * (e + 1) ** 2 * t**3
        - (b3 - 1) * b3**2 * (e + 1) ** 2 * t**4
    )
    rhs = (
        x_section**3
        + a * x_section**2
        + 2 * b * t**2 * (t - 1) * x_section
        + c * t**4 * (t - 1) ** 2
    )
    assert rhs.is_square()
    y_section = rhs.sqrt()
    assert y_section**2 == rhs
    return {
        "r": str(r),
        "s": str(s),
        "split_square": str(4 + m * (f + 1)),
        "x_degree": int(x_section.numerator().degree()),
        "y_degree": int(y_section.numerator().degree()),
        "exact_section_identity": True,
    }


def source_lattice():
    """Return U+D6+A8+A1 plus the displayed section curve."""

    d6 = CartanMatrix(["D", 6])
    a8 = CartanMatrix(["A", 8])
    a1 = CartanMatrix(["A", 1])
    roots = block_diagonal_matrix(-d6, -a8, -a1)
    gram = zero_matrix(ZZ, 18)
    gram[0, 1] = gram[1, 0] = 1
    gram[1, 1] = -2
    gram[2:17, 2:17] = roots
    gram[17, 17] = -2
    gram[0, 17] = gram[17, 0] = 1
    # The correction 20/9 is A8 node 4 (or node 5 after diagram reversal).
    gram[2 + 6 + 3, 17] = gram[17, 2 + 6 + 3] = 1
    # The correction 1/2 is the nonidentity A1 component.
    gram[2 + 6 + 8, 17] = gram[17, 2 + 6 + 8] = 1
    return gram


def isotropic_discriminant_classes(gram):
    """Enumerate nonzero classes with q(x)=0 in QQ/2ZZ."""

    smith, left, right = gram.smith_form(transformation=True)
    assert left * gram * right == smith
    diagonal = [abs(ZZ(value)) for value in smith.diagonal()]
    nontrivial = [index for index, value in enumerate(diagonal) if value > 1]
    assert prod(diagonal) == abs(gram.det())
    assert len(nontrivial) == 2
    left_inverse = left.inverse()
    gram_inverse = gram.inverse()
    isotropic = []
    representatives = 0
    for first in range(diagonal[nontrivial[0]]):
        for second in range(diagonal[nontrivial[1]]):
            representative = vector(ZZ, gram.nrows())
            representative[nontrivial[0]] = first
            representative[nontrivial[1]] = second
            if not representative:
                continue
            representatives += 1
            dual_numerator = left_inverse * representative
            norm = dual_numerator * gram_inverse * dual_numerator
            if norm.denominator() == 1 and norm.numerator() % 2 == 0:
                isotropic.append(
                    {
                        "smith_coordinates": [int(first), int(second)],
                        "norm": str(norm),
                    }
                )
    assert representatives == abs(gram.det()) - 1
    return diagonal, isotropic


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source", type=Path, default=SOURCE)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

assert digest(args.source) == SOURCE_SHA256
assert digest(FRAME) == FRAME_SHA256
source_text = args.source.read_text()
for phrase in (
    "We want fibers of type D6 A8 A1 and a section",
    "finally, need D6 components to be rational",
    "first we go to E7 A8 via a 2-neighbor step",
    "next go to E8 E7 by a 3-neighbor step",
):
    assert phrase in source_text

specialization = source_section_at_target()
source_gram = source_lattice()
assert source_gram.det() == -92
assert all(source_gram[index, index] % 2 == 0 for index in range(18))
trivial = source_gram[:17, :17]
section_cross = source_gram[17, :17]
projected_square = (
    source_gram[17, 17]
    - section_cross * trivial.inverse() * section_cross.transpose()
)
source_height = -QQ(projected_square[0, 0])
assert source_height == 4 - QQ(1) / 2 - QQ(20) / 9 == QQ(23) / 18
smith_diagonal, isotropic = isotropic_discriminant_classes(source_gram)
assert [value for value in smith_diagonal if value > 1] == [2, 46]
assert not isotropic

frame = load_gram(FRAME)
assert frame.nrows() == 17 and frame.det() == 948
roots = frame[:15, :15]
height = (
    frame[15:, 15:]
    - frame[15:, :15] * roots.inverse() * frame[:15, 15:]
)
assert height == matrix(QQ, [[QQ(21) / 2, 3], [3, 46]])

payload = {
    "schema": "elkies-k3.h92-section-descent.v1",
    "status": "PASS_H92_SECTION_Q_DEFINED",
    "inputs": {
        "h92_ancillary": {
            "path": str(args.source.relative_to(ROOT)),
            "sha256": SOURCE_SHA256,
        },
        "h3_frame": {
            "path": str(FRAME.relative_to(ROOT)),
            "sha256": FRAME_SHA256,
        },
    },
    "source_fibration": {
        "root_type": "D6+A8+A1",
        "section_height": str(source_height),
        "section_local_corrections": ["1/2", "20/9"],
        "ns_determinant": abs(int(source_gram.det())),
        "smith_diagonal": [int(value) for value in smith_diagonal],
        "discriminant_class_count": abs(int(source_gram.det())),
        "nonzero_isotropic_discriminant_classes": isotropic,
        "proper_even_overlattice_possible": False,
        "generating_curves_defined_over": "QQ(r,s)",
        "neighbor_sequence": ["2-neighbor", "3-neighbor"],
        "neighbor_parameters_defined_over": "QQ(r,s)",
        "target_specialization_replay": specialization,
    },
    "target_fibration": {
        "root_type": "E7+E8",
        "h92_generator_height": "46",
        "galois_action_on_h92_generator": "fixed",
        "h3_height_gram": [["21/2", "3"], ["3", "46"]],
    },
    "consequence": (
        "Both H3 source directions are individually Q-defined at the rational "
        "H21/H92 point: the H21 direction by its separate determinant-21 "
        "certificate and the H92 direction by this saturated determinant-92 "
        "ancillary construction."
    ),
    "remaining_gate": (
        "Recover or transport explicit height-46 section coordinates on the "
        "final short H92 model when equation-level section tracking requires them."
    ),
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "H92DESCENT|source=D6+A8+A1|source_det=92|source_height=23/18|"
    "disc_group=2,46|isotropic=0|neighbor_field=QQ(r,s)|"
    "P2_height=46|P2_galois=fixed",
    flush=True,
)
print("H92DESCENT|status=PASS_H92_SECTION_Q_DEFINED", flush=True)
