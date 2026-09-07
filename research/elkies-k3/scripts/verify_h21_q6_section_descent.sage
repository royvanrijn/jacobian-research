#!/usr/bin/env sage -python
"""Certify that the H21 height-21/2 section is defined over Q(r,s).

The Elkies--Kumar H21 ancillary construction starts from a split
A2+A6+E8 fibration.  Its fiber components, the old zero section, and the
displayed 3-neighbor parameter are all defined over ``QQ(r,s)``.  The
resulting trivial lattice has squarefree determinant 21, so it has no proper
integral saturation.  Consequently the new E7+E8 Mordell--Weil generator is
an individually rational divisor class, rather than an unordered conjugate
pair.

This settles the Galois action on the desired q=6 divisor class.  It does not
yet construct the rational function spanning its two-dimensional RR pencil.
"""

from sage.all import *

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/local/humbert-inputs/21/21.txt"
SOURCE_SHA256 = "e0d0ea5ae18502fc0b51cf1999ab4e8b5755a40bffe2b7418ace6891d40a71a6"
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
FRAME_SHA256 = "ba09ec834a7229e11e4ca687d187f663b6368c3e2fac9b5133bb1570e7031599"
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elkies-k3-h21-q6-section-descent.json"
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


def coordinates(value):
    return [int(entry) for entry in value]


parser = argparse.ArgumentParser()
parser.add_argument("--source", type=Path, default=SOURCE)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

assert digest(args.source) == SOURCE_SHA256
assert digest(FRAME) == FRAME_SHA256
source_text = args.source.read_text()
for pinned_phrase in (
    "We'll begin with A2 A6 E8.",
    "We put E8 at infinity, A6 at t=0, and A2 at t=1",
    "regular on 1/3 and  2/7 components",
    "Applying the Jacobian(cubic) formulas gives us E8,E7",
):
    assert pinned_phrase in source_text

# Replay the endpoint square conditions used to split the A6 and A2 fibers.
R = PolynomialRing(QQ, names=("r", "s", "t"))
r, s, t = R.gens()
e = (r - s) ** 2 * (r**2 - 1)
c = e**2 * (1 - t) + e**2 * r**2 * t
b = e * (1 - t) + e * r * s * t
assert c(t=0) == e**2 and c(t=1) == (e * r) ** 2
assert b(t=0) == e and b(t=1) == e * r * s
assert b(t=0) ** 2 / c(t=0) == 1
assert b(t=1) ** 2 / c(t=1) == s**2

# U+A2+A6+E8 is generated over QQ(r,s).  Its determinant is squarefree, so
# no proper integral overlattice can contain it with finite index.
source_root_determinants = {
    "A2": ZZ(CartanMatrix(["A", 2]).det()),
    "A6": ZZ(CartanMatrix(["A", 6]).det()),
    "E8": ZZ(CartanMatrix(["E", 8]).det()),
}
source_ns_determinant = prod(source_root_determinants.values())
assert source_root_determinants == {"A2": 3, "A6": 7, "E8": 1}
assert source_ns_determinant == 21 and ZZ(21).is_squarefree()
assert all(21 % (index**2) for index in range(2, 1 + ZZ(21).isqrt()))

# In the specialized H3 frame, coordinate 15 is the H21 generator.  Project
# it away from E7+E8 and recover its exact height 21/2.
frame = load_gram(FRAME)
assert frame.nrows() == 17 and frame.det() == 948
roots = frame[:15, :15]
height = (
    frame[15:, 15:]
    - frame[15:, :15] * roots.inverse() * frame[:15, 15:]
)
assert height == matrix(QQ, [[QQ(21) / 2, 3], [3, 46]])

NS = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -frame)
F = vector(ZZ, [1, 0] + [0] * 17)
O = vector(ZZ, [-1, 1] + [0] * 17)
twice_minuscule = (2, 3, 4, 6, 5, 4, 3)
minus_P1 = vector(
    ZZ,
    [5, 1]
    + [-value for value in twice_minuscule]
    + [0] * 8
    + [1, 0],
)
D = O + minus_P1 - F
assert minus_P1 * NS * minus_P1 == -2
assert D * NS * D == 0 and D * NS * F == 2

payload = {
    "schema": "elkies-k3.h21-q6-section-descent.v1",
    "status": "PASS_H21_SECTION_Q_DEFINED",
    "inputs": {
        "h21_ancillary": {
            "path": str(args.source.relative_to(ROOT)),
            "sha256": SOURCE_SHA256,
        },
        "h3_frame": {
            "path": str(FRAME.relative_to(ROOT)),
            "sha256": FRAME_SHA256,
        },
    },
    "source_fibration": {
        "root_type": "A2+A6+E8",
        "root_determinants": {
            key: int(value) for key, value in source_root_determinants.items()
        },
        "ns_determinant": int(source_ns_determinant),
        "squarefree": True,
        "proper_integral_overlattice_possible": False,
        "generating_curves_defined_over": "QQ(r,s)",
        "neighbor_parameter_defined_over": "QQ(r,s)",
    },
    "target_fibration": {
        "root_type": "E7+E8",
        "P1_height": "21/2",
        "galois_action_on_P1": "fixed",
        "galois_action_on_minus_P1": "fixed",
        "minus_P1_class": coordinates(minus_P1),
        "q6_divisor_class": coordinates(D),
        "q6_divisor_old_fiber_degree": int(D * NS * F),
        "q6_divisor_defined_over": "QQ(r,s)",
    },
    "consequence": (
        "At the rational H21/H92 specialization, P1, -P1, and "
        "D=O+(-P1)-F are individually Q-defined.  The oriented Hilbert-cover "
        "quadratic field does not exchange these K3 section classes."
    ),
    "remaining_gate": (
        "Construct an explicit QQ-valued basis of H0(D), equivalently the "
        "q=6 divisor function and its Weierstrass pencil."
    ),
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "H21Q6DESCENT|source=A2+A6+E8|source_det=21|squarefree=1|"
    "neighbor_field=QQ(r,s)|P1_height=21/2|P1_galois=fixed|"
    "minusP1_galois=fixed|D_galois=fixed",
    flush=True,
)
print("H21Q6DESCENT|status=PASS_H21_SECTION_Q_DEFINED", flush=True)
