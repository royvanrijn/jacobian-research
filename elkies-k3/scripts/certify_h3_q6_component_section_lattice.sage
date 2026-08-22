#!/usr/bin/env sage -python
"""Pin the two old E7 components which become low-height q=6 sections.

This is a lattice-side bridge for the equation-level component transport.  It
does not replace resolved-chart coordinates: its purpose is to state exactly
which curves those coordinates must represent.  In the first H3 q=6 chamber,
the old zero section, the E7_7 exceptional component, and the old E7 affine
component all have degree one against D.  With O chosen as child zero, E7_7
and E7_7-affine_E7 have the first 2-by-2 principal block of the predicted
Mordell--Weil height Gram.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, gcd, identity_matrix, matrix, pari, vector, xgcd


ROOT = Path(__file__).resolve().parents[2]
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h3-q6-component-sections.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_gram(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()]
         for line in path.read_text().splitlines()
         if line.strip() and not line.startswith("#")],
    )


def bezout_isotropic_mate(ns, fiber):
    current = ZZ(0)
    row = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(ns * fiber):
        if not value:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        row = [left * entry for entry in row]
        row[index] += right
        current = divisor
    assert abs(current) == 1
    if current == -1:
        row = [-entry for entry in row]
    mate = vector(ZZ, row)
    mate -= (mate * ns * mate // 2) * fiber
    assert mate * ns * mate == 0 and fiber * ns * mate == 1
    return mate


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()

frame = load_gram(FRAME)
assert frame.nrows() == 17 and frame.det() == 948
ns = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -frame)
fiber = vector(ZZ, [1, 0] + [0] * 17)
zero = vector(ZZ, [-1, 1] + [0] * 17)
simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(index == node) for index in range(17)])
    for node in range(15)
)

# Replay the recorded reduction rather than using only its endpoint.
raw = vector(ZZ, [3, 2] + [
    0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
])
divisor = vector(ZZ, raw)
reflection_nodes = (1, 2, 4, 3, 5, 4, 2, 6, 5, 4, 3, 1, 7, 6, 5, 4, 2, 3, 4, 5, 6, 7)
for node in reflection_nodes:
    pairing = divisor * ns * simple[node - 1]
    assert pairing == -1
    divisor += pairing * simple[node - 1]
assert divisor * ns * divisor == 0
assert gcd(tuple(ns * divisor)) == 1

highest_e7 = (2, 2, 3, 4, 3, 2, 1)
affine_e7 = fiber - sum(
    (coefficient * simple[index] for index, coefficient in enumerate(highest_e7)),
    vector(ZZ, [0] * 19),
)
e7_7 = simple[6]
for curve in (zero, e7_7, affine_e7):
    assert curve * ns * curve == -2
    assert divisor * ns * curve == 1

# Construct the child root projection directly from the exact primitive U
# split.  This avoids assuming a Kodaira label while retaining the already
# proved E8+E6 root lattice as an integral sublattice of NS.
mate = bezout_isotropic_mate(ns, divisor)
orthogonal = matrix(ZZ, [list(divisor * ns), list(mate * ns)]).right_kernel_matrix()
child = -(orthogonal * ns * orthogonal.transpose())
roots = matrix(ZZ, pari(child).qfminim(2)[2]).transpose().row_module().basis_matrix()
assert roots.rank() == 14
root_curves = roots * orthogonal
root_gram = root_curves * ns * root_curves.transpose()
assert abs(root_gram.det()) == 3
projection = identity_matrix(QQ, 19) - ns * root_curves.transpose() * root_gram.inverse() * root_curves


def shioda_projection(section):
    """Project a degree-one curve relative to the old zero child section."""
    correction = section * ns * zero + 2
    degree_zero = section - zero - correction * divisor
    assert degree_zero * ns * divisor == degree_zero * ns * zero == 0
    return vector(QQ, degree_zero) * projection


p = shioda_projection(e7_7)
a = shioda_projection(affine_e7)
q = p - a
height_gram = matrix(QQ, [
    [-p * ns * p, -p * ns * q],
    [-q * ns * p, -q * ns * q],
])
expected = matrix(QQ, [[QQ(8) / 3, QQ(1) / 3], [QQ(1) / 3, QQ(8) / 3]])
assert height_gram == expected

payload = {
    "schema": "elkies-k3.h3-q6-component-sections.v1",
    "status": "PASS_EXACT_Q6_COMPONENT_SECTION_IDENTITIES",
    "inputs": {"frame": {"path": str(FRAME.relative_to(ROOT)), "sha256": digest(FRAME)}},
    "divisor": {
        "square": int(divisor * ns * divisor),
        "primitive": True,
        "weyl_reflections": list(reflection_nodes),
    },
    "curves": {
        "child_zero": "old O",
        "first_section": "old E7_7 exceptional component",
        "second_pre_difference": "old affine E7 component",
        "second_section": "MW difference E7_7-affine_E7 relative to old O",
        "D_intersections": {
            "old_O": int(divisor * ns * zero),
            "E7_7": int(divisor * ns * e7_7),
            "affine_E7": int(divisor * ns * affine_e7),
        },
    },
    "root_lattice": {"rank": 14, "determinant": 3, "type": "E8+E6"},
    "height_gram": [[str(value) for value in row] for row in height_gram.rows()],
    "boundary": "This fixes the two resolved old E7 curves to transport. Exact child-equation coordinates and the third rank-three generator remain separate gates.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("H3Q6COMPONENTS|D.O=1|D.E7_7=1|D.affineE7=1|height=((8/3,1/3),(1/3,8/3))|status=PASS_EXACT_Q6_COMPONENT_SECTION_IDENTITIES", flush=True)
