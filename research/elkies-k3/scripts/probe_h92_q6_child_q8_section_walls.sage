#!/usr/bin/env sage -python
"""Exhaustively test old-fibre section walls of a marked degree-two class.

For a class ``D=(a,2,v)`` in the old-fibre frame
``U + (-L)``, every section has the form ``C=(k,1,x)`` with
``x^2=2k+2``.  Setting ``y=2*x-v`` gives the exact identity

    D.C = y^2/4 - 2.

Thus every negative section wall occurs among the vectors ``y`` in the fixed
coset ``v (mod 2L)`` of norm strictly less than 8.  The PARI short-vector
enumeration below is therefore complete for old-fibre sections.  For a
bisection ``C=(k,2,x)``, the companion identity
``(x-v)^2=2*(D.C+1)`` shows that a negative intersection would force
``D.C=-1`` and ``x=v``; then ``x^2=4*k+2`` conflicts with ``v^2=4*a``.
The script records that universal parity exclusion too.  Vertical component
or affine-component walls must still be checked separately.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, gcd, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TARGET = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-component-nef-physical-root-target.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-component-nef-section-walls.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
args.target = args.target.resolve()
args.output = args.output.resolve()

target = json.loads(args.target.read_text())
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
data = target["selected_q8"]["physical_old_fibre_coordinates"]
coordinates = vector(ZZ, data["coordinates"])
assert coordinates[1] == 2
a = coordinates[0]
v = vector(ZZ, coordinates[2:])
lattice = matrix(ZZ, data["positive_rank17_gram"])
assert lattice.is_positive_definite()
assert v * lattice * v == 4 * a
primitive = gcd(tuple(coordinates)) == 1
assert primitive

# qfminim(bound) returns every vector of norm at most ``bound``; include both
# signs because its returned list is only one representative of each pair.
short = pari(lattice).qfminim(7)
half = matrix(ZZ, short[2]).transpose().rows()
short_vectors = {tuple([0] * lattice.nrows())}
for row in half:
    short_vectors.add(tuple(row))
    short_vectors.add(tuple(-vector(ZZ, row)))
short_vectors = tuple(vector(ZZ, row) for row in short_vectors)
assert all(y * lattice * y <= 7 for y in short_vectors)

walls = []
coset_short_vectors = 0
for y in sorted(short_vectors, key=lambda row: (row * lattice * row, tuple(row))):
    if any((y[index] - v[index]) % 2 for index in range(lattice.nrows())):
        continue
    coset_short_vectors += 1
    x = vector(ZZ, (y + v) / 2)
    x_norm = x * lattice * x
    assert (x_norm - 2) % 2 == 0
    k = (x_norm - 2) // 2
    intersection = a + 2 * k - v * lattice * x
    assert QQ(y * lattice * y) / 4 - 2 == intersection
    if intersection < 0:
        walls.append({
            "short_coset_vector": list(map(int, y)),
            "norm": int(y * lattice * y),
            "section_rank17_coordinate": list(map(int, x)),
            "section_old_fibre_coordinate": [int(k), 1],
            "intersection": int(intersection),
        })

vertical_degrees = (
    target["selected_q8"]["E6"]["component_degrees"]
    + [target["selected_q8"]["E6"]["affine_component_degree"]]
    + target["selected_q8"]["E8"]["component_degrees"]
    + [target["selected_q8"]["E8"]["affine_component_degree"]]
)
vertical_nonnegative = all(value >= 0 for value in vertical_degrees)
nef = vertical_nonnegative and not walls

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-section-walls.v1",
    "status": (
        "PASS_PRIMITIVE_NEF_DEGREE_TWO_CLASS" if nef
        else "PASS_NO_NEGATIVE_OLD_SECTION_WALL" if not walls
        else "WITNESS_NEGATIVE_OLD_SECTION_WALL"
    ),
    "inputs": {"target": {"path": str(args.target.relative_to(ROOT)), "sha256": digest(args.target)}},
    "identity": "D.C=(2*x-v)^2/4-2 for D=(a,2,v), C=(k,1,x)",
    "enumeration": {
        "short_vector_bound": 7,
        "short_vectors_including_signs_and_zero": len(short_vectors),
        "short_vectors_in_required_coset": coset_short_vectors,
        "complete_for_negative_old_fibre_sections": True,
    },
    "bisection_parity": {
        "identity": "(x-v)^2=2*(D.C+1) for C=(k,2,x)",
        "negative_intersection_forces": "D.C=-1 and x=v",
        "contradiction": "x^2=4*k+2 whereas v^2=4*a",
        "negative_old_fibre_bisection_wall_impossible": True,
    },
    "nef_certificate": {
        "primitive": primitive,
        "vertical_component_degrees_including_affine": vertical_degrees,
        "vertical_nonnegative": vertical_nonnegative,
        "reduction": (
            "If an irreducible curve C has D.C<0, then C is a fixed component "
            "of the effective class D, so D-C is effective and 0<=F.(D-C)=2-F.C. "
            "Hence C is vertical, a section, or a bisection; the three cases are "
            "checked above."
        ),
        "conclusion": (
            "D is primitive nef isotropic and therefore defines a genus-one pencil "
            "whose generic member is a degree-two multisection of the old fibration."
            if nef else "No global nefness conclusion is made."
        ),
    },
    "negative_section_walls": walls,
    "boundary": "This is a lattice certificate only.  It does not construct the genus-one pencil equation, a bisection cover, an extension collision, or a generic rank claim.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDQ8SECTIONS|representative={}|short={}|coset_short={}|negative={}|vertical_nef={}|status={}".format(
        target["normalization"]["representative"], len(short_vectors), coset_short_vectors,
        len(walls), int(vertical_nonnegative), payload["status"],
    ),
    flush=True,
)
