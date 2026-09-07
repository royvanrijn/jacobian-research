#!/usr/bin/env sage -python
"""Reconcile the H3 q8 dominant, classifier-nef, and old-fibre chamber representatives.

This audit compares the q8 D13-dominant orbit hit, the independently certified
classifier nef representative, and the class used by the historical generic
RR ambient, all in the same original H3 Neron--Severi coordinates.

It reproduces two important facts:

* finite old-fibre root reduction of the dominant D13 hit gives the historical
  old-fibre-degree-18 ambient class;
* finite old-fibre root reduction of the classifier nef representative gives
  the old-fibre-degree-16 class found by the later horizontal/fibre chamber
  audit.

Thus degree 18 and degree 16 are different Weyl representatives of the q8
orbit.  The latter is the classifier-nef representative placed in the actual
old E7+E8 component chamber; this is a lattice statement only and does not
validate any later hand-translated q6^8 local RR compiler.
"""

import argparse
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
ORBITS = GEN / "elkies-k3-h3-q6-q8-orbits.json"
AMBIENT = GEN / "elkies-k3-h92-q8-generic-rr-ambient.json"
DEFAULT_CHAMBER = GEN / "zz-h92-q8-complete-chamber-reduction.json"


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--orbits", type=Path, default=ORBITS)
parser.add_argument("--ambient", type=Path, default=AMBIENT)
parser.add_argument("--chamber", type=Path, default=DEFAULT_CHAMBER)
args = parser.parse_args()

frame = load_gram(FRAME)
ns = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -frame)
F = vector(ZZ, [1, 0] + [0] * 17)
simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(index == node) for index in range(17)])
    for node in range(15)
)
zero = vector(ZZ, [0] * 19)
highest_e7 = (2, 2, 3, 4, 3, 2, 1)
highest_e8 = (2, 3, 4, 6, 5, 4, 3, 2)
affine_e7 = F - sum((c * simple[i] for i, c in enumerate(highest_e7)), zero)
affine_e8 = F - sum((c * simple[7+i] for i, c in enumerate(highest_e8)), zero)

orbits = json.loads(args.orbits.read_text())
ambient = json.loads(args.ambient.read_text())
assert orbits["status"] == "PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION"
assert ambient["status"] == "PASS_EXACT_Q8_GENERIC_RR_AMBIENT"

dominant = next(
    hit for hit in orbits["q8"]["d13_mw4_hits"]
    if hit["mw_projection"] == [0, -2, 0]
)
nef = orbits["q8"]["nef_representative"]
classes = {
    "dominant_raw": vector(ZZ, dominant["fiber_source_h3_ns"]),
    "classifier_nef": vector(ZZ, nef["fiber_source_h3_ns"]),
    "ambient_degree18": vector(ZZ, ambient["source_q8_lattice_class"]),
}
if args.chamber.exists():
    chamber = json.loads(args.chamber.read_text())
    classes["experimental_degree16"] = vector(ZZ, chamber["final_class"])


def finite_reduce(value):
    value = vector(ZZ, value)
    steps = []
    for _ in range(1000):
        negative = [
            (index, int(value * ns * root))
            for index, root in enumerate(simple)
            if value * ns * root < 0
        ]
        if not negative:
            return value, steps
        index, pairing = negative[0]
        value += pairing * simple[index]
        steps.append((index + 1, pairing))
    raise RuntimeError("finite-root reduction did not terminate")


def summary(name, value):
    return (
        f"Q8REP|name={name}|square={int(value*ns*value)}|oldF={int(value*ns*F)}|"
        f"simple={','.join(str(int(value*ns*root)) for root in simple)}|"
        f"affine={int(value*ns*affine_e7)},{int(value*ns*affine_e8)}|"
        f"vector={','.join(map(str, value))}"
    )


for name, value in classes.items():
    print(summary(name, value), flush=True)

reduced_dominant, dominant_steps = finite_reduce(classes["dominant_raw"])
reduced_nef, nef_steps = finite_reduce(classes["classifier_nef"])
assert reduced_dominant == classes["ambient_degree18"]
assert int(reduced_dominant * ns * F) == 18

if "experimental_degree16" in classes:
    assert reduced_nef == classes["experimental_degree16"]
assert int(reduced_nef * ns * F) == 16

print(
    "Q8REP_FINITEREDUCE|name=dominant_raw|"
    f"steps={len(dominant_steps)}|oldF=18|equals_ambient18=1|"
    f"vector={','.join(map(str, reduced_dominant))}",
    flush=True,
)
print(
    "Q8REP_FINITEREDUCE|name=classifier_nef|"
    f"steps={len(nef_steps)}|oldF=16|equals_degree16={int('experimental_degree16' in classes and reduced_nef == classes['experimental_degree16'])}|"
    f"vector={','.join(map(str, reduced_nef))}",
    flush=True,
)
print(
    "Q8REP_RESULT|dominant_to_degree18=1|classifier_nef_to_degree16=1|"
    "status=PASS_EXACT_Q8_REPRESENTATIVE_RECONCILIATION",
    flush=True,
)
