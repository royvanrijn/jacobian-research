#!/usr/bin/env sage -python
"""Count the closest distinct rootless rational-bisection orbit pairs.

For rootless bisection classes B_w and B_v in U+(-M),

    B_w.B_v = (w-v).M.(w-v)/2 - 2.

Section translations change w and v independently by 2M.  Consequently an
unordered pair of translation orbits has minimum intersection zero precisely
when the XOR of its M/2M masks is represented by a norm-four lattice vector.
This script counts that finite frontier without materializing its potentially
large list of pairs.  It is a lattice prioritization for the later equation
search only; equality of quadratic extensions and base-change heights still
require bisection equations.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

from sage.all import QuadraticForm, ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
GRAM = ROOT / "elkies-k3/data/lattice/short_vector_basis_gram.txt"
DEFAULT_ORBITS = ROOT / "artifacts/generated-results/elkies-k3-rootless-bisection-orbits.tsv"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-rootless-bisection-disjoint-frontier.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def parity_mask(value):
    result = 0
    for index, entry in enumerate(value):
        if ZZ(entry) % 2:
            result |= 1 << index
    return result


def quadratic_form(gram):
    coefficients = []
    for row in range(gram.nrows()):
        for column in range(row, gram.ncols()):
            coefficients.append(gram[row, row] // 2 if row == column else gram[row, column])
    return QuadraticForm(ZZ, gram.nrows(), coefficients)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--orbits", type=Path, default=DEFAULT_ORBITS)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
args.orbits = args.orbits.resolve()
args.output = args.output.resolve()

gram = load_matrix(GRAM)
assert gram.nrows() == gram.ncols() == 17 and gram.is_positive_definite()
assert pari(gram).qfminim(2)[0] == 0
rows = list(csv.DictReader(args.orbits.open(newline=""), delimiter="\t"))
assert rows and {"orbit_mask", "short_basis_w"} <= set(rows[0])
orbits = {int(row["orbit_mask"], 0) for row in rows}
assert len(orbits) == len(rows) == 39120

# qfminim uses substantial PARI stack even for the small norm-four shell on
# this rank-17 lattice; allocation is local to this checker.
pari.allocatemem(4 * 1024**3)
shells = quadratic_form(gram).short_vector_list_up_to_length(4, True)
assert not shells[1]
norm_four_masks = {parity_mask(vector(ZZ, item)) for item in shells[2]}
assert 0 not in norm_four_masks

# The difference mask determines an unordered orbit pair uniquely.  Count an
# edge only in its increasing mask orientation, avoiding materialization of
# the disjointness graph.
edge_count = 0
active_difference_masks = 0
for delta in norm_four_masks:
    contribution = sum(1 for orbit in orbits if orbit < (orbit ^ delta) and (orbit ^ delta) in orbits)
    if contribution:
        active_difference_masks += 1
        edge_count += contribution

payload = {
    "schema": "elkies-k3.rootless-bisection-disjoint-frontier.v1",
    "status": "PASS_EXACT_ROOTLESS_BISECTION_DISJOINT_FRONTIER",
    "inputs": {
        "short_basis_gram": {"path": str(GRAM.relative_to(ROOT)), "sha256": digest(GRAM)},
        "bisection_orbits": {"path": str(args.orbits.relative_to(ROOT)), "sha256": digest(args.orbits)},
    },
    "formula": {
        "intersection": "B_w.B_v=(w-v).M.(w-v)/2-2",
        "translation_action": "w,v may independently change by 2M",
        "minimum_intersection_zero": "the orbit-mask XOR has a norm-four representative",
    },
    "frontier": {
        "translation_orbit_count": len(orbits),
        "norm_four_unoriented_difference_masks": len(norm_four_masks),
        "active_difference_masks": active_difference_masks,
        "unordered_orbit_pairs_with_minimum_intersection_zero": edge_count,
        "pair_list": "not materialized; count obtained by exact mask-XOR membership",
    },
    "boundary": (
        "This is a lattice-only priority graph. It does not construct any bisection "
        "equation, assert equality of quadratic extensions, compute a base-change "
        "height matrix, or establish a rank claim."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "R17BISECTDISJOINT|orbits={}|norm4_masks={}|active_masks={}|pairs={}|"
    "status=PASS_EXACT_ROOTLESS_BISECTION_DISJOINT_FRONTIER".format(
        len(orbits), len(norm_four_masks), active_difference_masks, edge_count
    ),
    flush=True,
)
