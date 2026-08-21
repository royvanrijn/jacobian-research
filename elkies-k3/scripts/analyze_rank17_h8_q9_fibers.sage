#!/usr/bin/env sage
"""Classify the first fibers orthogonal to the intrinsic Humbert-8 class.

For ``h=(4,4,x)`` in ``U + (-R17)``, an orthogonal q=9 fiber must have
``(a,b)=(3,3)``: the alternative ``(1,9)`` is excluded by the exact dual
norm bound.  Completing the square centers its frame vector at ``2*x/3``;
the residual norm is two in the affine lattice ``ell*v=6``.  The close-vector
enumeration below recovers the thirteen representatives up to fiber sign,
then constructs every child and reports its complete root invariants.

This is deliberately a small affine-CVP calculation.  It does not call the
rank-17 PARI shell enumeration, which exceeds the configured PARI stack.
"""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from sage.all import (
    QQ,
    ZZ,
    block_diagonal_matrix,
    gcd,
    matrix,
    pari,
    vector,
    xgcd,
)
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]
U = matrix(ZZ, [[0, 1], [1, 0]])

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--write-artifact", action="store_true")
arguments = parser.parse_args()


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def bezout_vector(pairings):
    current = ZZ(0)
    coefficients = [ZZ(0)]*len(pairings)
    for index, pairing in enumerate(pairings):
        if pairing == 0:
            continue
        new_gcd, left, right = xgcd(current, ZZ(pairing))
        coefficients = [left*value for value in coefficients]
        coefficients[index] += right
        current = new_gcd
    assert abs(current) == 1
    return vector(
        ZZ, coefficients if current == 1 else [-value for value in coefficients]
    )


def neighbor(parent, qnorm, a, b, coordinates):
    ns = block_diagonal_matrix(U, -parent)
    fiber = vector(ZZ, [a, b]+list(coordinates))
    assert a*b == qnorm
    assert fiber*ns*fiber == 0
    assert gcd([abs(ZZ(value)) for value in ns*fiber]) == 1
    mate = bezout_vector(list(ns*fiber))
    mate -= ZZ(mate*ns*mate)//2*fiber
    complement = matrix(
        ZZ, [list(fiber*ns), list(mate*ns)]
    ).right_kernel_matrix()
    child = -(complement*ns*complement.transpose())
    assert child.det() == parent.det()
    return child


def root_data(frame):
    minimum = pari(frame).qfminim(2)
    count = ZZ(minimum[0])
    if not count:
        return (0, 0, 1)
    roots = matrix(ZZ, minimum[2]).transpose()
    basis = roots.row_module().basis_matrix()
    gram = basis*frame*basis.transpose()
    return (basis.rank(), count, abs(gram.det()))


rank17 = load_matrix(ROOT / "elkies-k3/data/lattice/rank17_gram.txt")
x = vector(ZZ, (-1, 0, -3, 0, 2, -2, 1, -2, 1, 1, 0, 1, 0, 0, -2, -2, 2))
ell = vector(ZZ, x*rank17/4)
assert ell in ZZ**17 and gcd(ell) == 1
assert ell*rank17.inverse()*ell == QQ(9)/4

# For (a,b)=(1,9), ell.v=10 while v.R.v=18, contradicting the
# exact Cauchy--Schwarz lower bound 10^2/(9/4)>18.
assert QQ(100)/(QQ(9)/4) > 18

kernel = matrix(ZZ, [list(ell)]).right_kernel_matrix()
assert kernel.nrows() == kernel.rank() == 16

# A small Bezout solution of ell.v=6.
particular = vector(ZZ, [0]*14+[-6, 0, 0])
assert ell*particular == 6
center = vector(QQ, 2*x/3)
assert ell*center == 6
assert center*rank17*center == 16

affine_lattice = IntegralLattice(rank17, basis=kernel)
close = affine_lattice.enumerate_close_vectors(center-particular)

# BEST_N close-vector enumeration is requested through 128 candidates.  Its
# exact distance spectrum has a gap after the thirteen norm-two residuals;
# this is the bounded affine-CVP certificate used here.
residuals = []
representatives = []
for _ in range(128):
    lattice_vector = vector(ZZ, next(close))
    representative = particular+lattice_vector
    distance = (
        (representative-center)*rank17*(representative-center)
    )
    residuals.append(distance)
    if distance == 2:
        representatives.append(representative)
assert Counter(residuals) == Counter({QQ(2): 13, QQ(4): 115})
representatives = sorted(set(map(tuple, representatives)))
assert len(representatives) == 13
assert all(
    ell*vector(ZZ, row) == 6
    and vector(ZZ, row)*rank17*vector(ZZ, row) == 18
    and max(map(abs, row)) <= 3
    for row in representatives
)

rows = []
for representative in representatives:
    child = neighbor(rank17, ZZ(9), ZZ(3), ZZ(3), vector(ZZ, representative))
    roots = root_data(child)
    rows.append((representative, roots))
    print(
        "R17H8Q9|"
        f"v={representative}|root_rank={roots[0]}|roots={roots[1]}|"
        f"rootdet={roots[2]}|MW={17-roots[0]}",
        flush=True,
    )

rootless = [row for row in rows if row[1][0] == 0]
best_rank = min(row[1][0] for row in rows)
best = [row for row in rows if row[1][0] == best_rank]
print(
    "R17H8Q9|"
    f"classes_up_to_sign={len(rows)}|rootless={len(rootless)}|"
    f"best_root_rank={best_rank}|best_count={len(best)}|"
    f"best_root_data={tuple(sorted(set(row[1] for row in best)))}|"
    "status=PASS_BOUNDED_AFFINE_CVP_CLASSIFICATION",
    flush=True,
)

if arguments.write_artifact:
    output = {
        "schema": "rank17-h8-orthogonal-q9-fibers-v1",
        "h": [4, 4]+list(map(int, x)),
        "ell": list(map(int, ell)),
        "factor_pair": [3, 3],
        "classes_up_to_sign": len(rows),
        "affine_cvp_first_candidate_distance_counts": {
            str(distance): count
            for distance, count in sorted(Counter(residuals).items())
        },
        "rows": [
            {
                "v": list(map(int, representative)),
                "root_rank": int(roots[0]),
                "roots": int(roots[1]),
                "rootdet": int(roots[2]),
                "MW": int(17-roots[0]),
            }
            for representative, roots in rows
        ],
        "rootless_count": len(rootless),
        "best_root_rank": int(best_rank),
        "reproduce": (
            "sage elkies-k3/scripts/analyze_rank17_h8_q9_fibers.sage "
            "--write-artifact"
        ),
    }
    output_path = (
        ROOT / "artifacts/generated-results/rank17-h8-orthogonal-q9-fibers.json"
    )
    encoded = json.dumps(output, indent=2, sort_keys=True, default=int)+"\n"
    output_path.write_text(encoded)
    digest = hashlib.sha256(encoded.encode()).hexdigest()
    print(
        "R17H8Q9|"
        f"artifact={output_path}|sha256={digest}|status=PASS_ARTIFACT_WRITE",
        flush=True,
    )
