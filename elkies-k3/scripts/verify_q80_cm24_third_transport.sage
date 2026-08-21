#!/usr/bin/env sage
"""Transport the first three q80 neighbor classes through the CM24 embedding.

This is the exact marking gate that distinguishes the two small CM24 chord
completions of the third q=12 pencil.  The generic q80 frame embeds primitively
in the rank-18 CM24 frame.  We transport each abstract neighbor class back to
the q80 basis, apply that embedding, and then express it in the sequential CM
child bases.  Root invariants of the specialized third child select the live
equation-level completion without relying on discriminant order alone.
"""

import csv
from pathlib import Path

from sage.all import (
    Matrix, QuadraticForm, QQ, ZZ, block_diagonal_matrix, gcd, matrix, vector,
    xgcd,
)


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "elkies-k3/data/fibrations"
U = matrix(ZZ, [[0, 1], [1, 0]])


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def isotropic_mate(ns, fiber):
    pairings = list(ns*fiber)
    current = ZZ(0)
    coefficients = [ZZ(0)]*len(pairings)
    for index, pairing in enumerate(pairings):
        if not pairing:
            continue
        new_gcd, left, right = xgcd(current, ZZ(pairing))
        coefficients = [left*value for value in coefficients]
        coefficients[index] += right
        current = new_gcd
    assert abs(current) == 1
    if current == -1:
        coefficients = [-value for value in coefficients]
    mate = vector(ZZ, coefficients)
    mate -= ZZ(mate*ns*mate)//2*fiber
    assert fiber*ns*mate == 1 and mate*ns*mate == 0
    return mate


def neighbor_from_fiber(parent, fiber):
    ns = block_diagonal_matrix(U, -parent)
    fiber = vector(ZZ, fiber)
    assert fiber*ns*fiber == 0
    assert gcd([abs(ZZ(value)) for value in ns*fiber]) == 1
    mate = isotropic_mate(ns, fiber)
    complement = matrix(
        ZZ, [list(fiber*ns), list(mate*ns)]
    ).right_kernel_matrix()
    transition = matrix(ZZ, [list(fiber), list(mate)]+complement.rows())
    assert abs(transition.det()) == 1
    child = -(complement*ns*complement.transpose())
    return child, transition


def root_invariants(gram):
    half = QuadraticForm(ZZ, gram).short_vector_list_up_to_length(
        2, up_to_sign_flag=True
    )[1]
    basis = matrix(ZZ, [list(row) for row in half]).row_module().basis_matrix()
    root_gram = basis*gram*basis.transpose()
    return basis.nrows(), 2*len(half), abs(ZZ(root_gram.det()))


generic = load_matrix(DATA / "kumar_q80_e6_d5_a3_mw3_frame.txt")
cm24 = matrix(
    ZZ,
    [
        (2,0,-1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1),
        (0,2,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
        (-1,-1,2,-1,0,0,0,0,0,0,0,0,-1,0,0,0,0,-1),
        (0,0,-1,2,-1,0,0,0,0,0,0,0,1,0,0,0,0,1),
        (0,0,0,-1,2,0,0,0,0,0,0,0,-1,1,0,0,3,0),
        (0,0,0,0,0,2,0,-1,0,0,0,0,0,0,0,0,0,0),
        (0,0,0,0,0,0,2,0,-1,0,0,0,0,0,0,0,0,0),
        (0,0,0,0,0,-1,0,2,-1,0,0,0,0,0,0,0,0,0),
        (0,0,0,0,0,0,-1,-1,2,0,0,-1,0,1,0,0,-1,0),
        (0,0,0,0,0,0,0,0,0,2,0,1,0,1,0,0,-1,0),
        (0,0,0,0,0,0,0,0,0,0,158,0,32,0,0,0,0,32),
        (0,0,0,0,0,0,0,0,-1,1,0,4,0,0,0,0,0,0),
        (1,0,-1,1,-1,0,0,0,0,0,32,0,8,-1,0,0,-3,6),
        (0,0,0,0,1,0,0,0,1,1,0,0,-1,4,0,-1,4,5),
        (0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,6,-4,18),
        (0,0,0,0,0,0,0,0,0,0,0,0,0,-1,6,18,-8,12),
        (0,0,0,0,3,0,0,0,-1,-1,0,0,-3,4,-4,-8,24,-8),
        (1,0,-1,1,0,0,0,0,0,0,32,0,6,5,18,12,-8,130),
    ],
)
embedding = matrix(ZZ, 17, 18)
embedding[0] = vector(
    ZZ, (-41,-2,-3,-42,-2,0,0,0,0,0,-16,0,78,10,-2,5,10,1)
)
for row in range(1, 11):
    embedding[row, row-1] = 1
embedding[11,11] = 1
for row, column in ((12,13), (13,14), (14,15), (15,16)):
    embedding[row, column] = 1
embedding[16,12] = 1
embedding[16,17] = -1
assert generic == embedding*cm24*embedding.transpose()

ns_embedding = matrix(ZZ, 19, 20)
ns_embedding[0,0] = ns_embedding[1,1] = 1
ns_embedding[2:,2:] = embedding

with (DATA / "kumar_q80_to_rootless_path.tsv").open() as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))


def row_fiber(row):
    return vector(
        ZZ,
        [ZZ(row["a"]), ZZ(row["b"])]
        + list(map(ZZ, row["v"].split(","))),
    )


generic_frame = generic
generic_total = Matrix.identity(ZZ, 19)
cm_frame = cm24
cm_total = Matrix.identity(ZZ, 20)
specialized_fibers = []
for row in rows[:3]:
    generic_fiber = row_fiber(row)
    old_generic = generic_fiber*generic_total
    old_cm = old_generic*ns_embedding
    specialized = old_cm*cm_total.inverse()
    assert all(value in ZZ for value in specialized)
    specialized = vector(ZZ, specialized)
    specialized_fibers.append(tuple(specialized))
    generic_frame, generic_transition = neighbor_from_fiber(
        generic_frame, generic_fiber
    )
    cm_frame, cm_transition = neighbor_from_fiber(cm_frame, specialized)
    generic_total = generic_transition*generic_total
    cm_total = cm_transition*cm_total

assert generic_frame.det() == 948 and cm_frame.det() == 24
third_roots = root_invariants(cm_frame)
assert third_roots == (15, 90, 392)  # 2A6 + 3A1
print(
    f"Q80CM24THIRDTRANSPORT|specialized_fibers={tuple(specialized_fibers)}",
    flush=True,
)
print(
    f"Q80CM24THIRDTRANSPORT|third_root_rank={third_roots[0]}|"
    f"roots={third_roots[1]}|rootdet={third_roots[2]}|"
    f"ADE=2A6+3A1|geometric_MW={18-third_roots[0]}|status=PASS",
    flush=True,
)
