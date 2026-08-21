#!/usr/bin/env sage
"""Verify the CM-43 q=60 factorization through q=8 and q=9.

The q=60 fiber is an automorphism of the Picard-rank-20, discriminant-43
Kumar surface.  In the Q79-marked glue-211 closure, a complete
W(E7)xW(E8) orbit classification finds a primitive q=8 neighbor for which
that fiber has coordinates (a,b)=(3,3), hence gives a second q=9 neighbor.
The intermediate frame already has roots E7+E8 and MW rank three.

For comparison, this script pins both the cheapest q=8 orbit having D9+E7
roots and the P3-marked orbit corresponding to the explicit source-level
Humbert-8 two-neighbor.  They need q=289 and q=2,466,464 respectively to
reach the same q=60 fiber, so neither is the arithmetically shortest
factorization of this automorphism.
"""

from pathlib import Path

from sage.all import *


BASE = Path(__file__).resolve().parents[1]


def load_gram(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ],
    )


U = matrix(ZZ, ((0, 1), (1, 0)))
source = load_gram(
    BASE / "data/fibrations/kumar_cm43_marked_e7e8_mw3_frame.txt"
)
assert source.nrows() == 18 and source.det() == 43
source_ns = block_diagonal_matrix(U, -source)

# This is the exact q=60 fiber transported to the CM-43 Kumar closure.
q60_fiber = vector(ZZ, (
    5, 12,
    0, 0, -1, -1, -1, -1, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0,
))
assert q60_fiber*source_ns*q60_fiber == 0


def bezout_vector_for_pairing(ns, fiber):
    current = ZZ(0)
    coefficients = [ZZ(0)]*ns.nrows()
    for index, value in enumerate(ns*fiber):
        if value == 0:
            continue
        new_gcd, old_scale, new_scale = xgcd(current, ZZ(value))
        coefficients = [old_scale*entry for entry in coefficients]
        coefficients[index] += new_scale
        current = new_gcd
    assert abs(current) == 1
    if current == -1:
        coefficients = [-entry for entry in coefficients]
    return vector(ZZ, coefficients)


def neighbor(frame, a, b, witness):
    ns = block_diagonal_matrix(U, -frame)
    fiber = vector(ZZ, [a, b]+list(witness))
    assert fiber*ns*fiber == 0
    mate = bezout_vector_for_pairing(ns, fiber)
    mate_square = ZZ(mate*ns*mate)
    assert mate_square % 2 == 0
    mate -= (mate_square//2)*fiber
    assert fiber*ns*mate == 1 and mate*ns*mate == 0
    kernel = matrix(ZZ, [list(fiber*ns), list(mate*ns)]).right_kernel_matrix()
    child = -(kernel*ns*kernel.transpose())
    basis = matrix(ZZ, [list(fiber), list(mate)]+[list(row) for row in kernel])
    assert abs(basis.det()) == 1
    assert basis*ns*basis.transpose() == block_diagonal_matrix(U, -child)
    assert child.is_positive_definite() and child.det() == 43
    return child, basis


def root_components(gram):
    """Return exact (rank, signed-root-count, determinant) components."""
    result = pari(gram).qfminim(2)
    signed_count = ZZ(result[0])
    half_roots = [
        vector(ZZ, column) for column in matrix(ZZ, result[2]).columns()
    ]
    assert 2*len(half_roots) == signed_count
    parents = list(range(len(half_roots)))

    def find(index):
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left, right):
        left = find(left)
        right = find(right)
        if left != right:
            parents[right] = left

    for left in range(len(half_roots)):
        for right in range(left):
            if half_roots[left]*gram*half_roots[right] != 0:
                union(left, right)
    groups = {}
    for index, root in enumerate(half_roots):
        groups.setdefault(find(index), []).append(root)
    answer = []
    for roots in groups.values():
        basis = matrix(ZZ, [list(root) for root in roots]).row_module().basis_matrix()
        root_gram = basis*gram*basis.transpose()
        answer.append((basis.rank(), 2*len(roots), abs(ZZ(root_gram.det()))))
    return sorted(answer)


# Arithmetic-optimal q=8 orbit.  The complete orbit certificate is in
# classify_kumar_cm43_q8_orbits.sage.
q8_witness = vector(ZZ, (
    156, -78, 0, 0, -78, 0, -78, 0, 0,
    0, 0, 0, 0, 0, 0, -1, -155, -32,
))
q8_child, q8_basis = neighbor(source, 2, 4, q8_witness)
assert root_components(q8_child) == [(7, 126, 2), (8, 240, 1)]

q60_in_q8 = vector(ZZ, q60_fiber*q8_basis.inverse())
expected_q60_in_q8 = vector(ZZ, (
    3, 3, 468, -1, -1, 467, 0, 0, 0, 0,
    0, 0, 0, 0, 206, -147, -101, -329, 30, 167,
))
assert q60_in_q8 == expected_q60_in_q8
assert q60_in_q8[0]*q60_in_q8[1] == 9

q9_child, q9_basis = neighbor(
    q8_child, q60_in_q8[0], q60_in_q8[1], q60_in_q8[2:]
)
assert root_components(q9_child) == [(7, 126, 2), (8, 240, 1)]
return_isometry = QuadraticForm(ZZ, q9_child).is_globally_equivalent_to(
    QuadraticForm(ZZ, source), return_matrix=True
)
assert return_isometry is not False and abs(return_isometry.det()) == 1
assert return_isometry.transpose()*q9_child*return_isometry == source


# Cheapest root-equivalent D9+E7 q=8 orbit.  Root invariants alone do not make
# it the source-level two-neighbor; the P3 marking below distinguishes that.
d9e7_witness = vector(ZZ, (
    422, -211, 0, 0, -211, 0, -211, 1, 1,
    2, 2, 2, 2, 2, 2, -1, -420, -87,
))
d9e7_child, d9e7_basis = neighbor(source, 2, 4, d9e7_witness)
assert root_components(d9e7_child) == [(7, 126, 2), (9, 144, 4)]
d9e7_q60 = vector(ZZ, q60_fiber*d9e7_basis.inverse())
expected_d9e7_q60 = vector(ZZ, (
    17, 17, 7174, -1, -1, 7173, -1, -34, -68, -68,
    -7319, -78, -1589, 2135, -85, -8264, -1538, 2958, -5983, -1521,
))
assert d9e7_q60 == expected_d9e7_q60
assert d9e7_q60[0]*d9e7_q60[1] == 289

# The source-level D9+E7 inverse neighbor is distinguished from the other
# nine root-equivalent q=8 orbits by horizontal projection +P3.  In this
# marked frame P3 is represented by (0^15,4,0,-1); subtracting it from the
# witness lands in the rational span of the E7+E8 root lattice.
source_d9e7_witness = vector(ZZ, (
    0, 0, 0, 0, 0, 0, 0, 1, 1,
    1, 1, 1, 1, 2, 3, 4, 0, -1,
))
source_d9e7_child, source_d9e7_basis = neighbor(
    source, 2, 4, source_d9e7_witness
)
assert root_components(source_d9e7_child) == [(7, 126, 2), (9, 144, 4)]
source_roots_result = pari(source).qfminim(2)
source_root_basis = matrix(
    ZZ,
    [list(column) for column in matrix(ZZ, source_roots_result[2]).columns()],
).row_module().basis_matrix()
assert source_root_basis.rank() == 15
height4_frame = vector(ZZ, [0]*15+[4, 0, -1])
assert source_d9e7_witness-height4_frame in source_root_basis.row_space(QQ)
source_d9e7_q60 = vector(ZZ, q60_fiber*source_d9e7_basis.inverse())
expected_source_d9e7_q60 = vector(ZZ, (
    56056, 44, 0, -1, -1, -1, -1, -1, -112112, -112112,
    -224224, 1, 101046, 121902, 101684, 223899, -380336, 203681,
    3564, 2288,
))
assert source_d9e7_q60 == expected_source_d9e7_q60
assert source_d9e7_q60[0]*source_d9e7_q60[1] == 2466464

print(
    "KUMARCM43Q8Q9|first_q=8|intermediate=E7+E8|MW=3"
    f"|first_witness={tuple(q8_witness)}",
    flush=True,
)
print(
    "KUMARCM43Q8Q9|second_q=9|second_ab=3,3|return=E7+E8/MW3"
    f"|q60_in_child={tuple(q60_in_q8)}|isometry=1",
    flush=True,
)
print(
    "KUMARCM43D9E7|first_q=8|intermediate=D9+E7|MW=2"
    f"|cheapest_q60_second_q=289|q60_in_child={tuple(d9e7_q60)}",
    flush=True,
)
print(
    "KUMARCM43D9E7|source_marking=horizontal_P3|first_q=8"
    f"|q60_second_q=2466464|q60_in_child={tuple(source_d9e7_q60)}",
    flush=True,
)
print("KUMARCM43Q8Q9|status=PASS", flush=True)
