#!/usr/bin/env sage -python
"""Certify two exact low-degree chambers from the H3-chain D13/MW4 frame.

The source is the root-adapted D13 frame reached by the exact H3 q6,q8
chain.  The dominant q=4 witness is already in the effective D13 chamber.
Besides checking the displayed components, this script proves full nefness:

* section intersections reduce to a closest-vector problem in the D13 root
  lattice and the rank-four MW quotient; the exact minimum is two;
* a negative bisection would force its frame coordinate to equal the fiber
  coordinate, contradicting the parity of its self-intersection.

The resulting child has root data (13,158,26), hence A12+A1 and MW rank four.
This is a lateral presentation, not a rank-growing step.

The same exact checks are then applied to the first rank-growing shell.  Its
preferred q=24, degree-two orbit has child root data (12,264,4), hence D12
and MW rank five.  In that case the MW-quotient projection alone excludes a
negative section.
"""

import itertools
from pathlib import Path

from sage.all import *


BASE = Path(__file__).resolve().parents[1]
FRAME = BASE / "data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt"
WITNESS = vector(ZZ, (
    3, -1, 1, 4, 8, 5, 7, 6, 5, 4, -1, -3, 3, -1, 1, -1, 0,
))
EXPECTED_MW_PROJECTION = vector(ZZ, (-1, 1, -1, 0))
Q24_WITNESS = vector(ZZ, (
    0, 5, 0, 1, 2, 1, 2, 2, 2, 2, 4, 8, 2, 0, -1, 1, 1,
))
Q24_MW_PROJECTION = vector(ZZ, (0, -1, 1, 1))


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


def bezout_vector_for_pairing(ns, fiber):
    pairings = tuple(ns * fiber)
    current = ZZ(0)
    result = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(pairings):
        if value == 0:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        result = [left * entry for entry in result]
        result[index] += right
        current = divisor
    assert abs(current) == 1
    if current == -1:
        result = [-entry for entry in result]
    return vector(ZZ, result)


def neighbor_frame(ns, fiber):
    mate = bezout_vector_for_pairing(ns, fiber)
    mate -= (mate * ns * mate // 2) * fiber
    assert fiber * ns * mate == 1 and mate * ns * mate == 0
    complement = matrix(
        ZZ, [list(fiber * ns), list(mate * ns)]
    ).right_kernel_matrix()
    full_basis = matrix(
        ZZ, [list(fiber), list(mate)] + [list(row) for row in complement]
    )
    assert abs(full_basis.det()) == 1
    return -(complement * ns * complement.transpose())


def closest_d13_squared(target, standard_simple_isometry):
    """Return the exact distance from ``target`` to the D13 root lattice."""
    standard_cartan = CartanMatrix(["D", 13])
    ambient_simple = matrix(ZZ, 13, 13)
    for index in range(12):
        ambient_simple[index, index] = 1
        ambient_simple[index, index + 1] = -1
    ambient_simple[12, 11] = 1
    ambient_simple[12, 12] = 1
    assert ambient_simple * ambient_simple.transpose() == standard_cartan

    standard_coordinates = standard_simple_isometry.inverse() * target
    ambient_target = ambient_simple.transpose() * standard_coordinates

    # D13 is {x in ZZ^13 : sum(x) even}.  A closest vector has every
    # coordinate in {floor(target_i),ceil(target_i)}: moving any farther
    # coordinate by two preserves parity and strictly reduces the distance.
    choices = []
    for value in ambient_target:
        lower, upper = floor(value), ceil(value)
        choices.append((lower,) if lower == upper else (lower, upper))
    minimum = None
    for coordinates in itertools.product(*choices):
        if sum(coordinates) % 2:
            continue
        difference = vector(QQ, coordinates) - ambient_target
        norm = difference * difference
        if minimum is None or norm < minimum:
            minimum = norm
    assert minimum is not None
    return minimum


frame = load_gram(FRAME)
assert frame.nrows() == 17 and frame.det() == 948
root = frame[:13, :13]
root_mw = frame[:13, 13:]
height = (
    frame[13:, 13:]
    - frame[13:, :13] * root.inverse() * root_mw
)
assert root.det() == 4
assert height == matrix(QQ, [
    [QQ(3)/4, QQ(1)/4, -QQ(1)/4, 0],
    [QQ(1)/4, QQ(11)/4, QQ(1)/4, 1],
    [-QQ(1)/4, QQ(1)/4, QQ(11)/4, -1],
    [0, 1, -1, 46],
])

ns = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -frame)
old_fiber = vector(ZZ, [1, 0] + [0] * 17)
old_zero = vector(ZZ, [-1, 1] + [0] * 17)
divisor = vector(ZZ, [2, 2] + list(WITNESS))
assert divisor * ns * divisor == 0
assert divisor * ns * old_fiber == 2
assert divisor * ns * old_zero == 0
assert gcd(tuple(ns * divisor)) == 1
assert WITNESS * frame * WITNESS == 8
assert vector(ZZ, WITNESS[13:]) == EXPECTED_MW_PROJECTION

# The effective D13 simple roots are the negatives of the first thirteen
# frame basis vectors.  The sole nonzero label is node six.
effective_simple = tuple(
    vector(ZZ, [0, 0] + [-ZZ(index == node) for index in range(17)])
    for node in range(13)
)
simple_pairings = tuple(divisor * ns * curve for curve in effective_simple)
assert simple_pairings == (0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)

# Recover the highest D13 root directly in this pinned Cartan numbering.
half_roots = matrix(ZZ, pari(root).qfminim(2)[2]).transpose().rows()
roots = tuple(half_roots) + tuple(-row for row in half_roots)
positive_coordinate_roots = [row for row in roots if all(value >= 0 for value in row)]
highest = max(positive_coordinate_roots, key=lambda row: sum(row))
assert tuple(highest) == (2, 2, 1, 1, 2, 1, 2, 2, 2, 2, 2, 2, 2)
affine = old_fiber + vector(ZZ, [0, 0] + list(highest) + [0] * 4)
assert affine * ns * affine == -2 and divisor * ns * affine == 1

# Complete the section CVP along the primitive D13 roots.  Writing a frame
# vector as (r,m), only quotient vectors with ||m-z/2||_H^2<2 can threaten a
# negative section.  qfminim enumerates that finite shifted parity shell.
scaled_height = (4 * height).change_ring(ZZ)
z = EXPECTED_MW_PROJECTION
short_data = pari(scaled_height).qfminim(31)
half_short = matrix(ZZ, short_data[2]).transpose().rows()
n_candidates = [vector(ZZ, [0] * 4)] + list(half_short) + [-row for row in half_short]
n_candidates = [
    row for row in n_candidates
    if all((row[index] - z[index]) % 2 == 0 for index in range(4))
]
assert len(n_candidates) == 8

standard_cartan = CartanMatrix(["D", 13])
standard_simple_isometry = matrix(ZZ, pari(standard_cartan).qfisom(pari(root)))
assert abs(standard_simple_isometry.det()) == 1
assert standard_simple_isometry.transpose() * root * standard_simple_isometry == standard_cartan

root_coordinate = vector(QQ, WITNESS[:13])
section_distances = []
for n in n_candidates:
    m = vector(QQ, (n + z) / 2)
    quotient_difference = m - vector(QQ, z) / 2
    root_center = (
        root_coordinate / 2
        - root.inverse() * root_mw * quotient_difference
    )
    root_distance = closest_d13_squared(root_center, standard_simple_isometry)
    quotient_distance = QQ(n * scaled_height * n) / 16
    section_distances.append(root_distance + quotient_distance)

assert sorted(section_distances) == [QQ(2), QQ(2)] + [QQ(3)] * 6
assert min(section_distances) == 2

# If a degree-two (-2)-curve C=[k,2,w] were fixed, then D.C<0 and
# ||w-v||^2=2(D.C+1) force w=v and D.C=-1.  But C^2=-2 would require
# v^2=4k+2, impossible for v^2=8.  Together with the section CVP and all
# vertical components this proves nefness.
assert (WITNESS * frame * WITNESS - 2) % 4 != 0

child = neighbor_frame(ns, divisor)
minimum = pari(child).qfminim(2)
child_roots = matrix(ZZ, minimum[2]).transpose()
child_root_basis = child_roots.row_module().basis_matrix()
child_root_gram = child_root_basis * child * child_root_basis.transpose()
root_data = (
    child_root_basis.rank(),
    ZZ(minimum[0]),
    abs(child_root_gram.det()),
)
assert child.det() == 948 and root_data == (13, 158, 26)

print(
    "H3D13Q4|source=D13/MW4|q=4|ab=2,2|old_degree=2|O=0|"
    "reflections=0|component_pairings={}|affine=1|section_cvp_min=2|"
    "bisection_parity=1|nef=1|child=A12+A1/MW4|root_data=13,158,26|"
    "status=PASS".format(",".join(map(str, simple_pairings))),
    flush=True,
)

# Preferred first rank-growing orbit: q=24 with factor order (12,2).  It is
# also already in the effective D13 chamber.
q24_divisor = vector(ZZ, [12, 2] + list(Q24_WITNESS))
assert q24_divisor * ns * q24_divisor == 0
assert q24_divisor * ns * old_fiber == 2
assert q24_divisor * ns * old_zero == 10
assert gcd(tuple(ns * q24_divisor)) == 1
assert Q24_WITNESS * frame * Q24_WITNESS == 48
assert vector(ZZ, Q24_WITNESS[13:]) == Q24_MW_PROJECTION
q24_simple_pairings = tuple(
    q24_divisor * ns * curve for curve in effective_simple
)
assert q24_simple_pairings == (0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
assert q24_divisor * ns * affine == 1

# A negative section would require a quotient vector m with
# ||m-z/2||_H^2<2.  The complete parity-filtered qfminim list is empty, so no
# root-coordinate completion can make a section negative.
q24_z = Q24_MW_PROJECTION
q24_n_candidates = [vector(ZZ, [0] * 4)] + list(half_short) + [
    -row for row in half_short
]
q24_n_candidates = [
    row for row in q24_n_candidates
    if all((row[index] - q24_z[index]) % 2 == 0 for index in range(4))
]
assert q24_n_candidates == []

# The same degree-two norm identity excludes a negative bisection.
assert (Q24_WITNESS * frame * Q24_WITNESS - 2) % 4 != 0

q24_child = neighbor_frame(ns, q24_divisor)
q24_minimum = pari(q24_child).qfminim(2)
q24_roots = matrix(ZZ, q24_minimum[2]).transpose()
q24_root_basis = q24_roots.row_module().basis_matrix()
q24_root_gram = q24_root_basis * q24_child * q24_root_basis.transpose()
q24_root_data = (
    q24_root_basis.rank(),
    ZZ(q24_minimum[0]),
    abs(q24_root_gram.det()),
)
assert q24_child.det() == 948 and q24_root_data == (12, 264, 4)

print(
    "H3D13Q24|source=D13/MW4|q=24|ab=12,2|old_degree=2|O=10|"
    "reflections=0|component_pairings={}|affine=1|"
    "section_quotient_ball_empty=1|bisection_parity=1|nef=1|"
    "child=D12/MW5|root_data=12,264,4|status=PASS".format(
        ",".join(map(str, q24_simple_pairings))
    ),
    flush=True,
)
