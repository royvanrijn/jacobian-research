#!/usr/bin/env sage
"""Recover the section profile of the preferred disc-43 A12+A3+A2 MW1 frame."""

from pathlib import Path

from sage.all import *


FRAME = Path("elkies-k3/data/fibrations/picard20_mw1_a12_a3_a2_frame.txt")
F = matrix(ZZ, [
    [ZZ(value) for value in line.split()]
    for line in FRAME.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
assert F.nrows() == 18 and F.det() == 43 and F.is_positive_definite()


def qform_from_gram(gram):
    coefficients = []
    for i in range(gram.nrows()):
        for j in range(i, gram.ncols()):
            coefficients.append(gram[i, i] // 2 if i == j else gram[i, j])
    return QuadraticForm(ZZ, gram.nrows(), coefficients)


def fractional_class(coordinates):
    return tuple(value - floor(value) for value in coordinates)


def mod_two(value):
    value = QQ(value)
    return value - 2 * floor(value / 2)


def class_order(point, exponent):
    zero = tuple(QQ(0) for _ in point)
    for order in divisors(exponent):
        if fractional_class(tuple(order * value for value in point)) == zero:
            return order
    raise AssertionError("class order exceeds component exponent")


half_roots = [
    vector(ZZ, root)
    for root in qform_from_gram(F).short_vector_list_up_to_length(2, True)[1]
]
roots = half_roots + [-root for root in half_roots]
graph = Graph()
graph.add_vertices(range(len(roots)))
for i in range(len(roots)):
    for j in range(i):
        if roots[i] * F * roots[j] != 0:
            graph.add_edge(i, j)
components = sorted(graph.connected_components(sort=False), key=len, reverse=True)
assert [len(component) for component in components] == [156, 12, 6]

component_bases = [
    matrix(ZZ, [roots[index] for index in component]).row_module().basis_matrix()
    for component in components
]
component_grams = [basis * F * basis.transpose() for basis in component_bases]
assert [
    (basis.rank(), abs(gram.det()))
    for basis, gram in zip(component_bases, component_grams)
] == [(12, 13), (3, 4), (2, 3)]

R = block_matrix([[basis] for basis in component_bases], subdivide=False)
GR = R * F * R.transpose()
C = (R * F).right_kernel_matrix()
GC = C * F * C.transpose()
assert R.rank() == 17 and C.rank() == 1 and GC == matrix(ZZ, [[6708]])

A = block_matrix([[R], [C]], subdivide=False)
index = abs(A.det())
assert index == 156
A_inverse = A.inverse()


def coset_key(point):
    return fractional_class(vector(QQ, point) * A_inverse)


zero = vector(ZZ, [0] * 18)
cosets = {coset_key(zero): zero}
queue = [zero]
head = 0
while head < len(queue) and len(cosets) < index:
    point = queue[head]
    head += 1
    for coordinate in range(18):
        unit = vector(ZZ, [0] * 18)
        unit[coordinate] = 1
        for sign in (1, -1):
            candidate = point + sign * unit
            key = coset_key(candidate)
            if key not in cosets:
                cosets[key] = candidate
                queue.append(candidate)
                if len(cosets) == index:
                    break
        if len(cosets) == index:
            break
assert len(cosets) == index

GR_inverse = GR.inverse()
GC_inverse = GC.inverse()


def project_mw(point):
    point = vector(QQ, point)
    return point - (point * F * R.transpose()) * GR_inverse * R


def coordinates_in_C(point):
    return (vector(QQ, point) * F * C.transpose()) * GC_inverse


target = vector(QQ, [QQ(1) / 156])
lifts = []
for representative in cosets.values():
    projected = coordinates_in_C(project_mw(representative))
    difference = target - projected
    if all(value.denominator() == 1 for value in difference):
        lifts.append(vector(ZZ, vector(QQ, representative) + difference * C))
assert len(lifts) == 1
lift = lifts[0]
assert project_mw(lift) == target * C
height = (target * C) * F * (target * C)
assert height == QQ(43) / 156

root_part = vector(QQ, lift) - target * C
root_coordinates = (root_part * F * R.transpose()) * GR_inverse
assert root_part == root_coordinates * R

data = []
left = 0
for basis, gram, exponent in zip(component_bases, component_grams, (13, 4, 3)):
    right = left + basis.nrows()
    coordinates = vector(QQ, root_coordinates[left:right])
    data.append((
        class_order(fractional_class(coordinates), exponent),
        mod_two(coordinates * gram * coordinates),
    ))
    left = right

assert data == [(13, QQ(4) / 13), (4, QQ(3) / 4), (3, QQ(2) / 3)]

# The A12 residue 4/13 is label 3 or 10, whose minimal correction is 30/13.
local_correction = QQ(30) / 13 + QQ(3) / 4 + QQ(2) / 3
zero_intersection = (height + local_correction - 4) / 2
assert local_correction == QQ(581) / 156
assert zero_intersection == 0

# root_det * height = 156 * 43/156 = 43 = |disc(NS)|, so the torsion group
# is trivial by Shioda's determinant formula.
assert 156 * height == 43

print("PICARD20MW1A12GLUE|frame_det=43|roots=A12+A3+A2|root_det=156", flush=True)
print("PICARD20MW1A12GLUE|height=43/156|torsion=trivial", flush=True)
print("PICARD20MW1A12GLUE|profile=A12:3_or_10,A3:1_or_3,A2:1_or_2|P.O=0", flush=True)
print("PICARD20MW1A12GLUE|expected_fibers=I13+I4+I3+4I1", flush=True)
print("PICARD20MW1A12GLUE|status=PASS", flush=True)
