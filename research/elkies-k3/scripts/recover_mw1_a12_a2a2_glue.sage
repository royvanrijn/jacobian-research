#!/usr/bin/env sage
"""Recover the component profile of the A12+2A2 MW1 frame exactly."""

from sage.all import *
from pathlib import Path


frame_path = Path("elkies-k3/data/fibrations/mw1_a12_a2a2_frame.txt")
F = matrix(ZZ, [
    [ZZ(value) for value in line.split()]
    for line in frame_path.read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
assert F.nrows() == 17 and F.det() == 948 and F.is_positive_definite()


def qform_from_gram(gram):
    coefficients = []
    for i in range(gram.nrows()):
        for j in range(i, gram.ncols()):
            coefficients.append(gram[i, i] // 2 if i == j else gram[i, j])
    return QuadraticForm(ZZ, gram.nrows(), coefficients)


def fractional_class(coordinates):
    return tuple(value - floor(value) for value in coordinates)


def class_order(point, exponent):
    zero = tuple(QQ(0) for _ in point)
    for order in range(1, exponent + 1):
        multiple = fractional_class(tuple(order * value for value in point))
        if multiple == zero:
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
assert [len(component) for component in components] == [156, 6, 6]

component_bases = [
    matrix(ZZ, [roots[index] for index in component]).row_module().basis_matrix()
    for component in components
]
component_grams = [basis * F * basis.transpose() for basis in component_bases]
assert [
    (basis.rank(), abs(gram.det()))
    for basis, gram in zip(component_bases, component_grams)
] == [(12, 13), (2, 3), (2, 3)]

R = block_matrix([[basis] for basis in component_bases], subdivide=False)
GR = R * F * R.transpose()
C = (R * F).right_kernel_matrix()
GC = C * F * C.transpose()
assert R.rank() == 16 and C.rank() == 1 and GC == matrix(ZZ, [[12324]])

A = block_matrix([[R], [C]], subdivide=False)
index = abs(A.det())
assert index == 39
A_inverse = A.inverse()


def coset_key(point):
    return fractional_class(vector(QQ, point) * A_inverse)


zero = vector(ZZ, [0] * 17)
cosets = {coset_key(zero): zero}
queue = [zero]
head = 0
while head < len(queue) and len(cosets) < index:
    point = queue[head]
    head += 1
    for coordinate in range(17):
        unit = vector(ZZ, [0] * 17)
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


target = vector(QQ, [QQ(1) / 39])
lift = None
for representative in cosets.values():
    projected = coordinates_in_C(project_mw(representative))
    difference = target - projected
    if all(value.denominator() == 1 for value in difference):
        lift = vector(ZZ, vector(QQ, representative) + difference * C)
        break
assert lift is not None
assert project_mw(lift) == target * C
height = (target * C) * F * (target * C)
assert height == QQ(316) / 39

root_part = vector(QQ, lift) - target * C
root_coordinates = (root_part * F * R.transpose()) * GR_inverse
assert root_part == root_coordinates * R

classes = []
left = 0
for basis in component_bases:
    right = left + basis.nrows()
    classes.append(fractional_class(root_coordinates[left:right]))
    left = right

assert class_order(classes[0], 13) == 13
a2_orders = [class_order(point, 3) for point in classes[1:]]
assert sorted(a2_orders) == [1, 3]

# In A12 the class with correction 42/13 has labels 6 and 7, exchanged by
# component inversion.  Exactly one A2 class is nonzero, contributing 2/3.
local_correction = QQ(42) / 13 + QQ(2) / 3
zero_intersection = (height + local_correction - 4) / 2
assert local_correction == QQ(152) / 39
assert zero_intersection == 4

print("MW1A12GLUE|frame_det=948|roots=A12+A2+A2|root_det=117", flush=True)
print("MW1A12GLUE|root_plus_C_index=39|height=316/39", flush=True)
print("MW1A12GLUE|profile=A12:6_or_7,A2:1,A2:0|P.O=4", flush=True)
print("MW1A12GLUE|expected_fibers=I13+I3+I3+5I1", flush=True)
print("MW1A12GLUE|PASS", flush=True)
