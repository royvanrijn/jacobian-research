#!/usr/bin/env sage
"""Recover the component profile of the E8+D7+A1 MW1 frame exactly."""

from sage.all import *
from pathlib import Path


frame_path = Path("elkies-k3/data/fibrations/mw1_e8_d7_a1_frame.txt")
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
assert [len(component) for component in components] == [240, 84, 2]

component_bases = [
    matrix(ZZ, [roots[index] for index in component]).row_module().basis_matrix()
    for component in components
]
component_grams = [basis * F * basis.transpose() for basis in component_bases]
assert [
    (basis.rank(), abs(gram.det()))
    for basis, gram in zip(component_bases, component_grams)
] == [(8, 1), (7, 4), (1, 2)]

R = block_matrix([[basis] for basis in component_bases], subdivide=False)
GR = R * F * R.transpose()
C = (R * F).right_kernel_matrix()
GC = C * F * C.transpose()
assert R.rank() == 16 and C.rank() == 1 and GC == matrix(ZZ, [[474]])

A = block_matrix([[R], [C]], subdivide=False)
assert abs(A.det()) == 2
A_inverse = A.inverse()


def coset_key(point):
    return fractional_class(vector(QQ, point) * A_inverse)


zero = vector(ZZ, [0] * 17)
cosets = {coset_key(zero): zero}
for coordinate in range(17):
    unit = vector(ZZ, [0] * 17)
    unit[coordinate] = 1
    if coset_key(unit) not in cosets:
        cosets[coset_key(unit)] = unit
        break
assert len(cosets) == 2

GR_inverse = GR.inverse()
GC_inverse = GC.inverse()


def project_mw(point):
    point = vector(QQ, point)
    return point - (point * F * R.transpose()) * GR_inverse * R


def coordinates_in_C(point):
    return (vector(QQ, point) * F * C.transpose()) * GC_inverse


target = vector(QQ, [QQ(1) / 2])
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
assert height == QQ(237) / 2

root_part = vector(QQ, lift) - target * C
root_coordinates = (root_part * F * R.transpose()) * GR_inverse
assert root_part == root_coordinates * R

classes = []
left = 0
for basis in component_bases:
    right = left + basis.nrows()
    classes.append(fractional_class(root_coordinates[left:right]))
    left = right

assert all(value == 0 for value in classes[0])
assert any(value != 0 for value in classes[1])
assert any(value != 0 for value in classes[2])
assert all((2 * value).denominator() == 1 for value in classes[1])
assert all((2 * value).denominator() == 1 for value in classes[2])

# The unique order-two D7 discriminant class is the vector class, of
# correction 1.  The nonzero A1 class has correction 1/2.
local_correction = QQ(1) + QQ(1) / 2
zero_intersection = (height + local_correction - 4) / 2
assert zero_intersection == 58

print("MW1GLUE|frame_det=948|roots=E8+D7+A1|root_det=8", flush=True)
print("MW1GLUE|root_plus_C_index=2|height=237/2", flush=True)
print("MW1GLUE|profile=E8:0,D7:vector,A1:1|P.O=58", flush=True)
print("MW1GLUE|expected_fibers=II*+I3*+I2+3I1", flush=True)
print("MW1GLUE|PASS", flush=True)
