from sage.all import *
from pathlib import Path


frame_path = Path("elkies-k3/data/fibrations/mw2_e6_d4_a2a2_a1_frame.txt")
F = matrix(ZZ, [
    [ZZ(value) for value in line.split()]
    for line in frame_path.read_text().splitlines()
    if line.strip()
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


def class_multiple(multiplier, point):
    return fractional_class(tuple(multiplier * value for value in point))


def class_order(point, exponent):
    zero = tuple(QQ(0) for _ in point)
    return next(
        order for order in range(1, exponent + 1)
        if class_multiple(order, point) == zero
    )


# Split the complete norm-2 shell into E6, D4, A2, A2, and A1.
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
assert [len(component) for component in components] == [72, 24, 6, 6, 2]

component_bases = [
    matrix(ZZ, [roots[index] for index in component]).row_module().basis_matrix()
    for component in components
]
component_data = [
    (basis.rank(), abs((basis * F * basis.transpose()).det()))
    for basis in component_bases
]
assert component_data == [(6, 3), (4, 4), (2, 3), (2, 3), (1, 2)]

R = block_matrix([[basis] for basis in component_bases], subdivide=False)
GR = R * F * R.transpose()
C = (R * F).right_kernel_matrix()
GC = C * F * C.transpose()
assert R.rank() == 15 and C.rank() == 2

# Enumerate frame/(root + orthogonal-intersection), then saturate its
# projection to recover the exact Mordell--Weil lattice.
A = block_matrix([[R], [C]], subdivide=False)
index = abs(A.det())
assert index == 36
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


projected_cosets = [
    (coordinates_in_C(project_mw(representative)), representative)
    for representative in cosets.values()
]
generators = [vector(QQ, (1, 0)), vector(QQ, (0, 1))]
generators += [coordinates for coordinates, _ in projected_cosets]
denominator = lcm(
    QQ(value).denominator() for row in generators for value in row
)
MW_integer = matrix(ZZ, [
    [ZZ(denominator * value) for value in row]
    for row in generators
]).row_module().basis_matrix()
MW_basis = MW_integer.change_ring(QQ) / denominator
H = MW_basis * GC * MW_basis.transpose()
assert H == matrix(QQ, [[QQ(3) / 2, -QQ(10) / 3], [-QQ(10) / 3, QQ(31) / 3]])
assert H.det() == QQ(79) / 18

# This integral change of basis gives the small Gram used for reconstruction.
target_scaled = matrix(ZZ, [[9, 2], [2, 18]])
transform = matrix(ZZ, [[-1, 0], [2, 1]])
target_basis = transform * MW_basis
assert target_basis * GC * target_basis.transpose() == target_scaled / 6

# Lift the reduced generators back to the frame.  Their root corrections,
# modulo each irreducible root lattice, are the reducible-fiber component
# classes.  The two A2 factors are subsequently ordered lexicographically by
# the resulting pair of labels.
lifts = []
root_classes = []
component_bounds = []
start = 0
for basis in component_bases:
    component_bounds.append((start, start + basis.nrows()))
    start += basis.nrows()

for target_vector in target_basis.rows():
    lift = None
    for projected, representative in projected_cosets:
        difference = target_vector - projected
        if all(QQ(value).denominator() == 1 for value in difference):
            lift = vector(QQ, representative) + difference * C
            break
    assert lift is not None
    assert all(QQ(value).denominator() == 1 for value in lift)
    lift = vector(ZZ, lift)
    assert project_mw(lift) == target_vector * C
    lifts.append(lift)

    root_part = vector(QQ, lift) - target_vector * C
    root_coordinates = (root_part * F * R.transpose()) * GR_inverse
    assert root_part == root_coordinates * R
    root_classes.append([
        fractional_class(root_coordinates[left:right])
        for left, right in component_bounds
    ])


def cyclic_labels(component, modulus):
    classes = [root_classes[row][component] for row in range(2)]
    generator = next(point for point in classes if class_order(point, modulus) == modulus)
    return tuple(
        next(k for k in range(modulus) if class_multiple(k, generator) == point)
        for point in classes
    )


e6_labels = cyclic_labels(0, 3)
a2_labels = sorted((cyclic_labels(2, 3), cyclic_labels(3, 3)))
a1_labels = cyclic_labels(4, 2)

# D4 has three triality-equivalent nonzero classes.  Name the first class d1;
# a distinct second class is d2 (and their sum would be d3).
d4_classes = [root_classes[row][1] for row in range(2)]
zero_d4 = tuple(QQ(0) for _ in d4_classes[0])
d4_names = []
named_d4 = {}
for point in d4_classes:
    if point == zero_d4:
        d4_names.append("0")
    elif point in named_d4:
        d4_names.append(named_d4[point])
    else:
        name = f"d{len(named_d4) + 1}"
        named_d4[point] = name
        d4_names.append(name)
assert all(class_order(point, 2) in (1, 2) for point in d4_classes)

profiles = tuple(
    (
        e6_labels[row],
        d4_names[row],
        a2_labels[0][row],
        a2_labels[1][row],
        a1_labels[row],
    )
    for row in range(2)
)
assert profiles == ((1, "0", 0, 1, 1), (2, "d1", 1, 0, 0))


def e6_correction(left, right):
    if left == 0 or right == 0:
        return QQ(0)
    return QQ(4) / 3 if left == right else QQ(2) / 3


def d4_correction(left, right):
    if left == "0" or right == "0":
        return QQ(0)
    return QQ(1) if left == right else QQ(1) / 2


def a2_correction(left, right):
    if left == 0 or right == 0:
        return QQ(0)
    return QQ(2) / 3 if left == right else QQ(1) / 3


def a1_correction(left, right):
    return QQ(1) / 2 if left == right == 1 else QQ(0)


def local_correction(left, right):
    return (
        e6_correction(left[0], right[0])
        + d4_correction(left[1], right[1])
        + a2_correction(left[2], right[2])
        + a2_correction(left[3], right[3])
        + a1_correction(left[4], right[4])
    )


zero_intersections = tuple(
    (target_scaled[row, row] / 6 + local_correction(profiles[row], profiles[row]) - 4) / 2
    for row in range(2)
)
assert all(value in ZZ and value >= 0 for value in zero_intersections)
pair_intersection = (
    2 + sum(zero_intersections)
    - local_correction(profiles[0], profiles[1])
    - QQ(target_scaled[0, 1]) / 6
)
assert pair_intersection in ZZ and pair_intersection >= 0
assert zero_intersections == (0, 1)
assert pair_intersection == 2

print("MW2GLUE|frame_det=948|roots=E6+D4+A2+A2+A1|root_det=216", flush=True)
print("MW2GLUE|root_plus_C_index=36|mw_det=79/18", flush=True)
print("MW2GLUE|height_gram=(1/6)*[9,2;2,18]", flush=True)
print(
    "MW2GLUE|profiles="
    + ";".join(
        f"P{row + 1}:" + ",".join(map(str, profile))
        for row, profile in enumerate(profiles)
    ),
    flush=True,
)
print(
    "MW2GLUE|zero_intersections="
    + ",".join(map(str, zero_intersections))
    + f"|pair_intersection={pair_intersection}",
    flush=True,
)
print("MW2GLUE|PASS", flush=True)
