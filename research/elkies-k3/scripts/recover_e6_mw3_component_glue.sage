from sage.all import *
from pathlib import Path


frame_path = Path("elkies-k3/data/fibrations/mw3_e6_a3a3_a1a1_frame.txt")
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
            coefficients.append(gram[i, i]//2 if i == j else gram[i, j])
    return QuadraticForm(ZZ, gram.nrows(), coefficients)


# Split the exact norm-2 shell into its irreducible root components.
half_roots = [
    vector(ZZ, root)
    for root in qform_from_gram(F).short_vector_list_up_to_length(2, True)[1]
]
roots = half_roots+[-root for root in half_roots]
graph = Graph()
graph.add_vertices(range(len(roots)))
for i in range(len(roots)):
    for j in range(i):
        if roots[i]*F*roots[j] != 0:
            graph.add_edge(i, j)
components = sorted(graph.connected_components(sort=False), key=len, reverse=True)
assert [len(component) for component in components] == [72, 12, 12, 2, 2]

component_bases = [
    matrix(ZZ, [roots[index] for index in component]).row_module().basis_matrix()
    for component in components
]
component_data = [
    (basis.rank(), abs((basis*F*basis.transpose()).det()))
    for basis in component_bases
]
assert component_data == [(6, 3), (3, 4), (3, 4), (1, 2), (1, 2)]

R = block_matrix([[basis] for basis in component_bases], subdivide=False)
GR = R*F*R.transpose()
C = (R*F).right_kernel_matrix()
GC = C*F*C.transpose()
assert R.rank() == 14 and C.rank() == 3

# Enumerate the exact frame/(root + orthogonal-intersection) glue and project
# every coset to the Mordell-Weil space.
A = block_matrix([[R], [C]], subdivide=False)
index = abs(A.det())
assert index == 96
A_inverse = A.inverse()


def coset_key(point):
    coordinates = vector(QQ, point)*A_inverse
    return tuple(value-floor(value) for value in coordinates)


zero = vector(ZZ, [0]*17)
cosets = {coset_key(zero): zero}
queue = [zero]
head = 0
while head < len(queue) and len(cosets) < index:
    point = queue[head]
    head += 1
    for coordinate in range(17):
        unit = vector(ZZ, [0]*17)
        unit[coordinate] = 1
        for sign in (1, -1):
            candidate = point+sign*unit
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
    return point-(point*F*R.transpose())*GR_inverse*R


def coordinates_in_C(point):
    return (vector(QQ, point)*F*C.transpose())*GC_inverse


projected_cosets = [
    (coordinates_in_C(project_mw(representative)), representative)
    for representative in cosets.values()
]
generators = [
    vector(QQ, (1, 0, 0)),
    vector(QQ, (0, 1, 0)),
    vector(QQ, (0, 0, 1)),
] + [coordinates for coordinates, _ in projected_cosets]
denominator = lcm(
    QQ(value).denominator() for row in generators for value in row
)
MW_integer = matrix(ZZ, [
    [ZZ(denominator*value) for value in row]
    for row in generators
]).row_module().basis_matrix()
MW_basis = MW_integer.change_ring(QQ)/denominator
H = MW_basis*GC*MW_basis.transpose()
assert H.det() == QQ(79)/16

# Recover the declared reduced basis without relying on an LLL convention.
target_scaled = matrix(ZZ, [
    [23, -10, -8],
    [-10, 23, 1],
    [-8, 1, 23],
])
scaled = (12*H).change_ring(ZZ)
minimum = pari(scaled).qfminim(23)
assert ZZ(minimum[0]) == 6 and ZZ(minimum[1]) == 23
representatives = list(matrix(ZZ, minimum[2]).columns())
minimal_vectors = representatives+[-point for point in representatives]
transforms = []
for first in minimal_vectors:
    for second in minimal_vectors:
        for third in minimal_vectors:
            transform = matrix(ZZ, [first, second, third])
            if transform*scaled*transform.transpose() == target_scaled:
                transforms.append(transform)
assert len(transforms) == 2 and transforms[0] == -transforms[1]
transform = next(item for item in transforms if item.det() == -1)
target_basis = transform*MW_basis
assert target_basis*GC*target_basis.transpose() == target_scaled/12

# Lift each reduced MW generator back to the frame.  The root correction of a
# lift, modulo the root lattice, is exactly its reducible-fiber component class.
lifts = []
root_classes = []
component_bounds = []
start = 0
for basis in component_bases:
    component_bounds.append((start, start+basis.nrows()))
    start += basis.nrows()

for target_vector in target_basis.rows():
    lift = None
    for projected, representative in projected_cosets:
        difference = target_vector-projected
        if all(QQ(value).denominator() == 1 for value in difference):
            lift = vector(QQ, representative)+difference*C
            break
    assert lift is not None
    assert all(QQ(value).denominator() == 1 for value in lift)
    lift = vector(ZZ, lift)
    assert project_mw(lift) == target_vector*C
    lifts.append(lift)

    root_part = vector(QQ, lift)-target_vector*C
    root_coordinates = (root_part*F*R.transpose())*GR_inverse
    assert root_part == root_coordinates*R
    root_classes.append([
        tuple(value-floor(value) for value in root_coordinates[left:right])
        for left, right in component_bounds
    ])


def class_multiple(multiplier, point):
    return tuple(
        multiplier*value-floor(multiplier*value)
        for value in point
    )


def class_order(point, modulus):
    zero_class = tuple(QQ(0) for _ in point)
    return next(
        order for order in range(1, modulus+1)
        if class_multiple(order, point) == zero_class
    )


# Put the A3 component on which P1 is trivial first.  The two A1 classes are
# identical on this basis, so their order is immaterial.
component_order = [0]
component_order += sorted(
    (1, 2), key=lambda index: class_order(root_classes[0][index], 4)
)
component_order += [3, 4]
moduli = (3, 4, 4, 2, 2)
profiles = [[] for _ in range(3)]
for output_index, component_index in enumerate(component_order):
    modulus = moduli[output_index]
    classes = [root_classes[row][component_index] for row in range(3)]
    generator = next(
        point for point in classes if class_order(point, modulus) == modulus
    )
    zero_class = tuple(QQ(0) for _ in generator)
    for row, point in enumerate(classes):
        label = next(
            multiplier for multiplier in range(modulus)
            if class_multiple(multiplier, generator) == point
        )
        profiles[row].append(label)
        assert class_multiple(modulus, point) == zero_class

profiles = tuple(tuple(profile) for profile in profiles)
assert profiles == (
    (1, 0, 1, 0, 0),
    (1, 1, 2, 1, 1),
    (2, 3, 0, 0, 0),
)

# The reduced height lattice has only +/- identity as automorphisms.  Thus the
# second 16-element Shioda-compatible label orbit is not a basis change of the
# actual frame glue.
target_minimum = pari(target_scaled).qfminim(23)
target_representatives = list(matrix(ZZ, target_minimum[2]).columns())
target_vectors = target_representatives+[-point for point in target_representatives]
automorphisms = []
for first in target_vectors:
    for second in target_vectors:
        for third in target_vectors:
            candidate = matrix(ZZ, [first, second, third])
            if candidate*target_scaled*candidate.transpose() == target_scaled:
                automorphisms.append(candidate)
assert automorphisms == [identity_matrix(ZZ, 3), -identity_matrix(ZZ, 3)] or automorphisms == [-identity_matrix(ZZ, 3), identity_matrix(ZZ, 3)]

print("E6GLUE|frame_det=948|roots=E6+A3+A3+A1+A1|root_det=192", flush=True)
print("E6GLUE|root_plus_C_index=96|mw_det=79/16", flush=True)
print("E6GLUE|height_gram=(1/12)*[23,-10,-8;-10,23,1;-8,1,23]", flush=True)
print(
    "E6GLUE|profiles="
    + ";".join(
        f"P{index+1}:"+",".join(map(str, profile))
        for index, profile in enumerate(profiles)
    ),
    flush=True,
)
print("E6GLUE|height_automorphisms=2|generators=+I,-I", flush=True)
print("E6GLUE|component_orbit=canonical_first_16|second_orbit=EXCLUDED_BY_FRAME_GLUE", flush=True)
print("E6GLUE|PASS", flush=True)
