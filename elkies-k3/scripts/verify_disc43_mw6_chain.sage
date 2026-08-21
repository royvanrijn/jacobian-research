#!/usr/bin/env sage
"""Replay the discriminant-43 MW3 -> MW5 -> MW5 -> MW6 chain."""

from sage.all import *
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "elkies-k3" / "data" / "fibrations"
U = matrix(ZZ, [[0, 1], [1, 0]])


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def bezout_vector(pairings):
    current = ZZ(0)
    coefficients = [ZZ(0)] * len(pairings)
    for index, pairing in enumerate(pairings):
        if pairing == 0:
            continue
        new_gcd, left, right = xgcd(current, ZZ(pairing))
        coefficients = [left * value for value in coefficients]
        coefficients[index] += right
        current = new_gcd
    assert abs(current) == 1
    return vector(ZZ, coefficients if current == 1 else [-x for x in coefficients])


def neighbor(parent, qnorm, a, b, coordinates):
    ns = block_diagonal_matrix(U, -parent)
    fiber = vector(ZZ, [a, b] + list(coordinates))
    assert a * b == qnorm
    assert coordinates * parent * coordinates == 2 * qnorm
    assert fiber * ns * fiber == 0
    assert gcd([abs(ZZ(value)) for value in ns * fiber]) == 1
    mate = bezout_vector(list(ns * fiber))
    assert fiber * ns * mate == 1
    mate -= ZZ(mate * ns * mate) // 2 * fiber
    assert mate * ns * mate == 0 and fiber * ns * mate == 1
    complement = matrix(
        ZZ, [list(fiber * ns), list(mate * ns)]
    ).right_kernel_matrix()
    child = -(complement * ns * complement.transpose())
    transport = matrix(ZZ, [list(fiber), list(mate)] + complement.rows())
    assert abs(transport.det()) == 1
    assert transport * ns * transport.transpose() == block_diagonal_matrix(U, -child)
    return child, transport


canonical = load_matrix(DATA / "picard20_e6_d4_a2a2_a1_mw3_frame.txt")
start = load_matrix(DATA / "disc43_mw3_frame.txt")
target = load_matrix(DATA / "disc43_mw6_a4a4a4_frame.txt")
assert canonical.det() == start.det() == target.det() == 43

# The two pinned MW3 frames use different integral complements of U.  PARI
# supplies an exact unimodular isometry, so the chain starts on the verified
# explicit Picard-20 frame rather than only an abstract genus mate.
initial_isometry = matrix(ZZ, pari(start).qfisom(pari(canonical)))
assert abs(initial_isometry.det()) == 1
assert initial_isometry.transpose() * canonical * initial_isometry == start

steps = (
    (4, 2, 2, (0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1, 0, -1, -1, 0, 0, 0)),
    (4, 1, 4, (-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
    (12, 3, 4, (0, -2, 0, 0, 1, -1, 0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0)),
)

frame = start
transports = []
expected_data = ((13, 78, 96), (13, 78, 96), (12, 60, 125))
for index, ((qnorm, a, b, values), expected) in enumerate(zip(steps, expected_data), 1):
    frame, transport = neighbor(frame, qnorm, a, b, vector(ZZ, values))
    transports.append(transport)
    roots = pari(frame).qfminim(2)
    count = ZZ(roots[0])
    root_vectors = matrix(ZZ, roots[2]).transpose()
    root_basis = root_vectors.row_module().basis_matrix()
    root_det = abs((root_basis * frame * root_basis.transpose()).det())
    actual = (root_basis.rank(), count, root_det)
    assert actual == expected
    print(
        f"DISC43MW6|step={index}|q={qnorm}|ab={a},{b}"
        f"|root_rank={actual[0]}|roots={actual[1]}|rootdet={actual[2]}"
        f"|MW={18-actual[0]}|transport_det={transport.det()}", flush=True,
    )
assert frame == target

# Identify the terminal root system exactly as A4^3.
root_result = pari(target).qfminim(2)
half_roots = [vector(ZZ, column) for column in matrix(ZZ, root_result[2]).columns()]
signed_roots = half_roots + [-root for root in half_roots]
graph = Graph()
graph.add_vertices(range(len(signed_roots)))
for left in range(len(signed_roots)):
    for right in range(left):
        if signed_roots[left] * target * signed_roots[right] != 0:
            graph.add_edge(left, right)
components = sorted(graph.connected_components(sort=False), key=len, reverse=True)
assert [len(component) for component in components] == [20, 20, 20]
component_bases = [
    matrix(ZZ, [signed_roots[index] for index in component]).row_module().basis_matrix()
    for component in components
]
assert all(basis.rank() == 4 for basis in component_bases)
assert all(abs((basis * target * basis.transpose()).det()) == 5 for basis in component_bases)

# Saturate the projection of target/(root + orthogonal intersection) to
# recover the exact rank-six Mordell-Weil height lattice and torsion.
R = block_matrix([[basis] for basis in component_bases], subdivide=False)
GR = R * target * R.transpose()
C = (R * target).right_kernel_matrix()
GC = C * target * C.transpose()
A = block_matrix([[R], [C]], subdivide=False)
index = abs(A.det())
assert index == 125
A_inverse = A.inverse()


def fractional_class(row):
    return tuple(QQ(value) - floor(QQ(value)) for value in row)


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
assert len(cosets) == index

GR_inverse = GR.inverse()
GC_inverse = GC.inverse()


def project_mw(point):
    point = vector(QQ, point)
    return point - (point * target * R.transpose()) * GR_inverse * R


def coordinates_in_C(point):
    point = vector(QQ, point)
    return (point * target * C.transpose()) * GC_inverse


projected = [coordinates_in_C(project_mw(point)) for point in cosets.values()]
assert sum(point == 0 for point in projected) == 1  # no torsion
generators = [
    vector(QQ, [ZZ(i == j) for i in range(6)]) for j in range(6)
] + projected
denominator = lcm(QQ(value).denominator() for row in generators for value in row)
integer_basis = matrix(ZZ, [
    [ZZ(denominator * value) for value in row] for row in generators
]).row_module().basis_matrix()
mw_basis = integer_basis.change_ring(QQ) / denominator
height = mw_basis * GC * mw_basis.transpose()
assert denominator == 5 and height.det() == QQ(43) / 125

scaled_height = (5 * height).change_ring(ZZ)
lll_transform = matrix(ZZ, pari(scaled_height).qflllgram(1))
reduced_height = lll_transform.transpose() * scaled_height * lll_transform
reported_height = matrix(ZZ, [
    [4, -2, 1, 1, -2, 2],
    [-2, 4, 1, -2, 0, -1],
    [1, 1, 6, -3, -1, 3],
    [1, -2, -3, 6, 0, -2],
    [-2, 0, -1, 0, 6, -2],
    [2, -1, 3, -2, -2, 8],
])
# LLL bases are not canonical across PARI versions.  Certify the pinned Gram
# by an exact integral isometry instead of requiring one particular LLL output.
height_isometry = matrix(ZZ, pari(reduced_height).qfisom(pari(reported_height)))
assert abs(height_isometry.det()) == 1
assert height_isometry.transpose() * reported_height * height_isometry == reduced_height

print("DISC43MW6|fibers=3I5+9I1|roots=3A4|torsion=0", flush=True)
print("DISC43MW6|height_gram=(1/5)*[4,-2,1,1,-2,2;-2,4,1,-2,0,-1;1,1,6,-3,-1,3;1,-2,-3,6,0,-2;-2,0,-1,0,6,-2;2,-1,3,-2,-2,8]", flush=True)
print("DISC43MW6|height_det=43/125|status=PASS", flush=True)
