#!/usr/bin/env sage
"""Classify every q=8 neighbor orbit of the CM-43 Kumar frame.

Raw norm-16 enumeration has 1,421,331,656 signed vectors, almost all related
by W(E7)xW(E8).  This verifier instead decomposes a vector into its saturated
MW projection and its dominant E7/E8 Dynkin labels.  The norm equation and
the root-discriminant congruence then give a complete finite orbit list.

For every orbit representative it constructs the (a,b)=(2,4) U-neighbor,
computes exact root invariants, and tests any A7+A4+A3+A2 candidate against
the pinned semistable q=8 frame by an exact integral quadratic-form isometry.
"""

from pathlib import Path
import argparse

from sage.all import *


BASE = Path(__file__).resolve().parents[1]

parser = argparse.ArgumentParser()
parser.add_argument(
    "--frame",
    default="data/fibrations/kumar_cm43_marked_e7e8_mw3_frame.txt",
    help="CM43 positive frame, relative to elkies-k3 unless absolute",
)
args = parser.parse_args()


def load_gram(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ],
    )


frame_path = Path(args.frame)
if not frame_path.is_absolute():
    frame_path = BASE / frame_path
frame = load_gram(frame_path)
target = load_gram(
    BASE / "data/fibrations/picard20_mw2_a7_a4_a3_a2_frame.txt"
)
assert frame.det() == target.det() == 43
rank = frame.nrows()
norm16_signed_count = ZZ(pari(frame).qfminim(16, 1)[0])
assert norm16_signed_count == 1421331656

# Recover all roots and a deterministic positive system.  A positive root is
# simple precisely when it is not a sum of two positive roots.
root_result = pari(frame).qfminim(2)
assert ZZ(root_result[0]) == 366
half_roots = [vector(ZZ, column) for column in matrix(ZZ, root_result[2]).columns()]
roots = half_roots + [-root for root in half_roots]

regular = None
for shift in range(1, 100):
    candidate = vector(ZZ, [
        (index+1)**2 + shift*(index+1) + shift**2
        for index in range(rank)
    ])
    if all(root * frame * candidate != 0 for root in roots):
        regular = candidate
        break
assert regular is not None
positive_roots = [root for root in roots if root * frame * regular > 0]
positive_set = {tuple(root) for root in positive_roots}
simple_roots = []
for root in positive_roots:
    decomposable = any(
        tuple(root-left) in positive_set for left in positive_roots
    )
    if not decomposable:
        simple_roots.append(root)
simple = matrix(ZZ, [list(root) for root in simple_roots])
assert simple.nrows() == simple.rank() == 15
cartan = simple * frame * simple.transpose()
assert set(cartan.diagonal()) == {2}
assert all(cartan[row, column] in (0, -1) for row in range(15)
           for column in range(15) if row != column)

# Connected Dynkin components.
parents = list(range(15))


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


for row in range(15):
    for column in range(row):
        if cartan[row, column]:
            union(row, column)
component_indices = {}
for index in range(15):
    component_indices.setdefault(find(index), []).append(index)
components = sorted(component_indices.values(), key=len)
assert tuple(map(len, components)) == (7, 8)
assert tuple(abs(cartan.matrix_from_rows_and_columns(indices, indices).det())
             for indices in components) == (2, 1)


def root_torsion_order(root_basis):
    smith, _, _ = root_basis.smith_form()
    diagonal = [abs(ZZ(smith[index, index]))
                for index in range(root_basis.nrows())]
    return prod(value for value in diagonal if value)


def saturated_mw_data(gram, root_basis):
    """Return a saturated height Gram and integral lifts, as in the CM audit."""
    mw_rank = gram.nrows()-root_basis.rank()
    orthogonal = (root_basis*gram).right_kernel_matrix()
    root_gram = root_basis*gram*root_basis.transpose()
    orthogonal_gram = orthogonal*gram*orthogonal.transpose()
    combined = root_basis.stack(orthogonal)
    index = abs(ZZ(combined.det()))
    inverse_combined = combined.inverse()

    def fractional_key(row):
        coordinates = vector(QQ, row)*inverse_combined
        return tuple(value-value.floor() for value in coordinates)

    zero = vector(ZZ, [0]*gram.nrows())
    representatives = {fractional_key(zero): zero}
    queue = [zero]
    head = 0
    while head < len(queue) and len(representatives) < index:
        row = queue[head]
        head += 1
        for coordinate in range(gram.nrows()):
            unit = vector(ZZ, [1 if i == coordinate else 0
                               for i in range(gram.nrows())])
            for sign in (1, -1):
                candidate = row+sign*unit
                key = fractional_key(candidate)
                if key not in representatives:
                    representatives[key] = candidate
                    queue.append(candidate)
                    if len(representatives) == index:
                        break
            if len(representatives) == index:
                break
    assert len(representatives) == index

    inverse_root = root_gram.inverse()
    inverse_orthogonal = orthogonal_gram.inverse()

    def project(row):
        row = vector(QQ, row)
        root_coordinates = row*gram*root_basis.transpose()*inverse_root
        root_free = row-root_coordinates*root_basis
        return root_free*gram*orthogonal.transpose()*inverse_orthogonal

    representative_rows = list(representatives.values())
    projected = [project(row) for row in representative_rows]
    generators = list(Matrix.identity(QQ, mw_rank).rows())+projected
    denominator = lcm([vector(QQ, row).denominator() for row in generators])
    integer_rows = matrix(ZZ, [
        [ZZ(denominator*value) for value in row] for row in generators
    ])
    mw_basis = integer_rows.row_module().basis_matrix().change_ring(QQ)/denominator
    lifts = []
    for basis_row in mw_basis.rows():
        for projected_row, representative_row in zip(projected, representative_rows):
            difference = vector(QQ, basis_row)-vector(QQ, projected_row)
            if all(value.denominator() == 1 for value in difference):
                lift = vector(QQ, representative_row)+difference*orthogonal
                assert all(value.denominator() == 1 for value in lift)
                lifts.append(vector(ZZ, lift))
                break
        else:
            raise RuntimeError("failed to lift saturated MW vector")
    lifts = matrix(ZZ, [list(row) for row in lifts])
    height = mw_basis*orthogonal_gram*mw_basis.transpose()
    assert height.det() == QQ(gram.det()*root_torsion_order(root_basis)**2)/root_gram.det()

    scale = lcm([QQ(value).denominator() for value in height.list()])
    integral = (scale*height).change_ring(ZZ)
    lll = matrix(ZZ, pari(integral).qflllgram())
    if abs(lll.det()) == 1:
        height = lll*height*lll.transpose()
        lifts = lll*lifts
    return height, lifts


height, mw_lifts = saturated_mw_data(frame, simple)
assert height.det() == QQ(43)/2
assert (2*height).change_ring(ZZ).is_positive_definite()

marked_projection = None
if frame_path.name == "kumar_cm43_marked_e7e8_mw3_frame.txt":
    marked_height = matrix(QQ, (
        (QQ(5)/2, -QQ(1)/2, -1),
        (-QQ(1)/2, QQ(5)/2, 0),
        (-1, 0, 4),
    ))
    marked_height4 = vector(ZZ, (0, 0, 1))
    marked_q79 = vector(ZZ, (4, -5, 1))
    height4_frame = vector(ZZ, [0]*15+[4, 0, -1])
    q79_frame = vector(ZZ, [0]*16+[1, 0])
    root_gram_inverse = (simple*frame*simple.transpose()).inverse()

    def internal_mw_coordinates(row):
        row = vector(QQ, row)
        root_coordinates = row*frame*simple.transpose()*root_gram_inverse
        projection = row-root_coordinates*simple
        pairings = projection*frame*mw_lifts.transpose()
        coordinates = pairings*height.inverse()
        assert all(value in ZZ for value in coordinates)
        return vector(ZZ, coordinates)

    marking_isometry = matrix(
        ZZ, pari(2*height).qfisom(pari(2*marked_height))
    )
    # PARI qfisom(A,B) returns C with C^t*B*C=A.
    assert (
        marking_isometry.transpose()*(2*marked_height)*marking_isometry
        == 2*height
    )

    def map_internal_to_marked(coordinates):
        return vector(ZZ, coordinates*marking_isometry.transpose())

    mapped_height4 = map_internal_to_marked(
        internal_mw_coordinates(height4_frame)
    )
    mapped_q79 = map_internal_to_marked(internal_mw_coordinates(q79_frame))
    if mapped_height4 == -marked_height4 and mapped_q79 == -marked_q79:
        marking_isometry = -marking_isometry
        mapped_height4 = -mapped_height4
        mapped_q79 = -mapped_q79
    assert mapped_height4 == marked_height4 and mapped_q79 == marked_q79

    def marked_projection(row):
        return tuple(map_internal_to_marked(internal_mw_coordinates(row)))

# Enumerate all MW vectors of height at most 16.
height_scale = lcm([QQ(value).denominator() for value in height.list()])
scaled_height = (height_scale*height).change_ring(ZZ)
mw_result = pari(scaled_height).qfminim(height_scale*16)
mw_vectors = {tuple([0]*height.nrows())}
for column in matrix(ZZ, mw_result[2]).columns():
    mw_vectors.add(tuple(column))
    mw_vectors.add(tuple(-column))


def dominant_weights(component, bound):
    component_cart = cartan.matrix_from_rows_and_columns(component, component)
    inverse = component_cart.inverse()
    size = len(component)
    weights = []

    def recurse(prefix, norm):
        index = len(prefix)
        if index == size:
            weights.append((tuple(prefix), norm))
            return
        value = 0
        while True:
            added = inverse[index, index]*value**2
            added += 2*value*sum(
                inverse[index, previous]*prefix[previous]
                for previous in range(index)
            )
            new_norm = norm+added
            if new_norm > bound:
                break
            recurse(prefix+[value], new_norm)
            value += 1

    recurse([], QQ(0))
    return weights


component_weights = [dominant_weights(component, QQ(16))
                     for component in components]
cartan_inverse = cartan.inverse()
dominant_vectors = {}
for mw_tuple in sorted(mw_vectors):
    mw = vector(ZZ, mw_tuple)
    mw_norm = mw*height*mw
    if mw_norm > 16:
        continue
    lift = mw*mw_lifts
    lift_labels = lift*frame*simple.transpose()
    for labels7, norm7 in component_weights[0]:
        for labels8, norm8 in component_weights[1]:
            if mw_norm+norm7+norm8 != 16:
                continue
            labels = vector(ZZ, [0]*15)
            for index, value in zip(components[0], labels7):
                labels[index] = value
            for index, value in zip(components[1], labels8):
                labels[index] = value
            root_coordinates = (labels-lift_labels)*cartan_inverse
            if not all(value in ZZ for value in root_coordinates):
                continue
            vector_value = lift+vector(ZZ, root_coordinates)*simple
            assert vector_value*frame*vector_value == 16
            assert vector_value*frame*simple.transpose() == labels
            dominant_vectors[tuple(vector_value)] = (
                mw_tuple, tuple(labels), mw_norm, norm7+norm8
            )


U = matrix(ZZ, ((0, 1), (1, 0)))
NS = block_diagonal_matrix(U, -frame)
q60_fiber = vector(ZZ, (
    5, 12,
    0, 0, -1, -1, -1, -1, -1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0,
))
assert q60_fiber*NS*q60_fiber == 0


def bezout_vector_for_pairing(fiber):
    pairings = list(NS*fiber)
    current = ZZ(0)
    coefficients = [ZZ(0)]*NS.nrows()
    for index, value in enumerate(pairings):
        if value == 0:
            continue
        new_gcd, old_scale, new_scale = xgcd(current, ZZ(value))
        coefficients = [old_scale*entry for entry in coefficients]
        coefficients[index] += new_scale
        current = new_gcd
    if abs(current) != 1:
        return None
    if current == -1:
        coefficients = [-entry for entry in coefficients]
    return vector(ZZ, coefficients)


def child_frame(witness):
    fiber = vector(ZZ, [2, 4]+list(witness))
    assert fiber*NS*fiber == 0
    mate = bezout_vector_for_pairing(fiber)
    if mate is None:
        return None
    mate_square = ZZ(mate*NS*mate)
    assert mate_square % 2 == 0
    mate -= (mate_square//2)*fiber
    assert fiber*NS*mate == 1 and mate*NS*mate == 0
    kernel = matrix(ZZ, [list(fiber*NS), list(mate*NS)]).right_kernel_matrix()
    child = -(kernel*NS*kernel.transpose())
    assert child.is_positive_definite() and child.det() == 43
    neighbor_basis = matrix(
        ZZ, [list(fiber), list(mate)] + [list(row) for row in kernel.rows()]
    )
    assert abs(neighbor_basis.det()) == 1
    assert neighbor_basis*NS*neighbor_basis.transpose() == block_diagonal_matrix(
        U, -child
    )
    return child, neighbor_basis


def root_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    if count == 0:
        return (0, 0, 1)
    root_rows = matrix(ZZ, result[2]).transpose()
    root_basis = root_rows.row_module().basis_matrix()
    root_gram = root_basis*gram*root_basis.transpose()
    return (root_basis.rank(), count, abs(ZZ(root_gram.det())))


histogram = {}
target_hits = []
factorizations = []
for witness_tuple, decomposition in sorted(dominant_vectors.items()):
    child_data = child_frame(vector(ZZ, witness_tuple))
    if child_data is None:
        continue
    child, neighbor_basis = child_data
    invariants = root_data(child)
    histogram[invariants] = histogram.get(invariants, 0)+1
    q60_in_child = q60_fiber*neighbor_basis.inverse()
    assert all(value in ZZ for value in q60_in_child)
    q60_in_child = vector(ZZ, q60_in_child)
    child_ns = block_diagonal_matrix(U, -child)
    assert q60_in_child*child_ns*q60_in_child == 0
    if q60_in_child[0] < 0 and q60_in_child[1] < 0:
        q60_in_child = -q60_in_child
    second_q = ZZ(q60_in_child[0]*q60_in_child[1])
    assert second_q >= 0
    factorizations.append((
        max(q60_in_child[0], q60_in_child[1]),
        second_q,
        q60_in_child[1],
        witness_tuple,
        invariants,
        tuple(q60_in_child),
    ))
    if invariants != (16, 94, 480):
        continue
    isometry = QuadraticForm(ZZ, child).is_globally_equivalent_to(
        QuadraticForm(ZZ, target), return_matrix=True
    )
    assert isometry is not False
    assert abs(isometry.det()) == 1
    assert isometry.transpose()*child*isometry == target
    target_hits.append((witness_tuple, decomposition))
    print(
        f"KUMARCM43Q8|target={len(target_hits)}|witness={witness_tuple}"
        f"|decomposition={decomposition}|isometric_to_semistable_q8=1",
        flush=True,
    )

factorizations.sort()
for index, record in enumerate(factorizations[:12], 1):
    maximum, second_q, degree, witness_tuple, invariants, coordinates = record
    print(
        f"KUMARCM43Q8FACTOR|rank={index}|first_q=8|second_q={second_q}"
        f"|second_ab={coordinates[0]},{coordinates[1]}|degree={degree}"
        f"|intermediate_roots={invariants}|first_witness={witness_tuple}"
        f"|first_mw_projection={dominant_vectors[witness_tuple][0]}"
        f"|first_marked_projection="
        f"{None if marked_projection is None else marked_projection(witness_tuple)}"
        f"|q60_in_child={coordinates}",
        flush=True,
    )

# The source-level Humbert-8 construction has a D9+E7 pre-neighbor.  Its
# root invariants are (rank,count,det)=(16,270,8).  Record the cheapest q=8
# orbit of that exact root type separately: it need not be the orbit giving
# the arithmetically cheapest factorization of the q=60 automorphism.
d9e7_factorizations = [
    record for record in factorizations if record[4] == (16, 270, 8)
]
assert len(d9e7_factorizations) == histogram.get((16, 270, 8), 0) == 10
for index, record in enumerate(d9e7_factorizations, 1):
    _, source_second_q, _, source_witness, _, source_coordinates = record
    print(
        f"KUMARCM43D9E7ORBIT|rank={index}|second_q={source_second_q}"
        f"|witness={source_witness}"
        f"|mw_projection={dominant_vectors[source_witness][0]}"
        f"|marked_projection="
        f"{None if marked_projection is None else marked_projection(source_witness)}"
        f"|q60_in_child={source_coordinates}",
        flush=True,
    )
d9e7_best = d9e7_factorizations[0]
_, d9e7_second_q, _, d9e7_witness, _, d9e7_coordinates = d9e7_best
print(
    f"KUMARCM43D9E7|first_q=8|orbits={len(d9e7_factorizations)}"
    f"|best_second_q={d9e7_second_q}"
    f"|second_ab={d9e7_coordinates[0]},{d9e7_coordinates[1]}"
    f"|first_witness={d9e7_witness}"
    f"|q60_in_child={d9e7_coordinates}|status=PASS",
    flush=True,
)

print(
    f"KUMARCM43Q8|frame={frame_path.name}"
    f"|raw_signed_shell={norm16_signed_count}"
    f"|mw_vectors={len(mw_vectors)}|dominant_orbits={len(dominant_vectors)}"
    f"|primitive_neighbors={sum(histogram.values())}"
    f"|root_histogram={tuple(sorted(histogram.items()))}"
    f"|target_hits={len(target_hits)}|best_factor_second_q={factorizations[0][1]}"
    f"|d9e7_orbits={len(d9e7_factorizations)}"
    f"|d9e7_best_second_q={d9e7_second_q}"
    f"|status=PASS",
    flush=True,
)
