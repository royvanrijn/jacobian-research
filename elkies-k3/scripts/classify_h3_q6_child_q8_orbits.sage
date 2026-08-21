#!/usr/bin/env sage
"""Classify the complete horizontal q=8 shell after the first H3 q=6 move.

The raw norm-16 shell in the E8+E6/MW3 child is far too large to enumerate
vector by vector.  This checker reconstructs that child from the labeled H3
frame, splits off its primitive E8+E6 root lattice, and works modulo its Weyl
group.  A horizontal vector is determined by:

* one of the finitely many nonzero vectors of projected MW norm at most 16;
* nonnegative dominant E6 and E8 Dynkin labels; and
* the exact root-discriminant integrality congruence.

This gives a complete orbit list.  The checker constructs every primitive
(a,b)=(4,2) neighbor, classifies its roots, and records the two D13/MW4
orbits together with all changes of basis needed to transport them back to
the original H3 Neron--Severi coordinates.
"""

import argparse
import json
from collections import Counter
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "elkies-k3"
SOURCE_FRAME = BASE / "data/fibrations/kumar_e7e8_mw2_frame_3.txt"
PINNED_D13_FRAME = (
    BASE / "data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt"
)

parser = argparse.ArgumentParser()
parser.add_argument(
    "--output",
    type=Path,
    default=ROOT / "artifacts/generated-results/elkies-k3-h3-q6-q8-orbits.json",
)
args = parser.parse_args()

U = matrix(ZZ, ((0, 1), (1, 0)))


def load_gram(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ],
    )


def bezout_vector_for_pairing(ns, fiber):
    pairings = list(ns * fiber)
    current = ZZ(0)
    result = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(pairings):
        if value == 0:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        result = [left * entry for entry in result]
        result[index] += right
        current = divisor
    if abs(current) != 1:
        return None
    if current == -1:
        result = [-entry for entry in result]
    return vector(ZZ, result)


def child_frame(ns, fiber, determinant):
    assert fiber * ns * fiber == 0
    mate = bezout_vector_for_pairing(ns, fiber)
    if mate is None:
        return None
    mate_square = ZZ(mate * ns * mate)
    assert mate_square % 2 == 0
    mate -= (mate_square // 2) * fiber
    assert mate * ns * mate == 0 and fiber * ns * mate == 1
    kernel = matrix(
        ZZ, [list(fiber * ns), list(mate * ns)]
    ).right_kernel_matrix()
    child = -(kernel * ns * kernel.transpose())
    assert child.is_positive_definite() and child.det() == determinant
    neighbor_basis = matrix(
        ZZ,
        [list(fiber), list(mate)] + [list(row) for row in kernel.rows()],
    )
    assert abs(neighbor_basis.det()) == 1
    assert (
        neighbor_basis * ns * neighbor_basis.transpose()
        == block_diagonal_matrix(U, -child)
    )
    return child, neighbor_basis


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    if count == 0:
        return (), matrix(ZZ, 0, gram.nrows()), (0, 0, 1)
    half = [
        vector(ZZ, column)
        for column in matrix(ZZ, result[2]).columns()
    ]
    roots = tuple(half + [-root for root in half])
    root_basis = matrix(ZZ, [list(root) for root in roots]).row_module().basis_matrix()
    root_gram = root_basis * gram * root_basis.transpose()
    return roots, root_basis, (
        root_basis.rank(),
        count,
        abs(ZZ(root_gram.det())),
    )


def deterministic_simple_roots(gram):
    roots, _, data = roots_and_data(gram)
    root_rank = data[0]
    regular = None
    for shift in range(1, 1000):
        candidate = vector(
            ZZ,
            [
                (index + 1) ** 2 + shift * (index + 1) + 1
                for index in range(gram.nrows())
            ],
        )
        if all(candidate * root != 0 for root in roots):
            regular = candidate
            break
    assert regular is not None
    positive = [root for root in roots if regular * root > 0]
    positive_set = {tuple(root) for root in positive}
    simple = []
    for root in positive:
        if not any(tuple(root - left) in positive_set for left in positive):
            simple.append(root)
    simple = matrix(ZZ, [list(root) for root in simple])
    assert simple.nrows() == simple.rank() == root_rank
    cartan = simple * gram * simple.transpose()
    assert set(cartan.diagonal()) == {2}
    assert all(
        cartan[row, column] in (0, -1)
        for row in range(root_rank)
        for column in range(root_rank)
        if row != column
    )
    return simple, cartan


def connected_components(cartan):
    unseen = set(range(cartan.nrows()))
    result = []
    while unseen:
        first = min(unseen)
        unseen.remove(first)
        todo = [first]
        component = []
        while todo:
            index = todo.pop()
            component.append(index)
            adjacent = [other for other in unseen if cartan[index, other]]
            for other in adjacent:
                unseen.remove(other)
                todo.append(other)
        result.append(tuple(sorted(component)))
    return tuple(sorted(result, key=lambda component: (len(component), component)))


def dominant_weights(cartan, component, bound):
    component_cartan = cartan.matrix_from_rows_and_columns(component, component)
    inverse = component_cartan.inverse()
    assert all(value >= 0 for value in inverse.list())
    weights = []

    def recurse(prefix, norm):
        index = len(prefix)
        if index == len(component):
            weights.append((tuple(prefix), norm))
            return
        value = 0
        while True:
            added = inverse[index, index] * value**2
            added += 2 * value * sum(
                inverse[index, previous] * prefix[previous]
                for previous in range(index)
            )
            new_norm = norm + added
            if new_norm > bound:
                break
            recurse(prefix + [value], new_norm)
            value += 1

    recurse([], QQ(0))
    return tuple(weights)


def exact_norm_count(gram, norm):
    at_norm = ZZ(pari(gram).qfminim(norm, 1)[0])
    below_norm = ZZ(pari(gram).qfminim(norm - 1, 1)[0])
    return at_norm - below_norm


def d13_root_adaptation(child):
    """Return a simple-root/LLL-quotient basis for a D13/MW4 frame."""
    _, root_basis, invariants = roots_and_data(child)
    assert invariants == (13, 312, 4)
    simple, cartan = deterministic_simple_roots(child)
    assert connected_components(cartan) == (tuple(range(13)),)
    assert cartan.det() == 4

    smith, smith_left, smith_right = root_basis.smith_form()
    assert smith == smith_left * root_basis * smith_right
    assert tuple(abs(smith[index, index]) for index in range(13)) == (1,) * 13
    completion = smith_right.inverse()
    assert abs(root_basis.stack(completion[13:]).det()) == 1
    adapted_basis = simple.stack(completion[13:])
    assert abs(adapted_basis.det()) == 1
    adapted = adapted_basis * child * adapted_basis.transpose()
    coupling = adapted[:13, 13:]
    tail = adapted[13:, 13:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling
    scale = lcm(entry.denominator() for entry in height.list())
    lll = matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram())
    assert abs(lll.det()) == 1
    # PARI returns a column change T with T^t*G*T reduced.  Frame basis
    # vectors are rows, so the corresponding row change is T^t.
    quotient_change = block_diagonal_matrix(
        identity_matrix(ZZ, 13), lll.transpose()
    )
    adapted_basis = quotient_change * adapted_basis
    adapted = adapted_basis * child * adapted_basis.transpose()
    assert adapted[:13, :13] == cartan
    coupling = adapted[:13, 13:]
    tail = adapted[13:, 13:]
    height = tail - coupling.transpose() * cartan.inverse() * coupling
    return simple, adapted_basis, adapted, height


def matrix_rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def rational_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


# Reconstruct the first q=6 child directly from the labeled H3 frame.
source = load_gram(SOURCE_FRAME)
assert source.nrows() == 17 and source.det() == 948
source_ns = block_diagonal_matrix(U, -source)
q6_witness = vector(
    ZZ,
    (0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
)
q6_fiber = vector(ZZ, [3, 2] + list(q6_witness))
q6_data = child_frame(source_ns, q6_fiber, determinant=948)
assert q6_data is not None
q6_child, q6_neighbor_basis = q6_data
q6_roots, q6_root_basis, q6_root_data = roots_and_data(q6_child)
assert q6_root_data == (14, 312, 3)

# These lifts are a saturated MW basis in the deterministic q=6 child
# coordinates above.  Stacking them below the primitive root basis is
# unimodular, so this also pins the exact root/MW coordinate convention.
mw_lifts = matrix(
    ZZ,
    [
        [-5, -4, -3, 0, 0, 0, 0, 0, 0, 0, 0, -4, 1, 0, -4, 2, -2],
        [-10, -8, -6, 0, 0, 0, 0, 0, 0, 0, 0, -8, 4, 1, -8, 5, -4],
        [-5, -4, -3, 0, 0, 0, 0, 0, 0, 0, 0, -3, 2, 0, -4, 2, -2],
    ],
)
root_mw_basis = q6_root_basis.stack(mw_lifts)
assert abs(root_mw_basis.det()) == 1
root_adapted = root_mw_basis * q6_child * root_mw_basis.transpose()
root_block = root_adapted[:14, :14]
root_coupling = root_adapted[:14, 14:]
mw_block = root_adapted[14:, 14:]
height = (
    mw_block
    - root_coupling.transpose() * root_block.inverse() * root_coupling
)
expected_height = matrix(
    QQ,
    (
        (QQ(8) / 3, QQ(1) / 3, -1),
        (QQ(1) / 3, QQ(8) / 3, 1),
        (-1, 1, 46),
    ),
)
assert height == expected_height and height.det() == 316

# Replace the primitive root basis by a deterministic E6+E8 simple system.
simple_in_root_block, cartan = deterministic_simple_roots(root_block)
assert abs(simple_in_root_block.det()) == 1
components = connected_components(cartan)
assert tuple(map(len, components)) == (6, 8)
assert tuple(
    abs(cartan.matrix_from_rows_and_columns(component, component).det())
    for component in components
) == (3, 1)
simple_change = block_diagonal_matrix(
    simple_in_root_block, identity_matrix(ZZ, 3)
)
simple_frame = simple_change * root_adapted * simple_change.transpose()
simple_root_block = simple_frame[:14, :14]
simple_coupling = simple_frame[:14, 14:]
simple_mw_block = simple_frame[14:, 14:]
assert simple_root_block == cartan
assert (
    simple_mw_block
    - simple_coupling.transpose()
    * simple_root_block.inverse()
    * simple_coupling
    == height
)

# Count the raw shell without storing its hundreds of millions of vectors.
raw_norm16 = exact_norm_count(simple_frame, 16)
root_norm16 = exact_norm_count(simple_root_block, 16)
horizontal_norm16 = raw_norm16 - root_norm16
assert (raw_norm16, root_norm16, horizontal_norm16) == (
    219758670,
    80751870,
    139006800,
)

# Enumerate the nonzero MW projections up to the global sign v ~ -v.
height_scale = lcm(entry.denominator() for entry in height.list())
scaled_height = (height_scale * height).change_ring(ZZ)
mw_result = pari(scaled_height).qfminim(height_scale * 16)
mw_vectors = {}
for column in matrix(ZZ, mw_result[2]).columns():
    for sign in (1, -1):
        value = sign * vector(ZZ, column)
        if value * height * value > 16:
            continue
        canonical = min(tuple(value), tuple(-value))
        mw_vectors[canonical] = vector(ZZ, canonical)
mw_vectors = tuple(
    sorted(mw_vectors.values(), key=lambda value: (value * height * value, tuple(value)))
)
assert len(mw_vectors) == 10

component_weights = tuple(
    dominant_weights(cartan, component, QQ(16))
    for component in components
)
assert tuple(len(weights) for weights in component_weights) == (45, 12)

# The dominant labels p satisfy
#   16 = p^t*A^-1*p + z^t*H*z,
# and x=A^-1*(p-B*z) must be integral.  Nonnegativity of A^-1 makes the
# recursive coordinate bounds exhaustive.
cartan_inverse = cartan.inverse()
dominant_orbits = {}
for mw in mw_vectors:
    mw_norm = mw * height * mw
    for first_labels, first_norm in component_weights[0]:
        for second_labels, second_norm in component_weights[1]:
            if mw_norm + first_norm + second_norm != 16:
                continue
            labels = vector(ZZ, [0] * 14)
            for index, value in zip(components[0], first_labels):
                labels[index] = value
            for index, value in zip(components[1], second_labels):
                labels[index] = value
            root_coordinates = cartan_inverse * (labels - simple_coupling * mw)
            if not all(value in ZZ for value in root_coordinates):
                continue
            witness = vector(ZZ, list(root_coordinates) + list(mw))
            assert witness * simple_frame * witness == 16
            assert (
                simple_root_block * root_coordinates + simple_coupling * mw
                == labels
            )
            dominant_orbits[tuple(witness)] = (tuple(mw), tuple(labels))
assert len(dominant_orbits) == 63

q8_ns = block_diagonal_matrix(U, -simple_frame)
histogram = Counter()
mw4_hits = []
nonprimitive = 0
for witness_tuple, (mw_tuple, labels_tuple) in sorted(dominant_orbits.items()):
    witness = vector(ZZ, witness_tuple)
    # This order has intersection two with the old fiber.
    q8_fiber = vector(ZZ, [4, 2] + list(witness))
    result = child_frame(q8_ns, q8_fiber, determinant=948)
    if result is None:
        nonprimitive += 1
        continue
    q8_child, q8_neighbor_basis = result
    _, _, invariants = roots_and_data(q8_child)
    histogram[invariants] += 1
    if invariants != (13, 312, 4):
        continue

    # Connected simply-laced rank 13 and determinant four identifies D13.
    (
        d13_simple,
        d13_root_adapted_basis,
        d13_root_adapted,
        d13_height,
    ) = d13_root_adaptation(q8_child)
    witness_root_adapted = witness * simple_change
    witness_q6_child = witness_root_adapted * root_mw_basis
    assert witness_q6_child * q6_child * witness_q6_child == 16
    fiber_q6_ns = vector(ZZ, [4, 2] + list(witness_q6_child))
    fiber_source_ns = fiber_q6_ns * q6_neighbor_basis
    assert fiber_source_ns * source_ns * fiber_source_ns == 0

    mw4_hits.append(
        {
            "mw_projection": list(map(int, mw_tuple)),
            "dominant_labels": list(map(int, labels_tuple)),
            "witness_simple_frame": list(map(int, witness)),
            "witness_root_adapted": list(map(int, witness_root_adapted)),
            "witness_q6_child": list(map(int, witness_q6_child)),
            "fiber_q6_ns": list(map(int, fiber_q6_ns)),
            "fiber_source_h3_ns": list(map(int, fiber_source_ns)),
            "child_root_data": [13, 312, 4],
            "child_type": "D13/MW4",
            "child_frame": matrix_rows(q8_child),
            "neighbor_basis_in_q6_ns": matrix_rows(q8_neighbor_basis),
            "d13_simple_roots": matrix_rows(d13_simple),
            "d13_root_adapted_basis_in_child": matrix_rows(
                d13_root_adapted_basis
            ),
            "d13_root_adapted_gram": matrix_rows(d13_root_adapted),
            "d13_mw_height": rational_rows(d13_height),
        }
    )

expected_histogram = {
    (15, 366, 2): 4,
    (14, 184, 12): 4,
    (15, 184, 32): 10,
    (14, 184, 16): 8,
    (14, 312, 3): 14,
    (14, 172, 24): 7,
    (16, 288, 21): 5,
    (14, 144, 48): 4,
    (15, 204, 18): 3,
    (13, 312, 4): 2,
}
assert nonprimitive == 2
assert dict(histogram) == expected_histogram
assert len(mw4_hits) == 2
assert {tuple(hit["mw_projection"]) for hit in mw4_hits} == {
    (-2, 0, 0),
    (0, -2, 0),
}

# Link the exhaustive orbit list to the independently chamber-certified nef
# representative.  Four simple-root reflections move it to the first
# dominant orbit above; its (4,2) presentation has old-fiber degree two.
nef_q8_witness_q6_child = vector(
    ZZ,
    (-5, -4, -3, 4, 5, 7, 10, 8, 6, 4, 2, -4, 2, -2, -4, 0, -2),
)
assert nef_q8_witness_q6_child * q6_child * nef_q8_witness_q6_child == 16
nef_root_adapted = nef_q8_witness_q6_child * root_mw_basis.inverse()
nef_simple = vector(ZZ, nef_root_adapted * simple_change.inverse())
assert tuple(nef_simple[-3:]) == (0, -2, 0)
nef_to_dominant = []
dominant = vector(ZZ, nef_simple)
while True:
    labels = vector(ZZ, dominant * simple_frame)[:14]
    negative = [index for index, value in enumerate(labels) if value < 0]
    if not negative:
        break
    index = negative[0]
    pairing = labels[index]
    dominant[index] -= pairing
    nef_to_dominant.append((index + 1, int(pairing)))
assert tuple(nef_to_dominant) == ((14, -1), (13, -1), (12, -1), (11, -1))
assert tuple(dominant) == tuple(mw4_hits[0]["witness_simple_frame"])

q6_child_ns = block_diagonal_matrix(U, -q6_child)
nef_q8_fiber = vector(ZZ, [4, 2] + list(nef_q8_witness_q6_child))
nef_q8_data = child_frame(q6_child_ns, nef_q8_fiber, determinant=948)
assert nef_q8_data is not None
nef_q8_child, nef_q8_neighbor_basis = nef_q8_data
assert roots_and_data(nef_q8_child)[2] == (13, 312, 4)
(
    nef_d13_simple,
    nef_d13_root_adapted_basis,
    nef_d13_root_adapted,
    nef_d13_height,
) = d13_root_adaptation(nef_q8_child)
nef_q8_fiber_source = nef_q8_fiber * q6_neighbor_basis
assert nef_q8_fiber_source * source_ns * nef_q8_fiber_source == 0
pinned_d13 = load_gram(PINNED_D13_FRAME)
assert pinned_d13 == matrix(ZZ, mw4_hits[0]["d13_root_adapted_gram"])

payload = {
    "status": "PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION",
    "source_frame": str(SOURCE_FRAME.relative_to(ROOT)),
    "q6": {
        "factor_order": [int(3), int(2)],
        "old_fiber_degree": 2,
        "witness": list(map(int, q6_witness)),
        "child_root_data": list(map(int, q6_root_data)),
        "child_type": "E8+E6/MW3",
        "neighbor_basis_in_source_ns": matrix_rows(q6_neighbor_basis),
        "root_mw_basis_in_child": matrix_rows(root_mw_basis),
        "root_adapted_gram": matrix_rows(root_adapted),
        "mw_height": rational_rows(height),
    },
    "q8": {
        "factor_order": [int(4), int(2)],
        "old_fiber_degree": 2,
        "raw_norm16_signed_vectors": int(raw_norm16),
        "root_span_norm16_signed_vectors": int(root_norm16),
        "horizontal_norm16_signed_vectors": int(horizontal_norm16),
        "simple_root_change_in_root_block": matrix_rows(simple_in_root_block),
        "simple_frame_gram": matrix_rows(simple_frame),
        "root_components": [list(map(int, component)) for component in components],
        "mw_projection_representatives": [
            {
                "vector": list(map(int, value)),
                "norm": str(value * height * value),
            }
            for value in mw_vectors
        ],
        "component_dominant_weight_counts": [45, 12],
        "dominant_orbits": len(dominant_orbits),
        "primitive_neighbors": sum(histogram.values()),
        "nonprimitive_orbits": nonprimitive,
        "root_histogram": [
            {
                "root_rank": int(invariants[0]),
                "root_count": int(invariants[1]),
                "root_determinant": int(invariants[2]),
                "mw_rank": 17 - int(invariants[0]),
                "orbit_count": int(count),
            }
            for invariants, count in sorted(histogram.items())
        ],
        "d13_mw4_hits": mw4_hits,
        "nef_representative": {
            "witness_q6_child": list(map(int, nef_q8_witness_q6_child)),
            "mw_projection": list(map(int, nef_simple[-3:])),
            "to_dominant_reflections": [
                [int(index), int(pairing)]
                for index, pairing in nef_to_dominant
            ],
            "fiber_q6_ns": list(map(int, nef_q8_fiber)),
            "fiber_source_h3_ns": list(map(int, nef_q8_fiber_source)),
            "child_frame": matrix_rows(nef_q8_child),
            "neighbor_basis_in_q6_ns": matrix_rows(nef_q8_neighbor_basis),
            "d13_simple_roots": matrix_rows(nef_d13_simple),
            "d13_root_adapted_basis_in_child": matrix_rows(
                nef_d13_root_adapted_basis
            ),
            "d13_root_adapted_gram": matrix_rows(nef_d13_root_adapted),
            "d13_mw_height": rational_rows(nef_d13_height),
        },
        "pinned_d13_frame": str(PINNED_D13_FRAME.relative_to(ROOT)),
    },
}

args.output.parent.mkdir(parents=True, exist_ok=True)
def json_default(value):
    if isinstance(value, Integer):
        return int(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


args.output.write_text(
    json.dumps(payload, indent=2, sort_keys=True, default=json_default) + "\n"
)

for index, hit in enumerate(mw4_hits, 1):
    print(
        "H3Q6Q8|hit={}|mw_projection={}|labels={}|witness={}|"
        "child=D13/MW4|root_data=13,312,4".format(
            index,
            tuple(hit["mw_projection"]),
            tuple(hit["dominant_labels"]),
            tuple(hit["witness_simple_frame"]),
        ),
        flush=True,
    )
print(
    "H3Q6Q8|raw_norm16={}|horizontal_norm16={}|mw_projections={}|"
    "dominant_orbits={}|primitive_neighbors={}|mw4_hits={}|"
    "status=PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION".format(
        raw_norm16,
        horizontal_norm16,
        len(mw_vectors),
        len(dominant_orbits),
        sum(histogram.values()),
        len(mw4_hits),
    ),
    flush=True,
)
