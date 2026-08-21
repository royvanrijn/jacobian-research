from sage.all import *
from pathlib import Path
from itertools import product
import hashlib


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "elkies-k3" / "data"
U = matrix(ZZ, [[0, 1], [1, 0]])


def load_matrix(path):
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            rows.append([ZZ(value) for value in line.split()])
    return matrix(ZZ, rows)


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
    if current == -1:
        coefficients = [-value for value in coefficients]
    return vector(ZZ, coefficients)


def reconstruct_step(parent, a, b, v):
    ns_parent = block_diagonal_matrix(U, -parent)
    fiber = vector(ZZ, [ZZ(a), ZZ(b)] + list(v))
    assert fiber * ns_parent * fiber == 0
    assert gcd([abs(ZZ(value)) for value in ns_parent * fiber]) == 1

    mate = bezout_vector(list(ns_parent * fiber))
    assert fiber * ns_parent * mate == 1
    mate_square = ZZ(mate * ns_parent * mate)
    assert mate_square % 2 == 0
    mate -= (mate_square // 2) * fiber
    assert mate * ns_parent * mate == 0 and fiber * ns_parent * mate == 1

    orthogonal = matrix(
        ZZ, [list(fiber * ns_parent), list(mate * ns_parent)]
    ).right_kernel_matrix()
    child = -(orthogonal * ns_parent * orthogonal.transpose())
    transport = matrix(ZZ, [list(fiber), list(mate)] + orthogonal.rows())
    assert transport.det() in (-1, 1)
    assert transport * ns_parent * transport.transpose() == block_diagonal_matrix(
        U, -child
    )
    return child, transport


STEPS = (
    {
        "name": "rank17_to_q25_mw7",
        "parent": DATA / "lattice" / "rank17_gram.txt",
        "child": DATA / "fibrations" / "q25_mw7_frame.txt",
        "a": 5,
        "b": 5,
        "v": (-1, 0, -4, 3, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, -3, 0, 0),
    },
    {
        "name": "q25_mw7_to_q25_mw4",
        "parent": DATA / "fibrations" / "q25_mw7_frame.txt",
        "child": DATA / "fibrations" / "q25_mw4_frame.txt",
        "a": 2,
        "b": 2,
        "v": (-1, -2, 1, 0, 1, 1, 2, -3, 0, -2, 0, 1, 0, 0, -1, 0, 0),
    },
    {
        "name": "q25_mw4_to_a5_d4_2a2_a1_mw3",
        "parent": DATA / "fibrations" / "q25_mw4_frame.txt",
        "child": DATA / "fibrations" / "mw3_a5_d4_a2a2_a1_frame.txt",
        "a": 2,
        "b": 2,
        "v": (-1, 0, 0, 2, 0, 2, -1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0),
    },
    {
        "name": "a5_d4_2a2_a1_mw3_to_e6_d4_2a2_a1_mw2",
        "parent": DATA / "fibrations" / "mw3_a5_d4_a2a2_a1_frame.txt",
        "child": DATA / "fibrations" / "mw2_e6_d4_a2a2_a1_frame.txt",
        "a": 2,
        "b": 2,
        "v": (0, -2, -2, 0, 1, 0, 3, 2, 4, 0, 3, 0, -1, -3, -1, -4, 0),
    },
)


transports = []
for index, step in enumerate(STEPS, 1):
    parent = load_matrix(step["parent"])
    expected_child = load_matrix(step["child"])
    child, transport = reconstruct_step(
        parent, step["a"], step["b"], vector(ZZ, step["v"])
    )
    assert child == expected_child
    transports.append(transport)
    print(
        f"MW2TRANSPORT|step={index}|name={step['name']}"
        f"|det={transport.det()}|child_exact=1",
        flush=True,
    )

standard_transport = transports[-1]
for transport in reversed(transports[:-1]):
    standard_transport *= transport

rank17_frame = load_matrix(STEPS[0]["parent"])
terminal_frame = load_matrix(STEPS[-1]["child"])
rank17_ns = block_diagonal_matrix(U, -rank17_frame)
terminal_ns = block_diagonal_matrix(U, -terminal_frame)
assert standard_transport.det() in (-1, 1)
assert standard_transport * rank17_ns * standard_transport.transpose() == terminal_ns


def qform_from_gram(gram):
    coefficients = []
    for row in range(gram.nrows()):
        for column in range(row, gram.ncols()):
            coefficients.append(
                gram[row, row] // 2 if row == column else gram[row, column]
            )
    return QuadraticForm(ZZ, gram.nrows(), coefficients)


def positive_root(root):
    return next(value for value in root if value != 0) > 0


half_roots = [
    vector(ZZ, root)
    for root in qform_from_gram(terminal_frame).short_vector_list_up_to_length(2, True)[1]
]
roots = half_roots + [-root for root in half_roots]
root_graph = Graph()
root_graph.add_vertices(range(len(roots)))
for left in range(len(roots)):
    for right in range(left):
        if roots[left] * terminal_frame * roots[right] != 0:
            root_graph.add_edge(left, right)
components = sorted(root_graph.connected_components(sort=False), key=len, reverse=True)
assert [len(component) for component in components] == [72, 24, 6, 6, 2]

# Lexicographic positivity is additive and selects one root from each +/- pair.
# Indecomposable positive roots give a deterministic simple-root system.
component_simple_bases = []
for component in components:
    component_roots = [roots[index] for index in component]
    positive = [root for root in component_roots if positive_root(root)]
    positive_set = {tuple(root) for root in positive}
    simple = []
    for root in positive:
        if not any(tuple(root - summand) in positive_set for summand in positive):
            simple.append(root)
    simple_basis = matrix(ZZ, simple)
    full_basis = matrix(ZZ, component_roots).row_module().basis_matrix()
    assert simple_basis.rank() == full_basis.rank()
    assert simple_basis.row_module() == full_basis.row_module()
    simple_gram = simple_basis * terminal_frame * simple_basis.transpose()
    assert all(simple_gram[i, i] == 2 for i in range(simple_gram.nrows()))
    assert all(
        simple_gram[i, j] in (0, -1)
        for i in range(simple_gram.nrows()) for j in range(i)
    )
    component_simple_bases.append(simple_basis)

component_ranks = tuple(basis.nrows() for basis in component_simple_bases)
assert component_ranks == (6, 4, 2, 2, 1)
R = block_matrix([[basis] for basis in component_simple_bases], subdivide=False)
GR = R * terminal_frame * R.transpose()
assert abs(GR.det()) == 216
C = (R * terminal_frame).right_kernel_matrix()
GC = C * terminal_frame * C.transpose()
assert C.rank() == 2


def fractional_class(coordinates):
    return tuple(QQ(value) - floor(QQ(value)) for value in coordinates)


def class_multiple(multiplier, point):
    return fractional_class(tuple(multiplier * value for value in point))


def class_order(point, exponent):
    zero = tuple(QQ(0) for _ in point)
    return next(
        order for order in range(1, exponent + 1)
        if class_multiple(order, point) == zero
    )


# Enumerate terminal_frame/(R + C), then saturate its MW projection.
A = block_matrix([[R], [C]], subdivide=False)
glue_index = abs(A.det())
assert glue_index == 36
A_inverse = A.inverse()


def coset_key(point):
    return fractional_class(vector(QQ, point) * A_inverse)


zero = vector(ZZ, [0] * 17)
cosets = {coset_key(zero): zero}
queue = [zero]
head = 0
while head < len(queue) and len(cosets) < glue_index:
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
                if len(cosets) == glue_index:
                    break
        if len(cosets) == glue_index:
            break
assert len(cosets) == glue_index

GR_inverse = GR.inverse()
GC_inverse = GC.inverse()


def project_mw(point):
    point = vector(QQ, point)
    return point - (point * terminal_frame * R.transpose()) * GR_inverse * R


def coordinates_in_C(point):
    return (vector(QQ, point) * terminal_frame * C.transpose()) * GC_inverse


projected_cosets = [
    (coordinates_in_C(project_mw(representative)), representative)
    for representative in cosets.values()
]
generators = [vector(QQ, (1, 0)), vector(QQ, (0, 1))]
generators += [coordinates for coordinates, _ in projected_cosets]
denominator = lcm(QQ(value).denominator() for row in generators for value in row)
MW_integer = matrix(ZZ, [
    [ZZ(denominator * value) for value in row] for row in generators
]).row_module().basis_matrix()
MW_basis = MW_integer.change_ring(QQ) / denominator
H = MW_basis * GC * MW_basis.transpose()
assert H.det() == QQ(79) / 18

target_scaled = matrix(ZZ, [[9, 2], [2, 18]])
scaled = (6 * H).change_ring(ZZ)
minimum = pari(scaled).qfminim(18)
representatives = list(matrix(ZZ, minimum[2]).columns())
vectors = representatives + [-point for point in representatives]
transforms = []
for first in vectors:
    if first * scaled * first != 9:
        continue
    for second in vectors:
        if second * scaled * second != 18:
            continue
        transform = matrix(ZZ, [first, second])
        if abs(transform.det()) == 1 and transform * scaled * transform.transpose() == target_scaled:
            transforms.append(transform)
assert len(transforms) == 2 and transforms[0] == -transforms[1]

component_bounds = []
start = 0
for basis in component_simple_bases:
    component_bounds.append((start, start + basis.nrows()))
    start += basis.nrows()


def lift_target_basis(transform):
    target_basis = transform * MW_basis
    lifts = []
    classes = []
    root_coordinates_all = []
    for target_vector in target_basis.rows():
        lift = None
        for projected, representative in projected_cosets:
            difference = target_vector - projected
            if all(QQ(value).denominator() == 1 for value in difference):
                lift = vector(QQ, representative) + difference * C
                break
        assert lift is not None and all(QQ(value).denominator() == 1 for value in lift)
        lift = vector(ZZ, lift)
        root_part = vector(QQ, lift) - target_vector * C
        root_coordinates = (root_part * terminal_frame * R.transpose()) * GR_inverse
        assert root_part == root_coordinates * R
        lifts.append(lift)
        classes.append([
            fractional_class(root_coordinates[left:right])
            for left, right in component_bounds
        ])
        root_coordinates_all.append(root_coordinates)
    return target_basis, lifts, classes, root_coordinates_all


def support_profile(classes):
    return tuple(
        tuple(class_order(classes[row][component], 6) for component in range(5))
        for row in range(2)
    )


candidates = []
for transform in transforms:
    target_basis, lifts, classes, root_coordinates_all = lift_target_basis(transform)
    orders = support_profile(classes)
    # Orders encode P1=(nonzero E6,0,one A2,A1) and
    # P2=(inverse E6,nonzero D4,the other A2,0), before ordering the A2s.
    assert orders[0][0:2] == (3, 1) and orders[0][4] == 2
    assert orders[1][0:2] == (3, 2) and orders[1][4] == 1
    assert sorted((orders[0][2], orders[0][3])) == [1, 3]
    assert sorted((orders[1][2], orders[1][3])) == [1, 3]
    assert all(orders[0][component] * orders[1][component] == 3 for component in (2, 3))
    candidates.append((target_basis, lifts, classes, root_coordinates_all))

# The two choices differ by simultaneous section negation.  Pin the one whose
# pair of raw integral lifts is lexicographically smaller; the opposite choice
# is obtained geometrically by (x,y) -> (x,-y).
target_basis, raw_lifts, root_classes, root_coordinates_all = min(
    candidates, key=lambda item: tuple(item[1][0]) + tuple(item[1][1])
)


def minimized_component_coordinates(point, gram, expected_norm):
    rank = gram.nrows()
    candidates = []
    for shift in product(range(-3, 4), repeat=rank):
        coordinates = vector(QQ, point) + vector(ZZ, shift)
        if coordinates * gram * coordinates == expected_norm:
            candidates.append(coordinates)
    assert candidates
    return min(candidates, key=lambda row: tuple(row))


expected_nonzero_norm = (QQ(4) / 3, QQ(1), QQ(2) / 3, QQ(2) / 3, QQ(1) / 2)
section_lifts = []
for row in range(2):
    minimized = []
    for component, gram in enumerate(
        basis * terminal_frame * basis.transpose() for basis in component_simple_bases
    ):
        point = root_classes[row][component]
        order = class_order(point, 6)
        expected = QQ(0) if order == 1 else expected_nonzero_norm[component]
        minimized.extend(minimized_component_coordinates(point, gram, expected))
    root_part = vector(QQ, minimized) * R
    lift = target_basis.row(row) * C + root_part
    assert all(QQ(value).denominator() == 1 for value in lift)
    section_lifts.append(vector(ZZ, lift))

section_zero_intersections = tuple(
    ZZ((lift * terminal_frame * lift - 4) / 2) for lift in section_lifts
)
assert section_zero_intersections == (0, 1)
section_pair_intersection = ZZ(
    2 + sum(section_zero_intersections)
    - section_lifts[0] * terminal_frame * section_lifts[1]
)
assert section_pair_intersection == 2

# Put the A2 met by P2 first, matching profiles
# P1=(1,0,0,1,1), P2=(2,d1,1,0,0).
a2_indices = sorted(
    (2, 3), key=lambda component: class_order(root_classes[0][component], 3)
)
component_order = (0, 1, a2_indices[0], a2_indices[1], 4)
ordered_root_bases = [component_simple_bases[index] for index in component_order]
ordered_R = block_matrix([[basis] for basis in ordered_root_bases], subdivide=False)

terminal_fiber = vector(ZZ, [1, 0] + [0] * 17)
terminal_zero = vector(ZZ, [-1, 1] + [0] * 17)
terminal_p1 = vector(ZZ, [section_zero_intersections[0] + 1, 1] + list(section_lifts[0]))
terminal_p2 = vector(ZZ, [section_zero_intersections[1] + 1, 1] + list(section_lifts[1]))
explicit_terminal_basis = matrix(
    ZZ,
    [terminal_fiber, terminal_zero]
    + [vector(ZZ, [0, 0] + list(root)) for root in ordered_R.rows()]
    + [terminal_p1, terminal_p2],
)
assert explicit_terminal_basis.det() in (-1, 1)
section_root_intersections = -matrix(ZZ, section_lifts) * terminal_frame * ordered_R.transpose()
assert section_root_intersections == matrix(ZZ, [
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
    [0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0],
])

explicit_transport = explicit_terminal_basis * standard_transport
assert explicit_transport.det() in (-1, 1)
explicit_gram = explicit_terminal_basis * terminal_ns * explicit_terminal_basis.transpose()
assert explicit_transport * rank17_ns * explicit_transport.transpose() == explicit_gram

explicit_basis_path = DATA / "fibrations" / "mw2_e6_d4_a2a2_a1_explicit_basis.txt"
transport_path = DATA / "fibrations" / "mw2_e6_d4_a2a2_a1_ns_transport_to_rank17.txt"
assert load_matrix(explicit_basis_path) == explicit_terminal_basis
assert load_matrix(transport_path) == explicit_transport
explicit_basis_hash = hashlib.sha256(explicit_basis_path.read_bytes()).hexdigest()
transport_hash = hashlib.sha256(transport_path.read_bytes()).hexdigest()

rank17_in_explicit = explicit_transport.inverse()
assert rank17_in_explicit in MatrixSpace(ZZ, 19)
rank17_fiber = vector(ZZ, rank17_in_explicit.row(0))
rank17_mate = vector(ZZ, rank17_in_explicit.row(1))
rank17_zero = rank17_mate - rank17_fiber
assert rank17_fiber * explicit_gram * rank17_fiber == 0
assert rank17_mate * explicit_gram * rank17_mate == 0
assert rank17_fiber * explicit_gram * rank17_mate == 1
assert rank17_zero * explicit_gram * rank17_zero == -2
assert rank17_fiber * explicit_gram * rank17_zero == 1

backward_stage_names = (
    "a5_d4_2a2_a1_mw3",
    "q25_mw4",
    "q25_mw7",
    "rank17",
)
backward_classes = []
partial_transport = explicit_terminal_basis
for stage_name, transition in zip(backward_stage_names, reversed(transports)):
    partial_transport *= transition
    inverse = partial_transport.inverse()
    assert inverse in MatrixSpace(ZZ, 19)
    old_fiber = vector(ZZ, inverse.row(0))
    old_mate = vector(ZZ, inverse.row(1))
    old_zero = old_mate - old_fiber
    assert old_fiber * explicit_gram * old_fiber == 0
    assert old_mate * explicit_gram * old_mate == 0
    assert old_fiber * explicit_gram * old_mate == 1
    assert old_zero * explicit_gram * old_zero == -2
    assert old_fiber * explicit_gram * old_zero == 1
    backward_classes.append((stage_name, old_fiber, old_mate, old_zero))
assert backward_classes[-1][1] == rank17_fiber
assert backward_classes[-1][2] == rank17_mate
assert backward_classes[-1][3] == rank17_zero


def reduce_against_explicit_curves(divisor):
    """Remove fixed (-2)-curve components visible in the explicit basis."""
    # Sage may return the same mutable vector when the ring is unchanged.
    reduced = vector(ZZ, list(divisor))
    reflections = []
    for _ in range(10000):
        pairings = reduced * explicit_gram
        negative = [index for index in range(1, 19) if pairings[index] < 0]
        if not negative:
            return reduced, reflections
        index = min(negative, key=lambda candidate: pairings[candidate])
        multiplicity = -ZZ(pairings[index])
        reduced[index] -= multiplicity
        reflections.append((index, multiplicity))
        assert reduced * explicit_gram * reduced == divisor * explicit_gram * divisor
    raise RuntimeError("explicit-curve fixed-component reduction did not terminate")

backward_path = DATA / "fibrations" / "mw2_e6_d4_a2a2_a1_inverse_neighbor_classes.tsv"
backward_rows = [
    line.split("\t")
    for line in backward_path.read_text().splitlines()
    if line.strip() and not line.startswith("#")
]
assert backward_rows[0] == ["stage", "class", "coordinates"]
recorded_backward = {
    (stage, class_name): vector(ZZ, coordinates.split(","))
    for stage, class_name, coordinates in backward_rows[1:]
}
for stage_name, old_fiber, old_mate, old_zero in backward_classes:
    assert recorded_backward[(stage_name, "fiber")] == old_fiber
    assert recorded_backward[(stage_name, "mate")] == old_mate
    assert recorded_backward[(stage_name, "zero")] == old_zero
backward_hash = hashlib.sha256(backward_path.read_bytes()).hexdigest()

print(
    f"MW2TRANSPORT|standard_composite=PASS|det={standard_transport.det()}"
    f"|max_entry={max(abs(value) for value in standard_transport.list())}",
    flush=True,
)
print(
    f"MW2TRANSPORT|explicit_basis=PASS|det={explicit_terminal_basis.det()}"
    f"|P1.O={section_zero_intersections[0]}|P2.O={section_zero_intersections[1]}"
    f"|P1.P2={section_pair_intersection}",
    flush=True,
)
print("MW2TRANSPORT|row_order=F,O,E6x6,D4x4,A2_P2x2,A2_P1x2,A1,P1,P2", flush=True)
print(f"MW2TRANSPORT|rank17_fiber_in_explicit={tuple(rank17_fiber)}", flush=True)
print(f"MW2TRANSPORT|rank17_mate_in_explicit={tuple(rank17_mate)}", flush=True)
print(f"MW2TRANSPORT|rank17_zero_in_explicit={tuple(rank17_zero)}", flush=True)
for stage_name, old_fiber, old_mate, old_zero in backward_classes:
    print(f"MW2TRANSPORT|backward_stage={stage_name}|fiber={tuple(old_fiber)}", flush=True)
    print(
        f"MW2TRANSPORT|backward_stage={stage_name}"
        f"|fiber_intersections={tuple(old_fiber * explicit_gram)}",
        flush=True,
    )
    print(f"MW2TRANSPORT|backward_stage={stage_name}|mate={tuple(old_mate)}", flush=True)
    print(f"MW2TRANSPORT|backward_stage={stage_name}|zero={tuple(old_zero)}", flush=True)
    reduced_fiber, fixed_reflections = reduce_against_explicit_curves(old_fiber)
    print(
        f"MW2TRANSPORT|backward_stage={stage_name}"
        f"|explicit_curve_reduced_fiber={tuple(reduced_fiber)}"
        f"|reflections={len(fixed_reflections)}"
        f"|intersections={tuple(reduced_fiber * explicit_gram)}",
        flush=True,
    )
print(f"MW2TRANSPORT|explicit_basis_sha256={explicit_basis_hash}", flush=True)
print(f"MW2TRANSPORT|transport_sha256={transport_hash}", flush=True)
print(f"MW2TRANSPORT|backward_classes_sha256={backward_hash}", flush=True)
print("MW2TRANSPORT|status=PASS", flush=True)
