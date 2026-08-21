#!/usr/bin/env sage
"""Certify the first low-degree neighbor from the labeled H3 Kumar frame.

The source is the E7+E8/MW2 frame with height Gram
[[21/2,3],[3,46]].  The complete constrained q=6 shell is enumerated by
``search_alternate_fibrations.sage``.  This checker takes its shortest
representative, chooses the factor order giving old-fiber degree two, reduces
it against the displayed effective E7 components, and identifies the movable
class as ``O+(-P1)-F``.  It also reconstructs the child frame and verifies its
E8+E6 root data.  It then recovers the saturated rank-three Mordell--Weil
height lattice, proves its optimal pole profile, and distinguishes it from the
older H2 q=60 E8+E6/MW3 frame at the level of the height lattice.  The final
stage certifies the exact q=8 continuation to D13/MW4, including its
fixed-zero-section reduction and a lattice-theoretic nefness proof.
"""

from itertools import combinations
from pathlib import Path

from sage.all import *


BASE = Path(__file__).resolve().parents[1]
FRAME = BASE / "data/fibrations/kumar_e7e8_mw2_frame_3.txt"
Q60_FRAME = BASE / "data/fibrations/kumar_q60_e8_e6_mw3_frame.txt"
Q8_D13_FRAME = (
    BASE / "data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt"
)


def load_gram(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ],
    )


frame = load_gram(FRAME)
assert frame.nrows() == 17 and frame.det() == 948
NS = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -frame)

F = vector(ZZ, [1, 0] + [0] * 17)
O = vector(ZZ, [-1, 1] + [0] * 17)
simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(index == node) for index in range(17)])
    for node in range(15)
)

# The search records (a,b)=(2,3).  Swapping the two U generators gives the
# same neighbor presentation with D.F=2, which is the useful geometric order.
witness = vector(
    ZZ,
    [0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
)
D_raw = vector(ZZ, [3, 2] + list(witness))
assert D_raw * NS * D_raw == 0
assert D_raw * NS * F == 2
assert gcd(tuple(NS * D_raw)) == 1

# Deterministic fixed-component reflections in the displayed E7 simple roots.
reflection_nodes = []
D = vector(ZZ, list(D_raw))
while True:
    negative = tuple(
        (node, ZZ(D * NS * curve))
        for node, curve in enumerate(simple)
        if D * NS * curve < 0
    )
    if not negative:
        break
    node, pairing = negative[0]
    reflection_nodes.append((node + 1, pairing))
    D = vector(ZZ, list(D + pairing * simple[node]))

expected_reflections = [
    1, 2, 4, 3, 5, 4, 2, 6, 5, 4, 3,
    1, 7, 6, 5, 4, 2, 3, 4, 5, 6, 7,
]
assert reflection_nodes == [(node, -1) for node in expected_reflections]

# The inverse section has the twice-minuscule E7 correction.  Its frame sign
# is the sign appearing in the q=6 witness; this makes its intersections with
# the effective simple components nonnegative.
twice_minuscule = (2, 3, 4, 6, 5, 4, 3)
minus_P1 = vector(
    ZZ,
    [5, 1]
    + [-value for value in twice_minuscule]
    + [0] * 8
    + [1, 0],
)
assert minus_P1 * NS * minus_P1 == -2
assert minus_P1 * NS * F == 1
assert minus_P1 * NS * O == 4
assert tuple(minus_P1 * NS * curve for curve in simple) == (
    0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0,
)

assert D == O + minus_P1 - F
assert D * NS * D == 0
assert D * NS * F == 2
assert D * NS * O == 1
assert D * NS * minus_P1 == 1
assert all(D * NS * curve >= 0 for curve in simple)

# Include the two affine identity components, not only the finite Dynkin
# simples.  The coefficients use Sage's pinned Cartan numbering.
highest_E7 = (2, 2, 3, 4, 3, 2, 1)
highest_E8 = (2, 3, 4, 6, 5, 4, 3, 2)
affine_E7 = F - sum((coefficient * simple[index]
                     for index, coefficient in enumerate(highest_E7)),
                    vector(ZZ, [0] * 19))
affine_E8 = F - sum((coefficient * simple[7 + index]
                     for index, coefficient in enumerate(highest_E8)),
                    vector(ZZ, [0] * 19))
assert affine_E7 * NS * affine_E7 == -2
assert affine_E8 * NS * affine_E8 == -2
assert (D * NS * affine_E7, D * NS * affine_E8) == (1, 2)


def bezout_vector_for_pairing(ns, fiber):
    pairings = tuple(ns * fiber)
    current = ZZ(0)
    result = [ZZ(0)] * 19
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


mate = bezout_vector_for_pairing(NS, D_raw)
assert D_raw * NS * mate == 1
mate -= (mate * NS * mate // 2) * D_raw
assert mate * NS * mate == 0 and D_raw * NS * mate == 1
orthogonal = matrix(ZZ, [list(D_raw * NS), list(mate * NS)]).right_kernel_matrix()
child = -(orthogonal * NS * orthogonal.transpose())
assert child.is_positive_definite() and child.det() == 948
minimum = pari(child).qfminim(2)
root_count = ZZ(minimum[0])
roots = matrix(ZZ, minimum[2]).transpose()
root_basis = roots.row_module().basis_matrix()
root_rank = root_basis.rank()
root_gram = root_basis * child * root_basis.transpose()
root_determinant = abs(root_gram.det())
assert (root_rank, root_count, root_determinant) == (14, 312, 3)


def root_torsion_order(basis):
    smith = basis.smith_form()[0]
    diagonal = [
        abs(ZZ(smith[index, index])) for index in range(basis.nrows())
    ]
    return prod(value for value in diagonal if value)


def projected_height_gram(gram, basis, lifts):
    root_form = basis * gram * basis.transpose()
    correction = (
        lifts * gram * basis.transpose() * root_form.inverse() * basis
    )
    projections = lifts - correction
    return projections * gram * projections.transpose()


# These three integral lifts are in the child coordinates constructed above.
# Their projected determinant equals det(frame)/det(E8+E6), so together with
# primitive roots they certify a saturated Mordell--Weil basis.
h3_lifts = matrix(ZZ, [
    [-5, -4, -3, 0, 0, 0, 0, 0, 0, 0, 0, -4, 1, 0, -4, 2, -2],
    [-10, -8, -6, 0, 0, 0, 0, 0, 0, 0, 0, -8, 4, 1, -8, 5, -4],
    [-5, -4, -3, 0, 0, 0, 0, 0, 0, 0, 0, -3, 2, 0, -4, 2, -2],
])
h3_height = projected_height_gram(child, root_basis, h3_lifts)
expected_h3_height = matrix(QQ, [
    [QQ(8) / 3, QQ(1) / 3, -1],
    [QQ(1) / 3, QQ(8) / 3, 1],
    [-1, 1, 46],
])
assert root_torsion_order(root_basis) == 1
assert h3_height == expected_h3_height
assert h3_height.det() == QQ(child.det()) / root_determinant == 316


def pole_order(height, coordinates):
    coordinates = vector(ZZ, coordinates)
    value = coordinates * height * coordinates
    fractional = value - value.floor()
    if fractional == 0:
        local_correction = QQ(0)
    else:
        assert fractional == QQ(2) / 3
        local_correction = QQ(4) / 3
    pole = (value - 4 + local_correction) / 2
    assert pole in ZZ
    return ZZ(pole)


def optimal_pole_basis(height):
    """Exhaustively minimize the largest P.O over unimodular MW bases."""
    scaled = (3 * height).change_ring(ZZ)
    form = QuadraticForm(ZZ, scaled)
    coordinate_upper = max(
        pole_order(height, row) for row in identity_matrix(ZZ, 3).rows()
    )
    for bound in range(coordinate_upper + 1):
        # Sage uses q(x)=x^t*(3H)*x/2.  The endpoint is exclusive.
        q_bound = ZZ(3 * (4 + 2 * bound) / 2) + 1
        short = []
        for shell in form.short_vector_list_up_to_length(
            q_bound, up_to_sign_flag=True
        ):
            for row in shell:
                row = vector(ZZ, row)
                if row and 0 <= pole_order(height, row) <= bound:
                    short.append(row)
        for rows in combinations(short, 3):
            basis = matrix(ZZ, [list(row) for row in rows])
            if abs(basis.det()) == 1:
                return (
                    bound,
                    tuple(pole_order(height, row) for row in basis.rows()),
                    basis,
                )
    raise RuntimeError("coordinate basis did not furnish a pole bound")


pole_bound, pole_profile, pole_basis = optimal_pole_basis(h3_height)
assert pole_bound == 21 and pole_profile == (0, 0, 21)
assert abs(pole_basis.det()) == 1

# Compare with the old H2 q=60 child using the same exact projection method.
# The scaled ternary forms are not isometric, so identical E8+E6 root data do
# not identify the two marked frames.
q60_frame = load_gram(Q60_FRAME)
q60_minimum = pari(q60_frame).qfminim(2)
q60_roots = matrix(ZZ, q60_minimum[2]).transpose()
q60_root_basis = q60_roots.row_module().basis_matrix()
q60_root_gram = q60_root_basis * q60_frame * q60_root_basis.transpose()
assert (
    q60_root_basis.rank(),
    ZZ(q60_minimum[0]),
    abs(q60_root_gram.det()),
    root_torsion_order(q60_root_basis),
) == (14, 312, 3, 1)
q60_lifts = matrix(ZZ, [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    [20, 16, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 14, 16, 17, 16, 8],
    [-25, -20, -15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -17, -20, -20, -19, -10],
])
q60_height = projected_height_gram(q60_frame, q60_root_basis, q60_lifts)
expected_q60_height = matrix(QQ, [
    [4, 0, 0],
    [0, QQ(20) / 3, 1],
    [0, 1, 12],
])
assert q60_height == expected_q60_height and q60_height.det() == 316
assert pari(3 * h3_height).qfisom(pari(3 * q60_height)) == 0


def all_roots_of_norm_two(gram):
    minimum = pari(gram).qfminim(2)
    half = [vector(ZZ, row) for row in matrix(ZZ, minimum[2]).transpose().rows()]
    return half + [-row for row in half]


def lex_positive(row):
    return next(value > 0 for value in row if value)


def deterministic_simple_roots(gram):
    roots = all_roots_of_norm_two(gram)
    positive = [row for row in roots if lex_positive(row)]
    positive_set = {tuple(row) for row in positive}
    simple_roots = []
    for row in positive:
        if not any(
            tuple(row - left) in positive_set
            for left in positive
            if left != row
        ):
            simple_roots.append(row)
    simple_roots = matrix(ZZ, [list(row) for row in simple_roots])
    assert simple_roots.nrows() == simple_roots.rank()
    return simple_roots, positive


def connected_components(cartan):
    unseen = set(range(cartan.nrows()))
    result = []
    while unseen:
        todo = [min(unseen)]
        unseen.remove(todo[0])
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


def highest_roots(gram, simple_roots, positive):
    cartan = simple_roots * gram * simple_roots.transpose()
    inverse_simple = simple_roots.pseudoinverse()
    result = []
    for component in connected_components(cartan):
        candidates = []
        for root in positive:
            coordinates = vector(QQ, root) * inverse_simple
            if not all(value in ZZ and value >= 0 for value in coordinates):
                continue
            support = tuple(
                index for index, value in enumerate(coordinates) if value
            )
            if support and all(index in component for index in support):
                candidates.append((sum(coordinates), vector(ZZ, root)))
        result.append(max(candidates, key=lambda item: item[0])[1])
    return tuple(result)


def neighbor_frame(ns, raw_fiber):
    mate = bezout_vector_for_pairing(ns, raw_fiber)
    assert raw_fiber * ns * mate == 1
    mate -= (mate * ns * mate // 2) * raw_fiber
    assert mate * ns * mate == 0 and raw_fiber * ns * mate == 1
    complement = matrix(
        ZZ, [list(raw_fiber * ns), list(mate * ns)]
    ).right_kernel_matrix()
    full_basis = matrix(
        ZZ,
        [list(raw_fiber), list(mate)]
        + [list(row) for row in complement.rows()],
    )
    assert abs(full_basis.det()) == 1
    return -(complement * ns * complement.transpose())


# The exact dominant q=8 orbit can be represented directly in the child
# coordinates constructed above.  In the deterministic E6+E8 chamber it has
# fundamental-weight pairings at one node in each component.
q8_witness = vector(ZZ, (
    -5, -4, -3, 4, 5, 7, 10, 8, 6, 4, 2, -4, 2, -2, -4, 0, -2,
))
assert q8_witness * child * q8_witness == 16
child_ns = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -child)
child_F = vector(ZZ, [1, 0] + [0] * 17)
child_O = vector(ZZ, [-1, 1] + [0] * 17)
q8_raw = vector(ZZ, [2, 4] + list(q8_witness))
assert q8_raw * child_ns * q8_raw == 0
assert q8_raw * child_ns * child_F == 4
assert q8_raw * child_ns * child_O == -2

# Reflecting once in O swaps the U factor order and exposes the movable
# degree-two representative.
q8_nef = q8_raw - 2 * child_O
assert q8_nef == vector(ZZ, [4, 2] + list(q8_witness))
assert q8_nef * child_ns * q8_nef == 0
assert q8_nef * child_ns * child_F == 2
assert q8_nef * child_ns * child_O == 2

child_simple, child_positive = deterministic_simple_roots(child)
assert child_simple.nrows() == 14
child_cartan = child_simple * child * child_simple.transpose()
assert tuple(sorted(map(len, connected_components(child_cartan)))) == (6, 8)

# Negate the deterministic positive roots to get the effective component
# chamber selected by the dominant q=8 witness.
effective_simple = tuple(
    vector(ZZ, [0, 0] + list(-row)) for row in child_simple.rows()
)
simple_pairings = tuple(q8_nef * child_ns * row for row in effective_simple)
assert simple_pairings == (
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
)
child_affine = tuple(
    child_F + vector(ZZ, [0, 0] + list(root))
    for root in highest_roots(child, child_simple, child_positive)
)
assert tuple(q8_nef * child_ns * row for row in child_affine) == (1, 0)
assert all(q8_nef * child_ns * row >= 0 for row in effective_simple + child_affine)

# Full nefness, not merely a displayed-component check.  A negative curve C
# has degree one or two because q8_nef.F=2.  Degree one is a section.  In the
# saturated MW quotient, q8_witness projects to -2 times the second reduced
# basis vector, so the center q8_witness/2 is the integral MW point -e2.  The
# nonzero MW minimum is 8/3>2; equality in the zero-MW coset would require the
# nonintegral frame vector q8_witness/2.  Hence no section is negative.
root_projection_matrix = (
    identity_matrix(QQ, 17)
    - child * root_basis.transpose() * root_gram.inverse() * root_basis
)
q8_projection = q8_witness * root_projection_matrix
h3_lift_projections = h3_lifts * root_projection_matrix
assert q8_projection == -2 * h3_lift_projections[1]
assert ZZ(pari(3 * h3_height).qfminim()[1]) == 8
assert any(value % 2 for value in q8_witness)

# If a degree-two (-2)-curve C=[k,2,w] were negative, then
# ||w-v||^2=2(D.C+1), forcing w=v and D.C=-1.  But C^2=-2 would then require
# v^2=4k+2, impossible because v^2=16.  Thus no negative bisection exists.
assert q8_witness * child * q8_witness == 16
assert (q8_witness * child * q8_witness - 2) % 4 != 0

q8_child = neighbor_frame(child_ns, q8_raw)
assert q8_child.det() == 948 and q8_child.is_positive_definite()
q8_minimum = pari(q8_child).qfminim(2)
q8_roots = matrix(ZZ, q8_minimum[2]).transpose()
q8_root_basis = q8_roots.row_module().basis_matrix()
q8_root_gram = q8_root_basis * q8_child * q8_root_basis.transpose()
assert (
    q8_root_basis.rank(),
    ZZ(q8_minimum[0]),
    abs(q8_root_gram.det()),
) == (13, 312, 4)

# Pin a compact root-adapted presentation for the next neighbor search.  The
# first 13 rows are the deterministic simple roots; the last four are a
# saturated MW basis shifted to minimal D13 discriminant-coset weights.
q8_simple, _ = deterministic_simple_roots(q8_child)
assert q8_simple.nrows() == 13
q8_adapted_lifts = matrix(ZZ, [
    [5, 4, -2, -4, -5, -4, -3, -2, -1, 3, 3, -6, -2, 0, 2, 1, 0],
    [-5, -4, 2, 4, 5, 4, 3, 2, 1, -2, -3, 6, 2, 0, -2, -1, 0],
    [-10, -8, 2, 4, 5, 4, 3, 2, 1, -7, -6, 9, 2, -1, -3, -3, 0],
    [450, 360, -450, -630, -900, -720, -540, -360, -180,
     343, 165, -530, -173, -3, 155, 166, 1],
])
q8_adapted_basis = q8_simple.stack(q8_adapted_lifts)
assert abs(q8_adapted_basis.det()) == 1
q8_adapted = q8_adapted_basis * q8_child * q8_adapted_basis.transpose()
q8_root_cartan = q8_adapted[:13, :13]
q8_internal_height = (
    q8_adapted[13:, 13:]
    - q8_adapted[13:, :13]
    * q8_root_cartan.inverse()
    * q8_adapted[:13, 13:]
)
assert q8_internal_height == matrix(QQ, [
    [QQ(3) / 4, QQ(1) / 4, QQ(1) / 4, 0],
    [QQ(1) / 4, QQ(11) / 4, -QQ(1) / 4, -1],
    [QQ(1) / 4, -QQ(1) / 4, QQ(11) / 4, -1],
    [0, -1, -1, 46],
])

# The full orbit classifier pins a second, smaller root-adapted presentation
# together with its source transport.  Check its primitive D13 block and MW
# quotient here; the classifier verifies the complete integral change of
# basis and writes it to the generated JSON certificate.
q8_pinned = load_gram(Q8_D13_FRAME)
assert q8_pinned.det() == 948
q8_pinned_root = q8_pinned[:13, :13]
assert q8_pinned_root.det() == 4
q8_height = (
    q8_pinned[13:, 13:]
    - q8_pinned[13:, :13]
    * q8_pinned_root.inverse()
    * q8_pinned[:13, 13:]
)
assert q8_height == matrix(QQ, [
    [QQ(3) / 4, QQ(1) / 4, -QQ(1) / 4, 0],
    [QQ(1) / 4, QQ(11) / 4, QQ(1) / 4, 1],
    [-QQ(1) / 4, QQ(1) / 4, QQ(11) / 4, -1],
    [0, 1, -1, 46],
])
assert q8_height.det() == q8_internal_height.det() == QQ(237)
assert pari(4 * q8_height).qfisom(pari(4 * q8_internal_height)) != 0

print(
    "H3Q6|source=E7+E8/MW2|q=6|ab=3,2|old_degree=2|"
    "reflections={}|nef=O+(-P1)-F|P1_height=21/2|P1.O=4".format(
        ",".join(map(str, expected_reflections))
    ),
    flush=True,
)
print(
    "H3Q6|nef_pairings=F:2,O:1,minusP1:1,affineE7:1,affineE8:2|"
    "child=E8+E6/MW3|root_data=14,312,3",
    flush=True,
)
print(
    "H3Q6|MW_height={}|optimal_PO={}|q60_MW_height={}|"
    "q60_height_isometric=0|status=PASS".format(
        tuple(tuple(value for value in row) for row in h3_height.rows()),
        ",".join(map(str, pole_profile)),
        tuple(tuple(value for value in row) for row in q60_height.rows()),
    ),
    flush=True,
)
print(
    "H3Q8|q=8|ab=2,4|raw_old_degree=4|raw_O=-2|"
    "fixed_reflections=O:-2|nef_ab=4,2|nef_old_degree=2|nef_O=2|"
    "component_pairings={}|affine_pairings=1,0|"
    "section_and_bisection_nef_proof=1|child=D13/MW4|"
    "root_data=13,312,4|MW_height={}|status=PASS".format(
        ",".join(map(str, simple_pairings))
        , tuple(tuple(value for value in row) for row in q8_height.rows())
    ),
    flush=True,
)
