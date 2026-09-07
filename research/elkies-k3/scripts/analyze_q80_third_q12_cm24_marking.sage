#!/usr/bin/env sage
"""Identify the horizontal section in the third q=12 pencil at CM24.

The first two q=4 pencils are already explicit.  This script transports the
pinned third fiber through their exact lattice frames after adjoining the
CM24 algebraic class.  Its Mordell--Weil projection has height three and
coordinates ``(2,-1)`` in the specialized second child's reduced MW basis.

The same second child has two explicit polynomial sections over
``K=QQ(sqrt(-6))``.  Their I2 component profiles identify the abstract MW
basis, proving that the horizontal point of the q12 trisection is, up to
sign, the polynomial section ``Q=P1+3*P2`` displayed below.  This is a
marking certificate, not yet an execution of the seven remaining vertical
Riemann--Roch gates or of the generic q12 pencil.
"""

from pathlib import Path

from sage.all import (
    EllipticCurve,
    PolynomialRing,
    QQ,
    QuadraticField,
    ZZ,
    block_diagonal_matrix,
    matrix,
    vector,
)


HERE = Path(__file__).resolve().parent

# The first loader supplies the exact CM24 extension of the q80 frame and
# saturated-MW helpers.  The second supplies the pinned q4,q4,q12 transports
# and deterministic chamber routines.
load(str(HERE / "classify_kumar_cm_frame_extensions.sage"))
special_q80 = q80_cm24
generic_to_special = q80_embedding24
load(str(HERE / "analyze_q80_second_neighbor_chamber.sage"))


def enhance_neighbor(transport, embedding, special_frame):
    """Transport one generic neighbor basis after a rank-one NS extension."""
    lifted = [
        vector(ZZ, list(row[:2]) + list(vector(ZZ, row[2:]) * embedding))
        for row in transport.rows()
    ]
    special_ns = block_diagonal_matrix(U, -special_frame)
    fiber, mate = lifted[:2]
    assert fiber * special_ns * fiber == 0
    assert mate * special_ns * mate == 0
    assert fiber * special_ns * mate == 1
    complement = matrix(
        ZZ, [list(fiber * special_ns), list(mate * special_ns)]
    ).right_kernel_matrix()
    basis = matrix(
        ZZ, [list(fiber), list(mate)] + [list(row) for row in complement]
    )
    assert abs(basis.det()) == 1
    child = -(complement * special_ns * complement.transpose())
    inverse_basis = basis.inverse()
    embedding_rows = []
    for row in lifted[2:]:
        coordinates = row * inverse_basis
        assert coordinates[0] == coordinates[1] == 0
        embedding_rows.append(list(coordinates[2:]))
    child_embedding = matrix(ZZ, embedding_rows)
    return child, child_embedding


# Follow the two pinned q4 moves in the generic and specialized lattices.
first_step = steps[0]
assert first_frame.det() == 948
special_first, first_embedding = enhance_neighbor(
    first_transport, generic_to_special, special_q80
)

second_step = steps[1]
second_frame, second_transport = neighbor(
    first_frame,
    ZZ(second_step["q"]),
    ZZ(second_step["a"]),
    ZZ(second_step["b"]),
    vector(ZZ, map(ZZ, second_step["v"].split(","))),
)
special_second, second_embedding = enhance_neighbor(
    second_transport, first_embedding, special_first
)
assert second_frame == second_embedding * special_second * second_embedding.transpose()

# Reduce the transported q12 divisor in the full CM24 chamber.
third_step = steps[2]
generic_third = vector(
    ZZ,
    [ZZ(third_step["a"]), ZZ(third_step["b"])]
    + list(map(ZZ, third_step["v"].split(","))),
)
generic_third_frame, third_transport = neighbor(
    second_frame,
    ZZ(third_step["q"]),
    ZZ(third_step["a"]),
    ZZ(third_step["b"]),
    vector(ZZ, map(ZZ, third_step["v"].split(","))),
)
special_third, third_embedding = enhance_neighbor(
    third_transport, second_embedding, special_second
)
assert (
    generic_third_frame
    == third_embedding * special_third * third_embedding.transpose()
)
assert root_components(special_third) == [
    (1, 2, 2),
    (1, 2, 2),
    (1, 2, 2),
    (6, 42, 7),
    (6, 42, 7),
]
third = vector(
    ZZ,
    list(generic_third[:2])
    + list(vector(ZZ, generic_third[2:]) * second_embedding),
)
special_ns = block_diagonal_matrix(U, -special_second)
simple, positive = deterministic_simple_roots(special_second)
fiber = vector(ZZ, [1, 0] + [0] * 18)
zero = vector(ZZ, [-1, 1] + [0] * 18)
curves = [("O", zero)] + [
    (f"R{index}", vector(ZZ, [0, 0] + list(root)))
    for index, root in enumerate(simple.rows(), 1)
]
for component_index, (_, root, _) in enumerate(
    highest_roots(special_second, simple, positive), 1
):
    curves.append(
        (
            f"Theta0_{component_index}",
            fiber - vector(ZZ, [0, 0] + list(root)),
        )
    )
reduced, reflection_sequence = chamber_reduce(third, curves, special_ns)
assert intersection(reduced, fiber, special_ns) == 3
assert intersection(reduced, zero, special_ns) == 0

root_gram = simple * special_second * simple.transpose()
special_frame_part = vector(QQ, reduced[2:])
special_root_coordinates = (
    special_frame_part
    * special_second
    * simple.transpose()
    * root_gram.inverse()
)
special_components = connected_components(root_gram)
special_affine_data = highest_roots(special_second, simple, positive)
special_vertical_components = tuple(
    (
        component,
        tuple(-special_root_coordinates[index] for index in component),
    )
    for component in special_components
)
special_affine_multiplicities = tuple(
    (
        component,
        tuple(coordinates[index] for index in component),
    )
    for component, _, coordinates in special_affine_data
)


def mw_projection(row):
    row = vector(QQ, row)
    return (
        row
        - row
        * special_second
        * simple.transpose()
        * root_gram.inverse()
        * simple
    )


reduced_mw = mw_projection(reduced[2:])
root_data = root_invariants(special_second)
assert root_components(special_second) == [
    (1, 2, 2),
    (1, 2, 2),
    (1, 2, 2),
    (6, 72, 3),
    (7, 84, 4),
]
_, mw_height, mw_lifts = mw_height_gram(
    special_second, root_data[3], return_lifts=True
)
assert mw_height == matrix(QQ, [[QQ(5) / 12, -QQ(1) / 6], [-QQ(1) / 6, QQ(2) / 3]])
optimal = optimal_section_pole_basis(special_second, mw_height, mw_lifts)
assert optimal == (0, (0, 0), ((1, 0), (0, 1)))
mw_basis_lifts = matrix(ZZ, optimal[2]) * mw_lifts
projected_basis = matrix(
    QQ, [list(mw_projection(row)) for row in mw_basis_lifts.rows()]
)
mw_coordinates = projected_basis.solve_left(reduced_mw)
assert mw_coordinates == vector(QQ, (2, -1))
assert reduced_mw * special_second * reduced_mw == 3

# The generic marked section has S.O=2, but its CM24 specialization has
# Q.O=0.  Thus the same nef divisor acquires two extra full old fibers:
# D=Q+2O+4F+R.  Unlike the orthogonal root projection, the correction R is
# integral and is the equation-level component gate needed below.
q_lift = vector(ZZ, mw_coordinates * mw_basis_lifts)
q_root_coordinates = (
    vector(QQ, q_lift)
    * special_second
    * simple.transpose()
    * root_gram.inverse()
)
effective_q_lifts = []
for shift in IntegralLattice(root_gram).enumerate_close_vectors(
    -q_root_coordinates
):
    candidate_lift = q_lift + vector(ZZ, shift) * simple
    candidate_norm = candidate_lift * special_second * candidate_lift
    if candidate_norm > 4 and effective_q_lifts:
        break
    if candidate_norm != 4:
        continue
    candidate_section = vector(ZZ, [1, 1] + list(candidate_lift))
    component_pairings = tuple(
        intersection(candidate_section, curve, special_ns)
        for _, curve in curves[1:]
    )
    if min(component_pairings) < 0:
        continue
    # Q=P1+3P2 meets the two conjugate quadratic I2 components and the
    # identity component of the rational I2 fiber.
    if component_pairings[:3] != (1, 1, 0):
        continue
    if component_pairings[len(simple.rows()):][:3] != (0, 0, 1):
        continue
    effective_q_lifts.append((candidate_lift, component_pairings))
assert len(effective_q_lifts) == 1
q_lift, q_component_pairings = effective_q_lifts[0]
q_section = vector(ZZ, [1, 1] + list(q_lift))
assert q_section * special_ns * q_section == -2
assert intersection(q_section, zero, special_ns) == 0
special_vertical_root = (
    reduced - q_section - 2 * zero - 4 * fiber
)
assert special_vertical_root[:2] == vector(ZZ, (0, 0))
special_integral_root_coordinates = simple.solve_left(
    vector(ZZ, special_vertical_root[2:])
)
assert all(value in ZZ for value in special_integral_root_coordinates)
special_integral_vertical_components = tuple(
    (
        component,
        tuple(
            -ZZ(special_integral_root_coordinates[index])
            for index in component
        ),
    )
    for component in special_components
)
assert tuple(coefficients for _, coefficients in special_integral_vertical_components) == (
    (0,),
    (1,),
    (0,),
    (1, 2, 3, 2, 1, 2),
    (2, 4, 3, 3, 5, 6, 3),
)
d7_component = special_components[-1]
d7_cartan = root_gram.matrix_from_rows_and_columns(
    d7_component, d7_component
)
d7_correction = vector(
    ZZ, special_integral_vertical_components[-1][1]
)
d7_boundary = d7_cartan * d7_correction

# The chosen abstract basis has A1 profiles (0,0,1) and (1,1,0).
component_data = root_component_data(special_second)
abstract_a1_profiles = []
for invariant, component_basis in component_data[:3]:
    assert invariant == (1, 2, 2)
    classes = tuple(
        fractional_root_class(special_second, lift, component_basis)
        for lift in mw_basis_lifts.rows()
    )
    abstract_a1_profiles.append(tuple(int(point != (QQ(0),)) for point in classes))
assert tuple(abstract_a1_profiles) == ((0, 1), (0, 1), (1, 0))

# Explicit CM24 second-child model and its two polynomial sections.
K = QuadraticField(-6, "s")
s = K.gen()
function_field = PolynomialRing(K, "W").fraction_field()
W = function_field.gen()
A2 = (
    -QQ(10460353203) / 64
    + QQ(1162261467) / 32 * W
    + QQ(129140163) / 8 * W**2
    + QQ(13286025) / 8 * W**3
    + QQ(59049) * W**4
    - 27 * W**6
)
B2 = (
    -QQ(5147278302366225) / 512
    - QQ(1303977169932777) / 256 * W
    - QQ(144886352214753) / 128 * W**2
    - QQ(4487491524087) / 32 * W**3
    - QQ(331244518095) / 32 * W**4
    - QQ(7360989291) / 16 * W**5
    - QQ(97253703) / 8 * W**6
    - QQ(177147) * W**7
    + 54 * W**9
)
curve = EllipticCurve(function_field, [0, 0, 0, A2, B2])
P1 = curve(
    3 * W**3 + QQ(6561) / 4 * W + QQ(59049) / 8,
    s
    * (
        QQ(2187) / 2 * W**3
        + QQ(177147) / 8 * W**2
        + QQ(4782969) / 16 * W
        + QQ(43046721) / 32
    ),
)
P2 = curve(
    3 * W**3 - QQ(6561) / 4 * W - QQ(59049) / 4,
    s
    * (
        QQ(2187) / 2 * W**3
        + QQ(295245) / 8 * W**2
        + QQ(1594323) / 4 * W
        + QQ(43046721) / 32
    ),
)
assert P1 and P2 and P1 != P2 and P1 != -P2


def to_polynomial(value, ring):
    value = function_field(value)
    assert value.denominator() == 1
    return ring(value.numerator().list())


polynomial_ring = PolynomialRing(K, "u")
u = polynomial_ring.gen()
polynomial_A2 = to_polynomial(A2, polynomial_ring)
polynomial_B2 = to_polynomial(B2, polynomial_ring)
i2_factors = (
    u + QQ(27) / 4,
    u**2 + QQ(27) / 2 * u + QQ(729) / 4,
)


def i2_profile(point):
    labels = []
    for factor in i2_factors:
        quotient = polynomial_ring.quotient(factor, "z")
        local_A = quotient(polynomial_A2)
        local_B = quotient(polynomial_B2)
        node_x = -3 * local_B / (2 * local_A)
        point_x = quotient(point[0].numerator().list()) / quotient(
            point[0].denominator().list()
        )
        labels.append(int(point_x == node_x))
    return tuple(labels)


# Tuple order is the rational I2 and then the conjugate quadratic pair.
assert i2_profile(P1) == (1, 1)
assert i2_profile(P2) == (1, 0)
assert i2_profile(P1 + P2) == (0, 1)

# Hence e1=+/-P2 and e2=-/+ (P1+P2).  The abstract vector (2,-1)
# is therefore +/- (P1+3P2).
Q12 = P1 + 3 * P2
expected_qx = (
    -QQ(8) / 27 * W**4
    + 22 * W**3
    - QQ(243) / 2 * W**2
    + 729 * W
    - QQ(492075) / 8
)
expected_qy = s * (
    QQ(16) / 243 * W**6
    - QQ(22) / 3 * W**5
    + QQ(333) / 2 * W**4
    - QQ(2025) / 4 * W**3
    + QQ(190269) / 4 * W**2
    - QQ(177147) / 16 * W
    + QQ(199290375) / 32
)
assert Q12 == curve(expected_qx, expected_qy)
assert i2_profile(Q12) == (0, 1)
assert Q12[0].denominator() == Q12[1].denominator() == 1
assert (Q12[0].numerator().degree(), Q12[1].numerator().degree()) == (4, 6)

print(
    "Q80THIRDQ12CM24|special_second=D7+E6+3A1|MW=2|"
    "height=((5/12,-1/6),(-1/6,2/3))|basis_PO=0,0",
    flush=True,
)
print(
    f"Q80THIRDQ12CM24|special_third_components={root_components(special_third)}|"
    f"special_third_root_rank={root_invariants(special_third)[0]}",
    flush=True,
)
print(
    f"Q80THIRDQ12CM24|old_vertical_components={special_vertical_components}|"
    f"old_affine_multiplicities={special_affine_multiplicities}",
    flush=True,
)
print(
    "Q80THIRDQ12CM24|decomposition=Q+2O+4F+root_correction|"
    f"Q_component_pairings={q_component_pairings}|"
    f"integral_vertical_components={special_integral_vertical_components}",
    flush=True,
)
print(
    f"Q80THIRDQ12CM24|D7_cartan={tuple(tuple(row) for row in d7_cartan.rows())}|"
    f"D7_correction={tuple(d7_correction)}|"
    f"D7_Cc={tuple(d7_boundary)}",
    flush=True,
)
print(
    f"Q80THIRDQ12CM24|third_DF=3|third_DO=0|"
    f"reflections={len(reflection_sequence)}|MW_coordinates={tuple(mw_coordinates)}|"
    "MW_height=3",
    flush=True,
)
print(
    "Q80THIRDQ12CM24|field=Q(sqrt(-6))|"
    "basis=e1:+/-P2,e2:-/+(P1+P2)|Q12=+/-(P1+3P2)|"
    "Q12_polynomial_degrees=4,6|I2_profile=quadratic_pair",
    flush=True,
)
print(
    f"Q80THIRDQ12CM24|Q12_x={expected_qx}|Q12_y={expected_qy}",
    flush=True,
)
# Compare the bounded alternate q=8 move in the same exact marking.  This is
# a genuinely different generic neighbor (E6+A7/MW4), not one of the A13
# return loops.  At CM24 we transport both its U embedding and its divisor
# into the full rank-20 chamber so that an equation-level use cannot confuse
# it with the CM43 q=8 class that collapses to the old fiber.
special_q8, q8_embedding = enhance_neighbor(
    alternate_transport, second_embedding, special_second
)
special_q8_raw = vector(
    ZZ,
    [2, 4]+list(vector(ZZ, alternate_q8)*second_embedding),
)
special_q8_reduced, special_q8_sequence = chamber_reduce(
    special_q8_raw, curves, special_ns
)
assert intersection(special_q8_reduced, fiber, special_ns) == 2
assert intersection(special_q8_reduced, zero, special_ns) == 0
assert all(
    intersection(special_q8_reduced, curve, special_ns) >= 0
    for _, curve in curves
)
special_q8_component_pairings = tuple(
    (name, intersection(special_q8_reduced, curve, special_ns))
    for name, curve in curves[1:]
)
special_q8_components = root_components(special_q8)
print(
    f"Q80SECONDCHILDQ8CM24|raw={tuple(special_q8_raw)}|"
    f"reduction={special_q8_sequence}|reduced={tuple(special_q8_reduced)}|"
    f"zero={tuple(zero)}|component_pairings={special_q8_component_pairings}|"
    f"components={special_q8_components}|"
    f"root_rank={root_invariants(special_q8)[0]}|"
    f"MW={18-root_invariants(special_q8)[0]}|D.F=2|D.O=0|"
    "rho=20|trivial_rank=18|shioda_tate_MW=2|"
    "status=PASS_SPECIALIZED_MARKING",
    flush=True,
)
print("Q80THIRDQ12CM24|status=PASS_MARKED_HORIZONTAL", flush=True)
