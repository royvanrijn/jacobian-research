#!/usr/bin/env sage
"""Mark the final q6 chord on the exact CM24 fifth child.

The productive fifth q4 class is the generic degree-47/pair23 neighbor.  Its
CM24 child has roots ``D4+3A3+3A1`` and MW rank two.  This checker applies the
same source-chamber Weyl word to the full specialized neighbor basis, then
transports the already certified generic q6 horizontal section into that
aligned CM24 child.

It reports the section's pole/component/MW data, proves the specialized
divisor nef, computes the specialized q6 child root data, and decomposes
``D-O-S`` into its local affine-component cycles.  This is the exact lattice
marking needed before constructing the saturated degree-two local module on
the explicit pair23 secant-gauge Jacobian; it does not yet construct that
function.
"""

from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, matrix, vector
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


HERE = Path(__file__).resolve().parent
load(str(HERE / "analyze_q80_fifth_q4_cm24_readiness.sage"))


alternate_q4_v = vector(ZZ, (
    -9, 8, -11, 10, -4, 0, 5, 1, -6,
    6, 1, -2, -1, -1, 1, 2, 0,
))
alternate_child, alternate_transport = neighbor(
    generic_fourth_frame, ZZ(4), ZZ(2), ZZ(2), alternate_q4_v
)
assert root_invariants(alternate_child)[:3] == (1, 2, 2)

source_ns = block_diagonal_matrix(U, -special_fourth)


def full_specialized_neighbor_basis(transport, embedding):
    """Return the full specialized U+frame basis, not only its embedding."""
    lifted = [
        vector(ZZ, list(row[:2]) + list(vector(ZZ, row[2:]) * embedding))
        for row in transport.rows()
    ]
    neighbor_fiber, neighbor_mate = lifted[:2]
    complement = matrix(
        ZZ,
        [
            list(neighbor_fiber * source_ns),
            list(neighbor_mate * source_ns),
        ],
    ).right_kernel_matrix()
    basis = matrix(
        ZZ,
        [list(neighbor_fiber), list(neighbor_mate)]
        + [list(row) for row in complement.rows()],
    )
    assert abs(basis.det()) == 1
    return basis


# Construct the full CM24 neighbor basis in the source fourth-child NS.  The
# first basis row is the raw specialization of the generic q4 fiber.
raw_basis = full_specialized_neighbor_basis(
    alternate_transport, fourth_embedding
)
special_raw_fiber = vector(
    ZZ, [2, 2] + list(alternate_q4_v * fourth_embedding)
)
assert raw_basis[0] == special_raw_fiber

# Reduce the specialized fiber and apply that identical Weyl word to every
# basis row.  This aligns the abstract special child with the explicit
# pair23 secant-gauge pencil rather than with the non-nef raw marking.
reduced_fiber, reflection_word = chamber_reduce(
    special_raw_fiber, new_curves, source_ns
)
curve_by_name = dict(new_curves)


def apply_source_word(row):
    row = vector(ZZ, row)
    for name, _ in reflection_word:
        curve = curve_by_name[name]
        row += intersection(row, curve, source_ns) * curve
    return row


aligned_basis = matrix(
    ZZ, [list(apply_source_word(row)) for row in raw_basis.rows()]
)
assert aligned_basis[0] == reduced_fiber
special_child_ns = (
    aligned_basis * source_ns * aligned_basis.transpose()
)
assert special_child_ns[:2, :2] == U
special_child = -special_child_ns[2:, 2:]
assert special_child_ns == block_diagonal_matrix(U, -special_child)
assert tuple(map(int, root_invariants(special_child)[:3])) == (16, 66, 2048)
aligned_inverse = aligned_basis.inverse().change_ring(ZZ)

# The final generic q6 chamber has D=O+S-F, with the horizontal section
# S=(5,1,v) in the generic fifth-child basis.
q6_v = vector(ZZ, (
    0, -2, 4, 2, -1, 2, 1, -1, 1,
    0, 1, -1, 1, 0, 0, 0, 0,
))
generic_section = vector(ZZ, [5, 1] + list(q6_v))
assert generic_section * block_diagonal_matrix(U, -alternate_child) * generic_section == -2
generic_divisor = vector(ZZ, [3, 2] + list(q6_v))
assert generic_divisor == vector(ZZ, [-1, 1] + [0] * 17) + generic_section - vector(ZZ, [1, 0] + [0] * 17)
assert generic_divisor * block_diagonal_matrix(U, -alternate_child) * generic_divisor == 0

# Move S back to the generic fourth child, specialize the frame coordinates,
# apply the CM chamber word, and express it in the aligned special child.
section_generic_source = vector(ZZ, generic_section * alternate_transport)
section_special_source = vector(
    ZZ,
    list(section_generic_source[:2])
    + list(vector(ZZ, section_generic_source[2:]) * fourth_embedding),
)
section_special_source = apply_source_word(section_special_source)
special_section_raw = vector(ZZ, section_special_source * aligned_inverse)

divisor_generic_source = vector(ZZ, generic_divisor * alternate_transport)
divisor_special_source = vector(
    ZZ,
    list(divisor_generic_source[:2])
    + list(vector(ZZ, divisor_generic_source[2:]) * fourth_embedding),
)
divisor_special_source = apply_source_word(divisor_special_source)
special_divisor_raw = vector(ZZ, divisor_special_source * aligned_inverse)

fiber = vector(ZZ, [1, 0] + [0] * 18)
zero = vector(ZZ, [-1, 1] + [0] * 18)
assert special_section_raw * special_child_ns * special_section_raw == -2
assert intersection(special_section_raw, fiber, special_child_ns) == 1

simple, positive = deterministic_simple_roots(special_child)
root_gram = simple * special_child * simple.transpose()
curves = [("O", zero)] + [
    (f"R{index}", vector(ZZ, [0, 0] + list(root)))
    for index, root in enumerate(simple.rows(), 1)
]
for component_index, (_, highest, _) in enumerate(
    highest_roots(special_child, simple, positive), 1
):
    curves.append((
        f"Theta0_{component_index}",
        fiber - vector(ZZ, [0, 0] + list(highest)),
    ))

# Specialization can add vertical fixed components to a generic section.
# Reflect them away inside the child fiber chamber; these reflections preserve
# the fiber degree and recover the actual effective CM section.
special_section = vector(ZZ, special_section_raw)
section_reflections = []
while True:
    for name, curve in curves[1:]:
        pairing = intersection(special_section, curve, special_child_ns)
        if pairing < 0:
            special_section += pairing * curve
            section_reflections.append((name, pairing))
            assert special_section * special_child_ns * special_section == -2
            break
    else:
        break
section_reflections = tuple(section_reflections)
component_pairings = tuple(
    (name, intersection(special_section, curve, special_child_ns))
    for name, curve in curves
)
assert all(value >= 0 for _, value in component_pairings)
section_pole = intersection(special_section, zero, special_child_ns)

frame_part = vector(QQ, special_section[2:])
root_coordinates = (
    frame_part * special_child * simple.transpose() * root_gram.inverse()
)
projection = frame_part - root_coordinates * simple
section_height = projection * special_child * projection

# The unspecialized relative ambient for the final local module is
# L=O+S+2F.  It has square six and is nonnegative on every displayed fiber
# component/section, so its K3 Riemann--Roch dimension is five.  The three
# disjoint root-cycle subtractions below cut it to the isotropic pencil.
ambient_divisor = zero + special_section + 2*fiber
assert ambient_divisor * special_child_ns * ambient_divisor == 6
ambient_pairings = tuple(
    (name, intersection(ambient_divisor, curve, special_child_ns))
    for name, curve in curves
)
assert all(value >= 0 for _, value in ambient_pairings)
ambient_frame = vector(ZZ, ambient_divisor[2:])
assert ambient_frame * special_child * ambient_frame == 6
ambient_lattice = IntegralLattice(special_child)
ambient_closest = vector(ZZ, next(ambient_lattice.enumerate_close_vectors(
    vector(QQ, ambient_frame)/2
)))
ambient_closest_distance = (
    ambient_closest-vector(QQ, ambient_frame)/2
) * special_child * (
    ambient_closest-vector(QQ, ambient_frame)/2
)
ambient_minimum_section_pairing = ambient_closest_distance-QQ(1)/2
assert ambient_minimum_section_pairing >= 0
# For a degree-two (-2)-curve the corresponding completed-square formula is
# L.C=||w-v||^2/2+2, so bisections are automatically nonnegative.  Any fixed
# (-2)-curve has old-fiber degree at most L.F=2; components, sections, and
# bisections therefore exhaust the walls and prove L nef.  Since L^2=6, K3
# Riemann--Roch gives h0(L)=2+L^2/2=5.
ambient_minimum_bisection_lower_bound = ZZ(2)

root_data = root_invariants(special_child)
_, mw_height, mw_lifts = mw_height_gram(
    special_child, root_data[3], return_lifts=True
)
optimal = optimal_section_pole_basis(special_child, mw_height, mw_lifts)
mw_basis_lifts = matrix(ZZ, optimal[2]) * mw_lifts


def project_mw(row):
    row = vector(QQ, row)
    return row - row * special_child * simple.transpose() * root_gram.inverse() * simple


projected_basis = matrix(
    QQ, [list(project_mw(row)) for row in mw_basis_lifts.rows()]
)
optimal_height = projected_basis * special_child * projected_basis.transpose()
mw_coordinates = projected_basis.solve_left(projection)
assert mw_coordinates * projected_basis == projection
basis_section_data = tuple(
    section_data_for_lift(
        special_child,
        optimal_height[index, index],
        vector(ZZ, mw_basis_lifts[index]),
    )
    for index in range(mw_basis_lifts.nrows())
)

divisor, divisor_reflections = chamber_reduce(
    special_divisor_raw, curves, special_child_ns
)
assert divisor * special_child_ns * divisor == 0
assert intersection(divisor, fiber, special_child_ns) == 2

# The old zero meets the reduced q6 fiber once, so it is a section of the
# final fibration.  This gives an integral U directly and pins the expected
# CM24 child root data for the equation-level chord classification.
assert intersection(divisor, zero, special_child_ns) == 1
q6_mate = divisor + zero
q6_complement = matrix(
    ZZ,
    [list(divisor * special_child_ns), list(q6_mate * special_child_ns)],
).right_kernel_matrix()
q6_basis = matrix(
    ZZ,
    [list(divisor), list(q6_mate)]
    + [list(row) for row in q6_complement.rows()],
)
assert abs(q6_basis.det()) == 1
q6_ns = q6_basis * special_child_ns * q6_basis.transpose()
assert q6_ns[:2, :2] == U
q6_frame = -q6_ns[2:, 2:]
q6_root_data = tuple(map(int, root_invariants(q6_frame)[:3]))

# Record how the CM-effective section and divisor differ vertically.  The
# generic identity is O+S-F; specialization and component reflections can add
# a root correction and change the whole-fiber coefficient.
vertical_difference = divisor - zero - special_section
assert vertical_difference[1] == 0
try:
    vertical_root_coordinates = simple.solve_left(
        vector(ZZ, vertical_difference[2:])
    )
except ValueError:
    vertical_root_coordinates = None
assert vertical_root_coordinates is not None

# Compute the least local affine-fiber copy needed to make the correction at
# each reducible fiber effective.  Their sum can exceed the global F
# coefficient: here it is three versus two, proving that D-O-S is not itself
# effective and that the constant function is absent from the saturated local
# module.
vertical_components = []
allocated_fibers = 0
for component_index, (component, _, highest_coordinates) in enumerate(
    highest_roots(special_child, simple, positive), 1
):
    coefficients = vector(ZZ, [
        vertical_root_coordinates[index] for index in component
    ])
    marks = vector(ZZ, [highest_coordinates[index] for index in component])
    fiber_copies = max(
        [ZZ(0)]
        + [(-coefficient+mark-1)//mark
           for coefficient, mark in zip(coefficients, marks)]
    )
    effective = coefficients+fiber_copies*marks
    assert all(value >= 0 for value in effective)
    component_gram = root_gram.matrix_from_rows_and_columns(component, component)
    component_root_data = tuple(map(int, root_invariants(component_gram)[:3]))
    vertical_components.append((
        component_index,
        tuple(map(lambda value: int(value)+1, component)),
        component_root_data,
        tuple(map(int, marks)),
        int(fiber_copies),
        tuple(map(int, coefficients)),
        tuple(map(int, effective)),
    ))
    allocated_fibers += fiber_copies
vertical_components = tuple(vertical_components)
assert vertical_difference[0] == 2
assert allocated_fibers == 3

# A negative fixed curve has degree at most two.  As in the generic proof,
# parity excludes bisections.  Sections are an exact CVP in the rank-18 CM
# frame: D.C=(w-v/2)^2-2 in the aligned U+(-M) coordinates.
divisor_frame = vector(ZZ, divisor[2:])
assert divisor_frame * special_child * divisor_frame == 4 * divisor[0]
lattice = IntegralLattice(special_child)
closest = vector(ZZ, next(lattice.enumerate_close_vectors(
    vector(QQ, divisor_frame) / 2
)))
closest_distance = (
    closest - vector(QQ, divisor_frame) / 2
) * special_child * (
    closest - vector(QQ, divisor_frame) / 2
)
minimum_section_pairing = closest_distance - 2

negative_sections = []
if minimum_section_pairing < 0:
    # Record every closest negative class at this distance; the finite loop
    # terminates as soon as the CVP distance increases.
    iterator = lattice.enumerate_close_vectors(vector(QQ, divisor_frame) / 2)
    for _ in range(65536):
        candidate = vector(ZZ, next(iterator))
        distance = (
            candidate - vector(QQ, divisor_frame) / 2
        ) * special_child * (
            candidate - vector(QQ, divisor_frame) / 2
        )
        if distance > closest_distance:
            break
        candidate_section = vector(ZZ, [
            (candidate * special_child * candidate - 2) // 2,
            1,
            *candidate,
        ])
        negative_sections.append((
            tuple(candidate_section),
            intersection(divisor, candidate_section, special_child_ns),
        ))

assert divisor_frame * special_child * divisor_frame % 4 == 0
negative_bisection_impossible = True

print(
    "Q80FINALQ6CM24|"
    f"reflection_count={len(reflection_word)}|"
    f"section_reflection_count={len(section_reflections)}|"
    f"special_child_roots={tuple(map(int, root_invariants(special_child)[:3]))}|"
    f"section={tuple(special_section)}|section_P.O={section_pole}|"
    f"component_pairings={component_pairings}|"
    f"section_height={section_height}|mw_coordinates={tuple(mw_coordinates)}|"
    f"ambient_L_square=6|ambient_h0=5|ambient_pairings={ambient_pairings}|"
    f"ambient_closest_distance={ambient_closest_distance}|"
    f"ambient_minimum_section_pairing={ambient_minimum_section_pairing}|"
    f"ambient_minimum_bisection_lower_bound={ambient_minimum_bisection_lower_bound}|"
    f"mw_height={tuple(tuple(row) for row in optimal_height.rows())}|"
    f"basis_section_data={basis_section_data}|"
    "status=PASS_MARKING",
    flush=True,
)
print(
    "Q80FINALQ6CM24|"
    f"divisor_reflections={divisor_reflections}|"
    f"divisor={tuple(divisor)}|"
    f"q6_CM24_root_data={q6_root_data}|q6_CM24_MW={18-q6_root_data[0]}|"
    f"vertical_fiber_coefficient={vertical_difference[0]}|"
    f"allocated_component_fibers={allocated_fibers}|"
    f"vertical_root_coordinates={None if vertical_root_coordinates is None else tuple(vertical_root_coordinates)}|"
    f"vertical_components={vertical_components}|"
    f"closest_section_vector={tuple(closest)}|"
    f"closest_distance={closest_distance}|"
    f"minimum_section_pairing={minimum_section_pairing}|"
    f"negative_sections={tuple(negative_sections)}|"
    f"negative_bisection_impossible={int(negative_bisection_impossible)}|"
    f"cm_nef={int(minimum_section_pairing >= 0)}|"
    "status=PASS_CM_FIXED_COMPONENT_AUDIT",
    flush=True,
)
