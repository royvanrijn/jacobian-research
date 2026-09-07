#!/usr/bin/env sage
"""Reduce the fourth q=12 class in the exact third-child chamber.

The first three neighbors are replayed by
``analyze_q80_second_neighbor_chamber.sage``.  This certificate constructs
the generic ``A5+A3+3A1/MW6`` third child, reduces the pinned fourth q=12
isotropic class against its zero section and all affine/simple fiber
components, and records its root and saturated Mordell--Weil projections.

This is an exact lattice/chamber planning certificate.  Nonnegativity on the
displayed old curves is necessary for the equation-level pencil, but is not
by itself a proof of nefness against every effective section or multisection.
"""

from pathlib import Path

from sage.all import (
    QQ,
    ZZ,
    QuadraticForm,
    block_diagonal_matrix,
    identity_matrix,
    lcm,
    matrix,
    vector,
)
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


HERE = Path(__file__).resolve().parent
load(str(HERE / "analyze_q80_second_neighbor_chamber.sage"))


# Construct the third child from the pinned third q12 witness.
third_child_frame, third_child_transport = neighbor(
    second_frame,
    ZZ(third["q"]),
    ZZ(third["a"]),
    ZZ(third["b"]),
    vector(ZZ, map(ZZ, third["v"].split(","))),
)
assert third_child_frame.det() == 948

third_child_ns = block_diagonal_matrix(U, -third_child_frame)
simple, positive = deterministic_simple_roots(third_child_frame)
assert simple.nrows() == 11
root_gram = simple * third_child_frame * simple.transpose()
components = connected_components(root_gram)
assert sorted(map(len, components)) == [1, 1, 1, 3, 5]

fiber = vector(ZZ, [1, 0] + [0] * 17)
zero = vector(ZZ, [-1, 1] + [0] * 17)
curves = [("O", zero)] + [
    (f"R{index}", vector(ZZ, [0, 0] + list(root)))
    for index, root in enumerate(simple.rows(), 1)
]
affine_data = highest_roots(third_child_frame, simple, positive)
for component_index, (_, root, _) in enumerate(affine_data, 1):
    curves.append(
        (
            f"Theta0_{component_index}",
            fiber - vector(ZZ, [0, 0] + list(root)),
        )
    )

fourth = steps[3]
raw = vector(
    ZZ,
    [ZZ(fourth["a"]), ZZ(fourth["b"])]
    + list(map(ZZ, fourth["v"].split(","))),
)
assert raw * third_child_ns * raw == 0
reduced, reflection_sequence = chamber_reduce(raw, curves, third_child_ns)
assert reduced * third_child_ns * reduced == 0
assert all(
    intersection(reduced, curve, third_child_ns) >= 0
    for _, curve in curves
)

frame_part = vector(QQ, reduced[2:])
root_coordinates = (
    frame_part
    * third_child_frame
    * simple.transpose()
    * root_gram.inverse()
)
root_projection = root_coordinates * simple
mw_projection = frame_part - root_projection
mw_norm = mw_projection * third_child_frame * mw_projection

# Saturated MW quotient, directly from the orthogonal projector.
projection = (
    identity_matrix(QQ, 17)
    - third_child_frame
    * simple.transpose()
    * root_gram.inverse()
    * simple
)
projection_denominator = lcm(value.denominator() for value in projection.list())
scaled_projection = (projection_denominator * projection).change_ring(ZZ)
mw_projected_integer = scaled_projection.row_module().basis_matrix()
assert mw_projected_integer.nrows() == 6
mw_basis = mw_projected_integer / projection_denominator
mw_height = mw_basis * third_child_frame * mw_basis.transpose()
height_denominator = lcm(value.denominator() for value in mw_height.list())
mw_transform = (
    (height_denominator * mw_height).change_ring(ZZ).LLL_gram().transpose()
)
assert abs(mw_transform.det()) == 1
mw_basis = mw_transform * mw_basis
mw_height = mw_basis * third_child_frame * mw_basis.transpose()
mw_coordinates = mw_basis.solve_left(mw_projection)
assert mw_coordinates * mw_basis == mw_projection


def integral_preimage(projected):
    """Return an integral frame row with the requested MW projection."""
    target = vector(QQ, projected) * projection_denominator
    assert all(value in ZZ for value in target)
    target = vector(ZZ, target)
    linear_map = scaled_projection.transpose()
    diagonal, left, right = linear_map.smith_form()
    transformed = left * target
    smith_coordinates = vector(ZZ, [0] * 17)
    for index in range(17):
        elementary_divisor = diagonal[index, index]
        if elementary_divisor:
            assert transformed[index] % elementary_divisor == 0
            smith_coordinates[index] = transformed[index] // elementary_divisor
        else:
            assert transformed[index] == 0
    preimage = right * smith_coordinates
    assert vector(QQ, preimage) * projection == vector(QQ, projected)
    return vector(ZZ, preimage)


# A complete wall check is smaller when centered at D.  For an irreducible
# (-2)-curve C=(a,m,l) fixed in |D|, effectiveness of D-C gives
# 1 <= m=C.F <= D.F=3.  Since C^2=-2 and D=(4,3,d), d^2=24,
#
#   D.C = (3/(2m))*||l-(m/3)d||^2 - 3/m.
#
# Therefore D.C<0 forces ||l-(m/3)d||^2<2.  Enumerate the rank-six MW
# projection inside that radius and lift each coset by an exact rank-eleven
# root-lattice CVP.  Every iterator is continued until the first vector at
# distance >=2, so the cap below is fail-closed rather than evidentiary.
divisor_frame = vector(ZZ, reduced[2:])
assert divisor_frame * third_child_frame * divisor_frame == 24
nearby_form = QuadraticForm(ZZ, (24 * mw_height).change_ring(ZZ))
nearby_shells = nearby_form.short_vector_list_up_to_length(
    216, up_to_sign_flag=True
)
nearby_numerators = {tuple([0] * 6)}
for shell in nearby_shells[1:]:
    for row in shell:
        row = vector(ZZ, row)
        nearby_numerators.add(tuple(row))
        nearby_numerators.add(tuple(-row))

root_lattice = IntegralLattice(root_gram)
cvp_cap = 256
cvp_truncated = 0
wall_candidates = {}
for degree in (1, 2, 3):
    center_three = degree * mw_coordinates
    center_frame = QQ(degree) / 3 * divisor_frame
    for numerator_tuple in sorted(nearby_numerators):
        numerator = vector(ZZ, numerator_tuple)
        if any(
            (numerator[index] + center_three[index]) % 3
            for index in range(6)
        ):
            continue
        coordinates = vector(
            ZZ,
            [
                (numerator[index] + center_three[index]) // 3
                for index in range(6)
            ],
        )
        mw_difference = vector(QQ, numerator) / 3
        mw_distance = mw_difference * mw_height * mw_difference
        if mw_distance >= 2:
            continue
        projected = coordinates * mw_basis
        raw_lift = integral_preimage(projected)
        root_difference = vector(QQ, raw_lift) - center_frame
        raw_root_coordinates = (
            root_difference
            * third_child_frame
            * simple.transpose()
            * root_gram.inverse()
        )
        iterator = root_lattice.enumerate_close_vectors(-raw_root_coordinates)
        exhausted_below_two = False
        for _ in range(cvp_cap):
            try:
                shift = vector(ZZ, next(iterator))
            except StopIteration:
                exhausted_below_two = True
                break
            lift = raw_lift + shift * simple
            difference = vector(QQ, lift) - center_frame
            distance = difference * third_child_frame * difference
            if distance >= 2:
                exhausted_below_two = True
                break
            norm = ZZ(lift * third_child_frame * lift)
            if (norm - 2) % (2 * degree):
                continue
            fiber_coefficient = (norm - 2) // (2 * degree)
            curve = vector(
                ZZ, [fiber_coefficient, degree] + list(lift)
            )
            assert curve * third_child_ns * curve == -2
            if any(
                intersection(curve, old_curve, third_child_ns) < 0
                for _, old_curve in curves
            ):
                continue
            pairing = intersection(reduced, curve, third_child_ns)
            wall_candidates[tuple(curve)] = (
                degree, tuple(coordinates), distance, pairing
            )
        if not exhausted_below_two:
            cvp_truncated += 1

assert cvp_truncated == 0
negative_walls = {
    curve: data for curve, data in wall_candidates.items() if data[-1] < 0
}
assert not negative_walls
minimum_wall_pairing = min(
    (data[-1] for data in wall_candidates.values()), default=None
)

# Recover the effective section representing the horizontal MW projection.
# Its pole count determines the generic decomposition of the trisection and
# measures the pole drop at the CM24 boundary.
raw_section_lift = integral_preimage(mw_projection)
raw_section_root_coordinates = (
    vector(QQ, raw_section_lift)
    * third_child_frame
    * simple.transpose()
    * root_gram.inverse()
)
section_iterator = root_lattice.enumerate_close_vectors(
    -raw_section_root_coordinates
)
effective_sections = []
first_effective_norm = None
for _ in range(1024):
    shift = vector(ZZ, next(section_iterator))
    lift = raw_section_lift + shift * simple
    norm = ZZ(lift * third_child_frame * lift)
    if first_effective_norm is not None and norm > first_effective_norm:
        break
    if norm < 4 or norm % 2:
        continue
    pole = norm // 2 - 2
    section = vector(ZZ, [pole + 1, 1] + list(lift))
    if any(
        intersection(section, old_curve, third_child_ns) < 0
        for _, old_curve in curves[1:]
    ):
        continue
    if first_effective_norm is None:
        first_effective_norm = norm
    effective_sections.append(section)
assert len(effective_sections) == 1
horizontal_section = effective_sections[0]
horizontal_pole = intersection(horizontal_section, zero, third_child_ns)
horizontal_pairing = intersection(reduced, horizontal_section, third_child_ns)
vertical_remainder = reduced - horizontal_section - 2*zero
fiber_shift = ZZ(vertical_remainder[0])
vertical_remainder -= fiber_shift*fiber
assert vertical_remainder[:2] == vector(ZZ, (0, 0))
integral_vertical_coordinates = simple.solve_left(vertical_remainder[2:])
assert all(value in ZZ for value in integral_vertical_coordinates)

component_pairings = tuple(
    (name, intersection(reduced, curve, third_child_ns))
    for name, curve in curves
)
vertical_components = tuple(
    (
        component,
        tuple(-root_coordinates[index] for index in component),
    )
    for component in components
)
affine_multiplicities = tuple(
    (
        component,
        tuple(coordinates[index] for index in component),
    )
    for component, _, coordinates in affine_data
)

print(
    "Q80FOURTHQ12CHAMBER|source=A5+A3+3A1/MW6|q=12|ab=3,4|"
    f"raw={tuple(raw)}|raw_D.F={intersection(raw, fiber, third_child_ns)}|"
    f"raw_D.O={intersection(raw, zero, third_child_ns)}",
    flush=True,
)
print(
    f"Q80FOURTHQ12CHAMBER|reflection_count={len(reflection_sequence)}|"
    f"reduction={reflection_sequence}|reduced={tuple(reduced)}|"
    f"D.F={intersection(reduced, fiber, third_child_ns)}|"
    f"D.O={intersection(reduced, zero, third_child_ns)}|"
    f"component_pairings={component_pairings}",
    flush=True,
)
print(
    f"Q80FOURTHQ12CHAMBER|root_coordinates={tuple(root_coordinates)}|"
    f"vertical_components={vertical_components}|"
    f"affine_multiplicities={affine_multiplicities}",
    flush=True,
)
print(
    f"Q80FOURTHQ12CHAMBER|mw_projection={tuple(mw_projection)}|"
    f"mw_norm={mw_norm}|mw_coordinates={tuple(mw_coordinates)}|"
    f"mw_height={tuple(tuple(row) for row in mw_height.rows())}|"
    f"mw_height_det={mw_height.det()}|"
    f"nearby_MW_numerators={len(nearby_numerators)}|"
    f"effective_wall_candidates={len(wall_candidates)}|"
    f"minimum_wall_pairing={minimum_wall_pairing}|"
    f"cvp_cap={cvp_cap}|cvp_truncated={cvp_truncated}|"
    "old_curves_nonnegative=1|all_effective_minus2_nonnegative=1|nef=1|"
    "status=PASS_NEF_CHAMBER",
    flush=True,
)
print(
    f"Q80FOURTHQ12CHAMBER|horizontal_section={tuple(horizontal_section)}|"
    f"section_P.O={horizontal_pole}|D.S={horizontal_pairing}|"
    f"decomposition=S+2O+{fiber_shift}F+root_correction|"
    f"integral_root_coordinates={tuple(integral_vertical_coordinates)}|"
    "generic_fiber_RR=<1,x,(y+yS)/(x-xS)>|"
    "status=PASS_HORIZONTAL_SECTION",
    flush=True,
)
