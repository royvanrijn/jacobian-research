#!/usr/bin/env sage
"""Audit fifth-q4 readiness on the compact CM24 fourth child.

The generic fifth neighbor starts from ``4A1/MW13``.  This script transports
its divisor to the Picard-rank-20 CM24 fourth child, performs the exact old
component reduction, and extracts its horizontal section in the saturated
rank-two MW quotient.  It then transports the already known polynomial
sections of the third child through the fourth neighbor and checks whether
the five nonconstant degree-one sections found on the compact moving cubic
span the required horizontal class.

This is a lattice/section readiness certificate.  It deliberately avoids a
global normalization or Jacobian conversion of the bidegree-(14,3) cubic.
"""

import sys
from itertools import product
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, matrix, vector
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


HERE = Path(__file__).resolve().parent
load(str(HERE / "analyze_q80_fourth_q12_cm24_marking.sage"))

# Preserve the third-child marking data before loading the equation search.
old_simple = simple
old_root_gram = root_gram
old_curves = tuple(curves)
old_mw_basis_lifts = mw_basis_lifts
old_projected_basis = projected_basis
old_special_ns = special_ns
old_fiber = fiber
old_zero = zero
fourth_reduced_fiber = reduced
fourth_reflection_sequence = reflection_sequence


def enhanced_basis(transport, embedding, source_frame):
    """Return the full specialized neighbor basis used by enhance_neighbor."""
    lifted = [
        vector(ZZ, list(row[:2]) + list(vector(ZZ, row[2:]) * embedding))
        for row in transport.rows()
    ]
    source_ns = block_diagonal_matrix(U, -source_frame)
    neighbor_fiber, neighbor_mate = lifted[:2]
    complement = matrix(
        ZZ,
        [list(neighbor_fiber * source_ns), list(neighbor_mate * source_ns)],
    ).right_kernel_matrix()
    basis = matrix(
        ZZ,
        [list(neighbor_fiber), list(neighbor_mate)]
        + [list(row) for row in complement.rows()],
    )
    assert abs(basis.det()) == 1
    assert -(complement * source_ns * complement.transpose()) == special_fourth
    return basis


raw_special_fourth_basis = enhanced_basis(
    fourth_transport, third_embedding, special_third
)

# The neighbor constructor uses the raw isotropic fiber.  The compact equation
# uses its chamber-reduced representative.  Apply the identical Weyl word to
# every basis row so that the abstract fourth-child frame is aligned with the
# geometric pencil and its explicit degree-one sections.
old_curve_by_name = dict(old_curves)


def apply_fourth_weyl_word(row):
    row = vector(ZZ, row)
    for name, _ in fourth_reflection_sequence:
        curve_row = old_curve_by_name[name]
        row += intersection(row, curve_row, old_special_ns)*curve_row
    return row


special_fourth_basis = matrix(
    ZZ, [list(apply_fourth_weyl_word(row)) for row in raw_special_fourth_basis.rows()]
)
assert special_fourth_basis[0] == fourth_reduced_fiber
assert (
    special_fourth_basis
    * old_special_ns
    * special_fourth_basis.transpose()
    == block_diagonal_matrix(U, -special_fourth)
)
special_fourth_basis_inverse = special_fourth_basis.inverse().change_ring(ZZ)
geometric_fourth_zero_old = (
    vector(ZZ, [-1, 1] + [0] * 18) * special_fourth_basis
)

fifth_step = steps[4]
generic_fifth = vector(
    ZZ,
    [ZZ(fifth_step["a"]), ZZ(fifth_step["b"])]
    + list(map(ZZ, fifth_step["v"].split(","))),
)
special_fifth = vector(
    ZZ,
    list(generic_fifth[:2])
    + list(vector(ZZ, generic_fifth[2:]) * fourth_embedding),
)

new_ns = block_diagonal_matrix(U, -special_fourth)
new_simple, new_positive = deterministic_simple_roots(special_fourth)
assert new_simple.nrows() == 16
new_root_gram = new_simple * special_fourth * new_simple.transpose()
new_fiber = vector(ZZ, [1, 0] + [0] * 18)
new_zero = vector(ZZ, [-1, 1] + [0] * 18)
new_curves = [("O", new_zero)] + [
    (f"R{index}", vector(ZZ, [0, 0] + list(root)))
    for index, root in enumerate(new_simple.rows(), 1)
]
for component_index, (_, root, _) in enumerate(
    highest_roots(special_fourth, new_simple, new_positive), 1
):
    new_curves.append(
        (
            f"Theta0_{component_index}",
            new_fiber-vector(ZZ, [0, 0] + list(root)),
        )
    )

new_reduced, new_reflections = chamber_reduce(
    special_fifth, new_curves, new_ns
)
assert new_reduced * new_ns * new_reduced == 0
assert all(
    intersection(new_reduced, curve, new_ns) >= 0
    for _, curve in new_curves
)

new_frame_part = vector(QQ, new_reduced[2:])
new_root_coordinates = (
    new_frame_part
    * special_fourth
    * new_simple.transpose()
    * new_root_gram.inverse()
)
new_mw_projection = new_frame_part-new_root_coordinates*new_simple
new_mw_norm = new_mw_projection*special_fourth*new_mw_projection
fifth_reduced_in_old_basis = vector(
    ZZ, new_reduced*special_fourth_basis
)

new_root_data = root_invariants(special_fourth)
assert new_root_data[0] == 16
_, new_mw_height, new_mw_lifts = mw_height_gram(
    special_fourth, new_root_data[3], return_lifts=True
)
new_optimal = optimal_section_pole_basis(
    special_fourth, new_mw_height, new_mw_lifts
)
new_mw_basis_lifts = matrix(ZZ, new_optimal[2]) * new_mw_lifts


def new_project_mw(row):
    row = vector(QQ, row)
    return (
        row
        - row
        * special_fourth
        * new_simple.transpose()
        * new_root_gram.inverse()
        * new_simple
    )


new_projected_basis = matrix(
    QQ, [list(new_project_mw(row)) for row in new_mw_basis_lifts.rows()]
)
new_optimal_height = (
    new_projected_basis
    * special_fourth
    * new_projected_basis.transpose()
)
new_mw_coordinates = new_projected_basis.solve_left(new_mw_projection)
assert new_mw_coordinates*new_projected_basis == new_mw_projection

# Recover the horizontal section's exact component corrections and pole count
# directly from the discriminant class.  This avoids a large rank-16 CVP
# enumeration; the readiness gate only needs the section profile.
new_raw_section_lift = vector(ZZ, new_mw_coordinates*new_mw_basis_lifts)
needed_corrections, needed_pole = section_data_for_lift(
    special_fourth, new_mw_norm, new_raw_section_lift
)

# Load the exact old polynomial-section group and identify the five x-classes
# whose restriction to the compact fourth pencil is nonconstant of degree one.
saved_arguments = sys.argv
search_script = HERE / "search_q80_third_child_polynomial_sections_gf73.sage"
sys.argv = [str(search_script), "--match-mw"]
load(str(search_script))
sys.argv = saved_arguments

basis_points = (points[first], points[second], points[third])
coordinate_lookup = {}
for coefficients in product(range(-6, 7), repeat=3):
    point = sum(
        (coefficient*basis_point for coefficient, basis_point in zip(coefficients, basis_points)),
        curve(0),
    )
    coordinate_lookup.setdefault(point, coefficients)

degree_one_x_keys = (
    (41, 12, 67, 72, 43),
    (24, 63, 53, 61, 41),
    (7, 47, 45, 69, 62),
    (12, 26, 13, 39, 49),
    (48, 60, 25, 12, 30),
)
candidate_rows = []
old_root_lattice = IntegralLattice(old_root_gram)
for x_key in degree_one_x_keys:
    indices = tuple(
        index for index, (section_x, _) in enumerate(two_node_candidates)
        if tuple(map(int, section_x.list())) == x_key
    )
    assert len(indices) == 2
    transported_choices = []
    sign_diagnostics = []
    for index in indices:
        point = points[index]
        assert point in coordinate_lookup
        old_coordinates = vector(ZZ, coordinate_lookup[point])
        old_pole = section_pole(point)
        old_raw_lift = vector(ZZ, old_coordinates*old_mw_basis_lifts)
        old_raw_root_coordinates = (
            vector(QQ, old_raw_lift)
            * special_third
            * old_simple.transpose()
            * old_root_gram.inverse()
        )
        iterator = old_root_lattice.enumerate_close_vectors(
            -old_raw_root_coordinates
        )
        old_sections = []
        target_norm = 2*(old_pole+2)
        for _ in range(8192):
            shift = vector(ZZ, next(iterator))
            lift = old_raw_lift+shift*old_simple
            norm = ZZ(lift*special_third*lift)
            if norm > target_norm:
                break
            if norm != target_norm:
                continue
            section_class = vector(ZZ, [old_pole+1, 1] + list(lift))
            if any(
                intersection(section_class, curve_row, old_special_ns) < 0
                for _, curve_row in old_curves[1:]
            ):
                continue
            old_sections.append(section_class)
        assert len(old_sections) == 1
        transported = old_sections[0]*special_fourth_basis_inverse
        assert all(value in ZZ for value in transported)
        transported = vector(ZZ, transported)
        new_degree = intersection(transported, new_fiber, new_ns)
        sign_diagnostics.append(
            (index, tuple(old_coordinates), old_pole, new_degree)
        )
        if new_degree == 1:
            transported_choices.append(
                (old_coordinates, transported)
            )
    print(
        "Q80FIFTHQ4CM24|old_section_diagnostic|"
        f"X={x_key}|signs={tuple(sign_diagnostics)}|"
        f"degree_one_choices={len(transported_choices)}",
        flush=True,
    )
    assert len(transported_choices) == 1
    old_coordinates, transported = transported_choices[0]
    transported_projection = new_project_mw(transported[2:])
    transported_coordinates = new_projected_basis.solve_left(
        transported_projection
    )
    assert transported_coordinates*new_projected_basis == transported_projection
    assert all(value in ZZ for value in transported_coordinates)
    candidate_rows.append(
        (
            x_key,
            tuple(old_coordinates),
            tuple(map(ZZ, transported)),
            intersection(transported, new_zero, new_ns),
            tuple(map(ZZ, transported_coordinates)),
        )
    )

coordinate_matrix = matrix(
    ZZ, [list(row[-1]) for row in candidate_rows]
)
section_span = coordinate_matrix.row_module()
assert vector(ZZ, new_mw_coordinates) in section_span
span_basis = section_span.basis_matrix()
span_coefficients = span_basis.solve_left(vector(ZZ, new_mw_coordinates))
assert all(value in ZZ for value in span_coefficients)

print(
    "Q80FIFTHQ4CM24|source=A5+A3+4A2/MW2|q=4|"
    f"raw={tuple(special_fifth)}|"
    f"raw_D.F={intersection(special_fifth, new_fiber, new_ns)}|"
    f"raw_D.O={intersection(special_fifth, new_zero, new_ns)}|"
    f"reflection_count={len(new_reflections)}|reduction={new_reflections}|"
    f"reduced={tuple(new_reduced)}|"
    f"D.F={intersection(new_reduced, new_fiber, new_ns)}|"
    f"D.O={intersection(new_reduced, new_zero, new_ns)}",
    flush=True,
)
print(
    "Q80FIFTHQ4CM24|"
    f"mw_norm={new_mw_norm}|mw_coordinates={tuple(new_mw_coordinates)}|"
    f"mw_height={tuple(tuple(row) for row in new_optimal_height.rows())}|"
    f"horizontal_component_corrections={needed_corrections}|"
    f"section_P.O={needed_pole}|"
    f"degree_one_section_rows={tuple(candidate_rows)}|"
    f"section_span_basis={tuple(tuple(row) for row in span_basis.rows())}|"
    f"target_span_coordinates={tuple(span_coefficients)}|"
    "existing_degree_one_sections_span_target=1|"
    "status=PASS_FIFTH_Q4_READINESS",
    flush=True,
)
print(
    "Q80FIFTHQ4CM24|"
    f"geometric_fourth_zero_in_old_basis={tuple(geometric_fourth_zero_old)}|"
    f"old_fiber_degree={intersection(geometric_fourth_zero_old, old_fiber, old_special_ns)}|"
    f"old_zero_pairing={intersection(geometric_fourth_zero_old, old_zero, old_special_ns)}|"
    "status=PASS_GEOMETRIC_ZERO_TRANSPORT",
    flush=True,
)
print(
    "Q80FIFTHQ4CM24|"
    f"fifth_reduced_in_old_basis={tuple(fifth_reduced_in_old_basis)}|"
    f"old_fiber_degree={intersection(fifth_reduced_in_old_basis, old_fiber, old_special_ns)}|"
    f"old_zero_pairing={intersection(fifth_reduced_in_old_basis, old_zero, old_special_ns)}|"
    "status=PASS_FIFTH_DIVISOR_PULLBACK",
    flush=True,
)
