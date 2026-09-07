#!/usr/bin/env sage
"""Transport the alternate final q6 divisor to the CM24 fifth child.

This lattice checker transports the new generic q6 witness through the CM24
specialization of its selected fifth-q4 class, reduces it in the component
chamber, and records its horizontal MW coordinates and section profile.  The
selected child has CM root data distinct from the separate pair14 marked
Jacobian, so matching the vertical class must precede any equation transfer.
"""

from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, matrix, vector


HERE = Path(__file__).resolve().parent
load(str(HERE / "analyze_q80_fifth_q4_cm24_readiness.sage"))

alternate_q4_v = vector(ZZ, (
    -9, 8, -11, 10, -4, 0, 5, 1, -6,
    6, 1, -2, -1, -1, 1, 2, 0,
))
alternate_generic, alternate_transport = neighbor(
    generic_fourth_frame, ZZ(4), ZZ(2), ZZ(2), alternate_q4_v
)
alternate_special, alternate_embedding = enhance_neighbor(
    alternate_transport, fourth_embedding, special_fourth
)
assert (
    alternate_generic
    == alternate_embedding
    * alternate_special
    * alternate_embedding.transpose()
)

special_ns = block_diagonal_matrix(U, -alternate_special)
simple, positive = deterministic_simple_roots(alternate_special)
root_gram = simple*alternate_special*simple.transpose()
root_data = root_invariants(alternate_special)

fiber = vector(ZZ, [1, 0]+[0]*18)
zero = vector(ZZ, [-1, 1]+[0]*18)
curves = [("O", zero)] + [
    (f"R{index}", vector(ZZ, [0, 0]+list(root)))
    for index, root in enumerate(simple.rows(), 1)
]
for component_index, (_, root, _) in enumerate(
    highest_roots(alternate_special, simple, positive), 1
):
    curves.append((
        f"Theta0_{component_index}",
        fiber-vector(ZZ, [0, 0]+list(root)),
    ))

q6_v = vector(ZZ, (
    0, -2, 4, 2, -1, 2, 1, -1, 1,
    0, 1, -1, 1, 0, 0, 0, 0,
))
special_q6 = vector(ZZ, [2, 3]+list(q6_v*alternate_embedding))
assert special_q6*special_ns*special_q6 == 0
reduced, reflections = chamber_reduce(special_q6, curves, special_ns)
assert reduced*special_ns*reduced == 0
assert all(intersection(reduced, curve, special_ns) >= 0 for _, curve in curves)

frame_part = vector(QQ, reduced[2:])
root_coordinates = (
    frame_part
    * alternate_special
    * simple.transpose()
    * root_gram.inverse()
)
mw_projection = frame_part-root_coordinates*simple
mw_norm = mw_projection*alternate_special*mw_projection

_, mw_height, mw_lifts = mw_height_gram(
    alternate_special, root_data[3], return_lifts=True
)
optimal = optimal_section_pole_basis(
    alternate_special, mw_height, mw_lifts
)
mw_basis_lifts = matrix(ZZ, optimal[2])*mw_lifts


def project_mw(row):
    row = vector(QQ, row)
    return (
        row
        - row
        * alternate_special
        * simple.transpose()
        * root_gram.inverse()
        * simple
    )


projected_basis = matrix(
    QQ, [list(project_mw(row)) for row in mw_basis_lifts.rows()]
)
optimal_height = (
    projected_basis
    * alternate_special
    * projected_basis.transpose()
)
coordinates = projected_basis.solve_left(mw_projection)
assert coordinates*projected_basis == mw_projection
raw_section_lift = vector(ZZ, coordinates*mw_basis_lifts)
corrections, pole = section_data_for_lift(
    alternate_special, mw_norm, raw_section_lift
)

print(
    "Q80ALTFINALQ6CM24|"
    f"special_roots={root_data[:3]}|special_MW={18-root_data[0]}|"
    f"raw={tuple(special_q6)}|raw_D.F={intersection(special_q6, fiber, special_ns)}|"
    f"raw_D.O={intersection(special_q6, zero, special_ns)}|"
    f"reflections={tuple(reflections)}|reduced={tuple(reduced)}|"
    f"D.F={intersection(reduced, fiber, special_ns)}|"
    f"D.O={intersection(reduced, zero, special_ns)}|"
    f"mw_norm={mw_norm}|mw_coordinates={tuple(coordinates)}|"
    f"mw_height={tuple(tuple(row) for row in optimal_height.rows())}|"
    f"horizontal_component_corrections={corrections}|section_P.O={pole}|"
    "status=PASS_CM24_READINESS",
    flush=True,
)
