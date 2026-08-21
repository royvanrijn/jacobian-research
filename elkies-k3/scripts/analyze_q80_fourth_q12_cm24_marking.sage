#!/usr/bin/env sage
"""Transport the fourth q12 trisection to the exact CM24 third child.

The generic fourth class is proved nef in
``analyze_q80_fourth_q12_chamber.sage``.  This script adjoins the CM24 class,
transports the same divisor through the first three exact neighbors, reduces
it in the specialized ``2A6+3A1`` chamber, and expresses its horizontal part
in a saturated rank-three CM Mordell--Weil basis.  It also records the
shortest effective section representative and the specialized fourth-child
root system.

This is a lattice/marking certificate.  It does not yet construct the fourth
equation-level pencil over the characteristic-zero CM coefficient field.
"""

from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, matrix, vector
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


HERE = Path(__file__).resolve().parent
load(str(HERE / "analyze_q80_third_q12_cm24_marking.sage"))


fourth_step = steps[3]
generic_fourth = vector(
    ZZ,
    [ZZ(fourth_step["a"]), ZZ(fourth_step["b"])]
    + list(map(ZZ, fourth_step["v"].split(","))),
)
generic_fourth_frame, fourth_transport = neighbor(
    generic_third_frame,
    ZZ(fourth_step["q"]),
    ZZ(fourth_step["a"]),
    ZZ(fourth_step["b"]),
    vector(ZZ, map(ZZ, fourth_step["v"].split(","))),
)
special_fourth, fourth_embedding = enhance_neighbor(
    fourth_transport, third_embedding, special_third
)
assert (
    generic_fourth_frame
    == fourth_embedding * special_fourth * fourth_embedding.transpose()
)

special_divisor = vector(
    ZZ,
    list(generic_fourth[:2])
    + list(vector(ZZ, generic_fourth[2:]) * third_embedding),
)
special_ns = block_diagonal_matrix(U, -special_third)
simple, positive = deterministic_simple_roots(special_third)
assert simple.nrows() == 15
fiber = vector(ZZ, [1, 0] + [0] * 18)
zero = vector(ZZ, [-1, 1] + [0] * 18)
curves = [("O", zero)] + [
    (f"R{index}", vector(ZZ, [0, 0] + list(root)))
    for index, root in enumerate(simple.rows(), 1)
]
for component_index, (_, root, _) in enumerate(
    highest_roots(special_third, simple, positive), 1
):
    curves.append(
        (
            f"Theta0_{component_index}",
            fiber - vector(ZZ, [0, 0] + list(root)),
        )
    )

reduced, reflection_sequence = chamber_reduce(
    special_divisor, curves, special_ns
)
assert reduced * special_ns * reduced == 0
assert all(
    intersection(reduced, curve, special_ns) >= 0
    for _, curve in curves
)

root_gram = simple * special_third * simple.transpose()
frame_part = vector(QQ, reduced[2:])
root_coordinates = (
    frame_part * special_third * simple.transpose() * root_gram.inverse()
)
mw_projection = frame_part - root_coordinates * simple
mw_norm = mw_projection * special_third * mw_projection

root_data = root_invariants(special_third)
_, mw_height, mw_lifts = mw_height_gram(
    special_third, root_data[3], return_lifts=True
)
optimal = optimal_section_pole_basis(special_third, mw_height, mw_lifts)
mw_basis_lifts = matrix(ZZ, optimal[2]) * mw_lifts


def project_mw(row):
    row = vector(QQ, row)
    return (
        row
        - row
        * special_third
        * simple.transpose()
        * root_gram.inverse()
        * simple
    )


projected_basis = matrix(
    QQ, [list(project_mw(row)) for row in mw_basis_lifts.rows()]
)
optimal_height = projected_basis * special_third * projected_basis.transpose()
mw_coordinates = projected_basis.solve_left(mw_projection)
assert mw_coordinates * projected_basis == mw_projection
optimal_pole_data = section_pole_data(
    special_third, optimal_height, mw_basis_lifts
)
assert tuple(pole for _, pole in optimal_pole_data) == optimal[1]

# Find the effective section in this MW class.  Root-CVP candidates are
# ordered by norm; component nonnegativity selects the unique section class.
raw_section_lift = vector(ZZ, mw_coordinates * mw_basis_lifts)
raw_root_coordinates = (
    vector(QQ, raw_section_lift)
    * special_third
    * simple.transpose()
    * root_gram.inverse()
)
root_lattice = IntegralLattice(root_gram)
effective_sections = []
iterator = root_lattice.enumerate_close_vectors(-raw_root_coordinates)
first_effective_norm = None
for _ in range(512):
    shift = vector(ZZ, next(iterator))
    lift = raw_section_lift + shift * simple
    norm = ZZ(lift * special_third * lift)
    if first_effective_norm is not None and norm > first_effective_norm:
        break
    if norm < 4 or norm % 2:
        continue
    pole = norm // 2 - 2
    section = vector(ZZ, [pole + 1, 1] + list(lift))
    component_pairings = tuple(
        intersection(section, curve, special_ns)
        for _, curve in curves[1:]
    )
    if min(component_pairings) < 0:
        continue
    if first_effective_norm is None:
        first_effective_norm = norm
    effective_sections.append((section, component_pairings))

assert effective_sections
assert len(effective_sections) == 1
section, section_component_pairings = effective_sections[0]
section_pole = intersection(section, zero, special_ns)
divisor_section_pairing = intersection(reduced, section, special_ns)
vertical_remainder = reduced - section - 2*zero
fiber_shift = ZZ(vertical_remainder[0])
vertical_remainder -= fiber_shift*fiber
assert vertical_remainder[:2] == vector(ZZ, (0, 0))
integral_vertical_coordinates = simple.solve_left(vertical_remainder[2:])
assert all(value in ZZ for value in integral_vertical_coordinates)
components = connected_components(root_gram)
vertical_components = tuple(
    (
        component,
        tuple(-ZZ(integral_vertical_coordinates[index]) for index in component),
    )
    for component in components
)


def chain_order(component):
    component = tuple(component)
    if len(component) == 1:
        return component
    neighbors = {
        index: tuple(
            other
            for other in component
            if other != index and root_gram[index, other] != 0
        )
        for index in component
    }
    endpoints = sorted(index for index in component if len(neighbors[index]) == 1)
    assert len(endpoints) == 2
    order = [endpoints[0]]
    previous = None
    while len(order) < len(component):
        choices = [
            index for index in neighbors[order[-1]] if index != previous
        ]
        assert len(choices) == 1
        previous, current = order[-1], choices[0]
        order.append(current)
    return tuple(order)


component_chains = tuple(chain_order(component) for component in components)
vertical_chain_coefficients = tuple(
    tuple(-ZZ(integral_vertical_coordinates[index]) for index in chain)
    for chain in component_chains
)
section_simple_pairings = section_component_pairings[:simple.nrows()]
section_chain_hits = tuple(
    tuple(section_simple_pairings[index] for index in chain)
    for chain in component_chains
)

special_child_components = root_components(special_fourth)
special_child_root_rank = sum(component[0] for component in special_child_components)
special_child_mw_rank = 18 - special_child_root_rank

print(
    "Q80FOURTHQ12CM24|source=2A6+3A1/MW3|q=12|ab=3,4|"
    f"raw={tuple(special_divisor)}|"
    f"raw_D.F={intersection(special_divisor, fiber, special_ns)}|"
    f"raw_D.O={intersection(special_divisor, zero, special_ns)}|"
    f"reflection_count={len(reflection_sequence)}|"
    f"reduction={reflection_sequence}",
    flush=True,
)
print(
    f"Q80FOURTHQ12CM24|reduced={tuple(reduced)}|"
    f"D.F={intersection(reduced, fiber, special_ns)}|"
    f"D.O={intersection(reduced, zero, special_ns)}|"
    f"root_coordinates={tuple(root_coordinates)}|"
    f"mw_norm={mw_norm}|mw_coordinates={tuple(mw_coordinates)}|"
    f"optimal_P.O={optimal[1]}|"
    f"optimal_height={tuple(tuple(row) for row in optimal_height.rows())}|"
    f"optimal_component_corrections={optimal_pole_data}",
    flush=True,
)
print(
    f"Q80FOURTHQ12CM24|section={tuple(section)}|"
    f"section_P.O={section_pole}|D.S={divisor_section_pairing}|"
    f"section_component_pairings={section_component_pairings}|"
    f"decomposition=Q+2O+{fiber_shift}F+root_correction|"
    f"integral_root_coordinates={tuple(integral_vertical_coordinates)}|"
    f"vertical_components={vertical_components}|"
    f"component_chains={component_chains}|"
    f"vertical_chain_coefficients={vertical_chain_coefficients}|"
    f"section_chain_hits={section_chain_hits}|"
    f"child_components={special_child_components}|"
    f"child_root_rank={special_child_root_rank}|child_MW={special_child_mw_rank}|"
    "status=PASS_CM24_MARKING",
    flush=True,
)
