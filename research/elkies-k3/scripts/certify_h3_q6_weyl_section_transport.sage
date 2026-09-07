#!/usr/bin/env sage -python
"""Transport the canonical H3 q=6 MW divisor basis through its Weyl record.

The canonical MW3 frame in ``analyze_h3_first_q6_chamber.sage`` is expressed
in a complement for the raw isotropic class.  The equation-level q=6 pencil,
however, uses its nef representative ``O+(-P1)-F``.  This script applies the
recorded reflections to *each section class*, constructs the new root
projection relative to the old zero section, and only then reads the old MW
words.  It prevents the tempting but invalid identification of coordinate
vectors from the two complements.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import IntegralLattice, QQ, ZZ, block_diagonal_matrix, gcd, identity_matrix, matrix, pari, vector, xgcd


ROOT = Path(__file__).resolve().parents[2]
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h3-q6-weyl-section-transport.json"
REFLECTIONS = (1, 2, 4, 3, 5, 4, 2, 6, 5, 4, 3, 1, 7, 6, 5, 4, 2, 3, 4, 5, 6, 7)
H3_LIFTS = matrix(ZZ, [
    [-5, -4, -3, 0, 0, 0, 0, 0, 0, 0, 0, -4, 1, 0, -4, 2, -2],
    [-10, -8, -6, 0, 0, 0, 0, 0, 0, 0, 0, -8, 4, 1, -8, 5, -4],
    [-5, -4, -3, 0, 0, 0, 0, 0, 0, 0, 0, -3, 2, 0, -4, 2, -2],
])
OLD_ZERO_ROOT_SHIFTS = matrix(ZZ, [
    [5, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3, -1, 4],
    [12, 10, 8, 0, 0, 0, 0, 0, 0, 0, 0, 6, -1, 9],
    [5, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 4],
])
EXPECTED_HEIGHT = matrix(QQ, [
    [QQ(8) / 3, QQ(1) / 3, -1],
    [QQ(1) / 3, QQ(8) / 3, 1],
    [-1, 1, 46],
])

exec(compile(CORE.read_text(), str(CORE), "exec"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


def isotropic_mate(ns, fiber):
    current = ZZ(0)
    data = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(ns * fiber):
        if not value:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        data = [left * entry for entry in data]
        data[index] += right
        current = divisor
    assert abs(current) == 1
    if current == -1:
        data = [-entry for entry in data]
    mate = vector(ZZ, data)
    mate -= (mate * ns * mate // 2) * fiber
    assert mate * ns * mate == 0 and mate * ns * fiber == 1
    return mate


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--old-zero-cvp", action="store_true", help="diagnose root-CVP section lifts relative to old O")
parser.add_argument("--cvp-cap", type=int, default=8192)
args = parser.parse_args()

frame = load_gram(FRAME)
ns = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -frame)
fiber = vector(ZZ, [1, 0] + [0] * 17)
zero = vector(ZZ, [-1, 1] + [0] * 17)
simple = tuple(vector(ZZ, [0, 0] + [ZZ(index == node) for index in range(17)]) for node in range(15))

raw_fiber = vector(ZZ, [3, 2] + [0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0])
reflection_roots = tuple(("E7_{}".format(node), simple[node - 1]) for node in REFLECTIONS)
nef_fiber, reflection_record = replay_weyl_reflections(
    raw_fiber, ns, reflection_roots, expected_pairings=(-1,) * len(REFLECTIONS)
)
assert nef_fiber == zero + vector(ZZ, [5, 1] + [-value for value in (2, 3, 4, 6, 5, 4, 3)] + [0] * 8 + [1, 0]) - fiber
assert nef_fiber * ns * nef_fiber == 0 and gcd(tuple(ns * nef_fiber)) == 1


def weyl_transport(class_value):
    return replay_weyl_reflections(class_value, ns, reflection_roots)[0]


def inverse_weyl_transport(class_value):
    return replay_weyl_reflections(class_value, ns, tuple(reversed(reflection_roots)))[0]


assert weyl_transport(raw_fiber) == nef_fiber
assert inverse_weyl_transport(nef_fiber) == raw_fiber
raw_mate = isotropic_mate(ns, raw_fiber)
raw_orthogonal = matrix(ZZ, [list(raw_fiber * ns), list(raw_mate * ns)]).right_kernel_matrix()
raw_transport = matrix(ZZ, [list(raw_fiber), list(raw_mate)] + [list(row) for row in raw_orthogonal.rows()])
assert abs(raw_transport.det()) == 1
raw_child = -(raw_orthogonal * ns * raw_orthogonal.transpose())
raw_roots = matrix(ZZ, pari(raw_child).qfminim(2)[2]).transpose().row_module().basis_matrix()
raw_root_gram = raw_roots * raw_child * raw_roots.transpose()
assert raw_roots.rank() == 14 and abs(raw_root_gram.det()) == 3

# These are the canonical saturated MW directions in the raw complement.
raw_zero = inverse_weyl_transport(zero)
assert raw_zero * ns * raw_fiber == 1 and raw_zero * ns * raw_zero == -2
raw_zero_coordinates = vector(ZZ, raw_zero * raw_transport.inverse())
assert raw_zero_coordinates[1] == 1
raw_zero_mw_lift = vector(ZZ, raw_zero_coordinates[2:])
print("H3Q6WEYLTRANSPORT|raw_zero_mw_lift={}".format(tuple(raw_zero_mw_lift)), flush=True)

if args.old_zero_cvp:
    root_lattice = IntegralLattice(raw_root_gram)
    old_zero_cvp_sections = []
    for index, h3_lift in enumerate(H3_LIFTS.rows(), 1):
        target_lift = raw_zero_mw_lift + vector(ZZ, h3_lift)
        root_coordinates = target_lift * raw_child * raw_roots.transpose() * raw_root_gram.inverse()
        iterator = root_lattice.enumerate_close_vectors(-root_coordinates)
        selected = []
        first_norm = None
        for _ in range(args.cvp_cap):
            shift = vector(ZZ, next(iterator))
            candidate_lift = target_lift + shift * raw_roots
            norm = ZZ(candidate_lift * raw_child * candidate_lift)
            if first_norm is not None and norm > first_norm:
                break
            if norm < 4 or norm % 2:
                continue
            pole = norm // 2 - 2
            candidate_raw = (pole + 1) * raw_fiber + raw_mate + candidate_lift * raw_orthogonal
            candidate = weyl_transport(candidate_raw)
            pairings = tuple(int(candidate * ns * root) for root in simple)
            if candidate * ns * candidate != -2 or candidate * ns * nef_fiber != 1 or min(pairings) < 0:
                continue
            if first_norm is None:
                first_norm = norm
            selected.append((candidate, pole, pairings, shift))
        print(
            "H3Q6WEYLTRANSPORT|old_zero_cvp={}|candidates={}".format(index, len(selected)),
            flush=True,
        )
        for candidate, pole, pairings, shift in selected:
            print(
                "H3Q6WEYLTRANSPORT|old_zero_lift={}|pole={}|old_degree={}|shift={}|pairings={}".format(
                    index, pole, candidate * ns * fiber, tuple(shift), pairings
                ),
                flush=True,
            )
        if len(selected) == 1:
            old_zero_cvp_sections.append(selected[0][0])
else:
    old_zero_cvp_sections = []
raw_sections = []
for lift in H3_LIFTS.rows():
    perpendicular = vector(ZZ, lift) * raw_orthogonal
    start = raw_zero + perpendicular
    section = start + ((-2 - start * ns * start) // 2) * raw_fiber
    assert section * ns * section == -2 and section * ns * raw_fiber == 1
    raw_sections.append(section)

# The raw child frame has a different zero.  Translate the published MW
# lifts by the raw-frame coordinates of old O, then apply these closest-root
# corrections.  The shifts were found by the optional root-CVP diagnostic
# below and are checked here as exact data, not trusted as a search result.
sections = []
effective_lift_data = []
for h3_lift, shift in zip(H3_LIFTS.rows(), OLD_ZERO_ROOT_SHIFTS.rows()):
    candidate_lift = raw_zero_mw_lift + vector(ZZ, h3_lift) + vector(ZZ, shift) * raw_roots
    norm = ZZ(candidate_lift * raw_child * candidate_lift)
    assert norm >= 4 and norm % 2 == 0
    pole = norm // 2 - 2
    candidate_raw = (pole + 1) * raw_fiber + raw_mate + candidate_lift * raw_orthogonal
    candidate = weyl_transport(candidate_raw)
    pairings = tuple(int(candidate * ns * root) for root in simple)
    assert candidate * ns * candidate == -2 and candidate * ns * nef_fiber == 1
    assert min(pairings) >= 0
    sections.append(candidate)
    effective_lift_data.append({
        "root_shift": [int(value) for value in shift],
        "pole_against_old_zero": int(pole),
        "old_fiber_degree": int(candidate * ns * fiber),
        "old_simple_component_pairings": list(pairings),
    })
assert [entry["pole_against_old_zero"] for entry in effective_lift_data] == [5, 3, 22]
assert [entry["old_fiber_degree"] for entry in effective_lift_data] == [10, 3, 44]

old_fiber_degrees = [int(section * ns * fiber) for section in sections]
assert old_fiber_degrees == [10, 3, 44]

# In the nef model, compute the vertical root lattice and the Shioda map with
# the actual child zero old O.  This is the same model used by the q=6 pencil.
nef_mate = isotropic_mate(ns, nef_fiber)
nef_orthogonal = matrix(ZZ, [list(nef_fiber * ns), list(nef_mate * ns)]).right_kernel_matrix()
nef_child = -(nef_orthogonal * ns * nef_orthogonal.transpose())
nef_roots = matrix(ZZ, pari(nef_child).qfminim(2)[2]).transpose().row_module().basis_matrix() * nef_orthogonal
root_gram = nef_roots * ns * nef_roots.transpose()
assert nef_roots.rank() == 14 and abs(root_gram.det()) == 3
projection = identity_matrix(QQ, 19) - ns * nef_roots.transpose() * root_gram.inverse() * nef_roots


def shioda(section):
    horizontal = section - zero - (section * ns * zero + 2) * nef_fiber
    assert horizontal * ns * nef_fiber == horizontal * ns * zero == 0
    return vector(QQ, horizontal) * projection


shioda_sections = [shioda(section) for section in sections]
height = matrix(QQ, [[-left * ns * right for right in shioda_sections] for left in shioda_sections])
assert height == EXPECTED_HEIGHT
if args.old_zero_cvp:
    if len(old_zero_cvp_sections) != 3:
        raise RuntimeError("old-zero CVP did not select all three section representatives")
    old_zero_cvp_height = matrix(QQ, [
        [-left * ns * right for right in [shioda(section) for section in old_zero_cvp_sections]]
        for left in [shioda(section) for section in old_zero_cvp_sections]
    ])
    print("H3Q6WEYLTRANSPORT|old_zero_cvp_height={}".format(old_zero_cvp_height), flush=True)
    assert old_zero_cvp_height == EXPECTED_HEIGHT

# Read the exact old MW *projection* after the Weyl transport. Root
# corrections and old-fibre degree are discarded by the source Shioda map, so
# these coordinates must not be evaluated as old-model section coordinates.
minus_p1 = vector(ZZ, [5, 1] + [-value for value in (2, 3, 4, 6, 5, 4, 3)] + [0] * 8 + [1, 0])
p2 = vector(ZZ, [22, 1] + [0] * 16 + [1])
source_roots = matrix(ZZ, [list(root) for root in simple])


def source_shioda(section):
    horizontal = section - zero - (section * ns * zero + 2) * fiber
    root_gram = source_roots * ns * source_roots.transpose()
    source_projection = identity_matrix(QQ, 19) - ns * source_roots.transpose() * root_gram.inverse() * source_roots
    return vector(QQ, horizontal) * source_projection


source_basis = matrix(QQ, [source_shioda(minus_p1), source_shioda(p2)])
source_pairing = source_basis * ns * source_basis.transpose()
words = []
for section in sections:
    value = source_shioda(section)
    coordinates = vector(QQ, [value * ns * source_basis[index] for index in range(2)]) * source_pairing.inverse()
    assert all(coordinate in ZZ for coordinate in coordinates)
    words.append([int(coordinate) for coordinate in coordinates])
assert words == [[4, 0], [1, 0], [22, -1]]

# For an old-fibre degree-d multisection S, the old MW point governing
# ``S|_E ~(d-1)O+P`` is obtained from S-dO, not from the degree-one Shioda
# expression above.  Its source-MW coordinates are the data needed to build
# the high-degree old Riemann--Roch ambient for an equation-level transport.
abel_jacobi_words = []
for section in sections:
    degree = section * ns * fiber
    degree_zero = section - degree * zero
    degree_zero -= (degree_zero * ns * zero) * fiber
    abel_jacobi = vector(QQ, degree_zero) * (
        identity_matrix(QQ, 19)
        - ns * source_roots.transpose()
        * (source_roots * ns * source_roots.transpose()).inverse()
        * source_roots
    )
    coordinates = vector(QQ, [abel_jacobi * ns * source_basis[index] for index in range(2)]) * source_pairing.inverse()
    assert all(coordinate in ZZ for coordinate in coordinates)
    abel_jacobi_words.append([int(coordinate) for coordinate in coordinates])
assert abel_jacobi_words == [[4, 0], [1, 0], [22, -1]]

# Spell out the vertical correction for the third class. Its old generic
# fibre restriction is ``43*O + (22*(-P1)-P2)``; the residual is genuinely
# old-vertical. These coefficients are the input for resolved-chart
# conditions in the eventual degree-44 divisor compiler.
third_shioda = 22 * source_basis[0] - source_basis[1]
assert all(entry in ZZ for entry in third_shioda)
third_generic_point = zero + third_shioda - (third_shioda * ns * third_shioda // 2) * fiber
assert third_generic_point * ns * third_generic_point == -2
assert third_generic_point * ns * fiber == 1
third_horizontal = (old_fiber_degrees[2] - 1) * zero + third_generic_point
vertical_correction = sections[2] - third_horizontal
assert vertical_correction * ns * fiber == 0
vertical_basis = matrix(ZZ, [list(root) for root in simple] + [list(fiber)])
vertical_coordinates = vector(QQ, vertical_basis.transpose().solve_right(vertical_correction))
assert all(entry in ZZ for entry in vertical_coordinates)
assert vertical_correction == vector(ZZ, vertical_coordinates) * vertical_basis
assert nef_fiber * ns * third_generic_point == 4769
assert nef_fiber * ns * third_horizontal == 4812
assert nef_fiber * ns * vertical_correction == -4811
assert nef_fiber * ns * sections[2] == 1
third_component_intersections = [int(sections[2] * ns * root) for root in simple[:7]]
third_horizontal_component_intersections = [int(third_horizontal * ns * root) for root in simple[:7]]
third_vertical_component_intersections = [int(vertical_correction * ns * root) for root in simple[:7]]
assert all(
    section_value == horizontal_value + vertical_value
    for section_value, horizontal_value, vertical_value in zip(
        third_component_intersections,
        third_horizontal_component_intersections,
        third_vertical_component_intersections,
    )
)



payload = {
    "schema": "elkies-k3.h3-q6-weyl-section-transport.v1",
    "status": "PASS_EXACT_Q6_WEYL_SECTION_TRANSPORT",
    "inputs": {"frame": {"path": str(FRAME.relative_to(ROOT)), "sha256": digest(FRAME)}},
    "raw_to_nef": {"reflections": list(reflection_record), "raw_old_fiber_degree": int(raw_fiber * ns * fiber), "nef_old_fiber_degree": int(nef_fiber * ns * fiber)},
    "old_zero_in_raw_child_frame": {"coordinates": [int(entry) for entry in raw_zero_coordinates], "mw_lift": [int(entry) for entry in raw_zero_mw_lift]},
    "child": {"root_lattice": "E8+E6", "root_rank": 14, "root_determinant": 3, "height_gram": [[str(value) for value in row] for row in height.rows()]},
    "effective_lifts": effective_lift_data,
    "transported_section_mw_projections": {"basis": ["-P1", "P2"], "coordinates": words, "formulas": ["4*(-P1)", "(-P1)", "22*(-P1)-P2"], "old_fiber_degrees": old_fiber_degrees},
    "old_abel_jacobi_points": {"basis": ["-P1", "P2"], "coordinates": abel_jacobi_words},
    "third_vertical_correction": {
        "old_generic_restriction": "43*O + (22*(-P1)-P2)",
        "old_generic_point_new_fiber_degree": int(nef_fiber * ns * third_generic_point),
        "old_horizontal_new_fiber_degree": int(nef_fiber * ns * third_horizontal),
        "transported_section_new_fiber_degree": int(nef_fiber * ns * sections[2]),
        "correction_new_fiber_degree": int(nef_fiber * ns * vertical_correction),
        "basis": ["old_E7_1", "old_E7_2", "old_E7_3", "old_E7_4", "old_E7_5", "old_E7_6", "old_E7_7", "old_E8_1", "old_E8_2", "old_E8_3", "old_E8_4", "old_E8_5", "old_E8_6", "old_E8_7", "old_E8_8", "old_F"],
        "coordinates": [int(entry) for entry in vertical_coordinates],
        "old_E7_component_intersections": {
            "transported_section": third_component_intersections,
            "horizontal_part": third_horizontal_component_intersections,
            "vertical_correction": third_vertical_component_intersections,
        },
    },
    "boundary": "These are exact marked NS divisor classes and their source MW-projection coordinates. The latter are not old-model section coordinates: in particular, the third marked divisor has old-fibre degree 44 and nontrivial old E7 fixed-component data. Converting it to minimized child coordinates requires resolved divisor transport, not naive scalar multiplication on the old Weierstrass equation.",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("H3Q6WEYLTRANSPORT|words=4(-P1),(-P1),22(-P1)-P2|gram=PASS|status=PASS_EXACT_Q6_WEYL_SECTION_TRANSPORT", flush=True)
