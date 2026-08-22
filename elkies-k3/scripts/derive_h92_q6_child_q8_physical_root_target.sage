#!/usr/bin/env sage -python
"""Derive a physical E8+E6 vertical target for an H92 child q=8 class.

``derive_h92_q6_child_q8_marking.sage`` first expresses the selected q=8
class in a convenient *canonical-zero* root basis.  The explicit q=6 child
uses the transported old zero instead.  This script changes the hyperbolic
plane to ``<F,O_old+F>`` and extracts the finite fibre root lattice relative
to that actual zero.  It then rewrites the exact marked divisor there.

The result is a lattice-component target only: the labels below are a pinned
simple-root order, not yet names of charts in a resolved II*/IV* equation.
In particular it neither derives the local quotient maps nor a q=8 pencil.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
Q8_ORBITS = ROOT / "artifacts/generated-results/elkies-k3-h3-q6-q8-orbits.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-physical-root-target.json"

REFLECTIONS = (1, 2, 4, 3, 5, 4, 2, 6, 5, 4, 3, 1, 7, 6, 5, 4, 2, 3, 4, 5, 6, 7)
E6_QF_INDICES = (0, 1, 2, 3, 12, 13)
# This deterministic Weyl chamber has V.E6_i=(-1,-1,0,0,0,0).  The entries
# are coordinates in the qfminim basis indexed by E6_QF_INDICES.
E6_SIMPLE_IN_QF_COORDINATES = (
    (-1, -1, -1, -1, 0, 0),
    (0, 0, 0, 1, 0, 0),
    (0, 0, 1, 0, 0, 0),
    (0, 1, 0, 0, 0, 0),
    (1, 0, 0, 0, 0, 1),
    (2, 1, 0, 0, -1, 1),
)
E6_CARTAN = matrix(ZZ, [
    [2, -1, 0, 0, 0, 0],
    [-1, 2, -1, 0, 0, 0],
    [0, -1, 2, -1, 0, -1],
    [0, 0, -1, 2, -1, 0],
    [0, 0, 0, -1, 2, 0],
    [0, 0, -1, 0, 0, 2],
])
E8_CARTAN = matrix(ZZ, [
    [2, 0, -1, 0, 0, 0, 0, 0],
    [0, 2, 0, -1, 0, 0, 0, 0],
    [-1, 0, 2, -1, 0, 0, 0, 0],
    [0, -1, -1, 2, -1, 0, 0, 0],
    [0, 0, 0, -1, 2, -1, 0, 0],
    [0, 0, 0, 0, -1, 2, -1, 0],
    [0, 0, 0, 0, 0, -1, 2, -1],
    [0, 0, 0, 0, 0, 0, -1, 2],
])


def highest_root_coefficients(cartan):
    """Return the coordinate vector of the dominant highest root."""
    half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
    roots = tuple(vector(ZZ, row) for row in half) + tuple(-vector(ZZ, row) for row in half)
    return max(
        (root for root in roots if all(value >= 0 for value in root)),
        key=lambda root: sum(root),
    )


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--marking", type=Path, default=MARKING)
parser.add_argument("--q8-orbits", type=Path, default=Q8_ORBITS)
parser.add_argument(
    "--representative", choices=("nef", "dominant-d13", "component-nef"), default="dominant-d13",
    help="use the dominant D13, abstract Weyl-nef, or q6-component-dominant q8 fibre",
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
for name in ("marking", "q8_orbits", "output"):
    setattr(args, name, getattr(args, name).resolve())

marking = json.loads(args.marking.read_text())
q8_orbits = json.loads(args.q8_orbits.read_text())
assert marking["status"] == "PASS_EXACT_Q6_CHILD_Q8_MARKING"
assert q8_orbits["status"] == "PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION"

ns = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -load_gram(FRAME))
simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(index == node) for index in range(17)])
    for node in range(15)
)


def reflect(class_value, nodes):
    result = vector(ZZ, class_value)
    for node in nodes:
        curve = simple[node - 1]
        result += (result * ns * curve) * curve
    return result


fibre = vector(ZZ, [3, 2] + [
    0, 0, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
])
old_zero = reflect(vector(ZZ, [-1, 1] + [0] * 17), tuple(reversed(REFLECTIONS)))
assert fibre * ns * fibre == 0
assert old_zero * ns * old_zero == -2 and old_zero * ns * fibre == 1

# ``old_zero+fibre`` is the isotropic mate appropriate to the explicit child
# origin.  Pari returns a lattice basis, not necessarily simple roots; the
# fixed E6 chamber below turns its E6 completion into actual roots.
old_zero_orthogonal = matrix(
    ZZ, [list(fibre * ns), list((old_zero + fibre) * ns)]
).right_kernel_matrix()
physical_child = -(old_zero_orthogonal * ns * old_zero_orthogonal.transpose())
qf_basis = matrix(
    ZZ, pari(physical_child).qfminim(2)[2]
).transpose().row_module().basis_matrix()
assert qf_basis.rank() == 14
physical_roots = qf_basis * old_zero_orthogonal
physical_frame = matrix(ZZ, [list(fibre), list(old_zero + fibre)] + [
    list(row) for row in old_zero_orthogonal.rows()
])
assert abs(physical_frame.det()) == 1

e6_qf_basis = matrix(ZZ, [list(qf_basis[index]) for index in E6_QF_INDICES])
e6_roots = matrix(ZZ, [
    vector(ZZ, coordinates) * e6_qf_basis
    for coordinates in E6_SIMPLE_IN_QF_COORDINATES
]) * old_zero_orthogonal
e8_roots = physical_roots[4:12, :]
assert -e6_roots * ns * e6_roots.transpose() == E6_CARTAN
assert -e8_roots * ns * e8_roots.transpose() == E8_CARTAN
assert e6_roots * ns * e8_roots.transpose() == matrix(ZZ, 6, 8)
assert all(root * ns * old_zero == 0 for root in e6_roots.rows())
assert all(root * ns * old_zero == 0 for root in e8_roots.rows())

if args.representative in {"nef", "component-nef"}:
    representative_record = q8_orbits["q8"]["nef_representative"]
    assert representative_record["mw_projection"] == [0, -2, 0]
    selected_q8 = vector(ZZ, representative_record["fiber_source_h3_ns"])
else:
    representative_record = next(
        hit for hit in q8_orbits["q8"]["d13_mw4_hits"]
        if hit["mw_projection"] == [0, -2, 0]
    )
    selected_q8 = vector(ZZ, representative_record["fiber_source_h3_ns"])

# The abstract q8 chamber representative still has negative intersections
# with components that are effective in the explicit q6 Weierstrass model.
# For the actual old-fibre equation, move it to the dominant chamber of this
# E6+E8 configuration before interpreting it as a pencil divisor.
effective_reflections = []
if args.representative == "component-nef":
    actual_roots = tuple(e6_roots.rows()) + tuple(e8_roots.rows())
    for unused in range(500):
        pairings = [int(selected_q8 * ns * root) for root in actual_roots]
        negative = [index for index, value in enumerate(pairings) if value < 0]
        if not negative:
            break
        index = negative[0]
        pairing = pairings[index]
        selected_q8 += pairing * actual_roots[index]
        effective_reflections.append((index, pairing))
    else:
        raise RuntimeError("actual E6+E8 Weyl reduction did not terminate")
    assert all(selected_q8 * ns * root >= 0 for root in actual_roots)
    assert selected_q8 * ns * selected_q8 == 0
    assert selected_q8 * ns * fibre == 2
selected_section = vector(ZZ, marking["selected_q8"][
    "ns_horizontal_vertical_decomposition"
]["section_class_in_source_h3_ns"])
assert selected_section * ns * selected_section == -2
assert selected_section * ns * fibre == 1
assert all(selected_section * ns * root == 0 for root in e6_roots.rows())
assert all(selected_section * ns * root == 0 for root in e8_roots.rows())

vertical = selected_q8 - old_zero - selected_section
physical_basis = matrix(QQ, [list(fibre)] + [
    list(root) for root in e6_roots.rows()
] + [list(root) for root in e8_roots.rows()])
coordinates = physical_basis.solve_left(vertical)
assert all(value in ZZ for value in coordinates)
e6_cycle = vector(ZZ, coordinates[1:7])
e8_cycle = vector(ZZ, coordinates[7:15])
e6_degrees = vector(ZZ, vertical * ns * e6_roots.transpose())
e8_degrees = vector(ZZ, vertical * ns * e8_roots.transpose())
assert e6_cycle * E6_CARTAN == -e6_degrees
assert e8_cycle * E8_CARTAN == -e8_degrees
e6_highest = highest_root_coefficients(E6_CARTAN)
e8_highest = highest_root_coefficients(E8_CARTAN)
e6_affine_degree = selected_q8 * ns * (fibre - e6_highest * e6_roots)
e8_affine_degree = selected_q8 * ns * (fibre - e8_highest * e8_roots)
assert e6_affine_degree == 2 - e6_highest * e6_degrees
assert e8_affine_degree == 2 - e8_highest * e8_degrees
assert selected_q8 == old_zero + selected_section + vertical
selected_q8_physical_coordinates = physical_frame.solve_left(selected_q8)
assert all(value in ZZ for value in selected_q8_physical_coordinates)
assert selected_q8_physical_coordinates[1] == 2
assert (
    selected_q8_physical_coordinates[2:]
    * physical_child
    * selected_q8_physical_coordinates[2:]
    == 4 * selected_q8_physical_coordinates[0]
)

payload = {
    "schema": "elkies-k3.h92-q6-child-q8-physical-root-target.v2",
    "status": "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET",
    "inputs": {
        "source_frame": {"path": str(FRAME.relative_to(ROOT)), "sha256": digest(FRAME)},
        "marking": {"path": str(args.marking.relative_to(ROOT)), "sha256": digest(args.marking)},
        "q8_orbits": {"path": str(args.q8_orbits.relative_to(ROOT)), "sha256": digest(args.q8_orbits)},
    },
    "normalization": {
        "representative": args.representative,
        "representative_source": (
            "q8.nef_representative then q6 effective-component Weyl reduction" if args.representative == "component-nef"
            else "q8.nef_representative" if args.representative == "nef"
            else "q8.d13_mw4_hits[mw_projection=(0,-2,0)]"
        ),
        "effective_component_reflections": [[int(index), int(pairing)] for index, pairing in effective_reflections],
        "zero": "transported old zero on the explicit q6 child",
        "standard_weierstrass_zero_alignment": {
            "standard_zero": "the infinity section of the displayed short Weierstrass child",
            "marked_chord_status": (
                "The existing <1,m> marking is for the divisor translated from "
                "the transported old zero to the standard Weierstrass zero."
            ),
            "consequence": (
                "This untransformed physical target cannot be paired directly with "
                "that standard-coordinate chord until the section-translation action "
                "is transported in the exact NS marking."
            ),
        },
        "root_lattice": "orthogonal complement of <F, old_zero+F>",
        "component_labels": "pinned lattice-simple roots only; a resolved-chart map remains to be derived",
    },
    "selected_q8": {
        "identity": "q8_fibre=old_zero+marked_section+vertical",
        "source_h3_ns_vector": list(map(int, selected_q8)),
        "physical_old_fibre_coordinates": {
            "basis": "(old_fibre, old_zero+old_fibre, rank-17 complement)",
            "coordinates": list(map(int, selected_q8_physical_coordinates)),
            "positive_rank17_gram": [
                [int(value) for value in row] for row in physical_child.rows()
            ],
        },
        "old_zero_and_marked_section_component_degrees": {
            "E6": [0] * 6,
            "E8": [0] * 8,
        },
        "vertical_fibre_coefficient": int(coordinates[0]),
        "E6": {
            "simple_root_vectors_in_source_h3_ns": [list(map(int, row)) for row in e6_roots.rows()],
            "cartan": [[int(value) for value in row] for row in E6_CARTAN.rows()],
            "vertical_cycle": list(map(int, e6_cycle)),
            "component_degrees": list(map(int, e6_degrees)),
            "highest_root_coefficients": list(map(int, e6_highest)),
            "affine_component_degree": int(e6_affine_degree),
        },
        "E8": {
            "simple_root_vectors_in_source_h3_ns": [list(map(int, row)) for row in e8_roots.rows()],
            "cartan": [[int(value) for value in row] for row in E8_CARTAN.rows()],
            "vertical_cycle": list(map(int, e8_cycle)),
            "component_degrees": list(map(int, e8_degrees)),
            "highest_root_coefficients": list(map(int, e8_highest)),
            "affine_component_degree": int(e8_affine_degree),
        },
    },
    "boundary": (
        "This is the exact physical-component lattice target. It does not "
        "construct resolved II*/IV* coefficient modules, their finite "
        "quotient maps, the standard-zero translation of this divisor, a global "
        "q8 pencil, a bisection, an extension "
        "collision, or generic rank 18 or 19."
    ),
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "H92Q6CHILDQ8ROOTS|representative={}|E6_cycle={}|E6_degrees={}|"
    "E8_cycle={}|E8_degrees={}|status=PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET".format(
        args.representative,
        ",".join(map(str, e6_cycle)), ",".join(map(str, e6_degrees)),
        ",".join(map(str, e8_cycle)), ",".join(map(str, e8_degrees)),
    ),
    flush=True,
)
