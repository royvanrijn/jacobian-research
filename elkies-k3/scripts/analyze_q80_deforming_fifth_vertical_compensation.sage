#!/usr/bin/env sage
"""Decompose the two deforming fifth-q4 CM24 divisor classes.

The bounded generic q4 search specializes to two chamber-reduced classes on
the compact CM24 fourth child (old-fiber degrees 43 and 47).  Both have the
same horizontal MW coordinate (1,0).  This checker reconstructs the unique
effective section in that MW class and writes each isotropic divisor as

    D = O + P + m F + sum(c_i Theta_i).

The vertical coefficients are the exact input needed for an equation-level
compensated projection.  This is a lattice decomposition certificate, not a
construction of the corresponding rational function.
"""

import json
import sys
from pathlib import Path

from sage.all import QQ, ZZ, matrix, vector
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

# The readiness checker owns a small argparse parser in a loaded dependency.
saved_argv = list(sys.argv)
try:
    sys.argv = [sys.argv[0]]
    load(str(HERE / "analyze_q80_fifth_q4_cm24_readiness.sage"))
finally:
    sys.argv = saved_argv


def effective_sections(mw_coordinates):
    """Recover every torsion lift of the requested free MW class."""
    mw_coordinates = vector(ZZ, mw_coordinates)
    raw_lift = vector(ZZ, mw_coordinates * new_mw_basis_lifts)
    projection = new_project_mw(raw_lift)
    height = projection * special_fourth * projection
    # The special fourth frame has order-three torsion.  The ordinary ADE root
    # lattice sees only one of the three lifts; its saturation sees all of
    # them while preserving the same free MW projection.
    saturated_roots = new_simple.row_module(ZZ).saturation().basis_matrix()
    assert saturated_roots.nrows() == new_simple.nrows()
    assert new_simple.row_module(ZZ).index_in(saturated_roots.row_module(ZZ)) == 3
    saturated_gram = saturated_roots * special_fourth * saturated_roots.transpose()
    raw_root_coordinates = (
        vector(QQ, raw_lift)
        * special_fourth
        * saturated_roots.transpose()
        * saturated_gram.inverse()
    )
    root_lattice = IntegralLattice(saturated_gram)
    iterator = root_lattice.enumerate_close_vectors(-raw_root_coordinates)
    matches = []
    for _ in range(32768):
        shift = vector(ZZ, next(iterator))
        lift = raw_lift + shift * saturated_roots
        norm = ZZ(lift * special_fourth * lift)
        if norm > 4:
            break
        if norm != 4:
            continue
        section = vector(ZZ, [1, 1] + list(lift))
        if any(
            intersection(section, curve, new_ns) < 0
            for _, curve in new_curves[1:]
        ):
            continue
        matches.append(section)
    assert len(matches) == 3, len(matches)
    for section in matches:
        assert section * new_ns * section == -2
        assert intersection(section, new_fiber, new_ns) == 1
        assert intersection(section, new_zero, new_ns) == 0
        assert new_project_mw(section[2:]) == projection
    return tuple(matches), height, ZZ(0)


sections, section_height, section_pole = effective_sections((1, 0))
affine_data = highest_roots(special_fourth, new_simple, new_positive)
assert sorted(len(component) for component, _, _ in affine_data) == [2, 2, 2, 2, 3, 5]
component_rows = root_component_data(special_fourth)


def discriminant_profile(section):
    """Return exact ADE discriminant classes, independent of chain orientation."""
    return tuple(
        fractional_root_class(special_fourth, section[2:], component_basis)
        for _, component_basis in component_rows
    )


def add_profiles(left, right):
    return tuple(
        tuple(
            value - value.floor()
            for value in (vector(QQ, left_row) + vector(QQ, right_row))
        )
        for left_row, right_row in zip(left, right)
    )


def profile_record(profile):
    return [[str(value) for value in component] for component in profile]


known_section_classes = tuple(vector(ZZ, row[2]) for row in candidate_rows)
known_section_labels = tuple(discriminant_profile(row) for row in known_section_classes)
pair_indices = ((0, 1), (1, 4), (2, 3))
pair_labels = {}
for left, right in pair_indices:
    pair_labels[left, right] = add_profiles(
        known_section_labels[left], known_section_labels[right]
    )
torsion_lift_labels = tuple(discriminant_profile(row) for row in sections)
torsion_lift_pair_matches = tuple(
    tuple(pair for pair, labels in pair_labels.items() if labels == target)
    for target in torsion_lift_labels
)
assert sorted(sum((list(row) for row in torsion_lift_pair_matches), [])) == sorted(pair_indices)

window_paths = sorted(
    ROOT.glob("artifacts/local/q80-q4-deforming-window-*.json")
)
assert len(window_paths) == 3
all_candidates = []
for path in window_paths:
    payload = json.loads(path.read_text())
    all_candidates.extend(payload["candidates"])

representatives = {}
for row in all_candidates:
    degree = int(row["old_fiber_degree"])
    representatives.setdefault(degree, row)
assert set(representatives) == {43, 47}

rows = []
for degree in sorted(representatives):
    source = representatives[degree]
    divisor = vector(ZZ, source["cm24_reduced"])
    assert divisor * new_ns * divisor == 0
    assert intersection(divisor, new_fiber, new_ns) == 2
    assert intersection(divisor, new_zero, new_ns) == 0
    choices = []
    for section_index, possible_section in enumerate(sections):
        possible_vertical = divisor - new_zero - possible_section
        try:
            possible_coordinates = new_simple.solve_left(
                vector(ZZ, possible_vertical[2:])
            )
        except ValueError:
            continue
        if not all(value in ZZ for value in possible_coordinates):
            continue
        choices.append(
            (section_index, possible_section, possible_vertical,
             vector(ZZ, possible_coordinates))
        )
    assert len(choices) == 1, (degree, len(choices))
    section_index, section, vertical, root_coordinates = choices[0]
    assert vertical[1] == 0
    assert root_coordinates * new_simple == vertical[2:]

    per_component = []
    allocated_fibers = 0
    for component, _, highest_coordinates in affine_data:
        coefficients = vector(ZZ, [root_coordinates[index] for index in component])
        marks = vector(ZZ, [highest_coordinates[index] for index in component])
        # Add the least whole copy of this reducible fiber which makes its
        # affine and nonidentity component multiplicities nonnegative.
        fiber_copies = max(
            [ZZ(0)]
            + [(-coefficient + mark - 1) // mark
               for coefficient, mark in zip(coefficients, marks)]
        )
        effective_coefficients = coefficients + fiber_copies * marks
        assert all(value >= 0 for value in effective_coefficients)
        allocated_fibers += fiber_copies
        per_component.append({
            "type": f"A{len(component)}",
            "simple_root_indices": list(map(int, component)),
            "highest_root_marks": list(map(int, marks)),
            "fiber_copies": int(fiber_copies),
            "affine_coefficient": int(fiber_copies),
            "coefficients": list(map(int, coefficients)),
            "effective_simple_coefficients": list(map(int, effective_coefficients)),
        })
    residual_fibers = ZZ(vertical[0]) - allocated_fibers
    assert residual_fibers >= 0

    record = {
        "old_fiber_degree": degree,
        "generic_neighbor_v": source["v"],
        "divisor": list(map(int, divisor)),
        "section": list(map(int, section)),
        "torsion_lift_index": int(section_index),
        "torsion_component_labels": profile_record(
            torsion_lift_labels[section_index]
        ),
        "explicit_section_pair_matches": [
            list(map(int, pair))
            for pair in torsion_lift_pair_matches[section_index]
        ],
        "known_section_pairings": [
            int(intersection(divisor, known_section, new_ns))
            for known_section in known_section_classes
        ],
        "fiber_coefficient": int(vertical[0]),
        "residual_free_fibers": int(residual_fibers),
        "simple_root_coefficients": list(map(int, root_coordinates)),
        "component_coefficients": per_component,
        "vertical_square": int(vertical * new_ns * vertical),
        "vertical_dot_section": int(vertical * new_ns * section),
        "vertical_dot_zero": int(vertical * new_ns * new_zero),
    }
    rows.append(record)
    print(
        "Q80DEFORMINGQ4VERTICAL|"
        f"old_degree={degree}|fiber_coefficient={record['fiber_coefficient']}|"
        f"residual_free_fibers={record['residual_free_fibers']}|"
        f"simple_root_coefficients={tuple(root_coordinates)}|"
        f"components={tuple((item['type'], item['fiber_copies'], tuple(item['effective_simple_coefficients'])) for item in per_component)}|"
        f"vertical_square={record['vertical_square']}|"
        f"vertical_dot_P={record['vertical_dot_section']}|"
        f"vertical_dot_O={record['vertical_dot_zero']}|status=PASS",
        flush=True,
    )

artifact = {
    "schema": "q80-deforming-fifth-q4-vertical-compensation-v1",
    "status": "exact_cm24_lattice_decomposition",
    "source_fibration": "A5+A3+4A2/MW2",
    "horizontal_mw_coordinates": [1, 0],
    "horizontal_height": str(section_height),
    "horizontal_P_dot_O": int(section_pole),
    "known_section_component_labels": [
        profile_record(labels) for labels in known_section_labels
    ],
    "torsion_lift_component_labels": [
        profile_record(labels) for labels in torsion_lift_labels
    ],
    "torsion_lift_pair_matches": [
        [list(map(int, pair)) for pair in matches]
        for matches in torsion_lift_pair_matches
    ],
    "rows": rows,
    "rank_claim": None,
    "reproduce": (
        "sage elkies-k3/scripts/"
        "analyze_q80_deforming_fifth_vertical_compensation.sage"
    ),
}
output = (
    ROOT / "artifacts/local/"
    "q80-deforming-fifth-q4-vertical-compensation.json"
)
output.write_text(
    json.dumps(artifact, indent=2, sort_keys=True, default=int) + "\n"
)
print(
    "Q80DEFORMINGQ4VERTICAL|"
    f"section_height={section_height}|section_P.O={section_pole}|"
    f"artifact={output}|status=PASS_EXACT_DECOMPOSITION",
    flush=True,
)
