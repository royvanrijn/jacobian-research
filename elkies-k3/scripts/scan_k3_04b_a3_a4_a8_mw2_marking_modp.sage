#!/usr/bin/env sage-python
"""Exhaust both pole-zero MW markings on the semistable MW2 source.

The determinant-500 source ``S2021`` has fibres ``I4+I5+I9`` and a physical
MW basis whose two generators both have pole order zero.  This script derives
their local correction classes directly from the root-adapted Gram matrix,
uses the local depth shared by the two inverse orientations of each
multiplicative component class, and exhausts the remaining degree-four X
coefficients over the finite field.  Polynomial Y
square roots and exact component depths are checked before two marked
generators are paired on the same fibre model.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import re
from pathlib import Path

from sage.all import GF, PolynomialRing, PowerSeriesRing, QQ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
DEFAULT_FIBRES = GEN / "elkies-k3-k3-04b86146cc6b284b-a3-a4-a8-mw2-fibre-ansatz-mod5-v1.json"
DEFAULT_SOURCES = GEN / "elkies-k3-k3-04b86146cc6b284b-prescribed-root-sources-large-a-v1.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-k3-04b86146cc6b284b-a3-a4-a8-mw2-marking-mod5-v1.json"
SOURCE_ID = "K3-04b86146cc6b284b-S2021"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def connected_components(gram):
    unseen = set(range(gram.nrows()))
    result = []
    while unseen:
        first = min(unseen)
        unseen.remove(first)
        todo = [first]
        component = []
        while todo:
            node = todo.pop()
            component.append(node)
            for other in tuple(unseen):
                if gram[node, other]:
                    unseen.remove(other)
                    todo.append(other)
        result.append(sorted(component))
    return result


def fibre_orders(root_type):
    orders = []
    for term in root_type.split("+"):
        match = re.fullmatch(r"(?:(\d+))?A(\d+)", term)
        if match is None:
            raise ValueError("source is not semistable")
        orders.extend([int(match.group(2)) + 1] * int(match.group(1) or 1))
    return orders


def formal_center(A, B, point, precision):
    field = A.base_ring()
    ring = A.parent()
    t = ring.gen()
    shifted_A = ring(A(t + point))
    shifted_B = ring(B(t + point))
    node = -field(3) * shifted_B[0] / (field(2) * shifted_A[0])
    series_ring = PowerSeriesRing(field, "s", default_prec=precision + 3)
    center = series_ring(node)
    series_A = series_ring(shifted_A)
    for unused in range(8):
        center = (center + (-series_A / 3) / center) / 2
    if (center**2 + series_A / 3).valuation() < precision + 1:
        raise ArithmeticError("formal center did not converge")
    return series_ring, center


def reversed_local(poly, weight, series_ring):
    u = series_ring.gen()
    return series_ring(
        sum(poly[index] * u ** (weight - index) for index in range(poly.degree() + 1))
    )


def polynomial_roots(right):
    if not right.is_square():
        return []
    positive = right.sqrt()
    return [positive] if not positive else [positive, -positive]


def serialize_poly(poly):
    return [int(value) for value in poly]


def depth_options(order, correction):
    return sorted(
        {
            min(depth, order - depth)
            for depth in range(order)
            if QQ(depth * (order - depth)) / order == correction
        }
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--fibres", type=Path, default=DEFAULT_FIBRES)
parser.add_argument("--sources", type=Path, default=DEFAULT_SOURCES)
parser.add_argument("--source-id", default=SOURCE_ID)
parser.add_argument("--quadratic-twist", type=int, default=1)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

fibres_path = arguments.fibres.resolve()
sources_path = arguments.sources.resolve()
output_path = arguments.output.resolve()
fibres = json.loads(fibres_path.read_text())
sources = json.loads(sources_path.read_text())
source_row = next(
    row for row in sources["sources"] if row["source_id"] == arguments.source_id
)
source = source_row["source"]
orders = fibre_orders(source["root_type"])
if orders != [4, 5, 9] or source["mw_rank_for_rho_19"] != 2:
    raise ValueError("marking scan requires the semistable A3+A4+A8/MW2 source")
if [row["pole_order"] for row in source["pole_audit"]["basis"]] != [0, 0]:
    raise ValueError("selected source no longer has a pole-zero MW basis")
if fibres["ansatz"]["normalized_reducible_supports"] != [
    "0:I4",
    "1:I5",
    "infinity:I9",
]:
    raise ValueError("fibre artifact does not have normalized I4+I5+I9 supports")
if not fibres["scan"]["exhausted"]:
    raise ValueError("marking scan requires an exhaustive fibre census")

root_rank = int(source["root_rank"])
root = matrix(QQ, source["root_adapted_gram"])[:root_rank, :root_rank]
components = connected_components(root)
component_by_rank = {len(component): component for component in components}
if set(component_by_rank) != {3, 4, 8}:
    raise ArithmeticError("unexpected root-component ranks")

generator_profiles = []
for basis_index, basis in enumerate(source["pole_audit"]["basis"]):
    labels = vector(QQ, basis["simple_root_pairings"])
    corrections = []
    options = []
    for order in orders:
        component = component_by_rank[order - 1]
        block = root.matrix_from_rows_and_columns(component, component)
        block_labels = vector(QQ, [labels[index] for index in component])
        correction = block_labels * block.inverse() * block_labels
        depths = depth_options(order, correction)
        if not depths:
            raise ArithmeticError("component correction has no I_n depth")
        corrections.append(correction)
        options.append(depths)
    generator_profiles.append(
        {
            "basis_index": basis_index,
            "height": source["mw_height_gram"][1 - basis_index][1 - basis_index],
            "component_corrections": [str(value) for value in corrections],
            "depth_profiles": [list(profile) for profile in itertools.product(*options)],
        }
    )

prime = int(fibres["prime"])
field = GF(prime)
twist = field(arguments.quadratic_twist)
if not twist:
    raise ValueError("quadratic twist must be nonzero")
ring = PolynomialRing(field, "t")
t = ring.gen()

accounting = {
    "affine_X_spaces": 0,
    "X_numerators_scanned": 0,
    "polynomial_Y_square_roots": 0,
    "marked_generator_sections": [0, 0],
    "models_with_both_generator_section_classes": 0,
    "models_with_marked_basis": 0,
    "marked_ordered_basis_pairs": 0,
    "component_matched_pair_candidates": 0,
    "pairs_meeting_singular_fibres": 0,
    "pairs_with_wrong_smooth_intersection": 0,
}
depth_histograms = [{}, {}]
records = []
for example_index, example in enumerate(fibres["examples"]):
    A = twist**2 * ring(example["A_coefficients_low_to_high"])
    B = twist**3 * ring(example["B_coefficients_low_to_high"])
    discriminant_core = 4 * A**3 + 27 * B**2
    zero_ring, zero_center = formal_center(A, B, field.zero(), 5)
    one_ring, one_center = formal_center(A, B, field.one(), 5)
    infinity_ring = PowerSeriesRing(field, "u", default_prec=12)
    infinity_A = reversed_local(A, 8, infinity_ring)
    infinity_B = reversed_local(B, 12, infinity_ring)
    infinity_center = infinity_ring(-field(3) * infinity_B[0] / (field(2) * infinity_A[0]))
    for unused in range(8):
        infinity_center = (
            infinity_center + (-infinity_A / 3) / infinity_center
        ) / 2
    if (infinity_center**2 + infinity_A / 3).valuation() < 10:
        raise ArithmeticError("infinity formal center did not converge")

    centers = (zero_center, one_center, infinity_center)
    series_rings = (zero_ring, one_ring, infinity_ring)

    def local_x(poly, support_index):
        if support_index == 0:
            return zero_ring(poly(t))
        if support_index == 1:
            return one_ring(poly(t + 1))
        return reversed_local(poly, 4, infinity_ring)

    def local_y(poly, support_index):
        if support_index == 0:
            return zero_ring(poly(t))
        if support_index == 1:
            return one_ring(poly(t + 1))
        return reversed_local(poly, 6, infinity_ring)

    sections_by_generator = [[], []]
    for generator_index, generator in enumerate(generator_profiles):
        seen_sections = set()
        for depths in generator["depth_profiles"]:
            rows = []
            targets = []
            for support_index, depth in enumerate(depths):
                for jet in range(depth):
                    rows.append(
                        [
                            local_x(t**degree, support_index)[jet]
                            for degree in range(5)
                        ]
                    )
                    targets.append(centers[support_index][jet])
            linear = matrix(field, rows) if rows else matrix(field, 0, 5)
            target = vector(field, targets)
            if linear.rank() != linear.augment(target.column()).rank():
                continue
            particular = linear.solve_right(target)
            kernel_basis = linear.right_kernel().basis()
            accounting["affine_X_spaces"] += 1
            for coefficients in itertools.product(field, repeat=len(kernel_basis)):
                solution = particular + sum(
                    (coefficient * basis for coefficient, basis in zip(coefficients, kernel_basis)),
                    vector(field, 5),
                )
                X = ring(list(solution))
                accounting["X_numerators_scanned"] += 1
                right = X**3 + A * X + B
                for Y in polynomial_roots(right):
                    accounting["polynomial_Y_square_roots"] += 1
                    actual_depths = [
                        int(
                            min(
                                (local_x(X, support_index) - centers[support_index]).valuation(),
                                local_y(Y, support_index).valuation(),
                            )
                        )
                        for support_index in range(3)
                    ]
                    key = ",".join(map(str, actual_depths))
                    histogram = depth_histograms[generator_index]
                    histogram[key] = histogram.get(key, 0) + 1
                    if actual_depths != depths:
                        continue
                    section_key = (tuple(X.list()), tuple(Y.list()), tuple(depths))
                    if section_key in seen_sections:
                        continue
                    seen_sections.add(section_key)
                    sections_by_generator[generator_index].append(
                        {
                            "X_coefficients_low_to_high": serialize_poly(X),
                            "Y_coefficients_low_to_high": serialize_poly(Y),
                            "component_depths_at_I4_I5_I9": depths,
                        }
                    )
        accounting["marked_generator_sections"][generator_index] += len(
            sections_by_generator[generator_index]
        )

    basis_pairs = []
    if sections_by_generator[0] and sections_by_generator[1]:
        accounting["models_with_both_generator_section_classes"] += 1
    for left_index, left in enumerate(sections_by_generator[0]):
        for right_index, right in enumerate(sections_by_generator[1]):
            accounting["component_matched_pair_candidates"] += 1
            left_X = ring(left["X_coefficients_low_to_high"])
            left_Y = ring(left["Y_coefficients_low_to_high"])
            right_X = ring(right["X_coefficients_low_to_high"])
            right_Y = ring(right["Y_coefficients_low_to_high"])
            common = (left_X - right_X).gcd(left_Y - right_Y)
            if common.gcd(discriminant_core).degree() != 0:
                accounting["pairs_meeting_singular_fibres"] += 1
                continue
            intersection = int(common.degree())
            if intersection != 1:
                accounting["pairs_with_wrong_smooth_intersection"] += 1
                continue
            basis_pairs.append(
                {
                    "left_generator_section_index": left_index,
                    "right_generator_section_index": right_index,
                    "intersection_on_smooth_fibres": intersection,
                    "mw_index_from_height_determinant": 1,
                }
            )
    pair_count = len(basis_pairs)
    if pair_count:
        accounting["models_with_marked_basis"] += 1
        accounting["marked_ordered_basis_pairs"] += pair_count
        records.append(
            {
                "example_index": example_index,
                "generator_sections": sections_by_generator,
                "marked_basis_pairs": basis_pairs,
                "ordered_basis_pair_count": pair_count,
            }
        )

status = (
    "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_WITH_MARKED_MW2_BASIS"
    if accounting["marked_ordered_basis_pairs"]
    else "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_EMPTY_MARKED_MW2_BASIS_LOCUS"
)
payload = {
    "schema": "elkies-k3.k3-04b-a3-a4-a8-mw2-marking-modp.v1",
    "status": status,
    "prime": prime,
    "quadratic_twist": int(twist),
    "quadratic_twist_square_class": "square" if twist.is_square() else "nonsquare",
    "source": {
        "surface_id": "K3-04b86146cc6b284b",
        "source_id": arguments.source_id,
        "source_gram_sha256": source["gram_sha256"],
        "root_type": source["root_type"],
        "mw_height_gram": source["mw_height_gram"],
        "minimum_basis_pole_profile": [0, 0],
        "generator_profiles": generator_profiles,
    },
    "scope": {
        "fibre_census_exhaustive": True,
        "fibre_models": len(fibres["examples"]),
        "all_component_class_local_depths": True,
        "all_remaining_degree_four_X_coefficients": True,
        "all_polynomial_Y_square_roots": True,
        "all_component_matched_pairs_tested_at_required_smooth_intersection": 1,
    },
    "accounting": accounting,
    "exact_depth_histograms_for_square_sections": [
        dict(sorted(histogram.items())) for histogram in depth_histograms
    ],
    "models": records,
    "inputs": {
        relative(fibres_path): digest(fibres_path),
        relative(sources_path): digest(sources_path),
    },
    "proof_boundary": {
        "proved": (
            "For every model in the exhaustive normalized fibre census, both "
            "local depths for every prescribed multiplicative component class "
            "and every remaining degree-four X coefficient were exhausted. "
            "Every retained ordered pair consists of polynomial sections with exact "
            "component depths in the two distinct lattice generator classes."
        ),
        "not_proved": (
            "A finite-field marked basis need not lift to characteristic zero or be "
            "defined over Q. The scan does not certify MW saturation, a rational "
            "parameterization, or a physical neighbour corridor."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/scan_k3_04b_a3_a4_a8_mw2_marking_modp.sage"
        + (
            f" --quadratic-twist {int(twist)}" if twist != 1 else ""
        )
        + (
            f" --output {relative(output_path)}"
            if output_path != DEFAULT_OUTPUT.resolve()
            else ""
        )
    ),
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != serialized:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    f"K304BMW2MARKING|p={prime}|twist={int(twist)}|models={len(fibres['examples'])}|"
    f"basis_pairs={accounting['marked_ordered_basis_pairs']}|"
    f"status={'PASS' if status.startswith('PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_WITH') else 'EMPTY'}",
    flush=True,
)
