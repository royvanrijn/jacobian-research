#!/usr/bin/env sage-python
"""Exhaust a pole-[0,1] MW2 marking on normalized I6+I11 models.

The determinant-384 source is the pinned default.  A different A5+A10/MW2
source can reuse the exhaustive fibre census; its component depths, physical
basis height, and required smooth intersection are derived from its lattice
data rather than supplied as search parameters.
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
DEFAULT_FIBRES = GEN / "elkies-k3-k3-6ce16abb9de3c7c5-a5-a10-mw2-fibre-ansatz-mod5-v1.json"
DEFAULT_SOURCES = GEN / "elkies-k3-k3-6ce16abb9de3c7c5-semistable-mw0-2-sources-large-a-partner1-v1.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-k3-6ce16abb9de3c7c5-a5-a10-mw2-marking-mod5-v1.json"
SOURCE_ID = "K3-6ce16abb9de3c7c5-S0008"
SURFACE_ID = "K3-6ce16abb9de3c7c5"
SCHEMA = "elkies-k3.k3-6ce-a5-a10-mw2-marking-modp.v1"


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def digest(path):
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
    return sorted(result, key=len)


def two_support_orders(root_type):
    ranks = []
    for term in root_type.split("+"):
        match = re.fullmatch(r"A(\d+)", term)
        if match is None:
            raise ValueError("source is not a two-A-component semistable frame")
        ranks.append(int(match.group(1)))
    if len(ranks) != 2 or sum(ranks) != 15:
        raise ValueError("source does not have two A components of total rank 15")
    return [rank + 1 for rank in sorted(ranks)]


def formal_center_at_zero(A, B, precision):
    field = A.base_ring()
    series_ring = PowerSeriesRing(field, "s", default_prec=precision + 3)
    node = -field(3) * B[0] / (field(2) * A[0])
    center = series_ring(node)
    series_A = series_ring(A)
    for unused in range(8):
        center = (center + (-series_A / 3) / center) / 2
    if (center**2 + series_A / 3).valuation() < precision + 1:
        raise ArithmeticError("zero formal center did not converge")
    return series_ring, center


def reversed_local(poly, weight, series_ring):
    u = series_ring.gen()
    return series_ring(
        sum(poly[index] * u ** (weight - index) for index in range(poly.degree() + 1))
    )


def infinity_center(A, B, precision):
    field = A.base_ring()
    series_ring = PowerSeriesRing(field, "u", default_prec=precision + 3)
    local_A = reversed_local(A, 8, series_ring)
    local_B = reversed_local(B, 12, series_ring)
    center = series_ring(-field(3) * local_B[0] / (field(2) * local_A[0]))
    for unused in range(8):
        center = (center + (-local_A / 3) / center) / 2
    if (center**2 + local_A / 3).valuation() < precision + 1:
        raise ArithmeticError("infinity formal center did not converge")
    return series_ring, center


def polynomial_roots(right):
    if not right.is_square():
        return []
    positive = right.sqrt()
    return [positive] if not positive else [positive, -positive]


def serialize(poly):
    return [int(value) for value in poly]


def affine_solutions(field, linear, target):
    if linear.rank() != linear.augment(target.column()).rank():
        return
    particular = linear.solve_right(target)
    kernel = linear.right_kernel().basis()
    for coefficients in itertools.product(field, repeat=len(kernel)):
        yield particular + sum(
            (coefficient * basis for coefficient, basis in zip(coefficients, kernel)),
            vector(field, linear.ncols()),
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fibres", type=Path, default=DEFAULT_FIBRES)
    parser.add_argument("--sources", type=Path, default=DEFAULT_SOURCES)
    parser.add_argument("--source-id", default=SOURCE_ID)
    parser.add_argument("--surface-id", default=SURFACE_ID)
    parser.add_argument("--schema", default=SCHEMA)
    parser.add_argument("--quadratic-twist", type=int, default=1)
    parser.add_argument(
        "--max-models", type=int, default=0,
        help="test only the first N fibre models; zero exhausts the census",
    )
    parser.add_argument(
        "--skip-models", type=int, default=0,
        help="skip this many fibre models before applying --max-models",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if arguments.max_models < 0 or arguments.skip_models < 0:
        parser.error("model bounds must be nonnegative")

    fibres_path = arguments.fibres.resolve()
    sources_path = arguments.sources.resolve()
    output_path = arguments.output.resolve()
    legacy_profile = (
        fibres_path == DEFAULT_FIBRES.resolve()
        and sources_path == DEFAULT_SOURCES.resolve()
        and arguments.source_id == SOURCE_ID
        and arguments.surface_id == SURFACE_ID
        and arguments.schema == SCHEMA
    )
    fibres = json.loads(fibres_path.read_text())
    sources = json.loads(sources_path.read_text())
    source = next(
        row["source"] for row in sources["sources"] if row["source_id"] == arguments.source_id
    )
    basis = source["pole_audit"]["basis"]
    orders = two_support_orders(source["root_type"])
    expected_supports = [f"0:I{orders[0]}", f"infinity:I{orders[1]}"]
    if fibres["ansatz"]["normalized_reducible_supports"] != expected_supports:
        raise ValueError("fibre census does not match the selected source supports")
    if not fibres["scan"]["exhausted"]:
        raise ValueError("marking scan requires an exhaustive fibre census")
    if not (
        source["mw_rank_for_rho_19"] == 2
        and [section["pole_order"] for section in basis] == [0, 1]
    ):
        raise ValueError("selected source no longer has the required MW2 basis")

    root = matrix(QQ, source["root_adapted_gram"])[:15, :15]
    components = connected_components(root)
    if list(map(len, components)) != [order - 1 for order in orders]:
        raise ArithmeticError("unexpected root components")
    depth_profiles = []
    for section in basis:
        labels = vector(QQ, section["simple_root_pairings"])
        depths = []
        for component in components:
            block = root.matrix_from_rows_and_columns(component, component)
            local = vector(QQ, [labels[index] for index in component])
            correction = local * block.inverse() * local
            order = len(component) + 1
            options = [
                depth
                for depth in range(order // 2 + 1)
                if QQ(depth * (order - depth)) / order == correction
            ]
            if len(options) != 1:
                raise ArithmeticError("ambiguous component depth")
            depths.append(options[0])
        depth_profiles.append(depths)
    labels = [vector(QQ, section["simple_root_pairings"]) for section in basis]
    component_cross = QQ(0)
    for component in components:
        block = root.matrix_from_rows_and_columns(component, component)
        local = [vector(QQ, [row[index] for index in component]) for row in labels]
        component_cross += local[0] * block.inverse() * local[1]
    height = matrix(QQ, source["mw_height_gram"])
    coordinates = matrix(QQ, [section["mw_quotient_coordinates"] for section in basis])
    if abs(coordinates.det()) != 1:
        raise ArithmeticError("selected sections do not form a primitive MW basis")
    basis_height = coordinates * height * coordinates.transpose()
    required_intersection = QQ(3) - component_cross - basis_height[0, 1]
    if required_intersection.denominator() != 1 or required_intersection < 0:
        raise ArithmeticError("lattice data gives an invalid smooth intersection")
    required_intersection = int(required_intersection)
    depth_key = (
        "component_depths_at_I6_I11"
        if legacy_profile
        else "component_depths_at_normalized_supports"
    )

    prime = int(fibres["prime"])
    field = GF(prime)
    twist = field(arguments.quadratic_twist)
    if not twist:
        raise ValueError("quadratic twist must be nonzero")
    ring = PolynomialRing(field, "t")
    t = ring.gen()

    accounting = {
        "pole_zero_X_numerators_scanned": 0,
        "pole_one_denominators_scanned": 0,
        "pole_one_X_numerators_scanned": 0,
        "marked_generator_sections": [0, 0],
        "models_with_both_generator_section_classes": 0,
        "component_matched_pair_candidates": 0,
        "pairs_meeting_singular_fibres": 0,
        "pairs_with_wrong_smooth_intersection": 0,
        "models_with_marked_basis": 0,
        "marked_ordered_basis_pairs": 0,
    }
    depth_histograms = [{}, {}]
    records = []
    available_examples = fibres["examples"][arguments.skip_models :]
    selected_examples = (
        available_examples[: arguments.max_models]
        if arguments.max_models
        else available_examples
    )
    for local_index, example in enumerate(selected_examples):
        example_index = arguments.skip_models + local_index
        A = twist**2 * ring(example["A_coefficients_low_to_high"])
        B = twist**3 * ring(example["B_coefficients_low_to_high"])
        discriminant_core = 4 * A**3 + 27 * B**2
        zero_ring, zero_center = formal_center_at_zero(A, B, 3)
        infinity_ring, inf_center = infinity_center(A, B, 6)

        sections = [[], []]

        # Pole-zero generator: constrain the lattice-prescribed jets at both supports.
        rows0 = []
        targets0 = []
        for support_index, depth in enumerate(depth_profiles[0]):
            for jet in range(depth):
                rows0.append(
                    [
                        (
                            zero_ring(t**degree)
                            if support_index == 0
                            else reversed_local(t**degree, 4, infinity_ring)
                        )[jet]
                        for degree in range(5)
                    ]
                )
                targets0.append(
                    (zero_center if support_index == 0 else inf_center)[jet]
                )
        linear0 = matrix(field, rows0) if rows0 else matrix(field, 0, 5)
        target0 = vector(field, targets0)
        for solution in affine_solutions(field, linear0, target0):
            accounting["pole_zero_X_numerators_scanned"] += 1
            X = ring(list(solution))
            for Y in polynomial_roots(X**3 + A * X + B):
                x_zero = zero_ring(X)
                y_zero = zero_ring(Y)
                x_inf = reversed_local(X, 4, infinity_ring)
                y_inf = reversed_local(Y, 6, infinity_ring)
                depths = [
                    int(min((x_zero - zero_center).valuation(), y_zero.valuation())),
                    int(min((x_inf - inf_center).valuation(), y_inf.valuation())),
                ]
                key = ",".join(map(str, depths))
                depth_histograms[0][key] = depth_histograms[0].get(key, 0) + 1
                if depths == depth_profiles[0]:
                    sections[0].append(
                        {
                            "X_coefficients_low_to_high": serialize(X),
                            "Y_coefficients_low_to_high": serialize(Y),
                            depth_key: depths,
                        }
                    )

        # Pole-one generator: exhaust every monic linear denominator away from t=0.
        for c0 in field:
            C = t + c0
            if not C(0):
                continue
            accounting["pole_one_denominators_scanned"] += 1
            local_C = reversed_local(C, 1, infinity_ring)
            zero_C = zero_ring(C)
            rows1 = []
            targets1 = []
            for support_index, depth in enumerate(depth_profiles[1]):
                for jet in range(depth):
                    rows1.append(
                        [
                            (
                                zero_ring(t**degree) / zero_C**2
                                if support_index == 0
                                else reversed_local(t**degree, 6, infinity_ring)
                                / local_C**2
                            )[jet]
                            for degree in range(7)
                        ]
                    )
                    targets1.append(
                        (zero_center if support_index == 0 else inf_center)[jet]
                    )
            linear1 = matrix(field, rows1) if rows1 else matrix(field, 0, 7)
            target1 = vector(field, targets1)
            for solution in affine_solutions(field, linear1, target1):
                accounting["pole_one_X_numerators_scanned"] += 1
                Xn = ring(list(solution))
                if Xn.gcd(C).degree() != 0:
                    continue
                right = Xn**3 + A * Xn * C**4 + B * C**6
                for Yn in polynomial_roots(right):
                    if Yn.gcd(C).degree() != 0:
                        continue
                    x_zero = zero_ring(Xn) / zero_ring(C) ** 2
                    y_zero = zero_ring(Yn) / zero_ring(C) ** 3
                    x_inf = reversed_local(Xn, 6, infinity_ring) / local_C**2
                    y_inf = reversed_local(Yn, 9, infinity_ring) / local_C**3
                    depths = [
                        int(min((x_zero - zero_center).valuation(), y_zero.valuation())),
                        int(min((x_inf - inf_center).valuation(), y_inf.valuation())),
                    ]
                    key = ",".join(map(str, depths))
                    depth_histograms[1][key] = depth_histograms[1].get(key, 0) + 1
                    if depths == depth_profiles[1]:
                        sections[1].append(
                            {
                                "C_coefficients_low_to_high": serialize(C),
                                "X_numerator_coefficients_low_to_high": serialize(Xn),
                                "Y_numerator_coefficients_low_to_high": serialize(Yn),
                                depth_key: depths,
                            }
                        )

        accounting["marked_generator_sections"][0] += len(sections[0])
        accounting["marked_generator_sections"][1] += len(sections[1])
        if sections[0] and sections[1]:
            accounting["models_with_both_generator_section_classes"] += 1
        pairs = []
        for left_index, left in enumerate(sections[0]):
            X0 = ring(left["X_coefficients_low_to_high"])
            Y0 = ring(left["Y_coefficients_low_to_high"])
            for right_index, right in enumerate(sections[1]):
                accounting["component_matched_pair_candidates"] += 1
                C = ring(right["C_coefficients_low_to_high"])
                X1 = ring(right["X_numerator_coefficients_low_to_high"])
                Y1 = ring(right["Y_numerator_coefficients_low_to_high"])
                common = (X0 * C**2 - X1).gcd(Y0 * C**3 - Y1)
                if common.gcd(discriminant_core).degree() != 0:
                    accounting["pairs_meeting_singular_fibres"] += 1
                    continue
                intersection = int(common.degree())
                if intersection != required_intersection:
                    accounting["pairs_with_wrong_smooth_intersection"] += 1
                    continue
                pairs.append(
                    {
                        "left_generator_section_index": left_index,
                        "right_generator_section_index": right_index,
                        "intersection_on_smooth_fibres": intersection,
                        "mw_index_from_height_determinant": 1,
                    }
                )
        if pairs:
            accounting["models_with_marked_basis"] += 1
            accounting["marked_ordered_basis_pairs"] += len(pairs)
            records.append(
                {
                    "example_index": example_index,
                    "generator_sections": sections,
                    "marked_basis_pairs": pairs,
                    "ordered_basis_pair_count": len(pairs),
                }
            )

    status = (
        (
            "PASS_BOUNDED_NORMALIZED_CHART_WITH_MARKED_MW2_BASIS"
            if accounting["marked_ordered_basis_pairs"]
            else "PASS_BOUNDED_NORMALIZED_CHART_EMPTY_MARKED_MW2_BASIS_LOCUS"
        )
        if arguments.max_models or arguments.skip_models
        else (
            "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_WITH_MARKED_MW2_BASIS"
            if accounting["marked_ordered_basis_pairs"]
            else "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_EMPTY_MARKED_MW2_BASIS_LOCUS"
        )
    )
    payload = {
        "schema": arguments.schema,
        "status": status,
        "prime": prime,
        "quadratic_twist": int(twist),
        "quadratic_twist_square_class": "square" if twist.is_square() else "nonsquare",
        "source": {
            "surface_id": arguments.surface_id,
            "source_id": arguments.source_id,
            "source_gram_sha256": source["gram_sha256"],
            "root_type": source["root_type"],
            "mw_height_gram": source["mw_height_gram"],
            "minimum_basis_pole_profile": [0, 1],
            depth_key: depth_profiles,
            "required_smooth_pair_intersection": required_intersection,
            **(
                {}
                if legacy_profile
                else {
                    "physical_basis_height_gram": [
                        [str(value) for value in row] for row in basis_height.rows()
                    ],
                    "component_cross_correction": str(component_cross),
                }
            ),
        },
        "scope": {
            "fibre_census_exhaustive": not arguments.max_models and not arguments.skip_models,
            "fibre_models": len(selected_examples),
            "all_pole_zero_degree_four_X_numerators": True,
            "all_monic_linear_pole_denominators_away_from_reducible_supports": True,
            "all_pole_one_degree_six_X_numerators": True,
            "all_polynomial_Y_square_roots": True,
            "all_component_matched_pairs_tested_at_required_smooth_intersection": (
                required_intersection
            ),
        }
        | (
            {"fibre_model_offset": arguments.skip_models}
            if arguments.skip_models
            else {}
        ),
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
                "On every model in the exhaustive normalized fibre census, the full "
                "pole-zero and pole-one section coefficient spaces, polynomial square "
                "roots, exact component depths, and required pair intersection are tested."
            ),
            "not_proved": (
                "A finite-field marked basis is not a characteristic-zero lift, rational "
                + (
                    "marking, primitive determinant-384 specialization, or neighbour route."
                    if legacy_profile
                    else "marking, primitive characteristic-zero specialization, or neighbour route."
                )
            ),
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/scan_k3_6ce_a5_a10_mw2_marking_modp.sage"
            + (f" --fibres {relative(fibres_path)}" if fibres_path != DEFAULT_FIBRES.resolve() else "")
            + (f" --sources {relative(sources_path)}" if sources_path != DEFAULT_SOURCES.resolve() else "")
            + (f" --source-id {arguments.source_id}" if arguments.source_id != SOURCE_ID else "")
            + (f" --surface-id {arguments.surface_id}" if arguments.surface_id != SURFACE_ID else "")
            + (f" --schema {arguments.schema}" if arguments.schema != SCHEMA else "")
            + (f" --quadratic-twist {int(twist)}" if twist != 1 else "")
            + (f" --max-models {arguments.max_models}" if arguments.max_models else "")
            + (f" --skip-models {arguments.skip_models}" if arguments.skip_models else "")
            + (f" --output {relative(output_path)}" if output_path != DEFAULT_OUTPUT.resolve() else "")
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
        "K36CEA5A10MARKING|"
        f"p={prime}|twist={int(twist)}|models={len(selected_examples)}|"
        f"sections={accounting['marked_generator_sections']}|"
        f"pairs={accounting['marked_ordered_basis_pairs']}|"
        f"status={'PASS' if accounting['marked_ordered_basis_pairs'] else 'EMPTY'}",
        flush=True,
    )


if __name__ == "__main__":
    main()
