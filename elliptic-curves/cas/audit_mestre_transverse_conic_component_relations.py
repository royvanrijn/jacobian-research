#!/usr/bin/env python3
"""Exact rank-neutrality audit at the seed of the second rational component.

The point ``s=-357/47, T=1`` is the normalized six-root seed
``(0,7,79,81,128,137)``.  This checker records exact group-law expressions
for the two new affine points in terms of the visible points and an
independent mod-3 finite-reduction quotient computation.  It is deliberately
only a specialization statement.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from alternate_quartic_covers import short_add
from mestre_root_tuples import SixRootMestreConstruction
from search_mestre_dsquare_four import rational_square_root
from search_mestre_root_tuple_scale import (
    primitive_visible_points,
    quartic_point_to_jacobian,
    quartic_value,
)
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate
from verify_mestre_transverse_two_section_conic_component import (
    component_coordinates,
    split_roots,
)


Q = Fraction
S = -Q(357, 47)
T = Q(1)
PRIME_BOUND = 499
EXPECTED_COEFFICIENTS = (
    Q(0), Q(0), Q(0), Q(-25796475290575965843), Q(7775320562002058135094691758)
)
EXPECTED_PIVOTS = [1, 2, 3, 4, 5, 7, 8, 9, 10]


def rational_text(value: Fraction) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def negate(point: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return point[0], -point[1]


def negative_sum(
    coefficients: tuple[Fraction, ...],
    points: tuple[tuple[Fraction, Fraction], ...],
    indices: tuple[int, ...],
) -> tuple[Fraction, Fraction] | None:
    answer = None
    for index in indices:
        answer = short_add(coefficients, answer, negate(points[index]))
    return answer


def replay() -> dict[str, object]:
    roots = split_roots(S)
    if roots != (Q(0), Q(1), Q(137, 7), Q(79, 7), Q(81, 7), Q(128, 7)):
        raise AssertionError("the normalized seed roots changed")
    construction = SixRootMestreConstruction(roots)
    coordinates = component_coordinates(S)
    coefficients = construction.primitive_jacobian_coefficients(T)
    if coefficients != EXPECTED_COEFFICIENTS:
        raise AssertionError("the specialization short model changed")
    quartic = construction.primitive_quartic_coefficients(T)
    visible = tuple(
        quartic_point_to_jacobian(construction, T, point)
        for point in primitive_visible_points(construction, T)
    )
    affine = []
    for intercept, slope in ((coordinates[4], coordinates[5]), (coordinates[6], coordinates[7])):
        x_value = intercept + slope * T
        y_value = rational_square_root(quartic_value(quartic, x_value))
        if y_value in (None, 0):
            raise AssertionError("an affine component point lost its ordinate")
        affine.append(quartic_point_to_jacobian(construction, T, (x_value, y_value)))

    # These are zero-based positions in primitive_visible_points.  They are
    # retained as positions because this is an exact audit of one fixed fibre.
    relation_indices = ((0, 3, 4, 6, 9), (0, 3, 4))
    if tuple(
        negative_sum(coefficients, visible, indices) for indices in relation_indices
    ) != tuple(affine):
        raise AssertionError("the exact visible-subgroup relations changed")

    visible_certificate = mod3_independence_certificate(
        coefficients, visible, prime_bound=PRIME_BOUND
    )
    augmented_certificate = mod3_independence_certificate(
        coefficients, (*visible, *affine), prime_bound=PRIME_BOUND
    )
    if (
        visible_certificate["combined_exact_rank_over_F3"] != 9
        or augmented_certificate["combined_exact_rank_over_F3"] != 9
        or visible_certificate["independent_subset_indices_one_based"] != EXPECTED_PIVOTS
        or augmented_certificate["independent_subset_indices_one_based"] != EXPECTED_PIVOTS
    ):
        raise AssertionError("the finite-reduction non-promotion certificate changed")
    return {
        "status": "exact second-component seed relation and non-promotion audit completed",
        "component_parameter": "s=-357/47",
        "mestre_parameter": "T=1",
        "normalized_root_seed": [0, 7, 79, 81, 128, 137],
        "short_weierstrass_coefficients": [rational_text(value) for value in coefficients],
        "exact_affine_relations": [
            "P1=-V1-V4-V5-V7-V10",
            "P2=-V1-V4-V5",
        ],
        "visible_mod3_rank": 9,
        "augmented_mod3_rank": 9,
        "finite_reduction_pivots": EXPECTED_PIVOTS,
        "conclusion": (
            "at the normalized seed specialization both affine points are exactly in the "
            "visible subgroup; this fibre cannot witness a new Mordell-Weil direction"
        ),
        "not_established": [
            "a generic relation on the second rational component",
            "saturation of the visible rank-nine sublattice",
            "pair intersection numbers, a Shioda Gram matrix, or a generic rank statement",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(replay(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
