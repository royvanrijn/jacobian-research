#!/usr/bin/env python3
"""Finite intersection and quotient audit for the Fermigier two-section curve.

At the local-continuation seed ``u=-3`` the two signed cubic ordinates meet at
one exact finite fibre.  The cubic ordinates are reconstructed solely from
the degree-six specialized quartics; this is a small univariate calculation,
not an expanded universal residual computation.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from icarm_curve245_mestre import fermigier_roots
from mestre_root_tuples import SixRootMestreConstruction
from nagao_1994 import (
    primitive_quartic_coefficients,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    quartic_value,
    short_jacobian_coefficients,
)
from probe_mestre_fermigier_two_section_local_continuation import (
    BASE_U,
    normalized_data,
    reconstructed_second_line,
)
from screen_mestre_fermigier_two_section_escape import square_root
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate


Q = Fraction
T_NODES = tuple(Q(value) for value in range(1, 8))
QUOTIENT_T = Q(1)
PRIME_BOUND = 499


def add(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    answer = [Q(0)] * max(len(left), len(right))
    for index, value in enumerate(left):
        answer[index] += value
    for index, value in enumerate(right):
        answer[index] += value
    return answer


def multiply(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    answer = [Q(0)] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] += left_value * right_value
    return answer


def interpolate(nodes: tuple[Fraction, ...], values: list[Fraction]) -> list[Fraction]:
    """Return the degree-at-most-six polynomial through seven exact values."""

    answer = [Q(0)]
    for index, (node, value) in enumerate(zip(nodes, values)):
        basis, denominator = [Q(1)], Q(1)
        for other_index, other in enumerate(nodes):
            if other_index == index:
                continue
            basis = multiply(basis, [-other, Q(1)])
            denominator *= node - other
        answer = add(answer, [value * coefficient / denominator for coefficient in basis])
    return answer


def evaluate(coefficients: list[Fraction], value: Fraction) -> Fraction:
    return sum((coefficient * value**degree for degree, coefficient in enumerate(coefficients)), Q(0))


def cubic_square_root(coefficients: list[Fraction]) -> list[Fraction]:
    """Recover the selected cubic square root of a degree-six polynomial."""

    if len(coefficients) != 7:
        raise ValueError("the specialized quartic must have degree six")
    cubic = [Q(0)] * 4
    cubic[3] = square_root(coefficients[6])
    cubic[2] = coefficients[5] / (2 * cubic[3])
    cubic[1] = (coefficients[4] - cubic[2] ** 2) / (2 * cubic[3])
    cubic[0] = (coefficients[3] - 2 * cubic[2] * cubic[1]) / (2 * cubic[3])
    if multiply(cubic, cubic) != coefficients:
        raise AssertionError("the triangular cubic square root did not close")
    return cubic


def seed_data() -> tuple[SixRootMestreConstruction, tuple[tuple[Fraction, Fraction], ...]]:
    u = BASE_U
    v, second_intercept, second_slope = reconstructed_second_line(u)
    first_intercept, first_slope = normalized_data(u, v)[4:6]
    source_roots = fermigier_roots(u, v)
    roots = tuple(
        (root - source_roots[0]) / (source_roots[1] - source_roots[0])
        for root in source_roots
    )
    return SixRootMestreConstruction(roots), (
        (Q(first_intercept), Q(first_slope)),
        (Q(second_intercept), Q(second_slope)),
    )


def replay() -> dict[str, object]:
    construction, lines = seed_data()
    cubic_ordinates = []
    for intercept, slope in lines:
        values = [
            quartic_value(
                primitive_quartic_coefficients(construction, parameter),
                intercept + slope * parameter,
            )
            for parameter in T_NODES
        ]
        cubic_ordinates.append(cubic_square_root(interpolate(T_NODES, values)))
    meeting_parameter = (lines[1][0] - lines[0][0]) / (lines[0][1] - lines[1][1])
    meeting_x = lines[0][0] + lines[0][1] * meeting_parameter
    meeting_y = evaluate(cubic_ordinates[0], meeting_parameter)
    if evaluate(cubic_ordinates[1], meeting_parameter) != meeting_y:
        raise AssertionError("the selected signed sections missed their finite intersection")

    quartic = primitive_quartic_coefficients(construction, QUOTIENT_T)
    visible = tuple(
        quartic_point_to_short_jacobian(construction, QUOTIENT_T, point)
        for point in primitive_visible_points(construction, QUOTIENT_T)
    )
    affine = tuple(
        quartic_point_to_short_jacobian(
            construction,
            QUOTIENT_T,
            (
                intercept + slope * QUOTIENT_T,
                square_root(quartic_value(quartic, intercept + slope * QUOTIENT_T)),
            ),
        )
        for intercept, slope in lines
    )
    coefficients = short_jacobian_coefficients(construction, QUOTIENT_T)
    baseline = mod3_independence_certificate(
        coefficients, visible, prime_bound=PRIME_BOUND
    )
    augmented = mod3_independence_certificate(
        coefficients, (*visible, *affine), prime_bound=PRIME_BOUND
    )
    if baseline["combined_exact_rank_over_F3"] != augmented["combined_exact_rank_over_F3"]:
        raise AssertionError("the seed unexpectedly escaped the visible finite quotient")
    return {
        "status": "Fermigier two-section seed intersection and finite-quotient audit completed",
        "component_parameter": "u=-3, v=-8/3",
        "normalized_roots": [str(value) for value in construction.roots],
        "signed_affine_lines": [[str(value) for value in line] for line in lines],
        "cubic_ordinates": [[str(value) for value in cubic] for cubic in cubic_ordinates],
        "cubic_square_identity_degree": 6,
        "square_identity_sample_count": len(T_NODES),
        "finite_intersection": {
            "T": str(meeting_parameter),
            "x": str(meeting_x),
            "y": str(meeting_y),
        },
        "finite_reduction": {
            "Mestre_parameter": str(QUOTIENT_T),
            "prime_bound": PRIME_BOUND,
            "visible_mod3_rank": baseline["combined_exact_rank_over_F3"],
            "augmented_mod3_rank": augmented["combined_exact_rank_over_F3"],
        },
        "conclusion": "the selected signed pair meets at the recorded finite fibre; this seed does not escape the visible mod-3 quotient",
        "not_established": [
            "the full pair intersection number including infinity contributions",
            "a Shioda Gram matrix or saturation",
            "independence from the existing generic rank-13 subgroup",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(replay(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
