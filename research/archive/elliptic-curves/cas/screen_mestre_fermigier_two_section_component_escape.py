#!/usr/bin/env python3
"""Bounded quotient-escape screen along the exact Fermigier component.

This tests the second reconstructed section only after adjoining the twelve
visible points and Fermigier's first section.  It is deliberately a negative
bounded screen, not a relation or saturation proof.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import gcd
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
    normalized_data,
    reconstructed_second_line,
)
from screen_mestre_fermigier_two_section_escape import square_root
from search_mestre_root_tuple_scale import finite_reduction_attempt
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate


Q = Fraction
DEFAULT_HEIGHT = 5
DEFAULT_T_BOUND = 4
DEFAULT_PRIME_BOUND = 151


def rational_parameters(height: int) -> tuple[Fraction, ...]:
    values: list[Fraction] = []
    for denominator in range(1, height + 1):
        for numerator in range(-height, height + 1):
            if gcd(numerator, denominator) != 1:
                continue
            value = Q(numerator, denominator)
            if value not in values and value not in (0, -2):
                values.append(value)
    return tuple(values)


def replay(
    *, height: int = DEFAULT_HEIGHT, t_bound: int = DEFAULT_T_BOUND,
    prime_bound: int = DEFAULT_PRIME_BOUND,
) -> dict[str, object]:
    if height < 1 or t_bound < 1:
        raise ValueError("height and t bound must be positive")
    records = []
    for u in rational_parameters(height):
        try:
            v, second_intercept, second_slope = reconstructed_second_line(u)
            first_intercept, first_slope = normalized_data(u, v)[4:6]
            source = fermigier_roots(u, v)
            roots = tuple(
                (root - source[0]) / (source[1] - source[0])
                for root in source
            )
            if len(set(roots)) != 6:
                continue
            construction = SixRootMestreConstruction(roots)
            for parameter in map(Q, range(1, t_bound + 1)):
                quartic = primitive_quartic_coefficients(construction, parameter)
                points = list(primitive_visible_points(construction, parameter))
                for intercept, slope in (
                    (first_intercept, first_slope),
                    (second_intercept, second_slope),
                ):
                    x_value = intercept + slope * parameter
                    points.append((x_value, square_root(quartic_value(quartic, x_value))))
                jacobian_points = tuple(
                    quartic_point_to_short_jacobian(construction, parameter, point)
                    for point in points
                )
                coefficients = short_jacobian_coefficients(construction, parameter)
                ranks = [
                    mod3_independence_certificate(
                        coefficients, subset, prime_bound=prime_bound
                    )["combined_exact_rank_over_F3"]
                    for subset in (
                        jacobian_points[:12], jacobian_points[:13], jacobian_points,
                    )
                ]
                mod2_baseline = finite_reduction_attempt(
                    coefficients, jacobian_points[:13], prime_bound=prime_bound
                )["combined_exact_rank_over_F2"]
                mod2_augmented = finite_reduction_attempt(
                    coefficients, jacobian_points, prime_bound=prime_bound
                )["combined_exact_rank_over_F2"]
                records.append(
                    {
                        "u": str(u), "T": str(parameter),
                        "visible_mod3_rank": ranks[0],
                        "visible_plus_first_mod3_rank": ranks[1],
                        "visible_plus_two_mod3_rank": ranks[2],
                        "visible_plus_first_mod2_rank": mod2_baseline,
                        "visible_plus_two_mod2_rank": mod2_augmented,
                    }
                )
        except (ArithmeticError, ValueError, ZeroDivisionError):
            continue
    first_gain = sum(
        record["visible_plus_first_mod3_rank"] > record["visible_mod3_rank"]
        for record in records
    )
    second_escape = sum(
        record["visible_plus_two_mod3_rank"] > record["visible_plus_first_mod3_rank"]
        for record in records
    )
    mod2_second_escape = sum(
        record["visible_plus_two_mod2_rank"]
        > record["visible_plus_first_mod2_rank"]
        for record in records
    )
    if second_escape or mod2_second_escape:
        raise AssertionError("the second section escaped; replace this non-escape artifact")
    return {
        "status": "bounded Fermigier two-section component quotient screen completed",
        "u_height": height,
        "T_values": [str(value) for value in range(1, t_bound + 1)],
        "prime_bound": prime_bound,
        "record_count": len(records),
        "first_section_positive_gain_count": first_gain,
        "second_section_escape_count": second_escape,
        "second_section_mod2_escape_count": mod2_second_escape,
        "records": records,
        "conclusion": "the second reconstructed section has no mod-2 or mod-3 quotient escape after the visible-plus-first baseline in this bounded grid",
        "not_established": [
            "a dependence relation for the second section",
            "saturation, a Shioda Gram matrix, or generic rank",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=DEFAULT_HEIGHT)
    parser.add_argument("--t-bound", type=int, default=DEFAULT_T_BOUND)
    parser.add_argument("--prime-bound", type=int, default=DEFAULT_PRIME_BOUND)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(
        replay(height=args.height, t_bound=args.t_bound, prime_bound=args.prime_bound),
        indent=2, sort_keys=True,
    ) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
