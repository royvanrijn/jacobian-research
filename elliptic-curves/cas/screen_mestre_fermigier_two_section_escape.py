#!/usr/bin/env python3
"""Bounded finite-quotient escape audit at the eight-companion Fermigier seed.

The normalized roots ``(0,8,58,77,85,102)`` are Fermigier's source roots at
``(u,v)=(-3,-8/3)`` after ``alpha_1 -> 0, alpha_2 -> 1``.  Fermigier's extra
line becomes the first displayed affine line.  This screen asks whether the
second displayed line gains a dimension in exact mod-2 or mod-3 finite
quotients after specializations ``T=1,...,40``.  Non-escape is only bounded
negative evidence, never a dependence or rank-upper-bound proof.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import isqrt
from pathlib import Path

from icarm_curve245_mestre import fermigier_extra_line, fermigier_roots
from mestre_root_tuples import SixRootMestreConstruction
from nagao_1994 import (
    primitive_visible_points,
    primitive_quartic_coefficients,
    quartic_point_to_short_jacobian,
    quartic_value,
    short_jacobian_coefficients,
)
from search_mestre_root_tuple_scale import finite_reduction_attempt
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate


Q = Fraction
ROOTS = (0, 8, 58, 77, 85, 102)
PARAMETERS = (Q(-3), -Q(8, 3))
EXTRA_LINE = (Q(-58, 11) * 8, -Q(17, 11))
SECOND_LINE = (Q(247, 44) * 8, -Q(3, 11))
DEFAULT_BOUND = 40
DEFAULT_PRIME_BOUND = 151


def square_root(value: Fraction) -> Fraction:
    value = Q(value)
    if value < 0 or isqrt(value.numerator) ** 2 != value.numerator or isqrt(value.denominator) ** 2 != value.denominator:
        raise AssertionError("a declared affine section did not specialize to a rational quartic point")
    return Q(isqrt(value.numerator), isqrt(value.denominator))


def verify_fermigier_alignment() -> None:
    source = fermigier_roots(*PARAMETERS)
    normalized = tuple(sorted((root - source[0]) / (source[1] - source[0]) for root in source))
    target = tuple(sorted(Q(root, ROOTS[1]) for root in ROOTS))
    if normalized != target:
        raise AssertionError("Fermigier root normalization changed")
    intercept, slope = fermigier_extra_line(*PARAMETERS)
    mapped = ((intercept - source[0]) / (source[1] - source[0]) * ROOTS[1], slope)
    if mapped != EXTRA_LINE:
        raise AssertionError("Fermigier extra line no longer matches the first affine line")


def specialized_points(construction: SixRootMestreConstruction, parameter: Fraction) -> tuple[tuple[Fraction, Fraction], ...]:
    points = list(primitive_visible_points(construction, parameter))
    coefficients = primitive_quartic_coefficients(construction, parameter)
    for intercept, slope in (EXTRA_LINE, SECOND_LINE):
        x_value = intercept + slope * parameter
        points.append((x_value, square_root(quartic_value(coefficients, x_value))))
    return tuple(points)


def replay(
    parameter_bound: int = DEFAULT_BOUND,
    prime_bound: int = DEFAULT_PRIME_BOUND,
) -> dict[str, object]:
    if parameter_bound < 1:
        raise ValueError("parameter bound must be positive")
    verify_fermigier_alignment()
    construction = SixRootMestreConstruction(tuple(map(Q, ROOTS)))
    records = []
    for integer_parameter in range(1, parameter_bound + 1):
        parameter = Q(integer_parameter)
        quartic_points = specialized_points(construction, parameter)
        jacobian_points = tuple(
            quartic_point_to_short_jacobian(construction, parameter, point)
            for point in quartic_points
        )
        coefficients = short_jacobian_coefficients(construction, parameter)
        mod2_base = finite_reduction_attempt(coefficients, jacobian_points[:-1], prime_bound=prime_bound)
        mod2_augmented = finite_reduction_attempt(coefficients, jacobian_points, prime_bound=prime_bound)
        mod3_base = mod3_independence_certificate(coefficients, jacobian_points[:-1], prime_bound=prime_bound)
        mod3_augmented = mod3_independence_certificate(coefficients, jacobian_points, prime_bound=prime_bound)
        if mod2_augmented["combined_exact_rank_over_F2"] > mod2_base["combined_exact_rank_over_F2"]:
            raise AssertionError("the candidate escaped the mod-2 baseline; update the pinned audit")
        if mod3_augmented["combined_exact_rank_over_F3"] > mod3_base["combined_exact_rank_over_F3"]:
            raise AssertionError("the candidate escaped the mod-3 baseline; update the pinned audit")
        records.append(
            {
                "T": str(parameter),
                "mod2_baseline_rank": mod2_base["combined_exact_rank_over_F2"],
                "mod2_augmented_rank": mod2_augmented["combined_exact_rank_over_F2"],
                "mod3_baseline_rank": mod3_base["combined_exact_rank_over_F3"],
                "mod3_augmented_rank": mod3_augmented["combined_exact_rank_over_F3"],
            }
        )
    return {
        "status": "bounded exact finite-quotient Fermigier two-section escape audit completed",
        "roots": list(ROOTS),
        "Fermigier_parameters": [str(value) for value in PARAMETERS],
        "Fermigier_extra_line": [str(value) for value in EXTRA_LINE],
        "candidate_second_line": [str(value) for value in SECOND_LINE],
        "specialization_parameters": [record["T"] for record in records],
        "prime_bound": prime_bound,
        "records": records,
        "mod2_escape_count": 0,
        "mod3_escape_count": 0,
        "conclusion": "the second line has no finite-quotient escape in this bounded specialization screen",
        "not_established": [
            "a dependence relation for the second line",
            "saturation, a Shioda Gram matrix, or a Mordell--Weil rank upper bound",
            "generic rank at least fourteen",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameter-bound", type=int, default=DEFAULT_BOUND)
    parser.add_argument("--prime-bound", type=int, default=DEFAULT_PRIME_BOUND)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(replay(args.parameter_bound, args.prime_bound), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
