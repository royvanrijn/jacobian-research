#!/usr/bin/env python3
"""Exact conic-rational two-section component through (0,7,79,81,128,137).

Let ``u=c1+425/7`` and put

    u = -12*(47*s+357)/(s^2-49).

The formulas below give a rational component of the recursive two-section
incidence on the rational conic that splits its six roots.  The checker never
expands a universal residual: it evaluates the compact recursion at enough
rational s-values to prove the cleared rational-function identities.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from mestre_root_tuples import SixRootMestreConstruction
from probe_mestre_two_section_local_continuation import Field, residuals
from search_mestre_root_tuple_scale import quartic_value


Q = Fraction
SAMPLE_VALUES = tuple(value for value in range(1923) if value not in (-7, 7))
ROOT_PARAMETER_SAMPLES = tuple(value for value in range(-20, 21) if value not in (-7, 7))


def parameter_u(s_value: int | Fraction) -> Fraction:
    s_value = Q(s_value)
    if s_value in (-7, 7):
        raise ValueError("the component parameter has a pole at s=+/-7")
    return -12 * (47 * s_value + 357) / (s_value**2 - 49)


def component_coordinates(s_value: int | Fraction) -> tuple[Fraction, ...]:
    s_value = Q(s_value)
    u = parameter_u(s_value)
    denominator = 65 * s_value**2 + 658 * s_value + 1813
    if denominator == 0:
        raise ValueError("the rational section chart has a pole")
    c1 = -Q(425, 7) + u
    c2 = Q(66335, 49) - Q(929, 21) * u + Q(13, 36) * u**2
    c3 = -Q(4501495, 343) + Q(93718, 147) * u - Q(2599, 252) * u**2 + Q(1, 18) * u**3
    c4 = Q(112212864, 2401) - Q(1029264, 343) * u + Q(63671, 882) * u**2 - Q(583, 756) * u**3 + Q(1, 324) * u**4
    x01 = (137 * s_value**2 + 1316 * s_value + 3283) / denominator
    x11 = -(47 * s_value**2 + 714 * s_value + 2303) / denominator
    x02 = (
        (4 * s_value + 21)
        * (9 * s_value + 35)
        * (137 * s_value**2 + 1316 * s_value + 3283)
        / (7 * (s_value - 7) * (s_value + 7) * denominator)
    )
    return c1, c2, c3, c4, x01, x11, x02, Q(0)


def leading_invariant(coordinates: tuple[Fraction, ...]) -> Fraction:
    c1, c2, c3, c4 = coordinates[:4]
    a1, a2, a3, a4 = c1 - 1, c2 - c1, c3 - c2, c4 - c3
    return 5 * a1**4 - 24 * a1**2 * a2 + 32 * a1 * a3 + 16 * a2**2 - 64 * a4


def split_roots(s_value: int | Fraction) -> tuple[Fraction, ...]:
    s_value = Q(s_value)
    u = parameter_u(s_value)
    w = 282 + s_value * u
    return (
        Q(0),
        Q(1),
        (411 - 7 * u) / 21,
        (474 - 7 * u) / 42,
        (1254 - 21 * u - w) / 84,
        (1254 - 21 * u + w) / 84,
    )


def monic_polynomial_from_roots(roots: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    coefficients = [Q(1)]
    for root in roots:
        updated = [Q(0)] * (len(coefficients) + 1)
        for degree, coefficient in enumerate(coefficients):
            updated[degree] -= root * coefficient
            updated[degree + 1] += coefficient
        coefficients = updated
    return tuple(coefficients)


def q_coefficients(coordinates: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    c1, c2, c3, c4 = coordinates[:4]
    return (Q(0), -c4, c4 - c3, c3 - c2, c2 - c1, c1 - 1, Q(1))


def verify_split_roots() -> None:
    """Check the degree-at-most-twelve cleared root-product identities."""

    for s_value in ROOT_PARAMETER_SAMPLES:
        roots = split_roots(s_value)
        if len(set(roots)) != 6:
            raise AssertionError(f"the split roots collided at s={s_value}")
        if monic_polynomial_from_roots(roots) != q_coefficients(component_coordinates(s_value)):
            raise AssertionError(f"the split-root identity failed at s={s_value}")


def evaluate_polynomial(coefficients: tuple[Fraction, ...], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def seed_finite_intersection() -> dict[str, object]:
    """Audit the selected signed sections at their common finite seed fibre."""

    s_value = -Q(357, 47)
    construction = SixRootMestreConstruction(split_roots(s_value))
    coordinates = component_coordinates(s_value)
    lines = ((coordinates[4], coordinates[5]), (coordinates[6], coordinates[7]))
    # Primitive quartic ordinates from the triangular recurrence using the
    # D-square root -576/7 at the seed.
    ordinates = (
        (Q(-186731), Q(287206437, 4225), Q(-17138926, 4225), Q(-691488, 4225)),
        (Q(0), Q(125398903, 4225), Q(0), Q(-343)),
    )
    for parameter in range(1, 8):
        quartic = construction.primitive_quartic_coefficients(Q(parameter))
        for (intercept, slope), ordinate in zip(lines, ordinates):
            x_value = intercept + slope * parameter
            y_value = evaluate_polynomial(ordinate, Q(parameter))
            if y_value**2 != quartic_value(quartic, x_value):
                raise AssertionError("a seed cubic ordinate missed the primitive quartic")
    parameter = Q(3973, 329)
    points = tuple(
        (
            intercept + slope * parameter,
            evaluate_polynomial(ordinate, parameter),
        )
        for (intercept, slope), ordinate in zip(lines, ordinates)
    )
    expected = (Q(4932, 455), -Q(107740485691272, 438652175))
    if points[0] != points[1] or points[0] != expected:
        raise AssertionError("the declared finite seed intersection changed")
    return {
        "base_parameter": "T=3973/329",
        "common_affine_quartic_point": ["4932/455", "-107740485691272/438652175"],
        "square_identities_checked_at_distinct_T_values": 7,
        "conclusion": "the selected signed sections meet at this finite seed fibre",
        "not_established": "the full section intersection number, infinity contributions, or Shioda corrections",
    }


def replay() -> dict[str, object]:
    for s_value in SAMPLE_VALUES:
        coordinates = component_coordinates(s_value)
        u = parameter_u(s_value)
        expected = Q(16, 441) * (7 * u - 432) ** 2
        if leading_invariant(coordinates) != expected:
            raise AssertionError(f"the leading-square identity failed at s={s_value}")
        values = residuals(coordinates, Field(tangent=False))
        if any(value.value for value in values):
            raise AssertionError(f"a recursive residual failed at s={s_value}")
    verify_split_roots()
    intersection = seed_finite_intersection()
    return {
        "status": "exact conic-rational two-section component identity verified",
        "parameter": "u=-12*(47*s+357)/(s^2-49), with u=c1+425/7",
        "open_parameter_locus": "s != +/-7 and 65*s^2+658*s+1813 != 0",
        "normalized_root_seed": [0, 7, 79, 81, 128, 137],
        "leading_invariant": "D=16*(7*u-432)^2/441",
        "admissible_exact_sample_count": len(SAMPLE_VALUES),
        "residual_degree_bound_after_common_denominator": 1920,
        "all_recursive_residuals_vanish": True,
        "split_six_root_parameterization": {
            "conic": "w^2=49*u^2-4284*u+79524, with w=282+s*u",
            "sample_count": len(ROOT_PARAMETER_SAMPLES),
            "cleared_coefficient_degree_bound": 12,
            "all_root_product_coefficients_match": True,
        },
        "seed_finite_intersection": intersection,
        "not_established": [
            "pair intersections, a Shioda Gram matrix, or saturation",
            "independence from a generic rank-13 subgroup",
            "generic rank at least 14",
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
