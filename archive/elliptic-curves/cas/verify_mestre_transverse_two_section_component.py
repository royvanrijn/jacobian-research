#!/usr/bin/env python3
"""Exact rational two-section component through roots (0,1,7,8,9,11).

Put z=c1+35 and d=30-z.  This script verifies, without expanding any
universal residual, that

 c1=-35+z,
 c2=455-77z/3+13z^2/36,
 c3=-2605+652z/3-217z^2/36+z^3/18,
 c4=5544-608z+449z^2/18-49z^3/108+z^4/324,

 x1=(z^2-66z+1098)/(3d),  slope1=(42-z)/d,
 x2=(z^2-69z+1188)/(6d),  slope2=0

satisfy M and both triples (E2,E1,E0).  The leading invariant is
16(z-36)^2/9, so the triangular recursion supplies rational cubic ordinates
on z != 30,36.  The verification checks 301 admissible z-values.  After
clearing d^60, each residual numerator has degree at most 300, from the
recorded total-degree bounds 40, 50, 60; therefore the check is an exact
identity certificate, not a bounded search.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from mestre_root_tuples import SixRootMestreConstruction
from probe_mestre_two_section_local_continuation import Field, residuals
from search_mestre_root_tuple_scale import quartic_value


Q = Fraction
SAMPLE_VALUES = tuple(value for value in range(303) if value not in (30, 36))
ROOT_PARAMETER_SAMPLES = tuple(
    value for value in range(-20, 21) if value not in (-1, 0, 1)
)


def component_coordinates(z: int | Fraction) -> tuple[Fraction, ...]:
    z = Q(z)
    denominator = 30 - z
    if denominator == 0:
        raise ValueError("z=30 is the declared pole")
    return (
        -35 + z,
        455 - Q(77, 3) * z + Q(13, 36) * z**2,
        -2605 + Q(652, 3) * z - Q(217, 36) * z**2 + Q(1, 18) * z**3,
        5544 - 608 * z + Q(449, 18) * z**2 - Q(49, 108) * z**3 + Q(1, 324) * z**4,
        (z**2 - 66 * z + 1098) / (3 * denominator),
        (42 - z) / denominator,
        (z**2 - 69 * z + 1188) / (6 * denominator),
        Q(0),
    )


def leading_invariant(coordinates: tuple[Fraction, ...]) -> Fraction:
    c1, c2, c3, c4 = coordinates[:4]
    a1, a2, a3, a4 = c1 - 1, c2 - c1, c3 - c2, c4 - c3
    return 5 * a1**4 - 24 * a1**2 * a2 + 32 * a1 * a3 + 16 * a2**2 - 64 * a4


def split_roots(parameter: int | Fraction) -> tuple[Fraction, ...]:
    """Return the six rational roots after parametrizing the residual conic."""

    parameter = Q(parameter)
    if parameter in (-1, 1):
        raise ValueError("the conic parametrization has a pole at r=+/-1")
    z = 12 * (parameter + 3) / (1 - parameter**2)
    w = 6 + parameter * z
    return (
        Q(0),
        Q(1),
        (33 - z) / 3,
        (42 - z) / 6,
        (102 - 3 * z - w) / 12,
        (102 - 3 * z + w) / 12,
    )


def monic_polynomial_from_roots(roots: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    """Return the ascending coefficients of the monic polynomial with roots."""

    coefficients = [Q(1)]
    for root in roots:
        updated = [Q(0)] * (len(coefficients) + 1)
        for degree, coefficient in enumerate(coefficients):
            updated[degree] -= root * coefficient
            updated[degree + 1] += coefficient
        coefficients = updated
    return tuple(coefficients)


def q_coefficients(coordinates: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    """Return q=X(X-1)(X^4+c1*X^3+c2*X^2+c3*X+c4), ascending."""

    c1, c2, c3, c4 = coordinates[:4]
    return (Q(0), -c4, c4 - c3, c3 - c2, c2 - c1, c1 - 1, Q(1))


def verify_split_root_parameterization() -> None:
    """Certify q equals the product of six displayed rational roots over Q(r).

    Every cleared coefficient has degree at most twelve in r: the roots have
    denominator 1-r^2 and degree-at-most-two numerators, while q has the
    displayed degree-four z-coefficients.  The 38 samples below consequently
    prove the rational-function identity, not merely a finite sample fact.
    """

    for parameter in ROOT_PARAMETER_SAMPLES:
        z = Q(12 * (parameter + 3), 1 - parameter**2)
        if z == 30:
            raise AssertionError("a declared component pole entered the samples")
        roots = split_roots(parameter)
        if len(set(roots)) != 6:
            raise AssertionError(f"the split-root sample collided at r={parameter}")
        if monic_polynomial_from_roots(roots) != q_coefficients(component_coordinates(z)):
            raise AssertionError(f"the split-root identity failed at r={parameter}")


def evaluate_polynomial(coefficients: tuple[Fraction, ...], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def seed_finite_intersection() -> dict[str, object]:
    """Audit the selected signed sections at their common finite seed fibre."""

    construction = SixRootMestreConstruction(tuple(map(Q, (0, 1, 7, 8, 9, 11))))
    lines = ((Q(61, 5), Q(7, 5)), (Q(33, 5), Q(0)))
    # These are the primitive quartic ordinates from the triangular recurrence
    # using the D-square root +48 at z=0.
    ordinates = (
        (Q(-196), Q(-2811, 25), Q(-434, 25), Q(-24, 25)),
        (Q(0), Q(-634, 25), Q(0), Q(1)),
    )
    for parameter in range(1, 8):
        quartic = construction.primitive_quartic_coefficients(Q(parameter))
        for (intercept, slope), ordinate in zip(lines, ordinates):
            x_value = intercept + slope * parameter
            y_value = evaluate_polynomial(ordinate, Q(parameter))
            if y_value**2 != quartic_value(quartic, x_value):
                raise AssertionError("a seed cubic ordinate missed the primitive quartic")
    # Both sides have degree at most six in T, so the seven checks above prove
    # their square identities.  The abscissae agree only at T=-4.
    parameter = Q(-4)
    points = tuple(
        (
            intercept + slope * parameter,
            evaluate_polynomial(ordinate, parameter),
        )
        for (intercept, slope), ordinate in zip(lines, ordinates)
    )
    if points[0] != points[1] or points[0] != (Q(33, 5), Q(936, 25)):
        raise AssertionError("the declared finite intersection changed")
    return {
        "base_parameter": "T=-4",
        "common_affine_quartic_point": ["33/5", "936/25"],
        "square_identities_checked_at_distinct_T_values": 7,
        "conclusion": "the selected signed sections meet at this finite seed fibre",
        "not_established": "the full section intersection number, infinity contributions, or Shioda corrections",
    }


def replay() -> dict[str, object]:
    for z in SAMPLE_VALUES:
        coordinates = component_coordinates(z)
        if leading_invariant(coordinates) != Q(16, 9) * (z - 36) ** 2:
            raise AssertionError(f"leading square identity failed at z={z}")
        values = residuals(coordinates, Field(tangent=False))
        if any(value.value != 0 for value in values):
            raise AssertionError(f"a recursive residual failed at z={z}")
    verify_split_root_parameterization()
    intersection = seed_finite_intersection()
    return {
        "status": "exact rational two-section component identity verified",
        "parameter": "z=c1+35",
        "open_parameter_locus": "z != 30,36",
        "normalized_root_seed": [0, 1, 7, 8, 9, 11],
        "leading_invariant": "D=16*(z-36)^2/9",
        "residual_degree_bounds_after_substitution": {
            "M": 20,
            "E2": 200,
            "E1": 250,
            "E0": 300,
        },
        "admissible_exact_sample_count": len(SAMPLE_VALUES),
        "all_recursive_residuals_vanish": True,
        "split_six_root_parameterization": {
            "conic": "w^2=z^2-36*z+36, with z=12*(r+3)/(1-r^2), w=6+r*z",
            "sample_count": len(ROOT_PARAMETER_SAMPLES),
            "cleared_coefficient_degree_bound": 12,
            "all_root_product_coefficients_match": True,
        },
        "seed_finite_intersection": intersection,
        "triangular_ordinate_conclusion": (
            "both cubic ordinates are rational on z != 30,36"
        ),
        "not_established": [
            "section intersections at infinity or a Shioda Gram matrix",
            "saturation or independence from the existing generic rank-13 subgroup",
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
