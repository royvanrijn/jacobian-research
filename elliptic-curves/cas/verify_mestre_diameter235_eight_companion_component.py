#!/usr/bin/env python3
"""Exact rational two-section component through the diameter-235 seed.

The local split-root continuation at ``(0,17,136,161,207,235)`` recognizes a
rational nodal-quartic normalization parameter ``p``.  This verifier keeps
the Mestre square root and the triangular ordinate equations recursive: it
never constructs the enormous eliminated residual polynomials.

Every coordinate denominator divides ``D=B^4*K``, where

    B=(p-66)(p+54)(p^2+18p+456),  K=3p^2+4p+1068.

Since the established recursive residual degree bound is 60, clearing
``D^60`` has degree at most 1080.  The 1,099 admissible rational evaluations
therefore certify the seven rational-function residual identities exactly.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from probe_mestre_diameter235_eight_companion import ROOT_PLANE_QUARTIC, quotient_ranks
from probe_mestre_two_section_local_continuation import Field, residuals


Q = Fraction
ROOTS_AT_SEED = (Q(0), Q(1), Q(8), Q(207, 17), Q(161, 17), Q(235, 17))
SEED_PARAMETER = Q(-294)
SEED_MODULI = (
    -Q(739, 17), Q(201815, 289), -Q(24125597, 4913), Q(62654760, 4913)
)
SEED_SECTIONS = ((-Q(68, 27), -Q(37, 27)), (Q(3067, 459), -Q(5, 27)))
RESIDUAL_TOTAL_DEGREE_BOUND = 60
COMMON_DENOMINATOR_DEGREE = 18
CLEARED_RESIDUAL_DEGREE_BOUND = (
    RESIDUAL_TOTAL_DEGREE_BOUND * COMMON_DENOMINATOR_DEGREE
)
SAMPLE_VALUES = tuple(Q(value) for value in range(-550, 551) if value not in {-54, 66})


def poly(value: Q, coefficients: tuple[int, ...]) -> Q:
    """Evaluate an ascending integer coefficient polynomial."""

    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def denominator_b(value: Q) -> Q:
    return (value - 66) * (value + 54) * (value * value + 18 * value + 456)


def denominator_k(value: Q) -> Q:
    return 3 * value * value + 4 * value + 1068


def component_coordinates(value: Q) -> tuple[Q, ...]:
    """The four moduli and the two affine section abscissae over Q(p)."""

    value = Q(value)
    b = denominator_b(value)
    k = denominator_k(value)
    c1 = -poly(value, (-7245936, -246096, 4704, -876, 79)) / (2 * b)
    c2 = 5 * poly(
        value,
        (
            14277841000704, 888999349248, -30000699648, 2783206656,
            -278543520, -14802624, 293232, -40272, 1849,
        ),
    ) / (16 * b**2)
    c3 = -poly(
        value,
        (
            -32402742765539106816, -2677382748374243328,
            191800563007140864, -2244816893422080, 657002839322880,
            128514275004672, -3626525084544, 276573028032,
            -9211864320, -761052280, 16418244, -1903778, 59329,
        ),
    ) / (16 * b**3)
    c4 = (
        25 * (value - 26) * (value - 6) * (value + 6) * (value + 14)
        * (value * value - 12 * value + 276) * k
        * (7 * value * value - 204 * value - 2628)
        * (29 * value**3 + 378 * value**2 + 3132 * value + 177336)
        * (37 * value**3 - 126 * value**2 + 8316 * value - 269352)
        / (64 * b**4)
    )
    x01 = -(
        29 * value**3 + 378 * value**2 + 3132 * value + 177336
    ) / (3 * (value - 66) * k)
    x11 = -(13 * value**2 + 204 * value + 1908) / (3 * k)
    x02 = poly(
        value,
        (-2491329312, -167485968, -6454080, -770760, 9870, -1053, 53),
    ) / (3 * b * k)
    x12 = -2 * (value - 6) * (value + 54) / (3 * k)
    return (c1, c2, c3, c4, x01, x11, x02, x12)


def split_roots(value: Q) -> tuple[Q, ...]:
    """The six rational roots of the normalized degree-six Mestre polynomial."""

    value = Q(value)
    b = denominator_b(value)
    k = denominator_k(value)
    r3 = 8 - (
        (value + 294) * (3 * value**3 - 314 * value**2 - 7356 * value - 161208)
        / (4 * b)
    )
    r6 = Q(235, 17) - (
        45 * (value + 294)
        * (value**3 - 118 * value**2 - 2292 * value - 57416)
        / (34 * b)
    )
    r4 = k * (7 * value**2 - 204 * value - 2628) / (2 * b)
    r5 = (value + 14) * (37 * value**3 - 126 * value**2 + 8316 * value - 269352) / (4 * b)
    return (Q(0), Q(1), r3, r4, r5, r6)


def monic_polynomial_from_roots(roots: tuple[Q, ...]) -> tuple[Q, ...]:
    coefficients = [Q(1)]
    for root in roots:
        updated = [Q(0)] * (len(coefficients) + 1)
        for degree, coefficient in enumerate(coefficients):
            updated[degree] -= root * coefficient
            updated[degree + 1] += coefficient
        coefficients = updated
    return tuple(coefficients)


def q_coefficients(coordinates: tuple[Q, ...]) -> tuple[Q, ...]:
    c1, c2, c3, c4 = coordinates[:4]
    return (Q(0), -c4, c4 - c3, c3 - c2, c2 - c1, c1 - 1, Q(1))


def leading_invariant(coordinates: tuple[Q, ...]) -> Q:
    c1, c2, c3, c4 = coordinates[:4]
    a1, a2, a3, a4 = c1 - 1, c2 - c1, c3 - c2, c4 - c3
    return 5 * a1**4 - 24 * a1**2 * a2 + 32 * a1 * a3 + 16 * a2**2 - 64 * a4


def leading_square_root(value: Q) -> Q:
    b = denominator_b(value)
    return (
        15 * (value - 26) * (value - 6) * (value + 6) * (value + 14)
        * denominator_k(value)
        / ((value - 66) * (value + 54) * (value * value + 18 * value + 456) ** 2)
    )


def root_plane_quartic(u: Q, v: Q) -> Q:
    return sum(Q(coefficient) * u**u_degree * v**v_degree for u_degree, v_degree, coefficient in ROOT_PLANE_QUARTIC)


def verify_component() -> dict[str, object]:
    if len(SAMPLE_VALUES) <= CLEARED_RESIDUAL_DEGREE_BOUND:
        raise AssertionError("too few exact samples for the degree certificate")
    if component_coordinates(SEED_PARAMETER) != (*SEED_MODULI, *SEED_SECTIONS[0], *SEED_SECTIONS[1]):
        raise AssertionError("the rational component missed its local seed")
    if split_roots(SEED_PARAMETER) != ROOTS_AT_SEED:
        raise AssertionError("the split-root parameterization missed its local seed")
    residuals_checked = 0
    split_open_samples = 0
    for value in SAMPLE_VALUES:
        coordinates = component_coordinates(value)
        roots = split_roots(value)
        split_open_samples += len(set(roots)) == 6
        if monic_polynomial_from_roots(roots) != q_coefficients(coordinates):
            raise AssertionError(f"the split-root identity failed at p={value}")
        if leading_invariant(coordinates) != leading_square_root(value) ** 2:
            raise AssertionError(f"the leading square identity failed at p={value}")
        u = roots[5] - Q(235, 17)
        v = roots[2] - 8
        if root_plane_quartic(u, v):
            raise AssertionError(f"the nodal quartic normalization failed at p={value}")
        values = residuals(coordinates, Field(tangent=False))
        if any(item.value for item in values):
            raise AssertionError(f"a recursive residual failed at p={value}")
        residuals_checked += len(values)
    quotient = quotient_ranks()
    if [
        (record["T"], record["visible_mod3_rank"], record["visible_plus_eight_companions_mod3_rank"])
        for record in quotient
    ] != [("1", 10, 10), ("2", 10, 11), ("3", 10, 11), ("-1", 10, 10)]:
        raise AssertionError("the seed finite-reduction escape profile changed")
    return {
        "status": "exact rational diameter-235 two-section component verified",
        "parameter": "p, with local seed p=-294",
        "normalized_root_seed": ["0", "1", "8", "161/17", "207/17", "235/17"],
        "selected_sections": [["-68/27", "-37/27"], ["3067/459", "-5/27"]],
        "split_root_parameterization": "all six roots lie in Q(p)",
        "leading_invariant": "D=(15*(p-26)*(p-6)*(p+6)*(p+14)*(3*p^2+4*p+1068)/((p-66)*(p+54)*(p^2+18*p+456)^2))^2",
        "root_plane_quartic": {
            "variables": "U=r6-235/17, V=r3-8",
            "terms": [list(term) for term in ROOT_PLANE_QUARTIC],
            "normalization": "a rational nodal-quartic parameter p",
        },
        "expanded_two_section_residual_materialized": False,
        "residual_identity_certificate": {
            "common_coordinate_denominator": "((p-66)*(p+54)*(p^2+18*p+456))^4*(3*p^2+4*p+1068)",
            "common_coordinate_denominator_degree": COMMON_DENOMINATOR_DEGREE,
            "residual_total_degree_bound": RESIDUAL_TOTAL_DEGREE_BOUND,
            "cleared_numerator_degree_bound": CLEARED_RESIDUAL_DEGREE_BOUND,
            "admissible_exact_sample_count": len(SAMPLE_VALUES),
            "split_six_root_open_samples": split_open_samples,
            "recursive_residual_evaluations": residuals_checked,
            "all_recursive_residuals_vanish": True,
        },
        "finite_reduction_visible_subgroup_escape_at_seed": {
            "p": "-294",
            "records": quotient,
            "conclusion": "at T=2 and T=3, each selected affine companion separately raises the exact mod-3 finite-reduction rank above the visible subgroup",
            "scope_limit": "this proves a specialization-level visible-subgroup escape, not generic independence, saturation, or a rank-14 assertion",
        },
        "not_established": [
            "generic Mordell-Weil independence or rank at least 14",
            "pair intersections away from the seed, a Shioda Gram matrix, or saturation",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(verify_component(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
