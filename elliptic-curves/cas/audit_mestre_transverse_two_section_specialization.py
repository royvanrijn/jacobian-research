#!/usr/bin/env python3
"""Finite-reduction quotient audit on the rational transverse component.

At r=8 and T=3, this script maps the twelve visible Mestre points and the
two affine points from the exact transverse component to the Jacobian.  It
then compares their exact mod-3 finite-reduction spans.  The two additional
points do not enlarge that *particular quotient span*.  This is an explicit
non-promotion gate, not a Mordell--Weil relation, saturation calculation,
height pairing, or generic-rank result.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from mestre_root_tuples import SixRootMestreConstruction
from search_mestre_dsquare_four import rational_square_root
from search_mestre_root_tuple_scale import (
    point_digest,
    primitive_visible_points,
    quartic_point_to_jacobian,
    quartic_value,
)
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate
from verify_mestre_transverse_two_section_component import (
    component_coordinates,
    split_roots,
)


Q = Fraction
R = Q(8)
T = Q(3)
PRIME_BOUND = 499

EXPECTED_VISIBLE = {
    "rank": 9,
    "pivots": [1, 2, 3, 4, 5, 7, 8, 9, 10],
    "primes": [29, 31, 61, 71, 83, 101, 109, 131, 139],
    "point_sha256": "92dd79f91f31219967cda2f9b1a46695ca8099ed0234cd2019e505edda2442b5",
    "subset_sha256": "dae01b8c494f3b3b1af111fa580b1361ed20634a13ee20431dd89b1a7eee91e7",
}
EXPECTED_AUGMENTED = {
    "rank": 9,
    "pivots": [1, 2, 3, 4, 5, 7, 8, 9, 10],
    "primes": [29, 31, 61, 71, 83, 101, 109, 131, 139],
    "point_sha256": "1aa222440fa702810ff5fd3e055c398abca3b0b7a87a01d8bc6fbd96292cfd2d",
    "subset_sha256": "dae01b8c494f3b3b1af111fa580b1361ed20634a13ee20431dd89b1a7eee91e7",
}


def certificate_summary(certificate: dict[str, object]) -> dict[str, object]:
    return {
        "rank": certificate["combined_exact_rank_over_F3"],
        "pivots": certificate["independent_subset_indices_one_based"],
        "primes": certificate["certificate_primes"],
        "point_sha256": certificate["point_sha256"],
        "subset_sha256": certificate["independent_subset_sha256"],
    }


def rational_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def replay() -> dict[str, object]:
    roots = split_roots(R)
    construction = SixRootMestreConstruction(roots)
    z = 12 * (R + 3) / (1 - R**2)
    coordinates = component_coordinates(z)
    x01, x11, x02, x12 = coordinates[4:]
    quartic = construction.primitive_quartic_coefficients(T)
    visible_quartic = primitive_visible_points(construction, T)
    visible = tuple(
        quartic_point_to_jacobian(construction, T, point)
        for point in visible_quartic
    )
    affine = []
    for intercept, slope in ((x01, x11), (x02, x12)):
        x_value = intercept + slope * T
        y_value = rational_square_root(quartic_value(quartic, x_value))
        if y_value in (None, 0):
            raise AssertionError("a component section lost its nonzero ordinate")
        affine.append((x_value, y_value))
    extra = tuple(
        quartic_point_to_jacobian(construction, T, point) for point in affine
    )
    coefficients = construction.primitive_jacobian_coefficients(T)
    visible_certificate = mod3_independence_certificate(
        coefficients, visible, prime_bound=PRIME_BOUND
    )
    augmented_certificate = mod3_independence_certificate(
        coefficients, (*visible, *extra), prime_bound=PRIME_BOUND
    )
    visible_summary = certificate_summary(visible_certificate)
    augmented_summary = certificate_summary(augmented_certificate)
    if visible_summary != EXPECTED_VISIBLE:
        raise AssertionError("the visible finite-reduction certificate changed")
    if augmented_summary != EXPECTED_AUGMENTED:
        raise AssertionError("the augmented finite-reduction certificate changed")
    return {
        "status": "finite-reduction non-promotion audit completed",
        "component_parameter": "r=8",
        "mestre_parameter": "T=3",
        "six_roots": [rational_text(value) for value in roots],
        "visible_point_count": len(visible),
        "affine_quartic_points": [
            [rational_text(value) for value in point] for point in affine
        ],
        "visible_certificate": visible_summary,
        "augmented_certificate": augmented_summary,
        "quotient_conclusion": (
            "the two affine points add no new direction to this combined exact "
            "mod-3 finite-reduction span"
        ),
        "not_established": [
            "a Mordell-Weil relation or dependence of either affine point",
            "saturation, pair intersections, a Shioda Gram matrix, or a height pairing",
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
