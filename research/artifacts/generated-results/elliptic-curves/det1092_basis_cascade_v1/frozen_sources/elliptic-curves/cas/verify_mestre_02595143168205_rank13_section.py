#!/usr/bin/env python3
"""Replay the generic identity and specialization certificate for a new section.

The companion Singular verifier proves over Q(T) that

    x/25 = -95/23 + (37/23)*(T/25)

has a cubic ordinate on the Mestre quartic for roots
``(0,25,95,143,168,205)``.  Here ``T`` is in the original integral-root
coordinate.  After the split-infinity base change

    T = (39146-u^2)/(2u),

this script specializes at ``u=197`` and gives an exact finite-reduction
certificate.  Twelve of the visible-plus-infinity sections and this new
section are independent, so specialization of a hypothetical generic
relation proves a rank lower bound of 13 over Q(u).  This is a lower bound,
not an upper bound for the generic or specialized Mordell--Weil rank.

The generic identity is independently reconstructed with exact polynomial
arithmetic: seven rational specializations interpolate every primitive
quartic coefficient (whose degree is at most six), and the displayed cubic
ordinate is then squared coefficient by coefficient.  No point search or
local artifact is read.
"""

from __future__ import annotations

from fractions import Fraction
import json
from typing import Any

from search_mestre_dsquare_four import (
    FAMILIES,
    base_parameter,
    known_jacobian_points,
    rational_square_root,
)
from search_mestre_root_tuple_scale import (
    point_digest,
    quartic_point_to_jacobian,
    quartic_value,
)
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate


Q = Fraction
FAMILY_INDEX = 2
PARAMETER_U = Q(197)
PARAMETER_T = Q(337, 394)
PRIME_BOUND = 499

EXPECTED_KNOWN = {
    "point_sha256": "9620340124a630387ebad9585474c06863a5060f9139b4b8bd2560060100777a",
    "rank": 12,
    "pivots": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13],
    "primes": [37, 41, 61, 67, 79, 83, 101, 103, 137, 139, 149],
    "subset_sha256": "fcf8da83bf4a02c1827c75f6e0160960e4b11a7e86ae8df71787e7f17f52eef8",
}
EXPECTED_AUGMENTED = {
    "point_sha256": "56b2ec76169bd880fdb1c046665b277494dd0664acd4f1a9fc43153ba0f31dc0",
    "rank": 13,
    "pivots": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14],
    "primes": [37, 41, 61, 67, 79, 83, 101, 103, 137, 139, 149, 163],
    "subset_sha256": "70c996b76fb8eb3adb7d90a8daa6fcddd189cb334f0b4a27fd185ebfd8b4e27e",
}
EXPECTED_ALL_COMPANIONS = {
    "point_sha256": "ce499539bec4d42d77590fade016eb858bed32c55993abc578ef0955b5c6714f",
    "rank": 13,
    "pivots": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14],
    "primes": [37, 41, 61, 67, 79, 83, 101, 103, 137, 139, 149, 163],
    "subset_sha256": "70c996b76fb8eb3adb7d90a8daa6fcddd189cb334f0b4a27fd185ebfd8b4e27e",
}


def trim(poly: list[Fraction]) -> list[Fraction]:
    """Return an ascending coefficient list without trailing zeroes."""

    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def polynomial_add(*polynomials: list[Fraction]) -> list[Fraction]:
    result = [Q(0)] * max(map(len, polynomials))
    for polynomial in polynomials:
        for index, coefficient in enumerate(polynomial):
            result[index] += coefficient
    return trim(result)


def polynomial_scale(poly: list[Fraction], scalar: Fraction) -> list[Fraction]:
    return trim([Q(scalar) * coefficient for coefficient in poly])


def polynomial_multiply(
    left: list[Fraction], right: list[Fraction]
) -> list[Fraction]:
    result = [Q(0)] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] += left_value * right_value
    return trim(result)


def polynomial_power(poly: list[Fraction], exponent: int) -> list[Fraction]:
    result = [Q(1)]
    for _ in range(exponent):
        result = polynomial_multiply(result, poly)
    return result


def interpolate(
    nodes: list[Fraction], values: list[Fraction]
) -> list[Fraction]:
    """Interpolate the supplied distinct rational nodes over Q."""

    if len(nodes) != len(values) or len(set(nodes)) != len(nodes):
        raise ValueError("interpolation nodes and values are inconsistent")
    result = [Q(0)]
    for index, value in enumerate(values):
        basis = [Q(1)]
        denominator = Q(1)
        for other in range(len(values)):
            if other == index:
                continue
            basis = polynomial_multiply(basis, [-nodes[other], Q(1)])
            denominator *= nodes[index] - nodes[other]
        result = polynomial_add(
            result, polynomial_scale(basis, Q(value) / denominator)
        )
    return trim(result)


def verify_generic_companion_identity() -> dict[str, Any]:
    """Reconstruct the primitive quartic and check one identity in Q[T]."""

    construction = FAMILIES[FAMILY_INDEX].construction
    nodes = [Q(parameter_t) for parameter_t in range(1, 8)]
    sampled = [
        construction.primitive_quartic_coefficients(parameter_t)
        for parameter_t in nodes
    ]
    reconstructed = tuple(
        interpolate(nodes, [row[index] for row in sampled]) for index in range(5)
    )
    expected = (
        [Q(3_140_781_450_625), 0, Q(371_623_025), 0, Q(-36_418), 0, Q(1)],
        [Q(-152_865_039_600), 0, Q(-8_356_594), 0, Q(424)],
        [Q(2_523_380_759), 0, Q(82_536), 0, Q(-2)],
        [Q(-16_830_158), 0, Q(-424)],
        [Q(39_146), 0, Q(1)],
    )
    assert_equal(reconstructed, expected, "interpolated primitive quartic")

    x_polynomial = [Q(-2_375, 23), Q(37, 23)]
    y_polynomial = [
        Q(8_295_400),
        Q(-71_029_947, 529),
        Q(356_162, 529),
        Q(-840, 529),
    ]
    right_hand_side = polynomial_add(
        *(
            polynomial_multiply(coefficient, polynomial_power(x_polynomial, power))
            for power, coefficient in enumerate(reconstructed)
        )
    )
    assert_equal(
        polynomial_multiply(y_polynomial, y_polynomial),
        right_hand_side,
        "generic companion square identity",
    )
    return {
        "x": "(-2375+37*T)/23",
        "y": "8295400-(71029947/529)*T+(356162/529)*T^2-(840/529)*T^3",
        "identity_verified_over": "Q[T]",
        "quartic_reconstruction": "exact degree-at-most-6 interpolation at T=1,...,7",
    }


def assert_equal(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise AssertionError(f"{label} changed: {actual!r} != {expected!r}")


def rational_text(value: Fraction) -> str:
    value = Q(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def certificate_summary(certificate: dict[str, Any]) -> dict[str, Any]:
    return {
        "point_sha256": certificate["point_sha256"],
        "rank": certificate["combined_exact_rank_over_F3"],
        "pivots": certificate["independent_subset_indices_one_based"],
        "primes": certificate["certificate_primes"],
        "subset_sha256": certificate["independent_subset_sha256"],
    }


def companion_point(
    intercept_normalized: Fraction,
    slope: Fraction,
    *,
    quartic_coefficients: tuple[Fraction, ...],
) -> tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]:
    family = FAMILIES[FAMILY_INDEX]
    # Original X is 25 times normalized X, while the slope is unchanged.
    quartic_x = 25 * Q(intercept_normalized) + Q(slope) * PARAMETER_T
    quartic_y = rational_square_root(quartic_value(quartic_coefficients, quartic_x))
    if quartic_y is None or quartic_y == 0:
        raise AssertionError("the exact companion specialization lost its ordinate")
    jacobian_point = quartic_point_to_jacobian(
        family.construction,
        PARAMETER_T,
        (quartic_x, quartic_y),
    )
    return (quartic_x, quartic_y), jacobian_point


def replay() -> dict[str, Any]:
    family = FAMILIES[FAMILY_INDEX]
    assert_equal(family.roots, (0, 25, 95, 143, 168, 205), "family roots")
    assert_equal(base_parameter(family, PARAMETER_U), PARAMETER_T, "base parameter T")
    generic_identity = verify_generic_companion_identity()

    quartic_coefficients = family.construction.primitive_quartic_coefficients(
        PARAMETER_T
    )
    jacobian_coefficients = family.construction.primitive_jacobian_coefficients(
        PARAMETER_T
    )
    known_points = known_jacobian_points(family, PARAMETER_U)
    assert_equal(len(known_points), 13, "visible-plus-infinity point count")
    assert_equal(point_digest(known_points), EXPECTED_KNOWN["point_sha256"], "known digest")
    known_certificate = mod3_independence_certificate(
        jacobian_coefficients,
        known_points,
        prime_bound=PRIME_BOUND,
    )
    assert_equal(certificate_summary(known_certificate), EXPECTED_KNOWN, "known certificate")

    section_data = (
        (Q(-95, 23), Q(37, 23)),
        (Q(-95, 23), Q(-37, 23)),
        (Q(583, 115), Q(13, 23)),
        (Q(583, 115), Q(-13, 23)),
        (Q(3444, 575), Q(7, 23)),
        (Q(3444, 575), Q(-7, 23)),
    )
    quartic_points = []
    jacobian_points = []
    for intercept, slope in section_data:
        quartic_point, jacobian_point = companion_point(
            intercept,
            slope,
            quartic_coefficients=quartic_coefficients,
        )
        quartic_points.append(quartic_point)
        jacobian_points.append(jacobian_point)

    assert_equal(quartic_points[0][0], Q(-923281, 9062), "new-section quartic x")
    assert_equal(
        quartic_points[0][1],
        Q(33087435768778776, 4044402317),
        "new-section quartic y",
    )
    section_ordinate = (
        Q(8_295_400)
        - Q(71_029_947, 529) * PARAMETER_T
        + Q(356_162, 529) * PARAMETER_T**2
        - Q(840, 529) * PARAMETER_T**3
    )
    assert_equal(
        quartic_points[0][1],
        section_ordinate,
        "new-section specialized generic ordinate",
    )

    augmented_certificate = mod3_independence_certificate(
        jacobian_coefficients,
        known_points + (jacobian_points[0],),
        prime_bound=PRIME_BOUND,
    )
    assert_equal(
        certificate_summary(augmented_certificate),
        EXPECTED_AUGMENTED,
        "one-section certificate",
    )

    all_certificate = mod3_independence_certificate(
        jacobian_coefficients,
        known_points + tuple(jacobian_points),
        prime_bound=PRIME_BOUND,
    )
    assert_equal(
        certificate_summary(all_certificate),
        EXPECTED_ALL_COMPANIONS,
        "all-companion certificate",
    )

    return {
        "status": "verified generic rank lower bound 13 after split-infinity base change",
        "family_roots": list(family.roots),
        "u": rational_text(PARAMETER_U),
        "T": rational_text(PARAMETER_T),
        "new_section_original_x": "(-2375+37*T)/23",
        "generic_companion_identity": generic_identity,
        "new_section_primitive_y": (
            "8295400-(71029947/529)*T+(356162/529)*T^2-(840/529)*T^3"
        ),
        "specialized_new_quartic_point": [
            rational_text(quartic_points[0][0]),
            rational_text(quartic_points[0][1]),
        ],
        "visible_plus_infinity_point_count": len(known_points),
        "visible_plus_infinity_rank_lower_bound": known_certificate[
            "combined_exact_rank_over_F3"
        ],
        "one_companion_rank_lower_bound": augmented_certificate[
            "combined_exact_rank_over_F3"
        ],
        "one_companion_pivots": augmented_certificate[
            "independent_subset_indices_one_based"
        ],
        "all_six_companions_rank_lower_bound": all_certificate[
            "combined_exact_rank_over_F3"
        ],
        "certificate_primes": augmented_certificate["certificate_primes"],
        "generic_rank_lower_bound_after_base_change": 13,
        "claim_limit": (
            "the six companions jointly certify only one additional direction at "
            "u=197; neither generic nor specialized rank is bounded above"
        ),
    }


def main() -> None:
    print(json.dumps(replay(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
