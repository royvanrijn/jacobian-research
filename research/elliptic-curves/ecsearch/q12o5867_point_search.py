#!/usr/bin/env python3
"""Exact helpers for direct point search on q12o5867 specializations."""

from __future__ import annotations

from fractions import Fraction
from math import comb, isqrt, lcm
from pathlib import Path
import re
import sys
from typing import Any, Sequence


ELLIPTIC_ROOT = Path(__file__).resolve().parents[1]
CAS = ELLIPTIC_ROOT / "cas"
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))

from elliptic_candidate_record import (  # noqa: E402
    WeierstrassChange,
    build_finite_quotient_certificate,
    is_on_weierstrass_curve,
    source_point_to_target,
    verify_finite_quotient_certificate,
    weierstrass_invariants,
)
from finite_quotient_escape import QuotientBlock, analyze_escape  # noqa: E402


Q = Fraction
RationalPoint = tuple[Fraction, Fraction]
RATPOINTS_X_PATTERN = re.compile(r"\((-?\d+)\s*:\s*(\d+)\)")


def completed_square_coefficients(
    model: Sequence[Fraction | int | str],
) -> tuple[int, int, int, int]:
    """Return low-to-high coefficients of ``Y^2=4x^3+b2*x^2+2b4*x+b6``."""

    coefficients = tuple(Q(value) for value in model)
    if len(coefficients) != 5 or any(value.denominator != 1 for value in coefficients):
        raise ValueError("completed-square search requires an integral model")
    invariants = weierstrass_invariants(coefficients)
    return (
        int(invariants["b6"]),
        int(2 * invariants["b4"]),
        int(invariants["b2"]),
        4,
    )


def evaluate_polynomial(coefficients: Sequence[int], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def affine_substitute_polynomial(
    coefficients_low_to_high: Sequence[Fraction | int],
    center: Fraction,
    scale: Fraction,
) -> tuple[Fraction, ...]:
    """Return the coefficients of ``f(center + scale*X)`` exactly."""

    center = Q(center)
    scale = Q(scale)
    if scale == 0:
        raise ValueError("an affine search chart needs nonzero scale")
    answer = [Q(0)] * len(coefficients_low_to_high)
    for degree, coefficient in enumerate(coefficients_low_to_high):
        coefficient = Q(coefficient)
        for x_degree in range(degree + 1):
            answer[x_degree] += (
                coefficient
                * comb(degree, x_degree)
                * center ** (degree - x_degree)
                * scale**x_degree
            )
    return tuple(answer)


def integral_square_scaled_coefficients(
    coefficients: Sequence[Fraction | int],
) -> tuple[tuple[int, ...], int]:
    """Clear coefficient denominators by the square ``D^2``.

    The transformed ratpoints equation is ``Z^2=D^2*f(X)`` with ``Z=D*Y``.
    """

    denominator = 1
    rational_coefficients = tuple(Q(value) for value in coefficients)
    for value in rational_coefficients:
        denominator = lcm(denominator, value.denominator)
    integral = tuple(int(value * denominator**2) for value in rational_coefficients)
    if any(
        Q(integer) != value * denominator**2
        for integer, value in zip(integral, rational_coefficients)
    ):
        raise AssertionError("square scaling failed to clear chart denominators")
    return integral, denominator


def rational_square_root(value: Fraction) -> Fraction | None:
    value = Q(value)
    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator**2 != value.numerator or denominator**2 != value.denominator:
        return None
    return Q(numerator, denominator)


def parse_ratpoints_abscissae(output: str) -> tuple[Fraction, ...]:
    """Parse quiet ``-y`` output, checking canonical reduced projective pairs."""

    answer = []
    seen = set()
    for line in output.splitlines():
        if not line.strip():
            continue
        match = RATPOINTS_X_PATTERN.fullmatch(line.strip())
        if match is None:
            raise ValueError(f"unexpected ratpoints output line: {line!r}")
        numerator, denominator = map(int, match.groups())
        if denominator == 0:
            if numerator != 1:
                raise ValueError("ratpoints emitted a noncanonical point at infinity")
            continue
        value = Q(numerator, denominator)
        if value not in seen:
            seen.add(value)
            answer.append(value)
    return tuple(answer)


def points_from_completed_square_abscissa(
    model: Sequence[Fraction | int | str], x_coordinate: Fraction
) -> tuple[RationalPoint, ...]:
    """Map an exact ratpoints abscissa back to the generalized minimal model."""

    coefficients = tuple(Q(value) for value in model)
    a1, _a2, a3, _a4, _a6 = coefficients
    square = evaluate_polynomial(completed_square_coefficients(coefficients), x_coordinate)
    root = rational_square_root(square)
    if root is None:
        raise AssertionError("a ratpoints abscissa failed exact square reconstruction")
    roots = (root,) if root == 0 else (root, -root)
    points = tuple(
        (x_coordinate, (ordinate - a1 * x_coordinate - a3) / 2)
        for ordinate in roots
    )
    if any(not is_on_weierstrass_curve(coefficients, point) for point in points):
        raise AssertionError("a completed-square point missed the generalized model")
    return points


def negate_point(
    model: Sequence[Fraction | int | str], point: RationalPoint
) -> RationalPoint:
    a1, _a2, a3, _a4, _a6 = (Q(value) for value in model)
    x_coordinate, y_coordinate = point
    return x_coordinate, -y_coordinate - a1 * x_coordinate - a3


def sign_key(
    model: Sequence[Fraction | int | str], point: RationalPoint
) -> tuple[tuple[int, int], tuple[int, int]]:
    def key(value: RationalPoint) -> tuple[tuple[int, int], tuple[int, int]]:
        return (
            (value[0].numerator, value[0].denominator),
            (value[1].numerator, value[1].denominator),
        )

    return min(key(point), key(negate_point(model, point)))


def novel_points_up_to_sign(
    model: Sequence[Fraction | int | str],
    baseline: Sequence[RationalPoint],
    searched: Sequence[RationalPoint],
) -> tuple[RationalPoint, ...]:
    seen = {sign_key(model, point) for point in baseline}
    answer = []
    for point in searched:
        key = sign_key(model, point)
        if key in seen:
            continue
        seen.add(key)
        answer.append(point)
    return tuple(answer)


def point_record(point: RationalPoint) -> list[str]:
    return [str(point[0]), str(point[1])]


def exact_escape_records(
    specialization_artifact: dict[str, Any],
    candidate_minimal_points: Sequence[RationalPoint],
    relation_primes: Sequence[int],
    reduction_prime_bound: int,
) -> dict[str, Any]:
    """Measure exact candidate escape from the specialized generic rank 17."""

    certificate_block = specialization_artifact["finite_quotient_independence"]
    certificate_model = tuple(
        Q(value) for value in certificate_block["certificate_short_model"]
    )
    baseline_points = tuple(
        (Q(point[0]), Q(point[1])) for point in certificate_block["points"]
    )
    if len(baseline_points) != 17:
        raise ValueError("the specialization artifact does not contain 17 baseline points")
    minimal_to_short = WeierstrassChange.from_values(
        certificate_block["minimal_to_certificate_short_change_u_r_s_t"]
    )
    candidate_points = tuple(
        source_point_to_target(point, minimal_to_short)
        for point in candidate_minimal_points
    )
    if any(
        not is_on_weierstrass_curve(certificate_model, point)
        for point in candidate_points
    ):
        raise AssertionError("a searched point missed the certificate short model")
    labels = tuple(f"searched-candidate-{index}" for index in range(len(candidate_points)))
    all_points = (*baseline_points, *candidate_points)
    attempts = []
    maximum_marginal_dimension = 0
    escape_basis_labels: list[str] = []
    if candidate_points:
        for relation_prime in relation_primes:
            certificate = build_finite_quotient_certificate(
                certificate_model,
                all_points,
                relation_prime=int(relation_prime),
                prime_bound=reduction_prime_bound,
            )
            verify_finite_quotient_certificate(
                certificate_model, all_points, certificate
            )
            blocks = tuple(
                QuotientBlock.build(
                    modulus=int(relation_prime),
                    rows=signature["rows"],
                    column_count=len(all_points),
                    source=f"good-reduction-p={signature['prime']}",
                )
                for signature in certificate["signatures"]
            )
            profile = analyze_escape(
                blocks, known_column_count=17, candidate_labels=labels
            )
            profile_record = profile.to_record()
            profile_record["finite_quotient_certificate"] = certificate
            attempts.append(profile_record)
            if profile.marginal_dimension > maximum_marginal_dimension:
                maximum_marginal_dimension = profile.marginal_dimension
                escape_basis_labels = list(profile.independent_escape_basis_labels)
    return {
        "known_column_count": 17,
        "candidate_count": len(candidate_points),
        "candidate_points_on_certificate_short_model": [
            point_record(point) for point in candidate_points
        ],
        "relation_prime_profiles": attempts,
        "maximum_marginal_dimension": maximum_marginal_dimension,
        "independent_escape_basis_labels": escape_basis_labels,
        "promotion_threshold": 15,
        "promotion_eligible": maximum_marginal_dimension >= 15,
        "promotion_status": (
            "ELIGIBLE_FOR_32_POINT_CERTIFICATION"
            if maximum_marginal_dimension >= 15
            else "NOT_PROMOTED_QUOTIENT_GAIN_BELOW_15"
        ),
        "claim_boundary": (
            "positive finite-quotient escape is exact; bounded non-escape is not "
            "a dependence proof"
        ),
    }
