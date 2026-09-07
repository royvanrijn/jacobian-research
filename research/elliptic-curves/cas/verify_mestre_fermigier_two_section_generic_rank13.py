#!/usr/bin/env python3
"""Certify thirteen independent sections on the Fermigier two-section curve.

The exact component identity supplies twelve visible Mestre sections, the
Fermigier line, and the reconstructed second line over ``Q(u)(T)``.  This
checker selects the first eleven visible points and both affine sections.

For a putative relation between these thirteen generic sections, every good
specialization has the same relation.  The thirteen fixed exact finite
quotients below make its coefficient vector zero modulo 3.  The generic
curve has no rational 3-torsion: at the smooth specialization ``u=-5,T=1``
its reduction modulo 19 has order 28.  Hence the relation is divisible by 3,
and repeated division proves it is zero.  This is a generic rank-at-least-13
certificate, not a saturation, Shioda, or rank-upper-bound computation.

All section coordinates are obtained from the compact root formulas and the
triangular square-root recurrence; no universal two-section residual is
expanded.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from mod2_reduction_independence import _reduce_rational, finite_curve_points
from screen_mestre_fermigier_two_section_height_triage import specialized_points
from search_mestre_root_tuple_scale_max200 import (
    gf_l_rank_and_pivots,
    mod_l_reduction_signature,
)


Q = Fraction
MODULUS = 3
TORSION_SPECIALIZATION = (Q(-5), Q(1))
TORSION_REDUCTION_PRIME = 19

# Each listed quotient raises the combined rank by one.  The first eleven
# columns are visible points in their fixed primitive order; columns twelve
# and thirteen are respectively Fermigier's line and the new second line.
PROBES = (
    (Q(-5), Q(1), 23),
    (Q(-5), Q(1), 41),
    (Q(-5), Q(1), 59),
    (Q(-5), Q(1), 67),
    (Q(-5), Q(1), 71),
    (Q(-5), Q(1), 73),
    (Q(-5), Q(1), 79),
    (Q(-5), Q(1), 83),
    (Q(-5), Q(1), 103),
    (Q(-5), Q(1), 109),
    (Q(-5), Q(2), 37),
    (Q(-3), Q(1), 41),
    (Q(-1, 2), Q(1), 7),
)
EXPECTED_GROUP_ORDERS = (30, 54, 72, 78, 84, 87, 81, 99, 117, 120, 48, 54, 12)


def rational_text(value: Fraction) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def basis_points(
    points: tuple[tuple[Fraction, Fraction], ...],
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Remove only the twelfth visible point from the fourteen-point list."""

    if len(points) != 14:
        raise AssertionError("the Fermigier component no longer supplied fourteen points")
    return points[:11] + points[12:]


def verify_point_equations(
    coefficients: tuple[Fraction, ...],
    points: tuple[tuple[Fraction, Fraction], ...],
) -> None:
    if len(coefficients) != 5 or any(coefficients[:3]):
        raise AssertionError("the expected short Weierstrass model changed")
    coefficient_a, coefficient_b = coefficients[3:]
    for x_value, y_value in points:
        if y_value**2 != x_value**3 + coefficient_a * x_value + coefficient_b:
            raise AssertionError("a selected generic section failed a specialization equation")


def torsion_exclusion() -> dict[str, object]:
    """Certify that the generic curve has no rational 3-torsion."""

    u, parameter = TORSION_SPECIALIZATION
    coefficients, points = specialized_points(u, parameter)
    verify_point_equations(coefficients, points)
    coefficient_a = _reduce_rational(coefficients[3], TORSION_REDUCTION_PRIME)
    coefficient_b = _reduce_rational(coefficients[4], TORSION_REDUCTION_PRIME)
    discriminant = -16 * (4 * coefficient_a**3 + 27 * coefficient_b**2)
    if discriminant % TORSION_REDUCTION_PRIME == 0:
        raise AssertionError("the torsion-exclusion reduction became bad")
    group_order = len(
        finite_curve_points(coefficient_a, coefficient_b, TORSION_REDUCTION_PRIME)
    )
    if group_order != 28 or group_order % MODULUS == 0:
        raise AssertionError("the generic 3-torsion exclusion changed")
    return {
        "u": rational_text(u),
        "T": rational_text(parameter),
        "reduction_prime": TORSION_REDUCTION_PRIME,
        "finite_group_order": group_order,
        "conclusion": "the specialization has no rational 3-torsion; smooth specialization excludes generic rational 3-torsion",
    }


def replay() -> dict[str, object]:
    rows: list[tuple[int, ...]] = []
    records = []
    for index, ((u, parameter, prime), expected_order) in enumerate(
        zip(PROBES, EXPECTED_GROUP_ORDERS), start=1
    ):
        coefficients, points = specialized_points(u, parameter)
        selected = basis_points(points)
        verify_point_equations(coefficients, selected)
        signature = mod_l_reduction_signature(
            coefficients, selected, prime, modulus=MODULUS
        )
        if signature.group_order != expected_order or signature.quotient_dimension != 1:
            raise AssertionError("a frozen Fermigier quotient probe changed")
        before_rank, _ = gf_l_rank_and_pivots(rows, len(selected), MODULUS)
        candidate_rows = [*rows, *signature.rows]
        after_rank, _ = gf_l_rank_and_pivots(
            candidate_rows, len(selected), MODULUS
        )
        if after_rank != index or after_rank != before_rank + 1:
            raise AssertionError("a quotient probe no longer raises rank by one")
        rows.extend(signature.rows)
        records.append(
            {
                "u": rational_text(u),
                "T": rational_text(parameter),
                "reduction_prime": prime,
                "finite_group_order": signature.group_order,
                "quotient_rows": [list(row) for row in signature.rows],
                "combined_rank_after_probe": after_rank,
            }
        )
    rank, pivots = gf_l_rank_and_pivots(rows, 13, MODULUS)
    if rank != 13 or pivots != tuple(range(13)):
        raise AssertionError("the Fermigier generic independence matrix changed")
    return {
        "status": "generic rank-at-least-13 finite-reduction certificate verified",
        "base_field": "Q(u)(T), on v=(u^2+u+2)/u",
        "basis": [
            "the first eleven primitive visible sections",
            "the Fermigier affine section",
            "the reconstructed second affine section",
        ],
        "section_count": 13,
        "descent_modulus": MODULUS,
        "torsion_exclusion": torsion_exclusion(),
        "quotient_probes": records,
        "combined_exact_rank_over_F3": rank,
        "conclusion": "the two affine sections, together with eleven visible sections, are independent over Q(u)(T); the Fermigier two-section component has generic Mordell-Weil rank at least 13",
        "not_established": [
            "generic rank at least 14 or independence from any separate rank-13 family",
            "saturation, pair intersections beyond the recorded finite meeting, or a Shioda Gram matrix",
            "a generic rank upper bound",
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
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
