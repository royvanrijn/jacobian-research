#!/usr/bin/env python3
"""Rejected mixed-label rank-13 entry point; retained helpers support its audit.

Sorting roots and choosing positive ordinate branches independently at each
specialization changes column identities. The old stacked matrix cannot prove
generic rank13. EC-MF2S13 records the corrected coherent lower bound11.
Use verify_mestre_parent_and_label_audits.sage --check for the current replay.
The original source and result are preserved in the cleanup archive.
"""

from __future__ import annotations

if __name__ == "__main__":
    raise SystemExit("REJECTED rank-13 certificate: mixed generic section labels. Use verify_mestre_parent_and_label_audits.sage --check; corrected lower bound11.")

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

# Historical probe roster. Sorted visible positions and positive extra signs
# are not fixed generic labels; the coherent audit corrects both separately.
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
    """Prevent import callers from reviving the rejected generic-rank claim."""
    raise RuntimeError(
        "REJECTED rank-13 certificate: mixed section labels; "
        "use verify_mestre_parent_and_label_audits.sage --check for lower bound11"
    )
