#!/usr/bin/env python3
"""Linear-abscissa sections on Nagao's rank-13 Mestre quartic.

Nagao singles out the section with abscissa ``(T+703)/15``.  The same
quartic has five companion sections whose abscissae are also linear in
``T``.  They do not increase the generic Mordell--Weil rank: their Jacobian
images satisfy the exact relations recorded below.

The separate verifier ``verify_nagao_linear_sections.py`` proves that these
six companion sections, together with the twelve Mestre sections, exhaust
the ansatz

``x = m*T+n,  y in Q[T],  degree(y) <= 3``.

This module only supplies exact formulas and specialization checks.  It does
not claim that the displayed sections generate the Mordell--Weil group.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

from nagao_1994 import (
    RANK13_CONSTRUCTION,
    primitive_quartic_coefficients,
    quartic_point_to_short_jacobian,
    quartic_value,
)


Q = Fraction


@dataclass(frozen=True)
class LinearQuarticSection:
    """A section ``x=m*T+n`` with a polynomial ordinate."""

    label: str
    slope: Fraction
    intercept: Fraction
    ordinate_coefficients: tuple[Fraction, ...]

    def point(self, parameter: Fraction) -> tuple[Fraction, Fraction]:
        parameter = Q(parameter)
        x = self.slope * parameter + self.intercept
        y = Q(0)
        for coefficient in reversed(self.ordinate_coefficients):
            y = y * parameter + coefficient
        coefficients = primitive_quartic_coefficients(
            RANK13_CONSTRUCTION, parameter
        )
        if y**2 != quartic_value(coefficients, x):
            raise AssertionError(f"linear section {self.label} failed exactly")
        return x, y

    def jacobian_point(self, parameter: Fraction) -> tuple[Fraction, Fraction]:
        return quartic_point_to_short_jacobian(
            RANK13_CONSTRUCTION, Q(parameter), self.point(parameter)
        )


LINEAR_COMPANION_SECTIONS = (
    LinearQuarticSection(
        "plus-1/15",
        Q(1, 15),
        Q(703, 15),
        (Q(2161725, 75), Q(900484, 75), Q(-844, 75), Q(-224, 75)),
    ),
    LinearQuarticSection(
        "minus-1/15",
        Q(-1, 15),
        Q(703, 15),
        (Q(2161725, 75), Q(-900484, 75), Q(-844, 75), Q(224, 75)),
    ),
    LinearQuarticSection(
        "plus-7/15",
        Q(7, 15),
        Q(928, 15),
        (Q(-33677700, 75), Q(-482659, 75), Q(2758, 75), Q(176, 75)),
    ),
    LinearQuarticSection(
        "minus-7/15",
        Q(-7, 15),
        Q(928, 15),
        (Q(33677700, 75), Q(-482659, 75), Q(-2758, 75), Q(176, 75)),
    ),
    LinearQuarticSection(
        "plus-5/3",
        Q(5, 3),
        Q(3628, 15),
        (Q(-915090000, 75), Q(-19161259, 75), Q(-125150, 75), Q(-400, 75)),
    ),
    LinearQuarticSection(
        "minus-5/3",
        Q(-5, 3),
        Q(3628, 15),
        (Q(915090000, 75), Q(-19161259, 75), Q(125150, 75), Q(-400, 75)),
    ),
)


# The relation basis is visible sections 0,...,10 in the ordering returned by
# primitive_visible_points, followed by Nagao's plus-1/15 section.  Visible
# section 11 is omitted because the twelve Mestre images already have a
# relation.  Coefficients below mean Q = sum_i coefficient_i * P_i.
COMPANION_RELATION_BASIS_VISIBLE_INDICES = tuple(range(11))
COMPANION_JACOBIAN_RELATIONS: dict[str, tuple[int, ...]] = {
    "minus-1/15": (0, 0, 0, 0, 0, 0, 1, -1, -1, 1, 0, -1),
    "plus-7/15": (0, 0, 0, 0, 0, -1, 0, -1, -1, 0, -1, -1),
    "minus-7/15": (1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, -1),
    "plus-5/3": (0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, -1),
    "minus-5/3": (1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, -1),
}


def omitted_companion_sections() -> tuple[LinearQuarticSection, ...]:
    """Return the five sections absent from ``rank13_known_quartic_points``."""

    return LINEAR_COMPANION_SECTIONS[1:]


def companion_abscissae(parameter: Fraction) -> tuple[Fraction, ...]:
    """Return all six exact companion abscissae at ``parameter``."""

    return tuple(section.point(parameter)[0] for section in LINEAR_COMPANION_SECTIONS)


def point_on_short_curve(
    coefficients: Sequence[Fraction], point: tuple[Fraction, Fraction]
) -> bool:
    """Check a point on an extended Weierstrass model exactly."""

    a1, a2, a3, a4, a6 = (Q(value) for value in coefficients)
    x, y = (Q(value) for value in point)
    return y**2 + a1 * x * y + a3 * y == x**3 + a2 * x**2 + a4 * x + a6
