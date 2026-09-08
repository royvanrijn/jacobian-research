#!/usr/bin/env python3
"""Exact generic data for Nagao's 1994 section-7 Mestre family.

Nagao's paper uses the roots ``(346,260,255,146,55,0)`` and writes the
special parameter as ``t=5081/94``.  The symmetric six-root constructor in
this repository uses ``T=2t``, hence ``T=5081/47`` at that specialization.

Besides the twelve Mestre points ``x=r_i +/- T``, the primitive quartic has
six further polynomial sections with linear abscissa.  The companion formulas
below are generic identities over ``Q(T)``.  Their Mordell--Weil relations and
the fact that the displayed list exhausts the linear-abscissa ansatz are
proved by ``verify_nagao_section7_linear_sections.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from mestre_root_tuples import SixRootMestreConstruction
from nagao_1994 import (
    PRIMARY_SOURCE,
    primitive_quartic_coefficients,
    quartic_point_to_short_jacobian,
    quartic_value,
)


Q = Fraction

SECTION7_ROOTS = (346, 260, 255, 146, 55, 0)
SECTION7_CONSTRUCTION = SixRootMestreConstruction(
    tuple(Q(root) for root in SECTION7_ROOTS)
)
SECTION7_PAPER_PARAMETER = Q(5081, 94)
SECTION7_CONSTRUCTOR_PARAMETER = 2 * SECTION7_PAPER_PARAMETER


def section7_primitive_quartic_coefficients(
    parameter: Fraction,
) -> tuple[Fraction, ...]:
    """Return the explicit primitive quartic coefficients in ascending x-order."""

    parameter = Q(parameter)
    coefficients = (
        9 * parameter**6
        - 910748 * parameter**4
        + 23718659440 * parameter**2
        + 557726319412900,
        18
        * (
            354 * parameter**4
            - 17901331 * parameter**2
            - 884640359570
        ),
        -3
        * (
            6 * parameter**4
            - 668642 * parameter**2
            - 52052853547
        ),
        -54 * (118 * parameter**2 + 11538729),
        9 * (parameter**2 + 96714),
    )
    generic = primitive_quartic_coefficients(SECTION7_CONSTRUCTION, parameter)
    if coefficients != generic:
        raise AssertionError("the explicit section-7 quartic normalization drifted")
    return coefficients


@dataclass(frozen=True)
class Section7LinearQuarticSection:
    """A section ``x=m*T+n`` with an ordinate in ``Q[T]``."""

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
        coefficients = section7_primitive_quartic_coefficients(parameter)
        if y**2 != quartic_value(coefficients, x):
            raise AssertionError(f"section-7 linear section {self.label} failed")
        return x, y

    def jacobian_point(self, parameter: Fraction) -> tuple[Fraction, Fraction]:
        return quartic_point_to_short_jacobian(
            SECTION7_CONSTRUCTION, Q(parameter), self.point(parameter)
        )


SECTION7_LINEAR_COMPANION_SECTIONS = (
    Section7LinearQuarticSection(
        "plus-7/27",
        Q(7, 27),
        Q(6920, 27),
        (Q(84770), Q(-18554923, 243), Q(-29974, 243), Q(680, 243)),
    ),
    Section7LinearQuarticSection(
        "minus-7/27",
        Q(-7, 27),
        Q(6920, 27),
        (Q(-84770), Q(-18554923, 243), Q(29974, 243), Q(680, 243)),
    ),
    Section7LinearQuarticSection(
        "plus-17/27",
        Q(17, 27),
        Q(5462, 27),
        (Q(5138284), Q(-6202747, 243), Q(-23222, 243), Q(440, 243)),
    ),
    Section7LinearQuarticSection(
        "minus-17/27",
        Q(-17, 27),
        Q(5462, 27),
        (Q(-5138284), Q(-6202747, 243), Q(23222, 243), Q(440, 243)),
    ),
    Section7LinearQuarticSection(
        "plus-43/27",
        Q(43, 27),
        Q(-4015, 27),
        (Q(94091525), Q(-236806588, 243), Q(756284, 243), Q(-1120, 243)),
    ),
    Section7LinearQuarticSection(
        "minus-43/27",
        Q(-43, 27),
        Q(-4015, 27),
        (Q(-94091525), Q(-236806588, 243), Q(-756284, 243), Q(-1120, 243)),
    ),
)


@dataclass(frozen=True)
class Section7QuadraticQuarticSection:
    """A section ``x=m*T^2+n*T+k`` with an ordinate in ``Q[T]``."""

    label: str
    quadratic_coefficient: Fraction
    linear_coefficient: Fraction
    constant_coefficient: Fraction
    ordinate_coefficients: tuple[Fraction, ...]

    def point(self, parameter: Fraction) -> tuple[Fraction, Fraction]:
        parameter = Q(parameter)
        x = (
            self.quadratic_coefficient * parameter**2
            + self.linear_coefficient * parameter
            + self.constant_coefficient
        )
        y = Q(0)
        for coefficient in reversed(self.ordinate_coefficients):
            y = y * parameter + coefficient
        coefficients = section7_primitive_quartic_coefficients(parameter)
        if y**2 != quartic_value(coefficients, x):
            raise AssertionError(f"section-7 quadratic section {self.label} failed")
        return x, y

    def jacobian_point(self, parameter: Fraction) -> tuple[Fraction, Fraction]:
        return quartic_point_to_short_jacobian(
            SECTION7_CONSTRUCTION, Q(parameter), self.point(parameter)
        )


SECTION7_QUADRATIC_COMPANION_SECTIONS = (
    Section7QuadraticQuarticSection(
        "quadratic-plus-56/5373",
        Q(56, 5373),
        Q(0),
        Q(1389190, 5373),
        (
            Q(0),
            Q(684218797630, 9623043),
            Q(0),
            Q(171853351, 9623043),
            Q(0),
            Q(3136, 9623043),
        ),
    ),
    Section7QuadraticQuarticSection(
        "quadratic-minus-22/5373",
        Q(-22, 5373),
        Q(0),
        Q(1389190, 5373),
        (
            Q(0),
            Q(638541742570, 9623043),
            Q(0),
            Q(-24743777, 9623043),
            Q(0),
            Q(484, 9623043),
        ),
    ),
    Section7QuadraticQuarticSection(
        "quadratic-minus-34/5373",
        Q(-34, 5373),
        Q(0),
        Q(1389190, 5373),
        (
            Q(0),
            Q(-631221199130, 9623043),
            Q(0),
            Q(-2763929, 9623043),
            Q(0),
            Q(1156, 9623043),
        ),
    ),
)


# The relation basis is visible sections 0,...,10 followed by plus-7/27.
# Each tuple records Q = sum_i c_i P_i.  The verifier proves these identities
# symbolically on the generic Jacobian, so bounded searches can safely treat
# the five companions as predeclared dependent sections.
SECTION7_RELATION_BASIS_VISIBLE_INDICES = tuple(range(11))
SECTION7_JACOBIAN_RELATIONS: dict[str, tuple[int, ...]] = {
    "visible-11": (-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0),
    "minus-7/27": (0, 0, -1, 1, 1, -1, 0, 0, 0, 0, 0, 1),
    "plus-17/27": (-1, 0, -1, 0, 0, -1, 0, -1, 0, 0, 0, 1),
    "minus-17/27": (0, -1, -1, 0, 0, -1, -1, 0, 0, 0, 0, 1),
    "plus-43/27": (-1, -1, -1, 0, 0, -1, -1, -1, 0, -1, -1, 1),
    "minus-43/27": (0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1),
}

SECTION7_QUADRATIC_JACOBIAN_RELATIONS: dict[str, tuple[int, ...]] = {
    "quadratic-plus-56/5373": (0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 1),
    "quadratic-minus-22/5373": (0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1),
    "quadratic-minus-34/5373": (
        -1,
        -1,
        -1,
        0,
        0,
        -1,
        -1,
        -1,
        0,
        0,
        0,
        1,
    ),
}


def section7_companion_abscissae(parameter: Fraction) -> tuple[Fraction, ...]:
    """Return all six generic companion abscissae at ``parameter``."""

    return tuple(
        section.point(Q(parameter))[0]
        for section in SECTION7_LINEAR_COMPANION_SECTIONS
    )


def section7_quadratic_companion_abscissae(
    parameter: Fraction,
) -> tuple[Fraction, ...]:
    """Return the three exact quadratic companion abscissae."""

    return tuple(
        section.point(Q(parameter))[0]
        for section in SECTION7_QUADRATIC_COMPANION_SECTIONS
    )
