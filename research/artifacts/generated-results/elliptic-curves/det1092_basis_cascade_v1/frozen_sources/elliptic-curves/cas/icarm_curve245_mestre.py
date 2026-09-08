#!/usr/bin/env python3
"""Exact Mestre parent data for ICARM curve #245.

The public commentary gives only ``T=5801/160``.  Exact comparison of
Fermigier's six-root formulas recovers the missing parameters
``(u,v)=(3/2,2)``.  Scaling the roots and ``T`` by 16 and translating the
least root to zero gives the integral presentation used here.

This module records identities and coordinates only.  It makes no rank claim
beyond those established by :mod:`verify_icarm_curve245_rank20`.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Sequence

from mestre_root_tuples import SixRootMestreConstruction


Q = Fraction


def _parameter(value: Any) -> Any:
    """Coerce ordinary inputs to Q while preserving compatible formal values."""

    try:
        return Q(value)
    except TypeError:
        return value

FERMIGIER_U = Q(3, 2)
FERMIGIER_V = Q(2)
NATIVE_PARAMETER = Q(5801, 160)
AFFINE_SCALE = Q(16)
CANONICAL_PARAMETER = AFFINE_SCALE * NATIVE_PARAMETER

NATIVE_ROOTS = tuple(
    map(Q, ("-375/16", "-269/16", "-31/16", "25/4", "219/16", "89/4"))
)
CANONICAL_ROOTS = (0, 106, 344, 475, 594, 731)

PUBLIC_MODEL = (
    Q(1),
    Q(-1),
    Q(1),
    Q(-25880411472355347134118026792),
    Q(1606663697747901005185875883284420820193259),
)

# Primitive short-Jacobian coefficients in ascending powers of the canonical
# parameter.  Odd coefficients vanish because T and -T give the same curve.
A_COEFFICIENTS = (
    -14557427332128769150107,
    0,
    -971401897768443192,
    0,
    11800027144800,
    0,
    -21555072,
    0,
    -432,
)
B_COEFFICIENTS = (
    648361633424537722009431313503306,
    0,
    82538787133870919273441032704,
    0,
    -1374614482522022028018828,
    0,
    19145150140414965408,
    0,
    -151279932590208,
    0,
    258660864,
    0,
    3456,
)

# Fermigier's thirteenth generic section on the primitive canonical quartic:
# x(T)=EXTRA_X[0]+EXTRA_X[1]*T and y(T)=sum EXTRA_Y[i]*T^i.
EXTRA_X = (Q(4558, 29), Q(7, 29))
EXTRA_Y = (Q(4801853), Q(-123478023, 841), Q(88438, 841), Q(792, 841))

# Sage's exact isomorphism from the primitive short Jacobian at T=5801/10 to
# the public minimal model, in (u,r,s,t) convention.
ANCHOR_SHORT_TO_PUBLIC_CHANGE = (Q(3, 25), Q(-9, 2500), Q(3, 50), Q(27, 31250))


def fermigier_roots(u: Fraction, v: Fraction) -> tuple[Fraction, ...]:
    """Return Fermigier's six labelled roots at exact rational ``(u,v)``."""

    u, v = _parameter(u), _parameter(v)
    alpha1 = (
        -v + v**2 - v**3 + v**4 + 2*u*v + v**2*u - 2*v**3*u
        + v**4*u + u**2 + u**2*v - 2*v**2*u**2 - 2*v**3*u**2
        + u**3*v**2 - u**4
    )
    alpha2 = (
        v**3*u**2 - 2*v**2*u**2 + v**2*u - 2*u**3*v**2
        - 2*u**3*v + u**2*v + u**4*v + 2*u*v - v**4 - u + v**2
        - u**3 + u**2 + u**4
    )
    alpha3 = (
        -v**4*u + 2*v**3*u + v**3*u**2 + v**2*u**2 + v**2*u
        - u**3*v**2 - 2*u**3*v - 2*u**2*v + u**4*v - v + v**3
        - 2*u**3 + u**2 + u**4
    )
    alpha4 = (
        v - v**2 + v**3 - v**4 + u - 2*u*v + v**2*u + 2*v**3*u
        - 2*u**2 - 2*u**2*v + v**2*u**2 + v**3*u**2 + u**3 - u**4*v
    )
    alpha5 = (
        v**4*u - 2*v**3*u - v**3*u**2 + v**2*u**2 - 2*v**2*u
        + u**3*v**2 + 2*u**3*v + u**2*v + v**4 - u - 2*v**3
        + v**2 + u**3 - u**4*v
    )
    alpha6 = (
        v - 2*v**2 + v**3 + u - 2*u*v - 2*v**2*u - v**4*u - u**2
        + u**2*v + v**2*u**2 + u**3 + 2*u**3*v + u**3*v**2 - u**4
    )
    return alpha1, alpha2, alpha3, alpha4, alpha5, alpha6


def fermigier_extra_line(u: Fraction, v: Fraction) -> tuple[Fraction, Fraction]:
    """Return Fermigier's extra affine abscissa ``A(u,v)+B(u,v)T``.

    The formula is the same source-coordinate line used in the exact
    Fermigier affine-section verifier.  Callers applying an affine change to
    the roots must translate the intercept and retain the slope.
    """

    u, v = _parameter(u), _parameter(v)
    denominator = u**2 + v**2 + 1
    intercept = (
        3*u**2*v + 3*u**3*v**2 + 2*u**4*v + u*v - 4*v**3*u
        - 3*v**2*u**2 + 3*v**3*u**2 - 4*u**3*v + v**3 + 2*v**4*u
        - 3*u**4 + u**3 - 3*v**4 + u + 3*v**2*u + v - u**4*v**3
        - 2*u**4*v**2 - u**3*v**4 - 4*u**3*v**3 + u**2*v**5
        - 2*u**2*v**4 + u*v**5 + u**5*v**2 + u**5*v - v**6 - u**6
        + 2*v**5 + 2*u**5
    ) / denominator
    slope = (-u**2 - v**2 + 2*u + 2*v + 1) / denominator
    return intercept, slope


def evaluate_polynomial(coefficients: Sequence[int | Fraction], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * Q(value) + Q(coefficient)
    return answer


def primitive_short_model(parameter: Fraction) -> tuple[Fraction, ...]:
    """Return ``[0,0,0,A(T),B(T)]`` in the canonical integral-root chart."""

    parameter = Q(parameter)
    return (
        Q(0), Q(0), Q(0),
        evaluate_polynomial(A_COEFFICIENTS, parameter),
        evaluate_polynomial(B_COEFFICIENTS, parameter),
    )


def extra_quartic_point(parameter: Fraction) -> tuple[Fraction, Fraction]:
    """Return Fermigier's exact extra generic point on the primitive quartic."""

    parameter = Q(parameter)
    return (
        EXTRA_X[0] + EXTRA_X[1] * parameter,
        evaluate_polynomial(EXTRA_Y, parameter),
    )


CONSTRUCTION = SixRootMestreConstruction(tuple(map(Q, CANONICAL_ROOTS)))
