#!/usr/bin/env python3
"""Formal local continuation of the unexplained Fermigier two-section seed.

The all-prime affine-section screen finds the normalized roots
``(0,8,58,77,85,102)`` with the two affine abscissae

    -58/11 - 17/11 T,     247/44 - 3/11 T.

The first is Fermigier's generic extra section at ``(u,v)=(-3,-8/3)``.
This script keeps that two-parameter root surface and locally imposes the
second section.  It works only with the recursive square-root residuals from
``probe_mestre_two_section_local_continuation``; no expanded elimination
residual is built.

The output is a formal Q[[t]] branch with ``u=-3+t`` when the transverse
three-by-three Jacobian in ``(v,a,b)`` is nonsingular.  Such a branch is local
evidence only: it is neither a global component equation nor a rank claim.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Sequence

from icarm_curve245_mestre import fermigier_extra_line, fermigier_roots
from probe_mestre_two_section_local_continuation import (
    Field,
    FormalSeries,
    residuals_from_jets,
    small_pade_models,
    solve_square_over_q,
)


Q = Fraction
BASE_U = Q(-3)
BASE_V = Q(-8, 3)
SECOND_SECTION = (Q(247, 44), Q(-3, 11))
SECOND_EQUATIONS = (4, 5, 6)


def normalized_data(u: object, v: object) -> tuple[object, ...]:
    """Return normalized moduli and Fermigier's extra line at ``(u,v)``.

    The arithmetic is deliberately generic: it accepts Q values and the
    truncated series class used for continuation.
    """

    roots = fermigier_roots(u, v)
    scale = roots[1] - roots[0]
    roots = tuple((root - roots[0]) / scale for root in roots)
    others = roots[2:]
    c1 = -sum(others)
    c2 = sum(others[left] * others[right] for left in range(4) for right in range(left + 1, 4))
    c3 = -sum(
        others[first] * others[second] * others[third]
        for first in range(4)
        for second in range(first + 1, 4)
        for third in range(second + 1, 4)
    )
    c4 = others[0] * others[1] * others[2] * others[3]
    intercept, slope = fermigier_extra_line(u, v)
    # The normalized Mestre parameter is T/scale, so a root-coordinate
    # normalization changes the intercept but leaves the line slope intact.
    return c1, c2, c3, c4, (intercept - fermigier_roots(u, v)[0]) / scale, slope


def base_coordinates() -> tuple[Fraction, ...]:
    """Return the labelled two-section point used by the all-prime screen."""

    data = normalized_data(BASE_U, BASE_V)
    return tuple(map(Q, data[:4])) + (Q(data[4]), Q(data[5])) + SECOND_SECTION


def reconstructed_second_line(u: object) -> tuple[object, object, object]:
    """The rational model recognized from the order-18 local branch."""

    v = (u**2 + u + 2) / u
    intercept = -(
        2*u**6 - u**5 + 4*u**4 - u**3 - 8*u**2 - 4*u - 16
    ) / (2 * (u + 2) * (u**2 + 2) * (u**2 - u + 4))
    slope = u / (u**2 + 2)
    return v, intercept, slope


def series_coordinates(
    coefficients: Sequence[Sequence[Fraction]], field: Field
) -> list[FormalSeries]:
    """Construct the eight recursive-residual inputs from local coordinates."""

    u = FormalSeries.seed(field, coefficients[0])
    v = FormalSeries.seed(field, coefficients[1])
    c1, c2, c3, c4, first_a, first_b = normalized_data(u, v)
    second_a = FormalSeries.seed(field, coefficients[2])
    second_b = FormalSeries.seed(field, coefficients[3])
    return [c1, c2, c3, c4, first_a, first_b, second_a, second_b]


def transverse_jacobian() -> list[list[Fraction]]:
    """Return d(E2,E1,E0)/d(v,a,b) at the labelled base point."""

    field = Field(series_order=1)
    base = base_coordinates()
    columns = []
    for variable in range(1, 4):
        coefficients = [[value, Q(0)] for value in (BASE_U, BASE_V, *SECOND_SECTION)]
        coefficients[variable][1] = Q(1)
        values = residuals_from_jets(series_coordinates(coefficients, field))
        columns.append([values[index].coefficients[1] for index in SECOND_EQUATIONS])
    return [[columns[column][row] for column in range(3)] for row in range(3)]


def formal_branch(order: int = 8, pade_degree: int = 5) -> dict[str, object]:
    """Lift the smooth local branch through the Fermigier seed to ``order``."""

    if order < 1:
        raise ValueError("formal order must be positive")
    field = Field(series_order=order)
    coefficients = [[Q(0)] * (order + 1) for _ in range(4)]
    coefficients[0][0], coefficients[0][1] = BASE_U, Q(1)
    coefficients[1][0] = BASE_V
    coefficients[2][0], coefficients[3][0] = SECOND_SECTION
    jacobian = transverse_jacobian()
    determinant = (
        jacobian[0][0] * (jacobian[1][1] * jacobian[2][2] - jacobian[1][2] * jacobian[2][1])
        - jacobian[0][1] * (jacobian[1][0] * jacobian[2][2] - jacobian[1][2] * jacobian[2][0])
        + jacobian[0][2] * (jacobian[1][0] * jacobian[2][1] - jacobian[1][1] * jacobian[2][0])
    )
    if not determinant:
        raise AssertionError("the Fermigier transverse Jacobian is singular")

    for degree in range(1, order + 1):
        values = residuals_from_jets(series_coordinates(coefficients, field))
        correction = solve_square_over_q(
            jacobian, [-values[index].coefficients[degree] for index in SECOND_EQUATIONS]
        )
        for variable, value in enumerate(correction, start=1):
            coefficients[variable][degree] = value
        checked = residuals_from_jets(series_coordinates(coefficients, field))
        if any(checked[index].coefficients[degree] for index in SECOND_EQUATIONS):
            raise AssertionError("formal transverse solve failed")
        if any(
            checked[index].coefficients[degree]
            for index in (0, 1, 2, 3)
        ):
            raise AssertionError("the Fermigier surface did not retain its first section")

    coordinate_series = {
        name: [str(value) for value in row]
        for name, row in zip(("u", "v", "second_intercept", "second_slope"), coefficients)
    }
    recognized = reconstructed_second_line(FormalSeries.seed(field, coefficients[0]))
    if any(
        expected.coefficients != tuple(actual)
        for expected, actual in zip(recognized, coefficients[1:])
    ):
        raise AssertionError("the reconstructed rational model disagrees with the formal lift")
    return {
        "base_parameters": {"u": str(BASE_U), "v": str(BASE_V)},
        "local_parameter": "u = -3 + t",
        "order": order,
        "transverse_variables": ["v", "second_intercept", "second_slope"],
        "transverse_jacobian_determinant": str(determinant),
        "coordinate_series": coordinate_series,
        "reconstructed_rational_model": {
            "v": "(u^2 + u + 2)/u",
            "second_intercept": "-(2u^6-u^5+4u^4-u^3-8u^2-4u-16)/(2(u+2)(u^2+2)(u^2-u+4))",
            "second_slope": "u/(u^2+2)",
            "matches_recursive_formal_lift_through_order": order,
            "not_yet_a_global_section_identity_certificate": True,
        },
        "recognized_low_degree_rational_series": small_pade_models(
            coordinate_series, max_degree=pade_degree
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=8)
    parser.add_argument("--pade-degree", type=int, default=5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = formal_branch(args.order, args.pade_degree)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
