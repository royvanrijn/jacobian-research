#!/usr/bin/env python3
"""Numerically triage rank on the exact Fermigier two-section component.

The component is kept in its compact rational-root and affine-line form.  At
each declared rational ``u`` and Mestre parameter ``T`` this program maps the
twelve visible points and the two selected affine points to the short
Jacobian, then asks PARI only for the numerical rank of their canonical-height
matrix.  This is an efficient way to locate a possible specialization for a
subsequent exact finite-reduction certificate; it is not a rank proof.

No universal two-section residual is expanded or materialized.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import gcd
import json
from pathlib import Path
import subprocess

from icarm_curve245_mestre import fermigier_roots
from mestre_root_tuples import SixRootMestreConstruction
from nagao_1994 import (
    primitive_quartic_coefficients,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    quartic_value,
    short_jacobian_coefficients,
)
from probe_mestre_fermigier_two_section_local_continuation import (
    normalized_data,
    reconstructed_second_line,
)
from screen_mestre_fermigier_two_section_escape import square_root


Q = Fraction
DEFAULT_HEIGHT = 10
DEFAULT_T_VALUES = (Q(1), Q(2), Q(3))


def rational_text(value: Fraction) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def gp_rational(value: Fraction) -> str:
    rendered = rational_text(value)
    return rendered if "/" not in rendered else f"({rendered})"


def rational_parameters(height: int) -> tuple[Fraction, ...]:
    """Enumerate reduced component parameters with numerator/denominator bound."""

    parameters: list[Fraction] = []
    for denominator in range(1, height + 1):
        for numerator in range(-height, height + 1):
            if gcd(numerator, denominator) != 1:
                continue
            parameter = Q(numerator, denominator)
            if parameter in (0, -2) or parameter in parameters:
                continue
            parameters.append(parameter)
    return tuple(parameters)


def specialized_points(
    u: Fraction, parameter: Fraction,
) -> tuple[tuple[Fraction, ...], tuple[tuple[Fraction, Fraction], ...]]:
    """Return the short model and all fourteen supplied points exactly."""

    v, second_intercept, second_slope = reconstructed_second_line(u)
    first_intercept, first_slope = normalized_data(u, v)[4:6]
    source_roots = fermigier_roots(u, v)
    roots = tuple(
        (root - source_roots[0]) / (source_roots[1] - source_roots[0])
        for root in source_roots
    )
    if len(set(roots)) != 6:
        raise ValueError("the specialization left the split-six-root open locus")
    construction = SixRootMestreConstruction(roots)
    quartic = primitive_quartic_coefficients(construction, parameter)
    points = [
        quartic_point_to_short_jacobian(construction, parameter, point)
        for point in primitive_visible_points(construction, parameter)
    ]
    for intercept, slope in (
        (first_intercept, first_slope),
        (second_intercept, second_slope),
    ):
        x_value = intercept + slope * parameter
        points.append(
            quartic_point_to_short_jacobian(
                construction,
                parameter,
                (x_value, square_root(quartic_value(quartic, x_value))),
            )
        )
    return short_jacobian_coefficients(construction, parameter), tuple(points)


def numerical_height_rank(
    coefficients: tuple[Fraction, ...],
    points: tuple[tuple[Fraction, Fraction], ...],
    *, timeout: float,
) -> int:
    """Ask PARI for the numerical rank of the canonical-height matrix."""

    curve = ",".join(gp_rational(value) for value in coefficients)
    group = ",".join(
        f"[{gp_rational(x_value)},{gp_rational(y_value)}]"
        for x_value, y_value in points
    )
    program = (
        f"E=ellinit([{curve}]);P=[{group}];"
        'print("HEIGHT_RANK_BEGIN");print(matrank(ellheightmatrix(E,P)));'
        'print("HEIGHT_RANK_END");quit\n'
    )
    completed = subprocess.run(
        ("gp", "-q"),
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    begin = lines.index("HEIGHT_RANK_BEGIN")
    if lines[begin + 2] != "HEIGHT_RANK_END":
        raise AssertionError("PARI height-rank delimiters changed")
    return int(lines[begin + 1])


def replay(
    *, height: int = DEFAULT_HEIGHT,
    parameters: tuple[Fraction, ...] = DEFAULT_T_VALUES,
    timeout: float = 20.0,
) -> dict[str, object]:
    if height < 1 or not parameters:
        raise ValueError("height and parameter list must be nonempty")
    records = []
    for u in rational_parameters(height):
        for parameter in parameters:
            try:
                coefficients, points = specialized_points(u, parameter)
                rank = numerical_height_rank(coefficients, points, timeout=timeout)
            except (ArithmeticError, subprocess.SubprocessError, ValueError):
                continue
            records.append(
                {"u": rational_text(u), "T": rational_text(parameter), "numerical_height_rank": rank}
            )
    if not records:
        raise AssertionError("the declared height panel had no smooth specializations")
    best_rank = max(record["numerical_height_rank"] for record in records)
    leaders = [
        {"u": record["u"], "T": record["T"]}
        for record in records
        if record["numerical_height_rank"] == best_rank
    ]
    return {
        "status": "bounded numerical height-rank triage completed",
        "component": "v=(u^2+u+2)/u with the two exact affine sections",
        "u_height": height,
        "T_values": [rational_text(value) for value in parameters],
        "smooth_specialization_count": len(records),
        "point_count_per_specialization": 14,
        "best_numerical_height_rank": best_rank,
        "best_specializations": leaders,
        "records": records,
        "conclusion": "the declared panel supplies no numerical height-rank-14 specialization for an immediate finite-reduction independence certificate",
        "not_established": [
            "an algebraic rank upper bound at any specialization",
            "a generic relation or generic independence for the second section",
            "saturation, intersections, or a Shioda Gram matrix",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=DEFAULT_HEIGHT)
    parser.add_argument("--t-values", default="1,2,3")
    parser.add_argument("--timeout", type=float, default=20.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    parameters = tuple(Q(value) for value in args.t_values.split(",") if value)
    rendered = json.dumps(
        replay(height=args.height, parameters=parameters, timeout=args.timeout),
        indent=2,
        sort_keys=True,
    ) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
