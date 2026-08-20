#!/usr/bin/env python3
"""Bounded direct search for a 31st point on ICARM curve 273, rank >= 30.

The curve has 30 independently certified rational points.  This script
searches exact birational charts of the rank-30 fiber for a 31st direction.  A line of slope ``m`` through a known point cuts the curve in two
additional rational points exactly when an explicit quartic ``D_P(m)`` is a
square.  Affine slope charts determined by the published points can reveal
points hidden far outside a uniform x-height box.  We also search affine
x-charts through pairs of published abscissas and integer-offset charts around
each published abscissa.

All PARI ``hyperellratpoints`` calls have individual time caps.  Returned
points are checked with Fraction arithmetic.  PARI height pairings are used
only to propose subgroup relations, which are replayed exactly.  A point not
so resolved is tested for an exact rank-31 finite-reduction certificate.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path
import platform
import re
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

from icarm_curve273 import (
    A as COEFFICIENT_A,
    B as COEFFICIENT_B,
    GENERAL_WEIERSTRASS_COEFFICIENTS,
    POINTS as PUBLISHED_POINTS,
    on_curve as point_on_general_curve,
    short_coefficients as short_weierstrass_coefficients,
    to_short as to_short_point,
)
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)
from pari_bridge import pari_version
from search_extra_points import gp_rational, gp_vector, run_gp


Q = Fraction
Polynomial = tuple[Fraction, ...]  # coefficients from constant term upward
RationalPoint = tuple[Fraction, Fraction]
GroupPoint = RationalPoint | None
DEFAULT_OUTPUT = Path(
    "artifacts/local/elliptic-curves/curve273-rank31-search.json"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_icarm_curve273_rank31.py"
)


@dataclass(frozen=True)
class QuarticChart:
    identifier: str
    kind: str
    polynomial: Polynomial
    center: Fraction
    scale: Fraction
    base_index: int | None
    seed_parameters: tuple[Fraction, ...]
    height_bound: int


def poly_add(left: Polynomial, right: Polynomial) -> Polynomial:
    answer = [Q(0)] * max(len(left), len(right))
    for index, value in enumerate(left):
        answer[index] += value
    for index, value in enumerate(right):
        answer[index] += value
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def poly_scale(polynomial: Polynomial, scalar: Fraction) -> Polynomial:
    return tuple(Q(scalar) * value for value in polynomial)


def poly_multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    answer = [Q(0)] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] += left_value * right_value
    return tuple(answer)


def poly_power(polynomial: Polynomial, exponent: int) -> Polynomial:
    answer: Polynomial = (Q(1),)
    for _ in range(exponent):
        answer = poly_multiply(answer, polynomial)
    return answer


def poly_evaluate(polynomial: Polynomial, value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(polynomial):
        answer = answer * value + coefficient
    return answer


def affine_substitute(
    polynomial: Polynomial, center: Fraction, scale: Fraction
) -> Polynomial:
    linear = (Q(center), Q(scale))
    answer: Polynomial = (Q(0),)
    power: Polynomial = (Q(1),)
    for coefficient in polynomial:
        answer = poly_add(answer, poly_scale(power, coefficient))
        power = poly_multiply(power, linear)
    return answer


def slope_discriminant(point: RationalPoint) -> Polynomial:
    """Return the quartic discriminant for lines through ``point``."""

    x_point, y_point = point
    m: Polynomial = (Q(0), Q(1))
    n: Polynomial = (y_point, -x_point)
    c2 = poly_add(poly_multiply(m, m), m)  # a1=1 and a2=0
    alpha = poly_add((x_point,), poly_scale(c2, Q(-1)))
    c1 = poly_add(
        poly_add(poly_scale(poly_multiply(m, n), Q(2)), n),
        (Q(-COEFFICIENT_A),),
    )
    discriminant = poly_add(
        poly_add(
            poly_multiply(alpha, alpha),
            poly_scale(poly_multiply(alpha, (x_point,)), Q(-4)),
        ),
        poly_scale(c1, Q(4)),
    )
    if len(discriminant) != 5 or discriminant[-1] != 1:
        raise AssertionError("the slope discriminant is not a monic quartic")
    return discriminant


def slope_to_points(
    base_point: RationalPoint, slope: Fraction, square_root: Fraction
) -> tuple[RationalPoint, RationalPoint]:
    """Map a rational point on the slope quartic back to two curve points."""

    x_point, y_point = base_point
    slope = Q(slope)
    square_root = Q(square_root)
    c2 = slope * slope + slope
    alpha = x_point - c2
    first_x = (-alpha + square_root) / 2
    second_x = (-alpha - square_root) / 2
    intercept = y_point - slope * x_point
    points = (
        (first_x, slope * first_x + intercept),
        (second_x, slope * second_x + intercept),
    )
    if not all(point_on_general_curve(point) for point in points):
        raise AssertionError("a slope-quartic image missed the elliptic curve")
    return points


def x_polynomial() -> Polynomial:
    """Return ``(2y+x)^2=4x^3+x^2+4Ax+4B``."""

    return (
        Q(4 * COEFFICIENT_B),
        Q(4 * COEFFICIENT_A),
        Q(1),
        Q(4),
    )


def x_chart_to_points(x_value: Fraction, square_root: Fraction) -> tuple[RationalPoint, RationalPoint]:
    points = (
        (x_value, (-x_value + square_root) / 2),
        (x_value, (-x_value - square_root) / 2),
    )
    if not all(point_on_general_curve(point) for point in points):
        raise AssertionError("an x-chart image missed the elliptic curve")
    return points


def point_negate(point: GroupPoint) -> GroupPoint:
    if point is None:
        return None
    x_value, y_value = point
    return x_value, -y_value - x_value


def point_add(left: GroupPoint, right: GroupPoint) -> GroupPoint:
    """Exact group law on ``y^2+x*y=x^3+A*x+B``."""

    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2 and y1 + y2 + x1 == 0:
        return None
    if left == right:
        denominator = 2 * y1 + x1
        if denominator == 0:
            return None
        slope = (3 * x1 * x1 + COEFFICIENT_A - y1) / denominator
    else:
        slope = (y2 - y1) / (x2 - x1)
    intercept = y1 - slope * x1
    x3 = slope * slope + slope - x1 - x2
    y3 = -(slope + 1) * x3 - intercept
    answer = (x3, y3)
    if not point_on_general_curve(answer):
        raise AssertionError("the exact group law returned an off-curve point")
    return answer


def point_multiply(point: GroupPoint, scalar: int) -> GroupPoint:
    if scalar < 0:
        return point_multiply(point_negate(point), -scalar)
    answer: GroupPoint = None
    addend = point
    while scalar:
        if scalar & 1:
            answer = point_add(answer, addend)
        addend = point_add(addend, addend)
        scalar >>= 1
    return answer


def exact_linear_combination(coefficients: Sequence[int]) -> GroupPoint:
    if len(coefficients) != len(PUBLISHED_POINTS):
        raise ValueError("a relation vector must have length 29")
    answer: GroupPoint = None
    for coefficient, point in zip(coefficients, PUBLISHED_POINTS):
        answer = point_add(answer, point_multiply(point, int(coefficient)))
    return answer


def polynomial_gp(polynomial: Polynomial) -> str:
    return "+".join(
        f"{gp_rational(coefficient)}*x^{degree}"
        for degree, coefficient in enumerate(polynomial)
        if coefficient
    ) or "0"


def secant_slope(base_index: int, point_index: int) -> Fraction:
    base = PUBLISHED_POINTS[base_index]
    point = PUBLISHED_POINTS[point_index]
    return (point[1] - base[1]) / (point[0] - base[0])


def build_charts(
    *,
    x_pair_height: int,
    x_offset_height: int,
    slope_offset_height: int,
    slope_pair_height: int,
    slope_pair_count: int,
) -> tuple[QuarticChart, ...]:
    charts: list[QuarticChart] = []
    x_poly = x_polynomial()

    # Every unordered pair gives a chart sending two public abscissas to 0,1.
    for left, right in combinations(range(len(PUBLISHED_POINTS)), 2):
        center = PUBLISHED_POINTS[left][0]
        scale = PUBLISHED_POINTS[right][0] - center
        charts.append(
            QuarticChart(
                identifier=f"xpair_p{left + 1:02d}_p{right + 1:02d}",
                kind="x_pair_affine",
                polynomial=affine_substitute(x_poly, center, scale),
                center=center,
                scale=scale,
                base_index=None,
                seed_parameters=(Q(0), Q(1)),
                height_bound=x_pair_height,
            )
        )

    for index, point in enumerate(PUBLISHED_POINTS):
        center = point[0]
        charts.append(
            QuarticChart(
                identifier=f"xoffset_p{index + 1:02d}",
                kind="x_integer_offset",
                polynomial=affine_substitute(x_poly, center, Q(1)),
                center=center,
                scale=Q(1),
                base_index=None,
                seed_parameters=(Q(0),),
                height_bound=x_offset_height,
            )
        )

    # Offset every oriented secant slope by a small rational u.
    for base_index in range(len(PUBLISHED_POINTS)):
        discriminant = slope_discriminant(PUBLISHED_POINTS[base_index])
        for point_index in range(len(PUBLISHED_POINTS)):
            if point_index == base_index:
                continue
            center = secant_slope(base_index, point_index)
            charts.append(
                QuarticChart(
                    identifier=(
                        f"soffset_p{base_index + 1:02d}_q{point_index + 1:02d}"
                    ),
                    kind="slope_integer_offset",
                    polynomial=affine_substitute(discriminant, center, Q(1)),
                    center=center,
                    scale=Q(1),
                    base_index=base_index,
                    seed_parameters=(Q(0),),
                    height_bound=slope_offset_height,
                )
            )

    # Rank all slope-pair charts by how many other public slopes they compress
    # at successively smaller heights, then by coefficient bit size.  The
    # scoring is deterministic and uses only public points.
    ranked_pairs: list[tuple[Any, QuarticChart]] = []
    thresholds = (10**6, 10**9, 10**12, 10**15, 10**18)
    for base_index in range(len(PUBLISHED_POINTS)):
        other_indices = tuple(
            index for index in range(len(PUBLISHED_POINTS)) if index != base_index
        )
        slopes = {
            index: secant_slope(base_index, index) for index in other_indices
        }
        discriminant = slope_discriminant(PUBLISHED_POINTS[base_index])
        for first, second in combinations(other_indices, 2):
            center = slopes[first]
            scale = slopes[second] - center
            parameters = tuple((slope - center) / scale for slope in slopes.values())
            counts = tuple(
                sum(
                    max(abs(value.numerator), value.denominator) <= threshold
                    for value in parameters
                )
                for threshold in thresholds
            )
            transformed = affine_substitute(discriminant, center, scale)
            coefficient_bits = max(
                max(abs(value.numerator).bit_length(), value.denominator.bit_length())
                for value in transformed
            )
            chart = QuarticChart(
                identifier=(
                    f"spair_p{base_index + 1:02d}_q{first + 1:02d}_"
                    f"q{second + 1:02d}"
                ),
                kind="slope_pair_affine",
                polynomial=transformed,
                center=center,
                scale=scale,
                base_index=base_index,
                seed_parameters=(Q(0), Q(1)),
                height_bound=slope_pair_height,
            )
            score = tuple(-count for count in reversed(counts)) + (
                coefficient_bits,
                chart.identifier,
            )
            ranked_pairs.append((score, chart))
    ranked_pairs.sort(key=lambda item: item[0])
    charts.extend(chart for _, chart in ranked_pairs[:slope_pair_count])
    return tuple(charts)


POINT_PATTERN = re.compile(
    r"^POINT\s+(\S+)\s+(-?\d+(?:/\d+)?)\s+(-?\d+(?:/\d+)?)$",
    re.MULTILINE,
)


def search_chart(
    chart: QuarticChart, *, timeout: float, stack_bytes: int
) -> tuple[tuple[tuple[Fraction, Fraction], ...], int, float]:
    program = "\n".join(
        (
            f"P={polynomial_gp(chart.polynomial)};",
            "gettime();",
            f"R=hyperellratpoints(P,{chart.height_bound});",
            'print("PARI_MILLISECONDS ",gettime());',
            (
                f'for(i=1,#R,print("POINT {chart.identifier} ",'
                "R[i][1],\" \",R[i][2]));"
            ),
            "quit",
        )
    ) + "\n"
    output, wall_seconds = run_gp(
        program, timeout=timeout, stack_bytes=stack_bytes
    )
    milliseconds_match = re.search(r"PARI_MILLISECONDS (\d+)", output)
    if milliseconds_match is None:
        raise AssertionError("PARI omitted its chart timing marker")
    points = tuple(
        (Q(parameter), Q(square_root))
        for identifier, parameter, square_root in POINT_PATTERN.findall(output)
        if identifier == chart.identifier
    )
    return points, int(milliseconds_match.group(1)), wall_seconds


def chart_images(
    chart: QuarticChart, quartic_points: Iterable[tuple[Fraction, Fraction]]
) -> tuple[RationalPoint, ...]:
    images: list[RationalPoint] = []
    for parameter, square_root in quartic_points:
        if poly_evaluate(chart.polynomial, parameter) != square_root * square_root:
            raise AssertionError("PARI returned a point off a search quartic")
        if parameter in chart.seed_parameters:
            continue
        coordinate = chart.center + chart.scale * parameter
        if chart.base_index is None:
            mapped = x_chart_to_points(coordinate, square_root)
        else:
            mapped = slope_to_points(
                PUBLISHED_POINTS[chart.base_index], coordinate, square_root
            )
        images.extend(mapped)
    return tuple(images)


def discover_relation(
    point: RationalPoint, *, timeout: float, stack_bytes: int
) -> tuple[int, ...] | None:
    curve = ",".join(gp_rational(value) for value in GENERAL_WEIERSTRASS_COEFFICIENTS)
    basis = ",".join(gp_vector(value) for value in PUBLISHED_POINTS)
    program = "\n".join(
        (
            "default(realprecision,140);",
            f"E=ellinit([{curve}]);",
            f"B=[{basis}];",
            "H=ellheightmatrix(E,B);",
            f"Q={gp_vector(point)};",
            "V=vector(#B,j,ellheight(E,B[j],Q))~;",
            "C=round(matsolve(H,V));",
            "S=[0];for(j=1,#B,S=elladd(E,S,ellmul(E,B[j],C[j])));",
            'print("RELATION ",Vec(C)," EXACT ",S==Q);',
            "quit",
        )
    ) + "\n"
    try:
        output, _ = run_gp(program, timeout=timeout, stack_bytes=stack_bytes)
    except (RuntimeError, subprocess.TimeoutExpired):
        return None
    match = re.search(r"^RELATION \[(.*?)\] EXACT ([01])$", output, re.MULTILINE)
    if match is None or match.group(2) != "1":
        return None
    coefficients = tuple(int(value.strip()) for value in match.group(1).split(","))
    if exact_linear_combination(coefficients) != point:
        raise AssertionError("PARI's proposed relation failed exact Python replay")
    return coefficients


def point_record(point: RationalPoint) -> dict[str, str]:
    return {"x": str(point[0]), "y": str(point[1])}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--x-pair-height", type=int, default=10_000)
    parser.add_argument("--x-offset-height", type=int, default=1_000_000)
    parser.add_argument("--slope-offset-height", type=int, default=10_000)
    parser.add_argument("--slope-pair-height", type=int, default=10_000)
    parser.add_argument("--slope-pair-count", type=int, default=400)
    parser.add_argument("--chart-timeout", type=float, default=3.0)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=500_000_000)
    parser.add_argument("--certificate-prime-bound", type=int, default=1000)
    args = parser.parse_args()
    if min(
        args.x_pair_height,
        args.x_offset_height,
        args.slope_offset_height,
        args.slope_pair_height,
        args.slope_pair_count,
    ) <= 0:
        raise SystemExit("all chart heights and counts must be positive")
    if args.chart_timeout <= 0 or args.relation_timeout <= 0:
        raise SystemExit("timeouts must be positive")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes must be at least 64MB")

    charts = build_charts(
        x_pair_height=args.x_pair_height,
        x_offset_height=args.x_offset_height,
        slope_offset_height=args.slope_offset_height,
        slope_pair_height=args.slope_pair_height,
        slope_pair_count=args.slope_pair_count,
    )
    print(
        f"R31|stage=start"
        f"|basis_rank=30"
        f"|basis_points={len(PUBLISHED_POINTS)}"
        f"|charts={len(charts)}"
        f"|x_pair_height={args.x_pair_height}"
        f"|x_offset_height={args.x_offset_height}"
        f"|slope_offset_height={args.slope_offset_height}"
        f"|slope_pair_height={args.slope_pair_height}"
        f"|slope_pair_count={args.slope_pair_count}",
        flush=True,
    )

    started = time.monotonic()
    completed: list[dict[str, Any]] = []
    timeouts: list[str] = []
    discovered: dict[RationalPoint, list[str]] = {}
    public_set = set(PUBLISHED_POINTS) | {
        point_negate(point) for point in PUBLISHED_POINTS
    }
    total_pari_milliseconds = 0
    for index, chart in enumerate(charts, 1):
        try:
            quartic_points, milliseconds, wall_seconds = search_chart(
                chart,
                timeout=args.chart_timeout,
                stack_bytes=args.stack_bytes,
            )
        except subprocess.TimeoutExpired:
            timeouts.append(chart.identifier)
            print(f"R31|timeout|chart={chart.identifier}", flush=True)
            continue
        images = chart_images(chart, quartic_points)
        unexpected = tuple(point for point in images if point not in public_set)
        for point in unexpected:
            discovered.setdefault(point, []).append(chart.identifier)
        total_pari_milliseconds += milliseconds
        completed.append(
            {
                "identifier": chart.identifier,
                "kind": chart.kind,
                "height_bound": chart.height_bound,
                "quartic_affine_point_count_including_signs": len(quartic_points),
                "nonseed_mapped_image_count_including_duplicates": len(images),
                "pari_milliseconds": milliseconds,
                "wall_seconds": wall_seconds,
            }
        )
        if index % 100 == 0 or unexpected:
            print(
                f"R31|progress|charts={index}/{len(charts)}"
                f"|timeouts={len(timeouts)}"
                f"|unexpected_images={len(discovered)}",
                flush=True,
            )

    candidate_records: list[dict[str, Any]] = []
    target_hit = False
    for point, sources in sorted(
        discovered.items(), key=lambda item: (item[0][0], item[0][1])
    ):
        relation = discover_relation(
            point, timeout=args.relation_timeout, stack_bytes=args.stack_bytes
        )
        record: dict[str, Any] = {
            **point_record(point),
            "source_charts": sorted(set(sources)),
            "exact_curve_membership_checked": True,
        }
        if relation is not None:
            record.update(
                {
                    "classification": "exactly_in_certified_rank30_subgroup",
                    "rank30_basis_relation": list(relation),
                    "exact_fraction_group_law_replay": True,
                }
            )
        else:
            augmented = tuple(to_short_point(value) for value in PUBLISHED_POINTS) + (
                to_short_point(point),
            )
            signatures = find_mod2_reduction_certificate(
                short_weierstrass_coefficients(),
                augmented,
                prime_bound=args.certificate_prime_bound,
            )
            binary_rank = combined_mod2_rank(signatures, len(augmented))
            record.update(
                {
                    "classification": (
                        "exact_independent_31st_point"
                        if binary_rank == 31
                        else "unresolved_after_rank30_relation_and_mod2_search"
                    ),
                    "augmented_mod2_rank": binary_rank,
                    "certificate_primes": [signature.prime for signature in signatures],
                    "certificate_prime_bound": args.certificate_prime_bound,
                }
            )
            target_hit |= binary_rank == 31
        candidate_records.append(record)

        print(
            f"R31|candidate"
            f"|classification={record['classification']}"
            f"|augmented_rank={record.get('augmented_mod2_rank', 'subgroup')}"
            f"|x={record['x']}"
            f"|y={record['y']}",
            flush=True,
        )

    chart_payload = [
        {
            "identifier": chart.identifier,
            "kind": chart.kind,
            "height_bound": chart.height_bound,
            "base_index": None if chart.base_index is None else chart.base_index + 1,
            "center": str(chart.center),
            "scale": str(chart.scale),
            "seed_parameters": [str(value) for value in chart.seed_parameters],
            "quartic_coefficients_constant_first": [
                str(value) for value in chart.polynomial
            ],
        }
        for chart in charts
    ]
    chart_digest = hashlib.sha256(
        json.dumps(chart_payload, separators=(",", ":")).encode()
    ).hexdigest()
    artifact = {
        "schema_version": 1,
        "artifact_kind": "bounded_elliptic_curve_rank31_point_search",
        "status": (
            "exact_rank31_target_hit"
            if target_hit
            else "bounded_search_no_certified_31st_point"
        ),
        "claim_scope": {
            "exact": (
                "chart construction, curve membership, subgroup relation replays, "
                "and any finite-reduction independence certificate"
            ),
            "bounded": (
                "PARI hyperellratpoints enumeration in exactly the completed "
                "quartic chart boxes; no rank upper bound"
            ),
            "search_basis": (
                "The search operates on the 30 independently certified rational "
                "points of ICARM curve 273 and looks for a 31st direction."
            ),
        },
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "parameters": {
            "x_pair_height": args.x_pair_height,
            "x_offset_height": args.x_offset_height,
            "slope_offset_height": args.slope_offset_height,
            "slope_pair_height": args.slope_pair_height,
            "slope_pair_count": args.slope_pair_count,
            "chart_timeout_seconds_each": args.chart_timeout,
            "relation_timeout_seconds_each": args.relation_timeout,
            "stack_bytes_each": args.stack_bytes,
            "certificate_prime_bound": args.certificate_prime_bound,
        },
        "chart_manifest": {
            "count": len(charts),
            "sha256": chart_digest,
            "kind_counts": {
                kind: sum(chart.kind == kind for chart in charts)
                for kind in sorted({chart.kind for chart in charts})
            },
            "manifest_stored_inline": False,
            "regeneration": "deterministic build_charts in the pinned script",
        },
        "search_result": {
            "completed_chart_count": len(completed),
            "timed_out_chart_count": len(timeouts),
            "timed_out_chart_identifiers": timeouts,
            "total_pari_milliseconds_completed": total_pari_milliseconds,
            "wall_seconds": time.monotonic() - started,
            "unique_nonpublic_images": len(discovered),
            "candidate_records": candidate_records,
            "certified_independent_31st_point_count": sum(
                record["classification"] == "exact_independent_31st_point"
                for record in candidate_records
            ),
            "rank31_target_hit": target_hit,
            "completed_chart_summaries": completed,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")
    print(
        f"R31|stage=done"
        f"|completed={len(completed)}/{len(charts)}"
        f"|timeouts={len(timeouts)}"
        f"|unique_nonbasis={len(discovered)}"
        f"|rank31_hit={str(target_hit).lower()}",
        flush=True,
    )


if __name__ == "__main__":
    main()
