#!/usr/bin/env python3
# <!-- status-consumer: EC-K3-ELKIES-2026-RESIDUAL-SELMER-GATE 7f8dffe58168acc8 -->
"""Search q12o5867 specializations through section-normalized slope quartics.

For a short Weierstrass curve ``y^2 = x^3 + A*x + B`` and a rational point
``P=(x0,y0)``, a line of slope ``m`` through ``P`` meets the curve at two
additional rational points precisely when

    w^2 = m^4 - 6*x0*m^2 + 8*y0*m - 3*x0^2 - 4*A.

Three secant slopes from ``P`` to other members of the certified generic
rank-17 basis are used to put known rational points of this quartic at
``z=0,1,infinity``. PARI ``hyperellratpoints`` or an explicitly supplied
ratpoints-compatible CPU/CUDA executable then searches a declared
projective-height box in ``z``. Every hit is mapped back and checked by exact
rational substitution. The common finite-reduction machinery measures its
image modulo the specialized generic subgroup.

This is a bounded point search, not a rank upper bound.  Every input must have
a completed unconditional residual-Selmer gate for its identical parameter
and minimal model.  A positive quotient escape is exact; failure to find or
detect an escape proves nothing.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import gcd, lcm
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence


REPOSITORY = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = REPOSITORY / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
sys.path.insert(0, str(ELLIPTIC_ROOT))
sys.path.insert(0, str(CAS))

from elliptic_candidate_record import (  # noqa: E402
    build_finite_quotient_certificate,
    is_on_weierstrass_curve,
    verify_finite_quotient_certificate,
)
from finite_quotient_escape import QuotientBlock, analyze_escape  # noqa: E402
from pari_bridge import pari_version  # noqa: E402
from elkies_residual_selmer_gate import require_gate_for_specialization  # noqa: E402
from ecsearch.rank_certification import add_rational_points  # noqa: E402
from ecsearch.q12o5867_point_search import (  # noqa: E402
    integral_square_scaled_coefficients,
    parse_ratpoints_abscissae,
    rational_square_root,
)


Q = Fraction
Point = tuple[Fraction, Fraction]
Matrix = tuple[int, int, int, int]
EXPECTED_STATUS = "PASS_EXACT_Q12O5867_SPECIALIZED_GENERIC_RANK17_LOWER_BOUND"


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def fraction_text(value: Fraction | int) -> str:
    value = Q(value)
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


def point_record(point: Point) -> list[str]:
    return [fraction_text(point[0]), fraction_text(point[1])]


def parse_point(record: Sequence[str]) -> Point:
    if len(record) != 2:
        raise ValueError("a point record must have two coordinates")
    return Q(record[0]), Q(record[1])


def parse_int_tuple(text: str) -> tuple[int, ...]:
    try:
        answer = tuple(int(value) for value in text.split(",") if value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from error
    if not answer:
        raise argparse.ArgumentTypeError("the integer list must be nonempty")
    return answer


def normalize_integer_matrix(values: Sequence[Fraction]) -> Matrix:
    denominator = lcm(*(value.denominator for value in values))
    integers = [int(value * denominator) for value in values]
    common = 0
    for value in integers:
        common = gcd(common, abs(value))
    if common:
        integers = [value // common for value in integers]
    first_nonzero = next(value for value in integers if value)
    if first_nonzero < 0:
        integers = [-value for value in integers]
    a, b, c, d = integers
    if a * d == b * c:
        raise ValueError("a slope chart must have nonzero determinant")
    return a, b, c, d


def slope_chart_matrix(s0: Fraction, s1: Fraction, sinfinity: Fraction) -> Matrix:
    """Return ``m=(a*z+b)/(c*z+d)`` sending 0,1,infinity as declared."""

    if len({s0, s1, sinfinity}) != 3:
        raise ValueError("the three slope anchors must be distinct")
    c = (s1 - s0) / (sinfinity - s1)
    a = sinfinity * c
    b = s0
    d = Q(1)
    matrix = normalize_integer_matrix((a, b, c, d))
    if mobius_value(matrix, Q(0)) != s0:
        raise AssertionError("the chart missed its z=0 anchor")
    if mobius_value(matrix, Q(1)) != s1:
        raise AssertionError("the chart missed its z=1 anchor")
    a_i, _b_i, c_i, _d_i = matrix
    if Q(a_i, c_i) != sinfinity:
        raise AssertionError("the chart missed its z=infinity anchor")
    return matrix


def mobius_value(matrix: Matrix, z_value: Fraction) -> Fraction:
    a, b, c, d = matrix
    denominator = c * z_value + d
    if denominator == 0:
        raise ZeroDivisionError("the chart value is infinite")
    return Q(a * z_value + b, denominator)


def inverse_mobius_value(matrix: Matrix, slope: Fraction) -> Fraction | None:
    a, b, c, d = matrix
    denominator = slope * c - a
    if denominator == 0:
        return None
    return Q(b - slope * d, denominator)


def rational_height_bits(value: Fraction) -> int:
    return max(abs(value.numerator).bit_length(), value.denominator.bit_length())


def in_search_box(value: Fraction, numerator_bound: int, denominator_bound: int) -> bool:
    return (
        abs(value.numerator) <= numerator_bound
        and value.denominator <= denominator_bound
    )


def section_slopes(points: Sequence[Point], base_index: int) -> dict[int, Fraction]:
    x0, y0 = points[base_index]
    answer = {}
    for index, (x_value, y_value) in enumerate(points):
        if index == base_index or x_value == x0:
            continue
        answer[index] = (y_value - y0) / (x_value - x0)
    return answer


def choose_charts(
    points: Sequence[Point],
    *,
    chart_count: int,
    charts_per_base: int,
    anchor_pool_size: int,
    numerator_bound: int,
    denominator_bound: int,
) -> tuple[dict[str, Any], ...]:
    """Choose deterministic charts that expose many known section slopes."""

    scored: list[tuple[tuple[Any, ...], dict[str, Any]]] = []
    for base_index in range(len(points)):
        slopes = section_slopes(points, base_index)
        pool = tuple(
            index
            for index, _slope in sorted(
                slopes.items(),
                key=lambda item: (rational_height_bits(item[1]), item[0]),
            )[:anchor_pool_size]
        )
        for triple in combinations(pool, 3):
            for ordered in permutations(triple):
                s0, s1, sinfinity = (slopes[index] for index in ordered)
                if len({s0, s1, sinfinity}) != 3:
                    continue
                matrix = slope_chart_matrix(s0, s1, sinfinity)
                visible: list[tuple[int, Fraction]] = []
                finite_coordinates: list[tuple[int, Fraction]] = []
                for section_index, slope in slopes.items():
                    coordinate = inverse_mobius_value(matrix, slope)
                    if coordinate is None:
                        continue
                    finite_coordinates.append((section_index, coordinate))
                    if in_search_box(coordinate, numerator_bound, denominator_bound):
                        visible.append((section_index, coordinate))
                bit_sizes = sorted(
                    rational_height_bits(coordinate)
                    for _index, coordinate in finite_coordinates
                )
                quality = (
                    -len(visible),
                    sum(bit_sizes),
                    max(bit_sizes, default=0),
                    base_index,
                    ordered,
                )
                scored.append(
                    (
                        quality,
                        {
                            "base_section_index": base_index,
                            "anchor_section_indices_z0_z1_zinfinity": list(ordered),
                            "matrix_a_b_c_d": list(matrix),
                            "known_finite_slope_coordinates": [
                                {
                                    "section_index": index,
                                    "z": fraction_text(coordinate),
                                    "inside_search_box": in_search_box(
                                        coordinate,
                                        numerator_bound,
                                        denominator_bound,
                                    ),
                                }
                                for index, coordinate in sorted(finite_coordinates)
                            ],
                            "known_finite_slopes_inside_search_box": len(visible),
                            "quality_key": [
                                -quality[0],
                                quality[1],
                                quality[2],
                                quality[3],
                                list(quality[4]),
                            ],
                        },
                    )
                )
    scored.sort(key=lambda item: item[0])
    selected: list[dict[str, Any]] = []
    base_counts: dict[int, int] = {}
    used_matrices: set[tuple[int, Matrix]] = set()
    used_anchor_sets: set[tuple[int, tuple[int, int, int]]] = set()
    for _quality, record in scored:
        base_index = int(record["base_section_index"])
        matrix = tuple(int(value) for value in record["matrix_a_b_c_d"])
        matrix_key = base_index, matrix  # type: ignore[assignment]
        anchor_key = (
            base_index,
            tuple(
                sorted(
                    int(value)
                    for value in record[
                        "anchor_section_indices_z0_z1_zinfinity"
                    ]
                )
            ),
        )
        if base_counts.get(base_index, 0) >= charts_per_base:
            continue
        if matrix_key in used_matrices:
            continue
        if anchor_key in used_anchor_sets:
            continue
        selected.append(record)
        base_counts[base_index] = base_counts.get(base_index, 0) + 1
        used_matrices.add(matrix_key)
        used_anchor_sets.add(anchor_key)
        if len(selected) == chart_count:
            break
    if len(selected) != chart_count:
        raise RuntimeError("not enough distinct-base slope charts were available")
    return tuple(selected)


def polynomial_add(left: Sequence[Fraction], right: Sequence[Fraction]) -> list[Fraction]:
    size = max(len(left), len(right))
    answer = [Q(0) for _ in range(size)]
    for index, value in enumerate(left):
        answer[index] += value
    for index, value in enumerate(right):
        answer[index] += value
    return answer


def polynomial_multiply(
    left: Sequence[Fraction], right: Sequence[Fraction]
) -> list[Fraction]:
    answer = [Q(0) for _ in range(len(left) + len(right) - 1)]
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] += left_value * right_value
    return answer


def polynomial_power(base: Sequence[Fraction], exponent: int) -> list[Fraction]:
    answer = [Q(1)]
    for _ in range(exponent):
        answer = polynomial_multiply(answer, base)
    return answer


def polynomial_scale(values: Sequence[Fraction], scale: Fraction) -> list[Fraction]:
    return [scale * value for value in values]


def slope_quartic_coefficients(
    coefficient_a: Fraction, base_point: Point, matrix: Matrix
) -> tuple[Fraction, ...]:
    """Return low-to-high coefficients of ``(cz+d)^4 D((az+b)/(cz+d))``."""

    x0, y0 = base_point
    a, b, c, d = matrix
    numerator = [Q(b), Q(a)]
    denominator = [Q(d), Q(c)]
    constant = -3 * x0**2 - 4 * coefficient_a
    terms = (
        polynomial_power(numerator, 4),
        polynomial_scale(
            polynomial_multiply(
                polynomial_power(numerator, 2), polynomial_power(denominator, 2)
            ),
            -6 * x0,
        ),
        polynomial_scale(
            polynomial_multiply(numerator, polynomial_power(denominator, 3)),
            8 * y0,
        ),
        polynomial_scale(polynomial_power(denominator, 4), constant),
    )
    answer: list[Fraction] = []
    for term in terms:
        answer = polynomial_add(answer, term)
    answer += [Q(0)] * (5 - len(answer))
    return tuple(answer[:5])


def polynomial_value(coefficients: Sequence[Fraction], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
    return answer


def gp_polynomial(coefficients: Sequence[Fraction]) -> str:
    return "+".join(
        f"({fraction_text(coefficient)})*x^{index}"
        for index, coefficient in enumerate(coefficients)
    )


def run_quartic_search(
    coefficients: Sequence[Fraction],
    *,
    numerator_bound: int,
    denominator_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[Point, ...], dict[str, Any]]:
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    height = f"[{numerator_bound},{denominator_bound}]"
    program = "\n".join(
        (
            f"Q={gp_polynomial(coefficients)};gettime();",
            f"R=hyperellratpoints(Q,{height});",
            'print("PARI_MILLISECONDS|",gettime());',
            'print("POINTS_BEGIN");print(R);print("POINTS_END");',
            "quit",
        )
    ) + "\n"
    started = time.monotonic()
    try:
        completed = subprocess.run(
            [executable, "-q", "-f", "-s", str(stack_bytes)],
            input=program,
            text=True,
            capture_output=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as error:
        return (), {
            "status": "timeout",
            "timeout_seconds": timeout,
            "wall_seconds": time.monotonic() - started,
            "stdout_tail": (error.stdout or "")[-1000:],
            "stderr_tail": (error.stderr or "")[-1000:],
            "program_sha256": sha256(program.encode()).hexdigest(),
        }
    wall_seconds = time.monotonic() - started
    fatal = [
        line for line in completed.stderr.splitlines() if "***" in line
    ]
    if completed.returncode != 0 or fatal:
        return (), {
            "status": "pari_error",
            "returncode": completed.returncode,
            "wall_seconds": wall_seconds,
            "error": " ".join(fatal or completed.stderr.splitlines())[:2000],
            "program_sha256": sha256(program.encode()).hexdigest(),
        }
    marker = re.search(
        r"POINTS_BEGIN\s*(.*?)\s*POINTS_END", completed.stdout, re.DOTALL
    )
    milliseconds = re.search(
        r"^PARI_MILLISECONDS\|(\d+)$", completed.stdout, re.MULTILINE
    )
    if marker is None or milliseconds is None:
        raise RuntimeError("PARI omitted slope-quartic output markers")
    points = tuple(
        (Q(x_value), Q(y_value))
        for x_value, y_value in re.findall(
            r"\[(-?\d+(?:/\d+)?),\s*(-?\d+(?:/\d+)?)\]", marker.group(1)
        )
    )
    if any(y_value**2 != polynomial_value(coefficients, x_value) for x_value, y_value in points):
        raise AssertionError("PARI returned a point off the exact slope quartic")
    return points, {
        "status": "completed",
        "height_specification": height,
        "wall_seconds": wall_seconds,
        "pari_milliseconds": int(milliseconds.group(1)),
        "signed_quartic_point_count": len(points),
        "program_sha256": sha256(program.encode()).hexdigest(),
    }


def quartic_points_from_abscissae(
    coefficients: Sequence[Fraction], abscissae: Sequence[Fraction]
) -> tuple[Point, ...]:
    points = []
    for x_value in abscissae:
        root = rational_square_root(polynomial_value(coefficients, x_value))
        if root is None:
            raise AssertionError("ratpoints returned an abscissa off the exact quartic")
        points.append((x_value, root))
        if root:
            points.append((x_value, -root))
    return tuple(points)


def run_ratpoints_quartic_search(
    coefficients: Sequence[Fraction],
    *,
    numerator_bound: int,
    denominator_bound: int,
    timeout: float,
    executable: Path,
    library: Path | None,
) -> tuple[tuple[Point, ...], dict[str, Any]]:
    """Search one normalized quartic with CPU or CUDA ratpoints exactly."""

    integral, ordinate_scale = integral_square_scaled_coefficients(coefficients)
    command = [
        str(executable.resolve()),
        " ".join(str(value) for value in integral),
        str(numerator_bound),
        "-du",
        str(denominator_bound),
        "-q",
        "-y",
    ]
    environment = os.environ.copy()
    if library is not None:
        old_library = environment.get("LD_LIBRARY_PATH")
        environment["LD_LIBRARY_PATH"] = (
            str(library.resolve())
            if not old_library
            else str(library.resolve()) + os.pathsep + old_library
        )
    started = time.monotonic()
    try:
        completed = subprocess.run(
            command,
            text=True,
            capture_output=True,
            env=environment,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as error:
        return (), {
            "status": "timeout",
            "timeout_seconds": timeout,
            "wall_seconds": time.monotonic() - started,
            "stdout_tail": (error.stdout or "")[-1000:],
            "stderr_tail": (error.stderr or "")[-1000:],
            "command": command,
        }
    wall_seconds = time.monotonic() - started
    if completed.returncode != 0:
        return (), {
            "status": "ratpoints_error",
            "returncode": completed.returncode,
            "wall_seconds": wall_seconds,
            "stderr_tail": completed.stderr[-2000:],
            "command": command,
        }
    abscissae = parse_ratpoints_abscissae(completed.stdout)
    points = quartic_points_from_abscissae(coefficients, abscissae)
    return points, {
        "status": "completed",
        "engine": "ratpoints-compatible executable",
        "command": command,
        "wall_seconds": wall_seconds,
        "integral_coefficients_low_to_high": [str(value) for value in integral],
        "ordinate_square_clearing_scale": str(ordinate_scale),
        "abscissa_count": len(abscissae),
        "signed_quartic_point_count": len(points),
        "stdout_sha256": sha256(completed.stdout.encode()).hexdigest(),
        "stderr_tail": completed.stderr[-2000:],
    }


def curve_points_from_quartic_point(
    coefficient_a: Fraction,
    coefficient_b: Fraction,
    base_point: Point,
    matrix: Matrix,
    quartic_point: Point,
) -> tuple[Point, Point] | None:
    z_value, transformed_ordinate = quartic_point
    a, b, c, d = matrix
    denominator = c * z_value + d
    if denominator == 0:
        return None
    slope = Q(a * z_value + b, denominator)
    ordinate = transformed_ordinate / denominator**2
    x0, y0 = base_point
    x_values = (
        (slope**2 - x0 + ordinate) / 2,
        (slope**2 - x0 - ordinate) / 2,
    )
    points = tuple(
        (x_value, slope * (x_value - x0) + y0) for x_value in x_values
    )
    model = (Q(0), Q(0), Q(0), coefficient_a, coefficient_b)
    if any(not is_on_weierstrass_curve(model, point) for point in points):
        raise AssertionError("a slope-quartic point failed exact curve transport")
    return points  # type: ignore[return-value]


def canonical_short_point(point: Point) -> Point:
    negative = (point[0], -point[1])

    def key(current: Point) -> tuple[tuple[int, int], tuple[int, int]]:
        return tuple(
            (coordinate.numerator, coordinate.denominator) for coordinate in current
        )  # type: ignore[return-value]

    return min((point, negative), key=key)


def signed_pair_relation_lookup(
    model: Sequence[Fraction], baseline: Sequence[Point]
) -> dict[Point, dict[str, int]]:
    """Recognize exact signed sums of two listed baseline sections up to sign."""

    lookup: dict[Point, dict[str, int]] = {}
    for left_index, left in enumerate(baseline):
        for right_index, right in enumerate(baseline):
            for left_sign in (1, -1):
                signed_left = left if left_sign == 1 else (left[0], -left[1])
                for right_sign in (1, -1):
                    signed_right = right if right_sign == 1 else (right[0], -right[1])
                    total = add_rational_points(model, signed_left, signed_right)
                    if total is None:
                        continue
                    lookup.setdefault(
                        canonical_short_point(total),
                        {
                            "left_section_index": left_index,
                            "left_sign": left_sign,
                            "right_section_index": right_index,
                            "right_sign": right_sign,
                        },
                    )
    return lookup


def exact_quotient_profiles(
    model: Sequence[Fraction],
    baseline: Sequence[Point],
    candidates: Sequence[Point],
    *,
    relation_primes: Sequence[int],
    reduction_prime_bound: int,
) -> dict[str, Any]:
    if not candidates:
        return {
            "candidate_count": 0,
            "maximum_marginal_dimension": 0,
            "independent_escape_basis_labels": [],
            "relation_prime_profiles": [],
        }
    all_points = (*baseline, *candidates)
    labels = tuple(f"candidate-{index:03d}" for index in range(len(candidates)))
    attempts = []
    maximum = 0
    basis_labels: list[str] = []
    for relation_prime in relation_primes:
        try:
            certificate = build_finite_quotient_certificate(
                model,
                all_points,
                relation_prime=relation_prime,
                prime_bound=reduction_prime_bound,
            )
            verify_finite_quotient_certificate(model, all_points, certificate)
            blocks = tuple(
                QuotientBlock.build(
                    modulus=relation_prime,
                    rows=signature["rows"],
                    column_count=len(all_points),
                    source=f"good-reduction-p={signature['prime']}",
                )
                for signature in certificate["signatures"]
            )
            profile = analyze_escape(
                blocks,
                known_column_count=len(baseline),
                candidate_labels=labels,
            )
            record = profile.to_record()
            record["finite_quotient_certificate"] = certificate
            attempts.append(record)
            if profile.marginal_dimension > maximum:
                maximum = profile.marginal_dimension
                basis_labels = list(profile.independent_escape_basis_labels)
        except (RuntimeError, ValueError, AssertionError) as error:
            attempts.append(
                {
                    "modulus": relation_prime,
                    "status": "bounded-certificate-error",
                    "exception": type(error).__name__,
                    "detail": str(error),
                }
            )
    return {
        "candidate_count": len(candidates),
        "candidate_labels": list(labels),
        "maximum_marginal_dimension": maximum,
        "independent_escape_basis_labels": basis_labels,
        "relation_prime_profiles": attempts,
        "promotion_threshold": 15,
        "promotion_eligible": maximum >= 15,
        "claim_boundary": (
            "positive finite-quotient escape is exact; bounded non-escape is not "
            "a dependence proof"
        ),
    }


def search_specialization(path: Path, args: argparse.Namespace) -> dict[str, Any]:
    raw = path.read_bytes()
    artifact = json.loads(raw)
    if artifact.get("status") != EXPECTED_STATUS:
        raise ValueError(f"{path} is not an exact-certified q12o5867 specialization")
    certificate = artifact["finite_quotient_independence"]
    model = tuple(Q(value) for value in certificate["certificate_short_model"])
    if model[:3] != (Q(0), Q(0), Q(0)):
        raise ValueError("the q12o5867 certificate model is no longer short")
    coefficient_a, coefficient_b = model[3], model[4]
    baseline = tuple(parse_point(record) for record in certificate["points"])
    if len(baseline) != 17:
        raise ValueError("the specialization does not contain 17 baseline sections")
    if any(not is_on_weierstrass_curve(model, point) for point in baseline):
        raise AssertionError("a baseline point missed the certificate model")

    charts = choose_charts(
        baseline,
        chart_count=args.charts_per_candidate,
        charts_per_base=args.charts_per_base,
        anchor_pool_size=args.anchor_pool_size,
        numerator_bound=args.numerator_bound,
        denominator_bound=args.denominator_bound,
    )
    known_x = {point[0] for point in baseline}
    discovered_by_x: dict[Fraction, dict[str, Any]] = {}
    chart_records = []
    for chart_index, chart in enumerate(charts):
        base_index = int(chart["base_section_index"])
        matrix = tuple(int(value) for value in chart["matrix_a_b_c_d"])
        coefficients = slope_quartic_coefficients(
            coefficient_a, baseline[base_index], matrix  # type: ignore[arg-type]
        )
        # The z=0 and z=1 anchors must remain exact rational points.
        for anchor in (Q(0), Q(1)):
            value = polynomial_value(coefficients, anchor)
            from math import isqrt

            numerator_root = isqrt(value.numerator)
            denominator_root = isqrt(value.denominator)
            if numerator_root**2 != value.numerator or denominator_root**2 != value.denominator:
                raise AssertionError("a normalized known slope is not on its quartic")
        if args.ratpoints is None:
            quartic_points, search_record = run_quartic_search(
                coefficients,
                numerator_bound=args.numerator_bound,
                denominator_bound=args.denominator_bound,
                timeout=args.chart_timeout,
                stack_bytes=args.stack_bytes,
            )
        else:
            quartic_points, search_record = run_ratpoints_quartic_search(
                coefficients,
                numerator_bound=args.numerator_bound,
                denominator_bound=args.denominator_bound,
                timeout=args.chart_timeout,
                executable=args.ratpoints,
                library=args.ratpoints_library,
            )
        transported_count = 0
        for quartic_point in quartic_points:
            pair = curve_points_from_quartic_point(
                coefficient_a,
                coefficient_b,
                baseline[base_index],
                matrix,  # type: ignore[arg-type]
                quartic_point,
            )
            if pair is None:
                continue
            transported_count += 2
            for point in pair:
                canonical = canonical_short_point(point)
                if canonical[0] in known_x:
                    continue
                discovered_by_x.setdefault(
                    canonical[0],
                    {
                        "point": canonical,
                        "first_chart_index": chart_index,
                        "source_quartic_point": point_record(quartic_point),
                    },
                )
        chart_records.append(
            {
                **chart,
                "quartic_coefficients_low_to_high": [
                    fraction_text(value) for value in coefficients
                ],
                "search": search_record,
                "exact_transported_curve_points_before_deduplication": transported_count,
            }
        )

    candidates = tuple(
        record["point"]
        for _x, record in sorted(
            discovered_by_x.items(),
            key=lambda item: (item[0].numerator, item[0].denominator),
        )
    )
    quotient = exact_quotient_profiles(
        model,
        baseline,
        candidates,
        relation_primes=args.relation_primes,
        reduction_prime_bound=args.reduction_prime_bound,
    )
    labels = quotient.get("candidate_labels", [])
    pair_relations = signed_pair_relation_lookup(model, baseline)
    candidate_records = []
    for label, point in zip(labels, candidates, strict=True):
        source = discovered_by_x[point[0]]
        candidate_records.append(
            {
                "label": label,
                "point_on_certificate_short_model": point_record(point),
                "exact_curve_membership_verified": True,
                "first_chart_index": source["first_chart_index"],
                "source_quartic_point": source["source_quartic_point"],
                "exact_signed_pair_relation_up_to_overall_sign": pair_relations.get(
                    canonical_short_point(point)
                ),
            }
        )
    return {
        "input": {
            "path": str(path),
            "sha256": sha256(raw).hexdigest(),
        },
        "parameter": artifact["parameter"],
        "certificate_short_model": [fraction_text(value) for value in model],
        "baseline_section_count": len(baseline),
        "charts": chart_records,
        "distinct_exact_points_beyond_the_listed_17": candidate_records,
        "exact_signed_pair_relation_count": sum(
            record["exact_signed_pair_relation_up_to_overall_sign"] is not None
            for record in candidate_records
        ),
        "all_reported_points_exactly_verified": True,
        "finite_quotient_escape": quotient,
        "status": (
            "EXACT_QUOTIENT_ESCAPE_FOUND"
            if quotient["maximum_marginal_dimension"] > 0
            else "COMPLETE_BOUNDED_SEARCH_NO_DETECTED_QUOTIENT_GAIN"
        ),
    }


def default_inputs() -> tuple[Path, ...]:
    root = REPOSITORY / "artifacts/local/elliptic-curves/q12o5867-specializations"
    preferred = (
        "q12o5867-specialization-m7801_1463.json",
        "q12o5867-specialization-m5954_7203.json",
        "q12o5867-specialization-601_418.json",
        "q12o5867-specialization-677_3402.json",
    )
    available = tuple(root / name for name in preferred if (root / name).is_file())
    if available:
        return available[:3]
    return tuple(sorted(root.glob("*.json")))[:3]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="*", type=Path)
    parser.add_argument(
        "--residual-selmer-gate",
        action="append",
        type=Path,
        required=True,
        help="repeat once per input, in input order",
    )
    parser.add_argument("--charts-per-candidate", type=int, default=2)
    parser.add_argument(
        "--charts-per-base",
        type=int,
        default=1,
        help="maximum diverse normalized slope charts retained for each base section",
    )
    parser.add_argument("--anchor-pool-size", type=int, default=8)
    parser.add_argument("--numerator-bound", type=int, default=200_000)
    parser.add_argument("--denominator-bound", type=int, default=2_000)
    parser.add_argument("--chart-timeout", type=float, default=120.0)
    parser.add_argument("--stack-bytes", type=int, default=256_000_000)
    parser.add_argument(
        "--ratpoints",
        type=Path,
        help="optional ratpoints-compatible CPU/CUDA executable replacing PARI",
    )
    parser.add_argument("--ratpoints-library", type=Path)
    parser.add_argument("--relation-primes", type=parse_int_tuple, default=(2, 3, 5))
    parser.add_argument("--reduction-prime-bound", type=int, default=500)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            REPOSITORY
            / "artifacts/generated-results/elliptic-curves"
            / "q12o5867-section-normalized-slope-slices.json"
        ),
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main() -> None:
    sys.set_int_max_str_digits(0)
    args = build_parser().parse_args()
    if args.charts_per_candidate < 1:
        raise SystemExit("--charts-per-candidate must be positive")
    if args.charts_per_base < 1:
        raise SystemExit("--charts-per-base must be positive")
    if args.charts_per_candidate > 17 * args.charts_per_base:
        raise SystemExit(
            "--charts-per-candidate exceeds 17 times --charts-per-base"
        )
    if not 3 <= args.anchor_pool_size <= 16:
        raise SystemExit("--anchor-pool-size must be between 3 and 16")
    if args.numerator_bound < 1 or args.denominator_bound < 1:
        raise SystemExit("search bounds must be positive")
    if args.denominator_bound > args.numerator_bound:
        raise SystemExit("PARI requires denominator bound <= numerator bound")
    if args.chart_timeout <= 0 or args.stack_bytes < 8_000_000:
        raise SystemExit("timeout and stack size must be positive")
    inputs = tuple(path.resolve() for path in (args.inputs or default_inputs()))
    if not inputs:
        raise SystemExit("no q12o5867 specialization artifacts were found")
    if len(args.residual_selmer_gate) != len(inputs):
        raise SystemExit("repeat --residual-selmer-gate exactly once per input")
    for path, gate_path in zip(inputs, args.residual_selmer_gate, strict=True):
        require_gate_for_specialization(gate_path, json.loads(path.read_text()))
    started = time.monotonic()
    results = [search_specialization(path, args) for path in inputs]
    command = " ".join(
        (
            ".venv/bin/python",
            "elliptic-curves/scripts/search_q12o5867_section_slope_slices.py",
            *(shlex.quote(str(path.relative_to(REPOSITORY))) for path in inputs),
            *(
                item
                for gate_path in args.residual_selmer_gate
                for item in (
                    "--residual-selmer-gate",
                    shlex.quote(str(gate_path)),
                )
            ),
            "--charts-per-candidate",
            str(args.charts_per_candidate),
            "--charts-per-base",
            str(args.charts_per_base),
            "--anchor-pool-size",
            str(args.anchor_pool_size),
            "--numerator-bound",
            str(args.numerator_bound),
            "--denominator-bound",
            str(args.denominator_bound),
            "--chart-timeout",
            str(args.chart_timeout),
            "--stack-bytes",
            str(args.stack_bytes),
            *(
                (
                    "--ratpoints",
                    shlex.quote(str(args.ratpoints)),
                    *(
                        (
                            "--ratpoints-library",
                            shlex.quote(str(args.ratpoints_library)),
                        )
                        if args.ratpoints_library is not None
                        else ()
                    ),
                )
                if args.ratpoints is not None
                else ()
            ),
            "--relation-primes",
            ",".join(str(value) for value in args.relation_primes),
            "--reduction-prime-bound",
            str(args.reduction_prime_bound),
            "--output",
            shlex.quote(str(args.output)),
            "--overwrite",
        )
    )
    artifact = {
        "schema": "elliptic-curves.q12o5867-section-normalized-slope-slices.v1",
        "status": (
            "EXACT_QUOTIENT_ESCAPE_FOUND"
            if any(
                result["finite_quotient_escape"]["maximum_marginal_dimension"] > 0
                for result in results
            )
            else "COMPLETE_BOUNDED_SEARCH_NO_DETECTED_QUOTIENT_GAIN"
        ),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "method": {
            "description": (
                "three generic-section secant slopes normalize each exact slope "
                "quartic at z=0,1,infinity before the selected point-search engine"
            ),
            "engine": (
                "PARI/GP hyperellratpoints"
                if args.ratpoints is None
                else "ratpoints-compatible executable"
            ),
            "pari_version": pari_version() if args.ratpoints is None else None,
            "ratpoints_executable": (
                str(args.ratpoints.resolve()) if args.ratpoints is not None else None
            ),
            "eclib_used": False,
            "charts_per_candidate": args.charts_per_candidate,
            "charts_per_base": args.charts_per_base,
            "anchor_pool_size": args.anchor_pool_size,
            "numerator_bound": args.numerator_bound,
            "denominator_bound": args.denominator_bound,
            "chart_timeout_seconds": args.chart_timeout,
            "stack_bytes": args.stack_bytes,
            "relation_primes": list(args.relation_primes),
            "reduction_prime_bound": args.reduction_prime_bound,
            "residual_selmer_gates": [
                {
                    "path": str(path.resolve()),
                    "sha256": sha256_file(path),
                }
                for path in args.residual_selmer_gate
            ],
        },
        "results": results,
        "maximum_exact_quotient_gain": max(
            result["finite_quotient_escape"]["maximum_marginal_dimension"]
            for result in results
        ),
        "wall_seconds": time.monotonic() - started,
        "script": {
            "path": str(Path(__file__).resolve().relative_to(REPOSITORY)),
            "sha256": sha256_file(Path(__file__).resolve()),
        },
        "reproducing_command": command,
        "claim_boundary": [
            "This is a complete enumeration only inside the declared chart boxes.",
            "A positive finite-quotient escape is exact.",
            "No detected escape is not a dependence result or rank upper bound.",
        ],
    }
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    mode = "w" if args.overwrite else "x"
    with output.open(mode) as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"status={artifact['status']}")
    print(f"candidates={len(results)}")
    print(f"maximum_exact_quotient_gain={artifact['maximum_exact_quotient_gain']}")
    for result in results:
        parameter = result["parameter"]["normalized_projective"]
        print(
            f"parameter={parameter} status={result['status']} "
            f"new_points={len(result['distinct_exact_points_beyond_the_listed_17'])} "
            "quotient_gain="
            f"{result['finite_quotient_escape']['maximum_marginal_dimension']}"
        )
    print(f"output={output}")


if __name__ == "__main__":
    main()
