#!/usr/bin/env python3
"""Skew-height and rational-chart point search on Nagao's ``u=42`` quartic.

The earlier uniform search enumerated quartic abscissas ``n/d`` with
``max(abs(n), d) <= 10^6``.  This script searches a different, explicitly
bounded region: ten denominator slabs with successively larger numerator
bounds, followed by 76 unimodular Mobius charts centred at the 38 nonvisible
abscissas in that checkpoint.  Every returned point is mapped to the short
Jacobian and checked exactly.

New images are compared with the small-prime-saturated 17-point basis.  PARI
height pairings are used only to *discover* integral relation vectors; each
reported relation is then replayed with exact Fraction group arithmetic.  The
height-matrix rank is also replayed at 72 and 120 decimal digits, and remains
numerical search evidence rather than a proof of independence.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path
import platform
import re
import shlex
import sys
from typing import Any, Iterable, Sequence

from ek_k3 import rational_to_string
from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK13_CONSTRUCTION,
    primitive_quartic_coefficients,
    quartic_point_to_short_jacobian,
    quartic_value,
    rank13_base_changed_short_jacobian_coefficients,
    rank13_base_parameter,
    rank13_known_quartic_points,
)
from pari_bridge import pari_version
from search_extra_points import (
    gp_rational,
    gp_vector,
    parse_point_vector,
    run_gp,
    signless_quartic_points,
)
from triage_nagao_rank13_finalists import (
    height_matrix_replay,
    point_digest,
    point_on_short_curve,
    quartic_gp_polynomial,
    stable_height_rank,
)


Q = Fraction
PARAMETER_U = 42
UNIFORM_CHECKPOINT_HEIGHT = 1_000_000
MOBIUS_HEIGHT = 50_000
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_u42_skew_height.py"
)


@dataclass(frozen=True)
class SearchBox:
    identifier: str
    numerator_bound: int
    denominator_lower: int
    denominator_upper: int

    def __post_init__(self) -> None:
        if self.numerator_bound < self.denominator_upper:
            raise ValueError("PARI requires the denominator bound not to exceed N")
        if not 1 <= self.denominator_lower <= self.denominator_upper:
            raise ValueError("invalid denominator interval")

    @property
    def gp_height(self) -> str:
        if self.denominator_lower == 1:
            return f"[{self.numerator_bound},{self.denominator_upper}]"
        return (
            f"[{self.numerator_bound},"
            f"[{self.denominator_lower},{self.denominator_upper}]]"
        )


# A deterministic staircase with roughly 10^11 numerator/denominator trials
# per slab.  The intervals are disjoint and cover every reduced denominator
# from 1 through 128000, but the numerator bound is slab-dependent.
SEARCH_BOXES = (
    SearchBox("d000001_000010", 10_000_000_000, 1, 10),
    SearchBox("d000011_000100", 1_000_000_000, 11, 100),
    SearchBox("d000101_001000", 100_000_000, 101, 1_000),
    SearchBox("d001001_002000", 100_000_000, 1_001, 2_000),
    SearchBox("d002001_004000", 50_000_000, 2_001, 4_000),
    SearchBox("d004001_008000", 25_000_000, 4_001, 8_000),
    SearchBox("d008001_016000", 12_500_000, 8_001, 16_000),
    SearchBox("d016001_032000", 6_250_000, 16_001, 32_000),
    SearchBox("d032001_064000", 3_125_000, 32_001, 64_000),
    SearchBox("d064001_128000", 1_562_500, 64_001, 128_000),
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def companion_section_x_values(parameter_t: Fraction) -> tuple[Fraction, ...]:
    """Five generic companion abscissas absent from the original 13 list."""

    parameter_t = Q(parameter_t)
    return (
        -parameter_t / 15 + Q(703, 15),
        7 * parameter_t / 15 + Q(928, 15),
        -7 * parameter_t / 15 + Q(928, 15),
        5 * parameter_t / 3 + Q(3628, 15),
        -5 * parameter_t / 3 + Q(3628, 15),
    )


def centered_unimodular_matrix(
    center: Fraction, shift: int = 0
) -> tuple[int, int, int, int]:
    """Return ``(a,b,c,d)`` with det 1 and ``b/d == center``.

    The chart is ``x=(a*X+b)/(c*X+d)``.  Adding ``shift*(b,d)`` to
    the first column gives another determinant-one orientation with the same
    centre at ``X=0``.
    """

    center = Q(center)
    b_value, d_value = center.numerator, center.denominator
    if b_value == 0:
        return (1, 0, shift, 1)
    a_value = pow(d_value, -1, abs(b_value))
    c_value = (a_value * d_value - 1) // b_value
    matrix = (
        a_value + shift * b_value,
        b_value,
        c_value + shift * d_value,
        d_value,
    )
    a_value, b_value, c_value, d_value = matrix
    if a_value * d_value - b_value * c_value != 1:
        raise AssertionError("the centred chart is not unimodular")
    return matrix


def transform_binary_quartic(
    coefficients: Sequence[Fraction], matrix: Sequence[int]
) -> tuple[Fraction, ...]:
    """Apply ``x=(aX+b)/(cX+d)``, including the fourth-power factor."""

    if len(coefficients) != 5 or len(matrix) != 4:
        raise ValueError("a quartic and a 2-by-2 matrix are required")
    a_value, b_value, c_value, d_value = (int(value) for value in matrix)
    if a_value * d_value - b_value * c_value == 0:
        raise ValueError("the chart matrix must be invertible")
    answer = [Q(0) for _ in range(5)]
    for power, coefficient in enumerate(coefficients):
        coefficient = Q(coefficient)
        for left_power in range(power + 1):
            left = (
                comb(power, left_power)
                * a_value**left_power
                * b_value ** (power - left_power)
            )
            for right_power in range(5 - power):
                right = (
                    comb(4 - power, right_power)
                    * c_value**right_power
                    * d_value ** (4 - power - right_power)
                )
                answer[left_power + right_power] += coefficient * left * right
    return tuple(answer)


def map_chart_point(
    point: tuple[Fraction, Fraction], matrix: Sequence[int]
) -> tuple[Fraction, Fraction] | None:
    """Map ``(X,Z)`` on a transformed quartic back to ``(x,z)``."""

    x_value, z_value = point
    a_value, b_value, c_value, d_value = (int(value) for value in matrix)
    denominator = c_value * x_value + d_value
    if denominator == 0:
        return None
    return (
        (a_value * x_value + b_value) / denominator,
        z_value / denominator**2,
    )


def short_add(
    coefficient_a: Fraction,
    first: tuple[Fraction, Fraction] | None,
    second: tuple[Fraction, Fraction] | None,
) -> tuple[Fraction, Fraction] | None:
    """Exact group addition on ``y^2=x^3+A*x+B`` (``None`` is infinity)."""

    if first is None:
        return second
    if second is None:
        return first
    x_first, y_first = first
    x_second, y_second = second
    if x_first == x_second and y_first == -y_second:
        return None
    if first == second:
        if y_first == 0:
            return None
        slope = (3 * x_first**2 + coefficient_a) / (2 * y_first)
    else:
        slope = (y_second - y_first) / (x_second - x_first)
    x_third = slope**2 - x_first - x_second
    y_third = slope * (x_first - x_third) - y_first
    return x_third, y_third


def short_multiply(
    coefficient_a: Fraction,
    point: tuple[Fraction, Fraction],
    multiplier: int,
) -> tuple[Fraction, Fraction] | None:
    """Exact double-and-add scalar multiplication on a short model."""

    if multiplier < 0:
        point = point[0], -point[1]
        multiplier = -multiplier
    answer = None
    addend: tuple[Fraction, Fraction] | None = point
    while multiplier:
        if multiplier & 1:
            answer = short_add(coefficient_a, answer, addend)
        addend = short_add(coefficient_a, addend, addend)
        multiplier >>= 1
    return answer


def exact_linear_combination(
    coefficient_a: Fraction,
    basis: Sequence[tuple[Fraction, Fraction]],
    relation: Sequence[int],
) -> tuple[Fraction, Fraction] | None:
    if len(basis) != len(relation):
        raise ValueError("the relation length does not match the basis")
    answer = None
    for point, multiplier in zip(basis, relation):
        answer = short_add(
            coefficient_a,
            answer,
            short_multiply(coefficient_a, point, int(multiplier)),
        )
    return answer


def run_skew_box(
    coefficients: Sequence[Fraction],
    search_box: SearchBox,
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], int, float]:
    program = "\n".join(
        (
            f"Q={quartic_gp_polynomial(coefficients)};",
            "gettime();",
            f"R=hyperellratpoints(Q,{search_box.gp_height});",
            'print("PARI_MILLISECONDS ",gettime());',
            'print("POINTS ",R);',
            "quit",
        )
    ) + "\n"
    output, wall_seconds = run_gp(
        program, timeout=timeout, stack_bytes=stack_bytes
    )
    match = re.search(r"PARI_MILLISECONDS (\d+)", output)
    if match is None or "POINTS " not in output:
        raise AssertionError("PARI omitted a skew-box result marker")
    points = parse_point_vector(output.split("POINTS ", 1)[1])
    return points, int(match.group(1)), wall_seconds


def run_mobius_charts(
    coefficients: Sequence[Fraction],
    charts: Sequence[tuple[str, tuple[int, int, int, int]]],
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[dict[str, tuple[tuple[Fraction, Fraction], ...]], dict[str, int], float]:
    commands: list[str] = []
    for index, (identifier, matrix) in enumerate(charts):
        transformed = transform_binary_quartic(coefficients, matrix)
        commands.extend(
            (
                f"Q={quartic_gp_polynomial(transformed)};",
                "gettime();",
                f"R=hyperellratpoints(Q,{height_bound});",
                f'print("CHART_{index}_BEGIN");',
                'print("PARI_MILLISECONDS ",gettime());',
                'print("POINTS ",R);',
                f'print("CHART_{index}_END");',
            )
        )
    commands.append("quit")
    output, wall_seconds = run_gp(
        "\n".join(commands) + "\n", timeout=timeout, stack_bytes=stack_bytes
    )
    point_results: dict[str, tuple[tuple[Fraction, Fraction], ...]] = {}
    millisecond_results: dict[str, int] = {}
    for index, (identifier, _) in enumerate(charts):
        begin = f"CHART_{index}_BEGIN\n"
        end = f"\nCHART_{index}_END"
        if begin not in output:
            raise AssertionError(f"PARI omitted {identifier}")
        block = output.split(begin, 1)[1].split(end, 1)[0]
        match = re.search(r"PARI_MILLISECONDS (\d+)", block)
        if match is None or "POINTS " not in block:
            raise AssertionError(f"PARI omitted output within {identifier}")
        point_results[identifier] = parse_point_vector(block.split("POINTS ", 1)[1])
        millisecond_results[identifier] = int(match.group(1))
    return point_results, millisecond_results, wall_seconds


def discover_relations(
    coefficients: Sequence[Fraction],
    basis: Sequence[tuple[Fraction, Fraction]],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[int, ...], ...]:
    curve = ",".join(gp_rational(Q(value)) for value in coefficients)
    commands = [
        "default(realprecision,120);",
        f"E=ellinit([{curve}]);",
        f"B=[{','.join(gp_vector(point) for point in basis)}];",
        "H=ellheightmatrix(E,B);",
    ]
    for index, point in enumerate(points):
        commands.extend(
            (
                f"Q={gp_vector(point)};",
                "V=vector(#B,j,ellheight(E,B[j],Q))~;",
                "C=round(matsolve(H,V));",
                "S=[0];for(j=1,#B,S=elladd(E,S,ellmul(E,B[j],C[j])));",
                f'print("RELATION_{index} ",Vec(C)," EXACT ",S==Q);',
            )
        )
    commands.append("quit")
    output, _ = run_gp(
        "\n".join(commands) + "\n", timeout=timeout, stack_bytes=stack_bytes
    )
    relations = []
    for index in range(len(points)):
        match = re.search(
            rf"^RELATION_{index} \[(.*?)\] EXACT ([01])$", output, re.MULTILINE
        )
        if match is None or match.group(2) != "1":
            raise AssertionError(f"PARI did not verify relation {index}")
        relation = tuple(int(value.strip()) for value in match.group(1).split(","))
        if len(relation) != len(basis):
            raise AssertionError("a relation has the wrong length")
        relations.append(relation)
    return tuple(relations)


def _record_point(
    quartic_point: tuple[Fraction, Fraction],
    jacobian_point: tuple[Fraction, Fraction],
    relation: Sequence[int],
    sources: Sequence[str],
) -> dict[str, Any]:
    return {
        "quartic_x": rational_to_string(quartic_point[0]),
        "quartic_z": rational_to_string(quartic_point[1]),
        "jacobian_x": rational_to_string(jacobian_point[0]),
        "jacobian_y": rational_to_string(jacobian_point[1]),
        "discovered_in": list(sources),
        "saturated_basis_relation": list(relation),
        "exact_quartic_membership_checked": True,
        "exact_jacobian_membership_checked": True,
        "exact_relation_replayed_with_fraction_group_law": True,
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts" / "generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--checkpoint-input",
        type=Path,
        default=generated / "elliptic_nagao_rank13_finalist_triage.json",
    )
    parser.add_argument(
        "--saturation-input",
        type=Path,
        default=generated / "elliptic_nagao_u42_height_10000000.json",
    )
    parser.add_argument("--box-timeout", type=float, default=20.0)
    parser.add_argument("--chart-timeout", type=float, default=25.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--stack-bytes", type=int, default=256_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_nagao_u42_skew_height.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 0 < args.box_timeout <= 30:
        raise SystemExit("--box-timeout must be in (0,30]")
    if not 0 < args.chart_timeout <= 30:
        raise SystemExit("--chart-timeout must be in (0,30]")
    if not 0 < args.height_timeout <= 60:
        raise SystemExit("--height-timeout must be in (0,60]")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes is too small")

    parameter_t = rank13_base_parameter(Q(PARAMETER_U))
    quartic_coefficients = primitive_quartic_coefficients(
        RANK13_CONSTRUCTION, parameter_t
    )
    jacobian_coefficients = rank13_base_changed_short_jacobian_coefficients(
        Q(PARAMETER_U)
    )
    if any(jacobian_coefficients[index] for index in range(3)):
        raise AssertionError("the exact relation verifier requires a short model")

    checkpoint_data = json.loads(args.checkpoint_input.read_text())
    checkpoint = next(
        record
        for record in checkpoint_data["candidates"]
        if int(record["parameter_u"]) == PARAMETER_U
    )
    escalated = checkpoint["escalated_bounded_search"]
    if int(escalated["quartic_naive_height_bound"]) != UNIFORM_CHECKPOINT_HEIGHT:
        raise AssertionError("the input is not the pinned height-10^6 checkpoint")
    checkpoint_records = escalated["new_point_records"]
    checkpoint_x = {point[0] for point in rank13_known_quartic_points(parameter_t)}
    checkpoint_x.update(Q(record["quartic_x"]) for record in checkpoint_records)
    if len(checkpoint_x) != int(escalated["distinct_quartic_x_values"]):
        raise AssertionError("the checkpoint abscissa count changed")

    companions = companion_section_x_values(parameter_t)
    if not set(companions).issubset(checkpoint_x):
        raise AssertionError("a generic companion section is absent from checkpoint")

    saturation_data = json.loads(args.saturation_input.read_text())
    saturated_basis = tuple(
        (Q(record["jacobian_x"]), Q(record["jacobian_y"]))
        for record in saturation_data["small_prime_saturation"]["saturated_basis"]
    )
    if len(saturated_basis) != 17:
        raise AssertionError("the pinned saturated basis no longer has 17 points")
    if any(
        not point_on_short_curve(jacobian_coefficients, point)
        for point in saturated_basis
    ):
        raise AssertionError("a saturated-basis point missed the exact curve")

    discovered: dict[Fraction, tuple[tuple[Fraction, Fraction], list[str]]] = {}
    box_records = []
    for search_box in SEARCH_BOXES:
        raw_points, pari_milliseconds, wall_seconds = run_skew_box(
            quartic_coefficients,
            search_box,
            timeout=args.box_timeout,
            stack_bytes=args.stack_bytes,
        )
        signless = signless_quartic_points(raw_points)
        new_x_values = []
        for point in signless:
            if point[1] ** 2 != quartic_value(quartic_coefficients, point[0]):
                raise AssertionError("PARI returned a point off the exact quartic")
            if point[0] in checkpoint_x:
                continue
            new_x_values.append(point[0])
            if point[0] not in discovered:
                discovered[point[0]] = (point, [f"skew:{search_box.identifier}"])
            else:
                discovered[point[0]][1].append(f"skew:{search_box.identifier}")
        box_records.append(
            {
                "id": search_box.identifier,
                "numerator_absolute_bound": search_box.numerator_bound,
                "denominator_lower_bound": search_box.denominator_lower,
                "denominator_upper_bound": search_box.denominator_upper,
                "signed_points_found": len(raw_points),
                "distinct_quartic_x_values": len(signless),
                "outside_uniform_checkpoint_x_values": [
                    rational_to_string(value) for value in new_x_values
                ],
                "pari_reported_milliseconds": pari_milliseconds,
                "wall_seconds": wall_seconds,
            }
        )
        print(
            f"{search_box.identifier}: x={len(signless)} "
            f"outside={len(new_x_values)}",
            flush=True,
        )

    charts = []
    for index, record in enumerate(checkpoint_records):
        center = Q(record["quartic_x"])
        for shift in (0, -1):
            identifier = f"checkpoint_{index:02d}_shift_{shift}"
            charts.append((identifier, centered_unimodular_matrix(center, shift)))
    chart_points, chart_milliseconds, chart_wall_seconds = run_mobius_charts(
        quartic_coefficients,
        charts,
        height_bound=MOBIUS_HEIGHT,
        timeout=args.chart_timeout,
        stack_bytes=args.stack_bytes,
    )
    chart_records = []
    chart_new_x: set[Fraction] = set()
    for identifier, matrix in charts:
        raw_points = chart_points[identifier]
        mapped_by_x: dict[Fraction, tuple[Fraction, Fraction]] = {}
        for transformed_point in signless_quartic_points(raw_points):
            mapped = map_chart_point(transformed_point, matrix)
            if mapped is None:
                continue
            if mapped[1] ** 2 != quartic_value(quartic_coefficients, mapped[0]):
                raise AssertionError("a chart point missed the original quartic")
            mapped_by_x.setdefault(mapped[0], mapped)
        outside = []
        for x_value, point in mapped_by_x.items():
            if x_value in checkpoint_x:
                continue
            outside.append(x_value)
            chart_new_x.add(x_value)
            if x_value not in discovered:
                discovered[x_value] = (point, [f"chart:{identifier}"])
            else:
                discovered[x_value][1].append(f"chart:{identifier}")
        chart_records.append(
            {
                "id": identifier,
                "matrix_a_b_c_d": list(matrix),
                "determinant": matrix[0] * matrix[3] - matrix[1] * matrix[2],
                "transformed_naive_height_bound": MOBIUS_HEIGHT,
                "signed_points_found": len(raw_points),
                "distinct_finite_original_x_values": len(mapped_by_x),
                "outside_uniform_checkpoint_x_values": [
                    rational_to_string(value) for value in outside
                ],
                "pari_reported_milliseconds": chart_milliseconds[identifier],
            }
        )
    print(
        f"mobius charts: charts={len(charts)} outside={len(chart_new_x)}",
        flush=True,
    )

    ordered_discoveries = tuple(discovered.values())
    quartic_points = tuple(item[0] for item in ordered_discoveries)
    jacobian_points = tuple(
        quartic_point_to_short_jacobian(
            RANK13_CONSTRUCTION, parameter_t, quartic_point
        )
        for quartic_point in quartic_points
    )
    if any(
        not point_on_short_curve(jacobian_coefficients, point)
        for point in jacobian_points
    ):
        raise AssertionError("a mapped point missed the exact Jacobian")

    relations = discover_relations(
        jacobian_coefficients,
        saturated_basis,
        jacobian_points,
        timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    coefficient_a = jacobian_coefficients[3]
    for point, relation in zip(jacobian_points, relations):
        if exact_linear_combination(coefficient_a, saturated_basis, relation) != point:
            raise AssertionError("an exact saturated-basis relation failed in Python")

    height_runs = height_matrix_replay(
        jacobian_coefficients,
        saturated_basis + jacobian_points,
        precisions=(72, 120),
        timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    stable_rank = stable_height_rank(height_runs)
    if stable_rank != 17:
        raise AssertionError("the pinned run unexpectedly changed numerical rank")

    point_records = tuple(
        _record_point(quartic_point, jacobian_point, relation, sources)
        for (quartic_point, sources), jacobian_point, relation in zip(
            ordered_discoveries, jacobian_points, relations
        )
    )
    script_path = Path(__file__).resolve()
    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    artifact = {
        "schema_version": 1,
        "status": "bounded_search_complete",
        "candidate": {
            "id": "nagao-u-42",
            "parameter_u": PARAMETER_U,
            "parameter_t": rational_to_string(parameter_t),
            "log_conductor": checkpoint["log_conductor"],
            "root_number": checkpoint["root_number"],
        },
        "primary_source": PRIMARY_SOURCE,
        "inputs": {
            str(args.checkpoint_input): sha256_file(args.checkpoint_input),
            str(args.saturation_input): sha256_file(args.saturation_input),
        },
        "parameters": {
            "uniform_checkpoint_height": UNIFORM_CHECKPOINT_HEIGHT,
            "skew_boxes": [
                {
                    "id": item.identifier,
                    "numerator_absolute_bound": item.numerator_bound,
                    "denominator_lower_bound": item.denominator_lower,
                    "denominator_upper_bound": item.denominator_upper,
                }
                for item in SEARCH_BOXES
            ],
            "mobius_chart_height": MOBIUS_HEIGHT,
            "mobius_center_count": len(checkpoint_records),
            "mobius_shifts": [0, -1],
            "box_timeout_seconds_each": args.box_timeout,
            "chart_batch_timeout_seconds": args.chart_timeout,
            "height_and_relation_timeout_seconds_each": args.height_timeout,
            "stack_bytes": args.stack_bytes,
        },
        "generic_companion_sections": {
            "count": len(companions),
            "quartic_x_values": [rational_to_string(value) for value in companions],
            "classification": (
                "generic companion sections, already present in the height-10^6 "
                "checkpoint; excluded from specialization-only discoveries"
            ),
        },
        "skew_search": {
            "scope": (
                "union of the ten explicit reduced n/d boxes; each record gives "
                "|n| and denominator bounds"
            ),
            "boxes": box_records,
        },
        "mobius_search": {
            "scope": (
                "76 explicit determinant-one charts x=(aX+b)/(cX+d), each "
                f"searched to transformed naive height {MOBIUS_HEIGHT}"
            ),
            "batch_wall_seconds": chart_wall_seconds,
            "charts": chart_records,
        },
        "outside_uniform_checkpoint": {
            "distinct_quartic_x_count": len(point_records),
            "distinct_jacobian_image_count": len(point_records),
            "all_mapped_and_checked_exactly": True,
            "all_exactly_in_saturated_17_point_span": True,
            "points": point_records,
        },
        "height_matrix_runs": height_runs,
        "stable_pool_numerical_rank": stable_rank,
        "saturated_basis": {
            "point_count": len(saturated_basis),
            "sha256": point_digest(saturated_basis),
            "provenance": args.saturation_input.name,
        },
        "interpretation": {
            "exact": (
                "all points returned outside the prior uniform box lie on the "
                "quartic and Jacobian and have displayed, exactly replayed integral "
                "relations in the pinned 17-point subgroup"
            ),
            "numerical": (
                "the augmented pool has stable height-matrix rank 17 at 72 and "
                "120 decimal digits"
            ),
            "not_claimed": (
                "this bounded result neither proves the full Mordell-Weil rank is "
                "17 nor excludes independent points outside the declared boxes/charts"
            ),
        },
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "platform": platform.platform(),
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "actual_command": command,
        "script_sha256": sha256_file(script_path),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        f"wrote {args.output}: outside={len(point_records)} stable_rank={stable_rank}",
        flush=True,
    )


if __name__ == "__main__":
    main()
