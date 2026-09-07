#!/usr/bin/env python3
"""Exact-point and bounded-rank triage for Nagao rank-13 finalists.

This script consumes the leakage-free integer-``u`` score/conductor scans and
replays every retained specialization on Nagao's primitive quartic.  For each
candidate it

* checks all thirteen displayed affine sections exactly;
* adds the exact Jacobian image obtained by taking the covariant-map limit at
  a split point at infinity;
* enumerates affine quartic points with PARI ``hyperellratpoints`` in a fixed
  naive-height box and maps every nonvisible point exactly;
* replays the baseline and augmented Neron--Tate height matrices at two
  precisions; and
* optionally runs strictly time-bounded PARI ``ellrank`` effort-zero,
  ``ellsaturation``, and analytic-rank probes.

The exact statements here are point membership and the bounded enumeration.
Height-matrix rank is numerical evidence, while equal ``ellrank`` bounds are
reported as a PARI software computation rather than a portable proof
certificate.  Root numbers and failed analytic probes are only prioritization
data.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

from ek_k3 import rational_to_string
from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK13_CONSTRUCTION,
    primitive_quartic_coefficients,
    quartic_point_to_short_jacobian,
    rank13_base_changed_short_jacobian_coefficients,
    rank13_base_parameter,
    rank13_known_quartic_points,
    rank13_leading_square,
)
from pari_bridge import pari_version
from search_extra_points import (
    gp_rational,
    gp_vector,
    parse_point_vector,
    parse_precisions,
    parse_vecsmall,
    run_gp,
    signless_quartic_points,
)


Q = Fraction
TARGET_LOG_CONDUCTOR = Decimal("182.72")
DEFAULT_RANK_U = (84, 42, 2, 11, 189, 171, 62, 50)
DEFAULT_SATURATION_U = (84, 42)
DEFAULT_ANALYTIC_U = (84, 2, 189)
DEFAULT_ESCALATE_U = (84, 42)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/triage_nagao_rank13_finalists.py"
)


@dataclass(frozen=True)
class Finalist:
    parameter_u: int
    parameter_t: Fraction
    score: str
    last_numerical_prime: int
    log_conductor: str
    root_number: int
    source_artifact: str

    @property
    def identifier(self) -> str:
        return f"nagao-u-{self.parameter_u}"

    @property
    def below_target(self) -> bool:
        return Decimal(self.log_conductor) < TARGET_LOG_CONDUCTOR


def parse_u_values(value: str) -> tuple[int, ...]:
    if value == "":
        return ()
    try:
        values = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("u values must be integers") from error
    if not values or any(parameter_u <= 0 for parameter_u in values):
        raise argparse.ArgumentTypeError("u values must be distinct positive integers")
    if len(set(values)) != len(values):
        raise argparse.ArgumentTypeError("u values must be distinct positive integers")
    return values


def _artifact_finalists(path: Path) -> dict[int, Finalist]:
    data = json.loads(path.read_text())
    records = data.get("final_conductor_candidates")
    if not isinstance(records, list):
        raise ValueError(f"{path} has no final_conductor_candidates list")
    answer: dict[int, Finalist] = {}
    for record in records:
        parameter_u = int(record["parameter_u"])
        parameter_t = Q(record["parameter_t"])
        if parameter_t != rank13_base_parameter(Q(parameter_u)):
            raise AssertionError(f"u={parameter_u} has an inconsistent T value")
        pari = record["pari"]
        candidate = Finalist(
            parameter_u=parameter_u,
            parameter_t=parameter_t,
            score=str(record["score"]),
            last_numerical_prime=int(record["last_numerical_prime"]),
            log_conductor=str(pari["log_conductor"]),
            root_number=int(pari["root_number"]),
            source_artifact=path.name,
        )
        if parameter_u in answer:
            raise ValueError(f"duplicate u={parameter_u} in {path}")
        answer[parameter_u] = candidate
    return answer


def load_finalists(
    primary_input: Path,
    additional_input: Path | None,
    additional_u: Sequence[int],
) -> tuple[Finalist, ...]:
    """Load every primary finalist plus named candidates from a larger scan."""

    primary = _artifact_finalists(primary_input)
    ordered = list(primary.values())
    seen = set(primary)
    if additional_u:
        if additional_input is None:
            raise ValueError("additional u values require an additional input")
        additional = _artifact_finalists(additional_input)
        for parameter_u in additional_u:
            if parameter_u not in additional:
                raise ValueError(
                    f"u={parameter_u} is absent from {additional_input.name}"
                )
            if parameter_u not in seen:
                ordered.append(additional[parameter_u])
                seen.add(parameter_u)
    return tuple(ordered)


def point_on_short_curve(
    coefficients: Sequence[Fraction], point: tuple[Fraction, Fraction]
) -> bool:
    if len(coefficients) != 5:
        raise ValueError("a Weierstrass coefficient vector has five entries")
    a1, a2, a3, a4, a6 = (Q(value) for value in coefficients)
    x_value, y_value = (Q(value) for value in point)
    return (
        y_value**2 + a1 * x_value * y_value + a3 * y_value
        == x_value**3 + a2 * x_value**2 + a4 * x_value + a6
    )


def split_infinity_jacobian_point(
    parameter_u: Fraction, *, sign: int = 1
) -> tuple[Fraction, Fraction]:
    """Return the exact covariant-map limit at a split quartic infinity.

    For ascending quartic coefficients ``(e,d,c,b,a)`` and square root ``s``
    of ``a``, the leading terms of the binary-quartic covariants give

    ``X=36*g0/a`` and ``Y=54*(a*g1-b*g0)/s^3``.

    Changing the chosen point at infinity negates ``Y``.
    """

    if sign not in (-1, 1):
        raise ValueError("the infinity sign must be -1 or 1")
    parameter_u = Q(parameter_u)
    parameter_t = rank13_base_parameter(parameter_u)
    _, d, c, b, a = primitive_quartic_coefficients(
        RANK13_CONSTRUCTION, parameter_t
    )
    if a == 0:
        raise ValueError("the quartic leading coefficient vanished")
    square_root = rank13_leading_square(parameter_u)
    if square_root**2 != a:
        raise AssertionError("the base-changed leading coefficient is not a square")
    g0 = b**2 / 16 - a * c / 6
    g1 = b * c / 12 - a * d / 2
    point = (
        36 * g0 / a,
        sign * 54 * (a * g1 - b * g0) / square_root**3,
    )
    coefficients = rank13_base_changed_short_jacobian_coefficients(parameter_u)
    if not point_on_short_curve(coefficients, point):
        raise AssertionError("the split-infinity covariant limit missed the Jacobian")
    return point


def point_digest(points: Iterable[tuple[Fraction, Fraction]]) -> str:
    text = "\n".join(
        f"{rational_to_string(x_value)},{rational_to_string(y_value)}"
        for x_value, y_value in points
    )
    return hashlib.sha256(text.encode()).hexdigest()


def quartic_gp_polynomial(coefficients: Sequence[Fraction]) -> str:
    return "+".join(
        f"{gp_rational(Q(coefficient))}*x^{power}"
        for power, coefficient in enumerate(coefficients)
    )


def bounded_quartic_points(
    parameter_t: Fraction,
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], float, int]:
    coefficients = primitive_quartic_coefficients(RANK13_CONSTRUCTION, parameter_t)
    program = "\n".join(
        (
            f"Q={quartic_gp_polynomial(coefficients)};",
            "gettime();",
            f"R=hyperellratpoints(Q,{height_bound});",
            'print("PARI_MILLISECONDS ",gettime());',
            'print("POINTS ",R);',
            "quit",
        )
    ) + "\n"
    output, wall_seconds = run_gp(program, timeout=timeout, stack_bytes=stack_bytes)
    match = re.search(r"PARI_MILLISECONDS (\d+)", output)
    if match is None or "POINTS " not in output:
        raise AssertionError("PARI omitted bounded-search output")
    points = parse_point_vector(output.split("POINTS ", 1)[1])
    return points, wall_seconds, int(match.group(1))


def exact_jacobian_record(
    quartic_point: tuple[Fraction, Fraction],
    jacobian_point: tuple[Fraction, Fraction],
) -> dict[str, Any]:
    return {
        "quartic_x": rational_to_string(quartic_point[0]),
        "quartic_z": rational_to_string(quartic_point[1]),
        "jacobian_x": rational_to_string(jacobian_point[0]),
        "jacobian_y": rational_to_string(jacobian_point[1]),
        "exact_quartic_membership_checked": True,
        "exact_jacobian_membership_checked": True,
    }


def height_matrix_replay(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    precisions: tuple[int, ...],
    timeout: float,
    stack_bytes: int,
) -> tuple[dict[str, Any], ...]:
    if not points:
        raise ValueError("height replay needs at least one point")
    if any(not point_on_short_curve(coefficients, point) for point in points):
        raise AssertionError("a height-replay point is not on the exact curve")
    curve = ",".join(gp_rational(Q(value)) for value in coefficients)
    point_vector = ",".join(gp_vector(point) for point in points)
    commands = [f"E=ellinit([{curve}]);", f"P=[{point_vector}];"]
    for precision in precisions:
        commands.extend(
            (
                f"default(realprecision,{precision});",
                "H=ellheightmatrix(E,P);",
                "IX=matindexrank(H);",
                "K=vecextract(P,IX[2]);",
                "HK=ellheightmatrix(E,K);",
                f'print("HEIGHT_{precision}_BEGIN");',
                "print(matrank(H));",
                "print(IX[2]);",
                "print(matdet(HK));",
                "EV=mateigen(HK,1)[1];print(vecmin(EV));print(vecmax(EV));",
                f'print("HEIGHT_{precision}_END");',
            )
        )
    commands.append("quit")
    output, wall_seconds = run_gp(
        "\n".join(commands) + "\n", timeout=timeout, stack_bytes=stack_bytes
    )
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    records = []
    for precision in precisions:
        start = lines.index(f"HEIGHT_{precision}_BEGIN") + 1
        end = lines.index(f"HEIGHT_{precision}_END")
        values = lines[start:end]
        records.append(
            {
                "decimal_precision": precision,
                "numerical_rank": int(values[0]),
                "subset_indices_one_based": parse_vecsmall(values[1]),
                "subset_height_determinant": values[2],
                "subset_smallest_eigenvalue": values[3],
                "subset_largest_eigenvalue": values[4],
                "gp_process_wall_seconds": wall_seconds,
            }
        )
    return tuple(records)


def stable_height_rank(runs: Sequence[dict[str, Any]]) -> int:
    ranks = {int(run["numerical_rank"]) for run in runs}
    subsets = {tuple(run["subset_indices_one_based"]) for run in runs}
    if len(ranks) != 1 or len(subsets) != 1:
        raise AssertionError("height rank or selected subset changed with precision")
    return next(iter(ranks))


def _bounded_gp_probe(
    program: str, *, timeout: float, stack_bytes: int
) -> tuple[str | None, float, str, str | None]:
    """Run GP and return output, elapsed, status, and a bounded error message."""

    executable = shutil.which("gp")
    if executable is None:
        return None, 0.0, "unavailable", "PARI/GP executable 'gp' was not found"
    started = time.monotonic()
    try:
        result = subprocess.run(
            [executable, "-q", "-s", str(stack_bytes)],
            input=program,
            text=True,
            capture_output=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return None, time.monotonic() - started, "timeout", None
    elapsed = time.monotonic() - started
    if result.returncode != 0 or "***" in result.stderr:
        message = " ".join(result.stderr.split())
        return None, elapsed, "pari_error", message[:1000]
    return result.stdout, elapsed, "completed", None


def ellrank_probe(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    curve = ",".join(gp_rational(Q(value)) for value in coefficients)
    point_vector = ",".join(gp_vector(point) for point in points)
    program = "\n".join(
        (
            "default(realprecision,80);",
            f"E=ellinit([{curve}]);",
            f"P=[{point_vector}];",
            'print("ON_CURVE ",vecsum(vector(#P,i,ellisoncurve(E,P[i]))));',
            "gettime();",
            "R=ellrank(E,0,P);",
            'print("PARI_MILLISECONDS ",gettime());',
            'print("RANK_LOWER ",R[1]);',
            'print("RANK_UPPER ",R[2]);',
            'print("RETURNED_POINTS ",#R[4]);',
            "quit",
        )
    ) + "\n"
    output, elapsed, status, error = _bounded_gp_probe(
        program, timeout=timeout, stack_bytes=stack_bytes
    )
    record: dict[str, Any] = {
        "status": status,
        "effort": 0,
        "supplied_numerically_independent_points": len(points),
        "timeout_seconds": timeout,
        "wall_seconds": elapsed,
        "pari_stack_bytes": stack_bytes,
    }
    if error is not None:
        record["error"] = error
    if output is not None:
        def integer(label: str) -> int:
            match = re.search(rf"^{label} (\d+)$", output, re.MULTILINE)
            if match is None:
                raise AssertionError(f"PARI omitted {label}")
            return int(match.group(1))

        record.update(
            {
                "exact_supplied_points_on_curve": integer("ON_CURVE"),
                "pari_milliseconds": integer("PARI_MILLISECONDS"),
                "lower_bound": integer("RANK_LOWER"),
                "upper_bound": integer("RANK_UPPER"),
                "returned_independent_points": integer("RETURNED_POINTS"),
            }
        )
        record["equal_returned_bounds"] = (
            record["lower_bound"] == record["upper_bound"]
        )
        record["interpretation"] = (
            "PARI effort-zero returned equal computational bounds; no portable "
            "exact independence or descent certificate is stored"
            if record["equal_returned_bounds"]
            else "PARI effort-zero returned a nontrivial rank interval"
        )
    elif status == "timeout":
        record["interpretation"] = "strict bounded probe timed out without rank bounds"
    else:
        record["interpretation"] = "probe failed without rank bounds"
    return record


def saturation_probe(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    prime_bound: int,
    timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    curve = ",".join(gp_rational(Q(value)) for value in coefficients)
    point_vector = ",".join(gp_vector(point) for point in points)
    program = "\n".join(
        (
            "default(realprecision,80);",
            f"E=ellinit([{curve}]);",
            f"P=[{point_vector}];",
            "H0=ellheightmatrix(E,P);",
            "gettime();",
            f"S=ellsaturation(E,P,{prime_bound});",
            'print("PARI_MILLISECONDS ",gettime());',
            'print("RETURNED_COUNT ",#S);',
            'print("ON_CURVE ",vecsum(vector(#S,i,ellisoncurve(E,S[i]))));',
            "H1=ellheightmatrix(E,S);",
            'print("ORIGINAL_DET ",matdet(H0));',
            'print("SATURATED_DET ",matdet(H1));',
            'print("DET_RATIO ",matdet(H0)/matdet(H1));',
            "quit",
        )
    ) + "\n"
    output, elapsed, status, error = _bounded_gp_probe(
        program, timeout=timeout, stack_bytes=stack_bytes
    )
    record: dict[str, Any] = {
        "status": status,
        "prime_bound_strict_upper_limit": prime_bound,
        "input_point_count": len(points),
        "timeout_seconds": timeout,
        "wall_seconds": elapsed,
        "pari_stack_bytes": stack_bytes,
        "scope_warning": (
            "PARI documents ellsaturation under a finite-index hypothesis; the "
            "full Mordell-Weil rank is unknown here, so this is only a subgroup "
            "basis/small-prime saturation computation"
        ),
    }
    if error is not None:
        record["error"] = error
    if output is not None:
        def value(label: str) -> str:
            match = re.search(rf"^{label} (.+)$", output, re.MULTILINE)
            if match is None:
                raise AssertionError(f"PARI omitted {label}")
            return match.group(1)

        record.update(
            {
                "pari_milliseconds": int(value("PARI_MILLISECONDS")),
                "returned_point_count": int(value("RETURNED_COUNT")),
                "exact_returned_points_on_curve": int(value("ON_CURVE")),
                "original_height_determinant": value("ORIGINAL_DET"),
                "returned_height_determinant": value("SATURATED_DET"),
                "height_determinant_ratio": value("DET_RATIO"),
            }
        )
    return record


def analytic_rank_probe(
    coefficients: Sequence[Fraction], *, timeout: float, stack_bytes: int
) -> dict[str, Any]:
    curve = ",".join(gp_rational(Q(value)) for value in coefficients)
    program = "\n".join(
        (
            "default(realprecision,80);",
            f"E=ellminimalmodel(ellinit([{curve}]));",
            "gettime();",
            "R=ellanalyticrank(E);",
            'print("PARI_MILLISECONDS ",gettime());',
            'print("ANALYTIC_RANK ",R[1]);',
            'print("LEADING_VALUE ",R[2]);',
            "quit",
        )
    ) + "\n"
    output, elapsed, status, error = _bounded_gp_probe(
        program, timeout=timeout, stack_bytes=stack_bytes
    )
    record: dict[str, Any] = {
        "status": status,
        "timeout_seconds": timeout,
        "wall_seconds": elapsed,
        "pari_stack_bytes": stack_bytes,
        "interpretation": (
            "strictly numerical L-series probe; never an algebraic-rank certificate"
        ),
    }
    if error is not None:
        record["error"] = error
    if output is not None:
        rank_match = re.search(r"^ANALYTIC_RANK (\d+)$", output, re.MULTILINE)
        value_match = re.search(r"^LEADING_VALUE (.+)$", output, re.MULTILINE)
        milliseconds = re.search(
            r"^PARI_MILLISECONDS (\d+)$", output, re.MULTILINE
        )
        if rank_match is None or value_match is None or milliseconds is None:
            raise AssertionError("PARI omitted analytic-rank output")
        record.update(
            {
                "pari_milliseconds": int(milliseconds.group(1)),
                "numerical_analytic_rank": int(rank_match.group(1)),
                "leading_value": value_match.group(1),
            }
        )
    return record


def exact_candidate_triage(
    candidate: Finalist,
    *,
    height_bound: int,
    precisions: tuple[int, ...],
    search_timeout: float,
    height_timeout: float,
    stack_bytes: int,
    escalated_height_bound: int | None = None,
    escalation_timeout: float | None = None,
) -> tuple[dict[str, Any], tuple[tuple[Fraction, Fraction], ...]]:
    parameter_u = Q(candidate.parameter_u)
    parameter_t = candidate.parameter_t
    coefficients = rank13_base_changed_short_jacobian_coefficients(parameter_u)

    affine_sections = rank13_known_quartic_points(parameter_t)
    if len(affine_sections) != 13 or len(set(affine_sections)) != 13:
        raise AssertionError("the thirteen specialized affine sections collided")
    affine_images = tuple(
        quartic_point_to_short_jacobian(RANK13_CONSTRUCTION, parameter_t, point)
        for point in affine_sections
    )
    if any(not point_on_short_curve(coefficients, point) for point in affine_images):
        raise AssertionError("an exact affine-section image missed the Jacobian")
    infinity = split_infinity_jacobian_point(parameter_u)
    baseline_points = affine_images + (infinity,)
    baseline_provenance = tuple(
        f"affine-section-{index}" for index in range(1, 14)
    ) + ("split-infinity",)

    baseline_runs = height_matrix_replay(
        coefficients,
        baseline_points,
        precisions=precisions,
        timeout=height_timeout,
        stack_bytes=stack_bytes,
    )
    baseline_rank = stable_height_rank(baseline_runs)

    raw_points, search_wall, pari_milliseconds = bounded_quartic_points(
        parameter_t,
        height_bound=height_bound,
        timeout=search_timeout,
        stack_bytes=stack_bytes,
    )
    signless = signless_quartic_points(raw_points)
    known_x = {point[0] for point in affine_sections}
    seen_jacobian_x = {point[0] for point in baseline_points}
    new_records: list[dict[str, Any]] = []
    new_images: list[tuple[Fraction, Fraction]] = []
    zero_ordinate_points: list[dict[str, str]] = []
    for quartic_point in signless:
        if quartic_point[0] in known_x:
            continue
        if quartic_point[1] == 0:
            zero_ordinate_points.append(
                {
                    "quartic_x": rational_to_string(quartic_point[0]),
                    "quartic_z": "0",
                }
            )
            continue
        jacobian_point = quartic_point_to_short_jacobian(
            RANK13_CONSTRUCTION, parameter_t, quartic_point
        )
        if not point_on_short_curve(coefficients, jacobian_point):
            raise AssertionError("a bounded point missed the exact Jacobian")
        if jacobian_point[0] in seen_jacobian_x:
            continue
        seen_jacobian_x.add(jacobian_point[0])
        new_records.append(exact_jacobian_record(quartic_point, jacobian_point))
        new_images.append(jacobian_point)

    pool_points = baseline_points + tuple(new_images)
    pool_provenance = baseline_provenance + tuple(
        f"bounded-search-{index}" for index in range(1, len(new_images) + 1)
    )
    pool_runs = height_matrix_replay(
        coefficients,
        pool_points,
        precisions=precisions,
        timeout=height_timeout,
        stack_bytes=stack_bytes,
    )
    pool_rank = stable_height_rank(pool_runs)
    subset_indices = pool_runs[-1]["subset_indices_one_based"]
    selected_points = tuple(pool_points[index - 1] for index in subset_indices)
    selected_records = [
        {
            "pool_index_one_based": index,
            "provenance": pool_provenance[index - 1],
            "jacobian_x": rational_to_string(pool_points[index - 1][0]),
            "jacobian_y": rational_to_string(pool_points[index - 1][1]),
            "exact_jacobian_membership_checked": True,
        }
        for index in subset_indices
    ]

    result = {
        "candidate_id": candidate.identifier,
        "parameter_u": candidate.parameter_u,
        "parameter_t": rational_to_string(parameter_t),
        "source_artifact": candidate.source_artifact,
        "score": candidate.score,
        "score_last_numerical_prime": candidate.last_numerical_prime,
        "log_conductor": candidate.log_conductor,
        "below_strict_log_conductor_target": candidate.below_target,
        "root_number": candidate.root_number,
        "exact_sections": {
            "affine_sections_checked": len(affine_sections),
            "distinct_affine_sections": len(set(affine_sections)),
            "exact_short_jacobian_images_checked": len(affine_images),
            "affine_section_sha256": point_digest(affine_sections),
            "affine_image_sha256": point_digest(affine_images),
            "split_infinity_point": {
                "jacobian_x": rational_to_string(infinity[0]),
                "jacobian_y": rational_to_string(infinity[1]),
                "exact_jacobian_membership_checked": True,
            },
        },
        "baseline_height_matrix_runs": list(baseline_runs),
        "stable_baseline_numerical_rank": baseline_rank,
        "bounded_search": {
            "quartic_naive_height_bound": height_bound,
            "scope": (
                "PARI hyperellratpoints affine rational points on the exact "
                "primitive quartic up to the declared naive height"
            ),
            "wall_seconds": search_wall,
            "pari_reported_milliseconds": pari_milliseconds,
            "signed_points_found": len(raw_points),
            "distinct_quartic_x_values": len(signless),
            "displayed_affine_section_x_values_found": sum(
                point[0] in known_x for point in signless
            ),
            "new_distinct_jacobian_images": len(new_images),
            "zero_ordinate_nonvisible_points_not_mapped": zero_ordinate_points,
            "all_nonzero_new_points_mapped_and_checked_exactly": True,
            "new_point_records": new_records,
        },
        "augmented_pool_point_count": len(pool_points),
        "augmented_height_matrix_runs": list(pool_runs),
        "stable_augmented_pool_numerical_rank": pool_rank,
        "numerical_rank_gain_over_baseline": pool_rank - baseline_rank,
        "explicit_numerically_independent_subset": selected_records,
        "height_status": (
            "stable across declared precisions and exact point membership; "
            "numerical independence evidence, not an exact certificate"
        ),
    }
    result["frontier_search_height_bound"] = height_bound
    result["frontier_stable_pool_numerical_rank"] = pool_rank
    result["frontier_explicit_numerically_independent_subset"] = selected_records

    if escalated_height_bound is not None:
        if escalated_height_bound <= height_bound:
            raise ValueError("the escalated height bound must exceed the base bound")
        if escalation_timeout is None or escalation_timeout <= 0:
            raise ValueError("an escalated search requires a positive timeout")
        raw_escalated, escalation_wall, escalation_milliseconds = (
            bounded_quartic_points(
                parameter_t,
                height_bound=escalated_height_bound,
                timeout=escalation_timeout,
                stack_bytes=stack_bytes,
            )
        )
        signless_escalated = signless_quartic_points(raw_escalated)
        escalated_seen_x = {point[0] for point in baseline_points}
        escalated_records: list[dict[str, Any]] = []
        escalated_images: list[tuple[Fraction, Fraction]] = []
        escalated_zero_ordinates: list[dict[str, str]] = []
        for quartic_point in signless_escalated:
            if quartic_point[0] in known_x:
                continue
            if quartic_point[1] == 0:
                escalated_zero_ordinates.append(
                    {
                        "quartic_x": rational_to_string(quartic_point[0]),
                        "quartic_z": "0",
                    }
                )
                continue
            jacobian_point = quartic_point_to_short_jacobian(
                RANK13_CONSTRUCTION, parameter_t, quartic_point
            )
            if not point_on_short_curve(coefficients, jacobian_point):
                raise AssertionError("an escalated point missed the exact Jacobian")
            if jacobian_point[0] in escalated_seen_x:
                continue
            escalated_seen_x.add(jacobian_point[0])
            escalated_records.append(
                exact_jacobian_record(quartic_point, jacobian_point)
            )
            escalated_images.append(jacobian_point)

        escalated_pool = baseline_points + tuple(escalated_images)
        escalated_provenance = baseline_provenance + tuple(
            f"escalated-bounded-search-{index}"
            for index in range(1, len(escalated_images) + 1)
        )
        escalated_runs = height_matrix_replay(
            coefficients,
            escalated_pool,
            precisions=precisions,
            timeout=height_timeout,
            stack_bytes=stack_bytes,
        )
        escalated_rank = stable_height_rank(escalated_runs)
        escalated_indices = escalated_runs[-1]["subset_indices_one_based"]
        escalated_selected_points = tuple(
            escalated_pool[index - 1] for index in escalated_indices
        )
        escalated_selected_records = [
            {
                "pool_index_one_based": index,
                "provenance": escalated_provenance[index - 1],
                "jacobian_x": rational_to_string(escalated_pool[index - 1][0]),
                "jacobian_y": rational_to_string(escalated_pool[index - 1][1]),
                "exact_jacobian_membership_checked": True,
            }
            for index in escalated_indices
        ]
        result["escalated_bounded_search"] = {
            "quartic_naive_height_bound": escalated_height_bound,
            "timeout_seconds": escalation_timeout,
            "scope": (
                "PARI hyperellratpoints affine rational points on the exact "
                "primitive quartic up to the escalated naive height"
            ),
            "wall_seconds": escalation_wall,
            "pari_reported_milliseconds": escalation_milliseconds,
            "signed_points_found": len(raw_escalated),
            "distinct_quartic_x_values": len(signless_escalated),
            "displayed_affine_section_x_values_found": sum(
                point[0] in known_x for point in signless_escalated
            ),
            "new_distinct_jacobian_images": len(escalated_images),
            "zero_ordinate_nonvisible_points_not_mapped": escalated_zero_ordinates,
            "all_nonzero_new_points_mapped_and_checked_exactly": True,
            "new_point_records": escalated_records,
        }
        result["escalated_augmented_pool_point_count"] = len(escalated_pool)
        result["escalated_augmented_height_matrix_runs"] = list(escalated_runs)
        result["stable_escalated_pool_numerical_rank"] = escalated_rank
        result["escalated_numerical_rank_gain_over_baseline"] = (
            escalated_rank - baseline_rank
        )
        result["escalated_explicit_numerically_independent_subset"] = (
            escalated_selected_records
        )
        result["frontier_search_height_bound"] = escalated_height_bound
        result["frontier_stable_pool_numerical_rank"] = escalated_rank
        result["frontier_explicit_numerically_independent_subset"] = (
            escalated_selected_records
        )
        selected_points = escalated_selected_points
    return result, selected_points


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_nagao_rank13_integer_u.json",
    )
    parser.add_argument(
        "--additional-input",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_nagao_rank13_integer_u2000.json",
    )
    parser.add_argument("--additional-u", type=parse_u_values, default=(1256,))
    parser.add_argument("--height-bound", type=int, default=50_000)
    parser.add_argument("--escalate-u", type=parse_u_values, default=DEFAULT_ESCALATE_U)
    parser.add_argument("--escalated-height-bound", type=int, default=1_000_000)
    parser.add_argument("--escalation-timeout", type=float, default=60.0)
    parser.add_argument("--precisions", type=parse_precisions, default=(72, 120))
    parser.add_argument("--search-timeout", type=float, default=30.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--rank-u", type=parse_u_values, default=DEFAULT_RANK_U)
    parser.add_argument("--rank-timeout", type=float, default=20.0)
    parser.add_argument(
        "--saturation-u", type=parse_u_values, default=DEFAULT_SATURATION_U
    )
    parser.add_argument("--saturation-bound", type=int, default=20)
    parser.add_argument("--saturation-timeout", type=float, default=5.0)
    parser.add_argument(
        "--analytic-u", type=parse_u_values, default=DEFAULT_ANALYTIC_U
    )
    parser.add_argument("--analytic-timeout", type=float, default=3.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--rank-stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_nagao_rank13_finalist_triage.json",
    )
    return parser


def _validate_requested_u(
    finalists: Sequence[Finalist], requested: Sequence[int], label: str
) -> None:
    available = {candidate.parameter_u for candidate in finalists}
    absent = [parameter_u for parameter_u in requested if parameter_u not in available]
    if absent:
        raise SystemExit(f"{label} contains absent finalist u values: {absent}")


def main() -> None:
    args = build_parser().parse_args()
    if args.height_bound <= 0:
        raise SystemExit("--height-bound must be positive")
    if args.escalated_height_bound <= args.height_bound:
        raise SystemExit("--escalated-height-bound must exceed --height-bound")
    if args.escalation_timeout <= 0 or args.escalation_timeout > 60:
        raise SystemExit("--escalation-timeout must be in (0,60]")
    if min(args.search_timeout, args.height_timeout, args.rank_timeout) <= 0:
        raise SystemExit("search, height, and rank timeouts must be positive")
    if args.rank_timeout > 30:
        raise SystemExit("--rank-timeout is intentionally capped at 30 seconds")
    if min(args.saturation_timeout, args.analytic_timeout) <= 0:
        raise SystemExit("saturation and analytic timeouts must be positive")
    if args.saturation_bound < 3:
        raise SystemExit("--saturation-bound must be at least 3")
    if min(args.stack_bytes, args.rank_stack_bytes) < 8_000_000:
        raise SystemExit("PARI stack bounds must be at least 8,000,000 bytes")

    finalists = load_finalists(args.input, args.additional_input, args.additional_u)
    _validate_requested_u(finalists, args.rank_u, "--rank-u")
    _validate_requested_u(finalists, args.saturation_u, "--saturation-u")
    _validate_requested_u(finalists, args.analytic_u, "--analytic-u")
    _validate_requested_u(finalists, args.escalate_u, "--escalate-u")

    results: list[dict[str, Any]] = []
    selected_by_u: dict[int, tuple[tuple[Fraction, Fraction], ...]] = {}
    coefficients_by_u: dict[int, tuple[Fraction, ...]] = {}
    for candidate in finalists:
        result, selected_points = exact_candidate_triage(
            candidate,
            height_bound=args.height_bound,
            precisions=args.precisions,
            search_timeout=args.search_timeout,
            height_timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
            escalated_height_bound=(
                args.escalated_height_bound
                if candidate.parameter_u in args.escalate_u
                else None
            ),
            escalation_timeout=args.escalation_timeout,
        )
        results.append(result)
        selected_by_u[candidate.parameter_u] = selected_points
        coefficients_by_u[candidate.parameter_u] = (
            rank13_base_changed_short_jacobian_coefficients(Q(candidate.parameter_u))
        )
        print(
            f"u={candidate.parameter_u} logN={candidate.log_conductor} "
            f"baseline={result['stable_baseline_numerical_rank']} "
            f"pool={result['frontier_stable_pool_numerical_rank']} "
            f"new={result['bounded_search']['new_distinct_jacobian_images']}",
            flush=True,
        )

    by_u = {int(result["parameter_u"]): result for result in results}
    for parameter_u in args.saturation_u:
        by_u[parameter_u]["small_prime_saturation_probe"] = saturation_probe(
            coefficients_by_u[parameter_u],
            selected_by_u[parameter_u],
            prime_bound=args.saturation_bound,
            timeout=args.saturation_timeout,
            stack_bytes=args.stack_bytes,
        )
    for parameter_u in args.rank_u:
        by_u[parameter_u]["pari_ellrank_effort_zero"] = ellrank_probe(
            coefficients_by_u[parameter_u],
            selected_by_u[parameter_u],
            timeout=args.rank_timeout,
            stack_bytes=args.rank_stack_bytes,
        )
        probe = by_u[parameter_u]["pari_ellrank_effort_zero"]
        print(
            f"u={parameter_u} ellrank={probe['status']} "
            f"bounds={probe.get('lower_bound')},{probe.get('upper_bound')}",
            flush=True,
        )
    for parameter_u in args.analytic_u:
        by_u[parameter_u]["numerical_analytic_rank_probe"] = analytic_rank_probe(
            coefficients_by_u[parameter_u],
            timeout=args.analytic_timeout,
            stack_bytes=args.rank_stack_bytes,
        )

    results.sort(
        key=lambda record: (
            -int(record["frontier_stable_pool_numerical_rank"]),
            Decimal(record["log_conductor"]),
            int(record["parameter_u"]),
        )
    )
    equal_bounds = []
    for result in results:
        probe = result.get("pari_ellrank_effort_zero", {})
        if probe.get("equal_returned_bounds"):
            equal_bounds.append(
                {
                    "parameter_u": result["parameter_u"],
                    "lower_bound": probe["lower_bound"],
                    "upper_bound": probe["upper_bound"],
                }
            )
    max_numerical_rank = max(
        int(result["frontier_stable_pool_numerical_rank"]) for result in results
    )
    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    script_path = Path(__file__).resolve()
    input_paths = [args.input]
    if args.additional_u:
        input_paths.append(args.additional_input)
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded exact-point and numerical-height triage; no rank-21 or "
            "rank-30 target hit is certified"
        ),
        "primary_source": PRIMARY_SOURCE,
        "inputs": [
            {
                "path": str(path),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for path in input_paths
        ],
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "certified_hits": [],
        },
        "method": {
            "exact_affine_sections_per_candidate": 13,
            "split_infinity_covariant_limit_added": True,
            "quartic_naive_height_bound": args.height_bound,
            "escalated_u": list(args.escalate_u),
            "escalated_quartic_naive_height_bound": args.escalated_height_bound,
            "escalation_timeout_seconds_per_candidate": args.escalation_timeout,
            "height_decimal_precisions": list(args.precisions),
            "ellrank_effort": 0,
            "ellrank_timeout_seconds_per_candidate": args.rank_timeout,
            "rank_probe_u": list(args.rank_u),
            "saturation_u": list(args.saturation_u),
            "saturation_prime_bound": args.saturation_bound,
            "analytic_probe_u": list(args.analytic_u),
        },
        "summary": {
            "candidate_count": len(results),
            "below_strict_log_conductor_target": sum(
                bool(result["below_strict_log_conductor_target"])
                for result in results
            ),
            "maximum_frontier_stable_pool_numerical_rank": max_numerical_rank,
            "candidates_at_maximum_numerical_rank": [
                int(result["parameter_u"])
                for result in results
                if int(result["frontier_stable_pool_numerical_rank"])
                == max_numerical_rank
            ],
            "pari_equal_rank_bounds": equal_bounds,
            "interpretation": (
                "exact point checks plus precision-stable numerical height rank; "
                "only explicitly stored PARI ellrank intervals are rank bounds"
            ),
        },
        "candidates": results,
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "parameters": {
            "search_timeout_seconds_per_candidate": args.search_timeout,
            "height_timeout_seconds_per_replay": args.height_timeout,
            "saturation_timeout_seconds_per_candidate": args.saturation_timeout,
            "analytic_timeout_seconds_per_candidate": args.analytic_timeout,
            "pari_stack_bytes": args.stack_bytes,
            "pari_rank_stack_bytes": args.rank_stack_bytes,
            "output": str(args.output),
        },
        "reproducing_command": command,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
