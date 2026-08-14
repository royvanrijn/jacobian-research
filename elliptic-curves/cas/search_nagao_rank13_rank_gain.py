#!/usr/bin/env python3
"""Rank-gain-driven rare-event search in Nagao's rank-13 base change.

This search deliberately does *not* select candidates by a Nagao prime score.
It starts with a deterministic rational ``u`` population, enumerates exact
points on the primitive quartic in increasing naive-height boxes, and spends
height-pairing and conductor work only on survivors of the preceding box.

Five linear sections which were absent from the original repository baseline
are treated as known dependent companions.  Their abscissae are removed from
the point-yield label before selection.  Thus the cheap label is the count of
exact, nonzero, non-displayed, non-companion quartic points, and the expensive
label is the precision-stable numerical Neron--Tate rank of their Jacobian
images together with the thirteen displayed sections and the split point at
infinity.

The final stage is fixed before its results are observed.  Conductors are
computed only after final-stage selection, so conductor data cannot leak into
candidate selection.  Exact point membership and bounded enumeration are
exact computations.  Height-matrix ranks remain numerical evidence and are
never reported as Mordell--Weil rank certificates.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from math import gcd, isqrt
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

from ek_k3 import rational_to_string
from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK13_BASE_CHANGE_CONSTANT,
    RANK13_CONSTRUCTION,
    primitive_quartic_coefficients,
    quartic_point_to_short_jacobian,
    quartic_value,
    rank13_base_changed_short_jacobian_coefficients,
    rank13_base_parameter,
    rank13_known_quartic_points,
)
from pari_bridge import pari_version
from search_extra_points import gp_rational, gp_vector, parse_point_vector
from triage_nagao_rank13_finalists import (
    point_digest,
    point_on_short_curve,
    split_infinity_jacobian_point,
)


Q = Fraction
TARGET_LOG_CONDUCTOR = Decimal("182.72")
DEFAULT_CENTERS = (Q(42), Q(118))
DEFAULT_ANCHORS = (Q(42), Q(84), Q(118))


@dataclass(frozen=True)
class ParameterCandidate:
    parameter_u: Fraction
    origins: tuple[str, ...]

    @property
    def identifier(self) -> str:
        return f"u-{self.parameter_u.numerator}-{self.parameter_u.denominator}"

    @property
    def parameter_t(self) -> Fraction:
        return rank13_base_parameter(self.parameter_u)


@dataclass(frozen=True)
class PointScreen:
    candidate: ParameterCandidate
    raw_points: tuple[tuple[Fraction, Fraction], ...]
    unexpected_points: tuple[tuple[Fraction, Fraction], ...]
    record: dict[str, Any]


def parse_rationals(value: str) -> tuple[Fraction, ...]:
    try:
        values = tuple(Q(part) for part in value.split(",") if part)
    except (ValueError, ZeroDivisionError) as error:
        raise argparse.ArgumentTypeError("values must be rational numbers") from error
    if not values or any(item <= 0 for item in values):
        raise argparse.ArgumentTypeError("provide distinct positive rationals")
    if len(set(values)) != len(values):
        raise argparse.ArgumentTypeError("rational values must be distinct")
    return values


def parse_precisions(value: str) -> tuple[int, ...]:
    try:
        values = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("precisions must be integers") from error
    if len(values) < 2 or values != tuple(sorted(set(values))):
        raise argparse.ArgumentTypeError(
            "provide at least two strictly increasing precisions"
        )
    if values[0] < 32:
        raise argparse.ArgumentTypeError("height precision must be at least 32")
    return values


def canonical_positive_u(parameter_u: Fraction) -> Fraction:
    """Quotient the exact ``u <-> 23550/u`` symmetry of the Jacobian.

    The two parameters give opposite ``T`` values.  Nagao's primitive quartic
    and its short Jacobian are even in ``T``.  We retain the representative
    below ``sqrt(23550)`` to avoid evaluating the same curve twice.
    """

    parameter_u = Q(parameter_u)
    if parameter_u <= 0:
        raise ValueError("u must be positive")
    if parameter_u**2 > RANK13_BASE_CHANGE_CONSTANT:
        parameter_u = Q(RANK13_BASE_CHANGE_CONSTANT, 1) / parameter_u
    return parameter_u


def generate_population(
    *,
    farey_denominator: int,
    mutation_denominator: int,
    mutation_numerator_radius: int,
    centers: Sequence[Fraction],
    anchors: Sequence[Fraction],
) -> tuple[ParameterCandidate, ...]:
    """Return a deterministic rational population with provenance labels."""

    if farey_denominator < 1 or mutation_denominator < 1:
        raise ValueError("denominator bounds must be positive")
    if mutation_numerator_radius < 0:
        raise ValueError("the mutation radius must be nonnegative")
    origins: dict[Fraction, set[str]] = {}

    def add(parameter_u: Fraction, origin: str) -> None:
        canonical = canonical_positive_u(parameter_u)
        origins.setdefault(canonical, set()).add(origin)

    # All reduced positive rationals n/d strictly below sqrt(23550).
    for denominator in range(1, farey_denominator + 1):
        maximum_numerator = isqrt(
            RANK13_BASE_CHANGE_CONSTANT * denominator**2 - 1
        )
        for numerator in range(1, maximum_numerator + 1):
            if gcd(numerator, denominator) == 1:
                add(Q(numerator, denominator), "farey-grid")

    # High-denominator local mutations around the two existing frontiers.
    for center in centers:
        center = Q(center)
        for denominator in range(1, mutation_denominator + 1):
            for numerator_delta in range(
                -mutation_numerator_radius, mutation_numerator_radius + 1
            ):
                mutated = center + Q(numerator_delta, denominator)
                if mutated > 0:
                    add(mutated, f"mutation-around-{rational_to_string(center)}")

    for anchor in anchors:
        add(Q(anchor), "forced-anchor")

    return tuple(
        ParameterCandidate(parameter_u, tuple(sorted(labels)))
        for parameter_u, labels in sorted(origins.items())
    )


def companion_section_x_values(parameter_t: Fraction) -> tuple[Fraction, ...]:
    """Return five omitted dependent linear-section abscissae.

    Exact substitution in Nagao's primitive quartic verifies these whenever a
    corresponding point is encountered.  They are excluded from the
    rare-event point-yield label because they are systematic, not exceptional.
    """

    parameter_t = Q(parameter_t)
    return (
        (-parameter_t + 703) / 15,
        7 * parameter_t / 15 + Q(928, 15),
        -7 * parameter_t / 15 + Q(928, 15),
        5 * parameter_t / 3 + Q(3628, 15),
        -5 * parameter_t / 3 + Q(3628, 15),
    )


def displayed_section_x_values(parameter_t: Fraction) -> tuple[Fraction, ...]:
    parameter_t = Q(parameter_t)
    values = tuple(
        root + sign * parameter_t
        for root in RANK13_CONSTRUCTION.roots
        for sign in (-1, 1)
    ) + ((parameter_t + 703) / 15,)
    return values


def _terminate_process_group(process: subprocess.Popen[str]) -> None:
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.communicate(timeout=0.75)
        return
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        return
    process.communicate()


def run_gp_capped(
    program: str, *, timeout: float, stack_bytes: int
) -> tuple[str, float]:
    """Run one foreground GP process with timeout and process-group cleanup."""

    if timeout <= 0 or stack_bytes < 8_000_000:
        raise ValueError("invalid GP process bounds")
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    started = time.monotonic()
    process = subprocess.Popen(
        [executable, "-q", "-s", str(stack_bytes)],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    try:
        stdout, stderr = process.communicate(program, timeout=timeout)
    except subprocess.TimeoutExpired:
        _terminate_process_group(process)
        raise
    finally:
        _terminate_process_group(process)
    elapsed = time.monotonic() - started
    if process.returncode != 0 or "***" in stderr:
        raise RuntimeError(f"PARI/GP failed: {' '.join(stderr.split())[:1000]}")
    return stdout, elapsed


def quartic_gp_polynomial(coefficients: Sequence[Fraction]) -> str:
    return "+".join(
        f"{gp_rational(Q(coefficient))}*x^{power}"
        for power, coefficient in enumerate(coefficients)
    )


def _screen_record(
    candidate: ParameterCandidate,
    raw_points: tuple[tuple[Fraction, Fraction], ...],
    *,
    height_bound: int,
) -> PointScreen:
    parameter_t = candidate.parameter_t
    coefficients = primitive_quartic_coefficients(RANK13_CONSTRUCTION, parameter_t)
    for x_value, z_value in raw_points:
        if z_value**2 != quartic_value(coefficients, x_value):
            raise AssertionError("PARI returned a point off the exact quartic")

    signless: list[tuple[Fraction, Fraction]] = []
    seen_x: set[Fraction] = set()
    for point in raw_points:
        if point[0] not in seen_x:
            signless.append(point)
            seen_x.add(point[0])

    displayed = set(displayed_section_x_values(parameter_t))
    companions = set(companion_section_x_values(parameter_t))
    unexpected = tuple(
        point
        for point in signless
        if point[1] != 0 and point[0] not in displayed and point[0] not in companions
    )
    record = {
        "candidate_id": candidate.identifier,
        "parameter_u": rational_to_string(candidate.parameter_u),
        "parameter_t": rational_to_string(parameter_t),
        "origins": list(candidate.origins),
        "quartic_naive_height_bound": height_bound,
        "signed_points_found": len(raw_points),
        "distinct_quartic_x_values": len(signless),
        "displayed_section_x_values_found": sum(
            point[0] in displayed for point in signless
        ),
        "dependent_companion_x_values_found": sum(
            point[0] in companions for point in signless
        ),
        "zero_ordinate_points_excluded": sum(point[1] == 0 for point in signless),
        "unexpected_nonzero_quartic_x_values": len(unexpected),
        "all_returned_points_checked_exactly": True,
        "unexpected_point_sha256": point_digest(unexpected),
    }
    return PointScreen(candidate, raw_points, unexpected, record)


def batch_point_screen(
    candidates: Sequence[ParameterCandidate],
    *,
    height_bound: int,
    batch_size: int,
    timeout_per_batch: float,
    stack_bytes: int,
) -> tuple[PointScreen, ...]:
    if height_bound <= 0 or batch_size <= 0:
        raise ValueError("height and batch bounds must be positive")
    answer: list[PointScreen] = []
    for start in range(0, len(candidates), batch_size):
        batch = candidates[start : start + batch_size]
        commands: list[str] = []
        for offset, candidate in enumerate(batch):
            coefficients = primitive_quartic_coefficients(
                RANK13_CONSTRUCTION, candidate.parameter_t
            )
            commands.extend(
                (
                    f"Q={quartic_gp_polynomial(coefficients)};",
                    f"R=hyperellratpoints(Q,{height_bound});",
                    f'print("ROW|{offset}|",R);',
                )
            )
        commands.append("quit")
        output, _ = run_gp_capped(
            "\n".join(commands) + "\n",
            timeout=timeout_per_batch,
            stack_bytes=stack_bytes,
        )
        parsed: dict[int, tuple[tuple[Fraction, Fraction], ...]] = {}
        for line in output.splitlines():
            if not line.startswith("ROW|"):
                continue
            _, offset_text, point_text = line.split("|", 2)
            parsed[int(offset_text)] = parse_point_vector(point_text)
        if set(parsed) != set(range(len(batch))):
            raise AssertionError("PARI omitted one or more point-screen rows")
        for offset, candidate in enumerate(batch):
            answer.append(
                _screen_record(
                    candidate, parsed[offset], height_bound=height_bound
                )
            )
    return tuple(answer)


def _candidate_sort_key(screen: PointScreen) -> tuple[int, int, Fraction]:
    return (
        -int(screen.record["unexpected_nonzero_quartic_x_values"]),
        -int(screen.record["distinct_quartic_x_values"]),
        screen.candidate.parameter_u,
    )


def select_point_survivors(
    screens: Sequence[PointScreen],
    *,
    keep_count: int,
    anchors: Sequence[Fraction],
    denominator_diversity: int,
) -> tuple[ParameterCandidate, ...]:
    """Select by declared point data, then force anchors and denominator strata."""

    if keep_count <= 0 or denominator_diversity < 0:
        raise ValueError("invalid survivor counts")
    ordered = sorted(screens, key=_candidate_sort_key)
    selected = {screen.candidate.parameter_u for screen in ordered[:keep_count]}
    by_denominator: dict[int, int] = {}
    for screen in ordered:
        denominator = screen.candidate.parameter_u.denominator
        if by_denominator.get(denominator, 0) >= denominator_diversity:
            continue
        if int(screen.record["unexpected_nonzero_quartic_x_values"]) == 0:
            continue
        selected.add(screen.candidate.parameter_u)
        by_denominator[denominator] = by_denominator.get(denominator, 0) + 1
    selected.update(canonical_positive_u(anchor) for anchor in anchors)
    by_u = {screen.candidate.parameter_u: screen.candidate for screen in screens}
    missing = selected - set(by_u)
    if missing:
        raise AssertionError(f"forced survivor absent from population: {missing}")
    return tuple(by_u[parameter_u] for parameter_u in sorted(selected))


def _parse_vecsmall(text: str) -> list[int]:
    match = re.search(r"Vecsmall\(\[(.*?)\]\)", text)
    if match is None:
        raise AssertionError("PARI omitted a subset index vector")
    if not match.group(1):
        return []
    return [int(part) for part in match.group(1).split(",")]


def height_rank_replay(
    coefficients: Sequence[Fraction],
    baseline: Sequence[tuple[Fraction, Fraction]],
    pool: Sequence[tuple[Fraction, Fraction]],
    *,
    precisions: Sequence[int],
    timeout: float,
    stack_bytes: int,
) -> tuple[dict[str, Any], ...]:
    if not baseline or len(pool) < len(baseline):
        raise ValueError("invalid height-replay point pools")
    if any(not point_on_short_curve(coefficients, point) for point in pool):
        raise AssertionError("an exact point missed the Jacobian")
    curve = ",".join(gp_rational(Q(value)) for value in coefficients)
    baseline_vector = ",".join(gp_vector(point) for point in baseline)
    pool_vector = ",".join(gp_vector(point) for point in pool)
    commands = [
        f"E=ellinit([{curve}]);",
        f"B=[{baseline_vector}];",
        f"P=[{pool_vector}];",
    ]
    for precision in precisions:
        commands.extend(
            (
                f"default(realprecision,{precision});",
                "HB=ellheightmatrix(E,B);HP=ellheightmatrix(E,P);",
                "IX=matindexrank(HP);K=vecextract(P,IX[2]);HK=ellheightmatrix(E,K);",
                f'print("HEIGHT_{precision}_BEGIN");',
                "print(matrank(HB));print(matrank(HP));print(IX[2]);",
                "print(matdet(HK));",
                "EV=mateigen(HK,1)[1];print(vecmin(EV));print(vecmax(EV));",
                f'print("HEIGHT_{precision}_END");',
            )
        )
    commands.append("quit")
    output, _ = run_gp_capped(
        "\n".join(commands) + "\n", timeout=timeout, stack_bytes=stack_bytes
    )
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    records: list[dict[str, Any]] = []
    for precision in precisions:
        start = lines.index(f"HEIGHT_{precision}_BEGIN") + 1
        end = lines.index(f"HEIGHT_{precision}_END")
        values = lines[start:end]
        records.append(
            {
                "decimal_precision": precision,
                "baseline_numerical_rank": int(values[0]),
                "pool_numerical_rank": int(values[1]),
                "subset_indices_one_based": _parse_vecsmall(values[2]),
                "subset_height_determinant": values[3],
                "subset_smallest_eigenvalue": values[4],
                "subset_largest_eigenvalue": values[5],
            }
        )
    return tuple(records)


def evaluate_height_rank(
    screen: PointScreen,
    *,
    precisions: Sequence[int],
    timeout: float,
    stack_bytes: int,
    store_selected_points: bool,
) -> tuple[dict[str, Any], tuple[tuple[Fraction, Fraction], ...]]:
    candidate = screen.candidate
    parameter_t = candidate.parameter_t
    coefficients = rank13_base_changed_short_jacobian_coefficients(
        candidate.parameter_u
    )
    sections = rank13_known_quartic_points(parameter_t)
    baseline = tuple(
        quartic_point_to_short_jacobian(RANK13_CONSTRUCTION, parameter_t, point)
        for point in sections
    ) + (split_infinity_jacobian_point(candidate.parameter_u),)
    seen_jacobian_x = {point[0] for point in baseline}
    new_images: list[tuple[Fraction, Fraction]] = []
    for quartic_point in screen.unexpected_points:
        jacobian_point = quartic_point_to_short_jacobian(
            RANK13_CONSTRUCTION, parameter_t, quartic_point
        )
        if not point_on_short_curve(coefficients, jacobian_point):
            raise AssertionError("an unexpected point missed the exact Jacobian")
        if jacobian_point[0] not in seen_jacobian_x:
            seen_jacobian_x.add(jacobian_point[0])
            new_images.append(jacobian_point)
    pool = baseline + tuple(new_images)
    runs = height_rank_replay(
        coefficients,
        baseline,
        pool,
        precisions=precisions,
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    stable = (
        len({run["baseline_numerical_rank"] for run in runs}) == 1
        and len({run["pool_numerical_rank"] for run in runs}) == 1
        and len({tuple(run["subset_indices_one_based"]) for run in runs}) == 1
    )
    final_run = runs[-1]
    indices = final_run["subset_indices_one_based"]
    selected = tuple(pool[index - 1] for index in indices)
    record = {
        **screen.record,
        "new_distinct_jacobian_images": len(new_images),
        "new_jacobian_image_sha256": point_digest(new_images),
        "height_matrix_runs": list(runs),
        "height_rank_stable_across_precisions": stable,
        "stable_baseline_numerical_rank": (
            int(final_run["baseline_numerical_rank"]) if stable else None
        ),
        "stable_pool_numerical_rank": (
            int(final_run["pool_numerical_rank"]) if stable else None
        ),
        "stable_numerical_rank_gain": (
            int(final_run["pool_numerical_rank"])
            - int(final_run["baseline_numerical_rank"])
            if stable
            else None
        ),
        "exact_pool_point_count": len(pool),
        "interpretation": (
            "all point memberships are exact; height rank is numerical evidence"
        ),
    }
    if store_selected_points:
        record["explicit_numerically_independent_subset"] = [
            {
                "pool_index_one_based": index,
                "jacobian_x": rational_to_string(pool[index - 1][0]),
                "jacobian_y": rational_to_string(pool[index - 1][1]),
                "exact_jacobian_membership_checked": True,
            }
            for index in indices
        ]
    else:
        record["selected_subset_sha256"] = point_digest(selected)
    return record, selected


def _height_sort_key(record: dict[str, Any]) -> tuple[int, int, int, Fraction]:
    rank = record.get("stable_pool_numerical_rank")
    gain = record.get("stable_numerical_rank_gain")
    return (
        -(int(rank) if rank is not None else -1),
        -(int(gain) if gain is not None else -1),
        -int(record["unexpected_nonzero_quartic_x_values"]),
        Q(record["parameter_u"]),
    )


def select_height_survivors(
    records: Sequence[dict[str, Any]],
    candidates: Sequence[ParameterCandidate],
    *,
    keep_count: int,
    anchors: Sequence[Fraction],
) -> tuple[ParameterCandidate, ...]:
    if keep_count <= 0:
        raise ValueError("keep_count must be positive")
    ordered = sorted(records, key=_height_sort_key)
    selected = {Q(record["parameter_u"]) for record in ordered[:keep_count]}
    selected.update(canonical_positive_u(anchor) for anchor in anchors)
    by_u = {candidate.parameter_u: candidate for candidate in candidates}
    missing = selected - set(by_u)
    if missing:
        raise AssertionError(f"height survivor absent from candidates: {missing}")
    return tuple(by_u[parameter_u] for parameter_u in sorted(selected))


def conductor_probe(
    coefficients: Sequence[Fraction], *, timeout: float, stack_bytes: int
) -> dict[str, Any]:
    curve = ",".join(gp_rational(Q(value)) for value in coefficients)
    program = "\n".join(
        (
            "default(realprecision,80);",
            f"E=ellminimalmodel(ellinit([{curve}]));G=ellglobalred(E);",
            'print("MODEL ",E.a1,"|",E.a2,"|",E.a3,"|",E.a4,"|",E.a6);',
            'print("CONDUCTOR ",G[1]);',
            'print("LOG_CONDUCTOR ",log(G[1]));',
            'print("MINIMAL_DISCRIMINANT ",E.disc);',
            'print("ROOT_NUMBER ",ellrootno(E));',
            "quit",
        )
    ) + "\n"
    started = time.monotonic()
    try:
        output, _ = run_gp_capped(
            program, timeout=timeout, stack_bytes=stack_bytes
        )
    except subprocess.TimeoutExpired:
        return {
            "status": "timeout",
            "timeout_seconds": timeout,
            "wall_seconds": time.monotonic() - started,
        }
    except RuntimeError as error:
        return {
            "status": "pari_error",
            "timeout_seconds": timeout,
            "error": str(error)[:1000],
        }

    def value(label: str) -> str:
        match = re.search(rf"^{label} (.+)$", output, re.MULTILINE)
        if match is None:
            raise AssertionError(f"PARI omitted {label}")
        return match.group(1)

    model = tuple(int(part) for part in value("MODEL").split("|"))
    log_conductor = value("LOG_CONDUCTOR")
    return {
        "status": "completed",
        "minimal_model": model,
        "conductor": int(value("CONDUCTOR")),
        "log_conductor": log_conductor,
        "below_strict_log_conductor_target": (
            Decimal(log_conductor) < TARGET_LOG_CONDUCTOR
        ),
        "minimal_discriminant": int(value("MINIMAL_DISCRIMINANT")),
        "root_number": int(value("ROOT_NUMBER")),
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts" / "generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--farey-denominator", type=int, default=12)
    parser.add_argument("--mutation-denominator", type=int, default=48)
    parser.add_argument("--mutation-numerator-radius", type=int, default=24)
    parser.add_argument("--centers", type=parse_rationals, default=DEFAULT_CENTERS)
    parser.add_argument("--anchors", type=parse_rationals, default=DEFAULT_ANCHORS)
    parser.add_argument("--screen-height", type=int, default=2_000)
    parser.add_argument("--screen-keep", type=int, default=384)
    parser.add_argument("--denominator-diversity", type=int, default=1)
    parser.add_argument("--rank-height", type=int, default=20_000)
    parser.add_argument("--height-keep", type=int, default=192)
    parser.add_argument("--rank-keep", type=int, default=16)
    parser.add_argument("--final-height", type=int, default=250_000)
    parser.add_argument("--escalation-height", type=int, default=1_000_000)
    parser.add_argument("--escalation-keep", type=int, default=4)
    parser.add_argument("--screen-precisions", type=parse_precisions, default=(48, 72))
    parser.add_argument("--final-precisions", type=parse_precisions, default=(72, 120))
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--batch-timeout", type=float, default=20.0)
    parser.add_argument("--height-timeout", type=float, default=20.0)
    parser.add_argument("--final-search-timeout", type=float, default=60.0)
    parser.add_argument("--escalation-search-timeout", type=float, default=90.0)
    parser.add_argument("--conductor-timeout", type=float, default=25.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_nagao_rank13_rank_gain_search.json",
    )
    return parser


def _validate_args(args: argparse.Namespace) -> None:
    positive_integers = (
        args.farey_denominator,
        args.mutation_denominator,
        args.screen_height,
        args.screen_keep,
        args.rank_height,
        args.height_keep,
        args.rank_keep,
        args.final_height,
        args.escalation_height,
        args.escalation_keep,
        args.batch_size,
    )
    if any(value <= 0 for value in positive_integers):
        raise SystemExit("all count, height, and denominator bounds must be positive")
    if args.mutation_numerator_radius < 0 or args.denominator_diversity < 0:
        raise SystemExit("mutation radius and diversity count must be nonnegative")
    if not (
        args.screen_height < args.rank_height < args.final_height
        < args.escalation_height
    ):
        raise SystemExit("search heights must be strictly increasing")
    if not (args.height_keep <= args.screen_keep and args.rank_keep <= args.height_keep):
        raise SystemExit("survivor counts must be nonincreasing")
    timeouts = (
        args.batch_timeout,
        args.height_timeout,
        args.final_search_timeout,
        args.escalation_search_timeout,
        args.conductor_timeout,
    )
    if any(value <= 0 or value > 120 for value in timeouts):
        raise SystemExit("every subprocess timeout must lie in (0,120]")
    if args.stack_bytes < 8_000_000:
        raise SystemExit("the PARI stack bound is too small")


def _compact_screen(screen: PointScreen) -> dict[str, Any]:
    return dict(screen.record)


def main() -> None:
    args = build_parser().parse_args()
    _validate_args(args)
    population = generate_population(
        farey_denominator=args.farey_denominator,
        mutation_denominator=args.mutation_denominator,
        mutation_numerator_radius=args.mutation_numerator_radius,
        centers=args.centers,
        anchors=args.anchors,
    )
    print(f"stage 0 population={len(population)}", flush=True)

    initial_screens = batch_point_screen(
        population,
        height_bound=args.screen_height,
        batch_size=args.batch_size,
        timeout_per_batch=args.batch_timeout,
        stack_bytes=args.stack_bytes,
    )
    point_survivors = select_point_survivors(
        initial_screens,
        keep_count=min(args.screen_keep, len(initial_screens)),
        anchors=args.anchors,
        denominator_diversity=args.denominator_diversity,
    )
    print(f"stage 0 survivors={len(point_survivors)}", flush=True)

    rank_box_screens = batch_point_screen(
        point_survivors,
        height_bound=args.rank_height,
        batch_size=args.batch_size,
        timeout_per_batch=args.batch_timeout,
        stack_bytes=args.stack_bytes,
    )
    height_candidates = select_point_survivors(
        rank_box_screens,
        keep_count=min(args.height_keep, len(rank_box_screens)),
        anchors=args.anchors,
        denominator_diversity=args.denominator_diversity,
    )
    rank_screen_by_u = {
        screen.candidate.parameter_u: screen for screen in rank_box_screens
    }
    height_records: list[dict[str, Any]] = []
    for index, candidate in enumerate(height_candidates, 1):
        record, _ = evaluate_height_rank(
            rank_screen_by_u[candidate.parameter_u],
            precisions=args.screen_precisions,
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
            store_selected_points=False,
        )
        height_records.append(record)
        print(
            f"height {index}/{len(height_candidates)} u={candidate.parameter_u} "
            f"rank={record['stable_pool_numerical_rank']} "
            f"unexpected={record['unexpected_nonzero_quartic_x_values']}",
            flush=True,
        )

    final_candidates = select_height_survivors(
        height_records,
        height_candidates,
        keep_count=min(args.rank_keep, len(height_records)),
        anchors=args.anchors,
    )
    print(f"final-box candidates={len(final_candidates)}", flush=True)
    final_screens = batch_point_screen(
        final_candidates,
        height_bound=args.final_height,
        batch_size=1,
        timeout_per_batch=args.final_search_timeout,
        stack_bytes=args.stack_bytes,
    )
    final_records: list[dict[str, Any]] = []
    final_selected: dict[Fraction, tuple[tuple[Fraction, Fraction], ...]] = {}
    for screen in final_screens:
        record, selected = evaluate_height_rank(
            screen,
            precisions=args.final_precisions,
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
            store_selected_points=True,
        )
        record["conductor_probe"] = conductor_probe(
            rank13_base_changed_short_jacobian_coefficients(
                screen.candidate.parameter_u
            ),
            timeout=args.conductor_timeout,
            stack_bytes=args.stack_bytes,
        )
        final_records.append(record)
        final_selected[screen.candidate.parameter_u] = selected
        print(
            f"final u={screen.candidate.parameter_u} "
            f"rank={record['stable_pool_numerical_rank']} "
            f"logN={record['conductor_probe'].get('log_conductor')}",
            flush=True,
        )

    escalation_candidates = select_height_survivors(
        final_records,
        final_candidates,
        keep_count=min(args.escalation_keep, len(final_records)),
        anchors=(),
    )
    escalation_screens = batch_point_screen(
        escalation_candidates,
        height_bound=args.escalation_height,
        batch_size=1,
        timeout_per_batch=args.escalation_search_timeout,
        stack_bytes=args.stack_bytes,
    )
    escalation_records: list[dict[str, Any]] = []
    final_by_u = {Q(record["parameter_u"]): record for record in final_records}
    for screen in escalation_screens:
        record, _ = evaluate_height_rank(
            screen,
            precisions=args.final_precisions,
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
            store_selected_points=True,
        )
        record["conductor_probe"] = final_by_u[
            screen.candidate.parameter_u
        ]["conductor_probe"]
        escalation_records.append(record)
        print(
            f"escalated u={screen.candidate.parameter_u} "
            f"rank={record['stable_pool_numerical_rank']}",
            flush=True,
        )

    all_frontier = escalation_records + [
        record
        for record in final_records
        if Q(record["parameter_u"])
        not in {Q(item["parameter_u"]) for item in escalation_records}
    ]
    all_frontier.sort(key=_height_sort_key)
    maximum_numerical_rank = max(
        (
            int(record["stable_pool_numerical_rank"])
            for record in all_frontier
            if record["stable_pool_numerical_rank"] is not None
        ),
        default=0,
    )
    numerical_target_leads = [
        {
            "parameter_u": record["parameter_u"],
            "stable_pool_numerical_rank": record["stable_pool_numerical_rank"],
            "log_conductor": record["conductor_probe"].get("log_conductor"),
            "below_strict_log_conductor_target": record["conductor_probe"].get(
                "below_strict_log_conductor_target"
            ),
        }
        for record in all_frontier
        if record["stable_pool_numerical_rank"] is not None
        and int(record["stable_pool_numerical_rank"]) >= 21
        and record["conductor_probe"].get("below_strict_log_conductor_target")
    ]

    script_path = Path(__file__).resolve()
    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    population_text = "\n".join(
        f"{candidate.identifier}|{','.join(candidate.origins)}"
        for candidate in population
    )
    artifact = {
        "schema_version": 1,
        "status": (
            "leakage-free staged exact-point and numerical-height search; "
            "numerical ranks are not exact Mordell-Weil certificates"
        ),
        "primary_source": PRIMARY_SOURCE,
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "certified_hits": [],
            "numerical_target_leads": numerical_target_leads,
        },
        "method": {
            "selection_objective": (
                "exact unexpected quartic-point yield followed by stable "
                "numerical height rank; no Nagao prime score"
            ),
            "known_dependent_companion_x_formulas": [
                "(-T+703)/15",
                "7*T/15+928/15",
                "-7*T/15+928/15",
                "5*T/3+3628/15",
                "-5*T/3+3628/15",
            ],
            "selection_order": [
                "initial exact point box",
                "rank-box exact point yield",
                "two-precision numerical height rank",
                "final exact point box and two-precision height rank",
                "post-selection conductor computation",
                "declared top-k escalation box",
            ],
            "conductor_not_used_for_selection": True,
            "u_symmetry_quotiented": "u <-> 23550/u, equivalently T <-> -T",
        },
        "population": {
            "count": len(population),
            "sha256": hashlib.sha256(population_text.encode()).hexdigest(),
            "farey_denominator_bound": args.farey_denominator,
            "mutation_centers": [rational_to_string(value) for value in args.centers],
            "mutation_denominator_bound": args.mutation_denominator,
            "mutation_numerator_radius": args.mutation_numerator_radius,
            "forced_anchors": [rational_to_string(value) for value in args.anchors],
        },
        "initial_point_screen": {
            "height_bound": args.screen_height,
            "population_count": len(initial_screens),
            "records": [_compact_screen(screen) for screen in initial_screens],
            "retained_candidate_ids": [
                candidate.identifier for candidate in point_survivors
            ],
        },
        "rank_box_point_screen": {
            "height_bound": args.rank_height,
            "population_count": len(rank_box_screens),
            "records": [_compact_screen(screen) for screen in rank_box_screens],
            "height_evaluated_candidate_ids": [
                candidate.identifier for candidate in height_candidates
            ],
        },
        "rank_box_height_evaluation": {
            "precisions": list(args.screen_precisions),
            "records": sorted(height_records, key=_height_sort_key),
            "retained_candidate_ids": [
                candidate.identifier for candidate in final_candidates
            ],
        },
        "final_box": {
            "height_bound": args.final_height,
            "precisions": list(args.final_precisions),
            "records": sorted(final_records, key=_height_sort_key),
        },
        "escalation_box": {
            "height_bound": args.escalation_height,
            "selection_count": args.escalation_keep,
            "precisions": list(args.final_precisions),
            "records": sorted(escalation_records, key=_height_sort_key),
        },
        "summary": {
            "maximum_stable_numerical_rank": maximum_numerical_rank,
            "frontier": all_frontier,
            "rank21_below_target_numerical_lead_count": len(
                numerical_target_leads
            ),
            "exact_rank_certificates": 0,
        },
        "parameters": {
            "screen_keep": args.screen_keep,
            "denominator_diversity": args.denominator_diversity,
            "height_keep": args.height_keep,
            "rank_keep": args.rank_keep,
            "batch_size": args.batch_size,
            "batch_timeout_seconds": args.batch_timeout,
            "height_timeout_seconds": args.height_timeout,
            "final_search_timeout_seconds": args.final_search_timeout,
            "escalation_search_timeout_seconds": args.escalation_search_timeout,
            "conductor_timeout_seconds": args.conductor_timeout,
            "pari_stack_bytes": args.stack_bytes,
            "output": str(args.output),
        },
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "reproducing_command": command,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}", flush=True)
    print(
        f"maximum stable numerical rank={maximum_numerical_rank}; "
        f"numerical rank-21/logN leads={len(numerical_target_leads)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
