#!/usr/bin/env python3
"""Deterministic Möbius-chart extension of Fermigier pair products.

The direct H=10^6 genus-3 search was stopped after its first two pairs each
exhausted a 20-second cap.  This fallback searches 24 exact projective charts
at chart height 50,000 for every one of the 220 published-direction pairs.

For k in +/-{2,4,8,12,16,19}, the charts are

    T = u + k,        T = 1/(u+k).

Every chart maps H(u)<=50,000 into projective height H(T)<=1,000,000, while
reaching both large-numerator and large-denominator regions beyond the prior
H(T)<=50,000 tranche.  The union is deterministic and materially deeper but
is not claimed to exhaust the full H(T)<=1,000,000 box.

Each transformed product is searched once with no retry.  Mapped parameters
must make both original quartic factors individually square and both forced
points are checked exactly.  Prior/generic fibers are excluded, conductor is
computed first, and only completed subtarget conductors receive H=50,000
rank triage and possible exact finite-reduction certification.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time
from typing import Any, Sequence

import sympy as sp

from ek_k3 import rational_to_string
from pari_bridge import pari_version
from search_fermigier_published_pair_fiber_products import (
    EXPECTED_LABELS,
    EXPECTED_PAIR_COUNT,
    EXPECTED_PRIOR_PARAMETER_COUNT,
    EXPECTED_PRIOR_PARAMETER_SHA256,
    EXPECTED_PUBLISHED_PREIMAGE_SHA256,
    PRIMARY_ARTIFACT,
    SPECIALIZATION_HEIGHT,
    T0,
    TARGET_LOG_CONDUCTOR,
    aggregate_candidates,
    classify_product_point,
    exact_slices,
    finalized_candidate_record,
    load_primary,
    pair_identifier,
    pair_population,
    polynomial_digest,
    prior_fermigier_parameters,
    published_accidentals,
    published_preimage_digest,
    rational_digest,
    sha256_bytes,
    sha256_file,
    triage_specialization,
    unique_product_parameters,
)
from search_fermigier_published_pair_fiber_products_h50000 import (
    EXPECTED_H50000_RESULT_SHA256,
    h50000_result_digest,
)
from search_fermigier_rank22_accidental_slices import (
    conductor_probe,
    poly_add,
    poly_multiply,
    search_polynomial,
)


CHART_HEIGHT = 50_000
GLOBAL_HEIGHT_CEILING = 1_000_000
SHIFTS = (-19, -16, -12, -8, -4, -2, 2, 4, 8, 12, 16, 19)
H50000_ARTIFACT = (
    "artifacts/generated-results/elliptic_fermigier_published_pair_fiber_products_h50000.json"
)
EXPECTED_TERMINAL_PRIOR_COUNT = 593
EXPECTED_TERMINAL_PRIOR_SHA256 = (
    "a4d06e4662d2e30c1a0f8873f91d8d348dae10f2abaffce88dcc0f480cfeede0"
)


def parse_precisions(value: str) -> tuple[int, ...]:
    values = tuple(int(part) for part in value.split(",") if part)
    if len(values) < 2 or tuple(sorted(set(values))) != values:
        raise argparse.ArgumentTypeError("provide increasing distinct precisions")
    return values


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chart-height", type=int, default=CHART_HEIGHT)
    parser.add_argument("--pair-timeout", type=float, default=15.0)
    parser.add_argument("--conductor-timeout", type=float, default=20.0)
    parser.add_argument("--specialization-timeout", type=float, default=30.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--saturation-timeout", type=float, default=60.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=2_000)
    parser.add_argument("--precisions", type=parse_precisions, default=(72, 120))
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_fermigier_published_pair_fiber_products_h1000000_charts.json",
    )
    return parser


def write_artifact(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")


def polynomial_power(
    polynomial: Sequence[Fraction], exponent: int
) -> tuple[Fraction, ...]:
    answer = (Fraction(1),)
    for _ in range(exponent):
        answer = poly_multiply(answer, polynomial)
    return answer


def mobius_quartic(
    coefficients: Sequence[Fraction], matrix: tuple[int, int, int, int]
) -> tuple[Fraction, ...]:
    """Return (c*u+d)^4 f((a*u+b)/(c*u+d)), ascending in u."""

    a_value, b_value, c_value, d_value = matrix
    numerator = (Fraction(b_value), Fraction(a_value))
    denominator = (Fraction(d_value), Fraction(c_value))
    terms = []
    for degree, coefficient in enumerate(coefficients):
        terms.append(
            tuple(
                coefficient * value
                for value in poly_multiply(
                    polynomial_power(numerator, degree),
                    polynomial_power(denominator, 4 - degree),
                )
            )
        )
    result = poly_add(*terms)
    if len(result) != 5:
        raise AssertionError("a declared Möbius chart lost quartic degree")
    return result


def chart_population() -> tuple[dict[str, Any], ...]:
    charts = []
    for kind in ("translate", "reciprocal-shift"):
        for shift in SHIFTS:
            matrix = (
                (1, shift, 0, 1)
                if kind == "translate"
                else (0, 1, 1, shift)
            )
            a_value, b_value, c_value, d_value = matrix
            if abs(a_value * d_value - b_value * c_value) != 1:
                raise AssertionError("a chart matrix was not unimodular")
            charts.append(
                {
                    "chart_id": f"{kind}-k{shift:+d}",
                    "kind": kind,
                    "shift": shift,
                    "matrix": list(matrix),
                    "mapped_projective_height_upper_bound": (
                        abs(shift) + 1
                    )
                    * CHART_HEIGHT,
                }
            )
    if len(charts) != 24 or max(
        chart["mapped_projective_height_upper_bound"] for chart in charts
    ) != GLOBAL_HEIGHT_CEILING:
        raise AssertionError("the deterministic chart population changed")
    return tuple(charts)


def map_chart_point(
    point: tuple[Fraction, Fraction], matrix: tuple[int, int, int, int]
) -> tuple[Fraction, Fraction] | None:
    u_value, chart_ordinate = map(Fraction, point)
    a_value, b_value, c_value, d_value = matrix
    denominator = c_value * u_value + d_value
    if denominator == 0:
        return None
    parameter = (a_value * u_value + b_value) / denominator
    product_ordinate = chart_ordinate / denominator**4
    return parameter, product_ordinate


def search_chart_pair(
    left: Any,
    right: Any,
    chart: dict[str, Any],
    *,
    timeout: float,
    stack_bytes: int,
    prior_parameters: set[Fraction],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    matrix = tuple(int(value) for value in chart["matrix"])
    left_chart = mobius_quartic(left.coefficients, matrix)
    right_chart = mobius_quartic(right.coefficients, matrix)
    product = poly_multiply(left_chart, right_chart)
    if len(product) != 9:
        raise AssertionError("a transformed pair product lost degree eight")
    points, search = search_polynomial(
        product,
        height_bound=CHART_HEIGHT,
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    incidences = []
    chart_poles = 0
    maximum_mapped_height = 0
    for chart_point in unique_product_parameters(points):
        mapped = map_chart_point(chart_point, matrix)
        if mapped is None:
            chart_poles += 1
            continue
        parameter, product_ordinate = mapped
        maximum_mapped_height = max(
            maximum_mapped_height,
            abs(parameter.numerator),
            parameter.denominator,
        )
        incidence = classify_product_point(
            left,
            right,
            (parameter, product_ordinate),
            prior_parameters,
        )
        incidence["chart_u"] = rational_to_string(chart_point[0])
        incidences.append(incidence)
    qualifying = [
        incidence
        for incidence in incidences
        if incidence["classification"] == "genuinely-new-double-forced-fiber"
    ]
    return (
        {
            "search": search,
            "transformed_product_polynomial_sha256": polynomial_digest(product),
            "chart_poles_excluded": chart_poles,
            "maximum_mapped_projective_height_observed": maximum_mapped_height,
            "incidences": incidences,
            "qualifying_new_parameter_count": len(
                {incidence["canonical_parameter_t"] for incidence in qualifying}
            ),
        },
        qualifying,
    )


def chart_result_digest(chart_rows: list[dict[str, Any]]) -> str:
    digest = hashlib.sha256()
    for chart in chart_rows:
        for row in chart["pair_searches"]:
            search = row["search"]
            digest.update(
                (
                    f"{chart['chart_id']}|{row['pair_id']}|"
                    f"{search['transformed_product_polynomial_sha256']}|"
                    f"{search['search']['status']}|"
                    f"{search['search'].get('signed_point_count')}|"
                    f"{search['qualifying_new_parameter_count']}\n"
                ).encode()
            )
    return digest.hexdigest()


def main() -> None:
    args = build_parser().parse_args()
    if args.chart_height != CHART_HEIGHT:
        raise SystemExit("the declared transformed-chart height is pinned at 50000")
    if min(
        args.pair_timeout,
        args.conductor_timeout,
        args.specialization_timeout,
        args.height_timeout,
        args.saturation_timeout,
    ) <= 0:
        raise SystemExit("all subprocess timeouts must be positive")
    if args.pair_timeout > 60:
        raise SystemExit("pair timeout may not exceed 60 seconds")

    root = Path(__file__).resolve().parents[2]
    artifact_dir = root / "artifacts" / "generated-results"
    primary_path = artifact_dir / PRIMARY_ARTIFACT
    primary, primary_raw = load_primary(primary_path)
    if published_preimage_digest(primary) != EXPECTED_PUBLISHED_PREIMAGE_SHA256:
        raise AssertionError("the exact published-preimage population changed")
    slices = exact_slices(published_accidentals(primary))
    pairs = pair_population(slices)
    charts = chart_population()

    h50000_path = root / H50000_ARTIFACT
    h50000 = json.loads(h50000_path.read_bytes())
    h50000_digest = h50000_result_digest(h50000["pair_searches"])
    if h50000_digest != EXPECTED_H50000_RESULT_SHA256:
        raise AssertionError("the exact H=50000 result digest changed")
    base_prior, base_prior_record = prior_fermigier_parameters(
        artifact_dir, args.output
    )
    if (
        len(base_prior) != EXPECTED_PRIOR_PARAMETER_COUNT
        or base_prior_record["prior_parameter_sha256"]
        != EXPECTED_PRIOR_PARAMETER_SHA256
    ):
        raise AssertionError("the exact base prior population changed")
    h50000_seen = {
        abs(Fraction(incidence["canonical_parameter_t"]))
        for row in h50000["pair_searches"]
        for incidence in row["search"]["incidences"]
    }
    prior_parameters = base_prior | h50000_seen
    prior_sha256 = rational_digest(sorted(prior_parameters))
    if (
        len(prior_parameters) != EXPECTED_TERMINAL_PRIOR_COUNT
        or prior_sha256 != EXPECTED_TERMINAL_PRIOR_SHA256
    ):
        raise AssertionError("the exact terminal prior population changed")

    script_path = Path(__file__).resolve()
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "in-progress deterministic transformed-chart pair screen",
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "hit": False,
        },
        "scope": {
            "chart_height": CHART_HEIGHT,
            "mapped_projective_height_ceiling": GLOBAL_HEIGHT_CEILING,
            "full_H1000000_exhaustive": False,
            "claim": (
                "deterministic 24-chart union reaches numerator- and denominator-heavy "
                "regions through height one million; it is not the full projective box"
            ),
        },
        "source": {
            "published_preimage_sha256": published_preimage_digest(primary),
            "published_accidental_labels": list(EXPECTED_LABELS),
            "primary_artifact_sha256_observed": sha256_bytes(primary_raw),
            "H50000_artifact_sha256_observed": sha256_file(h50000_path),
            "H50000_exact_pair_result_sha256": h50000_digest,
        },
        "prior_decontamination": {
            "base_prior_parameter_count": len(base_prior),
            "base_prior_parameter_sha256": base_prior_record[
                "prior_parameter_sha256"
            ],
            "H50000_seen_parameters": [
                rational_to_string(value) for value in sorted(h50000_seen)
            ],
            "terminal_prior_parameter_count": len(prior_parameters),
            "terminal_prior_parameter_sha256": prior_sha256,
        },
        "parameters": {
            "shifts": list(SHIFTS),
            "chart_count": len(charts),
            "pair_count_per_chart": EXPECTED_PAIR_COUNT,
            "declared_search_call_count": len(charts) * EXPECTED_PAIR_COUNT,
            "pair_timeout_seconds": args.pair_timeout,
            "conductor_timeout_seconds": args.conductor_timeout,
            "specialization_height": SPECIALIZATION_HEIGHT,
            "specialization_timeout_seconds": args.specialization_timeout,
            "height_timeout_seconds": args.height_timeout,
            "height_precisions": list(args.precisions),
            "stack_bytes": args.stack_bytes,
            "no_retries": True,
            "checkpoint_after_each_chart": True,
        },
        "chart_population": list(charts),
        "chart_searches": [],
        "candidates": [],
        "execution": {"phase": "chart-search-in-progress", "charts_completed": 0},
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
            "pari_gp": pari_version(),
        },
        "reproducing_command": " ".join(
            shlex.quote(part) for part in [sys.executable, *sys.argv]
        ),
        "script_sha256": sha256_file(script_path),
    }
    write_artifact(args.output, artifact)

    qualifying_all: list[tuple[str, dict[str, Any]]] = []
    started = time.monotonic()
    for chart in charts:
        chart_started = time.monotonic()
        chart_row = {
            **chart,
            "pair_searches": [],
        }
        for left, right in pairs:
            pair_id = pair_identifier(left, right)
            search, qualifying = search_chart_pair(
                left,
                right,
                chart,
                timeout=args.pair_timeout,
                stack_bytes=args.stack_bytes,
                prior_parameters=prior_parameters,
            )
            chart_row["pair_searches"].append(
                {
                    "pair_id": pair_id,
                    "left_source_label": left.accidental_label,
                    "right_source_label": right.accidental_label,
                    "search": search,
                }
            )
            qualifying_all.extend(
                (f"{chart['chart_id']}:{pair_id}", incidence)
                for incidence in qualifying
            )
        chart_row["wall_seconds"] = time.monotonic() - chart_started
        artifact["chart_searches"].append(chart_row)
        artifact["execution"].update(
            {
                "charts_completed": len(artifact["chart_searches"]),
                "last_chart_id": chart["chart_id"],
                "wall_seconds_so_far": time.monotonic() - started,
            }
        )
        write_artifact(args.output, artifact)

    aggregated = aggregate_candidates(qualifying_all)
    candidates = [
        finalized_candidate_record(candidate)
        for _, candidate in sorted(aggregated.items())
    ]
    candidates = [
        candidate
        for candidate in candidates
        if candidate["distinct_published_source_direction_count"] >= 2
        and candidate["distinct_group_pullback_classes_modulo_inverse"] >= 2
    ]
    artifact["candidates"] = candidates
    artifact["execution"]["phase"] = "conductor-first"
    write_artifact(args.output, artifact)
    for candidate in candidates:
        parameter = Fraction(candidate["parameter_t"])
        candidate["conductor_probe"] = conductor_probe(
            parameter,
            timeout=args.conductor_timeout,
            stack_bytes=args.stack_bytes,
        )
        write_artifact(args.output, artifact)
        if candidate["conductor_probe"].get("below_strict_log_conductor_target"):
            try:
                candidate["rank_triage"] = triage_specialization(
                    candidate,
                    search_timeout=args.specialization_timeout,
                    height_timeout=args.height_timeout,
                    precisions=args.precisions,
                    stack_bytes=args.stack_bytes,
                    saturation_timeout=args.saturation_timeout,
                    certificate_prime_bound=args.certificate_prime_bound,
                )
            except subprocess.TimeoutExpired as error:
                candidate["rank_triage"] = {
                    "status": "timeout-no-retry",
                    "error": str(error)[:1000],
                }
            except (RuntimeError, AssertionError, ValueError) as error:
                candidate["rank_triage"] = {
                    "status": "error-no-retry",
                    "error": str(error)[:1000],
                }
            write_artifact(args.output, artifact)

    all_pair_rows = [
        row
        for chart in artifact["chart_searches"]
        for row in chart["pair_searches"]
    ]
    classifications = Counter(
        incidence["classification"]
        for row in all_pair_rows
        for incidence in row["search"]["incidences"]
    )
    rank_records = [
        candidate["rank_triage"]
        for candidate in candidates
        if "rank_triage" in candidate
        and "full_pool_stable_numerical_rank" in candidate["rank_triage"]
    ]
    artifact["outcome"] = {
        "declared_search_call_count": len(charts) * EXPECTED_PAIR_COUNT,
        "search_calls_attempted": len(all_pair_rows),
        "search_calls_completed": sum(
            row["search"]["search"]["status"] == "completed"
            for row in all_pair_rows
        ),
        "search_calls_timed_out_or_errored": sum(
            row["search"]["search"]["status"] != "completed"
            for row in all_pair_rows
        ),
        "wall_seconds": time.monotonic() - started,
        "incidence_classification_counts": dict(sorted(classifications.items())),
        "maximum_mapped_projective_height_observed": max(
            (
                row["search"]["maximum_mapped_projective_height_observed"]
                for row in all_pair_rows
            ),
            default=0,
        ),
        "genuinely_new_double_forced_fibers": len(candidates),
        "completed_conductors": sum(
            candidate.get("conductor_probe", {}).get("status") == "completed"
            for candidate in candidates
        ),
        "subtarget_conductors": sum(
            candidate.get("conductor_probe", {}).get(
                "below_strict_log_conductor_target"
            )
            is True
            for candidate in candidates
        ),
        "rank_triage_count": len(rank_records),
        "maximum_stable_numerical_rank": max(
            (
                record["full_pool_stable_numerical_rank"]
                for record in rank_records
            ),
            default=None,
        ),
        "exact_chart_result_sha256": chart_result_digest(
            artifact["chart_searches"]
        ),
    }
    artifact["target"]["hit"] = any(
        candidate.get("rank_triage", {})
        .get("finite_reduction_attempt", {})
        .get("certified_algebraic_rank_lower_bound", 0)
        >= 21
        and candidate.get("conductor_probe", {}).get(
            "below_strict_log_conductor_target"
        )
        for candidate in candidates
    )
    if not artifact["target"]["hit"]:
        artifact["target"]["reason"] = (
            "no new deterministic-chart double-forced subtarget fiber received an exact rank-21 certificate"
        )
    artifact["status"] = (
        "completed deterministic 24-chart pair screen through mapped projective "
        "height one million; not a full H=1000000 exhaustion"
    )
    artifact["execution"]["phase"] = "complete"
    artifact["execution"]["wall_seconds"] = time.monotonic() - started
    artifact["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
    write_artifact(args.output, artifact)


if __name__ == "__main__":
    main()
