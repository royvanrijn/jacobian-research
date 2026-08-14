#!/usr/bin/env python3
"""Search alternate degree-two covers based at weight-two/three subgroup sums.

For a known point ``Q=(X_Q,Y_Q)`` on the integral short model, put

``t_Q(P) = (Y(P)+Y_Q)/(X(P)-X_Q)``.

This is the slope of the line through ``P`` and ``-Q``.  Dividing its
intersection cubic by the known root ``X_Q`` gives a quartic ``D_Q(t)``;
rational points on ``w^2=D_Q(t)`` map to pairs ``P_1+P_2=Q``.  Changing Q
changes the degree-two coordinate, so these are not merely affine charts of
the 29 public-base covers searched previously.

The finite tranche is pinned: form all 4060 positive subset sums of weights
two and three, retain the 64 with smallest public ``t_Q`` heights, and search
three offset plus four best cross-ratio charts per retained cover.
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
import subprocess
import sys
import time
from typing import Any, Iterable

from elkies_klagsbrun_rank29 import (
    PUBLISHED_POINTS,
    from_short_point,
    point_on_general_curve,
    point_on_short_curve,
    published_short_points,
    short_weierstrass_coefficients,
    to_short_point,
)
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)
from pari_bridge import pari_version
from search_elkies_klagsbrun_rank30 import (
    Polynomial,
    QuarticChart,
    RationalPoint,
    affine_substitute,
    discover_relation,
    point_add,
    point_negate,
    poly_add,
    poly_evaluate,
    poly_multiply,
    poly_scale,
    search_chart,
)


Q = Fraction
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/elliptic_elkies_klagsbrun_rank30_alternate_covers.json"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_elkies_klagsbrun_rank30_alternate_covers.py"
)


@dataclass(frozen=True)
class AlternateCover:
    subset_indices: tuple[int, ...]
    general_point: RationalPoint
    short_point: RationalPoint
    public_parameters: tuple[Fraction, ...]
    raw_height_score: tuple[int, ...]

    @property
    def identifier(self) -> str:
        indices = "_".join(f"{index + 1:02d}" for index in self.subset_indices)
        return f"w{len(self.subset_indices)}_{indices}"


@dataclass(frozen=True)
class AlternateChart:
    chart: QuarticChart
    cover: AlternateCover


def parameter_height(value: Fraction) -> int:
    return max(abs(value.numerator), value.denominator)


def short_alternate_discriminant(short_q: RationalPoint) -> Polynomial:
    """Return ``D_Q(t)`` for lines through ``-Q`` on the short model."""

    x_q, y_q = short_q
    coefficient_a = short_weierstrass_coefficients()[3]
    t: Polynomial = (Q(0), Q(1))
    intercept: Polynomial = (-y_q, -x_q)
    c2 = poly_multiply(t, t)
    alpha = poly_add((x_q,), poly_scale(c2, Q(-1)))
    c1 = poly_add(poly_scale(poly_multiply(t, intercept), Q(2)), (-coefficient_a,))
    discriminant = poly_add(
        poly_add(
            poly_multiply(alpha, alpha),
            poly_scale(poly_multiply(alpha, (x_q,)), Q(-4)),
        ),
        poly_scale(c1, Q(4)),
    )
    if len(discriminant) != 5 or discriminant[-1] != 1:
        raise AssertionError("an alternate discriminant is not a monic quartic")
    return discriminant


def alternate_parameter(short_q: RationalPoint, short_point: RationalPoint) -> Fraction:
    x_q, y_q = short_q
    x_value, y_value = short_point
    if x_value == x_q:
        raise ValueError("the alternate slope is vertical")
    return (y_value + y_q) / (x_value - x_q)


def alternate_to_short_points(
    short_q: RationalPoint, parameter: Fraction, square_root: Fraction
) -> tuple[RationalPoint, RationalPoint]:
    x_q, y_q = short_q
    parameter = Q(parameter)
    square_root = Q(square_root)
    alpha = x_q - parameter * parameter
    first_x = (-alpha + square_root) / 2
    second_x = (-alpha - square_root) / 2
    intercept = -y_q - parameter * x_q
    points = (
        (first_x, parameter * first_x + intercept),
        (second_x, parameter * second_x + intercept),
    )
    if not all(point_on_short_curve(point) for point in points):
        raise AssertionError("an alternate-cover image missed the short curve")
    return points


def build_all_covers() -> tuple[AlternateCover, ...]:
    short_public = published_short_points()
    covers: list[AlternateCover] = []
    for weight in (2, 3):
        for indices in combinations(range(len(PUBLISHED_POINTS)), weight):
            general_sum = None
            for index in indices:
                general_sum = point_add(general_sum, PUBLISHED_POINTS[index])
            if general_sum is None:
                raise AssertionError("a positive subset sum unexpectedly vanished")
            short_sum = to_short_point(general_sum)
            parameters = tuple(
                alternate_parameter(short_sum, point) for point in short_public
            )
            heights = tuple(sorted(parameter_height(value) for value in parameters))
            covers.append(
                AlternateCover(
                    subset_indices=indices,
                    general_point=general_sum,
                    short_point=short_sum,
                    public_parameters=parameters,
                    raw_height_score=heights[:5],
                )
            )
    covers.sort(
        key=lambda cover: (
            cover.raw_height_score,
            len(cover.subset_indices),
            cover.subset_indices,
        )
    )
    if len(covers) != 4060:
        raise AssertionError("the weight-two/three cover count changed")
    return tuple(covers)


def build_alternate_charts(
    *,
    cover_count: int,
    offset_count: int,
    cross_ratio_count: int,
    offset_height: int,
    cross_ratio_height: int,
) -> tuple[AlternateChart, ...]:
    selected = build_all_covers()[:cover_count]
    answers: list[AlternateChart] = []
    thresholds = (10**4, 10**6, 10**9, 10**12, 10**15, 10**18)
    for cover in selected:
        discriminant = short_alternate_discriminant(cover.short_point)
        distinct_parameters = tuple(
            sorted(
                set(cover.public_parameters),
                key=lambda value: (parameter_height(value), value),
            )
        )
        for offset_index, center in enumerate(distinct_parameters[:offset_count], 1):
            answers.append(
                AlternateChart(
                    chart=QuarticChart(
                        identifier=f"alt_{cover.identifier}_offset{offset_index}",
                        kind="alternate_cover_offset",
                        polynomial=affine_substitute(discriminant, center, Q(1)),
                        center=center,
                        scale=Q(1),
                        base_index=None,
                        seed_parameters=(Q(0),),
                        height_bound=offset_height,
                    ),
                    cover=cover,
                )
            )

        ranked: list[tuple[Any, QuarticChart]] = []
        # Cross-ratio scoring with every pair would spend most of the fixed
        # budget normalizing against enormous public parameters.  The eight
        # smallest raw t-heights are the declared normalization pool.
        normalization_pool_size = min(8, len(distinct_parameters))
        for first_index, second_index in combinations(range(normalization_pool_size), 2):
            center = distinct_parameters[first_index]
            scale = distinct_parameters[second_index] - center
            normalized = tuple(
                (value - center) / scale for value in cover.public_parameters
            )
            counts = tuple(
                sum(parameter_height(value) <= threshold for value in normalized)
                for threshold in thresholds
            )
            polynomial = affine_substitute(discriminant, center, scale)
            coefficient_bits = max(
                max(abs(value.numerator).bit_length(), value.denominator.bit_length())
                for value in polynomial
            )
            chart = QuarticChart(
                identifier=(
                    f"alt_{cover.identifier}_cross{first_index + 1:02d}_"
                    f"{second_index + 1:02d}"
                ),
                kind="alternate_cover_cross_ratio",
                polynomial=polynomial,
                center=center,
                scale=scale,
                base_index=None,
                seed_parameters=(Q(0), Q(1)),
                height_bound=cross_ratio_height,
            )
            score = tuple(-count for count in reversed(counts)) + (
                coefficient_bits,
                chart.identifier,
            )
            ranked.append((score, chart))
        ranked.sort(key=lambda item: item[0])
        answers.extend(
            AlternateChart(chart=chart, cover=cover)
            for _, chart in ranked[:cross_ratio_count]
        )
    return tuple(answers)


def alternate_chart_images(
    alternate_chart: AlternateChart,
    quartic_points: Iterable[tuple[Fraction, Fraction]],
) -> tuple[RationalPoint, ...]:
    chart = alternate_chart.chart
    images: list[RationalPoint] = []
    for local_parameter, square_root in quartic_points:
        if poly_evaluate(chart.polynomial, local_parameter) != square_root * square_root:
            raise AssertionError("PARI returned a point off an alternate quartic")
        if local_parameter in chart.seed_parameters:
            continue
        parameter = chart.center + chart.scale * local_parameter
        for short_point in alternate_to_short_points(
            alternate_chart.cover.short_point, parameter, square_root
        ):
            point = from_short_point(short_point)
            if not point_on_general_curve(point):
                raise AssertionError("inverse transport returned an off-curve point")
            images.append(point)
    return tuple(images)


def point_record(point: RationalPoint) -> dict[str, str]:
    return {"x": str(point[0]), "y": str(point[1])}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cover-count", type=int, default=64)
    parser.add_argument("--offset-count", type=int, default=3)
    parser.add_argument("--cross-ratio-count", type=int, default=4)
    parser.add_argument("--offset-height", type=int, default=100_000)
    parser.add_argument("--cross-ratio-height", type=int, default=50_000)
    parser.add_argument("--chart-timeout", type=float, default=8.0)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=500_000_000)
    parser.add_argument("--certificate-prime-bound", type=int, default=1000)
    args = parser.parse_args()
    if not 1 <= args.cover_count <= 4060:
        raise SystemExit("--cover-count must lie in [1,4060]")
    if min(
        args.offset_count,
        args.cross_ratio_count,
        args.offset_height,
        args.cross_ratio_height,
    ) <= 0:
        raise SystemExit("chart counts and heights must be positive")
    if args.chart_timeout <= 0 or args.relation_timeout <= 0:
        raise SystemExit("timeouts must be positive")

    charts = build_alternate_charts(
        cover_count=args.cover_count,
        offset_count=args.offset_count,
        cross_ratio_count=args.cross_ratio_count,
        offset_height=args.offset_height,
        cross_ratio_height=args.cross_ratio_height,
    )
    started = time.monotonic()
    public_set = set(PUBLISHED_POINTS) | {
        point_negate(point) for point in PUBLISHED_POINTS
    }
    discovered: dict[RationalPoint, list[str]] = {}
    completed: list[dict[str, Any]] = []
    timeouts: list[str] = []
    total_pari_milliseconds = 0
    for index, alternate_chart in enumerate(charts, 1):
        chart = alternate_chart.chart
        try:
            quartic_points, milliseconds, wall_seconds = search_chart(
                chart, timeout=args.chart_timeout, stack_bytes=args.stack_bytes
            )
        except subprocess.TimeoutExpired:
            timeouts.append(chart.identifier)
            print(f"timeout {chart.identifier}", flush=True)
            continue
        images = alternate_chart_images(alternate_chart, quartic_points)
        unexpected = tuple(point for point in images if point not in public_set)
        for point in unexpected:
            discovered.setdefault(point, []).append(chart.identifier)
        total_pari_milliseconds += milliseconds
        completed.append(
            {
                "identifier": chart.identifier,
                "cover": alternate_chart.cover.identifier,
                "kind": chart.kind,
                "height_bound": chart.height_bound,
                "quartic_affine_point_count_including_signs": len(quartic_points),
                "nonseed_mapped_image_count_including_duplicates": len(images),
                "pari_milliseconds": milliseconds,
                "wall_seconds": wall_seconds,
            }
        )
        if index % 50 == 0 or unexpected:
            print(
                f"alternate charts {index}/{len(charts)}; "
                f"timeouts={len(timeouts)}; unexpected_images={len(discovered)}",
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
                    "classification": "exactly_in_published_rank29_subgroup",
                    "published_basis_relation": list(relation),
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
            binary_rank = combined_mod2_rank(signatures, 30)
            record.update(
                {
                    "classification": (
                        "exact_independent_30th_point"
                        if binary_rank == 30
                        else "unresolved_after_relation_and_mod2_search"
                    ),
                    "augmented_mod2_rank": binary_rank,
                    "certificate_primes": [signature.prime for signature in signatures],
                    "certificate_prime_bound": args.certificate_prime_bound,
                }
            )
            target_hit |= binary_rank == 30
        candidate_records.append(record)

    selected_covers: list[AlternateCover] = []
    seen_cover_ids: set[str] = set()
    for alternate_chart in charts:
        if alternate_chart.cover.identifier not in seen_cover_ids:
            selected_covers.append(alternate_chart.cover)
            seen_cover_ids.add(alternate_chart.cover.identifier)
    manifest_payload = [
        {
            "identifier": alternate_chart.chart.identifier,
            "cover": alternate_chart.cover.identifier,
            "center": str(alternate_chart.chart.center),
            "scale": str(alternate_chart.chart.scale),
            "height": alternate_chart.chart.height_bound,
            "coefficients": [str(value) for value in alternate_chart.chart.polynomial],
        }
        for alternate_chart in charts
    ]
    manifest_sha256 = hashlib.sha256(
        json.dumps(manifest_payload, separators=(",", ":")).encode()
    ).hexdigest()
    artifact = {
        "schema_version": 1,
        "artifact_kind": "bounded_alternate_2cover_rank30_point_search",
        "status": (
            "exact_rank30_target_hit"
            if target_hit
            else "bounded_search_no_certified_30th_point"
        ),
        "claim_scope": {
            "exact": (
                "all 4060 subgroup sums, cover equations, chart maps, curve "
                "membership, relation replays, and finite-reduction certificates"
            ),
            "bounded": (
                "PARI enumeration only in the completed retained-cover chart boxes; "
                "no rank upper bound"
            ),
        },
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "selection": {
            "complete_weight_two_three_cover_count": 4060,
            "retained_cover_count": len(selected_covers),
            "ranking": (
                "lexicographic first five sorted naive heights of the 29 public "
                "t_Q coordinates, then subset weight and indices"
            ),
            "cross_ratio_normalization_pool": (
                "the eight distinct public t_Q coordinates of smallest naive height"
            ),
            "retained_covers": [
                {
                    "identifier": cover.identifier,
                    "subset_indices_one_based": [
                        index + 1 for index in cover.subset_indices
                    ],
                    "Q": point_record(cover.general_point),
                    "short_Q": point_record(cover.short_point),
                    "raw_height_score": list(cover.raw_height_score),
                }
                for cover in selected_covers
            ],
        },
        "parameters": {
            "cover_count": args.cover_count,
            "offset_charts_per_cover": args.offset_count,
            "cross_ratio_charts_per_cover": args.cross_ratio_count,
            "offset_height": args.offset_height,
            "cross_ratio_height": args.cross_ratio_height,
            "chart_timeout_seconds_each": args.chart_timeout,
            "relation_timeout_seconds_each": args.relation_timeout,
            "stack_bytes_each": args.stack_bytes,
            "certificate_prime_bound": args.certificate_prime_bound,
        },
        "chart_manifest": {
            "count": len(charts),
            "sha256": manifest_sha256,
            "kind_counts": {
                kind: sum(alternate_chart.chart.kind == kind for alternate_chart in charts)
                for kind in sorted(
                    {alternate_chart.chart.kind for alternate_chart in charts}
                )
            },
            "manifest_stored_inline": False,
            "regeneration": "deterministic build_alternate_charts in pinned script",
        },
        "search_result": {
            "completed_chart_count": len(completed),
            "timed_out_chart_count": len(timeouts),
            "timed_out_chart_identifiers": timeouts,
            "total_pari_milliseconds_completed": total_pari_milliseconds,
            "wall_seconds": time.monotonic() - started,
            "unique_nonpublic_images": len(discovered),
            "candidate_records": candidate_records,
            "certified_independent_30th_point_count": sum(
                record["classification"] == "exact_independent_30th_point"
                for record in candidate_records
            ),
            "rank30_target_hit": target_hit,
            "completed_chart_summaries": completed,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")
    print(
        f"completed={len(completed)}/{len(charts)} "
        f"unique_nonpublic={len(discovered)} target_hit={str(target_hit).lower()}"
    )


if __name__ == "__main__":
    main()
