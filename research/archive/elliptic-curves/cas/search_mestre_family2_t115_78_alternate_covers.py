#!/usr/bin/env python3
"""Bounded alternate-cover search on the low-conductor Mestre T=115/78 fiber.

The split-infinity specialization ``u=2553/13`` of the six-root family
``(0,25,95,143,168,205)`` has an exact mod-3 rank lower bound 17 and
``log(N)=140.0425...``.  This experiment reconstructs its certified pivot
basis from the cached three-chart point search and searches genuinely
different degree-two models attached to short nonzero subset sums.

Every returned point and every claimed rank gain is checked exactly.  The
bounded chart search is not a rank upper bound.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import itertools
import json
import random
from pathlib import Path
import sys
from typing import Any

from alternate_quartic_covers import alternate_cover, short_subset_sum
from search_mestre_dsquare_four import (
    FAMILIES,
    RELATION_PRIME_BOUND,
    atomic_json,
    base_parameter,
    candidate_identifier,
    known_jacobian_points,
    mod3_independence_certificate,
    point_digest,
    quartic_point_to_jacobian,
    quartic_value,
    rational_square_root,
    rational_text,
    run_ratpoints_chart,
)
from search_nagao_u135_alternate_covers import (
    CoverPlan,
    best_cross_ratio_charts,
    cover_parameters,
    point_record,
    projective_height,
    run_chart,
)


Q = Fraction
FAMILY_INDEX = 2
PARAMETER_U = Q(2553, 13)
PARAMETER_T = Q(115, 78)
EXPECTED_PIVOTS = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 19, 21)
EXPECTED_BASIS_SHA256 = "c6589bb0ff2a35cc802f64db3b19c4341fcf6639c2a37921e045b8bda904f59f"
STACK_BYTES = 512_000_000


def reconstruct_pool(
    raw_root: Path, parameter_u: Fraction, family_index: int = FAMILY_INDEX
) -> tuple[
    tuple[Fraction, ...],
    tuple[tuple[Fraction, Fraction], ...],
    dict[str, Any],
]:
    family = FAMILIES[family_index]
    parameter_u = Q(parameter_u)
    parameter_t = base_parameter(family, parameter_u)
    construction = family.construction
    coefficients = construction.primitive_jacobian_coefficients(parameter_t)
    raw_directory = raw_root / candidate_identifier(
        {
            "family_index": family_index,
            "numerator": parameter_u.numerator,
            "denominator": parameter_u.denominator,
        }
    )
    charts = tuple(
        run_ratpoints_chart(family, parameter_t, chart, raw_directory, 5.0)
        for chart in ("raw", "plus-T", "minus-T")
    )
    searched_x = sorted(
        {x_value for chart in charts for x_value in chart["finite_abscissae"]}
    )
    quartic = construction.primitive_quartic_coefficients(parameter_t)
    by_x: dict[Fraction, tuple[Fraction, Fraction]] = {}
    for point in known_jacobian_points(family, parameter_u):
        by_x.setdefault(point[0], point)
    for x_value in searched_x:
        square_root = rational_square_root(quartic_value(quartic, x_value))
        if square_root is None:
            raise AssertionError("a cached abscissa left the fixed quartic")
        if square_root == 0:
            continue
        point = quartic_point_to_jacobian(
            construction, parameter_t, (x_value, square_root)
        )
        by_x.setdefault(point[0], point)
    pool = tuple(by_x.values())
    certificate = mod3_independence_certificate(
        coefficients, pool, prime_bound=RELATION_PRIME_BOUND
    )
    pivots = tuple(certificate["independent_subset_indices_one_based"])
    basis = tuple(pool[index - 1] for index in pivots)
    if family_index == FAMILY_INDEX and parameter_u == PARAMETER_U:
        if parameter_t != PARAMETER_T or pivots != EXPECTED_PIVOTS:
            raise AssertionError("the pinned T=115/78 pivot set changed")
        if point_digest(basis) != EXPECTED_BASIS_SHA256:
            raise AssertionError("the pinned T=115/78 pivot basis changed")
    exact_rank = certificate["combined_exact_rank_over_F3"]
    if parameter_u == PARAMETER_U and exact_rank != 17:
        raise AssertionError("the pinned specialization lost rank 17")
    chart_records = [
        {key: value for key, value in chart.items() if key != "finite_abscissae"}
        for chart in charts
    ]
    return coefficients, basis, {
        "pool_point_count_modulo_inverse": len(pool),
        "pool_sha256": point_digest(pool),
        "basis_sha256": point_digest(basis),
        "parameter_T": rational_text(parameter_t),
        "pivot_indices_one_based": list(pivots),
        "certificate_primes": certificate["certificate_primes"],
        "certified_rank_lower_bound": exact_rank,
        "charts": chart_records,
    }


def build_plans(
    coefficients: tuple[Fraction, ...],
    basis: tuple[tuple[Fraction, Fraction], ...],
    *,
    minimum_weight: int,
    maximum_weight: int,
    charts_per_cover: int,
    cover_count: int,
) -> tuple[CoverPlan, ...]:
    plans: list[CoverPlan] = []
    for weight in range(minimum_weight, maximum_weight + 1):
        for indices in itertools.combinations(range(len(basis)), weight):
            base_point = short_subset_sum(coefficients, basis, indices)
            if base_point is None:
                raise AssertionError("a nonempty independent subset vanished")
            cover = alternate_cover(coefficients, base_point)
            charts = best_cross_ratio_charts(
                cover, basis, count=charts_per_cover
            )
            plans.append(CoverPlan(indices, cover, charts))
    plans.sort(key=lambda plan: plan.score)
    return tuple(plans[:cover_count])


def build_sampled_plans(
    coefficients: tuple[Fraction, ...],
    basis: tuple[tuple[Fraction, Fraction], ...],
    *,
    sample_count: int,
    charts_per_cover: int,
    cover_count: int,
) -> tuple[CoverPlan, ...]:
    """Rank a deterministic sample of higher-weight mod-2 classes."""

    generator = random.Random(1150782553)
    masks: set[int] = set()
    while len(masks) < sample_count:
        mask = generator.randrange(1, 1 << len(basis))
        if mask.bit_count() >= 3:
            masks.add(mask)
    retained: list[tuple[tuple[int, int, int], tuple[int, ...], Any]] = []
    for mask in sorted(masks):
        indices = tuple(index for index in range(len(basis)) if mask >> index & 1)
        base_point = short_subset_sum(coefficients, basis, indices)
        if base_point is None:
            raise AssertionError("a sampled independent subset vanished")
        cover = alternate_cover(coefficients, base_point)
        heights = tuple(
            projective_height(parameter).bit_length()
            for _, parameter in cover_parameters(cover, basis)
        )
        score = (max(heights), sum(heights), mask)
        retained.append((score, indices, cover))
    retained.sort(key=lambda item: item[0])
    plans = []
    for _, indices, cover in retained[:cover_count]:
        plans.append(
            CoverPlan(
                indices,
                cover,
                best_cross_ratio_charts(cover, basis, count=charts_per_cover),
            )
        )
    plans.sort(key=lambda plan: plan.score)
    return tuple(plans)


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-root",
        type=Path,
        default=root / "artifacts/local/elliptic-curves/mestre-dsquare-family2-expanded/ratpoints-raw",
    )
    parser.add_argument("--maximum-weight", type=int, default=1)
    parser.add_argument("--minimum-weight", type=int, default=1)
    parser.add_argument("--cover-count", type=int, default=8)
    parser.add_argument("--charts-per-cover", type=int, default=2)
    parser.add_argument("--height", type=int, default=50_000)
    parser.add_argument("--timeout", type=float, default=8.0)
    parser.add_argument("--numerator", type=int, default=2553)
    parser.add_argument("--denominator", type=int, default=13)
    parser.add_argument("--family-index", type=int, default=FAMILY_INDEX)
    parser.add_argument(
        "--sampled-mask-count",
        type=int,
        default=0,
        help="rank this many deterministic higher-weight mod-2 classes",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "artifacts/local/elliptic-curves/mestre-family2-t115-78-altcovers-pilot.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 1 <= args.minimum_weight <= args.maximum_weight <= 2:
        raise SystemExit("subset weights must satisfy 1 <= minimum <= maximum <= 2")
    if min(args.cover_count, args.charts_per_cover, args.height) <= 0:
        raise SystemExit("cover, chart, and height bounds must be positive")
    if not 0 <= args.family_index < len(FAMILIES):
        raise SystemExit("--family-index is outside the declared family range")
    parameter_u = Q(args.numerator, args.denominator)
    coefficients, basis, baseline = reconstruct_pool(
        args.raw_root, parameter_u, args.family_index
    )
    parameter_t = base_parameter(FAMILIES[args.family_index], parameter_u)
    baseline_rank = len(basis)
    if baseline_rank < 12:
        raise SystemExit("the reconstructed specialization has an implausibly small basis")
    if args.sampled_mask_count:
        plans = build_sampled_plans(
            coefficients,
            basis,
            sample_count=args.sampled_mask_count,
            charts_per_cover=args.charts_per_cover,
            cover_count=args.cover_count,
        )
    else:
        plans = build_plans(
            coefficients,
            basis,
            minimum_weight=args.minimum_weight,
            maximum_weight=args.maximum_weight,
            charts_per_cover=args.charts_per_cover,
            cover_count=args.cover_count,
        )
    basis_label = f"rank{baseline_rank}-basis"
    by_x: dict[Fraction, tuple[tuple[Fraction, Fraction], set[str]]] = {
        point[0]: (point, {basis_label}) for point in basis
    }
    runs: list[dict[str, Any]] = []
    for plan_index, plan in enumerate(plans, start=1):
        for chart_index, chart in enumerate(plan.charts, start=1):
            record, points = run_chart(
                plan,
                chart,
                stage=f"pilot-{plan_index}-{chart_index}",
                height_bound=args.height,
                timeout=args.timeout,
                stack_bytes=STACK_BYTES,
            )
            source = f"{plan.identifier}:chart-{chart_index}"
            novel = 0
            for point in points:
                if point[0] not in by_x:
                    by_x[point[0]] = (point, {source})
                    novel += 1
                else:
                    by_x[point[0]][1].add(source)
            record.update(
                {
                    "cover_subset_indices_one_based": [i + 1 for i in plan.subset_indices],
                    "chart_basis_indices_one_based": [i + 1 for i in chart.basis_indices],
                    "matrix": list(chart.matrix),
                    "novel_abscissae": novel,
                }
            )
            runs.append(record)
        print(
            f"cover={plan.identifier} completed={plan_index}/{len(plans)} "
            f"distinct_x={len(by_x)}",
            flush=True,
        )

    candidates = tuple(
        sorted(
            (entry[0] for x_value, entry in by_x.items() if basis_label not in entry[1]),
            key=lambda point: (projective_height(point[0]), point[0]),
        )
    )
    augmented = mod3_independence_certificate(
        coefficients, basis + candidates, prime_bound=RELATION_PRIME_BOUND
    )
    rank = augmented["combined_exact_rank_over_F3"]
    artifact = {
        "schema_version": 1,
        "status": "bounded alternate-cover pilot complete",
        "curve": {
            "family_index": args.family_index,
            "roots": list(FAMILIES[args.family_index].roots),
            "u": rational_text(parameter_u),
            "T": rational_text(parameter_t),
            "short_weierstrass_coefficients": [rational_text(value) for value in coefficients],
            "global_data": {
                "2553/13": {
                    "conductor": "6602174625938019235476464723256005161095499069821125031035190",
                    "log_conductor": "140.042504663184101558688341988343627345860523043619348958326",
                    "root_number": -1,
                },
                "202": {
                    "conductor": "449710538694747784031757572041467551513020085692209200945231300468333034",
                    "log_conductor": "164.986975545036178257262931955868682784724377350740842036578",
                    "root_number": -1,
                },
            }.get(rational_text(parameter_u)),
        },
        "baseline": baseline,
        "budget": {
            "maximum_subset_weight": args.maximum_weight,
            "minimum_subset_weight": args.minimum_weight,
            "selected_cover_count": len(plans),
            "charts_per_cover": args.charts_per_cover,
            "height_bound": args.height,
            "timeout_seconds_per_chart": args.timeout,
            "deterministic_higher_weight_mask_sample": args.sampled_mask_count,
        },
        "plans": [
            {
                "subset_indices_one_based": [i + 1 for i in plan.subset_indices],
                "score": list(plan.score[:3]),
            }
            for plan in plans
        ],
        "runs": runs,
        "results": {
            "nonbasis_candidate_count_modulo_inverse": len(candidates),
            "candidate_sha256": point_digest(candidates),
            "combined_exact_rank_over_F3": rank,
            "certified_new_directions": rank - baseline_rank,
            "rank21_target_achieved": rank >= 21,
            "independent_subset_indices_one_based": augmented[
                "independent_subset_indices_one_based"
            ],
            "certificate_primes": augmented["certificate_primes"],
            "candidate_points": [
                {
                    **point_record(point),
                    "sources": sorted(by_x[point[0]][1]),
                }
                for point in candidates
            ],
        },
        "interpretation": {
            "exact": "point membership and the mod-3 rank lower bound are exact",
            "bounded": "the declared cover charts only; no rank upper bound",
        },
    }
    atomic_json(args.output, artifact)
    print(
        f"complete candidates={len(candidates)} rank={rank} output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
