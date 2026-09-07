#!/usr/bin/env python3
"""One bounded skew-height/chart pass for Nagao ``u=135/2``.

This is the second application of the parameter-independent skew engine.  It
first replays the pinned uniform height-1,000,000 checkpoint, then searches
the same ten skew boxes used for ``u=42``.  To keep the chart budget no larger
than the first application, only the 38 largest-height unexpected checkpoint
abscissas are used, with the same two determinant-one orientations and the
same transformed height 50,000 (76 charts total).

All point memberships are exact.  Height ranks at 72 and 120 decimal digits
are search evidence.  If the augmented rank stays 17, PARI-proposed relations
are replayed using exact Fraction group arithmetic.  If it reaches at least
18, the exact selected point subset is stored for a subsequent finite-field
independence certificate.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import shlex
import sys
from typing import Any, Sequence

from ek_k3 import rational_to_string
from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK13_CONSTRUCTION,
    quartic_point_to_short_jacobian,
    quartic_value,
)
from nagao_skew_height import (
    MOBIUS_HEIGHT,
    SEARCH_BOXES,
    UniformCheckpoint,
    build_mobius_chart_plan,
    checkpoint_reference,
    classify_uniform_checkpoint,
    discover_relations,
    exact_linear_combination,
    load_rank17_target,
    map_chart_point,
    run_mobius_charts,
    run_skew_box,
)
from pari_bridge import pari_version
from search_extra_points import signless_quartic_points
from triage_nagao_rank13_finalists import (
    bounded_quartic_points,
    height_matrix_replay,
    point_digest,
    point_on_short_curve,
    stable_height_rank,
)


Q = Fraction
PARAMETER_U = Q(135, 2)
UNIFORM_HEIGHT = 1_000_000
CHART_CENTER_COUNT = 38
CHART_SHIFTS = (0, -1)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_u135_skew_height.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def projective_naive_height(point: tuple[Fraction, Fraction]) -> int:
    x_value = point[0]
    return max(abs(x_value.numerator), x_value.denominator)


def select_chart_centers(
    points: Sequence[tuple[Fraction, Fraction]], *, count: int
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Choose the largest projective-height centres deterministically."""

    if count <= 0:
        raise ValueError("the chart-centre count must be positive")
    unique: dict[Fraction, tuple[Fraction, Fraction]] = {}
    for point in points:
        unique.setdefault(Q(point[0]), (Q(point[0]), Q(point[1])))
    ordered = sorted(
        unique.values(),
        key=lambda point: (
            projective_naive_height(point),
            abs(point[0].numerator),
            point[0].denominator,
            point[0].numerator,
        ),
        reverse=True,
    )
    return tuple(ordered[:count])


def explicit_subset_records(
    pool: Sequence[tuple[Fraction, Fraction]], indices: Sequence[int]
) -> tuple[dict[str, Any], ...]:
    answer = []
    for index in indices:
        point = pool[index - 1]
        answer.append(
            {
                "pool_index_one_based": int(index),
                "jacobian_x": rational_to_string(point[0]),
                "jacobian_y": rational_to_string(point[1]),
                "exact_jacobian_membership_checked": True,
            }
        )
    return tuple(answer)


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts" / "generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--certificate-input",
        type=Path,
        default=generated / "elliptic_nagao_rank17_frontier_certificate.json",
    )
    parser.add_argument(
        "--rank-gain-input",
        type=Path,
        default=generated / "elliptic_nagao_rank13_rank_gain_search.json",
    )
    parser.add_argument("--uniform-timeout", type=float, default=90.0)
    parser.add_argument("--box-timeout", type=float, default=20.0)
    parser.add_argument("--chart-timeout", type=float, default=25.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_nagao_u135_skew_height.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 0 < args.uniform_timeout <= 90:
        raise SystemExit("--uniform-timeout must be in (0,90]")
    if not 0 < args.box_timeout <= 20:
        raise SystemExit("--box-timeout must be in (0,20]")
    if not 0 < args.chart_timeout <= 25:
        raise SystemExit("--chart-timeout must be in (0,25]")
    if not 0 < args.height_timeout <= 30:
        raise SystemExit("--height-timeout must be in (0,30]")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes is too small")

    target = load_rank17_target(args.certificate_input, PARAMETER_U)
    reference = checkpoint_reference(
        args.rank_gain_input, PARAMETER_U, height_bound=UNIFORM_HEIGHT
    )
    raw_uniform, uniform_wall, uniform_milliseconds = bounded_quartic_points(
        target.parameter_t,
        height_bound=UNIFORM_HEIGHT,
        timeout=args.uniform_timeout,
        stack_bytes=args.stack_bytes,
    )
    checkpoint = classify_uniform_checkpoint(
        target, raw_uniform, height_bound=UNIFORM_HEIGHT
    )
    if len(checkpoint.raw_signed_points) != reference.signed_point_count:
        raise AssertionError("the uniform signed-point count changed")
    if len(checkpoint.signless_points) != reference.distinct_x_count:
        raise AssertionError("the uniform distinct-x count changed")
    if len(checkpoint.unexpected_points) != reference.unexpected_x_count:
        raise AssertionError("the uniform unexpected-point count changed")
    if point_digest(checkpoint.unexpected_points) != reference.unexpected_point_sha256:
        raise AssertionError("the uniform unexpected-point digest changed")
    print(
        f"uniform H={UNIFORM_HEIGHT}: x={len(checkpoint.signless_points)} "
        f"unexpected={len(checkpoint.unexpected_points)}",
        flush=True,
    )

    checkpoint_x = {point[0] for point in checkpoint.signless_points}
    discovered: dict[Fraction, tuple[tuple[Fraction, Fraction], list[str]]] = {}
    box_records = []
    for search_box in SEARCH_BOXES:
        raw_points, pari_milliseconds, wall_seconds = run_skew_box(
            target.quartic_coefficients,
            search_box,
            timeout=args.box_timeout,
            stack_bytes=args.stack_bytes,
        )
        signless = signless_quartic_points(raw_points)
        outside = []
        for point in signless:
            if point[1] ** 2 != quartic_value(target.quartic_coefficients, point[0]):
                raise AssertionError("a skew-box point missed the exact quartic")
            if point[0] in checkpoint_x:
                continue
            outside.append(point[0])
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
                    rational_to_string(value) for value in outside
                ],
                "pari_reported_milliseconds": pari_milliseconds,
                "wall_seconds": wall_seconds,
            }
        )
        print(
            f"{search_box.identifier}: x={len(signless)} outside={len(outside)}",
            flush=True,
        )

    selected_centers = select_chart_centers(
        checkpoint.unexpected_points, count=CHART_CENTER_COUNT
    )
    limited_checkpoint = UniformCheckpoint(
        target=target,
        height_bound=checkpoint.height_bound,
        raw_signed_points=checkpoint.raw_signed_points,
        signless_points=checkpoint.signless_points,
        unexpected_points=selected_centers,
        displayed_x_count=checkpoint.displayed_x_count,
        companion_x_count=checkpoint.companion_x_count,
        zero_ordinate_count=checkpoint.zero_ordinate_count,
    )
    charts = build_mobius_chart_plan(limited_checkpoint, shifts=CHART_SHIFTS)
    if len(charts) != 76:
        raise AssertionError("the chart budget must equal the u=42 budget")
    chart_points, chart_milliseconds, chart_wall_seconds = run_mobius_charts(
        target.quartic_coefficients,
        charts,
        height_bound=MOBIUS_HEIGHT,
        timeout=args.chart_timeout,
        stack_bytes=args.stack_bytes,
    )
    chart_records = []
    for identifier, matrix in charts:
        raw_points = chart_points[identifier]
        mapped_by_x: dict[Fraction, tuple[Fraction, Fraction]] = {}
        for transformed_point in signless_quartic_points(raw_points):
            mapped = map_chart_point(transformed_point, matrix)
            if mapped is None:
                continue
            if mapped[1] ** 2 != quartic_value(
                target.quartic_coefficients, mapped[0]
            ):
                raise AssertionError("a chart point missed the original quartic")
            mapped_by_x.setdefault(mapped[0], mapped)
        outside = []
        for x_value, point in mapped_by_x.items():
            if x_value in checkpoint_x:
                continue
            outside.append(x_value)
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
        f"mobius charts: charts={len(charts)} "
        f"new_total_x={len(discovered)}",
        flush=True,
    )

    quartic_discoveries = tuple(discovered.values())
    quartic_points = tuple(item[0] for item in quartic_discoveries)
    mapped_points = tuple(
        quartic_point_to_short_jacobian(
            RANK13_CONSTRUCTION, target.parameter_t, quartic_point
        )
        for quartic_point in quartic_points
    )
    if any(
        not point_on_short_curve(target.jacobian_coefficients, point)
        for point in mapped_points
    ):
        raise AssertionError("a new image missed the exact Jacobian")

    basis_x = {point[0] for point in target.saturated_basis}
    distinct_images: list[tuple[Fraction, Fraction]] = []
    image_source_indices: list[int] = []
    seen_image_x = set(basis_x)
    for index, point in enumerate(mapped_points):
        if point[0] in seen_image_x:
            continue
        seen_image_x.add(point[0])
        distinct_images.append(point)
        image_source_indices.append(index)
    pool = target.saturated_basis + tuple(distinct_images)
    height_runs = height_matrix_replay(
        target.jacobian_coefficients,
        pool,
        precisions=(72, 120),
        timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    stable_rank = stable_height_rank(height_runs)
    subset_indices = height_runs[-1]["subset_indices_one_based"]
    selected_subset = explicit_subset_records(pool, subset_indices)

    relations: tuple[tuple[int, ...], ...] | None = None
    relation_status: str
    if stable_rank == target.certified_rank_lower_bound:
        relations = discover_relations(
            target.jacobian_coefficients,
            target.saturated_basis,
            tuple(distinct_images),
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
        )
        coefficient_a = target.jacobian_coefficients[3]
        for point, relation in zip(distinct_images, relations):
            if (
                exact_linear_combination(
                    coefficient_a, target.saturated_basis, relation
                )
                != point
            ):
                raise AssertionError("an exact basis relation failed")
        relation_status = "all distinct new images exactly in certified rank-17 span"
    else:
        relation_status = (
            "not attempted because the stable numerical pool rank exceeds the "
            "certified rank-17 basis"
        )

    point_records = []
    relation_by_source: dict[int, tuple[int, ...]] = {}
    if relations is not None:
        relation_by_source = dict(zip(image_source_indices, relations))
    for index, ((quartic_point, sources), jacobian_point) in enumerate(
        zip(quartic_discoveries, mapped_points)
    ):
        relation = relation_by_source.get(index)
        point_records.append(
            {
                "quartic_x": rational_to_string(quartic_point[0]),
                "quartic_z": rational_to_string(quartic_point[1]),
                "jacobian_x": rational_to_string(jacobian_point[0]),
                "jacobian_y": rational_to_string(jacobian_point[1]),
                "discovered_in": list(sources),
                "saturated_basis_relation": (
                    list(relation) if relation is not None else None
                ),
                "exact_quartic_membership_checked": True,
                "exact_jacobian_membership_checked": True,
                "exact_relation_replayed_with_fraction_group_law": (
                    relation is not None
                ),
            }
        )

    script_path = Path(__file__).resolve()
    actual_command = " ".join(
        shlex.quote(part) for part in [sys.executable, *sys.argv]
    )
    breakthrough = stable_rank >= 18
    artifact = {
        "schema_version": 1,
        "status": "bounded_search_complete",
        "candidate": {
            "id": target.identifier,
            "parameter_u": rational_to_string(target.parameter_u),
            "parameter_t": rational_to_string(target.parameter_t),
            "certified_algebraic_rank_lower_bound": (
                target.certified_rank_lower_bound
            ),
            "conductor": target.conductor,
            "log_conductor": target.log_conductor,
            "root_number": target.root_number,
        },
        "primary_source": PRIMARY_SOURCE,
        "inputs": {
            str(args.certificate_input): sha256_file(args.certificate_input),
            str(args.rank_gain_input): sha256_file(args.rank_gain_input),
        },
        "declared_budget": {
            "uniform_checkpoint_height": UNIFORM_HEIGHT,
            "uniform_checkpoint_timeout_seconds": args.uniform_timeout,
            "skew_box_count": len(SEARCH_BOXES),
            "skew_boxes_identical_to_u42": True,
            "skew_box_timeout_seconds_each": args.box_timeout,
            "chart_center_selection": (
                "38 largest projective-height unexpected H=10^6 abscissas"
            ),
            "chart_center_count": CHART_CENTER_COUNT,
            "chart_shifts": list(CHART_SHIFTS),
            "chart_count": len(charts),
            "chart_count_equal_to_u42": True,
            "chart_height": MOBIUS_HEIGHT,
            "chart_batch_timeout_seconds": args.chart_timeout,
            "height_and_relation_timeout_seconds_each": args.height_timeout,
            "stack_bytes": args.stack_bytes,
            "one_pass_no_retry": True,
        },
        "uniform_checkpoint_replay": {
            "signed_points_found": len(checkpoint.raw_signed_points),
            "distinct_quartic_x_values": len(checkpoint.signless_points),
            "unexpected_nonzero_quartic_x_values": len(
                checkpoint.unexpected_points
            ),
            "unexpected_point_sha256": point_digest(
                checkpoint.unexpected_points
            ),
            "reference_sha256_matched": True,
            "pari_reported_milliseconds": uniform_milliseconds,
            "wall_seconds": uniform_wall,
        },
        "selected_chart_centers": [
            {
                "quartic_x": rational_to_string(point[0]),
                "projective_naive_height": projective_naive_height(point),
            }
            for point in selected_centers
        ],
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
            "distinct_quartic_x_count": len(quartic_points),
            "distinct_new_jacobian_image_count": len(distinct_images),
            "all_mapped_and_checked_exactly": True,
            "relation_status": relation_status,
            "points": point_records,
        },
        "height_matrix_runs": height_runs,
        "stable_pool_numerical_rank": stable_rank,
        "stable_rank_at_least_18_observed": breakthrough,
        "explicit_stable_numerically_independent_subset": selected_subset,
        "subset_point_count": len(selected_subset),
        "interpretation": {
            "exact": (
                "all returned point memberships are exact; any displayed subgroup "
                "relations were replayed with exact Fraction group arithmetic"
            ),
            "numerical": (
                f"the augmented pool has stable height-matrix rank {stable_rank} "
                "at 72 and 120 decimal digits"
            ),
            "not_claimed": (
                "numerical height rank alone is not a Mordell-Weil independence "
                "certificate; a rank increase requires the planned finite-field "
                "certificate"
            ),
        },
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "platform": platform.platform(),
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "actual_command": actual_command,
        "script_sha256": sha256_file(script_path),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        f"wrote {args.output}: outside_x={len(quartic_points)} "
        f"new_images={len(distinct_images)} stable_rank={stable_rank}",
        flush=True,
    )


if __name__ == "__main__":
    main()
