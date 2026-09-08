#!/usr/bin/env python3
"""Deep skew and Mobius search on the Nagao ``T=6793/64`` quartic.

The input artifact certifies nineteen independent points on the short
Jacobian.  This standalone pass reconstructs the primitive quartic, replays
the old uniform height-one-million checkpoint, and then searches regions not
covered by that box:

* the ten denominator slabs of ``SEARCH_BOXES``;
* exact three-point cross-ratio charts selected by how strongly they compress
  all checkpoint abscissas; and
* centred determinant-one charts around the highest checkpoint and new skew
  abscissas.

The best-yielding pilot charts are escalated through two larger transformed
height bounds.  Each PARI invocation is a one-shot foreground process group
with a strict timeout and no retry.  Every returned point is mapped back to
the original quartic and Jacobian and checked over ``QQ``.  Height pairings
are selection evidence only.  A numerical rank gain triggers saturation and
an exact finite-reduction independence certificate.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import shlex
import statistics
import sys
from typing import Any, Iterable, Sequence

from alternate_quartic_covers import (
    mobius_preimage,
    point_on_short_curve,
    three_point_mobius_matrix,
)
from certify_nagao_rank21_t956 import (
    height_replay,
    point_digest,
    saturate_basis,
    signless_points,
)
from ek_k3 import rational_to_string
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK21_CONSTRUCTION,
    primitive_quartic_coefficients,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    quartic_value,
    short_jacobian_coefficients,
)
from search_nagao_rank21_t956_skew import (
    search_original_quartic,
    signature_records,
)
from search_nagao_u135_alternate_covers import relation_proposals
from search_nagao_u42_skew_height import (
    SEARCH_BOXES,
    SearchBox,
    centered_unimodular_matrix,
    map_chart_point,
    transform_binary_quartic,
)


Q = Fraction
PARAMETER_T = Q(6793, 64)
INPUT_BASIS_DIGEST = (
    "cdb7328683b523f49aac5efe3588631e82e9f18dd603b1329bc4e9f7c89e44dd"
)
INPUT_RANK = 19
TARGET_RANK = 21
CHECKPOINT_HEIGHT = 1_000_000
PILOT_HEIGHT = 50_000
ESCALATION_HEIGHT = 250_000
DEEP_HEIGHT = 1_000_000
CROSS_RATIO_CHART_COUNT = 48
CENTRED_CHECKPOINT_COUNT = 12
CENTRED_SKEW_COUNT = 8
CENTRED_SHIFTS = (0, -1, 1)
ESCALATION_CHART_COUNT = 16
DEEP_CHART_COUNT = 4
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank21_t6793_skew.py "
    "--output archive/elliptic-curves/artifacts/generated-results/elliptic_nagao_rank21_t6793_skew.json"
)


@dataclass(frozen=True)
class ChartPlan:
    identifier: str
    kind: str
    matrix: tuple[int, int, int, int]
    mean_log10_known_height: float
    median_log10_known_height: float
    maximum_log10_known_height: float

    @property
    def score(self) -> tuple[float, float, float, str]:
        return (
            self.mean_log10_known_height,
            self.median_log10_known_height,
            self.maximum_log10_known_height,
            self.identifier,
        )


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def projective_height(value: Fraction) -> int:
    value = Q(value)
    return max(abs(value.numerator), value.denominator)


def load_exact_basis(
    path: Path,
) -> tuple[tuple[Fraction, ...], tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    identifier = "unbiased-6793-64"
    certificate = data["finite_reduction_certificates"][identifier]
    if certificate["certified_algebraic_rank_lower_bound"] != INPUT_RANK:
        raise AssertionError("the input does not certify rank at least nineteen")
    if certificate["combined_exact_rank_over_F2"] != INPUT_RANK:
        raise AssertionError("the input finite-reduction rank changed")
    if certificate["saturated_point_sha256"] != INPUT_BASIS_DIGEST:
        raise AssertionError("the pinned saturated-basis digest changed")
    basis = tuple(
        (Q(record["jacobian_x"]), Q(record["jacobian_y"]))
        for record in certificate["saturation"]["saturated_basis"]
    )
    coefficients = short_jacobian_coefficients(RANK21_CONSTRUCTION, PARAMETER_T)
    if len(basis) != INPUT_RANK or point_digest(basis) != INPUT_BASIS_DIGEST:
        raise AssertionError("the exact input basis is incomplete")
    if any(not point_on_short_curve(coefficients, point) for point in basis):
        raise AssertionError("an input basis point is off the exact curve")
    conductor = data["final_conductor_replays"][identifier]
    if conductor["status"] != "completed" or conductor["root_number"] != -1:
        raise AssertionError("the exact conductor replay changed")
    return coefficients, basis, conductor


def chart_height_statistics(
    matrix: Sequence[int], known_x_values: Sequence[Fraction]
) -> tuple[float, float, float]:
    logarithms = [
        math.log10(projective_height(preimage))
        for value in known_x_values
        if (preimage := mobius_preimage(matrix, value)) is not None
    ]
    if not logarithms:
        raise AssertionError("a chart has no finite known preimages")
    return (
        sum(logarithms) / len(logarithms),
        statistics.median(logarithms),
        max(logarithms),
    )


def optimized_cross_ratio_charts(
    known_x_values: Sequence[Fraction], *, count: int
) -> tuple[ChartPlan, ...]:
    """Return the exact three-point charts best compressing known abscissas."""

    values = tuple(sorted(set(Q(value) for value in known_x_values)))
    if len(values) < 3 or count <= 0:
        raise ValueError("at least three known values and a positive count are required")
    candidates: list[ChartPlan] = []
    seen: set[tuple[int, int, int, int]] = set()
    # Cyclic orderings suffice: reversing zero and one merely applies s -> 1-s,
    # already represented by the same symmetric height box.
    for indices in itertools.combinations(range(len(values)), 3):
        selected = tuple(values[index] for index in indices)
        for rotation in range(3):
            images = selected[rotation:] + selected[:rotation]
            matrix = three_point_mobius_matrix(*images)
            if matrix in seen:
                continue
            seen.add(matrix)
            mean, median, maximum = chart_height_statistics(matrix, values)
            candidates.append(
                ChartPlan(
                    identifier=(
                        f"cross_{indices[0]:02d}_{indices[1]:02d}_{indices[2]:02d}"
                        f"_rot_{rotation}"
                    ),
                    kind="three_point_cross_ratio",
                    matrix=matrix,
                    mean_log10_known_height=mean,
                    median_log10_known_height=median,
                    maximum_log10_known_height=maximum,
                )
            )
    candidates.sort(key=lambda chart: chart.score)
    if len(candidates) < count:
        raise AssertionError("too few distinct optimized charts")
    return tuple(candidates[:count])


def centred_charts(
    x_values: Sequence[Fraction], *, prefix: str, count: int
) -> tuple[ChartPlan, ...]:
    values = sorted(
        set(Q(value) for value in x_values),
        key=lambda value: (-projective_height(value), value),
    )[:count]
    plans = []
    for index, center in enumerate(values):
        for shift in CENTRED_SHIFTS:
            matrix = centered_unimodular_matrix(center, shift)
            mean, median, maximum = chart_height_statistics(matrix, values)
            plans.append(
                ChartPlan(
                    identifier=f"{prefix}_{index:02d}_shift_{shift}",
                    kind="centred_unimodular",
                    matrix=matrix,
                    mean_log10_known_height=mean,
                    median_log10_known_height=median,
                    maximum_log10_known_height=maximum,
                )
            )
    return tuple(plans)


def deduplicate_chart_plans(plans: Iterable[ChartPlan]) -> tuple[ChartPlan, ...]:
    by_matrix: dict[tuple[int, int, int, int], ChartPlan] = {}
    for plan in plans:
        prior = by_matrix.get(plan.matrix)
        if prior is None or plan.score < prior.score:
            by_matrix[plan.matrix] = plan
    return tuple(sorted(by_matrix.values(), key=lambda plan: (plan.kind, plan.score)))


def map_run_points(
    quartic: Sequence[Fraction],
    raw_points: Iterable[tuple[Fraction, Fraction]],
    matrix: Sequence[int] | None,
) -> tuple[tuple[Fraction, Fraction], ...]:
    mapped: dict[Fraction, tuple[Fraction, Fraction]] = {}
    for raw in signless_points(raw_points):
        point = raw if matrix is None else map_chart_point(raw, matrix)
        if point is None:
            continue
        if point[1] ** 2 != quartic_value(quartic, point[0]):
            raise AssertionError("a returned point missed the original quartic")
        mapped.setdefault(point[0], point)
    return tuple(mapped.values())


def run_chart(
    quartic: Sequence[Fraction],
    plan: ChartPlan,
    *,
    stage: str,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    transformed = transform_binary_quartic(quartic, plan.matrix)
    raw, process = search_original_quartic(
        transformed,
        str(height_bound),
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    mapped = map_run_points(quartic, raw, plan.matrix)
    return mapped, {
        "id": plan.identifier,
        "kind": plan.kind,
        "stage": stage,
        "matrix_a_b_c_d": list(plan.matrix),
        "determinant": plan.matrix[0] * plan.matrix[3] - plan.matrix[1] * plan.matrix[2],
        "transformed_height_bound": height_bound,
        "mean_log10_known_projective_height": plan.mean_log10_known_height,
        "median_log10_known_projective_height": plan.median_log10_known_height,
        "maximum_log10_known_projective_height": plan.maximum_log10_known_height,
        **process,
        "distinct_finite_original_abscissa_count": len(mapped),
    }


def point_record(
    quartic_point: tuple[Fraction, Fraction],
    jacobian_point: tuple[Fraction, Fraction],
    sources: Sequence[str],
    relation: tuple[int, ...] | None,
) -> dict[str, Any]:
    return {
        "quartic_x": rational_to_string(quartic_point[0]),
        "quartic_z": rational_to_string(quartic_point[1]),
        "jacobian_x": rational_to_string(jacobian_point[0]),
        "jacobian_y": rational_to_string(jacobian_point[1]),
        "sources": list(sources),
        "basis_relation": list(relation) if relation is not None else None,
        "relation_replayed_exactly": relation is not None,
        "exact_quartic_and_jacobian_membership_checked": True,
    }


def exact_gain_certificate(
    coefficients: Sequence[Fraction],
    pool: Sequence[tuple[Fraction, Fraction]],
    height_runs: Sequence[dict[str, Any]],
    *,
    timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    stable_rank = int(height_runs[-1]["numerical_rank"])
    if stable_rank < 20:
        return {"status": "not_triggered", "stable_numerical_rank": stable_rank}
    selected = tuple(
        pool[index - 1]
        for index in height_runs[-1]["subset_indices_one_based"]
    )
    saturated, saturation = saturate_basis(
        coefficients, selected, timeout=timeout, stack_bytes=stack_bytes
    )
    signatures = find_mod2_reduction_certificate(
        coefficients, saturated, prime_bound=1_000
    )
    exact_rank = combined_mod2_rank(signatures, len(saturated))
    result: dict[str, Any] = {
        "status": "certified" if exact_rank == len(saturated) else "rank_deficient",
        "small_prime_saturation": saturation,
        "saturated_basis": [
            {
                "jacobian_x": rational_to_string(point[0]),
                "jacobian_y": rational_to_string(point[1]),
                "exact_membership_checked": True,
            }
            for point in saturated
        ],
        "saturated_basis_sha256": point_digest(saturated),
        "finite_reduction_prime_bound": 1_000,
        "finite_reduction_signatures": signature_records(signatures),
        "combined_exact_rank_over_F2": exact_rank,
    }
    if exact_rank == len(saturated):
        result["no_rational_2_torsion_prime"] = find_two_torsion_certificate_prime(
            coefficients
        )
        result["certified_algebraic_rank_lower_bound"] = exact_rank
        result["target_rank_21_achieved"] = exact_rank >= TARGET_RANK
    return result


def build_search(args: argparse.Namespace) -> dict[str, Any]:
    coefficients, basis, conductor = load_exact_basis(args.certificate_input)
    quartic = primitive_quartic_coefficients(RANK21_CONSTRUCTION, PARAMETER_T)
    baseline_signatures = find_mod2_reduction_certificate(
        coefficients, basis, prime_bound=1_000
    )
    if combined_mod2_rank(baseline_signatures, len(basis)) != INPUT_RANK:
        raise AssertionError("the input basis no longer independently certifies rank 19")

    checkpoint_raw, checkpoint_process = search_original_quartic(
        quartic,
        str(CHECKPOINT_HEIGHT),
        timeout=args.checkpoint_timeout,
        stack_bytes=args.stack_bytes,
    )
    if checkpoint_process["status"] != "completed":
        raise RuntimeError("the uniform checkpoint did not complete")
    checkpoint_points = map_run_points(quartic, checkpoint_raw, None)
    if len(checkpoint_points) != 31:
        raise AssertionError("the pinned uniform checkpoint no longer has 31 abscissas")
    visible = primitive_visible_points(RANK21_CONSTRUCTION, PARAMETER_T)
    visible_x = {point[0] for point in visible}
    if len(visible_x) != 12 or not visible_x.issubset(
        {point[0] for point in checkpoint_points}
    ):
        raise AssertionError("the uniform checkpoint lost a visible section")

    known_x = {point[0] for point in checkpoint_points}
    discovered: dict[Fraction, dict[str, Any]] = {}
    box_records = []
    skew_new_points: list[tuple[Fraction, Fraction]] = []
    for search_box in SEARCH_BOXES:
        raw, process = search_original_quartic(
            quartic,
            search_box.gp_height,
            timeout=args.box_timeout,
            stack_bytes=args.stack_bytes,
        )
        mapped = map_run_points(quartic, raw, None)
        outside = []
        for point in mapped:
            if point[0] in known_x:
                continue
            outside.append(point[0])
            skew_new_points.append(point)
            entry = discovered.setdefault(point[0], {"point": point, "sources": []})
            entry["sources"].append(f"skew:{search_box.identifier}")
        box_records.append(
            {
                "id": search_box.identifier,
                "numerator_absolute_bound": search_box.numerator_bound,
                "denominator_lower_bound": search_box.denominator_lower,
                "denominator_upper_bound": search_box.denominator_upper,
                **process,
                "outside_checkpoint_x_values": [
                    rational_to_string(value) for value in outside
                ],
            }
        )
        print(
            f"{search_box.identifier}: status={process['status']} outside={len(outside)}",
            flush=True,
        )

    checkpoint_nonvisible = [
        point[0] for point in checkpoint_points if point[0] not in visible_x
    ]
    plans = deduplicate_chart_plans(
        itertools.chain(
            optimized_cross_ratio_charts(
                tuple(known_x), count=args.cross_ratio_charts
            ),
            centred_charts(
                checkpoint_nonvisible,
                prefix="checkpoint",
                count=args.centred_checkpoint_count,
            ),
            centred_charts(
                [point[0] for point in skew_new_points],
                prefix="skew",
                count=args.centred_skew_count,
            ),
        )
    )
    chart_records: list[dict[str, Any]] = []
    pilot_yields: list[tuple[int, ChartPlan]] = []

    def absorb(
        plan: ChartPlan,
        stage: str,
        points: Sequence[tuple[Fraction, Fraction]],
        record: dict[str, Any],
    ) -> int:
        before = len(discovered)
        outside = []
        for point in points:
            if point[0] in known_x:
                continue
            outside.append(point[0])
            entry = discovered.setdefault(point[0], {"point": point, "sources": []})
            entry["sources"].append(f"chart:{plan.identifier}:{stage}")
        gained = len(discovered) - before
        record["outside_checkpoint_x_values"] = [
            rational_to_string(value) for value in outside
        ]
        record["new_global_abscissas"] = gained
        chart_records.append(record)
        return gained

    for index, plan in enumerate(plans, start=1):
        points, record = run_chart(
            quartic,
            plan,
            stage="pilot",
            height_bound=args.pilot_height,
            timeout=args.pilot_timeout,
            stack_bytes=args.stack_bytes,
        )
        gained = absorb(plan, "pilot", points, record)
        pilot_yields.append((gained, plan))
        print(
            f"pilot {index}/{len(plans)} {plan.identifier}: "
            f"status={record['status']} gain={gained}",
            flush=True,
        )

    pilot_yields.sort(key=lambda item: (-item[0], item[1].score))
    escalation_yields: list[tuple[int, ChartPlan]] = []
    for _, plan in pilot_yields[: args.escalation_charts]:
        points, record = run_chart(
            quartic,
            plan,
            stage="escalation",
            height_bound=args.escalation_height,
            timeout=args.escalation_timeout,
            stack_bytes=args.stack_bytes,
        )
        gained = absorb(plan, "escalation", points, record)
        escalation_yields.append((gained, plan))
    escalation_yields.sort(key=lambda item: (-item[0], item[1].score))
    for _, plan in escalation_yields[: args.deep_charts]:
        points, record = run_chart(
            quartic,
            plan,
            stage="deep",
            height_bound=args.deep_height,
            timeout=args.deep_timeout,
            stack_bytes=args.stack_bytes,
        )
        absorb(plan, "deep", points, record)

    basis_signs = {
        point
        for basis_point in basis
        for point in (basis_point, (basis_point[0], -basis_point[1]))
    }
    image_by_x: dict[Fraction, tuple[Fraction, Fraction]] = {}
    source_by_image_x: dict[Fraction, set[str]] = {}
    quartic_by_image_x: dict[Fraction, tuple[Fraction, Fraction]] = {}
    for entry in discovered.values():
        quartic_point = entry["point"]
        image = quartic_point_to_short_jacobian(
            RANK21_CONSTRUCTION, PARAMETER_T, quartic_point
        )
        if not point_on_short_curve(coefficients, image):
            raise AssertionError("a discovered image missed the exact Jacobian")
        if image in basis_signs:
            continue
        image_by_x.setdefault(image[0], image)
        quartic_by_image_x.setdefault(image[0], quartic_point)
        source_by_image_x.setdefault(image[0], set()).update(entry["sources"])
    new_images = tuple(
        image_by_x[x_value]
        for x_value in sorted(image_by_x, key=lambda value: (projective_height(value), value))
    )
    pool = basis + new_images
    height_runs = height_replay(
        coefficients,
        pool,
        timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    stable_rank = int(height_runs[-1]["numerical_rank"])
    proposals = relation_proposals(
        coefficients,
        basis,
        new_images,
        timeout=args.relation_timeout,
        stack_bytes=args.stack_bytes,
    )
    records = []
    for image, (relation, exact) in zip(new_images, proposals):
        records.append(
            point_record(
                quartic_by_image_x[image[0]],
                image,
                sorted(source_by_image_x[image[0]]),
                relation if exact else None,
            )
        )
    exact_gain = exact_gain_certificate(
        coefficients,
        pool,
        height_runs,
        timeout=args.saturation_timeout,
        stack_bytes=args.stack_bytes,
    )
    completed_boxes = sum(record["status"] == "completed" for record in box_records)
    completed_charts = sum(
        record["status"] == "completed" for record in chart_records
    )
    timed_out_charts = sum(record["status"] == "timeout" for record in chart_records)
    all_declared_runs_completed = (
        completed_boxes == len(box_records) and completed_charts == len(chart_records)
    )
    script_path = Path(__file__).resolve()
    return {
        "schema_version": 1,
        "status": (
            "bounded deep/skew/Mobius search complete"
            if all_declared_runs_completed
            else "bounded deep/skew/Mobius pass complete with recorded timeouts"
        ),
        "candidate": {
            "parameter_t": rational_to_string(PARAMETER_T),
            "conductor": conductor["conductor"],
            "log_conductor": conductor["log_conductor"],
            "root_number": conductor["root_number"],
            "certified_rank_lower_bound_before_search": INPUT_RANK,
            "target_rank": TARGET_RANK,
        },
        "primary_source": PRIMARY_SOURCE,
        "input": {
            "path": str(args.certificate_input),
            "sha256": sha256_file(args.certificate_input),
            "saturated_basis_sha256": point_digest(basis),
            "baseline_finite_reduction_primes": [
                signature.prime for signature in baseline_signatures
            ],
        },
        "uniform_checkpoint": {
            **checkpoint_process,
            "distinct_abscissa_count": len(checkpoint_points),
            "visible_abscissa_count": len(visible_x),
            "replayed_exactly": True,
        },
        "skew_staircase": {"boxes": box_records},
        "mobius_search": {
            "cross_ratio_selection": (
                "all triples and three cyclic orientations ranked by the mean, median, and maximum log10 projective heights of exact known preimages"
            ),
            "pilot_height": args.pilot_height,
            "pilot_chart_count": len(plans),
            "escalation_height": args.escalation_height,
            "escalation_chart_count": min(args.escalation_charts, len(plans)),
            "deep_height": args.deep_height,
            "deep_chart_count": min(args.deep_charts, len(escalation_yields)),
            "records": chart_records,
        },
        "new_point_analysis": {
            "outside_checkpoint_quartic_abscissa_count": len(discovered),
            "distinct_nonbasis_jacobian_sign_pairs": len(new_images),
            "exact_relations_in_certified_rank19_span": sum(
                exact for _, exact in proposals
            ),
            "unresolved_by_exact_relation_replay": sum(
                not exact for _, exact in proposals
            ),
            "point_sha256": point_digest(new_images),
            "records": records,
        },
        "height_selection": {
            "pool_point_count": len(pool),
            "runs": list(height_runs),
            "stable_numerical_rank": stable_rank,
            "stable_numerical_rank_gain": stable_rank - INPUT_RANK,
            "selection_is_not_certification": True,
        },
        "exact_rank_gain_attempt": exact_gain,
        "bounded_scope": {
            "checkpoint_timeout_seconds": args.checkpoint_timeout,
            "box_timeout_seconds_each": args.box_timeout,
            "pilot_timeout_seconds_each": args.pilot_timeout,
            "escalation_timeout_seconds_each": args.escalation_timeout,
            "deep_timeout_seconds_each": args.deep_timeout,
            "height_timeout_seconds": args.height_timeout,
            "relation_timeout_seconds": args.relation_timeout,
            "saturation_timeout_seconds_if_triggered": args.saturation_timeout,
            "one_pass_no_retry": True,
            "fresh_foreground_process_group_per_call": True,
            "completed_box_calls": completed_boxes,
            "completed_chart_calls": completed_charts,
            "timed_out_chart_calls": timed_out_charts,
            "all_declared_runs_completed": all_declared_runs_completed,
        },
        "interpretation": {
            "exact": (
                "all stored quartic and Jacobian points were checked over QQ; stored relations were replayed by exact Fraction group arithmetic"
            ),
            "numerical": "height rank was replayed at 72 and 120 decimal digits",
            "bounded": (
                "every declared box and chart was attempted once; timed-out calls are explicitly recorded and are not exhausted; no result is a rank upper bound or a complete descent"
            ),
        },
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "actual_command": " ".join(
            shlex.quote(part) for part in [sys.executable, *sys.argv]
        ),
        "script_sha256": sha256_file(script_path),
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts/generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--certificate-input",
        type=Path,
        default=generated / "elliptic_nagao_rank21_unbiased.json",
    )
    parser.add_argument("--checkpoint-timeout", type=float, default=60.0)
    parser.add_argument("--box-timeout", type=float, default=25.0)
    parser.add_argument("--pilot-timeout", type=float, default=8.0)
    parser.add_argument("--escalation-timeout", type=float, default=15.0)
    parser.add_argument("--deep-timeout", type=float, default=30.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--cross-ratio-charts", type=int, default=CROSS_RATIO_CHART_COUNT)
    parser.add_argument(
        "--centred-checkpoint-count", type=int, default=CENTRED_CHECKPOINT_COUNT
    )
    parser.add_argument("--centred-skew-count", type=int, default=CENTRED_SKEW_COUNT)
    parser.add_argument("--pilot-height", type=int, default=PILOT_HEIGHT)
    parser.add_argument("--escalation-height", type=int, default=ESCALATION_HEIGHT)
    parser.add_argument("--deep-height", type=int, default=DEEP_HEIGHT)
    parser.add_argument("--escalation-charts", type=int, default=ESCALATION_CHART_COUNT)
    parser.add_argument("--deep-charts", type=int, default=DEEP_CHART_COUNT)
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_nagao_rank21_t6793_skew.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    timeout_names = (
        "checkpoint_timeout",
        "box_timeout",
        "pilot_timeout",
        "escalation_timeout",
        "deep_timeout",
        "height_timeout",
        "relation_timeout",
        "saturation_timeout",
    )
    if any(not 0 < getattr(args, name) <= 60 for name in timeout_names):
        raise SystemExit("all PARI timeouts must lie in (0,60]")
    count_names = (
        "cross_ratio_charts",
        "centred_checkpoint_count",
        "centred_skew_count",
        "pilot_height",
        "escalation_height",
        "deep_height",
        "escalation_charts",
        "deep_charts",
    )
    if any(getattr(args, name) <= 0 for name in count_names):
        raise SystemExit("all search counts and heights must be positive")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes is too small")
    result = build_search(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"wrote {args.output}: outside="
        f"{result['new_point_analysis']['outside_checkpoint_quartic_abscissa_count']} "
        f"images={result['new_point_analysis']['distinct_nonbasis_jacobian_sign_pairs']} "
        f"stable_rank={result['height_selection']['stable_numerical_rank']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
