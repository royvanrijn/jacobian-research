#!/usr/bin/env python3
"""Search for a twenty-first direction on Nagao's section-7 rank-20 fiber.

The exact input certificate concerns the six-root construction
``(346,260,255,146,55,0)`` at constructor parameter ``T=5081/47``.  This pass
uses three complementary bounded searches:

* every one of the ``2^20-1`` nonzero classes represented by the certified
  basis is scored in Gray-code order, using a streaming frontier;
* the best and weight-diverse alternate degree-two covers receive optimized
  cross-ratio charts in three height stages; and
* the original primitive quartic receives a deeper uniform box, the standard
  ten skew denominator slabs, and deeper replays of the six productive
  cross-ratio charts from the input certificate.

All subprocess calls use fresh process groups, strict timeouts, one attempt,
and joined TERM/KILL cleanup.  Every returned point is mapped and checked over
``QQ``.  Height-pairing relations count only after exact ``Fraction`` group
law replay.  Any unresolved point is immediately sent to the exact finite-
reduction independence engine; only that engine can record a rank gain.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import heapq
import json
from pathlib import Path
import platform
import re
import shlex
import sys
import time
from typing import Any, Iterable, Sequence

from alternate_quartic_covers import alternate_cover, short_subset_sum
from certify_nagao_rank20_t5081 import (
    CONSTRUCTION,
    EXPECTED_CONDUCTOR,
EXPECTED_SATURATED_BASIS_SHA256,
    PARAMETER_T,
    ROOTS,
    exact_curve_data,
)
from certify_nagao_rank21_t956 import gp_curve, gp_point, point_digest, signless_points
from ek_k3 import rational_to_string
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)
from nagao_1994 import (
    PRIMARY_SOURCE,
    primitive_quartic_coefficients,
    quartic_point_to_short_jacobian,
    quartic_value,
)
from pari_bridge import pari_version
from search_nagao_rank21_t6793_skew import map_run_points
from search_nagao_rank21_t956_skew import run_gp_once, search_original_quartic
from search_nagao_u135_alternate_covers import (
    CoverPlan,
    best_cross_ratio_charts,
    point_record,
    projective_height,
)
from search_nagao_u42_skew_height import (
    SEARCH_BOXES,
    exact_linear_combination,
    map_chart_point,
    short_add as fast_short_add,
    transform_binary_quartic,
)
from triage_nagao_rank13_finalists import point_on_short_curve


Q = Fraction
INPUT_RANK = 20
TARGET_RANK = 21
FULL_FRONTIER_RETAIN = 64
GLOBAL_FRONTIER_SELECTED = 40
MAX_SELECTED_COVERS = 60
CHARTS_PER_COVER = 3
PILOT_HEIGHT = 50_000
ESCALATION_HEIGHT = 250_000
DEEP_HEIGHT = 500_000
ESCALATION_CHART_COUNT = 16
DEEP_CHART_COUNT = 4
DIRECT_UNIFORM_HEIGHT = 1_500_000
DIRECT_CHART_ESCALATION_HEIGHT = 250_000
DIRECT_CHART_DEEP_HEIGHT = 500_000
DIRECT_CHART_COUNT = 6
DIRECT_DEEP_CHART_COUNT = 2
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank20_t5081_direction.py "
    "--output artifacts/generated-results/elliptic_nagao_rank20_t5081_direction.json"
)
EXPECTED_INPUT_CERTIFICATE_SHA256 = (
    "466946076dc0c3fa02d0c5edd90b947d5ee3d10a4fb8cb16567049ab4380f88d"
)

# Nagao's twelve displayed Mestre sections are not the complete set of
# low-degree sections on this section-7 quartic.  These six further sections
# were recovered by an exact linear-abscissa ansatz.  The first contributes
# the twelfth generic direction; the other five specialize dependently.  In
# either case they are predeclared search seeds, never exceptional fixed-fiber
# discoveries.
GENERIC_COMPANION_SECTIONS = (
    ("plus-7/27", Q(7, 27), Q(6920, 27), (Q(84770), Q(-18554923, 243), Q(-29974, 243), Q(680, 243))),
    ("minus-7/27", Q(-7, 27), Q(6920, 27), (Q(-84770), Q(-18554923, 243), Q(29974, 243), Q(680, 243))),
    ("plus-17/27", Q(17, 27), Q(5462, 27), (Q(5138284), Q(-6202747, 243), Q(-23222, 243), Q(440, 243))),
    ("minus-17/27", Q(-17, 27), Q(5462, 27), (Q(-5138284), Q(-6202747, 243), Q(23222, 243), Q(440, 243))),
    ("plus-43/27", Q(43, 27), Q(-4015, 27), (Q(94091525), Q(-236806588, 243), Q(756284, 243), Q(-1120, 243))),
    ("minus-43/27", Q(-43, 27), Q(-4015, 27), (Q(-94091525), Q(-236806588, 243), Q(-756284, 243), Q(-1120, 243))),
)


@dataclass(frozen=True)
class FrontierEntry:
    maximum_bit_length: int
    sum_bit_lengths: int
    mask: int

    @property
    def score(self) -> tuple[int, int, int]:
        return self.maximum_bit_length, self.sum_bit_lengths, self.mask

    @property
    def subset_indices(self) -> tuple[int, ...]:
        return tuple(index for index in range(self.mask.bit_length()) if self.mask >> index & 1)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generic_companion_quartic_points(
    parameter: Fraction, quartic: Sequence[Fraction]
) -> tuple[tuple[str, tuple[Fraction, Fraction]], ...]:
    """Specialize and exactly verify the six predeclared companion sections."""

    records = []
    for label, slope, intercept, ordinate_coefficients in GENERIC_COMPANION_SECTIONS:
        x_value = slope * parameter + intercept
        y_value = Q(0)
        for coefficient in reversed(ordinate_coefficients):
            y_value = y_value * parameter + coefficient
        point = x_value, y_value
        if y_value**2 != quartic_value(quartic, x_value):
            raise AssertionError(f"generic companion {label} failed exact specialization")
        records.append((label, point))
    if len({point[0] for _, point in records}) != len(records):
        raise AssertionError("generic companion abscissas collided")
    return tuple(records)


def load_exact_basis(
    path: Path,
) -> tuple[tuple[Fraction, ...], tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    if sha256_file(path) != EXPECTED_INPUT_CERTIFICATE_SHA256:
        raise AssertionError("the pinned input certificate file changed")
    data = json.loads(path.read_text(encoding="utf-8"))
    candidate = data["candidate"]
    certificate = data["exact_rank_certificate"]
    if Q(candidate["constructor_parameter_T"]) != PARAMETER_T:
        raise AssertionError("the input constructor parameter changed")
    if tuple(candidate["roots"]) != ROOTS:
        raise AssertionError("the input root tuple changed")
    if int(candidate["conductor"]) != EXPECTED_CONDUCTOR:
        raise AssertionError("the exact conductor changed")
    if candidate["root_number"] != 1 or not candidate["below_strict_log_conductor_target"]:
        raise AssertionError("the input conductor target replay changed")
    if (
        certificate["certified_algebraic_rank_lower_bound"] != INPUT_RANK
        or certificate["combined_exact_rank_over_F2"] != INPUT_RANK
        or certificate["saturated_basis_sha256"]
        != EXPECTED_SATURATED_BASIS_SHA256
    ):
        raise AssertionError("the input rank-20 certificate changed")
    coefficients = tuple(Q(value) for value in candidate["short_weierstrass_coefficients"])
    expected_coefficients = exact_curve_data()[2]
    if coefficients != expected_coefficients:
        raise AssertionError("the input short model changed")
    basis = tuple(
        (Q(record["jacobian_x"]), Q(record["jacobian_y"]))
        for record in certificate["saturated_basis"]
    )
    if (
        len(basis) != INPUT_RANK
        or point_digest(basis) != EXPECTED_SATURATED_BASIS_SHA256
        or any(not point_on_short_curve(coefficients, point) for point in basis)
    ):
        raise AssertionError("the exact rank-20 basis is invalid")
    return coefficients, basis, data


def _push_frontier(
    heap: list[tuple[int, int, int, FrontierEntry]],
    entry: FrontierEntry,
    retain_count: int,
) -> None:
    item = (
        -entry.maximum_bit_length,
        -entry.sum_bit_lengths,
        -entry.mask,
        entry,
    )
    if len(heap) < retain_count:
        heapq.heappush(heap, item)
        return
    worst = heap[0][3]
    if entry.score < worst.score:
        heapq.heapreplace(heap, item)


def streaming_full_coset_frontier(
    coefficients: Sequence[Fraction],
    basis: Sequence[tuple[Fraction, Fraction]],
    *,
    retain_count: int = FULL_FRONTIER_RETAIN,
    progress_interval: int = 65_536,
) -> tuple[tuple[FrontierEntry, ...], tuple[FrontierEntry, ...], dict[str, Any]]:
    """Score all represented nonzero mod-2 classes with bounded memory."""

    if retain_count <= 0 or not basis:
        raise ValueError("a nonempty basis and positive retain count are required")
    basis_size = len(basis)
    coefficient_a = Q(coefficients[3])
    current = None
    previous_gray = 0
    frontier_heap: list[tuple[int, int, int, FrontierEntry]] = []
    weight_best: dict[int, FrontierEntry] = {}
    started = time.monotonic()
    total = (1 << basis_size) - 1
    for integer in range(1, total + 1):
        gray = integer ^ (integer >> 1)
        changed = gray ^ previous_gray
        if changed == 0 or changed & (changed - 1):
            raise AssertionError("the Gray-code invariant failed")
        index = changed.bit_length() - 1
        point = basis[index]
        if not (gray >> index & 1):
            point = point[0], -point[1]
        current = fast_short_add(coefficient_a, current, point)
        if current is None:
            raise AssertionError("a nonzero represented class summed to infinity")
        x_base, y_base = current
        bit_lengths = []
        for x_value, y_value in basis:
            if (x_value, y_value) == current:
                # The cover base point is the point at infinity in this
                # coordinate and therefore has no finite parameter to score.
                continue
            if x_value == x_base:
                if y_value != -y_base or y_base == 0:
                    raise AssertionError("an unexpected exceptional cover coordinate")
                parameter = -(3 * x_base**2 + coefficient_a) / (2 * y_base)
            else:
                parameter = (y_value + y_base) / (x_value - x_base)
            bit_lengths.append(projective_height(parameter).bit_length())
        if not bit_lengths:
            raise AssertionError("a cover has no finite certified parameters")
        entry = FrontierEntry(max(bit_lengths), sum(bit_lengths), gray)
        _push_frontier(frontier_heap, entry, retain_count)
        weight = gray.bit_count()
        prior = weight_best.get(weight)
        if prior is None or entry.score < prior.score:
            weight_best[weight] = entry
        previous_gray = gray
        if progress_interval and integer % progress_interval == 0:
            print(
                f"coset-score {integer}/{total} elapsed={time.monotonic()-started:.1f}s",
                flush=True,
            )
    frontier = tuple(sorted((item[3] for item in frontier_heap), key=lambda item: item.score))
    weights = tuple(weight_best[weight] for weight in sorted(weight_best))
    if len(frontier) != min(retain_count, total) or len(weights) != basis_size:
        raise AssertionError("the streaming frontier dimensions changed")
    return frontier, weights, {
        "status": "completed",
        "basis_size": basis_size,
        "nonzero_classes_scored": total,
        "streaming_frontier_retained": len(frontier),
        "weight_frontier_count": len(weights),
        "wall_seconds": time.monotonic() - started,
    }


def selected_frontier_entries(
    frontier: Sequence[FrontierEntry],
    weight_frontier: Sequence[FrontierEntry],
    *,
    global_count: int = GLOBAL_FRONTIER_SELECTED,
    maximum_count: int = MAX_SELECTED_COVERS,
) -> tuple[FrontierEntry, ...]:
    if global_count <= 0 or maximum_count < global_count:
        raise ValueError("invalid frontier selection counts")
    selected: dict[int, FrontierEntry] = {
        entry.mask: entry for entry in frontier[:global_count]
    }
    for entry in weight_frontier:
        selected.setdefault(entry.mask, entry)
    return tuple(sorted(selected.values(), key=lambda entry: entry.score)[:maximum_count])


def run_cover_chart(
    plan: CoverPlan,
    chart: Any,
    *,
    stage: str,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[dict[str, Any], tuple[tuple[Fraction, Fraction], ...]]:
    transformed = transform_binary_quartic(plan.cover.coefficients, chart.matrix)
    raw, process = search_original_quartic(
        transformed,
        str(height_bound),
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    mapped = []
    poles = 0
    for transformed_point in raw:
        cover_point = map_chart_point(transformed_point, chart.matrix)
        if cover_point is None:
            poles += 1
            continue
        curve_point = plan.cover.cover_point_to_curve(cover_point)
        if not point_on_short_curve(plan.cover.short_coefficients, curve_point):
            raise AssertionError("an alternate-cover point left the exact curve")
        mapped.append(curve_point)
    return {
        "stage": stage,
        "height_bound": height_bound,
        "status": process["status"],
        "signed_transformed_points": len(raw),
        "points_at_chart_pole": poles,
        "finite_exact_curve_points": len(mapped),
        **{key: value for key, value in process.items() if key not in {"status"}},
    }, tuple(mapped)


def exact_relation_proposals(
    coefficients: Sequence[Fraction],
    basis: Sequence[tuple[Fraction, Fraction]],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    timeout: float,
    stack_bytes: int,
    batch_size: int = 40,
) -> tuple[tuple[tuple[int, ...] | None, bool, str], ...]:
    """Propose relations in capped batches and replay every success exactly."""

    answers: list[tuple[tuple[int, ...] | None, bool, str]] = []
    for start in range(0, len(points), batch_size):
        batch = points[start : start + batch_size]
        commands = [
            "default(realprecision,120);",
            f"E=ellinit([{gp_curve(coefficients)}]);",
            f"B=[{','.join(gp_point(point) for point in basis)}];",
            "H=ellheightmatrix(E,B);",
        ]
        for index, point in enumerate(batch):
            commands.extend(
                (
                    f"Q={gp_point(point)};V=vector(#B,j,ellheight(E,B[j],Q))~;C=round(matsolve(H,V));",
                    "S=[0];for(j=1,#B,S=elladd(E,S,ellmul(E,B[j],C[j])));",
                    f'print("RELATION_{index} ",Vec(C)," EXACT ",S==Q);',
                )
            )
        commands.append("quit")
        output, process = run_gp_once(
            "\n".join(commands) + "\n", timeout=timeout, stack_bytes=stack_bytes
        )
        if output is None:
            answers.extend((None, False, process["status"]) for _ in batch)
            continue
        for index, point in enumerate(batch):
            match = re.search(
                rf"^RELATION_{index} \[(.*?)\] EXACT ([01])$",
                output,
                re.MULTILINE,
            )
            if match is None:
                answers.append((None, False, "missing_output"))
                continue
            relation = tuple(int(value.strip()) for value in match.group(1).split(","))
            exact = match.group(2) == "1"
            if exact and exact_linear_combination(Q(coefficients[3]), basis, relation) != point:
                raise AssertionError("a proposed relation failed exact Fraction replay")
            answers.append((relation, exact, "completed"))
    if len(answers) != len(points):
        raise AssertionError("relation batching lost a point")
    return tuple(answers)


def signature_records(signatures: Sequence[Any]) -> list[dict[str, Any]]:
    return [
        {
            "prime": signature.prime,
            "group_order": signature.group_order,
            "doubled_subgroup_order": signature.doubled_subgroup_order,
            "quotient_dimension": signature.quotient_dimension,
            "rows": [list(row) for row in signature.rows],
        }
        for signature in signatures
    ]


def build_search(args: argparse.Namespace) -> dict[str, Any]:
    coefficients, basis, certificate_data = load_exact_basis(args.certificate_input)
    baseline_signatures = find_mod2_reduction_certificate(
        coefficients, basis, prime_bound=500
    )
    baseline_rank = combined_mod2_rank(baseline_signatures, len(basis))
    if baseline_rank != INPUT_RANK:
        raise AssertionError("the input basis lost its exact finite-reduction rank")

    frontier, weights, scan = streaming_full_coset_frontier(
        coefficients,
        basis,
        retain_count=args.frontier_retain,
        progress_interval=args.progress_interval,
    )
    selected_entries = selected_frontier_entries(
        frontier,
        weights,
        global_count=args.global_frontier_selected,
        maximum_count=args.maximum_selected_covers,
    )
    plans = []
    entry_by_mask = {entry.mask: entry for entry in selected_entries}
    for entry in selected_entries:
        base_point = short_subset_sum(coefficients, basis, entry.subset_indices)
        if base_point is None:
            raise AssertionError("a selected nonzero cover class vanished")
        cover = alternate_cover(coefficients, base_point)
        plans.append(
            CoverPlan(
                entry.subset_indices,
                cover,
                best_cross_ratio_charts(cover, basis, count=CHARTS_PER_COVER),
            )
        )
    print(f"selected covers={len(plans)}", flush=True)

    discoveries: dict[tuple[Fraction, Fraction], set[str]] = {}
    cover_runs: list[dict[str, Any]] = []
    pilot_yields = []

    def absorb_cover(plan: CoverPlan, chart: Any, stage: str, record: dict[str, Any], points: Iterable[tuple[Fraction, Fraction]]) -> int:
        before = len(discoveries)
        source = f"alternate:{plan.identifier}:{stage}:{chart.basis_indices}"
        for point in points:
            discoveries.setdefault(point, set()).add(source)
        gained = len(discoveries) - before
        record.update(
            {
                "cover_id": plan.identifier,
                "cover_subset_indices_one_based": [index + 1 for index in plan.subset_indices],
                "normalizing_basis_indices_one_based": [index + 1 for index in chart.basis_indices],
                "matrix_a_b_c_d": list(chart.matrix),
                "new_global_exact_affine_points": gained,
            }
        )
        cover_runs.append(record)
        return gained

    for plan_index, plan in enumerate(plans, start=1):
        for chart in plan.charts:
            record, points = run_cover_chart(
                plan,
                chart,
                stage="pilot",
                height_bound=args.pilot_height,
                timeout=args.pilot_timeout,
                stack_bytes=args.stack_bytes,
            )
            gained = absorb_cover(plan, chart, "pilot", record, points)
            pilot_yields.append((gained, plan, chart))
        print(
            f"cover {plan_index}/{len(plans)} {plan.identifier} discoveries={len(discoveries)}",
            flush=True,
        )
    pilot_yields.sort(key=lambda item: (-item[0], item[1].score, item[2].score))
    escalation_yields = []
    for _, plan, chart in pilot_yields[: args.escalation_charts]:
        record, points = run_cover_chart(
            plan,
            chart,
            stage="escalation",
            height_bound=args.escalation_height,
            timeout=args.escalation_timeout,
            stack_bytes=args.stack_bytes,
        )
        gained = absorb_cover(plan, chart, "escalation", record, points)
        escalation_yields.append((gained, plan, chart))
    escalation_yields.sort(key=lambda item: (-item[0], item[1].score, item[2].score))
    for _, plan, chart in escalation_yields[: args.deep_charts]:
        record, points = run_cover_chart(
            plan,
            chart,
            stage="deep",
            height_bound=args.deep_height,
            timeout=args.deep_timeout,
            stack_bytes=args.stack_bytes,
        )
        absorb_cover(plan, chart, "deep", record, points)

    quartic = primitive_quartic_coefficients(CONSTRUCTION, PARAMETER_T)
    _, visible_jacobian, _ = exact_curve_data()
    companion_quartic = generic_companion_quartic_points(PARAMETER_T, quartic)
    companion_jacobian = tuple(
        (
            label,
            quartic_point_to_short_jacobian(CONSTRUCTION, PARAMETER_T, quartic_point),
        )
        for label, quartic_point in companion_quartic
    )
    if any(
        not point_on_short_curve(coefficients, point)
        for _, point in companion_jacobian
    ):
        raise AssertionError("a generic companion image missed the exact Jacobian")
    generic_seed_points = visible_jacobian + tuple(point for _, point in companion_jacobian)
    generic_seed_proposals = exact_relation_proposals(
        coefficients,
        basis,
        generic_seed_points,
        timeout=args.relation_timeout,
        stack_bytes=args.stack_bytes,
        batch_size=args.relation_batch_size,
    )
    if not all(exact for _, exact, _ in generic_seed_proposals):
        raise AssertionError("a predeclared generic seed was not replayed in the rank-20 basis")
    direct_runs: list[dict[str, Any]] = []
    direct_quartic_by_x: dict[Fraction, tuple[Fraction, Fraction]] = {}

    def absorb_direct(
        points: Iterable[tuple[Fraction, Fraction]], source: str, record: dict[str, Any]
    ) -> None:
        before = len(direct_quartic_by_x)
        for point in points:
            if point[1] ** 2 != quartic_value(quartic, point[0]):
                raise AssertionError("a direct point missed the primitive quartic")
            direct_quartic_by_x.setdefault(point[0], point)
        record["source"] = source
        record["new_global_direct_quartic_abscissas"] = len(direct_quartic_by_x) - before
        direct_runs.append(record)

    raw, process = search_original_quartic(
        quartic,
        str(args.direct_uniform_height),
        timeout=args.direct_uniform_timeout,
        stack_bytes=args.stack_bytes,
    )
    absorb_direct(signless_points(raw), "uniform", process)
    print(f"direct uniform status={process['status']} x={len(direct_quartic_by_x)}", flush=True)
    for box in SEARCH_BOXES:
        raw, process = search_original_quartic(
            quartic,
            box.gp_height,
            timeout=args.skew_timeout,
            stack_bytes=args.stack_bytes,
        )
        record = {
            "id": box.identifier,
            "numerator_absolute_bound": box.numerator_bound,
            "denominator_lower_bound": box.denominator_lower,
            "denominator_upper_bound": box.denominator_upper,
            **process,
        }
        absorb_direct(signless_points(raw), f"skew:{box.identifier}", record)
        print(f"direct {box.identifier} status={process['status']} x={len(direct_quartic_by_x)}", flush=True)

    productive = sorted(
        (
            record
            for record in certificate_data["point_search"]["chart_records"]
            if record["new_global_abscissas"] > 0
        ),
        key=lambda record: (-record["new_global_abscissas"], record["chart_id"]),
    )[:DIRECT_CHART_COUNT]
    direct_chart_yields = []
    for chart_record in productive:
        matrix = tuple(int(value) for value in chart_record["matrix"])
        transformed = transform_binary_quartic(quartic, matrix)
        raw, process = search_original_quartic(
            transformed,
            str(args.direct_chart_escalation_height),
            timeout=args.direct_chart_timeout,
            stack_bytes=args.stack_bytes,
        )
        mapped = map_run_points(quartic, raw, matrix)
        before = len(direct_quartic_by_x)
        record = {
            "id": chart_record["chart_id"],
            "stage": "escalation",
            "matrix_a_b_c_d": list(matrix),
            "height_bound": args.direct_chart_escalation_height,
            **process,
        }
        absorb_direct(mapped, f"direct-chart:{chart_record['chart_id']}:escalation", record)
        direct_chart_yields.append((len(direct_quartic_by_x) - before, chart_record, matrix))
    direct_chart_yields.sort(key=lambda item: (-item[0], item[1]["chart_id"]))
    for _, chart_record, matrix in direct_chart_yields[:DIRECT_DEEP_CHART_COUNT]:
        transformed = transform_binary_quartic(quartic, matrix)
        raw, process = search_original_quartic(
            transformed,
            str(args.direct_chart_deep_height),
            timeout=args.direct_deep_timeout,
            stack_bytes=args.stack_bytes,
        )
        mapped = map_run_points(quartic, raw, matrix)
        record = {
            "id": chart_record["chart_id"],
            "stage": "deep",
            "matrix_a_b_c_d": list(matrix),
            "height_bound": args.direct_chart_deep_height,
            **process,
        }
        absorb_direct(mapped, f"direct-chart:{chart_record['chart_id']}:deep", record)

    direct_image_sources: dict[tuple[Fraction, Fraction], set[str]] = {}
    for quartic_point in direct_quartic_by_x.values():
        if quartic_point[1] == 0:
            continue
        image = quartic_point_to_short_jacobian(CONSTRUCTION, PARAMETER_T, quartic_point)
        if not point_on_short_curve(coefficients, image):
            raise AssertionError("a direct quartic image missed the Jacobian")
        direct_image_sources.setdefault(image, set()).add("direct-original-quartic")
    for point, sources in direct_image_sources.items():
        discoveries.setdefault(point, set()).update(sources)

    predeclared_with_signs = {
        point
        for basis_point in basis + generic_seed_points
        for point in (basis_point, (basis_point[0], -basis_point[1]))
    }
    candidates = tuple(
        sorted(
            (point for point in discoveries if point not in predeclared_with_signs),
            key=lambda point: (projective_height(point[0]), projective_height(point[1]), point),
        )
    )
    proposals = exact_relation_proposals(
        coefficients,
        basis,
        candidates,
        timeout=args.relation_timeout,
        stack_bytes=args.stack_bytes,
        batch_size=args.relation_batch_size,
    )
    unresolved = tuple(
        point for point, (_, exact, _) in zip(candidates, proposals) if not exact
    )
    augmented_signatures = find_mod2_reduction_certificate(
        coefficients, basis + unresolved, prime_bound=1_000
    )
    augmented_rank = combined_mod2_rank(
        augmented_signatures, len(basis) + len(unresolved)
    )
    certified_gain = max(0, augmented_rank - baseline_rank)
    certified_rank = baseline_rank + certified_gain

    candidate_records = []
    for point, (relation, exact, status) in zip(candidates, proposals):
        candidate_records.append(
            {
                **point_record(point),
                "sources": sorted(discoveries[point]),
                "relation_process_status": status,
                "exact_relation_in_certified_rank20_subgroup": exact,
                "basis_relation": list(relation) if exact and relation is not None else None,
                "fraction_group_law_replay": exact,
            }
        )

    cover_records = []
    for plan in plans:
        mask = sum(1 << index for index in plan.subset_indices)
        entry = entry_by_mask[mask]
        cover_records.append(
            {
                "id": plan.identifier,
                "mask": hex(mask),
                "subset_weight": len(plan.subset_indices),
                "identity_score_maximum_known_t_projective_bit_length": entry.maximum_bit_length,
                "identity_score_sum_known_t_projective_bit_lengths": entry.sum_bit_lengths,
                "subset_indices_one_based": [index + 1 for index in plan.subset_indices],
                "base_point": point_record(plan.cover.base_point),
                "quartic_coefficients_ascending": [
                    rational_to_string(value) for value in plan.cover.coefficients
                ],
                "cross_ratio_charts": [
                    {
                        "normalizing_basis_indices_one_based": [index + 1 for index in chart.basis_indices],
                        "matrix_a_b_c_d": list(chart.matrix),
                        "mean_log10_known_projective_height": chart.mean_log_height,
                        "median_log10_known_projective_height": chart.median_log_height,
                        "maximum_log10_known_projective_height": chart.maximum_log_height,
                    }
                    for chart in plan.charts
                ],
            }
        )

    script_path = Path(__file__).resolve()
    return {
        "schema_version": 1,
        "status": "bounded section-7 rank-20 direction search complete",
        "candidate": {
            "constructor_parameter_T": rational_to_string(PARAMETER_T),
            "roots": list(ROOTS),
            "conductor": str(EXPECTED_CONDUCTOR),
            "log_conductor": certificate_data["candidate"]["log_conductor"],
            "root_number": certificate_data["candidate"]["root_number"],
            "certified_rank_lower_bound_before_search": baseline_rank,
            "target_rank": TARGET_RANK,
        },
        "input": {
            "path": str(args.certificate_input),
            "sha256": sha256_file(args.certificate_input),
            "saturated_basis_sha256": point_digest(basis),
        },
        "generic_seed_decontamination": {
            "displayed_mestre_section_count": len(visible_jacobian),
            "companion_section_count": len(companion_jacobian),
            "all_specialized_images_replayed_in_certified_rank20_basis": True,
            "predeclared_images_found_by_search_with_sign": sum(
                point in predeclared_with_signs for point in discoveries
            ),
            "companions": [
                {
                    "label": label,
                    "quartic_point": {
                        "x": rational_to_string(quartic_point[0]),
                        "y": rational_to_string(quartic_point[1]),
                    },
                    "jacobian_point": point_record(jacobian_point),
                    "basis_relation": list(proposal[0]) if proposal[0] is not None else None,
                    "fraction_group_law_replay": proposal[1],
                }
                for (label, quartic_point), (_, jacobian_point), proposal in zip(
                    companion_quartic,
                    companion_jacobian,
                    generic_seed_proposals[len(visible_jacobian) :],
                )
            ],
        },
        "full_mod2_class_scan": {
            **scan,
            "frontier": [entry.__dict__ | {"mask_hex": hex(entry.mask)} for entry in frontier],
            "weight_frontier": [entry.__dict__ | {"mask_hex": hex(entry.mask), "weight": entry.mask.bit_count()} for entry in weights],
            "global_frontier_selected": args.global_frontier_selected,
            "selected_cover_count_after_weight_diversification": len(plans),
        },
        "cover_plans": cover_records,
        "alternate_cover_search": {
            "pilot_height": args.pilot_height,
            "pilot_chart_count": len(plans) * CHARTS_PER_COVER,
            "escalation_height": args.escalation_height,
            "escalation_chart_count": min(args.escalation_charts, len(pilot_yields)),
            "deep_height": args.deep_height,
            "deep_chart_count": min(args.deep_charts, len(escalation_yields)),
            "runs": cover_runs,
        },
        "direct_primitive_quartic_search": {
            "uniform_height": args.direct_uniform_height,
            "standard_skew_box_count": len(SEARCH_BOXES),
            "productive_input_charts_deepened": len(productive),
            "distinct_quartic_abscissas": len(direct_quartic_by_x),
            "runs": direct_runs,
        },
        "results": {
            "distinct_exact_affine_curve_points": len(discoveries),
            "nonbasis_candidate_points": len(candidates),
            "candidate_point_sha256": point_digest(candidates),
            "exact_relations_in_certified_rank20_subgroup": sum(exact for _, exact, _ in proposals),
            "unresolved_by_exact_relation_replay": len(unresolved),
            "augmented_finite_reduction_signatures": signature_records(augmented_signatures),
            "combined_exact_finite_reduction_rank": augmented_rank,
            "certified_new_directions": certified_gain,
            "certified_rank_lower_bound_after_search": certified_rank,
            "target_rank_21_achieved": certified_rank >= TARGET_RANK,
            "candidate_points": candidate_records,
        },
        "bounded_scope": {
            "one_pass_no_retry": True,
            "fresh_foreground_process_group_per_call": True,
            "all_subprocess_timeouts_at_most_60_seconds": True,
            "negative_search_is_not_a_rank_upper_bound": True,
        },
        "primary_source": PRIMARY_SOURCE,
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "platform": platform.platform(),
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "actual_command": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "script_sha256": sha256_file(script_path),
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts/generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--certificate-input",
        type=Path,
        default=generated / "elliptic_nagao_rank20_t5081_rank20_certificate.json",
    )
    parser.add_argument("--frontier-retain", type=int, default=FULL_FRONTIER_RETAIN)
    parser.add_argument("--global-frontier-selected", type=int, default=GLOBAL_FRONTIER_SELECTED)
    parser.add_argument("--maximum-selected-covers", type=int, default=MAX_SELECTED_COVERS)
    parser.add_argument("--progress-interval", type=int, default=65_536)
    parser.add_argument("--pilot-height", type=int, default=PILOT_HEIGHT)
    parser.add_argument("--pilot-timeout", type=float, default=8.0)
    parser.add_argument("--escalation-height", type=int, default=ESCALATION_HEIGHT)
    parser.add_argument("--escalation-charts", type=int, default=ESCALATION_CHART_COUNT)
    parser.add_argument("--escalation-timeout", type=float, default=15.0)
    parser.add_argument("--deep-height", type=int, default=DEEP_HEIGHT)
    parser.add_argument("--deep-charts", type=int, default=DEEP_CHART_COUNT)
    parser.add_argument("--deep-timeout", type=float, default=40.0)
    parser.add_argument("--direct-uniform-height", type=int, default=DIRECT_UNIFORM_HEIGHT)
    parser.add_argument("--direct-uniform-timeout", type=float, default=60.0)
    parser.add_argument("--skew-timeout", type=float, default=25.0)
    parser.add_argument("--direct-chart-escalation-height", type=int, default=DIRECT_CHART_ESCALATION_HEIGHT)
    parser.add_argument("--direct-chart-timeout", type=float, default=15.0)
    parser.add_argument("--direct-chart-deep-height", type=int, default=DIRECT_CHART_DEEP_HEIGHT)
    parser.add_argument("--direct-deep-timeout", type=float, default=40.0)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--relation-batch-size", type=int, default=40)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_nagao_rank20_t5081_direction.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    timeout_names = (
        "pilot_timeout",
        "escalation_timeout",
        "deep_timeout",
        "direct_uniform_timeout",
        "skew_timeout",
        "direct_chart_timeout",
        "direct_deep_timeout",
        "relation_timeout",
    )
    if any(not 0 < getattr(args, name) <= 60 for name in timeout_names):
        raise SystemExit("all subprocess timeouts must lie in (0,60]")
    positive_names = (
        "frontier_retain",
        "global_frontier_selected",
        "maximum_selected_covers",
        "progress_interval",
        "pilot_height",
        "escalation_height",
        "escalation_charts",
        "deep_height",
        "deep_charts",
        "direct_uniform_height",
        "direct_chart_escalation_height",
        "direct_chart_deep_height",
        "relation_batch_size",
    )
    if any(getattr(args, name) <= 0 for name in positive_names):
        raise SystemExit("all counts and height bounds must be positive")
    if not args.global_frontier_selected <= args.frontier_retain:
        raise SystemExit("global selection exceeds the retained frontier")
    if args.maximum_selected_covers < args.global_frontier_selected:
        raise SystemExit("maximum cover count is below the global selection")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes is too small")
    result = build_search(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"wrote {args.output}: candidates={result['results']['nonbasis_candidate_points']} "
        f"unresolved={result['results']['unresolved_by_exact_relation_replay']} "
        f"certified_rank={result['results']['certified_rank_lower_bound_after_search']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
