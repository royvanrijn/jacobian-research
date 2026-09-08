#!/usr/bin/env python3
"""Bounded alternate degree-two chart search for Nagao ``u=135/2``.

The original Mestre quartic exposes one particular height on the curve.  This
search instead uses the exact degree-two models ``C_Q`` from
``alternate_quartic_covers``.  The base points ``Q`` are nonzero subset sums
of one or two members of the certified rank-17 basis.  Their finite-reduction
signatures show that the resulting ``Q mod 2E(Q)`` classes are distinct.

For each cover, three known cover parameters are sent to ``0,1,infinity``.
The remaining known parameters rank these cross-ratio charts by projective
height.  PARI's ``hyperellratpoints`` then searches a fixed small collection
of the best charts in three strictly bounded stages.  Every returned point is
mapped back and checked exactly.  Exact height-pairing proposals are replayed
with Fraction group arithmetic; any direction not so explained is tested by
the finite-reduction independence engine.

This is a bounded experiment, not a complete 2-descent or a rank upper bound.
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
import re
import shlex
import statistics
import sys
from typing import Any, Iterable, Sequence

from alternate_quartic_covers import (
    AlternateQuarticCover,
    alternate_cover,
    mobius_preimage,
    short_subset_sum,
    three_point_mobius_matrix,
)
from ek_k3 import rational_to_string
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)
from nagao_1994 import PRIMARY_SOURCE
from nagao_skew_height import load_rank17_target
from pari_bridge import pari_version
from search_extra_points import run_gp
from search_nagao_u42_skew_height import (
    exact_linear_combination,
    map_chart_point,
    run_mobius_charts,
    short_add as fast_short_add,
)
from triage_nagao_rank13_finalists import (
    gp_rational,
    gp_vector,
    point_digest,
    point_on_short_curve,
)


Q = Fraction
PARAMETER_U = Q(135, 2)
MAX_SUBSET_WEIGHT = 2
LOW_WEIGHT_COVER_COUNT = 12
FULL_COSET_TRANCHE_COUNT = 8
FULL_COSET_RETAIN_COUNT = 20
SELECTED_COVER_COUNT = LOW_WEIGHT_COVER_COUNT + FULL_COSET_TRANCHE_COUNT
CHARTS_PER_COVER = 3
PILOT_HEIGHT = 50_000
ESCALATION_HEIGHT = 250_000
ESCALATION_CHART_COUNT = 8
DEEP_HEIGHT = 1_000_000
DEEP_CHART_COUNT = 2
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_u135_alternate_covers.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def projective_height(value: Fraction) -> int:
    value = Q(value)
    return max(abs(value.numerator), value.denominator)


@dataclass(frozen=True)
class CrossRatioChart:
    basis_indices: tuple[int, int, int]
    matrix: tuple[int, int, int, int]
    mean_log_height: float
    median_log_height: float
    maximum_log_height: float

    @property
    def score(self) -> tuple[float, float, float, tuple[int, int, int]]:
        return (
            self.mean_log_height,
            self.median_log_height,
            self.maximum_log_height,
            self.basis_indices,
        )


@dataclass(frozen=True)
class CoverPlan:
    subset_indices: tuple[int, ...]
    cover: AlternateQuarticCover
    charts: tuple[CrossRatioChart, ...]

    @property
    def score(self) -> tuple[float, float, float, tuple[int, ...]]:
        best = self.charts[0]
        return (
            best.mean_log_height,
            best.median_log_height,
            best.maximum_log_height,
            self.subset_indices,
        )

    @property
    def identifier(self) -> str:
        return "q_" + "_".join(str(index + 1) for index in self.subset_indices)


def cover_parameters(
    cover: AlternateQuarticCover,
    basis: Sequence[tuple[Fraction, Fraction]],
) -> tuple[tuple[int, Fraction], ...]:
    """Return finite cover parameters of the certified basis points."""

    answer = []
    for index, point in enumerate(basis):
        try:
            parameter, _ = cover.curve_point_to_cover(point)
        except ValueError:
            # The chart base point Q itself is the only intended point at
            # infinity.  Pair sums are normally not individual basis points.
            continue
        answer.append((index, parameter))
    return tuple(answer)


def best_cross_ratio_charts(
    cover: AlternateQuarticCover,
    basis: Sequence[tuple[Fraction, Fraction]],
    *,
    count: int,
) -> tuple[CrossRatioChart, ...]:
    """Rank exact three-point charts by the other known coordinate heights."""

    if count <= 0:
        raise ValueError("the chart count must be positive")
    parameters = cover_parameters(cover, basis)
    candidates: list[CrossRatioChart] = []
    seen_matrices: set[tuple[int, int, int, int]] = set()
    for positions in itertools.combinations(range(len(parameters)), 3):
        selected = tuple(parameters[position] for position in positions)
        values = tuple(item[1] for item in selected)
        if len(set(values)) != 3:
            continue
        matrix = three_point_mobius_matrix(*values)
        if matrix in seen_matrices:
            continue
        seen_matrices.add(matrix)
        logarithms = []
        for _, parameter in parameters:
            preimage = mobius_preimage(matrix, parameter)
            if preimage is not None:
                logarithms.append(math.log10(projective_height(preimage)))
        if not logarithms:
            raise AssertionError("a chart has no finite known preimages")
        candidates.append(
            CrossRatioChart(
                tuple(item[0] for item in selected),
                matrix,
                sum(logarithms) / len(logarithms),
                statistics.median(logarithms),
                max(logarithms),
            )
        )
    candidates.sort(key=lambda chart: chart.score)
    if len(candidates) < count:
        raise AssertionError("too few distinct cross-ratio charts")
    return tuple(candidates[:count])


def enumerate_cover_plans(
    coefficients: Sequence[Fraction],
    basis: Sequence[tuple[Fraction, Fraction]],
    *,
    maximum_subset_weight: int = MAX_SUBSET_WEIGHT,
    charts_per_cover: int = CHARTS_PER_COVER,
) -> tuple[CoverPlan, ...]:
    """Build and deterministically rank all declared low-weight cover classes."""

    if maximum_subset_weight <= 0:
        raise ValueError("the maximum subset weight must be positive")
    plans = []
    for weight in range(1, maximum_subset_weight + 1):
        for indices in itertools.combinations(range(len(basis)), weight):
            base_point = short_subset_sum(coefficients, basis, indices)
            if base_point is None:
                raise AssertionError("a nonempty certified subset summed to infinity")
            cover = alternate_cover(coefficients, base_point)
            charts = best_cross_ratio_charts(
                cover, basis, count=charts_per_cover
            )
            plans.append(CoverPlan(indices, cover, charts))
    plans.sort(key=lambda plan: plan.score)
    return tuple(plans)


def mask_indices(mask: int, basis_size: int) -> tuple[int, ...]:
    if mask <= 0 or mask >= 1 << basis_size:
        raise ValueError("the subset mask is outside the nonzero basis range")
    return tuple(index for index in range(basis_size) if mask >> index & 1)


def full_coset_identity_frontier(
    coefficients: Sequence[Fraction],
    basis: Sequence[tuple[Fraction, Fraction]],
    *,
    retain_count: int,
) -> tuple[tuple[tuple[int, int, int], tuple[int, ...]], ...]:
    """Scan every known mod-2 class and retain the best identity coordinates.

    Gray-code order changes one basis bit at a time, so all ``2^17-1`` exact
    base points require one group addition apiece.  The score is entirely
    integral: first minimize the maximum projective bit length of the finite
    known ``t`` parameters, then their sum, then the subset mask.
    """

    if retain_count <= 0:
        raise ValueError("the retained full-coset count must be positive")
    basis_size = len(basis)
    if basis_size <= 0:
        raise ValueError("the basis must be nonempty")
    current = None
    previous_gray = 0
    scored: list[tuple[tuple[int, int, int], tuple[int, ...]]] = []
    for integer in range(1, 1 << basis_size):
        gray = integer ^ (integer >> 1)
        changed = gray ^ previous_gray
        if changed == 0 or changed & (changed - 1):
            raise AssertionError("consecutive Gray words must differ in one bit")
        index = changed.bit_length() - 1
        point = basis[index]
        if not (gray >> index & 1):
            point = point[0], -point[1]
        # The shared search group law omits the expensive point-on-curve cube
        # replay after every Gray step.  It is still exact Fraction arithmetic;
        # every retained point is validated when ``alternate_cover`` is built.
        current = fast_short_add(Q(coefficients[3]), current, point)
        if current is None:
            raise AssertionError("a nonempty certified subset summed to infinity")
        x_base, y_base = current
        coefficient_a = Q(coefficients[3])
        parameters = []
        for candidate in basis:
            x_value, y_value = candidate
            if candidate == current:
                continue
            if x_value == x_base:
                if y_value != -y_base or y_base == 0:
                    raise AssertionError("an unexpected exceptional chart point")
                parameter = -(3 * x_base**2 + coefficient_a) / (2 * y_base)
            else:
                parameter = (y_value + y_base) / (x_value - x_base)
            parameters.append(parameter)
        bit_lengths = tuple(
            projective_height(parameter).bit_length() for parameter in parameters
        )
        if not bit_lengths:
            raise AssertionError("a cover has no finite certified parameters")
        indices = mask_indices(gray, basis_size)
        scored.append(((max(bit_lengths), sum(bit_lengths), gray), indices))
        previous_gray = gray
    scored.sort(key=lambda item: item[0])
    return tuple(scored[:retain_count])


def relation_proposals(
    coefficients: Sequence[Fraction],
    basis: Sequence[tuple[Fraction, Fraction]],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[int, ...], bool], ...]:
    """Ask PARI for height-pairing relations and retain its exact checks."""

    if not points:
        return ()
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
    answer = []
    for index in range(len(points)):
        match = re.search(
            rf"^RELATION_{index} \[(.*?)\] EXACT ([01])$", output, re.MULTILINE
        )
        if match is None:
            raise AssertionError(f"PARI omitted relation proposal {index}")
        relation = tuple(int(value.strip()) for value in match.group(1).split(","))
        if len(relation) != len(basis):
            raise AssertionError("a relation proposal has the wrong length")
        exact = match.group(2) == "1"
        if exact:
            replay = exact_linear_combination(Q(coefficients[3]), basis, relation)
            if replay != points[index]:
                raise AssertionError("a PARI relation failed the Fraction replay")
        answer.append((relation, exact))
    return tuple(answer)


def point_record(point: tuple[Fraction, Fraction]) -> dict[str, str]:
    return {
        "curve_x": rational_to_string(point[0]),
        "curve_y": rational_to_string(point[1]),
    }


def run_chart(
    plan: CoverPlan,
    chart: CrossRatioChart,
    *,
    stage: str,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[dict[str, Any], tuple[tuple[Fraction, Fraction], ...]]:
    identifier = f"{plan.identifier}_{stage}"
    try:
        raw_by_chart, milliseconds, wall_seconds = run_mobius_charts(
            plan.cover.coefficients,
            ((identifier, chart.matrix),),
            height_bound=height_bound,
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
    except Exception as error:
        # A timeout or PARI failure is a bounded negative diagnostic.  Keep it
        # in the artifact and continue to the next independent chart.
        return (
            {
                "stage": stage,
                "height_bound": height_bound,
                "status": "failed_or_timed_out",
                "error_type": type(error).__name__,
                "error": str(error),
            },
            (),
        )

    curve_points: list[tuple[Fraction, Fraction]] = []
    pole_count = 0
    for transformed_point in raw_by_chart[identifier]:
        cover_point = map_chart_point(transformed_point, chart.matrix)
        if cover_point is None:
            pole_count += 1
            continue
        curve_point = plan.cover.cover_point_to_curve(cover_point)
        if not point_on_short_curve(plan.cover.short_coefficients, curve_point):
            raise AssertionError("a mapped alternate-cover point left E")
        curve_points.append(curve_point)
    return (
        {
            "stage": stage,
            "height_bound": height_bound,
            "status": "completed",
            "signed_transformed_points": len(raw_by_chart[identifier]),
            "points_at_chart_pole": pole_count,
            "finite_exact_curve_points": len(curve_points),
            "pari_reported_milliseconds": milliseconds[identifier],
            "wall_seconds": wall_seconds,
        },
        tuple(curve_points),
    )


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts/generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--certificate-input",
        type=Path,
        default=generated / "elliptic_nagao_rank17_frontier_certificate.json",
    )
    parser.add_argument("--pilot-timeout", type=float, default=8.0)
    parser.add_argument("--escalation-timeout", type=float, default=12.0)
    parser.add_argument("--deep-timeout", type=float, default=20.0)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_nagao_u135_alternate_covers.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    for name in (
        "pilot_timeout",
        "escalation_timeout",
        "deep_timeout",
        "relation_timeout",
    ):
        value = getattr(args, name)
        maximum = 60.0
        if not 0 < value <= maximum:
            raise SystemExit(f"--{name.replace('_', '-')} must be in (0,{maximum:g}]")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes is too small")

    target = load_rank17_target(args.certificate_input, PARAMETER_U)
    low_weight_plans = enumerate_cover_plans(
        target.jacobian_coefficients,
        target.saturated_basis,
        maximum_subset_weight=MAX_SUBSET_WEIGHT,
        charts_per_cover=CHARTS_PER_COVER,
    )
    expected_plan_count = sum(
        math.comb(len(target.saturated_basis), weight)
        for weight in range(1, MAX_SUBSET_WEIGHT + 1)
    )
    if len(low_weight_plans) != expected_plan_count:
        raise AssertionError("the low-weight cover plan count changed")
    full_coset_frontier = full_coset_identity_frontier(
        target.jacobian_coefficients,
        target.saturated_basis,
        retain_count=FULL_COSET_RETAIN_COUNT,
    )
    selected_by_indices: dict[tuple[int, ...], CoverPlan] = {
        plan.subset_indices: plan
        for plan in low_weight_plans[:LOW_WEIGHT_COVER_COUNT]
    }
    full_coset_selected_scores: dict[tuple[int, ...], tuple[int, int, int]] = {}
    for score, indices in full_coset_frontier:
        if indices in selected_by_indices:
            continue
        base_point = short_subset_sum(
            target.jacobian_coefficients, target.saturated_basis, indices
        )
        if base_point is None:
            raise AssertionError("a retained full-coset subset vanished")
        cover = alternate_cover(target.jacobian_coefficients, base_point)
        selected_by_indices[indices] = CoverPlan(
            indices,
            cover,
            best_cross_ratio_charts(
                cover, target.saturated_basis, count=CHARTS_PER_COVER
            ),
        )
        full_coset_selected_scores[indices] = score
        if len(full_coset_selected_scores) == FULL_COSET_TRANCHE_COUNT:
            break
    selected_plans = tuple(selected_by_indices.values())
    if len(selected_plans) != SELECTED_COVER_COUNT:
        raise AssertionError("the selected low-weight/full-coset union changed size")

    run_records: list[dict[str, Any]] = []
    discoveries: dict[tuple[Fraction, Fraction], set[str]] = {}
    chart_yields: list[tuple[int, CoverPlan, CrossRatioChart]] = []

    def absorb(
        plan: CoverPlan,
        chart: CrossRatioChart,
        stage: str,
        record: dict[str, Any],
        points: Iterable[tuple[Fraction, Fraction]],
    ) -> int:
        source = f"{plan.identifier}:{stage}:{chart.basis_indices}"
        before = len(discoveries)
        for point in points:
            discoveries.setdefault(point, set()).add(source)
        gained = len(discoveries) - before
        record.update(
            {
                "cover_id": plan.identifier,
                "cover_subset_indices_one_based": [
                    index + 1 for index in plan.subset_indices
                ],
                "normalizing_basis_indices_one_based": [
                    index + 1 for index in chart.basis_indices
                ],
                "matrix_a_b_c_d": list(chart.matrix),
                "new_global_exact_affine_points": gained,
            }
        )
        run_records.append(record)
        return gained

    for plan in selected_plans:
        for chart in plan.charts:
            record, points = run_chart(
                plan,
                chart,
                stage="pilot",
                height_bound=PILOT_HEIGHT,
                timeout=args.pilot_timeout,
                stack_bytes=args.stack_bytes,
            )
            gained = absorb(plan, chart, "pilot", record, points)
            chart_yields.append((gained, plan, chart))
        print(
            f"{plan.identifier}: pilot charts={len(plan.charts)} "
            f"global_points={len(discoveries)}",
            flush=True,
        )

    # Prefer charts that exposed new exact affine points, with the original
    # compression score as the deterministic tie breaker.
    chart_yields.sort(key=lambda item: (-item[0], item[1].score, item[2].score))
    escalated = chart_yields[:ESCALATION_CHART_COUNT]
    escalation_yields = []
    for _, plan, chart in escalated:
        record, points = run_chart(
            plan,
            chart,
            stage="escalation",
            height_bound=ESCALATION_HEIGHT,
            timeout=args.escalation_timeout,
            stack_bytes=args.stack_bytes,
        )
        gained = absorb(plan, chart, "escalation", record, points)
        escalation_yields.append((gained, plan, chart))

    escalation_yields.sort(
        key=lambda item: (-item[0], item[1].score, item[2].score)
    )
    for _, plan, chart in escalation_yields[:DEEP_CHART_COUNT]:
        record, points = run_chart(
            plan,
            chart,
            stage="deep",
            height_bound=DEEP_HEIGHT,
            timeout=args.deep_timeout,
            stack_bytes=args.stack_bytes,
        )
        absorb(plan, chart, "deep", record, points)

    basis_with_signs = {
        point
        for basis_point in target.saturated_basis
        for point in (basis_point, (basis_point[0], -basis_point[1]))
    }
    candidates = tuple(
        sorted(
            (point for point in discoveries if point not in basis_with_signs),
            key=lambda point: (
                projective_height(point[0]),
                projective_height(point[1]),
                point,
            ),
        )
    )
    proposals = relation_proposals(
        target.jacobian_coefficients,
        target.saturated_basis,
        candidates,
        timeout=args.relation_timeout,
        stack_bytes=args.stack_bytes,
    )
    unresolved = tuple(
        point for point, (_, exact) in zip(candidates, proposals) if not exact
    )

    finite_signatures = find_mod2_reduction_certificate(
        target.jacobian_coefficients,
        target.saturated_basis + unresolved,
        prime_bound=500,
    )
    finite_rank = combined_mod2_rank(
        finite_signatures, len(target.saturated_basis) + len(unresolved)
    )
    certified_gain = max(0, finite_rank - target.certified_rank_lower_bound)

    candidate_records = []
    for point, (relation, exact) in zip(candidates, proposals):
        record: dict[str, Any] = point_record(point)
        record.update(
            {
                "sources": sorted(discoveries[point]),
                "exact_relation_in_certified_rank17_subgroup": exact,
                "basis_relation": list(relation) if exact else None,
                "fraction_group_law_replay": exact,
            }
        )
        candidate_records.append(record)

    plan_records = []
    for plan in selected_plans:
        plan_records.append(
            {
                "id": plan.identifier,
                "subset_indices_one_based": [
                    index + 1 for index in plan.subset_indices
                ],
                "base_point": point_record(plan.cover.base_point),
                "quartic_coefficients_ascending": [
                    rational_to_string(value) for value in plan.cover.coefficients
                ],
                "forced_full_coset_identity_score": (
                    {
                        "maximum_known_t_projective_bit_length": (
                            full_coset_selected_scores[plan.subset_indices][0]
                        ),
                        "sum_known_t_projective_bit_lengths": (
                            full_coset_selected_scores[plan.subset_indices][1]
                        ),
                        "mask": hex(
                            full_coset_selected_scores[plan.subset_indices][2]
                        ),
                    }
                    if plan.subset_indices in full_coset_selected_scores
                    else None
                ),
                "charts": [
                    {
                        "normalizing_basis_indices_one_based": [
                            index + 1 for index in chart.basis_indices
                        ],
                        "matrix_a_b_c_d": list(chart.matrix),
                        "mean_log10_known_projective_height": (
                            chart.mean_log_height
                        ),
                        "median_log10_known_projective_height": (
                            chart.median_log_height
                        ),
                        "maximum_log10_known_projective_height": (
                            chart.maximum_log_height
                        ),
                    }
                    for chart in plan.charts
                ],
            }
        )

    script_path = Path(__file__).resolve()
    actual_command = " ".join(
        shlex.quote(part) for part in [sys.executable, *sys.argv]
    )
    artifact = {
        "schema_version": 1,
        "status": "bounded_alternate_cover_search_complete",
        "candidate": {
            "parameter_u": rational_to_string(target.parameter_u),
            "parameter_t": rational_to_string(target.parameter_t),
            "conductor": target.conductor,
            "log_conductor": target.log_conductor,
            "certified_rank_lower_bound_before_search": (
                target.certified_rank_lower_bound
            ),
        },
        "primary_source": PRIMARY_SOURCE,
        "input": {
            "path": str(args.certificate_input),
            "sha256": sha256_file(args.certificate_input),
        },
        "construction": {
            "cover_equation": (
                "v^2=t^4-6*x_Q*t^2-8*y_Q*t-3*x_Q^2-4*A"
            ),
            "map_to_curve": (
                "x=(t^2-x_Q+v)/2; y=t*(x-x_Q)-y_Q"
            ),
            "exact_round_trips_checked": True,
            "distinct_mod2_class_basis": (
                "the certified basis has full finite-reduction rank mod 2; "
                "therefore its distinct nonempty weight<=2 subsets define "
                "distinct nonzero Q classes modulo 2E(Q)"
            ),
        },
        "declared_budget": {
            "basis_size": len(target.saturated_basis),
            "maximum_cover_subset_weight": MAX_SUBSET_WEIGHT,
            "low_weight_cover_classes_cross_ratio_scored": len(low_weight_plans),
            "all_nonzero_known_mod2_classes_identity_scored": (
                (1 << len(target.saturated_basis)) - 1
            ),
            "full_coset_identity_frontier_retained": len(full_coset_frontier),
            "low_weight_cross_ratio_tranche": LOW_WEIGHT_COVER_COUNT,
            "full_coset_forced_tranche": FULL_COSET_TRANCHE_COUNT,
            "cover_classes_selected": len(selected_plans),
            "cross_ratio_charts_per_cover": CHARTS_PER_COVER,
            "pilot_chart_count": len(selected_plans) * CHARTS_PER_COVER,
            "pilot_height": PILOT_HEIGHT,
            "pilot_timeout_seconds_each": args.pilot_timeout,
            "escalation_chart_count": len(escalated),
            "escalation_height": ESCALATION_HEIGHT,
            "escalation_timeout_seconds_each": args.escalation_timeout,
            "deep_chart_count": min(DEEP_CHART_COUNT, len(escalation_yields)),
            "deep_height": DEEP_HEIGHT,
            "deep_timeout_seconds_each": args.deep_timeout,
            "relation_timeout_seconds": args.relation_timeout,
            "finite_reduction_prime_bound": 500,
            "stack_bytes_each": args.stack_bytes,
            "one_pass_no_retry": True,
        },
        "cover_plans": plan_records,
        "full_coset_identity_frontier": [
            {
                "maximum_known_t_projective_bit_length": score[0],
                "sum_known_t_projective_bit_lengths": score[1],
                "mask": hex(score[2]),
                "subset_indices_one_based": [index + 1 for index in indices],
                "selected_for_search": indices in full_coset_selected_scores,
            }
            for score, indices in full_coset_frontier
        ],
        "runs": run_records,
        "results": {
            "distinct_exact_affine_curve_points": len(discoveries),
            "nonbasis_candidate_points": len(candidates),
            "candidate_point_sha256": point_digest(candidates),
            "exact_relations_in_certified_rank17_subgroup": sum(
                exact for _, exact in proposals
            ),
            "unresolved_by_exact_relation_replay": len(unresolved),
            "combined_finite_reduction_rank": finite_rank,
            "certified_new_directions": certified_gain,
            "certified_rank_lower_bound_after_search": (
                target.certified_rank_lower_bound + certified_gain
            ),
            "candidate_points": candidate_records,
        },
        "interpretation": {
            "exact": (
                "all returned finite cover and elliptic-curve memberships are "
                "exact; displayed subgroup relations are exact group-law replays"
            ),
            "bounded": (
                "the declared finite chart/height boxes are exhausted, but this "
                "is not a complete 2-descent, rank upper bound, or exhaustive "
                "search of all degree-two models"
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
        f"wrote {args.output}: exact_points={len(discoveries)} "
        f"candidates={len(candidates)} unresolved={len(unresolved)} "
        f"finite_rank={finite_rank}",
        flush=True,
    )


if __name__ == "__main__":
    main()
