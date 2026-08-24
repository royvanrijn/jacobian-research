#!/usr/bin/env python3
"""Bounded all-class alternate-cover pass on Mestre T=2731/36.

The Delta=11/5 explicit formula narrowly fails to close this exact rank-15
fiber below analytic rank 17.  This fixed-fiber follow-up therefore scores all
``2^15-1`` nonzero classes of the certified saturated basis, combines the
best low-weight and full-coset classes, and searches three optimized
cross-ratio charts per selected cover in fixed H=50k/250k/1m tiers.

Every returned point is mapped and checked exactly.  Height-pairing relation
proposals are replayed with exact Fraction group arithmetic, and any
unresolved directions receive a fresh exact finite-reduction rank test.  This
is a bounded cover search, not a complete descent or rank upper bound.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shlex
import sys
from typing import Any, Iterable

from alternate_quartic_covers import alternate_cover, short_subset_sum
from ek_k3 import rational_to_string
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from pari_bridge import pari_version
from search_mestre_root_tuple_scale import point_digest, point_on_short_curve
from search_mestre_root_tuple_scale_max100 import stable_json_digest
from search_nagao_u135_alternate_covers import (
    CoverPlan,
    best_cross_ratio_charts,
    enumerate_cover_plans,
    full_coset_identity_frontier,
    point_record,
    projective_height,
    relation_proposals,
    run_chart,
)


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

ROOTS = (0, 7, 93, 154, 161, 191)
PARAMETER = Q(2731, 36)
CERTIFICATE_FILENAME = "elliptic_mestre_rank15_2731_36.json"
EXPECTED_CERTIFICATE_SHA256 = (
    "5f91987e9fd21887afbe0cd376e7b56844a37e0e70ade6fc713aaa3121e87c1a"
)
EXPECTED_CERTIFICATE_RESULT_SHA256 = (
    "4d76dc4adf65376027135a1c6023c4cf9040e5d19aa2f8b84a1b79c7b8877252"
)
EXPLICIT_FORMULA_FILENAME = (
    "elliptic_mestre_rank15_2731_36_explicit_formula_delta22.json"
)
MAX_SUBSET_WEIGHT = 2
LOW_WEIGHT_COVER_COUNT = 12
FULL_COSET_RETAIN_COUNT = 20
FULL_COSET_TRANCHE_COUNT = 8
SELECTED_COVER_COUNT = LOW_WEIGHT_COVER_COUNT + FULL_COSET_TRANCHE_COUNT
CHARTS_PER_COVER = 3
PILOT_HEIGHT = 50_000
PILOT_TIMEOUT = 8.0
ESCALATION_HEIGHT = 250_000
ESCALATION_CHART_COUNT = 8
ESCALATION_TIMEOUT = 12.0
DEEP_HEIGHT = 1_000_000
DEEP_CHART_COUNT = 2
DEEP_TIMEOUT = 20.0
RELATION_TIMEOUT = 60.0
STACK_BYTES = 512_000_000
FINITE_PRIME_BOUND = 500


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exclusive_write(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts/generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--certificate", type=Path, default=generated / CERTIFICATE_FILENAME
    )
    parser.add_argument(
        "--explicit-formula", type=Path,
        default=generated / EXPLICIT_FORMULA_FILENAME,
    )
    parser.add_argument("--pilot-timeout", type=float, default=PILOT_TIMEOUT)
    parser.add_argument("--escalation-timeout", type=float, default=ESCALATION_TIMEOUT)
    parser.add_argument("--deep-timeout", type=float, default=DEEP_TIMEOUT)
    parser.add_argument("--relation-timeout", type=float, default=RELATION_TIMEOUT)
    parser.add_argument("--stack-bytes", type=int, default=STACK_BYTES)
    parser.add_argument(
        "--output", type=Path,
        default=generated / "elliptic_mestre_rank15_2731_36_alternate_covers.json",
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if (
        args.pilot_timeout != PILOT_TIMEOUT
        or args.escalation_timeout != ESCALATION_TIMEOUT
        or args.deep_timeout != DEEP_TIMEOUT
        or args.relation_timeout != RELATION_TIMEOUT
        or args.stack_bytes != STACK_BYTES
    ):
        raise SystemExit("the fixed-cover resource bounds are pinned")


def main() -> None:
    args = build_parser().parse_args()
    validate_args(args)
    if args.output.exists():
        raise SystemExit("refusing to overwrite the fixed-cover artifact")
    if file_sha256(args.certificate) != EXPECTED_CERTIFICATE_SHA256:
        raise AssertionError("the exact rank-15 certificate changed")
    certificate = json.loads(args.certificate.read_text())
    if certificate["result_sha256"] != EXPECTED_CERTIFICATE_RESULT_SHA256:
        raise AssertionError("the rank-15 certificate result digest changed")
    formula = json.loads(args.explicit_formula.read_text())
    if formula["explicit_formula"]["strictly_below_17"]:
        raise AssertionError("the fixed fiber is already conditionally closed")
    if tuple(certificate["curve"]["roots"]) != ROOTS:
        raise AssertionError("the certified family changed")
    if Q(certificate["curve"]["parameter"]) != PARAMETER:
        raise AssertionError("the certified parameter changed")

    coefficients = tuple(Q(value) for value in certificate["curve"]["weierstrass_coefficients"])
    basis = tuple(
        (Q(record["x"]), Q(record["y"]))
        for record in certificate["saturated_basis"]
    )
    if len(basis) != 15 or any(
        not point_on_short_curve(coefficients, point) for point in basis
    ):
        raise AssertionError("the exact saturated basis changed")
    basis_signatures = find_mod2_reduction_certificate(
        coefficients, basis, prime_bound=FINITE_PRIME_BOUND
    )
    exact_basis_rank = combined_mod2_rank(basis_signatures, len(basis))
    if exact_basis_rank != 15:
        raise AssertionError("mod-2 reductions did not replay the rank-15 basis")
    two_torsion_prime = find_two_torsion_certificate_prime(coefficients)

    low_weight_plans = enumerate_cover_plans(
        coefficients,
        basis,
        maximum_subset_weight=MAX_SUBSET_WEIGHT,
        charts_per_cover=CHARTS_PER_COVER,
    )
    expected_low_weight = sum(
        math.comb(len(basis), weight)
        for weight in range(1, MAX_SUBSET_WEIGHT + 1)
    )
    if len(low_weight_plans) != expected_low_weight:
        raise AssertionError("the low-weight class count changed")
    full_frontier = full_coset_identity_frontier(
        coefficients, basis, retain_count=FULL_COSET_RETAIN_COUNT
    )
    selected_by_indices: dict[tuple[int, ...], CoverPlan] = {
        plan.subset_indices: plan
        for plan in low_weight_plans[:LOW_WEIGHT_COVER_COUNT]
    }
    full_scores: dict[tuple[int, ...], tuple[int, int, int]] = {}
    for score, indices in full_frontier:
        if indices in selected_by_indices:
            continue
        base_point = short_subset_sum(coefficients, basis, indices)
        if base_point is None:
            raise AssertionError("a retained nonzero class vanished")
        cover = alternate_cover(coefficients, base_point)
        selected_by_indices[indices] = CoverPlan(
            indices,
            cover,
            best_cross_ratio_charts(cover, basis, count=CHARTS_PER_COVER),
        )
        full_scores[indices] = score
        if len(full_scores) == FULL_COSET_TRANCHE_COUNT:
            break
    plans = tuple(selected_by_indices.values())
    if len(plans) != SELECTED_COVER_COUNT:
        raise AssertionError("the selected cover union changed size")

    run_records: list[dict[str, Any]] = []
    discoveries: dict[tuple[Fraction, Fraction], set[str]] = {}
    chart_yields: list[tuple[int, CoverPlan, Any]] = []

    def absorb(
        plan: CoverPlan,
        chart: Any,
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

    for plan in plans:
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
        for basis_point in basis
        for point in (basis_point, (basis_point[0], -basis_point[1]))
    }
    candidates = tuple(
        sorted(
            (point for point in discoveries if point not in basis_with_signs),
            key=lambda point: (
                projective_height(point[0]), projective_height(point[1]), point
            ),
        )
    )
    proposals = relation_proposals(
        coefficients,
        basis,
        candidates,
        timeout=args.relation_timeout,
        stack_bytes=args.stack_bytes,
    )
    unresolved = tuple(
        point for point, (_, exact) in zip(candidates, proposals) if not exact
    )
    augmented_signatures = find_mod2_reduction_certificate(
        coefficients, basis + unresolved, prime_bound=FINITE_PRIME_BOUND
    )
    augmented_rank = combined_mod2_rank(
        augmented_signatures, len(basis) + len(unresolved)
    )
    certified_gain = max(0, augmented_rank - exact_basis_rank)
    certified_rank = exact_basis_rank + certified_gain

    candidate_records = []
    for point, (relation, exact) in zip(candidates, proposals):
        row: dict[str, Any] = point_record(point)
        row.update(
            {
                "sources": sorted(discoveries[point]),
                "exact_relation_in_certified_rank15_subgroup": exact,
                "basis_relation": list(relation) if exact else None,
                "fraction_group_law_replay": exact,
            }
        )
        candidate_records.append(row)

    score_by_indices = {indices: score for score, indices in full_frontier}
    cover_records = []
    for plan in plans:
        score = score_by_indices.get(plan.subset_indices)
        cover_records.append(
            {
                "id": plan.identifier,
                "subset_indices_one_based": [
                    index + 1 for index in plan.subset_indices
                ],
                "selection_strata": [
                    label
                    for label, condition in (
                        ("low-weight-cross-ratio", plan in low_weight_plans[:LOW_WEIGHT_COVER_COUNT]),
                        ("full-coset-identity", plan.subset_indices in full_scores),
                    )
                    if condition
                ],
                "full_coset_identity_score": list(score) if score is not None else None,
                "base_point": point_record(plan.cover.base_point),
                "quartic_coefficients_ascending": [
                    rational_to_string(value) for value in plan.cover.coefficients
                ],
                "cross_ratio_charts": [
                    {
                        "basis_indices_one_based": [
                            index + 1 for index in chart.basis_indices
                        ],
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
    root = script_path.parents[2]
    helper_path = script_path.with_name("search_nagao_u135_alternate_covers.py")
    cover_engine = script_path.with_name("alternate_quartic_covers.py")
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "bounded fixed-fiber alternate-cover search complete",
        "curve": {
            "roots": list(ROOTS),
            "parameter": rational_to_string(PARAMETER),
            "short_weierstrass_coefficients": [
                rational_to_string(value) for value in coefficients
            ],
            "conductor": certificate["curve"]["conductor"],
            "log_conductor": certificate["curve"]["log_conductor"],
            "root_number": certificate["curve"]["root_number"],
        },
        "exact_rank15_baseline": {
            "saturated_basis_sha256": point_digest(basis),
            "two_torsion_certificate_prime": two_torsion_prime,
            "finite_reduction_primes": [
                signature.prime for signature in basis_signatures
            ],
            "combined_exact_rank_over_F2": exact_basis_rank,
            "certified_algebraic_rank_lower_bound": exact_basis_rank,
        },
        "explicit_formula_gate": {
            "artifact": args.explicit_formula.name,
            "artifact_sha256": file_sha256(args.explicit_formula),
            "conservative_upper": formula["explicit_formula"][
                "conservative_explicit_formula_upper"
            ],
            "strictly_below_17": False,
            "fixed_cover_search_required": True,
        },
        "declared_budget": {
            "all_nonzero_certified_mod2_classes_identity_scored": (1 << len(basis)) - 1,
            "low_weight_classes_built": len(low_weight_plans),
            "low_weight_classes_selected": LOW_WEIGHT_COVER_COUNT,
            "full_coset_frontier_retained": FULL_COSET_RETAIN_COUNT,
            "full_coset_classes_selected": FULL_COSET_TRANCHE_COUNT,
            "cover_classes_selected": len(plans),
            "charts_per_cover": CHARTS_PER_COVER,
            "pilot_chart_count": len(plans) * CHARTS_PER_COVER,
            "pilot_height": PILOT_HEIGHT,
            "pilot_timeout_seconds_each": PILOT_TIMEOUT,
            "escalation_chart_count": len(escalated),
            "escalation_height": ESCALATION_HEIGHT,
            "escalation_timeout_seconds_each": ESCALATION_TIMEOUT,
            "deep_chart_count": min(DEEP_CHART_COUNT, len(escalation_yields)),
            "deep_height": DEEP_HEIGHT,
            "deep_timeout_seconds_each": DEEP_TIMEOUT,
            "relation_timeout_seconds": RELATION_TIMEOUT,
            "finite_reduction_prime_bound": FINITE_PRIME_BOUND,
            "one_pass_no_retry": True,
        },
        "cover_plans": cover_records,
        "runs": run_records,
        "results": {
            "distinct_exact_affine_curve_points": len(discoveries),
            "nonbasis_candidate_points": len(candidates),
            "candidate_point_sha256": point_digest(candidates),
            "exact_relations_in_certified_rank15_subgroup": sum(
                exact for _, exact in proposals
            ),
            "unresolved_by_exact_relation_replay": len(unresolved),
            "combined_finite_reduction_rank": augmented_rank,
            "certified_new_directions": certified_gain,
            "certified_rank_lower_bound_after_search": certified_rank,
            "rank16_signal": certified_rank >= 16,
            "rank21_target_achieved": certified_rank >= 21,
            "candidate_points": candidate_records,
        },
        "interpretation": {
            "exact": (
                "the baseline independence, returned point membership, exact "
                "subgroup relations, and final finite-reduction rank are exact"
            ),
            "bounded": (
                "the declared chart boxes are exhausted; this is not a complete "
                "2-descent or rank upper bound"
            ),
        },
        "provenance": {
            "script": str(script_path.relative_to(root)),
            "script_sha256": file_sha256(script_path),
            "certificate": str(args.certificate.relative_to(root)),
            "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
            "cover_helper": str(helper_path.relative_to(root)),
            "cover_helper_sha256": file_sha256(helper_path),
            "cover_engine": str(cover_engine.relative_to(root)),
            "cover_engine_sha256": file_sha256(cover_engine),
            "command": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
            "temporary_files_removed": True,
            "no_retry": True,
        },
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "platform": platform.platform(),
        },
    }
    artifact["result_sha256"] = stable_json_digest(
        {
            "curve": artifact["curve"],
            "baseline": artifact["exact_rank15_baseline"],
            "formula": artifact["explicit_formula_gate"],
            "budget": artifact["declared_budget"],
            "plans": artifact["cover_plans"],
            "runs": artifact["runs"],
            "results": artifact["results"],
        }
    )
    exclusive_write(args.output, artifact)
    print(
        f"complete exact_points={len(discoveries)} candidates={len(candidates)} "
        f"unresolved={len(unresolved)} certified_rank={certified_rank} "
        f"output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
