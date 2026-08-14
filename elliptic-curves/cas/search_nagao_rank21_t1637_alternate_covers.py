#!/usr/bin/env python3
"""Bounded all-class alternate-cover search for Nagao ``T=1637/12``.

The exact finite-reduction-certified rank-16 basis is loaded from the pinned
standalone certificate.  Every one of its ``2^16-1`` nonzero known mod-2
classes is scored by the degree-two coordinates of the other basis points.
The best twenty classes receive three optimized cross-ratio charts and the
fixed 50k/250k/1m staged search used by the preceding frontier passes.

All returned points are mapped and checked with exact rational arithmetic.
Height-pairing relation proposals count only after exact Fraction group-law
replay.  Any unresolved points enter a new finite-reduction independence
certificate.  The declared boxes are bounded and are not a complete descent.
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
from typing import Any, Iterable

from alternate_quartic_covers import alternate_cover, short_subset_sum
from ek_k3 import rational_to_string
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)
from nagao_1994 import RANK21_CONSTRUCTION, short_jacobian_coefficients
from pari_bridge import pari_version
from search_nagao_u135_alternate_covers import (
    CHARTS_PER_COVER,
    DEEP_CHART_COUNT,
    DEEP_HEIGHT,
    ESCALATION_CHART_COUNT,
    ESCALATION_HEIGHT,
    PILOT_HEIGHT,
    CoverPlan,
    best_cross_ratio_charts,
    full_coset_identity_frontier,
    point_record,
    projective_height,
    relation_proposals,
    run_chart,
)
from triage_nagao_rank13_finalists import point_digest, point_on_short_curve
from verify_nagao_rank21_t1637 import PARAMETER_T


Q = Fraction
TARGET_RANK = 21
PARITY_COMPATIBLE_TARGET_RANK = 22
SELECTED_COVER_COUNT = 20
FULL_COSET_RETAIN_COUNT = 20
CERTIFICATE_SHA256 = (
    "0f71ad6f41e215bcf92cf0c19d6985d1d7ce64ea9d989388fdf600d138d9813e"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank21_t1637_alternate_covers.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_exact_basis(
    certificate_path: Path,
) -> tuple[tuple[Fraction, ...], tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    if sha256_file(certificate_path) != CERTIFICATE_SHA256:
        raise ValueError("the pinned rank-16 certificate hash changed")
    data = json.loads(certificate_path.read_text())
    if Q(data["candidate"]["parameter_t"]) != PARAMETER_T:
        raise AssertionError("the certificate parameter changed")
    coefficients = short_jacobian_coefficients(RANK21_CONSTRUCTION, PARAMETER_T)
    stored_coefficients = tuple(
        Q(value) for value in data["candidate"]["short_weierstrass_coefficients"]
    )
    if stored_coefficients != coefficients:
        raise AssertionError("the certificate short model changed")
    basis = tuple(
        (Q(record["jacobian_x"]), Q(record["jacobian_y"]))
        for record in data["exact_rank_certificate"]["saturated_basis"]
    )
    if len(basis) != 16 or any(
        not point_on_short_curve(coefficients, point) for point in basis
    ):
        raise AssertionError("the exact rank-16 basis is incomplete")
    if data["exact_rank_certificate"]["certified_algebraic_rank_lower_bound"] != 16:
        raise AssertionError("the input lacks its exact rank-16 claim")
    return coefficients, basis, data


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts/generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--certificate-input",
        type=Path,
        default=generated / "elliptic_nagao_rank21_t1637_rank16_certificate.json",
    )
    parser.add_argument("--pilot-timeout", type=float, default=8.0)
    parser.add_argument("--escalation-timeout", type=float, default=12.0)
    parser.add_argument("--deep-timeout", type=float, default=20.0)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_nagao_rank21_t1637_alternate_covers.json",
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
        if not 0 < getattr(args, name) <= 60:
            raise SystemExit(f"--{name.replace('_', '-')} must be in (0,60]")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes is too small")

    coefficients, basis, certificate = load_exact_basis(args.certificate_input)
    baseline_signatures = find_mod2_reduction_certificate(
        coefficients, basis, prime_bound=500
    )
    baseline_rank = combined_mod2_rank(baseline_signatures, len(basis))
    if baseline_rank != 16:
        raise AssertionError("the exact input basis no longer has finite rank 16")

    full_frontier = full_coset_identity_frontier(
        coefficients, basis, retain_count=FULL_COSET_RETAIN_COUNT
    )
    plans = []
    for _, subset_indices in full_frontier[:SELECTED_COVER_COUNT]:
        base_point = short_subset_sum(coefficients, basis, subset_indices)
        if base_point is None:
            raise AssertionError("a selected nonzero cover class vanished")
        cover = alternate_cover(coefficients, base_point)
        plans.append(
            CoverPlan(
                subset_indices,
                cover,
                best_cross_ratio_charts(cover, basis, count=CHARTS_PER_COVER),
            )
        )

    run_records: list[dict[str, Any]] = []
    discoveries: dict[tuple[Fraction, Fraction], set[str]] = {}
    chart_yields = []

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
            if not point_on_short_curve(coefficients, point):
                raise AssertionError("a mapped alternate-cover point missed the curve")
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
                projective_height(point[0]),
                projective_height(point[1]),
                point,
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
        coefficients, basis + unresolved, prime_bound=500
    )
    augmented_rank = combined_mod2_rank(
        augmented_signatures, len(basis) + len(unresolved)
    )
    certified_gain = max(0, augmented_rank - baseline_rank)
    certified_rank = baseline_rank + certified_gain

    candidate_records = []
    for point, (relation, exact) in zip(candidates, proposals):
        record: dict[str, Any] = point_record(point)
        record.update(
            {
                "sources": sorted(discoveries[point]),
                "exact_relation_in_certified_rank16_subgroup": exact,
                "basis_relation": list(relation) if exact else None,
                "fraction_group_law_replay": exact,
            }
        )
        candidate_records.append(record)

    score_by_indices = {indices: score for score, indices in full_frontier}
    cover_records = []
    for plan in plans:
        score = score_by_indices[plan.subset_indices]
        cover_records.append(
            {
                "id": plan.identifier,
                "mask": hex(score[2]),
                "identity_score_maximum_known_t_projective_bit_length": score[0],
                "identity_score_sum_known_t_projective_bit_lengths": score[1],
                "subset_indices_one_based": [
                    index + 1 for index in plan.subset_indices
                ],
                "base_point": point_record(plan.cover.base_point),
                "quartic_coefficients_ascending": [
                    rational_to_string(value) for value in plan.cover.coefficients
                ],
                "cross_ratio_charts": [
                    {
                        "normalizing_basis_indices_one_based": [
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
    artifact = {
        "schema_version": 1,
        "status": "bounded_alternate_cover_search_complete",
        "candidate": {
            "parameter_t": rational_to_string(PARAMETER_T),
            "minimal_model": certificate["candidate"]["minimal_model"],
            "conductor": certificate["candidate"]["conductor"],
            "log_conductor": certificate["candidate"]["log_conductor"],
            "root_number": certificate["candidate"]["root_number"],
            "certified_rank_lower_bound_before_search": baseline_rank,
            "target_rank": TARGET_RANK,
            "parity_heuristic_next_compatible_target_rank": (
                PARITY_COMPATIBLE_TARGET_RANK
            ),
        },
        "input": {
            "path": str(args.certificate_input),
            "sha256": sha256_file(args.certificate_input),
        },
        "declared_budget": {
            "all_nonzero_certified_mod2_classes_identity_scored": (
                (1 << len(basis)) - 1
            ),
            "cover_classes_selected": len(plans),
            "cross_ratio_charts_per_cover": CHARTS_PER_COVER,
            "pilot_chart_count": len(plans) * CHARTS_PER_COVER,
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
        "cover_plans": cover_records,
        "runs": run_records,
        "results": {
            "distinct_exact_affine_curve_points": len(discoveries),
            "nonbasis_candidate_points": len(candidates),
            "candidate_point_sha256": point_digest(candidates),
            "exact_relations_in_certified_rank16_subgroup": sum(
                exact for _, exact in proposals
            ),
            "unresolved_by_exact_relation_replay": len(unresolved),
            "augmented_finite_reduction_primes": [
                signature.prime for signature in augmented_signatures
            ],
            "combined_finite_reduction_rank": augmented_rank,
            "certified_new_directions": certified_gain,
            "certified_rank_lower_bound_after_search": certified_rank,
            "target_rank_21_achieved": certified_rank >= TARGET_RANK,
            "parity_heuristic_rank_22_reached": (
                certified_rank >= PARITY_COMPATIBLE_TARGET_RANK
            ),
            "candidate_points": candidate_records,
        },
        "interpretation": {
            "exact": (
                "all returned memberships and displayed subgroup relations are "
                "exact; any rank gain is certified by finite reductions"
            ),
            "bounded": (
                "the declared chart boxes are exhausted; this is not a rank "
                "upper bound or a complete search of all degree-two models"
            ),
            "parity": (
                "root number +1 is only heuristic evidence for even analytic "
                "rank and is not used in the algebraic certificate"
            ),
        },
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "platform": platform.platform(),
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "actual_command": " ".join(
            shlex.quote(part) for part in [sys.executable, *sys.argv]
        ),
        "script_sha256": sha256_file(script_path),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        f"wrote {args.output}: exact_points={len(discoveries)} "
        f"candidates={len(candidates)} unresolved={len(unresolved)} "
        f"certified_rank={certified_rank}",
        flush=True,
    )


if __name__ == "__main__":
    main()
