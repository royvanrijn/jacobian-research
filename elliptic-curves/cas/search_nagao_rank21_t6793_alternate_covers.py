#!/usr/bin/env python3
"""Bounded alternate-quartic-cover search for the ``T=6793/64`` rank-19 lead.

The input is the exact 19-point finite-reduction certificate embedded in the
unbiased Nagao rank-21-family search artifact.  Every one of the ``2^19-1``
nonzero classes represented by that basis is enumerated in Gray-code order
and scored by the projective bit lengths of its known degree-two coordinates.
The best 32 classes receive three deterministically optimized cross-ratio
charts.  PARI searches fixed pilot, escalation, and deep boxes with a strict
one-pass timeout policy.

Every returned point is mapped back to the short curve and checked using exact
``Fraction`` arithmetic.  Height-pairing relation proposals are accepted only
after exact group-law replay.  Any point not so proved to lie in the certified
rank-19 subgroup is immediately included in an exact finite-reduction
independence calculation.  This remains a bounded point search, not a rank
upper bound or a complete descent.
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
    CoverPlan,
    best_cross_ratio_charts,
    full_coset_identity_frontier,
    point_record,
    projective_height,
    relation_proposals,
    run_chart,
)
from triage_nagao_rank13_finalists import point_digest, point_on_short_curve


Q = Fraction
PARAMETER_T = Q(6793, 64)
INPUT_CANDIDATE_ID = "unbiased-6793-64"
INPUT_CERTIFIED_RANK = 19
TARGET_RANK = 21
SELECTED_COVER_COUNT = 32
FULL_COSET_RETAIN_COUNT = 32
CHARTS_PER_COVER = 3
PILOT_HEIGHT = 50_000
ESCALATION_HEIGHT = 250_000
ESCALATION_CHART_COUNT = 12
DEEP_HEIGHT = 1_000_000
DEEP_CHART_COUNT = 4
CERTIFICATE_SHA256 = "0df1ea176d66c0446de16774ac5129253ade2379d6b6af9a591c444c8822bb6e"
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank21_t6793_alternate_covers.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_exact_basis(
    certificate_path: Path,
) -> tuple[
    tuple[Fraction, ...],
    tuple[tuple[Fraction, Fraction], ...],
    dict[str, Any],
    dict[str, Any],
]:
    """Load and replay the pinned exact rank-19 certificate."""

    if sha256_file(certificate_path) != CERTIFICATE_SHA256:
        raise ValueError("the pinned rank-19 certificate hash changed")
    data = json.loads(certificate_path.read_text())
    certificate = data["exact_rank_certificate"]
    conductor = data["candidate"]
    if (
        certificate["certified_algebraic_rank_lower_bound"]
        != INPUT_CERTIFIED_RANK
        or certificate["combined_exact_rank_over_F2"] != INPUT_CERTIFIED_RANK
    ):
        raise AssertionError("the input no longer certifies rank at least 19")
    if Q(conductor["parameter_t"]) != PARAMETER_T:
        raise AssertionError("the pinned certificate parameter changed")
    if conductor["root_number"] != -1 or not conductor[
        "below_strict_log_conductor_target"
    ]:
        raise AssertionError("the pinned conductor replay changed")

    coefficients = short_jacobian_coefficients(RANK21_CONSTRUCTION, PARAMETER_T)
    if tuple(Q(value) for value in conductor["short_weierstrass_coefficients"]) != (
        coefficients
    ):
        raise AssertionError("the pinned short model changed")
    basis = tuple(
        (Q(record["jacobian_x"]), Q(record["jacobian_y"]))
        for record in certificate["saturated_basis"]
    )
    if len(basis) != INPUT_CERTIFIED_RANK:
        raise AssertionError("the pinned saturated basis has the wrong size")
    if any(not point_on_short_curve(coefficients, point) for point in basis):
        raise AssertionError("a pinned saturated-basis point left the short curve")
    return coefficients, basis, certificate, conductor


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts/generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--certificate-input",
        type=Path,
        default=(
            generated / "elliptic_nagao_rank21_t6793_rank19_certificate.json"
        ),
    )
    parser.add_argument("--pilot-timeout", type=float, default=8.0)
    parser.add_argument("--escalation-timeout", type=float, default=12.0)
    parser.add_argument("--deep-timeout", type=float, default=20.0)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            generated / "elliptic_nagao_rank21_t6793_alternate_covers.json"
        ),
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

    coefficients, basis, certificate, conductor = load_exact_basis(
        args.certificate_input
    )
    baseline_signatures = find_mod2_reduction_certificate(
        coefficients, basis, prime_bound=500
    )
    baseline_rank = combined_mod2_rank(baseline_signatures, len(basis))
    if baseline_rank != INPUT_CERTIFIED_RANK:
        raise AssertionError("the pinned basis did not replay at finite rank 19")

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
                best_cross_ratio_charts(
                    cover, basis, count=CHARTS_PER_COVER
                ),
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
                raise AssertionError("an absorbed point left the short curve")
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
    deepened = escalation_yields[:DEEP_CHART_COUNT]
    for _, plan, chart in deepened:
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
    if len(proposals) != len(candidates):
        raise AssertionError("not every candidate received a membership replay")
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
                "exact_relation_in_certified_rank19_subgroup": exact,
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
                        "maximum_log10_known_projective_height": (
                            chart.maximum_log_height
                        ),
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
            "candidate_id": INPUT_CANDIDATE_ID,
            "parameter_t": rational_to_string(PARAMETER_T),
            "minimal_model": conductor["minimal_model"],
            "conductor": conductor["conductor"],
            "log_conductor": conductor["log_conductor"],
            "root_number": conductor["root_number"],
            "certified_rank_lower_bound_before_search": baseline_rank,
            "target_rank": TARGET_RANK,
        },
        "input": {
            "path": str(args.certificate_input),
            "sha256": sha256_file(args.certificate_input),
            "embedded_saturated_basis_sha256": certificate[
                "saturated_basis_sha256"
            ],
        },
        "construction": {
            "cover_equation": "v^2=t^4-6*x_Q*t^2-8*y_Q*t-3*x_Q^2-4*A",
            "map_to_curve": "x=(t^2-x_Q+v)/2; y=t*(x-x_Q)-y_Q",
            "exact_curve_membership_checked_for_every_returned_point": True,
            "mod2_class_scope": (
                "the finite-reduction-independent rank-19 basis represents "
                "2^19 distinct classes; every nonzero represented class was "
                "identity-coordinate scored"
            ),
        },
        "declared_budget": {
            "basis_size": len(basis),
            "all_nonzero_certified_mod2_classes_identity_scored": (
                (1 << len(basis)) - 1
            ),
            "identity_frontier_retained": len(full_frontier),
            "cover_classes_selected": len(plans),
            "cross_ratio_charts_per_cover": CHARTS_PER_COVER,
            "pilot_chart_count": len(plans) * CHARTS_PER_COVER,
            "pilot_height": PILOT_HEIGHT,
            "pilot_timeout_seconds_each": args.pilot_timeout,
            "escalation_chart_count": len(escalated),
            "escalation_height": ESCALATION_HEIGHT,
            "escalation_timeout_seconds_each": args.escalation_timeout,
            "deep_chart_count": len(deepened),
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
            "exact_relations_in_certified_rank19_subgroup": sum(
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
