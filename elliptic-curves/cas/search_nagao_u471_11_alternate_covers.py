#!/usr/bin/env python3
"""Bounded alternate-cover search for the exact Nagao ``u=471/11`` lead.

The input is the independently certified rank-17 subgroup from
``elliptic_nagao_rank17_frontier_certificate.json``.  Its selected subgroup
is cross-checked against the exact height-1,000,000 point subset in the
rank-gain artifact before use.  All ``2^17-1`` nonzero known mod-2 classes are
then scored by their degree-two coordinates.  The best twenty classes receive
three cross-ratio charts and one fixed, strictly bounded search schedule.

Every returned point is mapped to the short curve and checked exactly.  PARI
height-pairing relation proposals are accepted only after exact Fraction
group-law replay; any unresolved points enter a fresh finite-reduction rank
certificate.  This is a bounded point search, not a complete 2-descent or a
rank upper bound.
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
from nagao_skew_height import checkpoint_reference, load_rank17_target
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


Q = Fraction
PARAMETER_U = Q(471, 11)
PARAMETER_T = Q(5579, 22)
UNIFORM_REFERENCE_HEIGHT = 1_000_000
TARGET_RANK = 21
SELECTED_COVER_COUNT = 20
FULL_COSET_RETAIN_COUNT = 20
CERTIFICATE_SHA256 = (
    "7378ce59c72974fe39e0e2a40c740f6a96e8dc555a1361b5aaeef67f4d9e0213"
)
RANK_GAIN_SHA256 = (
    "5f55e0f35368760c65dc2ee66da8edee9ca1153f88e70123e1844d4e2bade559"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_u471_11_alternate_covers.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _unique_certificate_record(data: dict[str, Any]) -> dict[str, Any]:
    matches = [
        record
        for record in data["certificates"]
        if Q(record["parameter_u"]) == PARAMETER_U
    ]
    if len(matches) != 1:
        raise ValueError("the rank-17 certificate has no unique u=471/11 record")
    return matches[0]


def _unique_rank_gain_record(data: dict[str, Any]) -> dict[str, Any]:
    matches = []
    for section_name in ("final_box", "escalation_box"):
        for record in data.get(section_name, {}).get("records", []):
            if (
                Q(record["parameter_u"]) == PARAMETER_U
                and int(record["quartic_naive_height_bound"])
                == UNIFORM_REFERENCE_HEIGHT
            ):
                matches.append(record)
    if len(matches) != 1:
        raise ValueError("the rank-gain artifact has no unique u=471/11 H=1e6 record")
    return matches[0]


def load_exact_inputs(
    certificate_path: Path,
    rank_gain_path: Path,
) -> tuple[Any, Any, str]:
    """Load the certified basis and verify its exact rank-gain lineage."""

    if sha256_file(certificate_path) != CERTIFICATE_SHA256:
        raise ValueError("the pinned rank-17 certificate hash changed")
    if sha256_file(rank_gain_path) != RANK_GAIN_SHA256:
        raise ValueError("the pinned rank-gain artifact hash changed")

    target = load_rank17_target(certificate_path, PARAMETER_U)
    if target.parameter_t != PARAMETER_T:
        raise AssertionError("the derived Nagao base parameter changed")
    if target.certified_rank_lower_bound != 17 or len(target.saturated_basis) != 17:
        raise AssertionError("the input does not contain a certified rank-17 basis")

    certificate_data = json.loads(certificate_path.read_text())
    certificate_record = _unique_certificate_record(certificate_data)
    rank_gain_data = json.loads(rank_gain_path.read_text())
    rank_gain_record = _unique_rank_gain_record(rank_gain_data)
    selected = tuple(
        (Q(record["jacobian_x"]), Q(record["jacobian_y"]))
        for record in rank_gain_record["explicit_numerically_independent_subset"]
    )
    if len(selected) != 17 or any(
        not point_on_short_curve(target.jacobian_coefficients, point)
        for point in selected
    ):
        raise AssertionError("the exact rank-gain subset is incomplete")
    selected_digest = point_digest(selected)
    if selected_digest != certificate_record["input_subset"]["sha256"]:
        raise AssertionError("the certificate does not derive from the pinned subset")
    if certificate_record["input_subset"]["source"] != rank_gain_path.name:
        raise AssertionError("the certificate records a different rank-gain source")

    checkpoint = checkpoint_reference(
        rank_gain_path,
        PARAMETER_U,
        height_bound=UNIFORM_REFERENCE_HEIGHT,
    )
    if checkpoint.stable_pool_numerical_rank != 17:
        raise AssertionError("the pinned uniform checkpoint rank changed")
    if checkpoint.signed_point_count != 100 or checkpoint.unexpected_x_count != 32:
        raise AssertionError("the pinned uniform checkpoint yield changed")
    return target, checkpoint, selected_digest


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts/generated-results"
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
    parser.add_argument("--pilot-timeout", type=float, default=8.0)
    parser.add_argument("--escalation-timeout", type=float, default=12.0)
    parser.add_argument("--deep-timeout", type=float, default=20.0)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_nagao_u471_11_alternate_covers.json",
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

    target, checkpoint, selected_digest = load_exact_inputs(
        args.certificate_input, args.rank_gain_input
    )
    coefficients = target.jacobian_coefficients
    basis = target.saturated_basis
    baseline_signatures = find_mod2_reduction_certificate(
        coefficients, basis, prime_bound=500
    )
    baseline_rank = combined_mod2_rank(baseline_signatures, len(basis))
    if baseline_rank != 17:
        raise AssertionError("the exact input basis no longer has finite rank 17")

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
                "exact_relation_in_certified_rank17_subgroup": exact,
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
            "parameter_u": rational_to_string(PARAMETER_U),
            "parameter_t": rational_to_string(PARAMETER_T),
            "short_weierstrass_coefficients": [
                rational_to_string(value) for value in coefficients
            ],
            "conductor": str(target.conductor),
            "log_conductor": target.log_conductor,
            "root_number": target.root_number,
            "certified_rank_lower_bound_before_search": baseline_rank,
            "target_rank": TARGET_RANK,
        },
        "inputs": {
            "rank17_certificate": {
                "path": str(args.certificate_input),
                "sha256": sha256_file(args.certificate_input),
            },
            "rank_gain_checkpoint": {
                "path": str(args.rank_gain_input),
                "sha256": sha256_file(args.rank_gain_input),
                "height_bound": checkpoint.height_bound,
                "signed_point_count": checkpoint.signed_point_count,
                "unexpected_x_count": checkpoint.unexpected_x_count,
                "unexpected_point_sha256": checkpoint.unexpected_point_sha256,
                "stable_numerical_rank": checkpoint.stable_pool_numerical_rank,
                "exact_selected_subset_sha256": selected_digest,
            },
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
            "exact_relations_in_certified_rank17_subgroup": sum(
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
