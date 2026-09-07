#!/usr/bin/env python3
"""One bounded alternate-cover rank-gain pass for Nagao ``u=75/2``.

An exact bivariate sieve supplied four quartic points.  Together with the
thirteen affine sections and the split point at infinity they form an
18-point pool whose stable numerical height rank is 15.  A small-prime PARI
saturation is used only to propose a better basis; exact finite reductions
then prove that the returned 15 points are independent.

All ``2^15-1`` nonzero classes represented by that exact basis are Gray-code
enumerated and scored by the projective bit heights of the other known points
in the associated degree-two coordinate.  The best twenty classes receive
three optimized cross-ratio charts and a strictly bounded staged point search.
Every result is checked exactly, replayed against the rank-15 subgroup when
possible, and included in a fresh finite-reduction rank test otherwise.
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
from extend_nagao_u42_frontier import saturate_exact_basis
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK13_CONSTRUCTION,
    primitive_quartic_coefficients,
    quartic_point_to_short_jacobian,
    quartic_value,
    rank13_base_changed_short_jacobian_coefficients,
    rank13_base_parameter,
    rank13_known_quartic_points,
)
from pari_bridge import minimal_curve_data, pari_version
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
from triage_nagao_rank13_finalists import (
    height_matrix_replay,
    point_digest,
    point_on_short_curve,
    split_infinity_jacobian_point,
    stable_height_rank,
)


Q = Fraction
PARAMETER_U = Q(75, 2)
PARAMETER_T = Q(1181, 4)
TARGET_RANK = 21
SIEVE_ROWS = (
    # (h, k, W), with x=T+k/h and z=W/(2*a^2*b^2*h^2), a/b=75/2.
    (2, -4049, 573371851297500),
    (2, -1651, 54597770955000),
    (6, -2755, 60839174182500),
    (8, -3703, 107230492946250),
)
EXPECTED_HEIGHT_SUBSET = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 16, 17)
FULL_COSET_RETAIN_COUNT = 20
SELECTED_COVER_COUNT = 20
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_u75_alternate_covers.py"
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_seed_pool() -> tuple[
    tuple[Fraction, ...],
    tuple[tuple[Fraction, Fraction], ...],
    tuple[dict[str, Any], ...],
]:
    """Construct the pinned eighteen exact points and their provenance."""

    if rank13_base_parameter(PARAMETER_U) != PARAMETER_T:
        raise AssertionError("the pinned base-change parameters disagree")
    coefficients = rank13_base_changed_short_jacobian_coefficients(PARAMETER_U)
    affine_quartic = rank13_known_quartic_points(PARAMETER_T)
    affine_images = tuple(
        quartic_point_to_short_jacobian(
            RANK13_CONSTRUCTION, PARAMETER_T, quartic_point
        )
        for quartic_point in affine_quartic
    )
    infinity = split_infinity_jacobian_point(PARAMETER_U)
    quartic_coefficients = primitive_quartic_coefficients(
        RANK13_CONSTRUCTION, PARAMETER_T
    )
    extras = []
    records = []
    for h_value, k_value, w_value in SIEVE_ROWS:
        x_value = PARAMETER_T + Q(k_value, h_value)
        z_value = Q(
            w_value,
            2
            * PARAMETER_U.numerator**2
            * PARAMETER_U.denominator**2
            * h_value**2,
        )
        if z_value**2 != quartic_value(quartic_coefficients, x_value):
            raise AssertionError("a pinned sieve row missed the exact quartic")
        image = quartic_point_to_short_jacobian(
            RANK13_CONSTRUCTION, PARAMETER_T, (x_value, z_value)
        )
        if not point_on_short_curve(coefficients, image):
            raise AssertionError("a pinned sieve image missed the Jacobian")
        extras.append(image)
        records.append(
            {
                "h": h_value,
                "k": k_value,
                "W": str(w_value),
                "quartic_x": rational_to_string(x_value),
                "quartic_z": rational_to_string(z_value),
                "jacobian_x": rational_to_string(image[0]),
                "jacobian_y": rational_to_string(image[1]),
                "exact_quartic_and_jacobian_membership_checked": True,
            }
        )
    pool = affine_images + (infinity,) + tuple(extras)
    if len(pool) != 18 or any(
        not point_on_short_curve(coefficients, point) for point in pool
    ):
        raise AssertionError("the exact seed pool is incomplete")
    return coefficients, pool, tuple(records)


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height-timeout", type=float, default=20.0)
    parser.add_argument("--saturation-timeout", type=float, default=20.0)
    parser.add_argument("--pilot-timeout", type=float, default=8.0)
    parser.add_argument("--escalation-timeout", type=float, default=12.0)
    parser.add_argument("--deep-timeout", type=float, default=20.0)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--conductor-timeout", type=float, default=10.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts/generated-results/elliptic_nagao_u75_alternate_covers.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    for name in (
        "height_timeout",
        "saturation_timeout",
        "pilot_timeout",
        "escalation_timeout",
        "deep_timeout",
        "relation_timeout",
        "conductor_timeout",
    ):
        if not 0 < getattr(args, name) <= 60:
            raise SystemExit(f"--{name.replace('_', '-')} must be in (0,60]")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes is too small")

    coefficients, pool, sieve_records = exact_seed_pool()
    height_runs = height_matrix_replay(
        coefficients,
        pool,
        precisions=(72, 120),
        timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    if stable_height_rank(height_runs) != 15:
        raise AssertionError("the seed-pool stable numerical rank changed")
    indices = tuple(height_runs[-1]["subset_indices_one_based"])
    if indices != EXPECTED_HEIGHT_SUBSET:
        raise AssertionError("the stable height subset changed")
    numerical_basis = tuple(pool[index - 1] for index in indices)
    saturated_basis, saturation = saturate_exact_basis(
        coefficients,
        numerical_basis,
        prime_bound=20,
        timeout=args.saturation_timeout,
        stack_bytes=args.stack_bytes,
    )
    if len(saturated_basis) != 15:
        raise AssertionError("small-prime saturation changed the basis length")
    basis_signatures = find_mod2_reduction_certificate(
        coefficients, saturated_basis, prime_bound=500
    )
    exact_basis_rank = combined_mod2_rank(basis_signatures, len(saturated_basis))
    if exact_basis_rank != 15:
        raise AssertionError("finite reductions did not certify the rank-15 basis")
    two_torsion_prime = find_two_torsion_certificate_prime(coefficients)

    full_frontier = full_coset_identity_frontier(
        coefficients, saturated_basis, retain_count=FULL_COSET_RETAIN_COUNT
    )
    plans = []
    for score, subset_indices in full_frontier[:SELECTED_COVER_COUNT]:
        base_point = short_subset_sum(
            coefficients, saturated_basis, subset_indices
        )
        if base_point is None:
            raise AssertionError("a selected exact cover class vanished")
        cover = alternate_cover(coefficients, base_point)
        plans.append(
            CoverPlan(
                subset_indices,
                cover,
                best_cross_ratio_charts(
                    cover, saturated_basis, count=CHARTS_PER_COVER
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
        for basis_point in saturated_basis
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
        saturated_basis,
        candidates,
        timeout=args.relation_timeout,
        stack_bytes=args.stack_bytes,
    )
    unresolved = tuple(
        point for point, (_, exact) in zip(candidates, proposals) if not exact
    )
    augmented_signatures = find_mod2_reduction_certificate(
        coefficients, saturated_basis + unresolved, prime_bound=500
    )
    augmented_rank = combined_mod2_rank(
        augmented_signatures, len(saturated_basis) + len(unresolved)
    )
    certified_gain = max(0, augmented_rank - exact_basis_rank)
    certified_rank = exact_basis_rank + certified_gain
    conductor = minimal_curve_data(
        coefficients, timeout=args.conductor_timeout, stack_bytes=args.stack_bytes
    )

    candidate_records = []
    for point, (relation, exact) in zip(candidates, proposals):
        record: dict[str, Any] = point_record(point)
        record.update(
            {
                "sources": sorted(discoveries[point]),
                "exact_relation_in_certified_rank15_subgroup": exact,
                "basis_relation": list(relation) if exact else None,
                "fraction_group_law_replay": exact,
            }
        )
        candidate_records.append(record)

    cover_records = []
    score_by_indices = {indices: score for score, indices in full_frontier}
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
            "minimal_model": list(conductor["minimal_model"]),
            "conductor": str(conductor["conductor"]),
            "log_conductor": conductor["log_conductor"],
            "root_number": conductor["root_number"],
            "target_rank": TARGET_RANK,
        },
        "primary_source": PRIMARY_SOURCE,
        "exact_seed_pool": {
            "point_count": len(pool),
            "point_sha256": point_digest(pool),
            "sieve_rows": list(sieve_records),
            "height_matrix_runs": list(height_runs),
            "stable_numerical_rank": 15,
            "selected_pool_indices_one_based": list(indices),
            "selection_is_numerical_not_certification": True,
        },
        "exact_rank15_certificate": {
            "small_prime_saturation": saturation,
            "saturated_basis_sha256": point_digest(saturated_basis),
            "two_torsion_certificate_prime": two_torsion_prime,
            "finite_reduction_primes": [
                signature.prime for signature in basis_signatures
            ],
            "combined_exact_rank_over_F2": exact_basis_rank,
            "certified_algebraic_rank_lower_bound": exact_basis_rank,
        },
        "declared_budget": {
            "all_nonzero_certified_mod2_classes_identity_scored": (
                (1 << len(saturated_basis)) - 1
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
            "exact_relations_in_certified_rank15_subgroup": sum(
                exact for _, exact in proposals
            ),
            "unresolved_by_exact_relation_replay": len(unresolved),
            "combined_finite_reduction_rank": augmented_rank,
            "certified_new_directions": certified_gain,
            "certified_rank_lower_bound_after_search": certified_rank,
            "target_rank_21_achieved": certified_rank >= TARGET_RANK,
            "candidate_points": candidate_records,
        },
        "interpretation": {
            "exact": (
                "the rank-15 baseline, every returned point membership, every "
                "displayed subgroup relation, and the final finite-reduction "
                "rank are exact"
            ),
            "bounded": (
                "the finite chart boxes are exhausted; this is not a rank upper "
                "bound or a complete search of degree-two models"
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
