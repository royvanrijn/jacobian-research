#!/usr/bin/env python3
"""Search extreme small-denominator slabs on the best section-7 covers.

This is a complementary continuation of
``search_nagao_rank20_t5081_direction.py``.  Its input pass exhausted all
``2^20-1`` represented mod-2 classes but searched each selected alternate
quartic only in balanced cross-ratio boxes.  Here the best identity-score
covers receive two much more skew boxes: numerator height ``10^10`` with
denominator at most 10, then numerator height ``10^9`` with denominator
11--100.  Every PARI process is a fresh process group with a strict timeout
and no retry.

All points are checked over ``QQ`` and decontaminated against the certified
rank-20 basis, the 18 predeclared generic seeds, and every point proved
dependent by the input direction pass.  Only exact Fraction relation replay
or exact finite-reduction signatures determine whether a new direction was
found.
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
from typing import Any, Iterable, Sequence

from alternate_quartic_covers import alternate_cover, short_subset_sum
from certify_nagao_rank20_t5081 import (
    CONSTRUCTION,
    EXPECTED_CONDUCTOR,
    PARAMETER_T,
    ROOTS,
    exact_curve_data,
)
from ek_k3 import rational_to_string
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)
from nagao_1994 import quartic_point_to_short_jacobian, quartic_value
from pari_bridge import pari_version
from search_nagao_rank20_t5081_direction import (
    EXPECTED_INPUT_CERTIFICATE_SHA256,
    INPUT_RANK,
    TARGET_RANK,
    exact_relation_proposals,
    generic_companion_quartic_points,
    load_exact_basis,
    point_digest,
    point_record,
    sha256_file,
    signature_records,
)
from search_nagao_rank21_t956_skew import search_original_quartic
from triage_nagao_rank13_finalists import point_on_short_curve


Q = Fraction
EXPECTED_DIRECTION_ARTIFACT_SHA256 = (
    "9462cfbd7a909f7922a47f1d4662452091500553801856a12fab00483d9cd997"
)
DEFAULT_COVER_COUNT = 8
SKEW_BOXES = (
    ("d000001_000010", "[10000000000,10]", 10_000_000_000, 1, 10),
    ("d000011_000100", "[1000000000,[11,100]]", 1_000_000_000, 11, 100),
)
GENERIC_QUADRATIC_SECTIONS = (
    (
        "quadratic-plus-56/5373",
        Q(56, 5373),
        Q(1389190, 5373),
        (Q(0), Q(684218797630, 9623043), Q(0), Q(171853351, 9623043), Q(0), Q(3136, 9623043)),
    ),
    (
        "quadratic-minus-22/5373",
        Q(-22, 5373),
        Q(1389190, 5373),
        (Q(0), Q(638541742570, 9623043), Q(0), Q(-24743777, 9623043), Q(0), Q(484, 9623043)),
    ),
    (
        "quadratic-minus-34/5373",
        Q(-34, 5373),
        Q(1389190, 5373),
        (Q(0), Q(-631221199130, 9623043), Q(0), Q(-2763929, 9623043), Q(0), Q(1156, 9623043)),
    ),
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank20_t5081_cover_skew.py "
    "--output artifacts/generated-results/elliptic_nagao_rank20_t5081_cover_skew.json"
)


def parse_curve_point(record: dict[str, Any]) -> tuple[Fraction, Fraction]:
    return Q(record["curve_x"]), Q(record["curve_y"])


def generic_quadratic_quartic_points(
    parameter: Fraction, quartic: Sequence[Fraction]
) -> tuple[tuple[str, tuple[Fraction, Fraction]], ...]:
    records = []
    for label, quadratic_coefficient, constant, ordinate_coefficients in GENERIC_QUADRATIC_SECTIONS:
        x_value = quadratic_coefficient * parameter**2 + constant
        y_value = Q(0)
        for coefficient in reversed(ordinate_coefficients):
            y_value = y_value * parameter + coefficient
        point = x_value, y_value
        if y_value**2 != quartic_value(quartic, x_value):
            raise AssertionError(f"generic quadratic section {label} failed exactly")
        records.append((label, point))
    return tuple(records)


def load_direction_artifact(path: Path) -> dict[str, Any]:
    if sha256_file(path) != EXPECTED_DIRECTION_ARTIFACT_SHA256:
        raise AssertionError("the pinned direction-search artifact changed")
    data = json.loads(path.read_text(encoding="utf-8"))
    scan = data["full_mod2_class_scan"]
    results = data["results"]
    if (
        scan["nonzero_classes_scored"] != 2**INPUT_RANK - 1
        or results["certified_rank_lower_bound_after_search"] != INPUT_RANK
        or results["unresolved_by_exact_relation_replay"] != 0
    ):
        raise AssertionError("the input direction-search conclusion changed")
    return data


def reconstruct_best_covers(
    coefficients: Sequence[Fraction],
    basis: Sequence[tuple[Fraction, Fraction]],
    direction_data: dict[str, Any],
    count: int,
) -> tuple[tuple[dict[str, Any], Any], ...]:
    records = direction_data["cover_plans"][:count]
    answer = []
    for record in records:
        subset = tuple(index - 1 for index in record["subset_indices_one_based"])
        base_point = short_subset_sum(coefficients, basis, subset)
        if base_point is None:
            raise AssertionError("a selected cover base point vanished")
        cover = alternate_cover(coefficients, base_point)
        expected_coefficients = tuple(Q(value) for value in record["quartic_coefficients_ascending"])
        if cover.coefficients != expected_coefficients:
            raise AssertionError("a pinned alternate quartic changed")
        answer.append((record, cover))
    if len(answer) != count:
        raise AssertionError("the requested cover count exceeds the pinned frontier")
    return tuple(answer)


def run_cover_skew(
    cover: Any,
    height_specification: str,
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[dict[str, Any], tuple[tuple[Fraction, Fraction], ...]]:
    raw, process = search_original_quartic(
        cover.coefficients,
        height_specification,
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    mapped = []
    for cover_point in raw:
        curve_point = cover.cover_point_to_curve(cover_point)
        if not point_on_short_curve(cover.short_coefficients, curve_point):
            raise AssertionError("a skew-cover point left the exact curve")
        mapped.append(curve_point)
    return process, tuple(mapped)


def build_search(args: argparse.Namespace) -> dict[str, Any]:
    coefficients, basis, _ = load_exact_basis(args.certificate_input)
    direction_data = load_direction_artifact(args.direction_input)
    cover_plans = reconstruct_best_covers(
        coefficients, basis, direction_data, args.cover_count
    )
    baseline_signatures = find_mod2_reduction_certificate(
        coefficients, basis, prime_bound=500
    )
    baseline_rank = combined_mod2_rank(baseline_signatures, len(basis))
    if baseline_rank != INPUT_RANK:
        raise AssertionError("the certified basis lost finite-reduction rank")

    prior_points = tuple(
        parse_curve_point(record)
        for record in direction_data["results"]["candidate_points"]
    )
    quartic, visible_jacobian, _ = exact_curve_data()
    companions = generic_companion_quartic_points(PARAMETER_T, quartic)
    companion_jacobian = tuple(
        quartic_point_to_short_jacobian(CONSTRUCTION, PARAMETER_T, point)
        for _, point in companions
    )
    quadratic_sections = generic_quadratic_quartic_points(PARAMETER_T, quartic)
    quadratic_jacobian = tuple(
        quartic_point_to_short_jacobian(CONSTRUCTION, PARAMETER_T, point)
        for _, point in quadratic_sections
    )
    if any(not point_on_short_curve(coefficients, point) for point in quadratic_jacobian):
        raise AssertionError("a generic quadratic image missed the exact Jacobian")
    quadratic_relations = exact_relation_proposals(
        coefficients,
        basis,
        quadratic_jacobian,
        timeout=args.relation_timeout,
        stack_bytes=args.stack_bytes,
        batch_size=args.relation_batch_size,
    )
    if not all(exact for _, exact, _ in quadratic_relations):
        raise AssertionError("a generic quadratic section did not replay in rank 20")
    known_points = (
        basis
        + visible_jacobian
        + companion_jacobian
        + quadratic_jacobian
        + prior_points
    )
    known_with_signs = {
        point
        for known in known_points
        for point in (known, (known[0], -known[1]))
    }

    discoveries: dict[tuple[Fraction, Fraction], set[str]] = {}
    runs = []
    for cover_index, (cover_record, cover) in enumerate(cover_plans, start=1):
        for box_id, height, numerator, denominator_lower, denominator_upper in SKEW_BOXES:
            process, mapped = run_cover_skew(
                cover,
                height,
                timeout=args.timeout,
                stack_bytes=args.stack_bytes,
            )
            before = len(discoveries)
            source = f"{cover_record['id']}:{box_id}"
            for point in mapped:
                discoveries.setdefault(point, set()).add(source)
            runs.append(
                {
                    "cover_id": cover_record["id"],
                    "cover_frontier_position_one_based": cover_index,
                    "cover_subset_indices_one_based": cover_record["subset_indices_one_based"],
                    "identity_score_maximum_known_t_projective_bit_length": cover_record[
                        "identity_score_maximum_known_t_projective_bit_length"
                    ],
                    "box_id": box_id,
                    "numerator_absolute_bound": numerator,
                    "denominator_lower_bound": denominator_lower,
                    "denominator_upper_bound": denominator_upper,
                    "mapped_exact_affine_points_with_sign": len(mapped),
                    "new_global_exact_affine_points_with_sign": len(discoveries) - before,
                    **process,
                }
            )
        print(
            f"cover-skew {cover_index}/{len(cover_plans)} "
            f"{cover_record['id']} discoveries={len(discoveries)}",
            flush=True,
        )

    prior_or_seed_hits = sum(point in known_with_signs for point in discoveries)
    candidates = tuple(
        sorted(
            (point for point in discoveries if point not in known_with_signs),
            key=lambda point: (
                max(abs(point[0].numerator), point[0].denominator),
                max(abs(point[1].numerator), point[1].denominator),
                point,
            ),
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
    signatures = find_mod2_reduction_certificate(
        coefficients, basis + unresolved, prime_bound=1_000
    )
    augmented_rank = combined_mod2_rank(signatures, len(basis) + len(unresolved))
    certified_gain = max(0, augmented_rank - baseline_rank)

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

    script_path = Path(__file__).resolve()
    return {
        "schema_version": 1,
        "status": "bounded section-7 alternate-cover skew search complete",
        "candidate": {
            "constructor_parameter_T": rational_to_string(PARAMETER_T),
            "roots": list(ROOTS),
            "conductor": str(EXPECTED_CONDUCTOR),
            "certified_rank_lower_bound_before_search": baseline_rank,
            "target_rank": TARGET_RANK,
        },
        "inputs": {
            "rank20_certificate": {
                "path": str(args.certificate_input),
                "sha256": sha256_file(args.certificate_input),
                "expected_sha256": EXPECTED_INPUT_CERTIFICATE_SHA256,
            },
            "direction_search": {
                "path": str(args.direction_input),
                "sha256": sha256_file(args.direction_input),
                "expected_sha256": EXPECTED_DIRECTION_ARTIFACT_SHA256,
                "prior_exact_dependent_candidate_count": len(prior_points),
            },
        },
        "search": {
            "cover_count": len(cover_plans),
            "boxes_per_cover": len(SKEW_BOXES),
            "declared_process_call_count": len(cover_plans) * len(SKEW_BOXES),
            "runs": runs,
        },
        "generic_quadratic_seed_decontamination": {
            "section_count": len(quadratic_sections),
            "generic_status": (
                "symbolically proven dependent over QQ(T); the fixed-fiber "
                "relations below are independently replayed over Fraction"
            ),
            "all_specialized_images_replayed_in_certified_rank20_basis": True,
            "already_present_among_prior_dependent_points_with_sign": sum(
                point in {
                    candidate
                    for prior in prior_points
                    for candidate in (prior, (prior[0], -prior[1]))
                }
                for point in quadratic_jacobian
            ),
            "sections": [
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
                for (label, quartic_point), jacobian_point, proposal in zip(
                    quadratic_sections, quadratic_jacobian, quadratic_relations
                )
            ],
        },
        "results": {
            "distinct_exact_affine_points_with_sign": len(discoveries),
            "prior_or_predeclared_seed_hits_with_sign": prior_or_seed_hits,
            "new_decontaminated_candidate_count": len(candidates),
            "candidate_point_sha256": point_digest(candidates),
            "exact_relations_in_certified_rank20_subgroup": sum(
                exact for _, exact, _ in proposals
            ),
            "unresolved_by_exact_relation_replay": len(unresolved),
            "augmented_finite_reduction_signatures": signature_records(signatures),
            "combined_exact_finite_reduction_rank": augmented_rank,
            "certified_new_directions": certified_gain,
            "certified_rank_lower_bound_after_search": baseline_rank + certified_gain,
            "target_rank_21_achieved": baseline_rank + certified_gain >= TARGET_RANK,
            "candidate_points": candidate_records,
        },
        "bounded_scope": {
            "one_pass_no_retry": True,
            "fresh_foreground_process_group_per_call": True,
            "all_subprocess_timeouts_at_most_60_seconds": True,
            "negative_search_is_not_a_rank_upper_bound": True,
        },
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
    generated = root / "artifacts" / "generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--certificate-input",
        type=Path,
        default=generated / "elliptic_nagao_rank20_t5081_rank20_certificate.json",
    )
    parser.add_argument(
        "--direction-input",
        type=Path,
        default=generated / "elliptic_nagao_rank20_t5081_direction.json",
    )
    parser.add_argument("--cover-count", type=int, default=DEFAULT_COVER_COUNT)
    parser.add_argument("--timeout", type=float, default=45.0)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--relation-batch-size", type=int, default=40)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_nagao_rank20_t5081_cover_skew.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 1 <= args.cover_count <= 20:
        raise SystemExit("--cover-count must lie in [1,20]")
    if not 0 < args.timeout <= 60 or not 0 < args.relation_timeout <= 60:
        raise SystemExit("all subprocess timeouts must lie in (0,60]")
    if args.relation_batch_size <= 0:
        raise SystemExit("--relation-batch-size must be positive")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes is too small")
    result = build_search(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"wrote {args.output}: new_candidates="
        f"{result['results']['new_decontaminated_candidate_count']} "
        f"unresolved={result['results']['unresolved_by_exact_relation_replay']} "
        f"certified_rank={result['results']['certified_rank_lower_bound_after_search']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
