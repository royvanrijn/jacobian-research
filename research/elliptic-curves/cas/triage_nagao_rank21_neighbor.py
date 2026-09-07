#!/usr/bin/env python3
"""Bounded rank triage for the low-conductor Nagao rank-21 neighbor.

The CRT/Gauss neighborhood search found ``T=6041/198`` with
``log(N)=170.765...`` in the six-root family used by Nagao's published
rank-21 curve.  This verifier injects the twelve exact visible Mestre points,
searches the primitive quartic with ``hyperellratpoints`` in one declared
height box, maps and deduplicates every returned nonzero point exactly, and
replays the resulting height matrix at two precisions.  One effort-zero PARI
``ellrank`` call is allowed, with a hard timeout of at most twenty seconds.

The published specialization is calibration only: its twelve visible points
are compared with Nagao's full printed set of twenty-one points.  Numerical
height rank is evidence, and a bounded search failure is never a rank bound.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import re
import shlex
import sys
from typing import Any, Sequence

from ek_k3 import rational_to_string
from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK21_CONSTRUCTION,
    RANK21_CONSTRUCTOR_PARAMETER,
    RANK21_PUBLISHED_MODEL,
    RANK21_PUBLISHED_POINTS,
    point_on_extended_weierstrass,
    primitive_quartic_coefficients,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    short_jacobian_coefficients,
)
from pari_bridge import pari_version
from search_extra_points import (
    parse_point_vector,
    parse_precisions,
    run_gp,
    signless_quartic_points,
)
from triage_nagao_rank13_finalists import (
    ellrank_probe,
    height_matrix_replay,
    point_digest,
    point_on_short_curve,
    quartic_gp_polynomial,
    stable_height_rank,
)


Q = Fraction
CANDIDATE_PARAMETER = Q(6041, 198)
TARGET_LOG_CONDUCTOR = Decimal("182.72")
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/triage_nagao_rank21_neighbor.py"
)


@dataclass(frozen=True)
class NeighborhoodRecord:
    candidate_id: str
    parameter: Fraction
    score: str
    score_last_numerical_prime: int
    conductor: int
    log_conductor: str
    root_number: int
    source_local_checks: tuple[dict[str, Any], ...]

    @property
    def below_target(self) -> bool:
        return Decimal(self.log_conductor) < TARGET_LOG_CONDUCTOR


def load_candidate(path: Path) -> NeighborhoodRecord:
    """Load and validate the pinned low-conductor neighborhood candidate."""

    data = json.loads(path.read_text())
    records = data.get("final_conductor_records")
    if not isinstance(records, list):
        raise ValueError("the input has no final_conductor_records list")
    matching = [
        record
        for record in records
        if Q(record["constructor_parameter"]) == CANDIDATE_PARAMETER
    ]
    if len(matching) != 1:
        raise ValueError("the input does not contain exactly one T=6041/198 record")
    record = matching[0]
    if record.get("rank_claim") is not None:
        raise AssertionError("the source unexpectedly promoted a rank claim")
    curve = record["curve"]
    candidate = NeighborhoodRecord(
        candidate_id=str(record["candidate_id"]),
        parameter=Q(record["constructor_parameter"]),
        score=str(record["score"]),
        score_last_numerical_prime=int(record["last_numerical_prime"]),
        conductor=int(curve["conductor"]),
        log_conductor=str(curve["log_conductor"]),
        root_number=int(curve["root_number"]),
        source_local_checks=tuple(record["local_checks"]),
    )
    if not candidate.below_target:
        raise AssertionError("the pinned candidate no longer lies below the target")
    if candidate.root_number != 1:
        raise AssertionError("the pinned candidate's root number changed")
    if [check["prime"] for check in candidate.source_local_checks] != [5, 7, 13, 17, 23]:
        raise AssertionError("the pinned local profile changed")
    return candidate


def exact_visible_seeds(
    parameter: Fraction,
) -> tuple[
    tuple[tuple[Fraction, Fraction], ...],
    tuple[tuple[Fraction, Fraction], ...],
    tuple[Fraction, ...],
]:
    """Return twelve checked quartic points and their exact Jacobian images."""

    parameter = Q(parameter)
    quartic_points = primitive_visible_points(RANK21_CONSTRUCTION, parameter)
    if len(quartic_points) != 12 or len(set(quartic_points)) != 12:
        raise AssertionError("the twelve visible points collided")
    if len({point[0] for point in quartic_points}) != 12:
        raise AssertionError("two visible points acquired the same abscissa")
    jacobian_points = tuple(
        quartic_point_to_short_jacobian(
            RANK21_CONSTRUCTION, parameter, quartic_point
        )
        for quartic_point in quartic_points
    )
    coefficients = short_jacobian_coefficients(RANK21_CONSTRUCTION, parameter)
    if any(not point_on_short_curve(coefficients, point) for point in jacobian_points):
        raise AssertionError("an exact visible image missed the Jacobian")
    if len({point[0] for point in jacobian_points}) != 12:
        raise AssertionError("two Jacobian seed points differ only by sign")
    return quartic_points, jacobian_points, coefficients


def bounded_quartic_points(
    parameter: Fraction,
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], float, int]:
    if height_bound < 1:
        raise ValueError("the quartic height bound must be positive")
    coefficients = primitive_quartic_coefficients(RANK21_CONSTRUCTION, parameter)
    program = "\n".join(
        (
            f"Q={quartic_gp_polynomial(coefficients)};",
            "gettime();",
            f"R=hyperellratpoints(Q,{height_bound});",
            'print("PARI_MILLISECONDS ",gettime());',
            'print("POINTS ",R);',
            "quit",
        )
    ) + "\n"
    output, wall_seconds = run_gp(
        program, timeout=timeout, stack_bytes=stack_bytes
    )
    milliseconds = re.search(r"PARI_MILLISECONDS (\d+)", output)
    if milliseconds is None or "POINTS " not in output:
        raise AssertionError("PARI omitted bounded-search output")
    return (
        parse_point_vector(output.split("POINTS ", 1)[1]),
        wall_seconds,
        int(milliseconds.group(1)),
    )


def map_and_deduplicate(
    parameter: Fraction,
    raw_points: Sequence[tuple[Fraction, Fraction]],
    seed_quartic: Sequence[tuple[Fraction, Fraction]],
    seed_jacobian: Sequence[tuple[Fraction, Fraction]],
) -> tuple[
    tuple[dict[str, Any], ...],
    tuple[tuple[Fraction, Fraction], ...],
    tuple[dict[str, str], ...],
]:
    """Map one point per quartic abscissa and retain new Jacobian signs once."""

    parameter = Q(parameter)
    signless = signless_quartic_points(tuple(raw_points))
    seed_quartic_x = {point[0] for point in seed_quartic}
    seen_jacobian_x = {point[0] for point in seed_jacobian}
    coefficients = short_jacobian_coefficients(RANK21_CONSTRUCTION, parameter)
    records: list[dict[str, Any]] = []
    new_images: list[tuple[Fraction, Fraction]] = []
    zero_ordinates: list[dict[str, str]] = []
    for quartic_point in signless:
        if quartic_point[1] == 0:
            zero_ordinates.append(
                {
                    "quartic_x": rational_to_string(quartic_point[0]),
                    "quartic_z": "0",
                }
            )
            continue
        jacobian_point = quartic_point_to_short_jacobian(
            RANK21_CONSTRUCTION, parameter, quartic_point
        )
        if not point_on_short_curve(coefficients, jacobian_point):
            raise AssertionError("a bounded quartic point missed the Jacobian")
        visible_x = quartic_point[0] in seed_quartic_x
        duplicate_jacobian_sign_pair = jacobian_point[0] in seen_jacobian_x
        records.append(
            {
                "quartic_x": rational_to_string(quartic_point[0]),
                "quartic_z": rational_to_string(quartic_point[1]),
                "jacobian_x": rational_to_string(jacobian_point[0]),
                "jacobian_y": rational_to_string(jacobian_point[1]),
                "visible_section_abscissa": visible_x,
                "duplicate_seed_or_prior_sign_pair": duplicate_jacobian_sign_pair,
                "exact_quartic_and_jacobian_membership_checked": True,
            }
        )
        if not duplicate_jacobian_sign_pair:
            seen_jacobian_x.add(jacobian_point[0])
            new_images.append(jacobian_point)
    return tuple(records), tuple(new_images), tuple(zero_ordinates)


def numerical_subset(
    points: Sequence[tuple[Fraction, Fraction]],
    height_runs: Sequence[dict[str, Any]],
) -> tuple[tuple[Fraction, Fraction], ...]:
    stable_height_rank(height_runs)
    indices = height_runs[-1]["subset_indices_one_based"]
    return tuple(points[index - 1] for index in indices)


def exact_point_records(
    points: Sequence[tuple[Fraction, Fraction]],
    provenance: Sequence[str],
    height_runs: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    indices = height_runs[-1]["subset_indices_one_based"]
    return [
        {
            "pool_index_one_based": index,
            "provenance": provenance[index - 1],
            "jacobian_x": rational_to_string(points[index - 1][0]),
            "jacobian_y": rational_to_string(points[index - 1][1]),
            "exact_jacobian_membership_checked": True,
        }
        for index in indices
    ]


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_nagao_rank21_neighborhood.json",
    )
    parser.add_argument("--height-bound", type=int, default=50_000)
    parser.add_argument("--precisions", type=parse_precisions, default=(72, 120))
    parser.add_argument("--search-timeout", type=float, default=30.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--rank-timeout", type=float, default=20.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--rank-stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_nagao_rank21_neighbor_triage.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.height_bound < 1:
        raise SystemExit("--height-bound must be positive")
    if min(args.search_timeout, args.height_timeout, args.rank_timeout) <= 0:
        raise SystemExit("all timeouts must be positive")
    if args.rank_timeout > 20:
        raise SystemExit("--rank-timeout is intentionally capped at 20 seconds")
    if min(args.stack_bytes, args.rank_stack_bytes) < 8_000_000:
        raise SystemExit("PARI stack bounds must be at least 8,000,000 bytes")

    candidate = load_candidate(args.input)
    seed_quartic, seed_jacobian, coefficients = exact_visible_seeds(
        candidate.parameter
    )
    baseline_runs = height_matrix_replay(
        coefficients,
        seed_jacobian,
        precisions=args.precisions,
        timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    baseline_rank = stable_height_rank(baseline_runs)

    raw_points, search_wall, search_milliseconds = bounded_quartic_points(
        candidate.parameter,
        height_bound=args.height_bound,
        timeout=args.search_timeout,
        stack_bytes=args.stack_bytes,
    )
    mapped_records, new_images, zero_ordinates = map_and_deduplicate(
        candidate.parameter, raw_points, seed_quartic, seed_jacobian
    )
    pool = seed_jacobian + new_images
    provenance = tuple(f"visible-section-{index}" for index in range(1, 13)) + tuple(
        f"bounded-search-{index}" for index in range(1, len(new_images) + 1)
    )
    augmented_runs = height_matrix_replay(
        coefficients,
        pool,
        precisions=args.precisions,
        timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    augmented_rank = stable_height_rank(augmented_runs)
    selected_points = numerical_subset(pool, augmented_runs)
    ellrank = ellrank_probe(
        coefficients,
        selected_points,
        timeout=args.rank_timeout,
        stack_bytes=args.rank_stack_bytes,
    )

    published_seed_quartic, published_seed_jacobian, published_coefficients = (
        exact_visible_seeds(RANK21_CONSTRUCTOR_PARAMETER)
    )
    published_seed_runs = height_matrix_replay(
        published_coefficients,
        published_seed_jacobian,
        precisions=args.precisions,
        timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    published_seed_rank = stable_height_rank(published_seed_runs)
    if not all(
        point_on_extended_weierstrass(RANK21_PUBLISHED_MODEL, point)
        for point in RANK21_PUBLISHED_POINTS
    ):
        raise AssertionError("a printed rank-21 point missed Nagao's model")
    published_full_runs = height_matrix_replay(
        RANK21_PUBLISHED_MODEL,
        RANK21_PUBLISHED_POINTS,
        precisions=args.precisions,
        timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    published_full_rank = stable_height_rank(published_full_runs)

    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded exact-point and precision-stable numerical-height triage; "
            "the effort-zero rank probe is strictly time bounded, no algebraic "
            "rank is inferred from the neighborhood, and no target hit is claimed"
        ),
        "primary_source": PRIMARY_SOURCE,
        "input": {
            "path": str(args.input),
            "sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        },
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "certified_hits": [],
        },
        "candidate": {
            "candidate_id": candidate.candidate_id,
            "constructor_parameter": rational_to_string(candidate.parameter),
            "source_score": candidate.score,
            "source_score_last_numerical_prime": candidate.score_last_numerical_prime,
            "conductor": candidate.conductor,
            "log_conductor": candidate.log_conductor,
            "below_strict_log_conductor_target": candidate.below_target,
            "root_number": candidate.root_number,
            "source_local_checks": list(candidate.source_local_checks),
            "rank_claim": None,
            "exact_visible_seeds": {
                "quartic_point_count": len(seed_quartic),
                "distinct_quartic_abscissas": len({point[0] for point in seed_quartic}),
                "jacobian_image_count": len(seed_jacobian),
                "distinct_jacobian_sign_pairs": len({point[0] for point in seed_jacobian}),
                "quartic_point_sha256": point_digest(seed_quartic),
                "jacobian_point_sha256": point_digest(seed_jacobian),
                "all_points_checked_exactly": True,
            },
            "baseline_height_matrix_runs": list(baseline_runs),
            "stable_baseline_numerical_rank": baseline_rank,
            "bounded_quartic_search": {
                "naive_height_bound": args.height_bound,
                "scope": (
                    "PARI hyperellratpoints affine rational points on the exact "
                    "primitive quartic up to the declared naive height"
                ),
                "wall_seconds": search_wall,
                "pari_reported_milliseconds": search_milliseconds,
                "signed_points_found": len(raw_points),
                "distinct_quartic_abscissas": len(signless_quartic_points(raw_points)),
                "mapped_nonzero_signless_points": len(mapped_records),
                "visible_seed_abscissas_returned": sum(
                    bool(record["visible_section_abscissa"])
                    for record in mapped_records
                ),
                "new_distinct_jacobian_sign_pairs": len(new_images),
                "zero_ordinate_points_not_mapped": list(zero_ordinates),
                "mapped_point_records": list(mapped_records),
                "exact_mapping_and_sign_pair_deduplication_checked": True,
            },
            "augmented_pool_point_count": len(pool),
            "augmented_height_matrix_runs": list(augmented_runs),
            "stable_augmented_numerical_rank": augmented_rank,
            "numerical_rank_gain_over_visible_seed_pool": augmented_rank - baseline_rank,
            "explicit_numerically_independent_subset": exact_point_records(
                pool, provenance, augmented_runs
            ),
            "pari_ellrank_effort_zero": ellrank,
            "interpretation": (
                "all stored points satisfy exact equations, but height-matrix "
                "rank is numerical and the H-bounded search supplies no upper "
                "bound on rational points or Mordell-Weil rank"
            ),
        },
        "published_record_calibration": {
            "constructor_parameter": rational_to_string(
                RANK21_CONSTRUCTOR_PARAMETER
            ),
            "visible_quartic_point_count": len(published_seed_quartic),
            "visible_jacobian_point_count": len(published_seed_jacobian),
            "visible_seed_height_matrix_runs": list(published_seed_runs),
            "stable_visible_seed_numerical_rank": published_seed_rank,
            "printed_point_count": len(RANK21_PUBLISHED_POINTS),
            "all_printed_points_checked_exactly_on_printed_model": True,
            "printed_point_height_matrix_runs": list(published_full_runs),
            "stable_printed_point_numerical_rank": published_full_rank,
            "independence_status": (
                "Nagao's published independence result is cited; the numerical "
                "height replay is calibration and not a replacement proof"
            ),
            "calibration_lesson": (
                "the published rank-21 specialization has only numerical rank 11 "
                "inside the same twelve visible seeds, so failure to find new "
                "low-height points at the neighbor cannot exclude hidden high rank"
            ),
        },
        "summary": {
            "candidate_stable_visible_seed_numerical_rank": baseline_rank,
            "candidate_signed_points_found_at_height_bound": len(raw_points),
            "candidate_new_jacobian_sign_pairs_at_height_bound": len(new_images),
            "candidate_stable_augmented_numerical_rank": augmented_rank,
            "ellrank_probe_status": ellrank["status"],
            "published_visible_seed_numerical_rank": published_seed_rank,
            "published_printed_point_numerical_rank": published_full_rank,
            "target_hit": False,
        },
        "bounds": {
            "quartic_naive_height_bound": args.height_bound,
            "search_timeout_seconds": args.search_timeout,
            "height_timeout_seconds_per_replay": args.height_timeout,
            "height_decimal_precisions": list(args.precisions),
            "ellrank_effort": 0,
            "ellrank_timeout_seconds": args.rank_timeout,
            "pari_stack_bytes": args.stack_bytes,
            "pari_rank_stack_bytes": args.rank_stack_bytes,
            "process_policy": (
                "every PARI invocation is a synchronous foreground subprocess; "
                "Python subprocess timeouts kill and wait for that child"
            ),
        },
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "invocation": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    print(
        f"T={candidate.parameter} visible_rank={baseline_rank} "
        f"signed_H{args.height_bound}={len(raw_points)} "
        f"new={len(new_images)} pool_rank={augmented_rank} "
        f"ellrank={ellrank['status']}"
    )
    print(
        f"published calibration: visible_rank={published_seed_rank}, "
        f"printed_point_rank={published_full_rank}"
    )


if __name__ == "__main__":
    main()
