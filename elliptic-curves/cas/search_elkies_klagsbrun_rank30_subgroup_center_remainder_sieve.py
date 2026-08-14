#!/usr/bin/env python3
"""Exhaust the remaining 713 weight-at-most-two rank-29 centers.

The predecessor selected the 128 smallest exact projective-bit-height centers
from all 841 doubles and signed pair sums of the published rank-29 basis.  This
standalone continuation takes the exact complementary 713 centers and searches
the identical disjoint denominator boxes

``x=x_Q+k/b^2, 50001<=b<=100000, 0<|k|<=16384``.

The predecessor artifact is hash-pinned and exact separation is reproved both
internally and against its 128 centers, all older public/companion boxes, and
the prior direct x charts.  The quadratic-residue sieve and terminal exact
integer-square checks are unchanged.  This is a bounded search, not a rank
upper bound.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, isqrt
from pathlib import Path
import platform
import shlex
import sys
import time
from typing import Any

from pari_bridge import pari_version
from search_elkies_klagsbrun_rank30_companion_center_sieve import (
    CompanionCenter,
    allowed_offset_mask,
    exact_nonoverlap_proof,
    load_companion_centers,
    normalized_abscissa,
)
from search_elkies_klagsbrun_rank30_denominator_sieve import (
    SIEVE_PRIMES,
    build_offset_residue_masks,
    classify_point,
    homogeneous_square_value,
    map_square_abscissa,
    positive_coprime_offset_count,
    small_companion_lookup,
)
from search_elkies_klagsbrun_rank30_subgroup_center_sieve import (
    EXPECTED_SOURCE_ENGINE_SHA256,
    SOURCE_ENGINE,
    exact_prior_center_separation,
    generate_subgroup_centers,
)


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]
PARENT_SCRIPT = (
    ROOT
    / "elliptic-curves/cas/"
    "search_elkies_klagsbrun_rank30_subgroup_center_sieve.py"
)
EXPECTED_PARENT_SCRIPT_SHA256 = (
    "ddcc4b314a21809c56c87caa55d13a5bc595f09b11b9b6765a454c8632d018ab"
)
PARENT_ARTIFACT = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_elkies_klagsbrun_rank30_subgroup_center_sieve.json"
)
EXPECTED_PARENT_ARTIFACT_SHA256 = (
    "a46d126cde3db7847bfcac75afb10f4cdd15c1ffd8a79983b899e78f2e14c5cd"
)
EXPECTED_FULL_POPULATION_SHA256 = (
    "079e43b25d8e61df848b02f9ba336a98eb66f68a8859262c0660110ae5ef4e0c"
)
EXPECTED_PARENT_SELECTION_SHA256 = (
    "584c037eda54b870bcbc7a50e42a01b0f08787aab6f4792b249cf6db6bb6441a"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_elkies_klagsbrun_rank30_subgroup_center_remainder_sieve.json"
)
FIRST_REMAINDER_INDEX = 128
EXPECTED_REMAINDER_COUNT = 713


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def selected_remainder() -> tuple[
    tuple[CompanionCenter, ...], tuple[CompanionCenter, ...], dict[str, Any]
]:
    population, first, selection = generate_subgroup_centers()
    if selection["full_population_sha256"] != EXPECTED_FULL_POPULATION_SHA256:
        raise AssertionError("the exact 841-center population changed")
    if selection["selected_sha256"] != EXPECTED_PARENT_SELECTION_SHA256:
        raise AssertionError("the predecessor 128-center selection changed")
    remainder = population[FIRST_REMAINDER_INDEX:]
    if len(first) != FIRST_REMAINDER_INDEX or len(remainder) != EXPECTED_REMAINDER_COUNT:
        raise AssertionError("the complementary center count changed")
    records = [
        {
            "population_index": FIRST_REMAINDER_INDEX + index + 1,
            "label": center.source_paths[0],
            "x": str(center.x),
            "representative_y": str(center.point[1]),
            "bit_height": center.bit_height,
            "published_basis_relation": list(center.relation),
        }
        for index, center in enumerate(remainder)
    ]
    return first, remainder, {
        "full_population_count": len(population),
        "full_population_sha256": EXPECTED_FULL_POPULATION_SHA256,
        "predecessor_selected_count": len(first),
        "predecessor_selection_sha256": EXPECTED_PARENT_SELECTION_SHA256,
        "remainder_count": len(remainder),
        "remainder_population_indices": [129, 841],
        "remainder_sha256": sha256(
            json.dumps(records, separators=(",", ":")).encode()
        ).hexdigest(),
        "selection_rule": "the exact complement of population indices 1..128",
        "records": records,
    }


def exact_cross_separation(
    left: tuple[CompanionCenter, ...],
    right: tuple[CompanionCenter, ...],
    *,
    denominator_min: int,
    offset_radius: int,
) -> dict[str, Any]:
    radius = Q(offset_radius, denominator_min**2)
    distance, left_index, right_index = min(
        (abs(a.x - b.x), i, j)
        for i, a in enumerate(left)
        for j, b in enumerate(right)
    )
    passed = distance > 2 * radius
    return {
        "passed": passed,
        "minimum_center_distance": str(distance),
        "remainder_center_index": left_index + 129,
        "predecessor_center_index": right_index + 1,
        "required_strict_lower_bound": str(2 * radius),
    }


def declared_count(
    centers: tuple[CompanionCenter, ...],
    denominator_min: int,
    denominator_max: int,
    signed_coprime_counts: dict[int, int],
) -> int:
    return sum(
        signed_coprime_counts[denominator]
        for denominator in range(denominator_min, denominator_max + 1)
        for center in centers
        if gcd(denominator, center.denominator_root) == 1
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--denominator-min", type=int, default=50_001)
    parser.add_argument("--denominator-max", type=int, default=100_000)
    parser.add_argument("--offset-radius", type=int, default=16_384)
    parser.add_argument("--wall-cap-seconds", type=float, default=1200.0)
    parser.add_argument("--progress-every-centers", type=int, default=25)
    parser.add_argument("--certificate-prime-bound", type=int, default=1000)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=500_000_000)
    args = parser.parse_args()
    if [args.denominator_min, args.denominator_max] != [50_001, 100_000]:
        raise SystemExit("this continuation is pinned to 50001<=b<=100000")
    if args.offset_radius != 16_384:
        raise SystemExit("this continuation is pinned to 0<|k|<=16384")
    if args.wall_cap_seconds <= 0 or args.progress_every_centers <= 0:
        raise SystemExit("resource and progress bounds must be positive")
    if args.relation_timeout <= 0 or args.stack_bytes < 64_000_000:
        raise SystemExit("invalid point-classifier bounds")
    if file_sha256(SOURCE_ENGINE) != EXPECTED_SOURCE_ENGINE_SHA256:
        raise AssertionError("the reusable denominator-sieve engine changed")
    if file_sha256(PARENT_SCRIPT) != EXPECTED_PARENT_SCRIPT_SHA256:
        raise AssertionError("the predecessor subgroup-center script changed")
    if file_sha256(PARENT_ARTIFACT) != EXPECTED_PARENT_ARTIFACT_SHA256:
        raise AssertionError("the predecessor bounded artifact changed")
    parent = json.loads(PARENT_ARTIFACT.read_text(encoding="utf-8"))
    if not parent["search_result"]["search_complete"]:
        raise AssertionError("the predecessor 128-center search is incomplete")

    started = time.monotonic()
    first, centers, population_manifest = selected_remainder()
    old_centers, _ = load_companion_centers()
    internal = exact_nonoverlap_proof(
        centers,
        denominator_min=args.denominator_min,
        offset_radius=args.offset_radius,
    )
    prior = exact_prior_center_separation(
        centers,
        old_centers,
        denominator_min=args.denominator_min,
        offset_radius=args.offset_radius,
    )
    predecessor = exact_cross_separation(
        centers,
        first,
        denominator_min=args.denominator_min,
        offset_radius=args.offset_radius,
    )
    if not internal["all_exact_nonoverlap_checks_passed"]:
        raise AssertionError("internal/direct-chart nonoverlap failed")
    if not prior["all_prior_center_boxes_disjoint"] or not predecessor["passed"]:
        raise AssertionError("cross-generation center-box nonoverlap failed")

    signed_coprime_counts = {
        denominator: 2
        * positive_coprime_offset_count(args.offset_radius, denominator)
        for denominator in range(args.denominator_min, args.denominator_max + 1)
    }
    expected = declared_count(
        centers,
        args.denominator_min,
        args.denominator_max,
        signed_coprime_counts,
    )
    residue_masks = build_offset_residue_masks(args.offset_radius)
    nonzero_mask = (1 << (2 * args.offset_radius + 1)) - 1
    nonzero_mask ^= 1 << args.offset_radius
    processed = 0
    survivors_before_primitivity = 0
    survivors = 0
    negatives = 0
    nonsquares = 0
    squares: list[dict[str, Any]] = []
    survivor_hasher = sha256()
    completed = 0
    wall_cap_reached = False
    cache_entries = 0

    for local_index, center in enumerate(centers):
        if time.monotonic() - started >= args.wall_cap_seconds:
            wall_cap_reached = True
            break
        tables: dict[int, tuple[int, ...]] = {}
        for prime in SIEVE_PRIMES:
            table = tuple(
                allowed_offset_mask(center, residue, prime, residue_masks)
                for residue in range(prime)
            )
            tables[prime] = table
            cache_entries += len(table)
        for denominator in range(args.denominator_min, args.denominator_max + 1):
            if gcd(denominator, center.denominator_root) != 1:
                continue
            processed += signed_coprime_counts[denominator]
            mask = nonzero_mask
            for prime in SIEVE_PRIMES:
                mask &= tables[prime][denominator % prime]
                if not mask:
                    break
            survivors_before_primitivity += mask.bit_count()
            while mask:
                low = mask & -mask
                bit_index = low.bit_length() - 1
                mask ^= low
                offset = bit_index - args.offset_radius
                if gcd(offset, denominator) != 1:
                    continue
                survivors += 1
                survivor_hasher.update(
                    f"{local_index + 129}|{denominator}|{offset}\n".encode()
                )
                numerator, root_denominator = normalized_abscissa(
                    center, denominator, offset
                )
                value = homogeneous_square_value(numerator, root_denominator)
                if value < 0:
                    negatives += 1
                    continue
                square_root = isqrt(value)
                if square_root * square_root != value:
                    nonsquares += 1
                    continue
                points = map_square_abscissa(
                    numerator, root_denominator, square_root
                )
                squares.append(
                    {
                        "population_center_index": local_index + 129,
                        "denominator": denominator,
                        "offset": offset,
                        "normalized_x_numerator": numerator,
                        "normalized_x_denominator_root": root_denominator,
                        "homogeneous_square_root": square_root,
                        "point": points[0],
                    }
                )
        completed = local_index + 1
        if completed % args.progress_every_centers == 0 or completed == len(centers):
            print(
                f"remainder centers {completed}/{len(centers)}; "
                f"primitive={processed}; survivors={survivors}; squares={len(squares)}",
                flush=True,
            )

    companion_lookup = None
    candidates: list[dict[str, Any]] = []
    for record in squares:
        point = record.pop("point")
        if companion_lookup is None:
            companion_lookup = small_companion_lookup()
        candidates.append(
            {
                **record,
                **classify_point(
                    point,
                    companion_lookup=companion_lookup,
                    certificate_prime_bound=args.certificate_prime_bound,
                    relation_timeout=args.relation_timeout,
                    stack_bytes=args.stack_bytes,
                ),
            }
        )
    target_hit = any(
        record["classification"] == "exact_independent_30th_point"
        for record in candidates
    )
    complete = not wall_cap_reached and completed == len(centers)
    status = (
        "exact_rank30_target_hit"
        if target_hit
        else "bounded_search_no_certified_30th_point"
        if complete
        else "bounded_search_incomplete_at_wall_cap"
    )
    script = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "artifact_kind": "exact_rank29_subgroup_center_remainder_sieve",
        "status": status,
        "target_hit": target_hit,
        "claim_scope": {
            "bounded": "the completed complementary 713 center boxes only",
            "one_pass_no_retries": True,
            "not_a_rank_upper_bound": True,
        },
        "source": {
            "parent_script": str(PARENT_SCRIPT.relative_to(ROOT)),
            "parent_script_sha256": EXPECTED_PARENT_SCRIPT_SHA256,
            "parent_artifact": str(PARENT_ARTIFACT.relative_to(ROOT)),
            "parent_artifact_sha256": EXPECTED_PARENT_ARTIFACT_SHA256,
        },
        "center_population": population_manifest,
        "parameters": {
            "denominator_interval": [args.denominator_min, args.denominator_max],
            "nonzero_offset_interval": [-args.offset_radius, args.offset_radius],
            "primitivity": "gcd(b,s)=gcd(k,b)=1",
            "sieve_primes": list(SIEVE_PRIMES),
            "wall_cap_seconds": args.wall_cap_seconds,
        },
        "exact_nonoverlap": {
            "remainder_internal_and_direct_charts": internal,
            "prior_public_and_32_companion_boxes": prior,
            "predecessor_128_center_boxes": predecessor,
        },
        "search_result": {
            "declared_primitive_candidate_count": expected,
            "processed_primitive_candidate_count": processed,
            "completed_center_count": completed,
            "search_complete": complete,
            "wall_cap_reached": wall_cap_reached,
            "modular_survivor_count_before_primitivity": survivors_before_primitivity,
            "modular_survivor_count_after_primitivity": survivors,
            "modular_survivor_manifest_sha256": survivor_hasher.hexdigest(),
            "negative_homogeneous_value_count_after_sieve": negatives,
            "exact_nonsquare_count_after_sieve": nonsquares,
            "exact_square_abscissa_count": len(squares),
            "candidate_records": candidates,
            "certified_independent_30th_point_count": sum(
                record["classification"] == "exact_independent_30th_point"
                for record in candidates
            ),
            "rank30_target_hit": target_hit,
            "streamed_residue_cache_entry_count": cache_entries,
            "wall_seconds": time.monotonic() - started,
        },
        "reproduction": {
            "command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "script": str(script.relative_to(ROOT)),
            "script_sha256": file_sha256(script),
            "python": platform.python_version(),
            "pari_gp": pari_version(),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, default=str) + "\n", encoding="utf-8"
    )
    print(f"wrote {args.output}")
    print(
        f"status={status}; centers={completed}/{len(centers)}; "
        f"processed={processed}/{expected}; survivors={survivors}; "
        f"squares={len(squares)}; target_hit={str(target_hit).lower()}",
        flush=True,
    )


if __name__ == "__main__":
    main()
