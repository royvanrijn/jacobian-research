#!/usr/bin/env python3
"""Exact denominator sieve around 128 fresh signed weight-3 centers.

Enumerate the global-sign quotient of all signed sums of three distinct
published rank-29 basis points.  Of the 14,616 exact abscissas, five already
occur among the 32 previously discovered higher-weight companions and are
excluded; none collide internally.  Before using any sieve outcome, retain
the 128 smallest exact projective-bit-height centers and exhaust

``x=x_Q+k/b^2, 50001<=b<=100000, 0<|k|<=16384``.

Exact rational separation proves disjointness from all public, companion,
weight-at-most-two, and prior direct-x boxes.  The modular sieve terminates in
an exact integer-square test and routes every point to the existing exact
rank-30 classifier.  This is a bounded search, not a rank upper bound.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
import json
from math import gcd, isqrt
from pathlib import Path
import platform
import shlex
import sys
import time
from typing import Any

from elkies_klagsbrun_rank29 import PUBLISHED_POINTS, point_on_general_curve
from pari_bridge import pari_version
from search_elkies_klagsbrun_rank30 import exact_linear_combination, point_negate
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
from search_elkies_klagsbrun_rank30_subgroup_center_remainder_sieve import (
    exact_cross_separation,
)
from search_elkies_klagsbrun_rank30_subgroup_center_sieve import (
    exact_prior_center_separation,
    generate_subgroup_centers,
)


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]
WEIGHT2_SCRIPT = (
    ROOT / "elliptic-curves/cas/search_elkies_klagsbrun_rank30_subgroup_center_sieve.py"
)
EXPECTED_WEIGHT2_SCRIPT_SHA256 = (
    "ddcc4b314a21809c56c87caa55d13a5bc595f09b11b9b6765a454c8632d018ab"
)
WEIGHT2_FIRST_ARTIFACT = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_elkies_klagsbrun_rank30_subgroup_center_sieve.json"
)
EXPECTED_WEIGHT2_FIRST_ARTIFACT_SHA256 = (
    "a46d126cde3db7847bfcac75afb10f4cdd15c1ffd8a79983b899e78f2e14c5cd"
)
WEIGHT2_REMAINDER_SCRIPT = (
    ROOT
    / "elliptic-curves/cas/"
    "search_elkies_klagsbrun_rank30_subgroup_center_remainder_sieve.py"
)
EXPECTED_WEIGHT2_REMAINDER_SCRIPT_SHA256 = (
    "164c7e708faea20e48797bfeb820f03f8ef6cf4cbcd96bced4e31f3f059dddab"
)
WEIGHT2_REMAINDER_ARTIFACT = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_elkies_klagsbrun_rank30_subgroup_center_remainder_sieve.json"
)
EXPECTED_WEIGHT2_REMAINDER_ARTIFACT_SHA256 = (
    "1f2253ac6e90c7ea44b17fc296d0efa088fe7d3fdb2c869c8f237cf5d3f53b5d"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_elkies_klagsbrun_rank30_weight3_center_sieve.json"
)
RAW_WEIGHT3_COUNT = 14_616
EXPECTED_PRIOR_OVERLAP_COUNT = 5
EXPECTED_POPULATION_COUNT = 14_611
SELECTED_COUNT = 128


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def generate_weight3_centers() -> tuple[
    tuple[CompanionCenter, ...], tuple[CompanionCenter, ...], dict[str, Any]
]:
    weight2, _, weight2_manifest = generate_subgroup_centers()
    old, _ = load_companion_centers()
    public_x = {point[0] for point in PUBLISHED_POINTS}
    weight2_x = {center.x for center in weight2}
    old_x = {center.x for center in old}
    by_x: dict[Fraction, CompanionCenter] = {}
    overlap_records: list[dict[str, Any]] = []
    raw_count = 0
    dimension = len(PUBLISHED_POINTS)
    for indices in combinations(range(dimension), 3):
        for signs in product((1, -1), repeat=2):
            raw_count += 1
            relation = [0] * dimension
            relation[indices[0]] = 1
            relation[indices[1]] = signs[0]
            relation[indices[2]] = signs[1]
            exact_relation = tuple(relation)
            point = exact_linear_combination(exact_relation)
            if point is None or not point_on_general_curve(point):
                raise AssertionError("a signed weight-3 center is not finite/on-curve")
            sources = []
            if point[0] in public_x:
                sources.append("public")
            if point[0] in old_x:
                sources.append("prior_32_companion")
            if point[0] in weight2_x:
                sources.append("weight_at_most_2")
            label = (
                f"p{indices[0] + 1:02d}_"
                f"{'p' if signs[0] > 0 else 'm'}{indices[1] + 1:02d}_"
                f"{'p' if signs[1] > 0 else 'm'}{indices[2] + 1:02d}"
            )
            if sources:
                overlap_records.append(
                    {
                        "label": label,
                        "x": str(point[0]),
                        "published_basis_relation": list(exact_relation),
                        "prior_sources": sources,
                    }
                )
                continue
            root = isqrt(point[0].denominator)
            if root * root != point[0].denominator:
                raise AssertionError("a weight-3 x denominator is not a square")
            center = CompanionCenter(
                x=point[0],
                point=point,
                relation=exact_relation,
                source_paths=(label,),
            )
            previous = by_x.get(center.x)
            if previous is not None:
                if previous.point not in (center.point, point_negate(center.point)):
                    raise AssertionError("equal weight-3 abscissas are not inverse")
                raise AssertionError("the new weight-3 population has an x collision")
            by_x[center.x] = center
    if raw_count != RAW_WEIGHT3_COUNT:
        raise AssertionError("the raw signed weight-3 count changed")
    if len(overlap_records) != EXPECTED_PRIOR_OVERLAP_COUNT:
        raise AssertionError("the prior weight-3 overlap count changed")
    population = tuple(
        sorted(
            by_x.values(),
            key=lambda center: (
                center.bit_height,
                center.x,
                center.point[1],
                center.relation,
            ),
        )
    )
    if len(population) != EXPECTED_POPULATION_COUNT:
        raise AssertionError("the decontaminated weight-3 population changed")
    selected = population[:SELECTED_COUNT]
    records = [
        {
            "population_index": index + 1,
            "label": center.source_paths[0],
            "x": str(center.x),
            "representative_y": str(center.point[1]),
            "bit_height": center.bit_height,
            "published_basis_relation": list(center.relation),
        }
        for index, center in enumerate(population)
    ]
    manifest = {
        "construction": (
            "all signed distinct-support weight-3 relations modulo global sign"
        ),
        "raw_count": raw_count,
        "prior_overlap_count": len(overlap_records),
        "prior_overlap_records": overlap_records,
        "decontaminated_population_count": len(population),
        "decontaminated_population_sha256": sha256(
            json.dumps(records, separators=(",", ":")).encode()
        ).hexdigest(),
        "selection_rule": (
            "first 128 by (exact projective bit height, x, y, relation), fixed "
            "before modular/square outcomes"
        ),
        "selected_count": len(selected),
        "selected_sha256": sha256(
            json.dumps(records[:SELECTED_COUNT], separators=(",", ":")).encode()
        ).hexdigest(),
        "selected_records": records[:SELECTED_COUNT],
        "weight_at_most_2_population_sha256": weight2_manifest[
            "full_population_sha256"
        ],
    }
    return population, selected, manifest


def signed_coprime_counts(
    denominator_min: int, denominator_max: int, radius: int
) -> dict[int, int]:
    return {
        denominator: 2 * positive_coprime_offset_count(radius, denominator)
        for denominator in range(denominator_min, denominator_max + 1)
    }


def declared_count(
    centers: tuple[CompanionCenter, ...], counts: dict[int, int]
) -> int:
    return sum(
        counts[denominator]
        for denominator in counts
        for center in centers
        if gcd(denominator, center.denominator_root) == 1
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--denominator-min", type=int, default=50_001)
    parser.add_argument("--denominator-max", type=int, default=100_000)
    parser.add_argument("--offset-radius", type=int, default=16_384)
    parser.add_argument("--wall-cap-seconds", type=float, default=600.0)
    parser.add_argument("--progress-every-centers", type=int, default=8)
    parser.add_argument("--certificate-prime-bound", type=int, default=1000)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=500_000_000)
    args = parser.parse_args()
    if [args.denominator_min, args.denominator_max] != [50_001, 100_000]:
        raise SystemExit("this lane is pinned to 50001<=b<=100000")
    if args.offset_radius != 16_384:
        raise SystemExit("this lane is pinned to 0<|k|<=16384")
    if args.wall_cap_seconds <= 0 or args.progress_every_centers <= 0:
        raise SystemExit("resource/progress bounds must be positive")
    for path, expected in (
        (WEIGHT2_SCRIPT, EXPECTED_WEIGHT2_SCRIPT_SHA256),
        (WEIGHT2_FIRST_ARTIFACT, EXPECTED_WEIGHT2_FIRST_ARTIFACT_SHA256),
        (WEIGHT2_REMAINDER_SCRIPT, EXPECTED_WEIGHT2_REMAINDER_SCRIPT_SHA256),
        (WEIGHT2_REMAINDER_ARTIFACT, EXPECTED_WEIGHT2_REMAINDER_ARTIFACT_SHA256),
    ):
        if file_sha256(path) != expected:
            raise AssertionError(f"a pinned predecessor changed: {path.name}")

    started = time.monotonic()
    _, centers, population_manifest = generate_weight3_centers()
    weight2, _, _ = generate_subgroup_centers()
    old, _ = load_companion_centers()
    internal = exact_nonoverlap_proof(
        centers,
        denominator_min=args.denominator_min,
        offset_radius=args.offset_radius,
    )
    old_boxes = exact_prior_center_separation(
        centers,
        old,
        denominator_min=args.denominator_min,
        offset_radius=args.offset_radius,
    )
    weight2_boxes = exact_cross_separation(
        centers,
        weight2,
        denominator_min=args.denominator_min,
        offset_radius=args.offset_radius,
    )
    if not internal["all_exact_nonoverlap_checks_passed"]:
        raise AssertionError("weight-3 internal/direct nonoverlap failed")
    if not old_boxes["all_prior_center_boxes_disjoint"]:
        raise AssertionError("weight-3 prior-center nonoverlap failed")
    if not weight2_boxes["passed"]:
        raise AssertionError("weight-3/weight-at-most-two separation failed")

    counts = signed_coprime_counts(
        args.denominator_min, args.denominator_max, args.offset_radius
    )
    expected = declared_count(centers, counts)
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
    for center_index, center in enumerate(centers):
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
            processed += counts[denominator]
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
                    f"{center_index + 1}|{denominator}|{offset}\n".encode()
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
                        "center_index": center_index + 1,
                        "denominator": denominator,
                        "offset": offset,
                        "normalized_x_numerator": numerator,
                        "normalized_x_denominator_root": root_denominator,
                        "homogeneous_square_root": square_root,
                        "point": points[0],
                    }
                )
        completed = center_index + 1
        if completed % args.progress_every_centers == 0 or completed == len(centers):
            print(
                f"weight3 centers {completed}/{len(centers)}; "
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
        "artifact_kind": "exact_rank29_weight3_subgroup_center_sieve",
        "status": status,
        "target_hit": target_hit,
        "claim_scope": {
            "bounded": "the completed 128 selected weight-3 center boxes only",
            "one_pass_no_retries": True,
            "selection_was_outcome_blind": True,
            "not_a_rank_upper_bound": True,
        },
        "source": {
            "weight2_first_artifact_sha256": EXPECTED_WEIGHT2_FIRST_ARTIFACT_SHA256,
            "weight2_remainder_artifact_sha256": EXPECTED_WEIGHT2_REMAINDER_ARTIFACT_SHA256,
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
            "weight3_internal_and_direct_charts": internal,
            "public_and_prior_32_companion_boxes": old_boxes,
            "all_841_weight_at_most_2_boxes": weight2_boxes,
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
