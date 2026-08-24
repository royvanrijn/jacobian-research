#!/usr/bin/env python3
"""Exact denominator sieve around 128 fresh rank-29 subgroup centers.

Generate all 29 doubles ``2 P_i`` and all 812 points ``P_i +/- P_j`` from
the published independent rank-29 subgroup.  Their 841 abscissas are exact
and distinct; none is public or among the 32 companion centers used by the
earlier sieve.  Select the 128 smallest exact projective-bit-height centers,
before seeing any sieve outcome, and exhaust

``x = x_Q + k/b^2, 50001 <= b <= 100000, 0 < |k| <= 16384``

under the usual denominator primitivity conditions.  Exact separation checks
prove that these boxes miss every completed public/companion center box and
the prior direct x charts.  Modular quadratic-residue masks are followed by
an exact integer-square test.  A returned point is immediately passed to the
existing finite-reduction rank-30 classifier.

This is a finite negative-search certificate for the declared boxes, not a
rank upper bound.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
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
from search_elkies_klagsbrun_rank30 import (
    exact_linear_combination,
    point_add,
    point_multiply,
    point_negate,
)
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


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_elkies_klagsbrun_rank30_subgroup_center_sieve.json"
)
SOURCE_ENGINE = (
    ROOT
    / "elliptic-curves"
    / "cas"
    / "search_elkies_klagsbrun_rank30_companion_center_sieve.py"
)
EXPECTED_SOURCE_ENGINE_SHA256 = (
    "89f63d500bdea8b3b6602d4321b2395f76740c0123e6583cb9adce8baf91e680"
)
PRIOR_PUBLIC_BASE = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_elkies_klagsbrun_rank30_denominator_sieve.json"
)
PRIOR_PUBLIC_EXTENSION = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_elkies_klagsbrun_rank30_denominator_sieve_b500000.json"
)
PRIOR_PUBLIC_ANNULUS = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_elkies_klagsbrun_rank30_offset_annulus.json"
)
PRIOR_COMPANION_BASE = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_elkies_klagsbrun_rank30_companion_center_sieve.json"
)
PRIOR_COMPANION_EXTENSION = (
    ROOT
    / "artifacts/generated-results/"
    "elliptic_elkies_klagsbrun_rank30_companion_center_sieve_b500000.json"
)
PRIOR_ARTIFACTS = (
    PRIOR_PUBLIC_BASE,
    PRIOR_PUBLIC_EXTENSION,
    PRIOR_PUBLIC_ANNULUS,
    PRIOR_COMPANION_BASE,
    PRIOR_COMPANION_EXTENSION,
)
GENERATED_CENTER_COUNT = 841
SELECTED_CENTER_COUNT = 128
PRIOR_DENOMINATOR_MIN = 3163
PRIOR_PUBLIC_OFFSET_MAX = 65_536
PRIOR_COMPANION_OFFSET_MAX = 16_384


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relation_point(relation: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    point = exact_linear_combination(relation)
    if point is None or not point_on_general_curve(point):
        raise AssertionError("a generated subgroup center is not a finite curve point")
    return point


def generate_subgroup_centers() -> tuple[
    tuple[CompanionCenter, ...], tuple[CompanionCenter, ...], dict[str, Any]
]:
    """Return the full 841-center population and its blind top-128 tranche."""

    old_centers, source_inventory = load_companion_centers()
    excluded_x = {point[0] for point in PUBLISHED_POINTS}
    excluded_x.update(center.x for center in old_centers)
    relation_records: list[tuple[tuple[int, ...], str]] = []
    dimension = len(PUBLISHED_POINTS)
    for index in range(dimension):
        relation = [0] * dimension
        relation[index] = 2
        relation_records.append((tuple(relation), f"double_p{index + 1:02d}"))
    for left, right in combinations(range(dimension), 2):
        for sign, suffix in ((1, "plus"), (-1, "minus")):
            relation = [0] * dimension
            relation[left] = 1
            relation[right] = sign
            relation_records.append(
                (
                    tuple(relation),
                    f"p{left + 1:02d}_{suffix}_p{right + 1:02d}",
                )
            )
    if len(relation_records) != GENERATED_CENTER_COUNT:
        raise AssertionError("the generated relation population changed")

    by_x: dict[Fraction, CompanionCenter] = {}
    for relation, label in relation_records:
        point = relation_point(relation)
        root = isqrt(point[0].denominator)
        if root * root != point[0].denominator:
            raise AssertionError("a generated x denominator is not a square")
        center = CompanionCenter(
            x=point[0],
            point=point,
            relation=relation,
            source_paths=(label,),
        )
        previous = by_x.get(center.x)
        if previous is not None:
            if previous.point not in (center.point, point_negate(center.point)):
                raise AssertionError("equal generated abscissas are not inverse")
            raise AssertionError("the 841 generated abscissas are no longer distinct")
        by_x[center.x] = center
    overlap = sorted(set(by_x).intersection(excluded_x))
    if overlap:
        raise AssertionError("a generated center overlaps a prior center")

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
    selected = population[:SELECTED_CENTER_COUNT]
    manifest = [
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
    selection = {
        "construction": "all 2P_i and P_i+/-P_j in the published rank-29 subgroup",
        "full_population_count": len(population),
        "full_population_sha256": sha256(
            json.dumps(manifest, separators=(",", ":")).encode()
        ).hexdigest(),
        "prior_public_and_companion_x_overlap_count": len(overlap),
        "selection_rule": (
            "first 128 by (exact projective bit height, x, y, relation); "
            "fixed before modular or square outcomes"
        ),
        "selected_count": len(selected),
        "selected_population_indices": list(range(1, len(selected) + 1)),
        "selected_sha256": sha256(
            json.dumps(manifest[: len(selected)], separators=(",", ":")).encode()
        ).hexdigest(),
        "selected_records": manifest[: len(selected)],
        "prior_companion_source_inventory": source_inventory,
    }
    return population, selected, selection


def exact_prior_center_separation(
    selected: tuple[CompanionCenter, ...],
    old_centers: tuple[CompanionCenter, ...],
    *,
    denominator_min: int,
    offset_radius: int,
) -> dict[str, Any]:
    new_radius = Q(offset_radius, denominator_min**2)
    prior_public_radius = Q(
        PRIOR_PUBLIC_OFFSET_MAX, PRIOR_DENOMINATOR_MIN**2
    )
    prior_companion_radius = Q(
        PRIOR_COMPANION_OFFSET_MAX, PRIOR_DENOMINATOR_MIN**2
    )
    public_witness = min(
        (abs(center.x - point[0]), center_index, public_index)
        for center_index, center in enumerate(selected)
        for public_index, point in enumerate(PUBLISHED_POINTS)
    )
    companion_witness = min(
        (abs(center.x - old.x), center_index, old_index)
        for center_index, center in enumerate(selected)
        for old_index, old in enumerate(old_centers)
    )
    public_pass = public_witness[0] > new_radius + prior_public_radius
    companion_pass = (
        companion_witness[0] > new_radius + prior_companion_radius
    )
    return {
        "new_box_maximum_radius": str(new_radius),
        "prior_public_box_maximum_radius": str(prior_public_radius),
        "prior_companion_box_maximum_radius": str(prior_companion_radius),
        "prior_public_boxes": {
            "passed": public_pass,
            "minimum_center_distance": str(public_witness[0]),
            "selected_center_index": public_witness[1] + 1,
            "public_center_index": public_witness[2] + 1,
            "required_strict_lower_bound": str(new_radius + prior_public_radius),
        },
        "prior_companion_boxes": {
            "passed": companion_pass,
            "minimum_center_distance": str(companion_witness[0]),
            "selected_center_index": companion_witness[1] + 1,
            "prior_companion_center_index": companion_witness[2] + 1,
            "required_strict_lower_bound": str(
                new_radius + prior_companion_radius
            ),
        },
        "all_prior_center_boxes_disjoint": public_pass and companion_pass,
    }


def declared_primitive_count(
    centers: tuple[CompanionCenter, ...],
    denominator_min: int,
    denominator_max: int,
    radius: int,
) -> int:
    return sum(
        2 * positive_coprime_offset_count(radius, denominator)
        for center in centers
        for denominator in range(denominator_min, denominator_max + 1)
        if gcd(denominator, center.denominator_root) == 1
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--denominator-min", type=int, default=50_001)
    parser.add_argument("--denominator-max", type=int, default=100_000)
    parser.add_argument("--offset-radius", type=int, default=16_384)
    parser.add_argument("--wall-cap-seconds", type=float, default=600.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=1000)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=500_000_000)
    parser.add_argument("--progress-every-centers", type=int, default=8)
    args = parser.parse_args()
    if args.denominator_min != 50_001 or args.denominator_max != 100_000:
        raise SystemExit("this lane is pinned to 50001<=b<=100000")
    if args.offset_radius != 16_384:
        raise SystemExit("this lane is pinned to 0<|k|<=16384")
    if args.wall_cap_seconds <= 0 or args.relation_timeout <= 0:
        raise SystemExit("all time caps must be positive")
    if args.stack_bytes < 64_000_000 or args.certificate_prime_bound < 3:
        raise SystemExit("invalid classifier resource bounds")
    if args.progress_every_centers <= 0:
        raise SystemExit("the progress interval must be positive")
    if file_sha256(SOURCE_ENGINE) != EXPECTED_SOURCE_ENGINE_SHA256:
        raise AssertionError("the reusable companion sieve engine changed")
    prior_inventory = [
        {
            "path": str(path.relative_to(ROOT)),
            "sha256": file_sha256(path),
        }
        for path in PRIOR_ARTIFACTS
    ]

    started = time.monotonic()
    population, centers, selection = generate_subgroup_centers()
    old_centers, _ = load_companion_centers()
    nonoverlap = exact_nonoverlap_proof(
        centers,
        denominator_min=args.denominator_min,
        offset_radius=args.offset_radius,
    )
    prior_separation = exact_prior_center_separation(
        centers,
        old_centers,
        denominator_min=args.denominator_min,
        offset_radius=args.offset_radius,
    )
    if not nonoverlap["all_exact_nonoverlap_checks_passed"]:
        raise AssertionError("the direct-chart/new-center nonoverlap gate failed")
    if not prior_separation["all_prior_center_boxes_disjoint"]:
        raise AssertionError("the prior-center-box nonoverlap gate failed")

    expected = declared_primitive_count(
        centers,
        args.denominator_min,
        args.denominator_max,
        args.offset_radius,
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
    completed_center_count = 0
    wall_cap_reached = False
    cache_entry_count = 0

    # Center-major order bounds memory: only one center's residue table is
    # retained at a time, rather than 128 copies of 32769-bit masks.
    for center_index, center in enumerate(centers):
        if time.monotonic() - started >= args.wall_cap_seconds:
            wall_cap_reached = True
            break
        mask_tables: dict[int, tuple[int, ...]] = {}
        for prime in SIEVE_PRIMES:
            table = tuple(
                allowed_offset_mask(center, residue, prime, residue_masks)
                for residue in range(prime)
            )
            mask_tables[prime] = table
            cache_entry_count += len(table)
        for denominator in range(args.denominator_min, args.denominator_max + 1):
            if gcd(denominator, center.denominator_root) != 1:
                continue
            processed += 2 * positive_coprime_offset_count(
                args.offset_radius, denominator
            )
            mask = nonzero_mask
            for prime in SIEVE_PRIMES:
                mask &= mask_tables[prime][denominator % prime]
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
        completed_center_count = center_index + 1
        if (
            completed_center_count % args.progress_every_centers == 0
            or completed_center_count == len(centers)
        ):
            print(
                f"centers {completed_center_count}/{len(centers)}; "
                f"primitive={processed}; survivors={survivors}; "
                f"squares={len(squares)}",
                flush=True,
            )

    companion_lookup = None
    candidate_records: list[dict[str, Any]] = []
    for record in squares:
        point = record.pop("point")
        if companion_lookup is None:
            companion_lookup = small_companion_lookup()
        candidate_records.append(
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
        for record in candidate_records
    )
    complete = not wall_cap_reached and completed_center_count == len(centers)
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
        "artifact_kind": "exact_rank29_subgroup_center_denominator_sieve",
        "status": status,
        "target_hit": target_hit,
        "claim_scope": {
            "bounded": (
                "the completed boxes around the selected 128 exact subgroup "
                "centers only; not a rank upper bound"
            ),
            "one_pass_no_retries": True,
            "selection_was_outcome_blind": True,
        },
        "source": {
            "reusable_engine": str(SOURCE_ENGINE.relative_to(ROOT)),
            "reusable_engine_sha256": EXPECTED_SOURCE_ENGINE_SHA256,
            "prior_artifacts": prior_inventory,
        },
        "center_population": selection,
        "parameters": {
            "denominator_interval": [args.denominator_min, args.denominator_max],
            "nonzero_offset_interval": [
                -args.offset_radius,
                args.offset_radius,
            ],
            "primitivity": "gcd(b,s)=gcd(k,b)=1",
            "sieve_primes": list(SIEVE_PRIMES),
            "wall_cap_seconds": args.wall_cap_seconds,
        },
        "exact_nonoverlap": {
            "new_centers_and_prior_direct_charts": nonoverlap,
            "prior_public_and_companion_center_boxes": prior_separation,
            "all_prior_centered_boxes_covered_by_exact_separation_bounds": True,
        },
        "search_result": {
            "declared_primitive_candidate_count": expected,
            "processed_primitive_candidate_count": processed,
            "completed_center_count": completed_center_count,
            "search_complete": complete,
            "wall_cap_reached": wall_cap_reached,
            "modular_survivor_count_before_primitivity": survivors_before_primitivity,
            "modular_survivor_count_after_primitivity": survivors,
            "modular_survivor_manifest_sha256": survivor_hasher.hexdigest(),
            "negative_homogeneous_value_count_after_sieve": negatives,
            "exact_nonsquare_count_after_sieve": nonsquares,
            "exact_square_abscissa_count": len(squares),
            "candidate_records": candidate_records,
            "certified_independent_30th_point_count": sum(
                record["classification"] == "exact_independent_30th_point"
                for record in candidate_records
            ),
            "rank30_target_hit": target_hit,
            "streamed_residue_cache_entry_count": cache_entry_count,
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
        f"status={status}; centers={completed_center_count}/{len(centers)}; "
        f"processed={processed}/{expected}; survivors={survivors}; "
        f"squares={len(squares)}; target_hit={str(target_hit).lower()}",
        flush=True,
    )


if __name__ == "__main__":
    main()
