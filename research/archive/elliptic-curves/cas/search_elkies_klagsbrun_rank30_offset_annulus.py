#!/usr/bin/env python3
"""Exact offset-annulus x-sieve on the rank-29 record curve.

For each of the 29 public abscissas ``x_i=u/s^2``, exhaust

``x=x_i+k/b^2, 3163<=b<=50000, 16385<=|k|<=65536,``

under ``gcd(b,s)=gcd(k,b)=1``.  The denominator range agrees with the first
denominator-normalized sieve, but the offset annulus is disjoint from its
``0<|k|<=16384`` box.  It is also disjoint from the later ``b>=50001``
extension.  Public centers are separated by more than twice the maximum
possible shift, checked exactly before the search.

The same exact quadratic-residue masks and terminal integer-square tests are
used.  Any square is mapped back to the curve and tested against the certified
rank-29 subgroup.  This is a finite negative-search certificate for the
declared annulus, not a rank upper bound.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt
import json
from pathlib import Path
import platform
import shlex
import sys
import time
from typing import Any

from elkies_klagsbrun_rank29 import PUBLISHED_POINTS
from pari_bridge import pari_version
from search_elkies_klagsbrun_rank30_denominator_sieve import (
    SIEVE_PRIMES,
    allowed_offset_mask,
    anchor_data,
    build_offset_residue_masks,
    classify_point,
    homogeneous_square_value,
    map_square_abscissa,
    normalized_abscissa,
    positive_coprime_offset_count,
    previous_x_chart_membership,
    small_companion_lookup,
)


Q = Fraction
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_elkies_klagsbrun_rank30_offset_annulus.json"
)
SOURCE_ENGINE = (
    ROOT
    / "elliptic-curves"
    / "cas"
    / "search_elkies_klagsbrun_rank30_denominator_sieve.py"
)
EXPECTED_SOURCE_ENGINE_SHA256 = (
    "73d46aa761e912feee18667a7ffea240db9a211ce9401ec73b8722fe54f09878"
)
PRIOR_BASE_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_elkies_klagsbrun_rank30_denominator_sieve.json"
)
EXPECTED_PRIOR_BASE_SHA256 = (
    "6f1d16c429b8d2df330eaa70e37161cc77e41a6bd1e641d0f2a90133406b70a9"
)
PRIOR_EXTENSION_ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_elkies_klagsbrun_rank30_denominator_sieve_b500000.json"
)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def annulus_mask(minimum: int, maximum: int) -> int:
    """Return the bit mask for ``minimum<=|k|<=maximum`` at index k+maximum."""

    if not 1 <= minimum <= maximum:
        raise ValueError("the annulus bounds must satisfy 1<=minimum<=maximum")
    width = maximum - minimum + 1
    negative = (1 << width) - 1
    positive = negative << (maximum + minimum)
    return negative | positive


def positive_coprime_annulus_count(minimum: int, maximum: int, denominator: int) -> int:
    return positive_coprime_offset_count(maximum, denominator) - positive_coprime_offset_count(
        minimum - 1, denominator
    )


def exact_center_separation(denominator_min: int, offset_max: int) -> dict[str, Any]:
    ordered = sorted((point[0], index + 1) for index, point in enumerate(PUBLISHED_POINTS))
    gap, left, right = min(
        (ordered[index + 1][0] - ordered[index][0], ordered[index][1], ordered[index + 1][1])
        for index in range(len(ordered) - 1)
    )
    twice_shift = Q(2 * offset_max, denominator_min**2)
    if gap <= twice_shift:
        raise AssertionError("public-center annuli are not pairwise separated")
    return {
        "minimum_public_center_gap": str(gap),
        "minimum_gap_center_indices": [left, right],
        "twice_maximum_annulus_shift": str(twice_shift),
        "pairwise_center_boxes_disjoint": True,
    }


def declared_count(
    denominator_min: int,
    denominator_max: int,
    offset_min: int,
    offset_max: int,
) -> int:
    answer = 0
    for denominator in range(denominator_min, denominator_max + 1):
        active = sum(
            gcd(denominator, anchor_data(index)[1]) == 1
            for index in range(len(PUBLISHED_POINTS))
        )
        answer += 2 * active * positive_coprime_annulus_count(
            offset_min, offset_max, denominator
        )
    return answer


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--denominator-min", type=int, default=3163)
    parser.add_argument("--denominator-max", type=int, default=50_000)
    parser.add_argument("--offset-min", type=int, default=16_385)
    parser.add_argument("--offset-max", type=int, default=65_536)
    parser.add_argument("--wall-cap-seconds", type=float, default=600.0)
    parser.add_argument("--progress-every", type=int, default=2500)
    parser.add_argument("--certificate-prime-bound", type=int, default=1000)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=500_000_000)
    args = parser.parse_args()
    if not 1 <= args.denominator_min <= args.denominator_max:
        raise SystemExit("the denominator interval must be positive and ordered")
    if not 1 <= args.offset_min <= args.offset_max:
        raise SystemExit("the offset annulus must be positive and ordered")
    if args.offset_min <= 16_384:
        raise SystemExit("the annulus must start outside the completed base box")
    if args.denominator_min != 3163 or args.denominator_max > 50_000:
        raise SystemExit("this lane is pinned to the completed base denominator range")
    if args.wall_cap_seconds <= 0 or args.progress_every <= 0:
        raise SystemExit("caps and progress interval must be positive")
    if file_sha256(SOURCE_ENGINE) != EXPECTED_SOURCE_ENGINE_SHA256:
        raise AssertionError("the source sieve engine changed")
    if file_sha256(PRIOR_BASE_ARTIFACT) != EXPECTED_PRIOR_BASE_SHA256:
        raise AssertionError("the completed base-box artifact changed")
    extension = json.loads(PRIOR_EXTENSION_ARTIFACT.read_text(encoding="utf-8"))
    if extension["parameters"]["denominator_interval"] != [50_001, 500_000]:
        raise AssertionError("the denominator extension boundary changed")

    started = time.monotonic()
    separation = exact_center_separation(args.denominator_min, args.offset_max)
    expected = declared_count(
        args.denominator_min,
        args.denominator_max,
        args.offset_min,
        args.offset_max,
    )
    residue_masks = build_offset_residue_masks(args.offset_max)
    active_annulus = annulus_mask(args.offset_min, args.offset_max)
    mask_cache: dict[tuple[int, int, int], int] = {}
    processed = 0
    survivors_before_primitivity = 0
    survivors = 0
    negative_values = 0
    nonsquares = 0
    squares: list[dict[str, Any]] = []
    survivor_hasher = sha256()
    completed_max = args.denominator_min - 1
    wall_cap_reached = False

    for denominator in range(args.denominator_min, args.denominator_max + 1):
        if time.monotonic() - started >= args.wall_cap_seconds:
            wall_cap_reached = True
            break
        active_anchor_count = sum(
            gcd(denominator, anchor_data(index)[1]) == 1
            for index in range(len(PUBLISHED_POINTS))
        )
        processed += 2 * active_anchor_count * positive_coprime_annulus_count(
            args.offset_min, args.offset_max, denominator
        )
        for anchor_index in range(len(PUBLISHED_POINTS)):
            _, anchor_root = anchor_data(anchor_index)
            if gcd(denominator, anchor_root) != 1:
                continue
            mask = active_annulus
            for prime in SIEVE_PRIMES:
                key = anchor_index, prime, denominator % prime
                allowed = mask_cache.get(key)
                if allowed is None:
                    allowed = allowed_offset_mask(
                        anchor_index, denominator, prime, residue_masks
                    )
                    mask_cache[key] = allowed
                mask &= allowed
                if not mask:
                    break
            survivors_before_primitivity += mask.bit_count()
            while mask:
                low = mask & -mask
                bit_index = low.bit_length() - 1
                mask ^= low
                offset = bit_index - args.offset_max
                if gcd(offset, denominator) != 1:
                    continue
                survivors += 1
                survivor_hasher.update(
                    f"{anchor_index + 1}|{denominator}|{offset}\n".encode()
                )
                numerator, root_denominator = normalized_abscissa(
                    anchor_index, denominator, offset
                )
                value = homogeneous_square_value(numerator, root_denominator)
                if value < 0:
                    negative_values += 1
                    continue
                square_root = isqrt(value)
                if square_root * square_root != value:
                    nonsquares += 1
                    continue
                points = map_square_abscissa(numerator, root_denominator, square_root)
                prior = previous_x_chart_membership(points[0][0])
                squares.append(
                    {
                        "anchor_index": anchor_index + 1,
                        "denominator": denominator,
                        "offset": offset,
                        "normalized_x_numerator": numerator,
                        "normalized_x_denominator_root": root_denominator,
                        "homogeneous_square_root": square_root,
                        "point": points[0],
                        "prior_direct_x_chart_membership": prior,
                    }
                )
        completed_max = denominator
        if (
            (denominator - args.denominator_min + 1) % args.progress_every == 0
            or denominator == args.denominator_max
        ):
            print(
                f"denominators through {denominator}/{args.denominator_max}; "
                f"primitive={processed}; survivors={survivors}; squares={len(squares)}",
                flush=True,
            )

    companion_lookup = None
    candidate_records = []
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
    complete = not wall_cap_reached and completed_max == args.denominator_max
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
        "artifact_kind": "exact_public_center_offset_annulus_sieve",
        "status": status,
        "target_hit": target_hit,
        "claim_scope": {
            "bounded": "the completed public-center offset annulus only",
            "not_a_rank_upper_bound": True,
            "one_pass_no_retries": True,
        },
        "source": {
            "engine": str(SOURCE_ENGINE.relative_to(ROOT)),
            "engine_sha256": EXPECTED_SOURCE_ENGINE_SHA256,
            "prior_base_artifact": str(PRIOR_BASE_ARTIFACT.relative_to(ROOT)),
            "prior_base_artifact_sha256": EXPECTED_PRIOR_BASE_SHA256,
            "prior_extension_artifact": str(PRIOR_EXTENSION_ARTIFACT.relative_to(ROOT)),
            "prior_extension_artifact_sha256": file_sha256(PRIOR_EXTENSION_ARTIFACT),
        },
        "parameters": {
            "anchor_count": len(PUBLISHED_POINTS),
            "denominator_interval": [args.denominator_min, args.denominator_max],
            "absolute_offset_interval": [args.offset_min, args.offset_max],
            "primitivity": "gcd(b,s)=gcd(k,b)=1",
            "sieve_primes": list(SIEVE_PRIMES),
            "wall_cap_seconds": args.wall_cap_seconds,
        },
        "nonoverlap": {
            **separation,
            "disjoint_from_prior_base_by_offset_magnitude": True,
            "disjoint_from_prior_extension_by_denominator": True,
            "reduced_offset_denominator_is_exactly_b_squared": True,
        },
        "search_result": {
            "declared_primitive_candidate_count": expected,
            "processed_primitive_candidate_count": processed,
            "completed_denominator_interval": (
                None
                if completed_max < args.denominator_min
                else [args.denominator_min, completed_max]
            ),
            "search_complete": complete,
            "wall_cap_reached": wall_cap_reached,
            "modular_survivor_count_before_primitivity": survivors_before_primitivity,
            "modular_survivor_count_after_primitivity": survivors,
            "modular_survivor_manifest_sha256": survivor_hasher.hexdigest(),
            "negative_homogeneous_value_count_after_sieve": negative_values,
            "exact_nonsquare_count_after_sieve": nonsquares,
            "exact_square_abscissa_count": len(squares),
            "candidate_records": candidate_records,
            "certified_independent_30th_point_count": sum(
                record["classification"] == "exact_independent_30th_point"
                for record in candidate_records
            ),
            "rank30_target_hit": target_hit,
            "allowed_mask_cache_entry_count": len(mask_cache),
            "wall_seconds": time.monotonic() - started,
        },
        "reproduction": {
            "command": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
            "script": str(script.relative_to(ROOT)),
            "script_sha256": file_sha256(script),
            "python": platform.python_version(),
            "pari_gp": pari_version(),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")
    print(
        f"status={status}; processed={processed}/{expected}; survivors={survivors}; "
        f"squares={len(squares)}; target_hit={str(target_hit).lower()}",
        flush=True,
    )


if __name__ == "__main__":
    main()
