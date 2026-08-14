#!/usr/bin/env python3
"""One-pass denominator sieve around exact nonpublic rank-29 companions.

The prior direct and alternate rank-30 artifacts contain no nonpublic points.
The higher-weight-cover artifact contains exact points in the published
rank-29 subgroup.  This script pins those artifacts, replays every recorded
relation over ``QQ``, removes all 29 public abscissas, and deduplicates the
remaining exact x-coordinates.  The resulting companion centers ``u/s^2``
define the boxes

``x = u/s^2 + k/b^2``

with the same denominator-normalized primitivity conditions as the public-
center sieve, but with a disjoint center set.

Disjointness is proved before searching.  Exact center separation keeps these
boxes disjoint from one another and from the public-center boxes.  For every
center and every prior deep x-pair chart, ``Fraction.limit_denominator`` finds
the nearest rational with denominator at most the old height bound; its exact
distance from the center parameter must exceed the entire new perturbation
radius.  The x-offset charts are excluded by an exact absolute-distance bound.

The search itself is one pass with no retries.  Modular quadratic-residue
bitsets eliminate impossible homogeneous square values, every survivor gets
an exact integer square test, and any point outside the cheap known-companion
table is sent immediately to the finite-reduction rank-30 engine before a
time-capped, exactly replayed subgroup-relation proposal.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import gcd, isqrt
from pathlib import Path
import platform
import time
from typing import Any

from elkies_klagsbrun_rank29 import PUBLISHED_POINTS, point_on_general_curve
from pari_bridge import pari_version
from search_elkies_klagsbrun_rank30 import exact_linear_combination, point_negate
from search_elkies_klagsbrun_rank30_denominator_sieve import (
    COEFFICIENT_A,
    COEFFICIENT_B,
    PREVIOUS_DEEP_X_OFFSET_HEIGHT,
    PREVIOUS_DEEP_X_PAIR_HEIGHT,
    SIEVE_PRIMES,
    build_offset_residue_masks,
    classify_point,
    homogeneous_square_value,
    map_square_abscissa,
    positive_coprime_offset_count,
    small_companion_lookup,
)


Q = Fraction
RationalPoint = tuple[Fraction, Fraction]
REPOSITORY = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    REPOSITORY
    / "artifacts/generated-results/"
    "elliptic_elkies_klagsbrun_rank30_companion_center_sieve.json"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/"
    "search_elkies_klagsbrun_rank30_companion_center_sieve.py"
)
SOURCE_ARTIFACT_HASHES = {
    "artifacts/generated-results/elliptic_elkies_klagsbrun_rank30_search.json": (
        "bca9853e5f05ea90d29d2a8eeb2707a3744c73063e834fd0c99cfdd9cbd08005"
    ),
    "artifacts/generated-results/elliptic_elkies_klagsbrun_rank30_search_deep.json": (
        "5f07458bdf57fb84d09abb830f3f21edecd470f0e7d1965c2e9f9a9d21efbe54"
    ),
    "artifacts/generated-results/elliptic_elkies_klagsbrun_rank30_alternate_covers.json": (
        "81d0b502d9cb550b2dcf5b15aa7bdfa2febf16ed9a5cad6af2060e2e294f36cb"
    ),
    "artifacts/generated-results/elliptic_elkies_klagsbrun_rank30_higher_weight_covers.json": (
        "c63af8aa2c6a972c34ecabc9acf0b8f564912a8e8b4956d065b33ab222d769e0"
    ),
}


@dataclass(frozen=True)
class CompanionCenter:
    x: Fraction
    point: RationalPoint
    relation: tuple[int, ...]
    source_paths: tuple[str, ...]

    @property
    def numerator(self) -> int:
        return self.x.numerator

    @property
    def denominator_root(self) -> int:
        root = isqrt(self.x.denominator)
        if root * root != self.x.denominator:
            raise AssertionError("a companion x denominator is not a square")
        return root

    @property
    def bit_height(self) -> int:
        return max(abs(self.x.numerator).bit_length(), self.x.denominator.bit_length())


def source_point_records(payload: dict[str, Any]) -> tuple[dict[str, Any], ...]:
    result = payload["search_result"]
    records = result.get("point_records")
    if records is None:
        records = result.get("candidate_records", [])
    return tuple(records)


def load_companion_centers() -> tuple[tuple[CompanionCenter, ...], list[dict[str, Any]]]:
    """Load, hash-check, and exactly replay the pinned source artifacts."""

    public_x = {point[0] for point in PUBLISHED_POINTS}
    by_x: dict[Fraction, list[tuple[RationalPoint, tuple[int, ...], str]]] = {}
    inventory: list[dict[str, Any]] = []
    for relative_path, expected_hash in SOURCE_ARTIFACT_HASHES.items():
        path = REPOSITORY / relative_path
        raw = path.read_bytes()
        actual_hash = sha256(raw).hexdigest()
        if actual_hash != expected_hash:
            raise AssertionError(f"pinned source artifact changed: {relative_path}")
        payload = json.loads(raw)
        records = source_point_records(payload)
        exact_relation_count = 0
        nonpublic_record_count = 0
        for record in records:
            point = Q(record["x"]), Q(record["y"])
            if not point_on_general_curve(point):
                raise AssertionError("a source artifact point is off the curve")
            relation_raw = record.get("published_basis_relation")
            if relation_raw is None:
                # A certified independent point would already solve the goal;
                # none occurs in the pinned negative artifacts.
                if record.get("classification") == "exact_independent_30th_point":
                    raise AssertionError("a pinned source already contains rank 30")
                continue
            relation = tuple(int(value) for value in relation_raw)
            if len(relation) != len(PUBLISHED_POINTS):
                raise AssertionError("a source relation has the wrong length")
            if exact_linear_combination(relation) != point:
                raise AssertionError("a source subgroup relation failed exact replay")
            exact_relation_count += 1
            if point[0] in public_x:
                continue
            nonpublic_record_count += 1
            by_x.setdefault(point[0], []).append((point, relation, relative_path))
        inventory.append(
            {
                "path": relative_path,
                "sha256": actual_hash,
                "status": payload["status"],
                "point_record_count": len(records),
                "exact_relation_replay_count": exact_relation_count,
                "nonpublic_relation_record_count": nonpublic_record_count,
            }
        )

    centers: list[CompanionCenter] = []
    for x_value, entries in by_x.items():
        entries.sort(key=lambda item: (item[0][1], item[1], item[2]))
        point, relation, _ = entries[0]
        for other_point, _, _ in entries[1:]:
            if other_point != point and other_point != point_negate(point):
                raise AssertionError("equal source abscissas are not inverse points")
        centers.append(
            CompanionCenter(
                x=x_value,
                point=point,
                relation=relation,
                source_paths=tuple(sorted({entry[2] for entry in entries})),
            )
        )
    centers.sort(key=lambda center: (center.bit_height, center.x, center.point[1]))
    if len(centers) != 32:
        raise AssertionError("the pinned nonpublic companion-center count changed")
    return tuple(centers), inventory


def normalized_abscissa(
    center: CompanionCenter, denominator: int, offset: int
) -> tuple[int, int]:
    if denominator <= 0 or offset == 0:
        raise ValueError("the denominator must be positive and the offset nonzero")
    root = center.denominator_root
    if gcd(denominator, root) != 1 or gcd(offset, denominator) != 1:
        raise ValueError("the center, denominator, and offset are not primitive")
    root_denominator = root * denominator
    numerator = center.numerator * denominator**2 + offset * root**2
    if gcd(numerator, root_denominator) != 1:
        raise AssertionError("a companion-centered abscissa is imprimitive")
    return numerator, root_denominator


def allowed_offset_mask(
    center: CompanionCenter,
    denominator: int,
    prime: int,
    residue_masks: dict[int, tuple[int, ...]],
) -> int:
    u_value = center.numerator
    anchor_root = center.denominator_root
    b_mod = denominator % prime
    s_mod = anchor_root % prime
    b2 = b_mod * b_mod % prime
    s2 = s_mod * s_mod % prime
    d_mod = s_mod * b_mod % prime
    d2 = d_mod * d_mod % prime
    d4 = d2 * d2 % prime
    d6 = d4 * d2 % prime
    residues = {value * value % prime for value in range(prime)}
    answer = 0
    for offset_residue in range(prime):
        a_mod = (u_value * b2 + offset_residue * s2) % prime
        value = (
            4 * a_mod * a_mod % prime * a_mod
            + a_mod * a_mod * d2
            + 4 * COEFFICIENT_A * a_mod * d4
            + 4 * COEFFICIENT_B * d6
        ) % prime
        if value in residues:
            answer |= residue_masks[prime][offset_residue]
    return answer


def exact_nonoverlap_proof(
    centers: tuple[CompanionCenter, ...],
    *,
    denominator_min: int,
    offset_radius: int,
) -> dict[str, Any]:
    """Prove the complete boxes miss public boxes and prior direct x charts."""

    public_x = tuple(point[0] for point in PUBLISHED_POINTS)
    radius = Q(offset_radius, denominator_min**2)
    public_witness = min(
        (abs(center.x - x_value), center_index, public_index)
        for center_index, center in enumerate(centers)
        for public_index, x_value in enumerate(public_x)
    )
    center_witness = min(
        (abs(left.x - right.x), left_index, right_index)
        for left_index, left in enumerate(centers)
        for right_index, right in enumerate(centers[left_index + 1 :], left_index + 1)
    )
    public_boxes_disjoint = public_witness[0] > 2 * radius
    companion_boxes_pairwise_disjoint = center_witness[0] > 2 * radius
    x_offset_disjoint = (
        public_witness[0] - radius > PREVIOUS_DEEP_X_OFFSET_HEIGHT
    )

    minimum_ratio: Fraction | None = None
    minimum_record: dict[str, Any] | None = None
    failure_records: list[dict[str, Any]] = []
    pair_count = 0
    for center_index, center in enumerate(centers):
        for left_index, right_index in combinations(range(len(public_x)), 2):
            scale = public_x[right_index] - public_x[left_index]
            center_parameter = (center.x - public_x[left_index]) / scale
            nearest = center_parameter.limit_denominator(
                PREVIOUS_DEEP_X_PAIR_HEIGHT
            )
            distance = abs(center_parameter - nearest)
            perturbation = radius / abs(scale)
            pair_count += 1
            record = {
                "center_index": center_index + 1,
                "public_pair_indices": [left_index + 1, right_index + 1],
                "center_parameter": str(center_parameter),
                "nearest_denominator_bounded_rational": str(nearest),
                "exact_distance": str(distance),
                "maximum_parameter_perturbation": str(perturbation),
            }
            if distance <= perturbation:
                failure_records.append(record)
                continue
            ratio = distance / perturbation
            if minimum_ratio is None or ratio < minimum_ratio:
                minimum_ratio = ratio
                minimum_record = {**record, "exact_clearance_ratio": str(ratio)}
    x_pair_disjoint = not failure_records
    if minimum_ratio is None:
        raise AssertionError("the x-pair nonoverlap audit was empty")
    all_passed = (
        public_boxes_disjoint
        and companion_boxes_pairwise_disjoint
        and x_offset_disjoint
        and x_pair_disjoint
    )
    return {
        "maximum_absolute_x_perturbation": str(radius),
        "public_center_box_separation": {
            "passed": public_boxes_disjoint,
            "minimum_center_distance": str(public_witness[0]),
            "companion_center_index": public_witness[1] + 1,
            "public_center_index": public_witness[2] + 1,
            "required_strict_lower_bound": str(2 * radius),
        },
        "companion_box_pairwise_separation": {
            "passed": companion_boxes_pairwise_disjoint,
            "minimum_center_distance": str(center_witness[0]),
            "center_indices": [center_witness[1] + 1, center_witness[2] + 1],
            "required_strict_lower_bound": str(2 * radius),
        },
        "prior_deep_x_offset_charts": {
            "passed": x_offset_disjoint,
            "chart_parameter_height_bound": PREVIOUS_DEEP_X_OFFSET_HEIGHT,
            "minimum_absolute_parameter_lower_bound": str(
                public_witness[0] - radius
            ),
        },
        "prior_deep_x_pair_charts": {
            "passed": x_pair_disjoint,
            "chart_count_per_center": 406,
            "exact_center_chart_pair_count": pair_count,
            "denominator_superset_bound": PREVIOUS_DEEP_X_PAIR_HEIGHT,
            "failure_count": len(failure_records),
            "failure_records": failure_records,
            "minimum_clearance_record": minimum_record,
        },
        "all_exact_nonoverlap_checks_passed": all_passed,
    }


def declared_primitive_count(
    centers: tuple[CompanionCenter, ...],
    denominator_min: int,
    denominator_max: int,
    radius: int,
) -> int:
    answer = 0
    for denominator in range(denominator_min, denominator_max + 1):
        active = sum(
            gcd(denominator, center.denominator_root) == 1 for center in centers
        )
        answer += active * 2 * positive_coprime_offset_count(radius, denominator)
    return answer


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--denominator-min", type=int, default=3163)
    parser.add_argument("--denominator-max", type=int, default=50_000)
    parser.add_argument("--offset-radius", type=int, default=16_384)
    parser.add_argument("--wall-cap-seconds", type=float, default=240.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=1000)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=500_000_000)
    parser.add_argument("--progress-every", type=int, default=1000)
    args = parser.parse_args()
    if not 1 <= args.denominator_min <= args.denominator_max:
        raise SystemExit("the denominator interval must be positive and ordered")
    if args.offset_radius <= 0:
        raise SystemExit("--offset-radius must be positive")
    if args.wall_cap_seconds <= 0 or args.relation_timeout <= 0:
        raise SystemExit("all time caps must be positive")
    if args.certificate_prime_bound < 3:
        raise SystemExit("--certificate-prime-bound must be at least 3")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes must be at least 64MB")
    if args.progress_every <= 0:
        raise SystemExit("--progress-every must be positive")

    started = time.monotonic()
    centers, source_inventory = load_companion_centers()
    nonoverlap = exact_nonoverlap_proof(
        centers,
        denominator_min=args.denominator_min,
        offset_radius=args.offset_radius,
    )
    if not nonoverlap["all_exact_nonoverlap_checks_passed"]:
        raise SystemExit("the exact nonoverlap gate failed; no search was run")

    declared_count = declared_primitive_count(
        centers,
        args.denominator_min,
        args.denominator_max,
        args.offset_radius,
    )
    residue_masks = build_offset_residue_masks(args.offset_radius)
    nonzero_mask = (1 << (2 * args.offset_radius + 1)) - 1
    nonzero_mask ^= 1 << args.offset_radius
    mask_cache: dict[tuple[int, int, int], int] = {}
    processed_count = 0
    survivors_before_primitivity = 0
    survivor_count = 0
    negative_count = 0
    nonsquare_count = 0
    square_records: list[dict[str, Any]] = []
    survivor_hasher = sha256()
    completed_denominator_max = args.denominator_min - 1
    wall_cap_reached = False

    for denominator in range(args.denominator_min, args.denominator_max + 1):
        if time.monotonic() - started >= args.wall_cap_seconds:
            wall_cap_reached = True
            break
        active = sum(
            gcd(denominator, center.denominator_root) == 1 for center in centers
        )
        processed_count += (
            active
            * 2
            * positive_coprime_offset_count(args.offset_radius, denominator)
        )
        for center_index, center in enumerate(centers):
            if gcd(denominator, center.denominator_root) != 1:
                continue
            mask = nonzero_mask
            for prime in SIEVE_PRIMES:
                key = center_index, prime, denominator % prime
                allowed = mask_cache.get(key)
                if allowed is None:
                    allowed = allowed_offset_mask(
                        center, denominator, prime, residue_masks
                    )
                    mask_cache[key] = allowed
                mask &= allowed
                if mask == 0:
                    break
            survivors_before_primitivity += mask.bit_count()
            while mask:
                low_bit = mask & -mask
                bit_index = low_bit.bit_length() - 1
                mask ^= low_bit
                offset = bit_index - args.offset_radius
                if gcd(offset, denominator) != 1:
                    continue
                survivor_count += 1
                survivor_hasher.update(
                    f"{center_index + 1}|{denominator}|{offset}\n".encode()
                )
                numerator, root_denominator = normalized_abscissa(
                    center, denominator, offset
                )
                value = homogeneous_square_value(numerator, root_denominator)
                if value < 0:
                    negative_count += 1
                    continue
                square_root = isqrt(value)
                if square_root * square_root != value:
                    nonsquare_count += 1
                    continue
                points = map_square_abscissa(
                    numerator, root_denominator, square_root
                )
                square_records.append(
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
        completed_denominator_max = denominator
        if (
            (denominator - args.denominator_min + 1) % args.progress_every == 0
            or denominator == args.denominator_max
        ):
            print(
                f"denominators through {denominator}/{args.denominator_max}; "
                f"primitive={processed_count}; survivors={survivor_count}; "
                f"squares={len(square_records)}",
                flush=True,
            )

    companion_lookup = None
    candidate_records: list[dict[str, Any]] = []
    for square in square_records:
        point = square.pop("point")
        if companion_lookup is None:
            companion_lookup = small_companion_lookup()
        classification = classify_point(
            point,
            companion_lookup=companion_lookup,
            certificate_prime_bound=args.certificate_prime_bound,
            relation_timeout=args.relation_timeout,
            stack_bytes=args.stack_bytes,
        )
        candidate_records.append({**square, **classification})

    target_hit = any(
        record["classification"] == "exact_independent_30th_point"
        for record in candidate_records
    )
    complete = (
        not wall_cap_reached
        and completed_denominator_max == args.denominator_max
    )
    if target_hit:
        status = "exact_rank30_target_hit"
    elif complete:
        status = "bounded_search_no_certified_30th_point"
    else:
        status = "bounded_search_incomplete_at_wall_cap"

    center_manifest = [
        {
            "index": index + 1,
            "x": str(center.x),
            "x_numerator": center.numerator,
            "x_denominator_root": center.denominator_root,
            "representative_y": str(center.point[1]),
            "published_basis_relation": list(center.relation),
            "source_paths": list(center.source_paths),
        }
        for index, center in enumerate(centers)
    ]
    center_digest = sha256(
        json.dumps(center_manifest, separators=(",", ":")).encode()
    ).hexdigest()
    artifact = {
        "schema_version": 1,
        "artifact_kind": "exact_companion_center_denominator_x_sieve",
        "status": status,
        "claim_scope": {
            "exact": (
                "source hashes, source subgroup relations, center extraction, "
                "complete-region nonoverlap, modular exclusions, integer square "
                "tests, and any returned point classification"
            ),
            "bounded": (
                "the completed companion-center denominator boxes only; this is "
                "not a rank upper bound"
            ),
            "one_pass": "one deterministic pass, no retries or adaptive broadening",
        },
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "source_inventory": source_inventory,
        "center_manifest": {
            "count": len(centers),
            "sha256": center_digest,
            "public_x_excluded": True,
            "exact_relation_replay_for_every_center": True,
            "records": center_manifest,
        },
        "parameters": {
            "denominator_interval": [args.denominator_min, args.denominator_max],
            "nonzero_offset_interval": [
                -args.offset_radius,
                args.offset_radius,
            ],
            "primitivity": "gcd(b,s)=gcd(k,b)=1",
            "sieve_primes": list(SIEVE_PRIMES),
            "wall_cap_seconds": args.wall_cap_seconds,
            "certificate_prime_bound": args.certificate_prime_bound,
            "relation_timeout_seconds_each": args.relation_timeout,
            "stack_bytes_each": args.stack_bytes,
        },
        "exact_nonoverlap_proof": nonoverlap,
        "search_result": {
            "one_pass_no_retry": True,
            "declared_primitive_candidate_count": declared_count,
            "processed_primitive_candidate_count": processed_count,
            "completed_denominator_interval": (
                None
                if completed_denominator_max < args.denominator_min
                else [args.denominator_min, completed_denominator_max]
            ),
            "search_complete": complete,
            "wall_cap_reached": wall_cap_reached,
            "modular_survivor_count_before_primitivity": survivors_before_primitivity,
            "modular_survivor_count_after_primitivity": survivor_count,
            "modular_survivor_manifest_sha256": survivor_hasher.hexdigest(),
            "negative_homogeneous_value_count_after_sieve": negative_count,
            "exact_nonsquare_count_after_sieve": nonsquare_count,
            "exact_square_abscissa_count": len(square_records),
            "candidate_records": candidate_records,
            "certified_independent_30th_point_count": sum(
                record["classification"] == "exact_independent_30th_point"
                for record in candidate_records
            ),
            "rank30_target_hit": target_hit,
            "allowed_mask_cache_entry_count": len(mask_cache),
            "wall_seconds": time.monotonic() - started,
        },
    }
    output = args.output
    if not output.is_absolute():
        output = REPOSITORY / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {output}")
    print(
        f"status={status}; processed={processed_count}/{declared_count}; "
        f"survivors={survivor_count}; squares={len(square_records)}; "
        f"target_hit={str(target_hit).lower()}",
        flush=True,
    )


if __name__ == "__main__":
    main()
