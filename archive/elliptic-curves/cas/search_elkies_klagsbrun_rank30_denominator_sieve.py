#!/usr/bin/env python3
"""Direct denominator-normalized x-sieve on the rank-29 record curve.

For the integral generalized Weierstrass model

``y^2 + x*y = x^3 + A*x + B``

put ``z=2*y+x``.  A rational point has ``x=a/d^2`` in lowest terms, and
then

``(d^3*z)^2 = 4*a^3 + a^2*d^2 + 4*A*a*d^4 + 4*B*d^6``.

The earlier searches enumerate bounded-height parameters on x and slope
charts.  This standalone lane instead fixes each of the 29 public abscissas
``u/s^2`` and exhausts the skew boxes

``x = u/s^2 + k/b^2,  b_min <= b <= b_max, 0 < |k| <= K,``

with ``gcd(b,s)=gcd(k,b)=1``.  Thus the new offset is reduced and its
denominator is exactly ``b^2``.  The default ``b_min=3163`` makes that
denominator strictly larger than the previous deep x-offset height bound
``10^7``.

Quadratic-residue conditions modulo a pinned prime list are applied as exact
bitset sieves.  Every surviving primitive pair receives an exact integer
square test and every square maps back to the curve with Fraction arithmetic.
Small, explicit subgroup companions are recognized exactly.  Any other point
is sent first to the existing finite-reduction mod-2 engine and only then to a
time-capped numerical relation proposal whose result is replayed exactly.

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
import time
from typing import Any, Iterable

from elkies_klagsbrun_rank29 import (
    COEFFICIENT_A,
    COEFFICIENT_B,
    PUBLISHED_POINTS,
    point_on_general_curve,
    short_weierstrass_coefficients,
    to_short_point,
)
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)
from pari_bridge import pari_version
from search_elkies_klagsbrun_rank30 import (
    discover_relation,
    point_add,
    point_multiply,
    point_negate,
)


Q = Fraction
RationalPoint = tuple[Fraction, Fraction]
REPOSITORY = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    REPOSITORY
    / "artifacts/generated-results/"
    "elliptic_elkies_klagsbrun_rank30_denominator_sieve.json"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/"
    "search_elkies_klagsbrun_rank30_denominator_sieve.py"
)

# The conditions are valid at bad as well as good reduction primes: an
# integer square remains a square modulo every prime.  Odd primes are used so
# the residue density is useful.
SIEVE_PRIMES = (
    3,
    5,
    7,
    11,
    13,
    17,
    19,
    23,
    29,
    31,
    37,
    41,
    43,
    47,
    53,
    59,
    61,
    67,
    71,
    73,
    79,
    83,
    89,
    97,
)
PREVIOUS_DEEP_X_OFFSET_HEIGHT = 10_000_000
PREVIOUS_DEEP_X_PAIR_HEIGHT = 50_000


def rational_height(value: Fraction) -> int:
    value = Q(value)
    return max(abs(value.numerator), value.denominator)


def anchor_data(index: int) -> tuple[int, int]:
    """Return ``(u,s)`` for the reduced public abscissa ``u/s^2``."""

    x_value = PUBLISHED_POINTS[index][0]
    denominator_root = isqrt(x_value.denominator)
    if denominator_root * denominator_root != x_value.denominator:
        raise AssertionError("a public x denominator is not a square")
    return x_value.numerator, denominator_root


def normalized_abscissa(index: int, denominator: int, offset: int) -> tuple[int, int]:
    """Return coprime ``(a,d)`` with ``x=a/d^2`` in the declared box."""

    if denominator <= 0 or offset == 0:
        raise ValueError("the denominator must be positive and the offset nonzero")
    u_value, anchor_root = anchor_data(index)
    if gcd(denominator, anchor_root) != 1 or gcd(offset, denominator) != 1:
        raise ValueError("the anchor, denominator, and offset are not primitive")
    root_denominator = anchor_root * denominator
    numerator = (
        u_value * denominator * denominator
        + offset * anchor_root * anchor_root
    )
    if gcd(numerator, root_denominator) != 1:
        raise AssertionError("the normalized x pair is unexpectedly imprimitive")
    return numerator, root_denominator


def homogeneous_square_value(numerator: int, root_denominator: int) -> int:
    """Return ``d^6*(2*y+x)^2`` for ``x=numerator/d^2``."""

    a_value = int(numerator)
    d_value = int(root_denominator)
    d2 = d_value * d_value
    d4 = d2 * d2
    d6 = d4 * d2
    return (
        4 * a_value**3
        + a_value * a_value * d2
        + 4 * COEFFICIENT_A * a_value * d4
        + 4 * COEFFICIENT_B * d6
    )


def map_square_abscissa(
    numerator: int, root_denominator: int, square_root: int
) -> tuple[RationalPoint, RationalPoint]:
    """Map an exact homogeneous square to the two inverse curve points."""

    x_value = Q(numerator, root_denominator**2)
    z_value = Q(square_root, root_denominator**3)
    points = (
        (x_value, (-x_value + z_value) / 2),
        (x_value, (-x_value - z_value) / 2),
    )
    if point_negate(points[0]) != points[1]:
        raise AssertionError("the two square-root images are not inverses")
    if not all(point_on_general_curve(point) for point in points):
        raise AssertionError("a homogeneous square mapped off the curve")
    return points


def distinct_prime_factors(value: int) -> tuple[int, ...]:
    answer: list[int] = []
    remaining = int(value)
    divisor = 2
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            answer.append(divisor)
            while remaining % divisor == 0:
                remaining //= divisor
        divisor += 1 if divisor == 2 else 2
    if remaining > 1:
        answer.append(remaining)
    return tuple(answer)


def positive_coprime_offset_count(radius: int, denominator: int) -> int:
    """Count ``1 <= k <= radius`` with ``gcd(k,denominator)=1`` exactly."""

    factors = distinct_prime_factors(denominator)
    answer = 0
    for mask in range(1 << len(factors)):
        product = 1
        parity = 0
        for index, prime in enumerate(factors):
            if (mask >> index) & 1:
                product *= prime
                parity ^= 1
        contribution = radius // product
        answer += -contribution if parity else contribution
    return answer


def build_offset_residue_masks(
    radius: int, primes: Iterable[int] = SIEVE_PRIMES
) -> dict[int, tuple[int, ...]]:
    """Return bit masks selecting each offset residue in ``[-radius,radius]``."""

    length = 2 * radius + 1
    answer: dict[int, tuple[int, ...]] = {}
    for prime in primes:
        masks: list[int] = []
        for residue in range(prime):
            mask = 0
            first_index = (residue + radius) % prime
            for index in range(first_index, length, prime):
                mask |= 1 << index
            masks.append(mask)
        answer[prime] = tuple(masks)
    return answer


def allowed_offset_mask(
    anchor_index: int,
    denominator: int,
    prime: int,
    residue_masks: dict[int, tuple[int, ...]],
) -> int:
    """Return offsets whose homogeneous value is a square modulo ``prime``."""

    u_value, anchor_root = anchor_data(anchor_index)
    b_mod = denominator % prime
    s_mod = anchor_root % prime
    b2 = b_mod * b_mod % prime
    s2 = s_mod * s_mod % prime
    d_mod = s_mod * b_mod % prime
    d2 = d_mod * d_mod % prime
    d4 = d2 * d2 % prime
    d6 = d4 * d2 % prime
    square_residues = {value * value % prime for value in range(prime)}
    answer = 0
    for offset_residue in range(prime):
        a_mod = (u_value * b2 + offset_residue * s2) % prime
        value = (
            4 * a_mod * a_mod % prime * a_mod
            + a_mod * a_mod * d2
            + 4 * COEFFICIENT_A * a_mod * d4
            + 4 * COEFFICIENT_B * d6
        ) % prime
        if value in square_residues:
            answer |= residue_masks[prime][offset_residue]
    return answer


def previous_x_chart_membership(x_value: Fraction) -> tuple[str, int] | None:
    """Return the first prior deep direct-x chart covering ``x_value``."""

    for index, point in enumerate(PUBLISHED_POINTS):
        parameter = x_value - point[0]
        height = rational_height(parameter)
        if height <= PREVIOUS_DEEP_X_OFFSET_HEIGHT:
            return f"xoffset_p{index + 1:02d}", height
    for left, right in combinations(range(len(PUBLISHED_POINTS)), 2):
        center = PUBLISHED_POINTS[left][0]
        scale = PUBLISHED_POINTS[right][0] - center
        parameter = (x_value - center) / scale
        height = rational_height(parameter)
        if height <= PREVIOUS_DEEP_X_PAIR_HEIGHT:
            return f"xpair_p{left + 1:02d}_p{right + 1:02d}", height
    return None


def overlap_calibration(
    denominator_min: int, denominator_max: int, radius: int
) -> dict[str, Any]:
    """Probe a fixed boundary/midpoint grid against the prior direct-x charts."""

    denominator_samples = sorted(
        {
            denominator_min,
            (denominator_min + denominator_max) // 2,
            denominator_max,
        }
    )
    offset_samples = sorted({-radius, -1, 1, radius})
    tested = 0
    overlaps: list[dict[str, Any]] = []
    for denominator in denominator_samples:
        for anchor_index in range(len(PUBLISHED_POINTS)):
            _, anchor_root = anchor_data(anchor_index)
            if gcd(denominator, anchor_root) != 1:
                continue
            for offset in offset_samples:
                if gcd(offset, denominator) != 1:
                    continue
                numerator, root_denominator = normalized_abscissa(
                    anchor_index, denominator, offset
                )
                x_value = Q(numerator, root_denominator**2)
                tested += 1
                membership = previous_x_chart_membership(x_value)
                if membership is not None:
                    identifier, height = membership
                    overlaps.append(
                        {
                            "anchor_index": anchor_index + 1,
                            "denominator": denominator,
                            "offset": offset,
                            "chart_identifier": identifier,
                            "chart_parameter_height": height,
                        }
                    )
    return {
        "denominator_samples": denominator_samples,
        "offset_samples": offset_samples,
        "tested_primitive_sample_count": tested,
        "prior_direct_x_chart_overlap_count": len(overlaps),
        "overlap_fraction": 0 if tested == 0 else len(overlaps) / tested,
        "overlap_records": overlaps,
    }


def relation_vector(index: int, coefficient: int) -> tuple[int, ...]:
    answer = [0] * len(PUBLISHED_POINTS)
    answer[index] = coefficient
    return tuple(answer)


def small_companion_lookup() -> dict[RationalPoint, tuple[int, ...]]:
    """Build exact +/- multiples and signed pair sums of the public points."""

    lookup: dict[RationalPoint, tuple[int, ...]] = {}

    def store(point: RationalPoint, coefficients: tuple[int, ...]) -> None:
        lookup.setdefault(point, coefficients)
        inverse = point_negate(point)
        if inverse is None:
            raise AssertionError("a nonzero companion negated to infinity")
        lookup.setdefault(inverse, tuple(-value for value in coefficients))

    for index, point in enumerate(PUBLISHED_POINTS):
        for scalar in range(1, 5):
            multiple = point_multiply(point, scalar)
            if multiple is None:
                raise AssertionError("a public point has unexpected small torsion")
            store(multiple, relation_vector(index, scalar))
    for left, right in combinations(range(len(PUBLISHED_POINTS)), 2):
        for right_sign in (1, -1):
            right_point = (
                PUBLISHED_POINTS[right]
                if right_sign == 1
                else point_negate(PUBLISHED_POINTS[right])
            )
            companion = point_add(PUBLISHED_POINTS[left], right_point)
            if companion is None:
                raise AssertionError("two distinct public points cancelled")
            coefficients = [0] * len(PUBLISHED_POINTS)
            coefficients[left] = 1
            coefficients[right] = right_sign
            store(companion, tuple(coefficients))
    return lookup


def point_record(point: RationalPoint) -> dict[str, str]:
    return {"x": str(point[0]), "y": str(point[1])}


def classify_point(
    point: RationalPoint,
    *,
    companion_lookup: dict[RationalPoint, tuple[int, ...]],
    certificate_prime_bound: int,
    relation_timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    """Classify exactly, with the rank-30 mod-2 test before relation discovery."""

    base: dict[str, Any] = {
        **point_record(point),
        "exact_curve_membership_checked": True,
    }
    companion = companion_lookup.get(point)
    if companion is not None:
        return {
            **base,
            "classification": "exact_small_companion_in_rank29_subgroup",
            "published_basis_relation": list(companion),
            "exact_fraction_group_law_construction": True,
        }

    augmented = tuple(to_short_point(value) for value in PUBLISHED_POINTS) + (
        to_short_point(point),
    )
    signatures = find_mod2_reduction_certificate(
        short_weierstrass_coefficients(),
        augmented,
        prime_bound=certificate_prime_bound,
    )
    binary_rank = combined_mod2_rank(signatures, len(augmented))
    base.update(
        {
            "augmented_mod2_rank": binary_rank,
            "certificate_primes": [signature.prime for signature in signatures],
            "certificate_prime_bound": certificate_prime_bound,
            "mod2_tested_before_general_relation_search": True,
        }
    )
    if binary_rank == 30:
        return {**base, "classification": "exact_independent_30th_point"}

    relation = discover_relation(
        point, timeout=relation_timeout, stack_bytes=stack_bytes
    )
    if relation is not None:
        return {
            **base,
            "classification": "exactly_in_published_rank29_subgroup",
            "published_basis_relation": list(relation),
            "exact_fraction_group_law_replay": True,
        }
    return {
        **base,
        "classification": "unresolved_after_mod2_and_relation_search",
    }


def declared_primitive_count(
    denominator_min: int, denominator_max: int, radius: int
) -> int:
    answer = 0
    for denominator in range(denominator_min, denominator_max + 1):
        active_anchor_count = sum(
            gcd(denominator, anchor_data(index)[1]) == 1
            for index in range(len(PUBLISHED_POINTS))
        )
        answer += (
            active_anchor_count
            * 2
            * positive_coprime_offset_count(radius, denominator)
        )
    return answer


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--denominator-min", type=int, default=3163)
    parser.add_argument("--denominator-max", type=int, default=50_000)
    parser.add_argument("--offset-radius", type=int, default=16_384)
    parser.add_argument("--wall-cap-seconds", type=float, default=240.0)
    parser.add_argument("--overlap-stop-fraction", type=float, default=0.5)
    parser.add_argument("--certificate-prime-bound", type=int, default=1000)
    parser.add_argument("--relation-timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=500_000_000)
    parser.add_argument("--progress-every", type=int, default=1000)
    args = parser.parse_args()
    if not 1 <= args.denominator_min <= args.denominator_max:
        raise SystemExit("the denominator interval must be positive and ordered")
    if args.offset_radius <= 0:
        raise SystemExit("--offset-radius must be positive")
    if args.denominator_min**2 <= PREVIOUS_DEEP_X_OFFSET_HEIGHT:
        raise SystemExit(
            "--denominator-min must force b^2 above the prior deep x-offset bound"
        )
    if args.wall_cap_seconds <= 0 or args.relation_timeout <= 0:
        raise SystemExit("all time caps must be positive")
    if not 0 <= args.overlap_stop_fraction <= 1:
        raise SystemExit("--overlap-stop-fraction must lie in [0,1]")
    if args.certificate_prime_bound < 3:
        raise SystemExit("--certificate-prime-bound must be at least 3")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes must be at least 64MB")
    if args.progress_every <= 0:
        raise SystemExit("--progress-every must be positive")

    started = time.monotonic()
    calibration = overlap_calibration(
        args.denominator_min, args.denominator_max, args.offset_radius
    )
    stopped_at_overlap_gate = (
        calibration["overlap_fraction"] >= args.overlap_stop_fraction
        and calibration["tested_primitive_sample_count"] > 0
    )
    declared_count = declared_primitive_count(
        args.denominator_min, args.denominator_max, args.offset_radius
    )

    residue_masks = build_offset_residue_masks(args.offset_radius)
    nonzero_mask = (1 << (2 * args.offset_radius + 1)) - 1
    nonzero_mask ^= 1 << args.offset_radius
    mask_cache: dict[tuple[int, int, int], int] = {}
    processed_primitive_count = 0
    modular_survivor_before_primitivity = 0
    modular_survivor_count = 0
    negative_after_sieve = 0
    exact_nonsquare_after_sieve = 0
    square_abscissas: list[dict[str, Any]] = []
    survivor_hasher = sha256()
    completed_denominator_max = args.denominator_min - 1
    wall_cap_reached = False

    if not stopped_at_overlap_gate:
        for denominator in range(args.denominator_min, args.denominator_max + 1):
            if time.monotonic() - started >= args.wall_cap_seconds:
                wall_cap_reached = True
                break
            active_anchor_count = sum(
                gcd(denominator, anchor_data(index)[1]) == 1
                for index in range(len(PUBLISHED_POINTS))
            )
            processed_primitive_count += (
                active_anchor_count
                * 2
                * positive_coprime_offset_count(args.offset_radius, denominator)
            )
            for anchor_index in range(len(PUBLISHED_POINTS)):
                _, anchor_root = anchor_data(anchor_index)
                if gcd(denominator, anchor_root) != 1:
                    continue
                mask = nonzero_mask
                for prime in SIEVE_PRIMES:
                    key = (anchor_index, prime, denominator % prime)
                    allowed = mask_cache.get(key)
                    if allowed is None:
                        allowed = allowed_offset_mask(
                            anchor_index,
                            denominator,
                            prime,
                            residue_masks,
                        )
                        mask_cache[key] = allowed
                    mask &= allowed
                    if mask == 0:
                        break
                modular_survivor_before_primitivity += mask.bit_count()
                while mask:
                    low_bit = mask & -mask
                    bit_index = low_bit.bit_length() - 1
                    mask ^= low_bit
                    offset = bit_index - args.offset_radius
                    if gcd(offset, denominator) != 1:
                        continue
                    modular_survivor_count += 1
                    survivor_hasher.update(
                        f"{anchor_index + 1}|{denominator}|{offset}\n".encode()
                    )
                    numerator, root_denominator = normalized_abscissa(
                        anchor_index, denominator, offset
                    )
                    value = homogeneous_square_value(numerator, root_denominator)
                    if value < 0:
                        negative_after_sieve += 1
                        continue
                    square_root = isqrt(value)
                    if square_root * square_root != value:
                        exact_nonsquare_after_sieve += 1
                        continue
                    points = map_square_abscissa(
                        numerator, root_denominator, square_root
                    )
                    x_value = points[0][0]
                    prior_membership = previous_x_chart_membership(x_value)
                    square_abscissas.append(
                        {
                            "anchor_index": anchor_index + 1,
                            "denominator": denominator,
                            "offset": offset,
                            "normalized_x_numerator": numerator,
                            "normalized_x_denominator_root": root_denominator,
                            "homogeneous_square_root": square_root,
                            "points": points,
                            "prior_direct_x_chart_membership": (
                                None
                                if prior_membership is None
                                else {
                                    "identifier": prior_membership[0],
                                    "parameter_height": prior_membership[1],
                                }
                            ),
                        }
                    )
            completed_denominator_max = denominator
            if (
                (denominator - args.denominator_min + 1) % args.progress_every == 0
                or denominator == args.denominator_max
            ):
                print(
                    f"denominators through {denominator}/{args.denominator_max}; "
                    f"primitive={processed_primitive_count}; "
                    f"sieve_survivors={modular_survivor_count}; "
                    f"squares={len(square_abscissas)}",
                    flush=True,
                )

    companion_lookup: dict[RationalPoint, tuple[int, ...]] | None = None
    candidate_records: list[dict[str, Any]] = []
    seen_x: set[Fraction] = set()
    for square in square_abscissas:
        point = square.pop("points")[0]
        if point[0] in seen_x:
            raise AssertionError("the declared primitive boxes produced duplicate x")
        seen_x.add(point[0])
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
    search_complete = (
        not stopped_at_overlap_gate
        and not wall_cap_reached
        and completed_denominator_max == args.denominator_max
    )
    if target_hit:
        status = "exact_rank30_target_hit"
    elif stopped_at_overlap_gate:
        status = "stopped_at_prior_chart_overlap_gate"
    elif not search_complete:
        status = "bounded_search_incomplete_at_wall_cap"
    else:
        status = "bounded_search_no_certified_30th_point"

    anchor_manifest = [
        {
            "index": index + 1,
            "x_numerator": anchor_data(index)[0],
            "x_denominator_root": anchor_data(index)[1],
        }
        for index in range(len(PUBLISHED_POINTS))
    ]
    artifact = {
        "schema_version": 1,
        "artifact_kind": "exact_denominator_normalized_x_residue_sieve",
        "status": status,
        "claim_scope": {
            "exact": (
                "every primitive pair in the completed declared boxes is either "
                "excluded by a necessary modular square condition or receives an "
                "exact integer square test; all returned points and relations are "
                "checked with exact arithmetic"
            ),
            "bounded": (
                "the completed denominator prefix and pinned offset radius only; "
                "this is not a rank upper bound"
            ),
            "new_region": (
                "each reduced offset denominator is b^2>10^7, outside the prior "
                "deep x-integer-offset parameter-height box"
            ),
        },
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "parameters": {
            "anchor_count": len(PUBLISHED_POINTS),
            "denominator_interval": [
                args.denominator_min,
                args.denominator_max,
            ],
            "nonzero_offset_interval": [
                -args.offset_radius,
                args.offset_radius,
            ],
            "primitivity": "gcd(b,s)=gcd(k,b)=1",
            "sieve_primes": list(SIEVE_PRIMES),
            "wall_cap_seconds": args.wall_cap_seconds,
            "relation_timeout_seconds_each": args.relation_timeout,
            "stack_bytes_each": args.stack_bytes,
            "certificate_prime_bound": args.certificate_prime_bound,
            "previous_deep_x_offset_height": PREVIOUS_DEEP_X_OFFSET_HEIGHT,
            "previous_deep_x_pair_height": PREVIOUS_DEEP_X_PAIR_HEIGHT,
            "overlap_stop_fraction": args.overlap_stop_fraction,
        },
        "anchor_manifest": anchor_manifest,
        "overlap_calibration": calibration,
        "search_result": {
            "declared_primitive_candidate_count": declared_count,
            "processed_primitive_candidate_count": processed_primitive_count,
            "completed_denominator_interval": (
                None
                if completed_denominator_max < args.denominator_min
                else [args.denominator_min, completed_denominator_max]
            ),
            "search_complete": search_complete,
            "wall_cap_reached": wall_cap_reached,
            "stopped_at_overlap_gate": stopped_at_overlap_gate,
            "modular_survivor_count_before_primitivity": (
                modular_survivor_before_primitivity
            ),
            "modular_survivor_count_after_primitivity": modular_survivor_count,
            "modular_survivor_manifest_sha256": survivor_hasher.hexdigest(),
            "negative_homogeneous_value_count_after_sieve": negative_after_sieve,
            "exact_nonsquare_count_after_sieve": exact_nonsquare_after_sieve,
            "exact_square_abscissa_count": len(square_abscissas),
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
        f"status={status}; processed={processed_primitive_count}/{declared_count}; "
        f"survivors={modular_survivor_count}; squares={len(square_abscissas)}; "
        f"target_hit={str(target_hit).lower()}",
        flush=True,
    )


if __name__ == "__main__":
    main()
