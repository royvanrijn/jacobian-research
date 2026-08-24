#!/usr/bin/env python3
"""Compressed p-adic CRT search for the Mestre tuple (0,4,30,31,39,46).

This lane is intentionally independent of ``scan_mestre_0430313946.cpp`` and
of the global root-tuple driver.  It derives and checks the exact primitive
discriminant geometry, exhausts all roots modulo ``p^4`` for every prime at
most 199, and compresses complete sibling sets into maximal p-adic balls.
Five automatically discovered clean multiplicative unions are then combined
by exact CRT.  A declared box in each exact Gauss-reduced two-dimensional
lattice is exhausted, with the active global scan box removed exactly up to
the even ``T -> -T`` symmetry.

The bounded population is ranked without conductor leakage by a small-prime
radical upper proxy and two held good-prime trace bands.  A predeclared union
of those rankings receives one capped exact conductor computation each.  All
sub-target fibers receive capped quartic point searches at heights 50k, 250k,
and 1m.  Numerical height ranks are triage only; an exact finite-reduction
independence attempt is triggered immediately if a stable rank reaches 21.
Every subprocess is foreground, single-attempt, and process-group capped.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import gcd, log
from pathlib import Path
import platform
import shlex
import sys
from typing import Any, Iterable, Sequence

from crt_lattice import crt_pair, gauss_reduce
from ek_k3 import primes_up_to, rational_to_string
from mestre_root_tuples import SixRootMestreConstruction
from multiple_root_lifting import (
    RootBall,
    affine_variable_coefficients,
    all_roots_mod_prime_power,
    fixed_divisor_valuation,
)
from search_mestre_root_tuple_scale import (
    CappedProcessTimeout,
    bounded_quartic_points,
    canonical_signless_points,
    capped_minimal_curve_data,
    finite_reduction_attempt,
    height_matrix_replay,
    numerical_subset,
    pari_version_capped,
    point_digest,
    primitive_visible_points,
    quartic_point_to_jacobian,
)


Q = Fraction
ROOTS = (0, 4, 30, 31, 39, 46)
TARGET_LOG_CONDUCTOR = Decimal("182.72")
LOCAL_PRIME_BOUND = 199
LOCAL_TARGET_EXPONENT = 4
LOCAL_ROOT_CAP = 100_000
SELECTED_POWER_PRIMES = (5, 11, 13, 37, 43)
HELD_TRACE_BANDS = (
    (101, 103, 107, 109, 113, 127, 131, 137, 139, 149),
    (151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199),
)
ACTIVE_NUMERATOR_BOUND = 100_000
ACTIVE_DENOMINATOR_BOUND = 5_000
DEFAULT_COEFFICIENT_RADIUS = 24
DEFAULT_CONDUCTOR_KEEP = 20
DEFAULT_SUBTARGET_TRIAGE_CAP = 4
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/elliptic_mestre_0430313946_power_crt.json"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_mestre_0430313946_power_crt.py"
)

# Ascending coefficients of the primitive short Jacobian y^2=x^3+A(T)x+B(T).
# The script verifies these against the general exact construction at enough
# parameters to determine polynomials of the declared degrees uniquely.
A_COEFFICIENTS = (
    -4_840_084_292_061_123,
    0,
    -611_317_954_173_024,
    0,
    1_593_386_668_512,
    0,
    178_536_960,
    0,
    -4_762_800,
)
B_COEFFICIENTS = (
    -84_220_770_768_445_491_421_122,
    0,
    56_408_643_410_668_151_386_896,
    0,
    -193_116_296_662_521_546_636,
    0,
    930_870_844_049_215_872,
    0,
    -2_013_991_695_596_160,
    0,
    -224_956_569_600,
    0,
    4_000_752_000,
)

# The irreducible-over-Q degree-16 factor of the coefficient-primitive
# discriminant core, in ascending order.  Odd coefficients vanish.
Q16_COEFFICIENTS = (
    4_790_673_696_756_029_110_247_045_369_088,
    0,
    7_919_240_759_869_795_553_035_334_023_164,
    0,
    -1_054_875_902_344_473_896_117_429_109_363,
    0,
    5_682_498_165_802_947_031_258_092_192,
    0,
    -8_703_761_868_038_350_177_904_133,
    0,
    -3_102_186_834_817_565_916_456,
    0,
    15_786_627_460_019_717_808,
    0,
    -5_837_288_789_113_600,
    0,
    404_860_108_800,
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def integer_digest(items: Iterable[Sequence[int]]) -> str:
    text = "\n".join(",".join(str(value) for value in item) for item in items)
    return hashlib.sha256(text.encode()).hexdigest()


def polynomial_value(coefficients: Sequence[int], value: int, modulus: int | None = None) -> int:
    answer = 0
    for coefficient in reversed(coefficients):
        answer = answer * value + coefficient
        if modulus is not None:
            answer %= modulus
    return answer


def polynomial_add(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    answer = [0] * max(len(left), len(right))
    for index in range(len(answer)):
        answer[index] = (left[index] if index < len(left) else 0) + (
            right[index] if index < len(right) else 0
        )
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def polynomial_scale(coefficients: Sequence[int], scale: int) -> tuple[int, ...]:
    return tuple(scale * coefficient for coefficient in coefficients)


def polynomial_multiply(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    answer = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] += left_value * right_value
    return tuple(answer)


def polynomial_power(coefficients: Sequence[int], exponent: int) -> tuple[int, ...]:
    answer = (1,)
    factor = tuple(coefficients)
    while exponent:
        if exponent & 1:
            answer = polynomial_multiply(answer, factor)
        factor = polynomial_multiply(factor, factor)
        exponent //= 2
    return answer


def homogenized_value(coefficients: Sequence[int], numerator: int, denominator: int) -> int:
    degree = len(coefficients) - 1
    return sum(
        coefficient * numerator**index * denominator ** (degree - index)
        for index, coefficient in enumerate(coefficients)
    )


def integer_valuation(value: int, prime: int) -> int:
    if value == 0:
        raise ValueError("the valuation of zero is not finite")
    value = abs(value)
    answer = 0
    while value % prime == 0:
        value //= prime
        answer += 1
    return answer


def derive_exact_geometry(
    construction: SixRootMestreConstruction,
) -> tuple[tuple[int, ...], dict[str, Any]]:
    discriminant = construction.primitive_discriminant_polynomial
    if any(value.denominator != 1 for value in discriminant):
        raise AssertionError("the primitive quartic discriminant ceased to be integral")
    integral = tuple(value.numerator for value in discriminant)
    coefficient_content = gcd(*map(abs, integral))
    core = tuple(value // coefficient_content for value in integral)
    if gcd(*map(abs, core)) != 1 or len(core) != 21:
        raise AssertionError("the discriminant core lost primitivity or degree 20")

    collision_factor = polynomial_multiply(
        polynomial_power((-21, 2), 2), polynomial_power((21, 2), 2)
    )
    expected_core = polynomial_multiply(collision_factor, Q16_COEFFICIENTS)
    if core != expected_core:
        raise AssertionError("the exact discriminant factorization changed")

    # Independent short-Weierstrass replay.  The 9 (respectively 13) exact
    # values determine A (respectively B) within their proven degree bounds.
    for parameter in range(1, 14):
        coefficients = construction.primitive_jacobian_coefficients(Q(parameter))
        if coefficients[:3] != (0, 0, 0):
            raise AssertionError("the primitive Jacobian is no longer short")
        if parameter <= 9 and coefficients[3] != polynomial_value(A_COEFFICIENTS, parameter):
            raise AssertionError("the pinned degree-eight A polynomial changed")
        if coefficients[4] != polynomial_value(B_COEFFICIENTS, parameter):
            raise AssertionError("the pinned degree-twelve B polynomial changed")
    short_delta = polynomial_scale(
        polynomial_add(
            polynomial_scale(polynomial_power(A_COEFFICIENTS, 3), 4),
            polynomial_scale(polynomial_power(B_COEFFICIENTS, 2), 27),
        ),
        -16,
    )
    short_delta_constant = 4_499_817_235_200
    if short_delta != polynomial_scale(core, short_delta_constant):
        raise AssertionError("the short discriminant/core identity changed")

    # Prove the homogeneous fixed divisor is exactly 28.  Exhaustion of all
    # primitive pairs modulo 4 and modulo 7 proves the lower divisibility.
    lower_certificates = {}
    for modulus in (4, 7):
        primitive_pairs = [
            (a_value, b_value)
            for a_value in range(modulus)
            for b_value in range(modulus)
            if gcd(gcd(a_value, b_value), modulus) == 1
        ]
        if not all(homogenized_value(core, *pair) % modulus == 0 for pair in primitive_pairs):
            raise AssertionError("the homogeneous fixed-divisor lower bound failed")
        lower_certificates[str(modulus)] = len(primitive_pairs)
    sample_pairs = ((0, 1), (1, 0), (1, 1), (1, 2), (2, 1))
    sample_gcd = gcd(
        *(abs(homogenized_value(core, *pair)) for pair in sample_pairs)
    )
    if sample_gcd != 28:
        raise AssertionError("the homogeneous fixed-divisor upper certificate changed")

    return core, {
        "status": "exact polynomial identities",
        "roots": list(ROOTS),
        "quartic_condition": str(construction.quartic_condition),
        "primitive_quartic_content": str(construction.quartic_content),
        "primitive_quartic_square_scale": str(construction.quartic_square_scale),
        "primitive_quartic_discriminant_degree": len(integral) - 1,
        "primitive_quartic_discriminant_coefficient_content": coefficient_content,
        "coefficient_primitive_discriminant_core_ascending": list(core),
        "factorization": {
            "text": "(2*T-21)^2*(2*T+21)^2*Q16(T)",
            "linear_collision_parameters": ["-21/2", "21/2"],
            "linear_factor_multiplicities": [2, 2],
            "Q16_ascending": list(Q16_COEFFICIENTS),
            "factorization_checked_by_exact_multiplication": True,
        },
        "short_jacobian": {
            "A_ascending": list(A_COEFFICIENTS),
            "B_ascending": list(B_COEFFICIENTS),
            "weierstrass_discriminant_core_multiplier": short_delta_constant,
            "core_multiplier_factorization": "2^8*3^15*5^2*7^2",
            "exact_polynomial_identity_checked": True,
        },
        "primitive_pair_homogeneous_fixed_divisor": {
            "value": 28,
            "factorization": "2^2*7",
            "lower_bound_exhaustive_residue_certificates": lower_certificates,
            "upper_bound_sample_pairs": [list(pair) for pair in sample_pairs],
            "gcd_of_upper_bound_sample_values": sample_gcd,
            "interpretation": (
                "every coprime (a,b) has 28 dividing the coefficient-primitive "
                "homogeneous core; the displayed finite sample proves no larger "
                "universal divisor"
            ),
        },
    }


def character_trace(coefficient_a: int, coefficient_b: int, prime: int) -> int:
    character_sum = 0
    for x_value in range(prime):
        rhs = (x_value**3 + coefficient_a * x_value + coefficient_b) % prime
        if rhs:
            character_sum += 1 if pow(rhs, (prime - 1) // 2, prime) == 1 else -1
    return -character_sum


def reduction_at_residue(residue: int, prime: int) -> dict[str, Any]:
    coefficient_a = polynomial_value(A_COEFFICIENTS, residue, prime)
    coefficient_b = polynomial_value(B_COEFFICIENTS, residue, prime)
    discriminant_core = (4 * coefficient_a**3 + 27 * coefficient_b**2) % prime
    if discriminant_core:
        trace = character_trace(coefficient_a, coefficient_b, prime)
        return {"kind": "good", "trace": trace, "point_count": prime + 1 - trace}
    if prime < 5 or coefficient_a == 0:
        return {
            "kind": "not-certified-clean-multiplicative",
            "c4_unit": False,
        }
    trace = character_trace(coefficient_a, coefficient_b, prime)
    if trace not in (-1, 1):
        raise AssertionError("a unit-c4 singular cubic did not have multiplicative trace")
    return {
        "kind": "split multiplicative" if trace == 1 else "nonsplit multiplicative",
        "c4_unit": True,
        "bad_euler_trace": trace,
    }


def ball_record(core: Sequence[int], ball: RootBall) -> dict[str, Any]:
    forced = fixed_divisor_valuation(
        affine_variable_coefficients(core, ball.residue, ball.modulus), ball.prime
    )
    reduction = reduction_at_residue(ball.residue % ball.prime, ball.prime)
    return {
        "residue": ball.residue,
        "modulus_exponent": ball.exponent,
        "modulus": ball.modulus,
        "forced_core_valuation": forced,
        "compression_gain_exponent": forced - ball.exponent,
        "reduction": reduction,
    }


def discover_local_profiles(
    core: Sequence[int],
) -> tuple[list[dict[str, Any]], dict[int, tuple[dict[str, Any], ...]]]:
    profiles: list[dict[str, Any]] = []
    selected: dict[int, tuple[dict[str, Any], ...]] = {}
    for prime in primes_up_to(LOCAL_PRIME_BOUND):
        result = all_roots_mod_prime_power(
            core,
            prime,
            LOCAL_TARGET_EXPONENT,
            max_roots=LOCAL_ROOT_CAP,
        )
        balls = tuple(ball_record(core, ball) for ball in result.maximal_balls())
        efficient_clean = tuple(
            ball
            for ball in balls
            if ball["modulus_exponent"] < LOCAL_TARGET_EXPONENT
            and ball["forced_core_valuation"] >= LOCAL_TARGET_EXPONENT
            and ball["reduction"].get("c4_unit") is True
        )
        profiles.append(
            {
                "prime": prime,
                "root_counts_mod_p_through_p4": list(result.level_counts),
                "roots_mod_p4": len(result.roots),
                "candidate_digits_checked": result.candidate_digits_checked,
                "maximal_balls": list(balls),
                "efficient_clean_multiplicative_balls": len(efficient_clean),
            }
        )
        if prime in SELECTED_POWER_PRIMES:
            if not efficient_clean:
                raise AssertionError(f"selected prime {prime} lost every efficient clean ball")
            selected[prime] = efficient_clean
    if tuple(sorted(selected)) != SELECTED_POWER_PRIMES:
        raise AssertionError("a selected power prime was not profiled")
    return profiles, selected


def combine_selected_balls(
    groups: dict[int, tuple[dict[str, Any], ...]],
) -> tuple[dict[str, Any], ...]:
    classes: list[dict[str, Any]] = []
    ordered_groups = [groups[prime] for prime in SELECTED_POWER_PRIMES]
    for choices in product(*ordered_groups):
        residue, modulus = 0, 1
        for choice in choices:
            residue, modulus = crt_pair(
                residue, modulus, choice["residue"], choice["modulus"]
            )
        classes.append(
            {
                "residue": residue,
                "modulus": modulus,
                "choices": tuple(choices),
            }
        )
    keys = [(item["residue"], item["modulus"]) for item in classes]
    if len(keys) != len(set(keys)):
        raise AssertionError("distinct local ball choices collapsed under CRT")
    return tuple(classes)


def matching_choice(
    numerator: int,
    denominator: int,
    group: Sequence[dict[str, Any]],
) -> dict[str, Any] | None:
    matches = [
        choice
        for choice in group
        if (numerator - choice["residue"] * denominator) % choice["modulus"] == 0
    ]
    if len(matches) > 1:
        raise AssertionError("selected p-adic balls overlap")
    return matches[0] if matches else None


def enumerate_lattice_population(
    classes: Sequence[dict[str, Any]],
    groups: dict[int, tuple[dict[str, Any], ...]],
    *,
    coefficient_radius: int,
) -> tuple[list[tuple[int, int]], dict[str, Any]]:
    representatives: dict[tuple[int, int], set[int]] = {}
    bounded_vectors = 0
    rejected_zero_or_infinite = 0
    for class_index, crt_class in enumerate(classes):
        residue = crt_class["residue"]
        modulus = crt_class["modulus"]
        first, second = gauss_reduce((modulus, 0), (residue, 1))
        if abs(first[0] * second[1] - first[1] * second[0]) != modulus:
            raise AssertionError("Gauss reduction changed the CRT lattice determinant")
        for left in range(-coefficient_radius, coefficient_radius + 1):
            for right in range(-coefficient_radius, coefficient_radius + 1):
                if left == 0 and right == 0:
                    continue
                bounded_vectors += 1
                numerator = left * first[0] + right * second[0]
                denominator = left * first[1] + right * second[1]
                if denominator == 0 or gcd(denominator, modulus) != 1:
                    rejected_zero_or_infinite += 1
                    continue
                common = gcd(abs(numerator), abs(denominator))
                numerator //= common
                denominator //= common
                if denominator < 0:
                    numerator = -numerator
                    denominator = -denominator
                # The family, discriminant, and selected ball unions are all
                # stable under T -> -T, so retain one exact sign quotient.
                numerator = abs(numerator)
                if numerator == 0 or gcd(numerator, denominator) != 1:
                    rejected_zero_or_infinite += 1
                    continue
                if any(
                    matching_choice(numerator, denominator, groups[prime]) is None
                    for prime in SELECTED_POWER_PRIMES
                ):
                    raise AssertionError("a reduced representative lost its p-adic union")
                representatives.setdefault((numerator, denominator), set()).add(class_index)

    ordered = sorted(
        representatives,
        key=lambda item: (max(item[0], item[1]), item[0], item[1]),
    )
    active_overlap = [
        item
        for item in ordered
        if item[0] <= ACTIVE_NUMERATOR_BOUND and item[1] <= ACTIVE_DENOMINATOR_BOUND
    ]
    calibration_present = (5, 1) in representatives
    outside = [item for item in ordered if item not in set(active_overlap) and item != (5, 1)]
    return outside, {
        "gauss_coefficient_radius": coefficient_radius,
        "bounded_coefficient_vectors_visited": bounded_vectors,
        "zero_infinite_or_nonunit_denominator_rejections": rejected_zero_or_infinite,
        "unique_sign_quotiented_representatives": len(ordered),
        "active_global_box": {
            "canonical_absolute_numerator_bound": ACTIVE_NUMERATOR_BOUND,
            "denominator_bound": ACTIVE_DENOMINATOR_BOUND,
            "symmetry": "T -> -T",
            "exact_overlap_count": len(active_overlap),
            "exact_overlap_parameters": [
                rational_to_string(Q(numerator, denominator))
                for numerator, denominator in active_overlap
            ],
            "overlap_excluded_before_scoring": True,
        },
        "T_equals_5_calibration_present_in_lattice_population": calibration_present,
        "T_equals_5_excluded_unconditionally": True,
        "genuinely_outside_active_box": len(outside),
        "outside_population_sha256": integer_digest(outside),
        "minimum_outside_height": min(max(item) for item in outside),
    }


def trace_score(trace: int, point_count: int, prime: int) -> float:
    return (2 - trace) / point_count * log(prime)


def build_good_trace_tables() -> dict[int, tuple[dict[str, int] | None, ...]]:
    tables: dict[int, tuple[dict[str, int] | None, ...]] = {}
    for prime in (*HELD_TRACE_BANDS[0], *HELD_TRACE_BANDS[1]):
        symbols: list[dict[str, int] | None] = []
        for residue in range(prime):
            local = reduction_at_residue(residue, prime)
            symbols.append(local if local["kind"] == "good" else None)
        coefficient_a = A_COEFFICIENTS[-1] % prime
        coefficient_b = B_COEFFICIENTS[-1] % prime
        discriminant_core = (4 * coefficient_a**3 + 27 * coefficient_b**2) % prime
        if discriminant_core:
            trace = character_trace(coefficient_a, coefficient_b, prime)
            symbols.append(
                {"kind": "good", "trace": trace, "point_count": prime + 1 - trace}
            )
        else:
            symbols.append(None)
        tables[prime] = tuple(symbols)
    return tables


def held_trace_data(
    numerator: int,
    denominator: int,
    tables: dict[int, tuple[dict[str, int] | None, ...]],
) -> dict[str, Any]:
    records = []
    band_values = []
    for band in HELD_TRACE_BANDS:
        value = 0.0
        used = 0
        skipped_bad = 0
        skipped_infinity_bad = 0
        for prime in band:
            symbol_index = prime
            if denominator % prime:
                symbol_index = numerator * pow(denominator, -1, prime) % prime
            local = tables[prime][symbol_index]
            if local is None:
                skipped_bad += 1
                if symbol_index == prime:
                    skipped_infinity_bad += 1
                continue
            value += trace_score(local["trace"], local["point_count"], prime)
            used += 1
            records.append({"prime": prime, "trace": local["trace"]})
        band_values.append(
            {
                "primes": list(band),
                "score": f"{value:.12f}",
                "good_primes_used": used,
                "bad_primes_skipped": skipped_bad,
                "bad_infinity_symbols_skipped": skipped_infinity_bad,
            }
        )
    total = sum(float(item["score"]) for item in band_values)
    return {
        "bands": band_values,
        "total_score": f"{total:.12f}",
        "traces": records,
    }


def score_population(
    population: Sequence[tuple[int, int]],
    core: Sequence[int],
    groups: dict[int, tuple[dict[str, Any], ...]],
) -> list[dict[str, Any]]:
    trace_tables = build_good_trace_tables()
    small_primes = tuple(primes_up_to(LOCAL_PRIME_BOUND))
    records: list[dict[str, Any]] = []
    for numerator, denominator in population:
        core_value = homogenized_value(core, numerator, denominator)
        if core_value == 0:
            # These are exactly the singular rational collision parameters.
            continue
        valuations = []
        radical_proxy = log(abs(core_value))
        for prime in small_primes:
            valuation = integer_valuation(core_value, prime)
            if valuation:
                valuations.append({"prime": prime, "valuation": valuation})
                radical_proxy -= (valuation - 1) * log(prime)
        forced_choices = []
        for prime in SELECTED_POWER_PRIMES:
            choice = matching_choice(numerator, denominator, groups[prime])
            if choice is None:
                raise AssertionError("a scored candidate lost a selected local union")
            actual = integer_valuation(core_value, prime)
            if actual < choice["forced_core_valuation"]:
                raise AssertionError("a selected p-adic valuation guarantee failed")
            forced_choices.append(
                {
                    "prime": prime,
                    "residue": choice["residue"],
                    "modulus": choice["modulus"],
                    "forced_core_valuation": choice["forced_core_valuation"],
                    "actual_core_valuation": actual,
                    "reduction": choice["reduction"]["kind"],
                }
            )
        held = held_trace_data(numerator, denominator, trace_tables)
        held_value = float(held["total_score"])
        records.append(
            {
                "t": rational_to_string(Q(numerator, denominator)),
                "numerator": numerator,
                "denominator": denominator,
                "height": max(numerator, denominator),
                "log_abs_homogeneous_core": f"{log(abs(core_value)):.12f}",
                "bounded_radical_log_proxy": f"{radical_proxy:.12f}",
                "small_prime_valuations_through_199": valuations,
                "selected_power_constraints": forced_choices,
                "held_good_prime_trace": held,
                "combined_proxy": f"{-radical_proxy + 6.0 * held_value:.12f}",
            }
        )
    return records


def select_conductor_population(
    records: Sequence[dict[str, Any]], keep: int
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if keep < 1:
        raise ValueError("the conductor keep count must be positive")
    by_radical = sorted(
        records,
        key=lambda item: (
            float(item["bounded_radical_log_proxy"]),
            item["height"],
            item["numerator"],
            item["denominator"],
        ),
    )
    by_combined = sorted(
        records,
        key=lambda item: (
            -float(item["combined_proxy"]),
            item["height"],
            item["numerator"],
            item["denominator"],
        ),
    )
    # Pure trace maximization over the entire lattice box selects enormous
    # heights and is not conductor-relevant.  Freeze a conductor-blind radical
    # frontier first, then use the held band only within that frontier.
    all_by_held = sorted(
        records,
        key=lambda item: (
            -float(item["held_good_prime_trace"]["total_score"]),
            item["height"],
            item["numerator"],
            item["denominator"],
        ),
    )
    held_frontier = by_radical[:512]
    by_held = sorted(
        held_frontier,
        key=lambda item: (
            -float(item["held_good_prime_trace"]["total_score"]),
            item["height"],
            item["numerator"],
            item["denominator"],
        ),
    )
    radical_rank = {item["t"]: index for index, item in enumerate(by_radical, 1)}
    combined_rank = {item["t"]: index for index, item in enumerate(by_combined, 1)}
    held_rank = {item["t"]: index for index, item in enumerate(all_by_held, 1)}
    selected: list[dict[str, Any]] = []
    seen: set[str] = set()
    lanes = (
        ("bounded-radical", by_radical[:10]),
        ("combined-radical-trace", by_combined[:8]),
        ("held-trace", by_held[:4]),
    )
    for lane, items in lanes:
        for item in items:
            if item["t"] in seen:
                continue
            copy = dict(item)
            copy["selection"] = {
                "first_lane": lane,
                "bounded_radical_rank": radical_rank[item["t"]],
                "combined_proxy_rank": combined_rank[item["t"]],
                "held_trace_rank": held_rank[item["t"]],
            }
            selected.append(copy)
            seen.add(item["t"])
            if len(selected) == keep:
                break
        if len(selected) == keep:
            break
    if len(selected) < keep:
        for item in by_combined:
            if item["t"] in seen:
                continue
            copy = dict(item)
            copy["selection"] = {
                "first_lane": "combined-fill",
                "bounded_radical_rank": radical_rank[item["t"]],
                "combined_proxy_rank": combined_rank[item["t"]],
                "held_trace_rank": held_rank[item["t"]],
            }
            selected.append(copy)
            seen.add(item["t"])
            if len(selected) == keep:
                break
    return selected, {
        "population_closed_before_conductor_calls": True,
        "lanes": [
            {"name": "bounded-radical", "predeclared_take": 10},
            {"name": "combined-radical-trace", "predeclared_take": 8},
            {
                "name": "held-trace",
                "predeclared_take": 4,
                "predeclared_bounded_radical_frontier": 512,
            },
        ],
        "deduplicated_keep": keep,
        "selected_parameters": [item["t"] for item in selected],
        "selection_sha256": hashlib.sha256(
            "\n".join(item["t"] for item in selected).encode()
        ).hexdigest(),
    }


def conductor_phase(
    construction: SixRootMestreConstruction,
    selected: list[dict[str, Any]],
    *,
    timeout: float,
    stack_bytes: int,
) -> None:
    for index, record in enumerate(selected, 1):
        print(f"conductor {index}/{len(selected)} T={record['t']}", flush=True)
        parameter = Q(record["numerator"], record["denominator"])
        try:
            conductor = capped_minimal_curve_data(
                construction.primitive_jacobian_coefficients(parameter),
                timeout=timeout,
                stack_bytes=stack_bytes,
            )
            exact_log = Decimal(conductor["conductor"]).ln()
            record["conductor_phase"] = {
                "status": "completed exact PARI minimal model and conductor",
                **conductor,
                "decimal_log_recomputed_from_exact_conductor": str(exact_log),
                "below_strict_log_conductor_target": exact_log < TARGET_LOG_CONDUCTOR,
            }
        except CappedProcessTimeout:
            record["conductor_phase"] = {
                "status": "timeout",
                "timeout_seconds": timeout,
                "retry_count": 0,
            }
        except Exception as error:
            record["conductor_phase"] = {
                "status": "error",
                "error": str(error)[:1000],
                "retry_count": 0,
            }


def stable_height_triage(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[dict[str, Any], ...], tuple[tuple[Fraction, Fraction], ...]]:
    height = height_matrix_replay(
        coefficients,
        points,
        precisions=(72, 120),
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    return height, numerical_subset(points, height)


def triage_subtarget(
    construction: SixRootMestreConstruction,
    record: dict[str, Any],
    *,
    point_timeout: float,
    height_timeout: float,
    stack_bytes: int,
    certificate_prime_bound: int,
) -> None:
    parameter = Q(record["numerator"], record["denominator"])
    coefficients = construction.primitive_jacobian_coefficients(parameter)
    quartic_coefficients = construction.primitive_quartic_coefficients(parameter)
    visible_quartic = primitive_visible_points(construction, parameter)
    visible_jacobian = tuple(
        quartic_point_to_jacobian(construction, parameter, point)
        for point in visible_quartic
    )
    triage: dict[str, Any] = {
        "status": "started exact point and stable numerical-height triage",
        "visible_quartic_points": len(visible_quartic),
        "visible_quartic_point_sha256": point_digest(visible_quartic),
        "visible_jacobian_points": len(visible_jacobian),
        "visible_jacobian_point_sha256": point_digest(visible_jacobian),
        "numerical_rank_is_not_a_rank_certificate": True,
        "searches": [],
    }
    try:
        visible_height, visible_subset = stable_height_triage(
            coefficients,
            visible_jacobian,
            timeout=height_timeout,
            stack_bytes=stack_bytes,
        )
        triage["visible_height_runs"] = list(visible_height)
        triage["visible_stable_numerical_rank"] = visible_height[-1]["numerical_rank"]
        accumulated_quartic = canonical_signless_points(visible_quartic)
        maximum_rank = visible_height[-1]["numerical_rank"]
        certificate = None
        if maximum_rank >= 21:
            certificate = finite_reduction_attempt(
                coefficients, visible_subset, prime_bound=certificate_prime_bound
            )

        for height_bound in (50_000, 250_000, 1_000_000):
            search_record: dict[str, Any] = {
                "height_bound": height_bound,
                "timeout_seconds": point_timeout,
            }
            try:
                found = bounded_quartic_points(
                    quartic_coefficients,
                    height_bound=height_bound,
                    timeout=point_timeout,
                    stack_bytes=stack_bytes,
                )
                canonical_found = canonical_signless_points(found)
                accumulated_quartic = canonical_signless_points(
                    (*accumulated_quartic, *canonical_found)
                )
                jacobian_points = tuple(
                    quartic_point_to_jacobian(construction, parameter, point)
                    for point in accumulated_quartic
                )
                height, subset = stable_height_triage(
                    coefficients,
                    jacobian_points,
                    timeout=height_timeout,
                    stack_bytes=stack_bytes,
                )
                stable_rank = height[-1]["numerical_rank"]
                maximum_rank = max(maximum_rank, stable_rank)
                search_record.update(
                    {
                        "status": "completed",
                        "raw_affine_points_parsed": len(found),
                        "canonical_signless_points_found": len(canonical_found),
                        "accumulated_distinct_abscissas": len(accumulated_quartic),
                        "accumulated_quartic_point_sha256": point_digest(accumulated_quartic),
                        "height_runs": list(height),
                        "stable_numerical_rank": stable_rank,
                    }
                )
                if stable_rank >= 21 and certificate is None:
                    certificate = finite_reduction_attempt(
                        coefficients, subset, prime_bound=certificate_prime_bound
                    )
            except CappedProcessTimeout:
                search_record.update(
                    {"status": "timeout", "retry_count": 0}
                )
            except Exception as error:
                search_record.update(
                    {"status": "error", "error": str(error)[:1000], "retry_count": 0}
                )
            triage["searches"].append(search_record)
        triage["maximum_stable_numerical_rank"] = maximum_rank
        triage["finite_reduction_certificate_attempt"] = certificate
        triage["status"] = "completed declared bounded point searches"
    except CappedProcessTimeout:
        triage.update(
            {
                "status": "visible-height-timeout",
                "timeout_seconds": height_timeout,
                "retry_count": 0,
            }
        )
    except Exception as error:
        triage.update(
            {"status": "error", "error": str(error)[:1000], "retry_count": 0}
        )
    record["point_rank_triage"] = triage


def calibration_record(
    construction: SixRootMestreConstruction,
    *,
    conductor_timeout: float,
    height_timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    parameter = Q(5)
    coefficients = construction.primitive_jacobian_coefficients(parameter)
    conductor = capped_minimal_curve_data(
        coefficients, timeout=conductor_timeout, stack_bytes=stack_bytes
    )
    quartic = primitive_visible_points(construction, parameter)
    jacobian = tuple(
        quartic_point_to_jacobian(construction, parameter, point) for point in quartic
    )
    height, _ = stable_height_triage(
        coefficients,
        jacobian,
        timeout=height_timeout,
        stack_bytes=stack_bytes,
    )
    return {
        "parameter": "5",
        "status": "pinned calibration and unconditional search exclusion",
        "excluded_from_lattice_and_conductor_candidate_population": True,
        "conductor": conductor,
        "visible_stable_numerical_rank": height[-1]["numerical_rank"],
        "visible_height_runs": list(height),
        "expected_prior_log_conductor_prefix": "79.729318123910",
        "expected_prior_visible_rank": 10,
        "prior_values_recovered": (
            conductor["log_conductor"].startswith("79.729318123910")
            and height[-1]["numerical_rank"] == 10
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coefficient-radius", type=int, default=DEFAULT_COEFFICIENT_RADIUS)
    parser.add_argument("--conductor-keep", type=int, default=DEFAULT_CONDUCTOR_KEEP)
    parser.add_argument("--subtarget-triage-cap", type=int, default=DEFAULT_SUBTARGET_TRIAGE_CAP)
    parser.add_argument("--conductor-timeout", type=float, default=20.0)
    parser.add_argument("--point-timeout", type=float, default=60.0)
    parser.add_argument("--height-timeout", type=float, default=40.0)
    parser.add_argument("--stack-bytes", type=int, default=256_000_000)
    parser.add_argument("--certificate-prime-bound", type=int, default=500)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not (1 <= args.coefficient_radius <= 40):
        raise SystemExit("--coefficient-radius must lie in [1,40]")
    if not (1 <= args.conductor_keep <= 24):
        raise SystemExit("--conductor-keep must lie in [1,24]")
    if not (0 <= args.subtarget_triage_cap <= 6):
        raise SystemExit("--subtarget-triage-cap must lie in [0,6]")
    if not (1 <= args.conductor_timeout <= 30):
        raise SystemExit("--conductor-timeout must lie in [1,30]")
    if not (1 <= args.point_timeout <= 60 and 1 <= args.height_timeout <= 60):
        raise SystemExit("point and height timeouts must lie in [1,60]")
    if not (32_000_000 <= args.stack_bytes <= 512_000_000):
        raise SystemExit("--stack-bytes must lie in [32000000,512000000]")
    if not (50 <= args.certificate_prime_bound <= 1000):
        raise SystemExit("--certificate-prime-bound must lie in [50,1000]")

    construction = SixRootMestreConstruction(tuple(Q(root) for root in ROOTS))
    core, geometry = derive_exact_geometry(construction)
    print("discovering complete p^4 root profiles through p=199", flush=True)
    local_profiles, groups = discover_local_profiles(core)
    classes = combine_selected_balls(groups)
    population, lattice = enumerate_lattice_population(
        classes, groups, coefficient_radius=args.coefficient_radius
    )
    print(
        f"scoring {len(population)} representatives outside the active box",
        flush=True,
    )
    scored = score_population(population, core, groups)
    scored_digest = hashlib.sha256(
        "\n".join(
            f"{item['t']}|{item['bounded_radical_log_proxy']}|"
            f"{item['held_good_prime_trace']['total_score']}"
            for item in sorted(scored, key=lambda candidate: candidate["t"])
        ).encode()
    ).hexdigest()
    selected, selection = select_conductor_population(scored, args.conductor_keep)
    conductor_phase(
        construction,
        selected,
        timeout=args.conductor_timeout,
        stack_bytes=args.stack_bytes,
    )

    subtarget = [
        record
        for record in selected
        if record.get("conductor_phase", {}).get("below_strict_log_conductor_target")
    ]
    triaged = subtarget[: args.subtarget_triage_cap]
    for index, record in enumerate(triaged, 1):
        print(f"point triage {index}/{len(triaged)} T={record['t']}", flush=True)
        triage_subtarget(
            construction,
            record,
            point_timeout=args.point_timeout,
            height_timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
            certificate_prime_bound=args.certificate_prime_bound,
        )

    calibration = calibration_record(
        construction,
        conductor_timeout=args.conductor_timeout,
        height_timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    if not calibration["prior_values_recovered"]:
        raise AssertionError("the T=5 calibration changed")

    target_hits = []
    for record in selected:
        triage = record.get("point_rank_triage", {})
        certificate = triage.get("finite_reduction_certificate_attempt")
        if (
            record.get("conductor_phase", {}).get("below_strict_log_conductor_target")
            and certificate is not None
            and certificate.get("certified_algebraic_rank_lower_bound", 0) >= 21
        ):
            target_hits.append(
                {
                    "t": record["t"],
                    "conductor": record["conductor_phase"]["conductor"],
                    "certified_rank_lower_bound": certificate[
                        "certified_algebraic_rank_lower_bound"
                    ],
                }
            )

    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[2]
    output = args.output if args.output.is_absolute() else repo_root / args.output
    modulus_distribution = Counter(item["modulus"] for item in classes)
    result_payload = {
        "geometry_core_sha256": integer_digest((core,)),
        "local_profile_sha256": hashlib.sha256(
            json.dumps(local_profiles, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest(),
        "crt_class_sha256": integer_digest(
            (item["residue"], item["modulus"]) for item in classes
        ),
        "outside_population_sha256": lattice["outside_population_sha256"],
        "scored_population_sha256": scored_digest,
        "conductor_selection_sha256": selection["selection_sha256"],
        "target_hits": target_hits,
    }
    result_sha256 = hashlib.sha256(
        json.dumps(result_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    artifact = {
        "schema_version": 1,
        "status": (
            "exact local/CRT experiment; rank conclusions require exact certificates"
        ),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_less_than": str(TARGET_LOG_CONDUCTOR),
            "hits": target_hits,
        },
        "parameters": {
            "roots": list(ROOTS),
            "local_prime_bound": LOCAL_PRIME_BOUND,
            "local_target_exponent": LOCAL_TARGET_EXPONENT,
            "local_root_cap": LOCAL_ROOT_CAP,
            "selected_power_primes": list(SELECTED_POWER_PRIMES),
            "held_trace_bands": [list(band) for band in HELD_TRACE_BANDS],
            "coefficient_radius": args.coefficient_radius,
            "conductor_keep": args.conductor_keep,
            "subtarget_triage_cap": args.subtarget_triage_cap,
            "timeouts_seconds": {
                "conductor": args.conductor_timeout,
                "point_search_each": args.point_timeout,
                "height_matrix_each": args.height_timeout,
            },
            "stack_bytes": args.stack_bytes,
            "certificate_prime_bound": args.certificate_prime_bound,
        },
        "geometry": geometry,
        "p_adic_discovery": {
            "status": "complete p^4 digit lifting; no truncated root sets",
            "prime_count": len(local_profiles),
            "profiles": local_profiles,
            "selected_clean_unions": {
                str(prime): list(groups[prime]) for prime in SELECTED_POWER_PRIMES
            },
            "selection_rule": (
                "maximal ball exponent below 4, forced core valuation at least 4, "
                "and c4 a unit in the singular reduction"
            ),
        },
        "crt_lattice": {
            "status": "exact CRT and exact two-dimensional Gauss reduction",
            "crt_classes": len(classes),
            "crt_modulus_distribution": {
                str(modulus): count
                for modulus, count in sorted(modulus_distribution.items())
            },
            **lattice,
        },
        "ranking": {
            "status": "conductor-blind bounded radical and held-trace ranking",
            "bounded_radical_proxy_definition": (
                "log|P_h(a,b)| - sum_{p<=199}(v_p(P_h(a,b))-1)log(p); "
                "unknown large-prime cofactor is treated as squarefree"
            ),
            "held_trace_score_definition": (
                "sum over held good primes of (2-a_p)/(p+1-a_p)*log(p)"
            ),
            "scored_nonsingular_population": len(scored),
            "scored_population_sha256": scored_digest,
            "selection": selection,
        },
        "calibration": calibration,
        "conductor_first_records": selected,
        "outcome": {
            "conductor_completed": sum(
                record["conductor_phase"]["status"].startswith("completed")
                for record in selected
            ),
            "conductor_timeouts": sum(
                record["conductor_phase"]["status"] == "timeout"
                for record in selected
            ),
            "conductor_errors": sum(
                record["conductor_phase"]["status"] == "error"
                for record in selected
            ),
            "subtarget_conductors": len(subtarget),
            "subtarget_triaged": len(triaged),
            "subtarget_untriaged_due_predeclared_cap": max(0, len(subtarget) - len(triaged)),
            "maximum_stable_numerical_rank": max(
                (
                    record.get("point_rank_triage", {}).get(
                        "maximum_stable_numerical_rank", 0
                    )
                    for record in selected
                ),
                default=0,
            ),
            "target_hits": target_hits,
            "no_rank_claim_from_numerical_height_or_bounded_point_search": True,
        },
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version_capped(),
        },
        "provenance": {
            "script": str(script_path.relative_to(repo_root)),
            "script_sha256": sha256_file(script_path),
            "reproducing_command": REPRODUCING_COMMAND,
            "active_cpp_scanner_read_or_invoked": False,
            "global_root_tuple_driver_read_or_invoked": False,
            "subprocess_policy": (
                "one foreground process group per declared call; hard kill on timeout; "
                "no retries and no background/orphan processes"
            ),
        },
        "result_sha256": result_sha256,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {output}", flush=True)
    print(f"result_sha256={result_sha256}", flush=True)


if __name__ == "__main__":
    main()
