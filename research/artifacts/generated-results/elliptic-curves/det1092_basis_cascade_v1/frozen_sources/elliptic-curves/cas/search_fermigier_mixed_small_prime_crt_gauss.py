#!/usr/bin/env python3
"""Bounded mixed-small-prime CRT--Gauss search in the Fermigier family.

The complete discovery lane exhausts every four-prime subset of the declared
local targets in a radius-4 reduced-basis box.  The ICARM-282 profile receives
a separate radius-64 calibration lane.  Both affine and infinity charts of
P^1(Q_p) are included.  Exact radical proxies select a compact point-search
prefix; the point search itself is only the declared integer-offset slice.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction as Q
import hashlib
import heapq
from itertools import combinations, product
import json
from math import gcd, isqrt
from pathlib import Path
import platform
import shutil
import sys
from typing import Any, Iterable, Sequence

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "elliptic-curves"), str(Path(__file__).resolve().parent)]

from ecsearch.fermigier import (  # noqa: E402
    FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS as H_COEFFICIENTS,
    evaluate_polynomial,
)
from ecsearch.fermigier_rank import (  # noqa: E402
    section_and_point_cloud_differences,
    specialize_fermigier_rank_sections,
)
from ecsearch.rank_certification import select_independent_subset  # noqa: E402
from multiple_root_lifting import (  # noqa: E402
    PrimePowerRootResult,
    all_roots_mod_prime_power,
)
from search_fermigier_high_power_crt_gauss import (  # noqa: E402
    load_prior,
    primes_up_to,
)
from search_mestre_02557104116148_power_root_crt import (  # noqa: E402
    ProjectiveBall,
    combine_rows,
    kernel_basis,
)


LOCAL_TARGETS = (
    (5, 3),
    (7, 5),
    (11, 4),
    (13, 3),
    (17, 4),
    (19, 4),
    (23, 2),
    (29, 3),
    (31, 2),
)
PROFILE_SIZE = 4
DISCOVERY_RADIUS = 4
CALIBRATION_PROFILE = (11, 13, 23, 31)
CALIBRATION_RADIUS = 64
HEIGHT_CAP = 100_000
FEATURE_KEEP = 256
POINT_KEEP = 32
POINT_OFFSET_BOUND = 5_000
TRIAL_PRIME_BOUND = 2_000
ICARM_282_U = Q(11671, 42)
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/fermigier_mixed_small_prime_crt_gauss_v1.json"
EXPECTED_DEFAULT_POPULATION = {
    "projective_classes": 36_352,
    "projective_class_sha256": "b613b62a4f1ad1807a408f0cb61cf03faf7f69208442f896412ea97aefa9d665",
    "bounded_vectors": 10_062_080,
    "invalid_vectors": 2_233_984,
    "above_height_cap": 649_972,
    "unique_parameters": 963_684,
    "parameter_sha256": "b96a68752dfbd7cdaf4030fc4c4656adedb9300a3a6772a0baa18eb514da1e75",
    "icarm_282_recovered": True,
}
EXPECTED_DEFAULT_RETENTION_COUNTS = {9: 1, 10: 6, 11: 37, 12: 337}


@dataclass(frozen=True)
class CandidateFeature:
    numerator: int
    denominator: int
    homogeneous_h: int
    target_valuations: tuple[tuple[int, int], ...]
    known_radical_upper: int
    known_repeated_divisor: int

    @property
    def parameter(self) -> Q:
        return Q(self.numerator, self.denominator)

    @property
    def height(self) -> int:
        return max(self.numerator, self.denominator)


def rational_text(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def sequence_digest(lines: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for line in lines:
        digest.update((line + "\n").encode())
    return digest.hexdigest()


def valuation(value: int, prime: int) -> int:
    value = abs(value)
    exponent = 0
    while value and value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def homogeneous_h(numerator: int, denominator: int) -> int:
    degree = len(H_COEFFICIENTS) - 1
    value = H_COEFFICIENTS[degree]
    denominator_power = denominator
    for index in range(degree - 1, -1, -1):
        value = value * numerator + H_COEFFICIENTS[index] * denominator_power
        denominator_power *= denominator
    return abs(value)


def _filtered_infinity_result(prime: int, exponent: int) -> PrimePowerRootResult:
    result = all_roots_mod_prime_power(
        tuple(reversed(H_COEFFICIENTS)), prime, exponent, max_roots=100_000
    )
    roots = tuple(root for root in result.roots if root % prime == 0)
    return PrimePowerRootResult(
        prime,
        exponent,
        prime**exponent,
        roots,
        result.level_counts,
        result.candidate_digits_checked,
    )


def _expand_balls(balls: Sequence[ProjectiveBall], exponent: int) -> tuple[int, ...]:
    covered: set[int] = set()
    for ball in balls:
        covered.update(
            ball.residue + digit * ball.modulus
            for digit in range(ball.prime ** (exponent - ball.exponent))
        )
    return tuple(sorted(covered))


def discover_groups() -> tuple[dict[int, tuple[ProjectiveBall, ...]], list[dict[str, Any]]]:
    groups: dict[int, tuple[ProjectiveBall, ...]] = {}
    records = []
    for prime, exponent in LOCAL_TARGETS:
        affine = all_roots_mod_prime_power(
            H_COEFFICIENTS, prime, exponent, max_roots=100_000
        )
        infinity = _filtered_infinity_result(prime, exponent)
        affine_balls = tuple(
            ProjectiveBall(prime, exponent, "affine", ball.exponent, ball.residue)
            for ball in affine.maximal_balls()
        )
        infinity_balls = tuple(
            ProjectiveBall(prime, exponent, "infinity", ball.exponent, ball.residue)
            for ball in infinity.maximal_balls()
        )
        if _expand_balls(affine_balls, exponent) != affine.roots:
            raise AssertionError("affine root-ball compression changed the root set")
        if _expand_balls(infinity_balls, exponent) != infinity.roots:
            raise AssertionError("infinity root-ball compression changed the root set")
        combined = tuple(sorted((*affine_balls, *infinity_balls)))
        if not combined:
            raise AssertionError(f"p={prime} has no projective target branch")
        groups[prime] = combined
        records.append(
            {
                "prime": prime,
                "target_h_valuation": exponent,
                "affine_level_root_counts": list(affine.level_counts),
                "affine_target_root_count": len(affine.roots),
                "infinity_target_root_count": len(infinity.roots),
                "maximal_balls": [
                    {
                        "chart": ball.chart,
                        "residue": ball.residue,
                        "exponent": ball.exponent,
                        "modulus": ball.modulus,
                    }
                    for ball in combined
                ],
            }
        )
    return groups, records


def _normalized_vector(
    first: tuple[int, int], second: tuple[int, int], left: int, right: int, modulus: int
) -> tuple[int, int] | None:
    numerator = left * first[0] + right * second[0]
    denominator = left * first[1] + right * second[1]
    if denominator == 0:
        return None
    common = gcd(abs(numerator), abs(denominator))
    if gcd(common, modulus) != 1:
        return None
    numerator //= common
    denominator //= common
    if denominator < 0:
        numerator = -numerator
        denominator = -denominator
    numerator = abs(numerator)
    if numerator == 0 or gcd(numerator, denominator) != 1:
        return None
    return numerator, denominator


def enumerate_population(
    groups: dict[int, tuple[ProjectiveBall, ...]], *, height_cap: int
) -> tuple[dict[tuple[int, int], int], list[tuple[int, ...]], dict[str, Any]]:
    primes = tuple(prime for prime, _ in LOCAL_TARGETS)
    profiles = list(combinations(primes, PROFILE_SIZE))
    population: dict[tuple[int, int], int] = {}
    counters = Counter()
    class_digest = hashlib.sha256()
    for profile_index, profile in enumerate(profiles):
        radius = CALIBRATION_RADIUS if profile == CALIBRATION_PROFILE else DISCOVERY_RADIUS
        for choices in product(*(groups[prime] for prime in profile)):
            counters["projective_classes"] += 1
            coefficient_a, coefficient_b, modulus = combine_rows(choices)
            class_digest.update(
                (
                    f"{profile}|{coefficient_a}|{coefficient_b}|{modulus}|"
                    + "/".join(choice.label for choice in choices)
                    + "\n"
                ).encode()
            )
            first, second = kernel_basis(coefficient_a, coefficient_b, modulus)
            for left in range(-radius, radius + 1):
                for right in range(-radius, radius + 1):
                    if left == 0 and right == 0:
                        continue
                    counters["bounded_vectors"] += 1
                    pair = _normalized_vector(first, second, left, right, modulus)
                    if pair is None:
                        counters["invalid_vectors"] += 1
                        continue
                    if max(pair) > height_cap:
                        counters["above_height_cap"] += 1
                        continue
                    population[pair] = population.get(pair, 0) | (1 << profile_index)
    ordered = sorted(population, key=lambda pair: (max(pair), pair[0], pair[1]))
    return population, profiles, {
        "profile_count": len(profiles),
        "projective_classes": counters["projective_classes"],
        "projective_class_sha256": class_digest.hexdigest(),
        "bounded_vectors": counters["bounded_vectors"],
        "invalid_vectors": counters["invalid_vectors"],
        "above_height_cap": counters["above_height_cap"],
        "unique_parameters": len(population),
        "parameter_sha256": sequence_digest(f"{a}/{b}" for a, b in ordered),
        "icarm_282_recovered": (ICARM_282_U.numerator, ICARM_282_U.denominator) in population,
    }


def feature(pair: tuple[int, int]) -> CandidateFeature | None:
    numerator, denominator = pair
    value = homogeneous_h(numerator, denominator)
    if value == 0:
        return None
    local = tuple((prime, valuation(value, prime)) for prime, _ in LOCAL_TARGETS)
    repeated = 1
    for prime, exponent in local:
        if exponent > 1:
            repeated *= prime ** (exponent - 1)
    return CandidateFeature(numerator, denominator, value, local, value // repeated, repeated)


def _bounded_smallest(heap: list[tuple[int, int, int, tuple[int, int]]], key: tuple[int, int, int], pair: tuple[int, int], limit: int) -> None:
    item = (-key[0], -key[1], -key[2], pair)
    if len(heap) < limit:
        heapq.heappush(heap, item)
    elif item > heap[0]:
        heapq.heapreplace(heap, item)


def select_feature_union(
    population: dict[tuple[int, int], int], excluded: set[Q]
) -> tuple[dict[tuple[int, int], CandidateFeature], dict[str, Any]]:
    radical_heap: list[tuple[int, int, int, tuple[int, int]]] = []
    height_heap: list[tuple[int, int, int, tuple[int, int]]] = []
    powerful_heap: list[tuple[int, int, int, tuple[int, int]]] = []
    singular = 0
    fresh = 0
    for pair in population:
        parameter = Q(*pair)
        if parameter in excluded or parameter == ICARM_282_U:
            continue
        fresh += 1
        row = feature(pair)
        if row is None:
            singular += 1
            continue
        _bounded_smallest(
            radical_heap,
            (row.known_radical_upper, row.height, row.numerator),
            pair,
            FEATURE_KEEP,
        )
        _bounded_smallest(
            height_heap,
            (row.height, row.known_radical_upper, row.numerator),
            pair,
            FEATURE_KEEP // 4,
        )
        # Negating the powerful part turns largest-powerful into smallest-key.
        _bounded_smallest(
            powerful_heap,
            (-row.known_repeated_divisor, row.known_radical_upper, row.height),
            pair,
            FEATURE_KEEP // 4,
        )
    selected_pairs = {item[3] for heap in (radical_heap, height_heap, powerful_heap) for item in heap}
    rows = {pair: feature(pair) for pair in selected_pairs}
    if any(row is None for row in rows.values()):
        raise AssertionError("a singular candidate survived feature selection")
    return rows, {
        "fresh_parameters_before_singularity": fresh,
        "singular_parameters": singular,
        "feature_union_count": len(rows),
        "lanes": {
            "smallest_known_radical_upper": FEATURE_KEEP,
            "lowest_height": FEATURE_KEEP // 4,
            "largest_known_powerful_part": FEATURE_KEEP // 4,
        },
    }


def trial_factor_record(row: CandidateFeature, trial_primes: Sequence[int]) -> dict[str, Any]:
    remaining = row.homogeneous_h
    factors: dict[str, int] = {}
    repeated = 1
    for prime in trial_primes:
        if remaining % prime:
            continue
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        factors[str(prime)] = exponent
        repeated *= prime ** (exponent - 1)
    return {
        "exact_trial_factorization": factors,
        "exact_repeated_prime_divisor": str(repeated),
        "exact_radical_upper_integer": str(row.homogeneous_h // repeated),
        "unfactored_cofactor_bits": remaining.bit_length(),
    }


def baseline_retention_rank(parameter: Q) -> int:
    """Return the exact mod-5 rank seen among the twelve generic directions."""

    specialization = specialize_fermigier_rank_sections(parameter)
    selected, _ = select_independent_subset(
        specialization.canonical_model,
        specialization.section_differences,
        relation_prime=5,
        maximum_reduction_prime=1_000,
    )
    return len(selected)


def _square_root(value: Q) -> Q | None:
    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator * numerator != value.numerator or denominator * denominator != value.denominator:
        return None
    return Q(numerator, denominator)


def point_screen(parameter: Q, offset_bound: int) -> dict[str, Any]:
    specialization = specialize_fermigier_rank_sections(parameter)
    shift = 2 * parameter
    known_x = {point[0] for point in specialization.quartic_points}
    found: dict[Q, Q] = {}
    for sign in (-1, 1):
        for offset in range(-offset_bound, offset_bound + 1):
            x_value = sign * shift + offset
            if x_value in known_x or x_value in found:
                continue
            root = _square_root(evaluate_polynomial(specialization.quartic_model.quartic, x_value))
            if root is not None:
                found[x_value] = root
    searched_points = tuple(found.items())
    cloud = section_and_point_cloud_differences(specialization, searched_points)
    attempts = []
    for relation_prime in (5, 3, 7):
        selected, certificate = select_independent_subset(
            specialization.canonical_model,
            cloud,
            relation_prime=relation_prime,
            maximum_reduction_prime=2_000,
        )
        attempts.append((len(selected), selected, certificate))
        if len(selected) == len(cloud):
            break
    certified_lower, selected, certificate = max(
        attempts, key=lambda item: (item[0], -item[2].relation_prime)
    )
    certificate_json = certificate.to_json_object()
    return {
        "integer_offset_slice": f"x=+/-2u+n, |n|<={offset_bound}",
        "extra_quartic_points": [
            [rational_text(x_value), rational_text(y_value)]
            for x_value, y_value in searched_points
        ],
        "extra_point_count": len(searched_points),
        "deduplicated_jacobian_difference_count": len(cloud),
        "selected_subset_contains_all_twelve_baseline_indices": set(range(12)).issubset(selected),
        "selected_indices": list(selected),
        "certified_rank_lower_bound": certified_lower,
        "finite_reduction_certificate": certificate_json,
        "boundary": "The integer-offset slice is bounded and not a complete quartic point search.",
    }


def profile_support(mask: int, profiles: Sequence[tuple[int, ...]]) -> list[list[int]]:
    return [list(profile) for index, profile in enumerate(profiles) if mask & (1 << index)]


def stable_result_digest(artifact: dict[str, Any]) -> str:
    stable = {key: value for key, value in artifact.items() if key not in {"generated_at_utc", "result_sha256"}}
    return hashlib.sha256(json.dumps(stable, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def run(*, height_cap: int, point_keep: int, offset_bound: int) -> dict[str, Any]:
    groups, local_records = discover_groups()
    population, profiles, population_record = enumerate_population(groups, height_cap=height_cap)
    if height_cap == HEIGHT_CAP:
        observed = {key: population_record[key] for key in EXPECTED_DEFAULT_POPULATION}
        if observed != EXPECTED_DEFAULT_POPULATION:
            raise AssertionError(f"the default closed population changed: {observed}")
    prior, prior_record = load_prior(ROOT)
    selected, selection_record = select_feature_union(population, prior)
    trial_primes = primes_up_to(TRIAL_PRIME_BOUND)
    trial_rows = []
    retention_counts = Counter()
    for pair, row in selected.items():
        assert row is not None
        retention_rank = baseline_retention_rank(row.parameter)
        retention_counts[retention_rank] += 1
        if retention_rank < 12:
            continue
        trial = trial_factor_record(row, trial_primes)
        trial_rows.append((int(trial["exact_radical_upper_integer"]), row.height, pair, row, trial))
    trial_rows.sort(key=lambda item: (item[0], item[1], item[2]))
    if height_cap == HEIGHT_CAP and dict(retention_counts) != EXPECTED_DEFAULT_RETENTION_COUNTS:
        raise AssertionError("the exact baseline-retention census changed")
    promoted = trial_rows[:point_keep]
    selected_records = []
    for _, _, pair, row, trial in promoted:
        parameter = row.parameter
        selected_records.append(
            {
                "adapter_u": rational_text(parameter),
                "literal_shift_s": rational_text(2 * parameter),
                "height": row.height,
                "profile_support": profile_support(population[pair], profiles),
                "target_prime_valuations": {str(p): e for p, e in row.target_valuations if e},
                "radical": trial,
                "baseline_mod5_rank_12_gate": True,
                "point_screen": point_screen(parameter, offset_bound),
                "conductor": (
                    {"status": "not_run", "reason": "PARI/GP is not installed"}
                    if shutil.which("gp") is None
                    else {"status": "not_run", "reason": "this bounded lane is point-first; run exact conductor promotion separately"}
                ),
            }
        )
    calibration = feature((ICARM_282_U.numerator, ICARM_282_U.denominator))
    if calibration is None:
        raise AssertionError("ICARM 282 became singular")
    artifact: dict[str, Any] = {
        "schema": "elliptic-curves.fermigier-mixed-small-prime-crt-gauss.v1",
        "status": "complete bounded projective CRT--Gauss population and integer-offset point screen",
        "claim_level": "bounded computation; only finite-reduction rank lower bounds are claims",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "coordinate_normalization": {
            "search_coordinate": "canonical adapter u=s/2",
            "literal_quartic_shift": "s=2u",
            "sign_quotient": "u and -u give the same canonical coefficients",
        },
        "local_targets": local_records,
        "population": {
            **population_record,
            "profile_size": PROFILE_SIZE,
            "profiles": [list(profile) for profile in profiles],
            "discovery_radius": DISCOVERY_RADIUS,
            "calibration_profile": list(CALIBRATION_PROFILE),
            "calibration_radius": CALIBRATION_RADIUS,
            "height_cap": height_cap,
        },
        "prior_exclusion": prior_record,
        "selection": {
            **selection_record,
            "baseline_retention_gate": {
                "method": "exact finite-reduction subset rank at relation prime 5 through good prime 1000",
                "observed_rank_counts": {
                    str(rank): count for rank, count in sorted(retention_counts.items())
                },
                "rank_12_survivors": retention_counts[12],
            },
            "trial_prime_bound": TRIAL_PRIME_BOUND,
            "point_screen_keep": point_keep,
            "point_offset_bound": offset_bound,
        },
        "icarm_282_calibration": {
            "adapter_u": rational_text(ICARM_282_U),
            "recovered_by_profile": list(CALIBRATION_PROFILE),
            "inside_height_cap": (ICARM_282_U.numerator, ICARM_282_U.denominator) in population,
            "profile_support": profile_support(
                population.get((ICARM_282_U.numerator, ICARM_282_U.denominator), 0), profiles
            ),
            "target_prime_valuations": {str(p): e for p, e in calibration.target_valuations if e},
            "excluded_from_novel_selection": True,
        },
        "selected": selected_records,
        "outcome": {
            "point_screened_candidates": len(selected_records),
            "candidates_with_extra_integer_offset_points": sum(
                bool(row["point_screen"]["extra_point_count"]) for row in selected_records
            ),
            "maximum_certified_rank_lower_bound": max(
                (row["point_screen"]["certified_rank_lower_bound"] for row in selected_records),
                default=0,
            ),
            "exact_conductors_completed": 0,
            "target_met": False,
            "interpretation": "No candidate is promoted without a complete quartic search and exact conductor replay.",
        },
        "software": {"python": platform.python_version(), "pari_gp_available": shutil.which("gp") is not None},
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 "
            "elliptic-curves/cas/search_fermigier_mixed_small_prime_crt_gauss.py"
        ),
    }
    artifact["result_sha256"] = stable_result_digest(artifact)
    return artifact


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height-cap", type=int, default=HEIGHT_CAP)
    parser.add_argument("--point-keep", type=int, default=POINT_KEEP)
    parser.add_argument("--offset-bound", type=int, default=POINT_OFFSET_BOUND)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.height_cap <= HEIGHT_CAP:
        raise SystemExit(f"height cap must lie in [1,{HEIGHT_CAP}]")
    if not 1 <= args.point_keep <= POINT_KEEP or not 1 <= args.offset_bound <= POINT_OFFSET_BOUND:
        raise SystemExit("point keep or offset bound exceeds the pinned maximum")
    artifact = run(
        height_cap=args.height_cap,
        point_keep=args.point_keep,
        offset_bound=args.offset_bound,
    )
    if args.check:
        expected = json.loads(args.output.read_text())
        if artifact["result_sha256"] != expected["result_sha256"]:
            raise AssertionError(
                f"result digest changed: {artifact['result_sha256']} != {expected['result_sha256']}"
            )
        print(f"PASS {artifact['result_sha256']}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps(artifact["outcome"], sort_keys=True))
    print(f"result_sha256={artifact['result_sha256']}")


if __name__ == "__main__":
    main()
