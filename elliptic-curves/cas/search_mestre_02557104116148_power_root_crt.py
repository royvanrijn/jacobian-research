#!/usr/bin/env python3
"""Projective p-adic power-root CRT search around the rank-16 T=62/35 fiber.

This lane is orthogonal to both completed direct searches for the Mestre root
tuple ``(0,25,57,104,116,148)``.  It starts from the exact content-free
degree-20 homogeneous discriminant and the factorization of its value at
``T=62/35``.  For every prime power in that factorization it exhausts *both*
projective charts, lifts every root digit by digit, and compresses the final
root set into disjoint maximal p-adic balls.  In particular, denominator-
divisible branches at 2, 3, and 5 are audited rather than silently discarded.

One ball per prime is combined by exact CRT as a homogeneous linear
congruence in ``(a,b)``.  Exact Gauss reduction then exhausts a pinned radius-4
coefficient box in every resulting lattice.  The complete earlier direct box,
the completed Farey annulus around 62/35, the frozen prior parameter census,
and a height cap are removed before any scoring.  Every retained candidate is
replayed against the projective root balls and its forced valuations are
checked on the exact homogeneous discriminant.

The closed candidate population receives exact radical features.  A fixed
160-element union of radical, powerful-part, height, and sparse-branch lanes
is chosen before conductor calls.  Conductors are completed before staged
quartic searches at H=5k, 50k, 250k, and 1m.  Stable numerical rank at least
15 immediately triggers an exact mod-3 finite-reduction certificate attempt.
Every external call is foreground, process-group capped, and single-attempt.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import gcd
import os
from pathlib import Path
import platform
import shlex
import sys
import time
from typing import Any, Iterable, Sequence

from crt_lattice import crt_pair, gauss_reduce
from mestre_root_tuples import SixRootMestreConstruction
from multiple_root_lifting import PrimePowerRootResult, all_roots_mod_prime_power
import search_mestre_02557104116148_direct_rational as direct
import search_mestre_rank14_pair_rational_frontier as engine
from search_mestre_root_tuple_scale import capped_gp, sha256_file
from search_mestre_root_tuple_scale_max100 import stable_json_digest


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

ROOTS = direct.ROOTS
FAMILY = direct.FAMILY
FAMILY_LABEL = direct.FAMILY_LABEL
ANCHOR = Q(62, 35)
TARGET_LOG_CONDUCTOR = Decimal("182.72")
TARGET_EXPONENTS = ((2, 4), (3, 8), (5, 7), (7, 2), (11, 3), (13, 4), (19, 3))
COEFFICIENT_RADIUS = 4
HEIGHT_CAP = 500_000
FINITE_REDUCTION_TRIGGER = 15
STACK_BYTES = 512_000_000
SELECTION_QUOTAS = {
    "smallest-exact-radical-upper-bound": 64,
    "largest-exact-known-powerful-part": 48,
    "lowest-projective-height": 24,
    "sparse-projective-branch": 24,
}
POINT_STAGES = (
    ("H5000", 5_000, None),
    ("H50000", 50_000, 32),
    ("H250000", 250_000, 8),
    ("H1000000", 1_000_000, 2),
)
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/elliptic_mestre_02557104116148_power_root_crt.json"
)

EXPECTED_DIRECT_SCRIPT_SHA256 = (
    "84fc344ed54f69b7a9d08e0635b08dac242014d87a1a5b5848c5b728f92ab05a"
)
FROZEN_DIRECT_ARTIFACT = "elliptic_mestre_02557104116148_direct_rational.json"
EXPECTED_DIRECT_ARTIFACT_SHA256 = (
    "4874478c553c81ed69fffb49738b5975900a26a17d96f4dca9203a8244e75db6"
)
EXPECTED_DIRECT_RESULT_SHA256 = (
    "c8d231506669c58fbb74b7e9d19b742a15564881f72b12c02d42fc9b1dadb687"
)
FROZEN_NEIGHBORHOOD_ARTIFACT = (
    "elliptic_mestre_02557104116148_t62_35_neighborhood.json"
)
EXPECTED_NEIGHBORHOOD_ARTIFACT_SHA256 = (
    "7953c61b51b27e07840d2dbaaa6016b77d3f443c8b519fb746860f56b6b67af1"
)
EXPECTED_NEIGHBORHOOD_RESULT_SHA256 = (
    "3338871a5a3a7d068e560e6596a3aeb0f10496c27d50bebd061ef806df5fa8c1"
)

ANCHOR_HOMOGENEOUS_DISCRIMINANT = (
    75_885_862_903_905_294_360_264_364_834_669_303_493_203_945_234_915_387_804_237_001_188_014_423_750_000
)
ANCHOR_FACTORIZATION = (
    (2, 4),
    (3, 8),
    (5, 7),
    (7, 2),
    (11, 3),
    (13, 4),
    (19, 3),
    (13_928_891, 1),
    (29_282_567, 1),
    (177_560_677_333_566_672_196_854_931_654_377_865_260_7, 1),
)

# Counts include all roots in the affine chart.  Infinity counts retain only
# u=b/a divisible by p, making the two projective charts a disjoint partition
# of primitive (a,b) modulo every p^k.
EXPECTED_LOCAL_PROFILES = {
    2: {
        "affine_counts": (2, 4, 4, 8),
        "affine_balls": ((1, 0),),
        "infinity_counts": (1, 2, 4, 8),
        "infinity_balls": ((1, 0),),
    },
    3: {
        "affine_counts": (3, 6, 18, 54, 162, 486, 1458, 4374),
        "affine_balls": ((1, 1), (1, 2)),
        "infinity_counts": (1, 3, 0, 0, 0, 0, 0, 0),
        "infinity_balls": (),
    },
    5: {
        "affine_counts": (5, 25, 25, 125, 375, 875, 1500),
        "affine_balls": ((3, 25), (3, 100), (4, 160), (4, 465)),
        "infinity_counts": (1, 5, 25, 125, 250, 1250, 2500),
        "infinity_balls": ((3, 45), (3, 55), (3, 70), (3, 80)),
    },
    7: {
        "affine_counts": (3, 21),
        "affine_balls": ((1, 0), (1, 3), (1, 4)),
        "infinity_counts": (1, 7),
        "infinity_balls": ((1, 0),),
    },
    11: {
        "affine_counts": (3, 33, 242),
        "affine_balls": ((1, 2), (1, 9)),
        "infinity_counts": (1, 11, 22),
        "infinity_balls": ((2, 22), (2, 99)),
    },
    13: {
        "affine_counts": (6, 78, 1014, 6422),
        "affine_balls": (
            (1, 4), (1, 9), (2, 19), (2, 32), (2, 36), (2, 49),
            (2, 68), (2, 84), (2, 85), (2, 101), (2, 120),
            (2, 133), (2, 137), (2, 150),
        ),
        "infinity_counts": (0, 0, 0, 0),
        "infinity_balls": (),
    },
    19: {
        "affine_counts": (9, 171, 722),
        "affine_balls": ((1, 8), (1, 11)),
        "infinity_counts": (0, 0, 0),
        "infinity_balls": (),
    },
}


@dataclass(frozen=True, order=True)
class ProjectiveBall:
    prime: int
    target_exponent: int
    chart: str
    exponent: int
    residue: int

    @property
    def modulus(self) -> int:
        return self.prime**self.exponent

    @property
    def label(self) -> str:
        variable = "t" if self.chart == "affine" else "u"
        return f"p{self.prime}:{self.chart}:{variable}={self.residue}(mod {self.modulus})"

    @property
    def row(self) -> tuple[int, int]:
        # C*a + D*b == 0 modulo p^j.
        if self.chart == "affine":
            return 1, -self.residue
        if self.chart == "infinity":
            return -self.residue, 1
        raise AssertionError("unknown projective chart")


def sequence_digest(values: Iterable[Any]) -> str:
    return hashlib.sha256("\n".join(map(str, values)).encode()).hexdigest()


def pair_digest(values: Iterable[tuple[int, int]]) -> str:
    return hashlib.sha256(
        "".join(f"{left}\t{right}\n" for left, right in values).encode()
    ).hexdigest()


def homogeneous_value(
    coefficients: Sequence[int], numerator: int, denominator: int
) -> int:
    degree = len(coefficients) - 1
    return sum(
        coefficient * numerator**power * denominator ** (degree - power)
        for power, coefficient in enumerate(coefficients)
    )


def integer_valuation(value: int, prime: int) -> int:
    if value == 0:
        raise ValueError("zero has no finite p-adic valuation")
    value = abs(value)
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def derive_discriminant(
    construction: SixRootMestreConstruction,
) -> tuple[tuple[int, ...], int]:
    raw = construction.primitive_discriminant_polynomial
    if len(raw) != 21 or any(value.denominator != 1 for value in raw):
        raise AssertionError("the primitive discriminant lost integral degree 20")
    content = engine.polynomial_content(raw)
    core = tuple(value.numerator // content for value in raw)
    if content != 11_664 or gcd(*(abs(value) for value in core if value)) != 1:
        raise AssertionError("the content-free discriminant normalization changed")
    for parameter in (Q(1), ANCHOR, Q(581, 58), Q(-17, 3)):
        if construction.primitive_jacobian_coefficients(parameter) != direct.family_coefficients(
            parameter
        ):
            raise AssertionError("the pinned short-family formula changed")
    return core, content


def exact_anchor_factorization(
    core: Sequence[int], *, timeout: float, stack_bytes: int
) -> dict[str, Any]:
    value = homogeneous_value(core, ANCHOR.numerator, ANCHOR.denominator)
    if value != ANCHOR_HOMOGENEOUS_DISCRIMINANT:
        raise AssertionError("the exact homogeneous discriminant at 62/35 changed")
    product_value = 1
    for prime, exponent in ANCHOR_FACTORIZATION:
        product_value *= prime**exponent
    if product_value != value:
        raise AssertionError("the pinned anchor factorization product is incomplete")
    primes = [prime for prime, _ in ANCHOR_FACTORIZATION]
    program = "\n".join(
        (
            f"P=[{','.join(map(str, primes))}];",
            'print("PRIME_BEGIN");',
            "for(i=1,#P,print(P[i],\" \",isprime(P[i])));",
            'print("PRIME_END");',
            "quit",
        )
    ) + "\n"
    stdout, _ = capped_gp(program, timeout=timeout, stack_bytes=stack_bytes)
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    start = lines.index("PRIME_BEGIN") + 1
    end = lines.index("PRIME_END")
    primality = []
    for line in lines[start:end]:
        prime_text, flag_text = line.split()
        primality.append((int(prime_text), int(flag_text)))
    if primality != [(prime, 1) for prime in primes]:
        raise AssertionError("PARI did not certify every factor base prime")
    return {
        "parameter_T": str(ANCHOR),
        "absolute_content_free_homogeneous_discriminant": str(value),
        "factorization": [
            {"prime": str(prime), "exponent": exponent, "PARI_isprime": True}
            for prime, exponent in ANCHOR_FACTORIZATION
        ],
        "exact_product_replay": True,
        "all_factor_bases_certified_prime_by_one_capped_PARi_batch": True,
    }


def filtered_infinity_result(
    reversed_core: Sequence[int], prime: int, exponent: int
) -> tuple[PrimePowerRootResult, tuple[int, ...]]:
    histories = []
    final = None
    for level in range(1, exponent + 1):
        result = all_roots_mod_prime_power(
            reversed_core, prime, level, max_roots=100_000
        )
        histories.append(sum(root % prime == 0 for root in result.roots))
        final = result
    assert final is not None
    roots = tuple(root for root in final.roots if root % prime == 0)
    return (
        PrimePowerRootResult(
            prime, exponent, prime**exponent, roots, tuple(histories),
            final.candidate_digits_checked,
        ),
        tuple(histories),
    )


def expand_balls(
    balls: Sequence[ProjectiveBall], target_exponent: int
) -> tuple[int, ...]:
    covered: set[int] = set()
    for ball in balls:
        covered.update(
            ball.residue + digit * ball.modulus
            for digit in range(ball.prime ** (target_exponent - ball.exponent))
        )
    return tuple(sorted(covered))


def discover_projective_root_balls(
    core: Sequence[int],
) -> tuple[dict[int, tuple[ProjectiveBall, ...]], list[dict[str, Any]]]:
    groups: dict[int, tuple[ProjectiveBall, ...]] = {}
    profiles = []
    reversed_core = tuple(reversed(core))
    for prime, exponent in TARGET_EXPONENTS:
        affine = all_roots_mod_prime_power(core, prime, exponent, max_roots=100_000)
        infinity, infinity_history = filtered_infinity_result(
            reversed_core, prime, exponent
        )
        affine_balls = tuple(
            ProjectiveBall(prime, exponent, "affine", ball.exponent, ball.residue)
            for ball in affine.maximal_balls()
        )
        infinity_balls = tuple(
            ProjectiveBall(prime, exponent, "infinity", ball.exponent, ball.residue)
            for ball in infinity.maximal_balls()
        )
        expected = EXPECTED_LOCAL_PROFILES[prime]
        observed = {
            "affine_counts": tuple(affine.level_counts),
            "affine_balls": tuple((ball.exponent, ball.residue) for ball in affine_balls),
            "infinity_counts": infinity_history,
            "infinity_balls": tuple(
                (ball.exponent, ball.residue) for ball in infinity_balls
            ),
        }
        if observed != expected:
            raise AssertionError(f"the complete projective p={prime} lift tree changed")
        if expand_balls(affine_balls, exponent) != affine.roots:
            raise AssertionError("affine maximal-ball coverage replay failed")
        if expand_balls(infinity_balls, exponent) != infinity.roots:
            raise AssertionError("infinity maximal-ball coverage replay failed")
        combined = tuple(sorted((*affine_balls, *infinity_balls)))
        if not combined:
            raise AssertionError(f"p={prime} lost every target valuation branch")
        groups[prime] = combined
        profiles.append(
            {
                "prime": prime,
                "target_valuation": exponent,
                "charts_partition": (
                    "affine b-unit t=a/b; infinity a-unit u=b/a with u=0 mod p"
                ),
                "affine": {
                    "level_root_counts": list(affine.level_counts),
                    "target_roots": len(affine.roots),
                    "target_root_sha256": sequence_digest(affine.roots),
                    "maximal_balls": [ball_record(ball) for ball in affine_balls],
                },
                "infinity": {
                    "level_root_counts_after_u_divisible_by_p_filter": list(
                        infinity_history
                    ),
                    "target_roots": len(infinity.roots),
                    "target_root_sha256": sequence_digest(infinity.roots),
                    "maximal_balls": [ball_record(ball) for ball in infinity_balls],
                },
                "complete_target_root_count": len(affine.roots) + len(infinity.roots),
                "complete_maximal_ball_count": len(combined),
                "coverage_replayed_exactly": True,
            }
        )
    if tuple(groups) != tuple(prime for prime, _ in TARGET_EXPONENTS):
        raise AssertionError("a target prime escaped projective profiling")
    return groups, profiles


def ball_record(ball: ProjectiveBall) -> dict[str, Any]:
    return {
        "label": ball.label,
        "chart": ball.chart,
        "residue": ball.residue,
        "modulus_exponent": ball.exponent,
        "modulus": ball.modulus,
        "forced_discriminant_valuation_at_least": ball.target_exponent,
        "target_leaf_count_covered": ball.prime ** (
            ball.target_exponent - ball.exponent
        ),
        "homogeneous_linear_row_C_D": list(ball.row),
    }


def combine_rows(choices: Sequence[ProjectiveBall]) -> tuple[int, int, int]:
    coefficient_a, coefficient_b, modulus = 0, 0, 1
    for choice in choices:
        row_a, row_b = choice.row
        next_a, next_modulus = crt_pair(
            coefficient_a, modulus, row_a, choice.modulus
        )
        next_b, check_modulus = crt_pair(
            coefficient_b, modulus, row_b, choice.modulus
        )
        if next_modulus != check_modulus:
            raise AssertionError("CRT row moduli diverged")
        coefficient_a, coefficient_b, modulus = next_a, next_b, next_modulus
    if gcd(gcd(coefficient_a, coefficient_b), modulus) != 1:
        raise AssertionError("a combined projective row ceased to be primitive")
    return coefficient_a, coefficient_b, modulus


def kernel_basis(coefficient_a: int, coefficient_b: int, modulus: int):
    divisor = gcd(coefficient_a, modulus)
    if gcd(divisor, coefficient_b) != 1:
        raise AssertionError("the homogeneous kernel lost its unit complement")
    quotient = modulus // divisor
    residue = 0
    if quotient > 1:
        residue = (
            -coefficient_b * pow(coefficient_a // divisor, -1, quotient)
        ) % quotient
    first, second = gauss_reduce((quotient, 0), (residue, divisor))
    determinant = abs(first[0] * second[1] - first[1] * second[0])
    if determinant != modulus:
        raise AssertionError("Gauss reduction changed the kernel determinant")
    for vector in (first, second):
        if (coefficient_a * vector[0] + coefficient_b * vector[1]) % modulus:
            raise AssertionError("a reduced basis vector left the CRT kernel")
    return first, second


def prior_parameters(generated: Path) -> tuple[set[Fraction], dict[str, Any]]:
    answer: set[Fraction] = set()
    records = []
    for name, expected_sha in sorted(direct.PRIOR_ARTIFACT_SHA256.items()):
        path = generated / name
        actual_sha = sha256_file(path)
        if actual_sha != expected_sha:
            raise AssertionError(f"frozen prior parameter artifact changed: {name}")
        local: set[Fraction] = set()
        direct._collect_parameter_t(json.loads(path.read_text()), local)
        answer.update(local)
        records.append(
            {"artifact": name, "sha256": actual_sha, "unique_absolute_T": len(local)}
        )
    if len(answer) != direct.EXPECTED_PRIOR_PARAMETER_T_ABSOLUTE_COUNT:
        raise AssertionError("the all-height frozen prior parameter census changed")
    ordered = sorted(answer)
    return answer, {
        "method": "recursive exact Fraction extraction of every parameter_t",
        "source_artifacts": records,
        "unique_absolute_parameter_T": len(answer),
        "canonical_parameter_lines_sha256": sequence_digest(ordered),
    }


def in_direct_box(numerator: int, denominator: int) -> bool:
    return numerator <= 30_000 and denominator <= 1_000


def in_completed_annulus(numerator: int, denominator: int) -> bool:
    return (
        1_001 <= denominator <= 5_000
        and 7 * numerator >= 11 * denominator
        and 35 * numerator <= 69 * denominator
    )


def matching_ball(
    numerator: int, denominator: int, balls: Sequence[ProjectiveBall]
) -> ProjectiveBall:
    prime = balls[0].prime
    chart = "affine" if denominator % prime else "infinity"
    matches = []
    for ball in balls:
        if ball.chart != chart:
            continue
        if chart == "affine":
            value = numerator * pow(denominator, -1, ball.modulus) % ball.modulus
        else:
            value = denominator * pow(numerator, -1, ball.modulus) % ball.modulus
        if value == ball.residue:
            matches.append(ball)
    if len(matches) != 1:
        raise AssertionError(
            f"a primitive candidate matched {len(matches)} p={prime} maximal balls"
        )
    return matches[0]


def enumerate_population(
    groups: dict[int, tuple[ProjectiveBall, ...]],
    prior: set[Fraction],
    *,
    coefficient_radius: int,
    height_cap: int,
) -> tuple[list[tuple[int, int]], dict[str, Any]]:
    ordered_groups = [groups[prime] for prime, _ in TARGET_EXPONENTS]
    branch_count = 1
    for group in ordered_groups:
        branch_count *= len(group)
    if branch_count != 14_336:
        raise AssertionError("the complete projective branch product changed")

    representatives: set[tuple[int, int]] = set()
    counters = Counter()
    class_digest = hashlib.sha256()
    raw_retained = 0
    for class_index, choices in enumerate(product(*ordered_groups)):
        coefficient_a, coefficient_b, modulus = combine_rows(choices)
        class_digest.update(
            (
                f"{class_index}|{coefficient_a}|{coefficient_b}|{modulus}|"
                f"{'/'.join(choice.label for choice in choices)}\n"
            ).encode()
        )
        first, second = kernel_basis(coefficient_a, coefficient_b, modulus)
        for left in range(-coefficient_radius, coefficient_radius + 1):
            for right in range(-coefficient_radius, coefficient_radius + 1):
                if left == 0 and right == 0:
                    continue
                counters["bounded_vectors_visited"] += 1
                numerator = left * first[0] + right * second[0]
                denominator = left * first[1] + right * second[1]
                if denominator == 0:
                    counters["zero_denominator"] += 1
                    continue
                common = gcd(abs(numerator), abs(denominator))
                if gcd(common, modulus) != 1:
                    counters["nonunit_common_divisor"] += 1
                    continue
                numerator //= common
                denominator //= common
                if denominator < 0:
                    numerator = -numerator
                    denominator = -denominator
                numerator = abs(numerator)
                if numerator == 0 or gcd(numerator, denominator) != 1:
                    counters["zero_or_nonprimitive_after_quotient"] += 1
                    continue
                if max(numerator, denominator) > height_cap:
                    counters["above_height_cap"] += 1
                    continue
                if in_direct_box(numerator, denominator):
                    counters["inside_completed_direct_box"] += 1
                    continue
                if in_completed_annulus(numerator, denominator):
                    counters["inside_completed_farey_annulus"] += 1
                    continue
                parameter = Q(numerator, denominator)
                if parameter in prior:
                    counters["frozen_prior_parameter_T"] += 1
                    continue
                raw_retained += 1
                representatives.add((numerator, denominator))

    ordered = sorted(
        representatives, key=lambda pair: (max(pair), pair[1], pair[0])
    )
    if not ordered:
        raise AssertionError("the p-adic lattice population became empty")
    return ordered, {
        "complete_projective_branch_combinations": branch_count,
        "branch_class_rows_sha256": class_digest.hexdigest(),
        "gauss_coefficient_radius": coefficient_radius,
        "height_cap": height_cap,
        "bounded_vectors_visited": counters["bounded_vectors_visited"],
        "zero_denominator": counters["zero_denominator"],
        "nonunit_common_divisor": counters["nonunit_common_divisor"],
        "zero_or_nonprimitive_after_quotient": counters[
            "zero_or_nonprimitive_after_quotient"
        ],
        "above_height_cap": counters["above_height_cap"],
        "inside_completed_direct_box": counters["inside_completed_direct_box"],
        "inside_completed_farey_annulus": counters[
            "inside_completed_farey_annulus"
        ],
        "frozen_prior_parameter_T": counters["frozen_prior_parameter_T"],
        "raw_retained_vector_instances": raw_retained,
        "unique_reduced_positive_parameters": len(ordered),
        "candidate_pair_lines_sha256": pair_digest(ordered),
        "minimum_projective_height": max(ordered[0]),
        "maximum_projective_height": max(max(pair) for pair in ordered),
        "completed_direct_box_excluded_exactly": True,
        "completed_farey_annulus_excluded_exactly": True,
        "T_sign_quotient": "absolute numerator because A and B are even",
    }


def feature_population(
    population: Sequence[tuple[int, int]],
    groups: dict[int, tuple[ProjectiveBall, ...]],
    core: Sequence[int],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    records = []
    signature_counts = Counter()
    digest = hashlib.sha256()
    targets = dict(TARGET_EXPONENTS)
    for index, (numerator, denominator) in enumerate(population, 1):
        value = homogeneous_value(core, numerator, denominator)
        if value == 0:
            raise AssertionError("a singular parameter entered the p-adic population")
        balls = tuple(
            matching_ball(numerator, denominator, groups[prime])
            for prime, _ in TARGET_EXPONENTS
        )
        valuations = {
            prime: integer_valuation(value, prime) for prime, _ in TARGET_EXPONENTS
        }
        if any(valuations[prime] < targets[prime] for prime in targets):
            raise AssertionError("a CRT candidate lost an exact forced valuation")
        signature = tuple(ball.label for ball in balls)
        signature_text = "|".join(signature)
        signature_counts[signature_text] += 1
        feature = engine.discriminant_feature(core, numerator, denominator)
        record = {
            "family_index": 0,
            "family_label": FAMILY_LABEL,
            "numerator": numerator,
            "denominator": denominator,
            "parameter": str(Q(numerator, denominator)),
            "projective_height": max(numerator, denominator),
            "projective_branch_signature": list(signature),
            "exact_forced_prime_valuations": {
                str(prime): valuations[prime] for prime, _ in TARGET_EXPONENTS
            },
            "discriminant_feature": feature,
        }
        records.append(record)
        digest.update(
            (
                f"{record['parameter']}|{signature_text}|"
                f"{feature['absolute_homogeneous_discriminant']}|"
                f"{feature['combined_radical_upper_bound']}|"
                f"{feature['known_powerful_part']}\n"
            ).encode()
        )
        if index % 5_000 == 0:
            print(f"exact radical features {index}/{len(population)}", flush=True)
    for record in records:
        signature_text = "|".join(record["projective_branch_signature"])
        record["projective_branch_population"] = signature_counts[signature_text]
    return records, {
        "population": len(records),
        "exact_singular_rejections": 0,
        "trial_division_prime_bound": engine.TRIAL_DIVISION_LIMIT,
        "distinct_projective_branch_signatures": len(signature_counts),
        "minimum_branch_population": min(signature_counts.values()),
        "maximum_branch_population": max(signature_counts.values()),
        "exact_feature_population_sha256": digest.hexdigest(),
    }


def select_conductor_population(
    records: Sequence[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    orders = {
        "smallest-exact-radical-upper-bound": sorted(
            records,
            key=lambda row: (
                int(row["discriminant_feature"]["combined_radical_upper_bound"]),
                row["projective_height"], row["denominator"], row["numerator"],
            ),
        ),
        "largest-exact-known-powerful-part": sorted(
            records,
            key=lambda row: (
                -int(row["discriminant_feature"]["known_powerful_part"]),
                int(row["discriminant_feature"]["combined_radical_upper_bound"]),
                row["projective_height"], row["denominator"], row["numerator"],
            ),
        ),
        "lowest-projective-height": sorted(
            records,
            key=lambda row: (
                row["projective_height"], row["denominator"], row["numerator"],
            ),
        ),
    }
    best_by_signature: dict[tuple[str, ...], dict[str, Any]] = {}
    for row in orders["smallest-exact-radical-upper-bound"]:
        signature = tuple(row["projective_branch_signature"])
        best_by_signature.setdefault(signature, row)
    orders["sparse-projective-branch"] = sorted(
        best_by_signature.values(),
        key=lambda row: (
            row["projective_branch_population"],
            int(row["discriminant_feature"]["combined_radical_upper_bound"]),
            row["projective_height"], row["denominator"], row["numerator"],
        ),
    )

    selected: dict[str, dict[str, Any]] = {}
    reasons: dict[str, set[str]] = {}
    for label, quota in SELECTION_QUOTAS.items():
        added = 0
        for row in orders[label]:
            parameter = row["parameter"]
            if parameter in selected:
                reasons[parameter].add(label)
                continue
            selected[parameter] = row
            reasons[parameter] = {label}
            added += 1
            if added == quota:
                break
        if added != quota:
            raise AssertionError(f"the fixed conductor quota failed: {label}")
    answer = []
    for parameter, row in selected.items():
        copy = dict(row)
        copy["conductor_selection_strata"] = sorted(reasons[parameter])
        answer.append(copy)
    answer.sort(key=lambda row: (row["projective_height"], row["denominator"], row["numerator"]))
    digest = hashlib.sha256()
    for row in answer:
        digest.update(
            (
                f"{row['parameter']}|{','.join(row['conductor_selection_strata'])}|"
                f"{'/'.join(row['projective_branch_signature'])}\n"
            ).encode()
        )
    return answer, {
        "population_closed_before_any_conductor_point_or_rank_call": True,
        "selection_uses_conductor_point_or_rank_data": False,
        "exact_radical_features_closed_before_selection": True,
        "novel_candidates_added_per_stratum": SELECTION_QUOTAS,
        "selected_population": len(answer),
        "selected_population_sha256": digest.hexdigest(),
    }


def stage_rank_key(record: dict[str, Any], prior_stage: str) -> tuple[Any, ...]:
    stage = record["point_stages"][prior_stage]
    conductor = record["conductor_phase"]
    return (
        -int(stage["stable_numerical_rank"]),
        not conductor.get("below_strict_log_conductor_target_numerically", False),
        Decimal(conductor["log_conductor"]),
        int(record["discriminant_feature"]["combined_radical_upper_bound"]),
        record["projective_height"], record["denominator"], record["numerator"],
    )


def frozen_scope_audit(generated: Path) -> dict[str, Any]:
    direct_path = generated / FROZEN_DIRECT_ARTIFACT
    neighborhood_path = generated / FROZEN_NEIGHBORHOOD_ARTIFACT
    if sha256_file(direct_path) != EXPECTED_DIRECT_ARTIFACT_SHA256:
        raise AssertionError("the completed direct-box artifact changed")
    if sha256_file(neighborhood_path) != EXPECTED_NEIGHBORHOOD_ARTIFACT_SHA256:
        raise AssertionError("the completed Farey-annulus artifact changed")
    direct_data = json.loads(direct_path.read_text())
    neighborhood_data = json.loads(neighborhood_path.read_text())
    if direct_data["result_sha256"] != EXPECTED_DIRECT_RESULT_SHA256:
        raise AssertionError("the direct-box result digest changed")
    if neighborhood_data["result_sha256"] != EXPECTED_NEIGHBORHOOD_RESULT_SHA256:
        raise AssertionError("the Farey-annulus result digest changed")
    return {
        "completed_direct_box": {
            "artifact": FROZEN_DIRECT_ARTIFACT,
            "sha256": EXPECTED_DIRECT_ARTIFACT_SHA256,
            "numerator_interval": [1, 30_000],
            "denominator_interval": [1, 1_000],
        },
        "completed_farey_annulus": {
            "artifact": FROZEN_NEIGHBORHOOD_ARTIFACT,
            "sha256": EXPECTED_NEIGHBORHOOD_ARTIFACT_SHA256,
            "denominator_interval": [1_001, 5_000],
            "parameter_interval": ["11/7", "69/35"],
        },
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection-only", action="store_true")
    parser.add_argument("--coefficient-radius", type=int, default=COEFFICIENT_RADIUS)
    parser.add_argument("--height-cap", type=int, default=HEIGHT_CAP)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--primality-timeout", type=float, default=20.0)
    parser.add_argument("--conductor-timeout", type=float, default=12.0)
    parser.add_argument("--h5000-timeout", type=float, default=15.0)
    parser.add_argument("--h50000-timeout", type=float, default=20.0)
    parser.add_argument("--h250000-timeout", type=float, default=30.0)
    parser.add_argument("--h1000000-timeout", type=float, default=45.0)
    parser.add_argument("--height-timeout", type=float, default=25.0)
    parser.add_argument("--ellrank-timeout", type=float, default=12.0)
    parser.add_argument("--mapping-cap", type=int, default=512)
    parser.add_argument("--certificate-prime-bound", type=int, default=1000)
    parser.add_argument("--stack-bytes", type=int, default=STACK_BYTES)
    parser.add_argument("--output", type=Path, default=root / DEFAULT_OUTPUT)
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.coefficient_radius != COEFFICIENT_RADIUS or args.height_cap != HEIGHT_CAP:
        raise SystemExit("the Gauss radius and height cap are pinned at 4 and 500000")
    if not 1 <= args.workers <= 8:
        raise SystemExit("workers must lie in [1,8]")
    timeouts = (
        args.primality_timeout, args.conductor_timeout, args.h5000_timeout,
        args.h50000_timeout, args.h250000_timeout, args.h1000000_timeout,
        args.height_timeout, args.ellrank_timeout,
    )
    if min(timeouts) <= 0 or max(timeouts) > 60:
        raise SystemExit("all foreground subprocess caps must lie in (0,60]")
    if not 32 <= args.mapping_cap <= 1024:
        raise SystemExit("mapping cap must lie in [32,1024]")
    if not 499 <= args.certificate_prime_bound <= 2000:
        raise SystemExit("certificate prime bound must lie in [499,2000]")
    if not 32_000_000 <= args.stack_bytes <= 1_000_000_000:
        raise SystemExit("stack bytes must lie in [32000000,1000000000]")


def exclusive_write(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")


def main() -> None:
    args = build_parser().parse_args()
    validate_args(args)
    if args.output.exists():
        raise SystemExit("refusing to overwrite the p-adic power-root artifact")
    started = time.monotonic()
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    generated = root / "artifacts/generated-results"
    direct_script = script_path.with_name(
        "search_mestre_02557104116148_direct_rational.py"
    )
    if sha256_file(direct_script) != EXPECTED_DIRECT_SCRIPT_SHA256:
        raise AssertionError("the frozen direct-family API changed")
    frozen_scopes = frozen_scope_audit(generated)
    prior, prior_audit = prior_parameters(generated)

    construction = SixRootMestreConstruction(tuple(Q(root_value) for root_value in ROOTS))
    if construction.quartic_condition or construction.is_reflection_symmetric:
        raise AssertionError("the selected family geometry changed")
    core, content = derive_discriminant(construction)
    anchor_factorization = exact_anchor_factorization(
        core, timeout=args.primality_timeout, stack_bytes=args.stack_bytes
    )
    groups, profiles = discover_projective_root_balls(core)
    population, population_audit = enumerate_population(
        groups, prior,
        coefficient_radius=args.coefficient_radius,
        height_cap=args.height_cap,
    )
    print(
        "complete p-adic CRT/Gauss population closed: "
        f"branches={population_audit['complete_projective_branch_combinations']} "
        f"unique={len(population)}",
        flush=True,
    )
    records, feature_audit = feature_population(population, groups, core)
    selected, selection = select_conductor_population(records)
    print(
        f"exact radical population closed; conductor selection={len(selected)}",
        flush=True,
    )

    direct.configure_engine(prior)
    engine.FINITE_REDUCTION_TRIGGER = FINITE_REDUCTION_TRIGGER
    common: dict[str, Any] = {
        "schema_version": 1,
        "status": (
            "completed selection-only projective p-adic power-root CRT scan"
            if args.selection_only
            else "in-progress conductor-first projective p-adic power-root CRT scan"
        ),
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "alternative_rank_at_least": 30,
            "hits": [],
        },
        "scope": {
            "family_roots": list(ROOTS),
            "anchor_T": str(ANCHOR),
            "direct_parameter": "positive reduced T=a/b modulo exact T<->-T symmetry",
            "orthogonal_to_prior_searches": frozen_scopes,
            "coefficient_radius": args.coefficient_radius,
            "projective_height_cap": args.height_cap,
        },
        "family": {
            "label": FAMILY_LABEL,
            "A_coefficients_ascending": list(direct.A_COEFFICIENTS),
            "B_coefficients_ascending": list(direct.B_COEFFICIENTS),
            "content_free_discriminant_coefficients_ascending": list(core),
            "removed_discriminant_polynomial_content": content,
        },
        "anchor_exact_factorization": anchor_factorization,
        "complete_projective_p_adic_profiles": profiles,
        "crt_branch_product": {
            "one_maximal_ball_chosen_per_prime": True,
            "prime_order": [prime for prime, _ in TARGET_EXPONENTS],
            "ball_counts": {str(prime): len(groups[prime]) for prime, _ in TARGET_EXPONENTS},
            "complete_combination_count": population_audit[
                "complete_projective_branch_combinations"
            ],
            "homogeneous_linear_Ca_plus_Db_congruences": True,
            "exact_two_dimensional_Gauss_reduction": True,
        },
        "prior_parameter_exclusion_audit": prior_audit,
        "lattice_population_audit": population_audit,
        "exact_radical_feature_screen": feature_audit,
        "conductor_selection": selection,
        "selected_records": selected,
        "parameters": {
            key: value for key, value in vars(args).items()
            if key not in {"output", "selection_only"}
        },
        "provenance": {
            "script_path": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "frozen_direct_family_script_sha256": EXPECTED_DIRECT_SCRIPT_SHA256,
            "reproducing_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "external_calls_use_foreground_process_groups": True,
            "same_stage_retries": 0,
        },
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "timings": {"pre_conductor_wall_seconds": time.monotonic() - started},
    }
    if args.selection_only:
        common["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
        common["result_sha256"] = stable_json_digest(
            {
                "scope": common["scope"],
                "family": common["family"],
                "factorization": anchor_factorization,
                "profiles": profiles,
                "prior": prior_audit,
                "population": population_audit,
                "features": feature_audit,
                "selection": selection,
                "records": selected,
            }
        )
        exclusive_write(args.output, common)
        return

    conductor_started = time.monotonic()
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(
                engine.conductor_worker,
                0, row["numerator"], row["denominator"],
                args.conductor_timeout, args.stack_bytes,
            )
            for row in selected
        ]
        for position, (row, future) in enumerate(zip(selected, futures), start=1):
            row["conductor_phase"] = future.result()
            if position % 20 == 0:
                print(f"conductors {position}/{len(selected)}", flush=True)
    eligible = [
        row for row in selected
        if row["conductor_phase"]["status"].startswith("completed")
    ]

    stage_timeouts = {
        "H5000": args.h5000_timeout,
        "H50000": args.h50000_timeout,
        "H250000": args.h250000_timeout,
        "H1000000": args.h1000000_timeout,
    }
    current = eligible
    for stage_index, (name, height, keep) in enumerate(POINT_STAGES):
        if stage_index:
            prior_name = POINT_STAGES[stage_index - 1][0]
            current = [
                row for row in current
                if row.get("point_stages", {}).get(prior_name, {}).get("status")
                == "completed"
            ]
            current.sort(key=lambda row: stage_rank_key(row, prior_name))
            current = current[: int(keep)]
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = [
                executor.submit(
                    engine.point_stage_worker,
                    0, row["numerator"], row["denominator"],
                    height, stage_timeouts[name], args.height_timeout,
                    args.ellrank_timeout, args.stack_bytes, args.mapping_cap,
                    args.certificate_prime_bound,
                )
                for row in current
            ]
            for row, future in zip(current, futures):
                row.setdefault("point_stages", {})[name] = future.result()
        maximum = max(
            (
                row["point_stages"][name].get("stable_numerical_rank", -1)
                for row in current
            ),
            default=-1,
        )
        print(f"{name} attempted={len(current)} max_rank={maximum}", flush=True)
        if maximum >= FINITE_REDUCTION_TRIGGER:
            leaders = [
                row["parameter"] for row in current
                if row["point_stages"][name].get("stable_numerical_rank", -1)
                == maximum
            ]
            print(f"EARLY_SIGNAL {name} stable_rank={maximum} T={leaders}", flush=True)

    completed_stages = [
        stage
        for row in selected
        for stage in row.get("point_stages", {}).values()
        if stage.get("status") == "completed"
    ]
    maximum_rank = max(
        (int(stage["stable_numerical_rank"]) for stage in completed_stages),
        default=None,
    )
    finite_attempts = []
    target_hits = []
    for row in selected:
        for stage_name, stage in row.get("point_stages", {}).items():
            certificate = stage.get("finite_reduction_attempt", {})
            certified = certificate.get("certified_algebraic_rank_lower_bound")
            if certified is None:
                continue
            finite_attempts.append(
                {
                    "parameter": row["parameter"],
                    "stage": stage_name,
                    "certified_rank_lower_bound": certified,
                    "point_sha256": certificate.get("point_sha256"),
                }
            )
            below = row["conductor_phase"].get(
                "below_strict_log_conductor_target_numerically", False
            )
            if certified >= 30 or (certified >= 21 and below):
                target_hits.append(
                    {
                        "parameter": row["parameter"],
                        "stage": stage_name,
                        "certified_rank_lower_bound": certified,
                        "conductor": row["conductor_phase"]["conductor"],
                        "log_conductor": row["conductor_phase"]["log_conductor"],
                    }
                )

    common["status"] = "completed fixed projective p-adic power-root CRT scan"
    common["target"]["hits"] = target_hits
    common["conductor_first_screen"] = {
        "population_closed_before_any_conductor_point_or_rank_call": True,
        "selected_population": len(selected),
        "completed": len(eligible),
        "timeouts": sum(
            row["conductor_phase"]["status"] == "timeout" for row in selected
        ),
        "errors": sum(
            row["conductor_phase"]["status"] == "error" for row in selected
        ),
        "subtarget": sum(
            row["conductor_phase"].get(
                "below_strict_log_conductor_target_numerically"
            ) is True
            for row in selected
        ),
    }
    common["point_search_protocol"] = {
        "stages": [
            {
                "name": name,
                "height_bound": height,
                "keep_after_previous_stage": keep,
                "attempted": sum(
                    name in row.get("point_stages", {}) for row in selected
                ),
            }
            for name, height, keep in POINT_STAGES
        ],
        "increasing_heights_are_not_retries": True,
        "same_height_retries": 0,
        "finite_reduction_trigger_stable_rank": FINITE_REDUCTION_TRIGGER,
        "finite_reduction_attempts": finite_attempts,
        "maximum_stable_numerical_rank": maximum_rank,
        "completed_stage_calls": len(completed_stages),
    }
    common["timings"].update(
        {
            "conductor_and_point_wall_seconds": time.monotonic() - conductor_started,
            "total_wall_seconds": time.monotonic() - started,
        }
    )
    common["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    common["provenance"]["owned_processes_remaining"] = 0
    common["result_sha256"] = stable_json_digest(
        {
            "scope": common["scope"],
            "family": common["family"],
            "factorization": anchor_factorization,
            "profiles": profiles,
            "prior": prior_audit,
            "population": population_audit,
            "features": feature_audit,
            "selection": selection,
            "records": selected,
            "conductor": common["conductor_first_screen"],
            "points": common["point_search_protocol"],
            "target": common["target"],
        }
    )
    exclusive_write(args.output, common)
    print(
        f"complete max_rank={maximum_rank} target_hits={len(target_hits)} "
        f"output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
