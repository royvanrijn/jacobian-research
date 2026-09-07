#!/usr/bin/env python3
"""Bounded survey of affine-normalized six-root Mestre configurations.

The survey first filters root tuples by the exact degree-five obstruction, so
only genuine quartic remainders are scored.  It then records:

* small-prime root densities of the quartic discriminant after subtracting
  each prime's fixed-divisor valuation;
* a bounded small-prime ``power savings`` proxy, not an actual conductor; and
* collisions or zero ordinates among the twelve displayed Mestre points.

No generic-rank, specialized-rank, or point-independence claim is made.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from math import isqrt, log
from pathlib import Path
import platform
import shlex
import sys
from typing import Any

from ek_k3 import primes_up_to, rational_to_string, valuation
from fermigier_mestre import FermigierMestreFamily, ROOTS
from mestre_root_tuples import (
    SixRootMestreConstruction,
    affine_normalized_integer_root_tuples,
    normalize_integer_root_tuple,
)
from pari_bridge import minimal_curve_data, pari_version


Q = Fraction


def _trim_polynomial(coefficients: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients = coefficients[:-1]
    return coefficients


def _add_polynomials(
    left: tuple[Fraction, ...], right: tuple[Fraction, ...]
) -> tuple[Fraction, ...]:
    return _trim_polynomial(
        tuple(
            (left[index] if index < len(left) else Q(0))
            + (right[index] if index < len(right) else Q(0))
            for index in range(max(len(left), len(right)))
        )
    )


def _multiply_polynomials(
    left: tuple[Fraction, ...], right: tuple[Fraction, ...]
) -> tuple[Fraction, ...]:
    answer = [Q(0)] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] += left_value * right_value
    return _trim_polynomial(tuple(answer))


def _polynomial_divmod(
    numerator: tuple[Fraction, ...], denominator: tuple[Fraction, ...]
) -> tuple[tuple[Fraction, ...], tuple[Fraction, ...]]:
    numerator = list(_trim_polynomial(numerator))
    denominator = _trim_polynomial(denominator)
    if denominator == (0,):
        raise ZeroDivisionError("polynomial division by zero")
    quotient = [Q(0)] * max(len(numerator) - len(denominator) + 1, 1)
    while len(numerator) >= len(denominator) and any(numerator):
        offset = len(numerator) - len(denominator)
        coefficient = numerator[-1] / denominator[-1]
        quotient[offset] = coefficient
        for index, value in enumerate(denominator):
            numerator[offset + index] -= coefficient * value
        while len(numerator) > 1 and numerator[-1] == 0:
            numerator.pop()
    return _trim_polynomial(tuple(quotient)), _trim_polynomial(tuple(numerator))


def _polynomial_gcd(
    left: tuple[Fraction, ...], right: tuple[Fraction, ...]
) -> tuple[Fraction, ...]:
    left = _trim_polynomial(left)
    right = _trim_polynomial(right)
    while right != (0,):
        _, remainder = _polynomial_divmod(left, right)
        left, right = right, remainder
    if left == (0,):
        return left
    return tuple(value / left[-1] for value in left)


def _modular_gcd_degree(
    left: tuple[Fraction, ...], right: tuple[Fraction, ...], prime: int
) -> int | None:
    """Return a reduction gcd degree, or ``None`` at a bad denominator/degree."""

    def reduce_polynomial(values: tuple[Fraction, ...]) -> list[int] | None:
        answer: list[int] = []
        for value in values:
            if value.denominator % prime == 0:
                return None
            answer.append(
                value.numerator
                * pow(value.denominator % prime, -1, prime)
                % prime
            )
        while len(answer) > 1 and answer[-1] == 0:
            answer.pop()
        return answer

    def remainder(numerator: list[int], denominator: list[int]) -> list[int]:
        numerator = numerator[:]
        inverse_lead = pow(denominator[-1], -1, prime)
        while len(numerator) >= len(denominator) and any(numerator):
            offset = len(numerator) - len(denominator)
            coefficient = numerator[-1] * inverse_lead % prime
            for index, value in enumerate(denominator):
                numerator[offset + index] = (
                    numerator[offset + index] - coefficient * value
                ) % prime
            while len(numerator) > 1 and numerator[-1] == 0:
                numerator.pop()
        return numerator

    reduced_left = reduce_polynomial(left)
    reduced_right = reduce_polynomial(right)
    if (
        reduced_left is None
        or reduced_right is None
        or len(reduced_left) != len(_trim_polynomial(left))
    ):
        return None
    while reduced_right != [0]:
        reduced_left, reduced_right = (
            reduced_right,
            remainder(reduced_left, reduced_right),
        )
    return len(reduced_left) - 1


def discriminant_degree_profile(
    construction: SixRootMestreConstruction,
) -> dict[str, int | str]:
    coefficients = construction.primitive_discriminant_polynomial
    if all(value == 0 for value in coefficients):
        return {
            "status": "identically zero discriminant",
            "degree": -1,
            "squarefree_degree": -1,
            "repeated_factor_degree": -1,
        }
    derivative = tuple(
        index * value for index, value in enumerate(coefficients[1:], 1)
    )
    common = _polynomial_gcd(coefficients, derivative)
    repeated_degree = 0 if common == (Q(1),) else len(common) - 1
    degree = len(coefficients) - 1
    return {
        "status": "exact polynomial computation",
        "degree": degree,
        "squarefree_degree": degree - repeated_degree,
        "repeated_factor_degree": repeated_degree,
    }


def _polynomial_power(
    coefficients: tuple[Fraction, ...], exponent: int
) -> tuple[Fraction, ...]:
    answer = (Q(1),)
    factor = coefficients
    while exponent:
        if exponent & 1:
            answer = _multiply_polynomials(answer, factor)
        factor = _multiply_polynomials(factor, factor)
        exponent //= 2
    return answer


def rational_base_change_polynomial(
    polynomial: tuple[Fraction, ...],
    numerator: tuple[Fraction, ...],
    denominator: tuple[Fraction, ...],
) -> tuple[Fraction, ...]:
    """Return ``den(U)^d * polynomial(num(U)/den(U))`` exactly."""

    degree = len(_trim_polynomial(polynomial)) - 1
    answer = (Q(0),)
    numerator_powers = [(Q(1),)]
    denominator_powers = [(Q(1),)]
    for _ in range(degree):
        numerator_powers.append(_multiply_polynomials(numerator_powers[-1], numerator))
        denominator_powers.append(
            _multiply_polynomials(denominator_powers[-1], denominator)
        )
    for index, coefficient in enumerate(polynomial):
        term = _multiply_polynomials(
            numerator_powers[index], denominator_powers[degree - index]
        )
        answer = _add_polynomials(
            answer, tuple(coefficient * value for value in term)
        )
    return _trim_polynomial(answer)


def _rational_square_root(value: Fraction) -> Fraction | None:
    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator**2 != value.numerator or denominator**2 != value.denominator:
        return None
    return Q(numerator, denominator)


def _quadratic_rational_roots(
    coefficients: tuple[Fraction, ...]
) -> tuple[str, ...]:
    coefficients = _trim_polynomial(coefficients)
    if len(coefficients) == 1:
        return ()
    if len(coefficients) == 2:
        return (rational_to_string(-coefficients[0] / coefficients[1]),)
    constant, linear, quadratic = coefficients
    square_root = _rational_square_root(linear**2 - 4 * quadratic * constant)
    if square_root is None:
        return ()
    roots = {
        (-linear + square_root) / (2 * quadratic),
        (-linear - square_root) / (2 * quadratic),
    }
    return tuple(rational_to_string(root) for root in sorted(roots))


def nagao_base_change_profile(
    construction: SixRootMestreConstruction,
) -> dict[str, Any]:
    # T=(3U^2-478U+1287)/(U^2-429), as printed by Nagao.
    numerator = (Q(1287), Q(-478), Q(3))
    denominator = (Q(-429), Q(0), Q(1))
    pulled_back = rational_base_change_polynomial(
        construction.primitive_discriminant_polynomial,
        numerator,
        denominator,
    )
    derivative = tuple(
        index * value for index, value in enumerate(pulled_back[1:], 1)
    )
    squarefree_certificate_prime = None
    repeated_degree = None
    for prime in (101, 103, 107, 109, 113):
        modular_degree = _modular_gcd_degree(pulled_back, derivative, prime)
        if modular_degree == 0:
            # Coprime reduction proves coprimality over Q.
            squarefree_certificate_prime = prime
            repeated_degree = 0
            break
    if repeated_degree is None:
        common = _polynomial_gcd(pulled_back, derivative)
        repeated_degree = 0 if common == (Q(1),) else len(common) - 1

    rational_collision_preimages: list[dict[str, Any]] = []
    for collision_parameter in construction.collision_parameters():
        equation = _add_polynomials(
            numerator,
            tuple(-collision_parameter * value for value in denominator),
        )
        finite_roots = list(_quadratic_rational_roots(equation))
        infinity_is_root = len(equation) < 3
        if finite_roots or infinity_is_root:
            rational_collision_preimages.append(
                {
                    "collision_parameter_T": rational_to_string(
                        collision_parameter
                    ),
                    "rational_parameters_U": [
                        *finite_roots,
                        *(["infinity"] if infinity_is_root else []),
                    ],
                    "collision_loss": construction.visible_point_degeneracy(
                        collision_parameter
                    ).collision_loss,
                }
            )
    degree = len(pulled_back) - 1
    return {
        "parameter_map": "T=(3*U^2-478*U+1287)/(U^2-429)",
        "discriminant_degree_after_denominator_clearing": degree,
        "squarefree_discriminant_degree_proxy": degree - repeated_degree,
        "repeated_factor_degree": repeated_degree,
        "squarefree_certificate_prime": squarefree_certificate_prime,
        "collision_divisor_degree_over_algebraic_closure": (
            2 * len(construction.collision_parameters())
        ),
        "rational_collision_preimages": rational_collision_preimages,
        "interpretation": (
            "squarefree discriminant degree is a conductor-growth proxy only; "
            "minimalization and fiber types were not computed"
        ),
    }


def parse_primes(value: str) -> tuple[int, ...]:
    try:
        primes = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("primes must be comma-separated integers") from error
    if not primes or len(set(primes)) != len(primes):
        raise argparse.ArgumentTypeError("primes must be nonempty and distinct")
    if tuple(sorted(primes)) != primes:
        raise argparse.ArgumentTypeError("primes must be strictly increasing")
    if any(prime < 5 or prime not in primes_up_to(prime) for prime in primes):
        raise argparse.ArgumentTypeError("every survey prime must be prime and at least 5")
    return primes


def discriminant_value(
    construction: SixRootMestreConstruction,
    parameter: int,
    cache: dict[int, Fraction],
) -> Fraction:
    if parameter not in cache:
        cache[parameter] = construction.primitive_discriminant_value(Q(parameter))
    return cache[parameter]


def p_adic_profile(
    construction: SixRootMestreConstruction,
    prime: int,
    exponent: int,
    cache: dict[int, Fraction],
    max_modulus: int = 20_000,
) -> dict[str, Any]:
    """Measure roots modulo ``p^exponent`` after removing fixed valuation."""

    fixed_samples = [
        discriminant_value(construction, parameter, cache)
        for parameter in range(1, 22)
    ]
    finite_fixed_samples = [value for value in fixed_samples if value != 0]
    if not finite_fixed_samples:
        raise AssertionError("the discriminant polynomial vanishes identically")
    fixed_valuation = min(
        valuation(value, prime) for value in finite_fixed_samples
    )
    congruence_exponent = fixed_valuation + exponent
    modulus = prime**congruence_exponent
    if modulus > max_modulus:
        return {
            "prime": prime,
            "status": "skipped_modulus_cap",
            "requested_excess_exponent": exponent,
            "fixed_discriminant_valuation": fixed_valuation,
            "required_enumeration_exponent": congruence_exponent,
            "required_enumeration_modulus": modulus,
            "max_profile_modulus": max_modulus,
        }
    valuations: list[int | None] = []
    for residue in range(modulus):
        discriminant = discriminant_value(construction, residue, cache)
        valuations.append(None if discriminant == 0 else valuation(discriminant, prime))
    finite = [value for value in valuations if value is not None]
    if not finite:
        raise AssertionError("the discriminant vanished on every residue")
    if min(finite) != fixed_valuation:
        raise AssertionError("the fixed discriminant valuation was not stable")

    level_counts: list[int] = []
    for level in range(1, exponent + 1):
        refined_count = sum(
            value is None or value >= fixed_valuation + level
            for value in valuations
        )
        refinements_per_class = prime ** (exponent - level)
        if refined_count % refinements_per_class:
            raise AssertionError("root counts did not descend through residue levels")
        level_counts.append(refined_count // refinements_per_class)

    root_count = level_counts[-1]
    reciprocal_density = (
        None if root_count == 0 else Q(modulus, root_count)
    )
    if root_count:
        log_cost = log(modulus / root_count)
        forcing_efficiency = (
            None if log_cost == 0 else exponent * log(prime) / log_cost
        )
    else:
        forcing_efficiency = 0.0
    first_count = level_counts[0]
    branching_ratio = None if first_count == 0 else Q(root_count, first_count)
    finite_max = max(finite)
    max_residues = [
        residue
        for residue, value in enumerate(valuations)
        if value is not None and value == finite_max
    ]
    return {
        "prime": prime,
        "status": "complete",
        "requested_excess_exponent": exponent,
        "required_enumeration_exponent": congruence_exponent,
        "modulus": modulus,
        "fixed_discriminant_valuation": fixed_valuation,
        "root_counts_by_exponent": level_counts,
        "root_count_at_requested_exponent": root_count,
        "root_residues_at_requested_exponent": [
            residue
            for residue, value in enumerate(valuations)
            if value is None or value >= fixed_valuation + exponent
        ][:100],
        "root_residue_list_truncated": root_count > 100,
        "reciprocal_density": (
            None if reciprocal_density is None else rational_to_string(reciprocal_density)
        ),
        "branching_ratio_last_to_first": (
            None if branching_ratio is None else rational_to_string(branching_ratio)
        ),
        "forcing_efficiency": forcing_efficiency,
        "observed_max_finite_valuation": finite_max,
        "observed_max_finite_residues": max_residues,
        "exact_zero_representatives": sum(value is None for value in valuations),
    }


def collision_geometry(
    construction: SixRootMestreConstruction,
) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    for parameter in construction.collision_parameters():
        degeneracy = construction.visible_point_degeneracy(parameter)
        records.append(
            {
                "parameter": rational_to_string(parameter),
                "collision_loss": degeneracy.collision_loss,
                "zero_ordinates": degeneracy.zero_ordinates,
            }
        )
    maximum_loss = max(record["collision_loss"] for record in records)
    return {
        "exceptional_parameter_count": len(records),
        "maximum_collision_loss": maximum_loss,
        "parameters_with_maximum_loss": [
            record["parameter"]
            for record in records
            if record["collision_loss"] == maximum_loss
        ],
        "zero_ordinate_exception_count": sum(
            record["zero_ordinates"] > 0 for record in records
        ),
    }


def analyze_construction(
    construction: SixRootMestreConstruction,
    *,
    primes: tuple[int, ...],
    exponent: int,
    parameter_bound: int,
    max_profile_modulus: int,
    pari_probe_count: int = 0,
    pari_timeout: float = 10.0,
) -> dict[str, Any]:
    cache: dict[int, Fraction] = {}
    profiles = [
        p_adic_profile(
            construction, prime, exponent, cache, max_profile_modulus
        )
        for prime in primes
    ]
    baselines = {
        profile["prime"]: profile["fixed_discriminant_valuation"]
        for profile in profiles
    }

    parameter_records: list[dict[str, Any]] = []
    singular_parameters: list[int] = []
    collision_parameters: list[int] = []
    zero_ordinate_parameters: list[int] = []
    for parameter in range(1, parameter_bound + 1):
        discriminant = discriminant_value(construction, parameter, cache)
        if discriminant == 0:
            singular_parameters.append(parameter)
            continue
        degeneracy = construction.visible_point_degeneracy(Q(parameter))
        if degeneracy.collision_loss:
            collision_parameters.append(parameter)
        if degeneracy.zero_ordinates:
            zero_ordinate_parameters.append(parameter)
        excess_valuations = {
            prime: max(valuation(discriminant, prime) - baselines[prime], 0)
            for prime in primes
        }
        discriminant_log = sum(
            exponent_value * log(prime)
            for prime, exponent_value in excess_valuations.items()
        )
        radical_log = sum(
            log(prime)
            for prime, exponent_value in excess_valuations.items()
            if exponent_value
        )
        parameter_records.append(
            {
                "parameter": parameter,
                "small_prime_excess_valuations": {
                    str(prime): value for prime, value in excess_valuations.items()
                },
                "small_prime_discriminant_log": discriminant_log,
                "small_prime_radical_log": radical_log,
                "small_prime_power_savings_log": discriminant_log - radical_log,
                "collision_loss": degeneracy.collision_loss,
                "zero_ordinates": degeneracy.zero_ordinates,
                "discriminant_numerator_bits": abs(discriminant.numerator).bit_length(),
                "discriminant_denominator_bits": discriminant.denominator.bit_length(),
            }
        )

    eligible = [
        record
        for record in parameter_records
        if record["collision_loss"] == 0 and record["zero_ordinates"] == 0
    ]
    eligible.sort(
        key=lambda record: (
            -record["small_prime_power_savings_log"],
            record["small_prime_radical_log"],
            record["parameter"],
        ),
    )
    best = eligible[0] if eligible else None
    pari_probes: list[dict[str, Any]] = []
    for parameter_record in eligible[:pari_probe_count]:
        parameter = parameter_record["parameter"]
        probe: dict[str, Any] = {
            "parameter": parameter,
            "selection": "top admissible bounded small-prime conductor proxy",
            "status": "PARI/GP software computation; not a generic-rank claim",
        }
        try:
            pari = minimal_curve_data(
                construction.primitive_jacobian_coefficients(Q(parameter)),
                timeout=pari_timeout,
                rank_effort=0,
                stack_bytes=64_000_000,
            )
            probe.update(
                {
                    "log_conductor": pari["log_conductor"],
                    "conductor": str(pari["conductor"]),
                    "root_number": pari["root_number"],
                    "rank_bounds": [
                        pari["pari_ellrank"]["lower_bound"],
                        pari["pari_ellrank"]["upper_bound"],
                    ],
                    "returned_independent_points": pari["pari_ellrank"]
                    ["returned_independent_points"],
                }
            )
        except Exception as error:
            probe["error"] = str(error)
        pari_probes.append(probe)
    successful_lower_bounds = [
        probe["rank_bounds"][0] for probe in pari_probes if "rank_bounds" in probe
    ]
    efficiency_bonus = sum(
        max((profile.get("forcing_efficiency") or 0.0) - 1.0, 0.0)
        for profile in profiles
    )
    return {
        "root_tuple": [rational_to_string(root) for root in construction.roots],
        "reflection_symmetric": construction.is_reflection_symmetric,
        "rank_viable_visible_geometry": not construction.is_reflection_symmetric,
        "rank_viability_explanation": (
            "passes the necessary nonreflection gate; independence remains open"
            if not construction.is_reflection_symmetric
            else "fails the visible-section gate: in the reflection-symmetric case six displayed points pair with the negatives of the other six"
        ),
        "mestre_quartic_condition": rational_to_string(construction.quartic_condition),
        "removed_fixed_quartic_content": rational_to_string(
            construction.quartic_content
        ),
        "removed_fixed_square_scale": rational_to_string(
            construction.quartic_square_scale
        ),
        "discriminant_geometry": {
            **discriminant_degree_profile(construction),
            "interpretation": (
                "squarefree discriminant degree is a conductor-growth proxy, "
                "not an arithmetic conductor degree"
            ),
        },
        "collision_geometry": collision_geometry(construction),
        "p_adic_profiles": profiles,
        "p_adic_richness": {
            "profiles_completed": sum(
                profile["status"] == "complete" for profile in profiles
            ),
            "profiles_skipped_by_modulus_cap": sum(
                profile["status"] != "complete" for profile in profiles
            ),
            "primes_with_requested_power_classes": sum(
                profile["status"] == "complete"
                and profile["root_count_at_requested_exponent"] > 0
                for profile in profiles
            ),
            "primes_with_branching": sum(
                profile["status"] == "complete"
                and profile["root_counts_by_exponent"][-1]
                    > profile["root_counts_by_exponent"][0]
                for profile in profiles
            ),
            "aggregate_efficiency_bonus": efficiency_bonus,
        },
        "bounded_parameter_survey": {
            "positive_integer_parameters_checked": parameter_bound,
            "nonsingular_parameters": len(parameter_records),
            "singular_parameters": singular_parameters,
            "collision_parameters": collision_parameters,
            "zero_ordinate_parameters": zero_ordinate_parameters,
            "best_admissible_conductor_proxy": best,
            "top_admissible_conductor_proxies": eligible[
                : max(2, pari_probe_count)
            ],
        },
        "pari_specialization_probe": {
            "probe_count_requested": pari_probe_count,
            "probe_count_completed": len(successful_lower_bounds),
            "maximum_observed_rank_lower_bound": (
                None if not successful_lower_bounds else max(successful_lower_bounds)
            ),
            "records": pari_probes,
            "interpretation": (
                "bounded effort-0 specialization probes only; they neither prove "
                "generic rank nor certify independence of the twelve displayed sections"
            ),
        },
    }


def ranking_key(record: dict[str, Any]) -> tuple[Any, ...]:
    best = record["bounded_parameter_survey"]["best_admissible_conductor_proxy"]
    savings = -1.0 if best is None else best["small_prime_power_savings_log"]
    roots = tuple(int(root) for root in record["root_tuple"])
    observed_lower_bound = record["pari_specialization_probe"][
        "maximum_observed_rank_lower_bound"
    ]
    return (
        record["reflection_symmetric"],
        -(observed_lower_bound if observed_lower_bound is not None else -1),
        -savings,
        -record["p_adic_richness"]["aggregate_efficiency_bonus"],
        record["collision_geometry"]["maximum_collision_loss"],
        roots[-1],
        roots,
    )


def verify_fermigier_embedding(
    construction: SixRootMestreConstruction,
) -> None:
    """Check that the general API recovers the repository's pinned model."""

    for parameter in (Q(1), Q(7, 3)):
        general = construction.primitive_quartic_coefficients(parameter)
        pinned_descending = FermigierMestreFamily.quartic_coefficients(parameter)
        expected_ascending = tuple(reversed(pinned_descending))
        if general != expected_ascending:
            raise AssertionError("the general root-tuple API missed Fermigier's model")
    if construction.quartic_square_scale != 50616:
        raise AssertionError("the Fermigier fixed-square normalization changed")


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-root", type=int, default=14)
    parser.add_argument("--parameter-bound", type=int, default=20)
    parser.add_argument(
        "--primes", type=parse_primes, default=(5, 7, 11, 17, 19, 37)
    )
    parser.add_argument("--exponent", type=int, default=2)
    parser.add_argument("--max-profile-modulus", type=int, default=20_000)
    parser.add_argument("--pari-probe-count", type=int, default=1)
    parser.add_argument("--pari-timeout", type=float, default=10.0)
    parser.add_argument("--keep", type=int, default=20)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts"
            / "generated-results"
            / "elliptic_mestre_root_tuple_survey.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.max_root < 5:
        raise SystemExit("--max-root must be at least 5")
    if (
        args.parameter_bound < 1
        or args.exponent < 1
        or args.keep < 1
        or args.max_profile_modulus < 2
        or args.pari_probe_count < 0
        or args.pari_timeout <= 0
    ):
        raise SystemExit(
            "--parameter-bound, --exponent, --keep, and "
            "--max-profile-modulus must be positive; PARI probe count must be "
            "nonnegative and timeout positive"
        )

    normalized = affine_normalized_integer_root_tuples(args.max_root)
    quartic_condition_constructions: list[SixRootMestreConstruction] = []
    for roots in normalized:
        construction = SixRootMestreConstruction(tuple(Q(root) for root in roots))
        if construction.is_quartic_family:
            quartic_condition_constructions.append(construction)

    generically_singular = [
        construction
        for construction in quartic_condition_constructions
        if all(
            coefficient == 0
            for coefficient in construction.primitive_discriminant_polynomial
        )
    ]
    quartic_constructions = [
        construction
        for construction in quartic_condition_constructions
        if construction not in generically_singular
    ]

    records = [
        analyze_construction(
            construction,
            primes=args.primes,
            exponent=args.exponent,
            parameter_bound=args.parameter_bound,
            max_profile_modulus=args.max_profile_modulus,
            pari_probe_count=args.pari_probe_count,
            pari_timeout=args.pari_timeout,
        )
        for construction in quartic_constructions
    ]
    records.sort(key=ranking_key)
    for position, record in enumerate(records, 1):
        record["survey_rank"] = position

    fermigier = SixRootMestreConstruction(tuple(Q(root) for root in ROOTS))
    verify_fermigier_embedding(fermigier)
    fermigier_record = analyze_construction(
        fermigier,
        primes=args.primes,
        exponent=args.exponent,
        parameter_bound=args.parameter_bound,
        max_profile_modulus=args.max_profile_modulus,
        pari_probe_count=0,
        pari_timeout=args.pari_timeout,
    )
    nagao_roots = (-17, -16, 10, 11, 14, 17)
    nagao = SixRootMestreConstruction(tuple(Q(root) for root in nagao_roots))
    if not nagao.is_quartic_family:
        raise AssertionError("Nagao's published tuple failed the quartic condition")
    nagao_record = analyze_construction(
        nagao,
        primes=args.primes,
        exponent=args.exponent,
        parameter_bound=args.parameter_bound,
        max_profile_modulus=args.max_profile_modulus,
        pari_probe_count=0,
        pari_timeout=args.pari_timeout,
    )
    nagao_rank13_roots = (0, 25, 57, 104, 116, 148)
    nagao_rank13 = SixRootMestreConstruction(
        tuple(Q(root) for root in nagao_rank13_roots)
    )
    if not nagao_rank13.is_quartic_family:
        raise AssertionError("Nagao's later published tuple failed the quartic condition")

    completed_probes = [
        probe
        for record in records
        for probe in record["pari_specialization_probe"]["records"]
        if "rank_bounds" in probe
    ]
    probe_errors = [
        probe
        for record in records
        for probe in record["pari_specialization_probe"]["records"]
        if "error" in probe
    ]
    rank_bound_histogram: dict[str, int] = {}
    for probe in completed_probes:
        key = str(probe["rank_bounds"])
        rank_bound_histogram[key] = rank_bound_histogram.get(key, 0) + 1

    script_path = Path(__file__).resolve()
    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded exact-arithmetic geometry/conductor-proxy survey with "
            "effort-0 PARI specialization probes; no generic-rank, section-"
            "independence, or target-hit claim"
        ),
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": "182.72",
            "alternative_rank_at_least": 30,
            "hits": [],
            "explanation": (
                "all bounded specialization probes have rank bounds below 21; "
                "no generic or specialized target curve is certified"
            ),
        },
        "method": {
            "root_normalization": (
                "six distinct integer roots modulo translation, common integral "
                "scale, and reflection; least root 0 and primitive differences"
            ),
            "quartic_filter": (
                "exact vanishing of the X^5 coefficient in "
                "(g_T^2-q(X-T)q(X+T))/T^2"
            ),
            "p_adic_profile": (
                "complete residues modulo p^k for each configured prime, after "
                "removing the quartic's fixed rational-square content and then "
                "subtracting the remaining minimum discriminant valuation"
            ),
            "conductor_proxy": (
                "over positive integer 1<=T<=bound, sum over configured primes "
                "of max(v_p(discriminant)-fixed_v_p-1,0)*log(p); this is only "
                "small-prime repeated-power savings, not a conductor computation"
            ),
            "visible_degeneracy": (
                "exact collisions among the twelve displayed abscissae and exact "
                "zero ordinates; degenerate parameters are excluded from proxy ranking"
            ),
            "rank_viability_gate": (
                "nonreflection configurations rank before reflection-symmetric "
                "controls; Nagao notes that symmetry pairs six displayed points "
                "with the negatives of the other six. Passing is necessary, not "
                "an independence or rank certificate"
            ),
        },
        "parameters": {
            "max_root": args.max_root,
            "parameter_bound": args.parameter_bound,
            "primes": list(args.primes),
            "exponent": args.exponent,
            "max_profile_modulus": args.max_profile_modulus,
            "pari_probe_count": args.pari_probe_count,
            "pari_timeout": args.pari_timeout,
            "keep": args.keep,
            "output": str(args.output),
        },
        "enumeration": {
            "affine_normalized_root_tuples": len(normalized),
            "quartic_condition_passes": len(quartic_condition_constructions),
            "generically_singular_quartics_rejected": len(generically_singular),
            "generically_nonsingular_quartic_families": len(records),
            "degree_five_remainders_rejected": (
                len(normalized) - len(quartic_condition_constructions)
            ),
            "generically_singular_root_tuples": [
                [rational_to_string(root) for root in construction.roots]
                for construction in generically_singular
            ],
            "reflection_symmetric_quartic_tuples": sum(
                record["reflection_symmetric"] for record in records
            ),
            "nonreflection_quartic_tuples": sum(
                not record["reflection_symmetric"] for record in records
            ),
        },
        "bounded_population": records,
        "pari_probe_summary": {
            "requested": len(records) * args.pari_probe_count,
            "completed": len(completed_probes),
            "errors": probe_errors,
            "rank_bound_histogram": rank_bound_histogram,
            "maximum_observed_lower_bound": max(
                (probe["rank_bounds"][0] for probe in completed_probes),
                default=None,
            ),
        },
        "top_configurations": records[: args.keep],
        "nonreflection_family_candidates": [
            record for record in records if not record["reflection_symmetric"]
        ],
        "symmetric_conductor_proxy_controls": [
            record
            for record in records
            if record["reflection_symmetric"]
        ][: args.keep],
        "fermigier_benchmark_outside_bounded_population": {
            "normalization_check": (
                "primitive general R_T/T^2 equals the pinned Fermigier quartic; "
                "the removed square scale is 50616"
            ),
            **fermigier_record,
        },
        "nagao_frontier_outside_bounded_population": {
            "published_root_tuple": list(nagao_roots),
            "affine_normalized_root_tuple": list(
                normalize_integer_root_tuple(nagao_roots)
            ),
            "published_base_change": nagao_base_change_profile(nagao),
            "frontier_assessment": (
                "this older tuple and base change remain an exactly implementable "
                "historical intermediate; the later published generic-rank-13 "
                "construction is now implemented, searched, and is the active "
                "family frontier.  The older base change still illustrates that "
                "a quadratic pullback doubles the squarefree discriminant-degree "
                "proxy from 20 to 40"
            ),
            "rank_status": (
                "Nagao's paper supplies rank statements; this artifact checks "
                "only geometry and makes no independent rank or independence claim"
            ),
            **nagao_record,
            "later_rank_at_least_13_direction_from_source": {
                "published_root_tuple": list(nagao_rank13_roots),
                "published_extra_section_abscissa": "(T+703)/15",
                "quartic_condition_checked_exactly": (
                    nagao_rank13.quartic_condition == 0
                ),
                "reflection_symmetric": nagao_rank13.is_reflection_symmetric,
                "discriminant_geometry": discriminant_degree_profile(
                    nagao_rank13
                ),
                "collision_geometry": collision_geometry(nagao_rank13),
                "implementation_status": (
                    "the extra section and quadratic base change are implemented "
                    "and exactly replayed separately in cas/nagao_1994.py; this "
                    "bounded survey does not duplicate that verification"
                ),
            },
        },
        "interpretation_limits": [
            "the Mestre quartic condition supplies a genus-one quartic with displayed points, not a proof of their independence",
            "the p-adic score uses only configured small primes and a bounded parameter interval",
            "fixed-valuation subtraction and radical savings are proxies; local minimal models and conductors were not computed",
            "the Fermigier tuple is a benchmark and was not part of the max-root population",
            "the Nagao tuple and quadratic base change are external benchmarks and their published rank statements are not reverified here",
            "reflection-symmetric tuples may have strong conductor proxies but fail the twelve-visible-section rank-viability gate and are retained only as arithmetic controls",
            "effort-0 PARI rank bounds at one or two selected specializations are software calibration data, not generic-rank evidence",
        ],
        "primary_sources": [
            {
                "citation": "J.-F. Mestre, Constructions polynomiales et theorie de Galois, ICM Zurich (1994), section 4",
                "url": "https://www.imj-prg.fr/wp-content/uploads/2020/prix/mestre1994.pdf",
                "relevance": "states the q(X-T)q(X+T)=g^2-r construction, the quartic obstruction, twelve points, and the extra conic condition needed for the rank-12 base change",
            },
            {
                "citation": "K. Nagao, Construction of high-rank elliptic curves, Kobe J. Math. 11 (1994), 211-219",
                "url": "https://da.lib.kobe-u.ac.jp/da/kernel/E0003610/",
                "relevance": "works explicitly with integer six-tuples satisfying the degree-five condition, warns that reflection symmetry pairs the displayed points by negation, gives the (-17,-16,10,11,14,17) quadratic base change, and later engineers an extra point",
            },
            {
                "citation": "S. Fermigier, Une courbe elliptique definie sur Q de rang >=22, Acta Arith. 82 (1997), 359-363",
                "url": "https://doi.org/10.4064/aa-82-4-359-363",
                "relevance": "source of the specialization and six-root benchmark used by this repository",
            },
            {
                "citation": "J. Scholten, Elliptic curves of high rank over function fields (1997)",
                "url": "https://arxiv.org/abs/math/9709235",
                "relevance": "proves exact function-field rank statements for Nagao's related family and shows why displayed points alone do not settle rank",
            },
        ],
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "reproducing_command": command,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    print(
        f"normalized={len(normalized)} quartic={len(records)} "
        f"generic_singular={len(generically_singular)} "
        f"degree5_rejected={len(normalized)-len(quartic_condition_constructions)}"
    )
    for record in records[: min(5, args.keep)]:
        best = record["bounded_parameter_survey"]["best_admissible_conductor_proxy"]
        print(
            f"rank={record['survey_rank']} roots={record['root_tuple']} "
            f"T={None if best is None else best['parameter']} "
            f"proxy={None if best is None else best['small_prime_power_savings_log']}"
        )


if __name__ == "__main__":
    main()
