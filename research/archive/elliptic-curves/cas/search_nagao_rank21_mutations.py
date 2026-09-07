#!/usr/bin/env python3
"""Search conductor-feasible rational mutations of the Nagao rank-21 family.

The starting curve is the specialization ``T=5777/32``.  Its actual local
residues are ``(1,4,9,3)`` modulo ``(5,7,11,13)``.  A second, separately
labelled branch uses ``(1,4,9,10)``.  Both branches force unusually large
primitive-discriminant valuations and split multiplicative reduction at all
four primes.

The search is deliberately leakage-free:

* rational mutations are generated exactly in declared CRT lattice slices;
* an inexpensive radical proxy removes candidates that cannot plausibly meet
  the conductor target;
* the first quartic box ranks candidates only by exact, decontaminated point
  yield (visible sections and repeated Jacobian sign-pairs are removed);
* larger boxes and numerical height matrices are used only on retained
  candidates; and
* exact conductor and finite-reduction independence replays are attempted on
  the deepest leaders.

Numerical height rank remains triage evidence.  A finite-reduction certificate
is the only route by which this script records an algebraic rank lower bound.
Every PARI invocation is synchronous and has a hard subprocess timeout.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import gcd, log
from pathlib import Path
import platform
import shlex
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

from crt_lattice import crt_pair, gauss_reduce
from ek_k3 import legendre_symbol, primes_up_to, rational_to_string
from extend_nagao_u42_frontier import saturate_exact_basis
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from multiple_root_lifting import (
    affine_variable_coefficients,
    fixed_divisor_valuation,
)
from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK21_CONSTRUCTION,
    primitive_visible_points,
    short_jacobian_coefficients,
)
from pari_bridge import minimal_curve_data, pari_version
from search_nagao_rank21_neighborhood import (
    DISCRIMINANT_POLYNOMIAL,
    LocalCondition,
    classify_condition,
    discriminant_valuation_at_rational,
)
from triage_nagao_rank13_finalists import (
    height_matrix_replay,
    point_digest,
    stable_height_rank,
)
from triage_nagao_rank21_neighbor import (
    bounded_quartic_points,
    exact_visible_seeds,
    map_and_deduplicate,
)


Q = Fraction
TARGET_LOG_CONDUCTOR = Decimal("182.72")
LEAD_PARAMETER = Q(5777, 32)
CRT_PRIMES = (5, 7, 11, 13)
CRT_MODULUS = 5 * 7 * 11 * 13
TRIAL_PRIME_BOUND = 1_000
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank21_mutations.py"
)


@dataclass(frozen=True)
class MutationProfile:
    label: str
    residues: tuple[int, int, int, int]

    @property
    def crt_residue(self) -> int:
        residue, modulus = 0, 1
        for prime, local_residue in zip(CRT_PRIMES, self.residues):
            residue, modulus = crt_pair(
                residue, modulus, local_residue, prime
            )
        if modulus != CRT_MODULUS:
            raise AssertionError("the profile CRT modulus changed")
        return residue


LEAD_PROFILE = MutationProfile("lead-exact-p13-3", (1, 4, 9, 3))
REQUESTED_PROFILE = MutationProfile("alternate-p13-10", (1, 4, 9, 10))
DEFAULT_PROFILES = (LEAD_PROFILE, REQUESTED_PROFILE)


@dataclass(frozen=True)
class MutationRegion:
    label: str
    minimum_denominator: int
    maximum_denominator: int
    offset_radius: int
    rational_radius: Fraction
    exclude_numerator_height_at_most: int | None = None


# The first slice is outside the earlier |a|<=10,000 scan.  The second extends
# the denominator from 200 through 1,000 while staying near the lead in R.
DEFAULT_REGIONS = (
    MutationRegion("new-large-numerator", 1, 200, 30, Q(1_000), 10_000),
    MutationRegion("new-denominator-extension", 201, 1_000, 3, Q(100), None),
)


@dataclass(frozen=True)
class MutationCandidate:
    parameter: Fraction
    profile: MutationProfile
    region_label: str
    radical_proxy: dict[str, Any]
    rank_friendly_residues: tuple[tuple[int, int], ...] = ()

    @property
    def identifier(self) -> str:
        sign = "m" if self.parameter.numerator < 0 else "p"
        return (
            f"{self.profile.label}-{self.region_label}-{sign}"
            f"{abs(self.parameter.numerator)}-{self.parameter.denominator}"
        )


@dataclass(frozen=True)
class PointSearchResult:
    candidate: MutationCandidate
    height_bound: int
    signed_point_count: int
    signless_point_count: int
    visible_abscissas_returned: int
    new_images: tuple[tuple[Fraction, Fraction], ...]
    zero_ordinate_count: int
    wall_seconds: float
    pari_milliseconds: int
    status: str
    error: str | None = None


def rational_residue(parameter: Fraction, prime: int) -> int:
    parameter = Q(parameter)
    if parameter.denominator % prime == 0:
        raise ValueError("the denominator is not a local unit")
    return parameter.numerator * pow(parameter.denominator, -1, prime) % prime


def profile_contains(profile: MutationProfile, parameter: Fraction) -> bool:
    return all(
        rational_residue(parameter, prime) == residue
        for prime, residue in zip(CRT_PRIMES, profile.residues)
    )


def homogenized_discriminant(parameter: Fraction) -> int:
    parameter = Q(parameter)
    numerator, denominator = parameter.numerator, parameter.denominator
    degree = len(DISCRIMINANT_POLYNOMIAL) - 1
    value = sum(
        coefficient * numerator**power * denominator ** (degree - power)
        for power, coefficient in enumerate(DISCRIMINANT_POLYNOMIAL)
    )
    if value == 0:
        raise ValueError("the mutation gives a singular specialization")
    return value


def conductor_radical_proxy(
    parameter: Fraction, *, trial_prime_bound: int = TRIAL_PRIME_BOUND
) -> dict[str, Any]:
    """Return a one-sided small-prime approximation to ``log(rad(Delta))``.

    Repeated factors through the declared bound are removed.  An un-factored
    cofactor is treated as squarefree, so the value is an upper proxy for the
    logarithm of the radical of this particular discriminant polynomial.  It
    is not an upper bound for the conductor because additive local exponents
    can exceed one.
    """

    discriminant = abs(homogenized_discriminant(parameter))
    remaining = discriminant
    repeated_savings = 0.0
    valuations: list[list[int]] = []
    for prime in primes_up_to(trial_prime_bound):
        valuation = 0
        while remaining % prime == 0:
            valuation += 1
            remaining //= prime
        if valuation:
            valuations.append([prime, valuation])
            repeated_savings += (valuation - 1) * log(prime)
    raw_log = log(discriminant)
    return {
        "trial_prime_bound": trial_prime_bound,
        "raw_log_absolute_homogenized_discriminant": raw_log,
        "known_repeated_prime_log_savings": repeated_savings,
        "log_radical_upper_proxy": raw_log - repeated_savings,
        "small_prime_valuations": valuations,
        "unfactored_cofactor_decimal_digits": len(str(remaining)),
    }


def nearest_congruent_numerator(
    target: Fraction, denominator: int, crt_residue: int
) -> tuple[int, int]:
    """Return ``a0,k`` with ``a0+k*M`` nearest to ``target*denominator``."""

    if denominator <= 0 or gcd(denominator, CRT_MODULUS) != 1:
        raise ValueError("the denominator must be a positive CRT unit")
    base = crt_residue * denominator % CRT_MODULUS
    displacement = target * denominator - base
    # Exact nearest-integer rounding, with deterministic ties toward +infinity.
    quotient, remainder = divmod(displacement.numerator, displacement.denominator * CRT_MODULUS)
    if 2 * remainder >= displacement.denominator * CRT_MODULUS:
        quotient += 1
    return base, quotient


def enumerate_mutations(
    *,
    profiles: Sequence[MutationProfile] = DEFAULT_PROFILES,
    regions: Sequence[MutationRegion] = DEFAULT_REGIONS,
    proxy_limit: Decimal = TARGET_LOG_CONDUCTOR,
) -> tuple[MutationCandidate, ...]:
    """Enumerate and exactly deduplicate the declared CRT lattice slices."""

    retained: dict[Fraction, MutationCandidate] = {}
    for profile in profiles:
        crt_residue = profile.crt_residue
        for region in regions:
            if (
                region.minimum_denominator < 1
                or region.maximum_denominator < region.minimum_denominator
                or region.offset_radius < 0
                or region.rational_radius < 0
            ):
                raise ValueError("invalid mutation region")
            for denominator in range(
                region.minimum_denominator, region.maximum_denominator + 1
            ):
                if gcd(denominator, CRT_MODULUS) != 1:
                    continue
                base, center = nearest_congruent_numerator(
                    LEAD_PARAMETER, denominator, crt_residue
                )
                for offset in range(-region.offset_radius, region.offset_radius + 1):
                    parameter = Q(
                        base + (center + offset) * CRT_MODULUS, denominator
                    )
                    if parameter == 0 or abs(parameter - LEAD_PARAMETER) > region.rational_radius:
                        continue
                    if not (
                        region.minimum_denominator
                        <= parameter.denominator
                        <= region.maximum_denominator
                    ):
                        continue
                    excluded_height = region.exclude_numerator_height_at_most
                    if (
                        excluded_height is not None
                        and abs(parameter.numerator) <= excluded_height
                    ):
                        continue
                    if not profile_contains(profile, parameter):
                        raise AssertionError("primitive reduction lost the CRT profile")
                    proxy = conductor_radical_proxy(parameter)
                    if Decimal(str(proxy["log_radical_upper_proxy"])) >= proxy_limit:
                        continue
                    candidate = MutationCandidate(
                        parameter=parameter,
                        profile=profile,
                        region_label=region.label,
                        radical_proxy=proxy,
                    )
                    previous = retained.get(parameter)
                    if previous is None or (
                        candidate.region_label,
                        candidate.profile.label,
                    ) < (previous.region_label, previous.profile.label):
                        retained[parameter] = candidate
    return tuple(
        sorted(
            retained.values(),
            key=lambda candidate: (
                candidate.radical_proxy["log_radical_upper_proxy"],
                max(abs(candidate.parameter.numerator), candidate.parameter.denominator),
                candidate.identifier,
            ),
        )
    )


def best_good_prime_residues(prime: int) -> tuple[dict[str, Any], ...]:
    """Return every nonsingular residue attaining the smallest trace ``a_p``.

    At fixed ``p`` the Nagao summand used in this repository is monotone in
    ``-a_p``.  Thus these are exactly the most rank-friendly local symbols at
    the prime.  This computes the lookup table from the family itself rather
    than pinning residue choices by hand.
    """

    if prime <= 3 or prime not in primes_up_to(prime):
        raise ValueError("the lookup modulus must be a prime greater than three")
    records: list[dict[str, Any]] = []
    for residue in range(prime):
        # T=0 is removable in the polynomial family, but the direct constructor
        # has T^2 in a denominator.  T=p is the same local specialization.
        representative = residue or prime
        coefficients = short_jacobian_coefficients(
            RANK21_CONSTRUCTION, Q(representative)
        )
        coefficient_a = (
            coefficients[3].numerator
            * pow(coefficients[3].denominator, -1, prime)
            % prime
        )
        coefficient_b = (
            coefficients[4].numerator
            * pow(coefficients[4].denominator, -1, prime)
            % prime
        )
        if (4 * coefficient_a**3 + 27 * coefficient_b**2) % prime == 0:
            continue
        trace = -sum(
            legendre_symbol(
                (x_value**3 + coefficient_a * x_value + coefficient_b) % prime,
                prime,
            )
            for x_value in range(prime)
        )
        records.append(
            {
                "prime": prime,
                "residue": residue,
                "ellap": trace,
                "nagao_summand": (
                    (2 - trace) / (prime + 1 - trace) * log(prime)
                ),
            }
        )
    best_trace = min(record["ellap"] for record in records)
    return tuple(
        record for record in records if record["ellap"] == best_trace
    )


def enumerate_rank_friendly_lattice_mutations(
    *,
    profile: MutationProfile = LEAD_PROFILE,
    good_primes: tuple[int, int] = (31, 43),
    coefficient_radius: int = 40,
    proxy_limit: Decimal = TARGET_LOG_CONDUCTOR,
) -> tuple[MutationCandidate, ...]:
    """Combine best local symbols and enumerate reduced CRT lattice boxes."""

    if coefficient_radius < 1 or len(set(good_primes)) != len(good_primes):
        raise ValueError("invalid rank-friendly lattice bounds")
    lookup = tuple(best_good_prime_residues(prime) for prime in good_primes)
    retained: dict[Fraction, MutationCandidate] = {}
    for local_symbols in product(*lookup):
        residue, modulus = profile.crt_residue, CRT_MODULUS
        symbol_pairs = tuple(
            (int(symbol["prime"]), int(symbol["residue"]))
            for symbol in local_symbols
        )
        for prime, local_residue in symbol_pairs:
            residue, modulus = crt_pair(
                residue, modulus, local_residue, prime
            )
        basis = gauss_reduce((modulus, 0), (residue, 1))
        for left in range(-coefficient_radius, coefficient_radius + 1):
            for right in range(-coefficient_radius, coefficient_radius + 1):
                if left == 0 and right == 0:
                    continue
                numerator = left * basis[0][0] + right * basis[1][0]
                denominator = left * basis[0][1] + right * basis[1][1]
                if denominator == 0:
                    continue
                parameter = Q(numerator, denominator)
                if parameter == 0 or gcd(parameter.denominator, modulus) != 1:
                    continue
                if not profile_contains(profile, parameter):
                    # A nonprimitive lattice vector may have common content
                    # divisible by the CRT modulus.  Dividing out that content
                    # need not preserve its residue; it belongs to another
                    # projective class and is outside this declared slice.
                    continue
                if any(
                    rational_residue(parameter, prime) != local_residue
                    for prime, local_residue in symbol_pairs
                ):
                    continue
                try:
                    proxy = conductor_radical_proxy(parameter)
                except ValueError:
                    continue
                if Decimal(str(proxy["log_radical_upper_proxy"])) >= proxy_limit:
                    continue
                label = "rank-friendly-" + "-".join(
                    f"p{prime}r{local_residue}"
                    for prime, local_residue in symbol_pairs
                )
                candidate = MutationCandidate(
                    parameter=parameter,
                    profile=profile,
                    region_label=label,
                    radical_proxy=proxy,
                    rank_friendly_residues=symbol_pairs,
                )
                previous = retained.get(parameter)
                if previous is None or candidate.identifier < previous.identifier:
                    retained[parameter] = candidate
    return tuple(
        sorted(
            retained.values(),
            key=lambda candidate: (
                candidate.radical_proxy["log_radical_upper_proxy"],
                max(abs(candidate.parameter.numerator), candidate.parameter.denominator),
                candidate.identifier,
            ),
        )
    )


def search_points(
    candidate: MutationCandidate,
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> PointSearchResult:
    started = time.monotonic()
    try:
        seed_quartic, seed_jacobian, _ = exact_visible_seeds(candidate.parameter)
        raw_points, wall_seconds, milliseconds = bounded_quartic_points(
            candidate.parameter,
            height_bound=height_bound,
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        mapped, new_images, zero_ordinates = map_and_deduplicate(
            candidate.parameter, raw_points, seed_quartic, seed_jacobian
        )
        return PointSearchResult(
            candidate=candidate,
            height_bound=height_bound,
            signed_point_count=len(raw_points),
            signless_point_count=len({point[0] for point in raw_points}),
            visible_abscissas_returned=sum(
                bool(record["visible_section_abscissa"]) for record in mapped
            ),
            new_images=new_images,
            zero_ordinate_count=len(zero_ordinates),
            wall_seconds=wall_seconds,
            pari_milliseconds=milliseconds,
            status="completed",
        )
    except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
        return PointSearchResult(
            candidate=candidate,
            height_bound=height_bound,
            signed_point_count=0,
            signless_point_count=0,
            visible_abscissas_returned=0,
            new_images=(),
            zero_ordinate_count=0,
            wall_seconds=time.monotonic() - started,
            pari_milliseconds=0,
            status="timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
            error=str(error)[:500],
        )


def result_priority(result: PointSearchResult) -> tuple[Any, ...]:
    return (
        result.status != "completed",
        -len(result.new_images),
        -result.signed_point_count,
        result.candidate.radical_proxy["log_radical_upper_proxy"],
        max(abs(result.candidate.parameter.numerator), result.candidate.parameter.denominator),
        result.candidate.identifier,
    )


def result_record(result: PointSearchResult, *, include_points: bool) -> dict[str, Any]:
    record: dict[str, Any] = {
        "candidate_id": result.candidate.identifier,
        "constructor_parameter": rational_to_string(result.candidate.parameter),
        "profile": result.candidate.profile.label,
        "profile_residues_mod_5_7_11_13": list(result.candidate.profile.residues),
        "profile_crt_residue_mod_5005": result.candidate.profile.crt_residue,
        "region": result.candidate.region_label,
        "rank_friendly_good_prime_residues": [
            list(pair) for pair in result.candidate.rank_friendly_residues
        ],
        "radical_proxy": result.candidate.radical_proxy,
        "height_bound": result.height_bound,
        "status": result.status,
        "signed_points_found": result.signed_point_count,
        "distinct_quartic_abscissas": result.signless_point_count,
        "visible_abscissas_returned": result.visible_abscissas_returned,
        "new_distinct_jacobian_sign_pairs": len(result.new_images),
        "zero_ordinate_points_not_mapped": result.zero_ordinate_count,
        "wall_seconds": result.wall_seconds,
        "pari_reported_milliseconds": result.pari_milliseconds,
    }
    if result.error is not None:
        record["error"] = result.error
    if include_points:
        record["new_exact_jacobian_images"] = [
            {
                "jacobian_x": rational_to_string(point[0]),
                "jacobian_y": rational_to_string(point[1]),
                "exact_membership_checked": True,
            }
            for point in result.new_images
        ]
        record["new_exact_jacobian_image_sha256"] = point_digest(result.new_images)
    return record


def numerical_rank_record(
    result: PointSearchResult,
    *,
    precisions: tuple[int, ...],
    timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    if result.status != "completed":
        return {"status": "not-run", "reason": "point search did not complete"}
    _, seed_jacobian, coefficients = exact_visible_seeds(result.candidate.parameter)
    pool = seed_jacobian + result.new_images
    try:
        runs = height_matrix_replay(
            coefficients,
            pool,
            precisions=precisions,
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        rank = stable_height_rank(runs)
        indices = tuple(runs[-1]["subset_indices_one_based"])
        selected = tuple(pool[index - 1] for index in indices)
        return {
            "status": "completed",
            "stable_numerical_rank": rank,
            "precision_runs": list(runs),
            "selected_subset_indices_one_based": list(indices),
            "selected_exact_points": selected,
            "selected_point_sha256": point_digest(selected),
        }
    except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
        return {
            "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
            "error": str(error)[:500],
        }


def rank_priority(
    result: PointSearchResult, rank_record: dict[str, Any] | None
) -> tuple[Any, ...]:
    stable_rank = (
        int(rank_record["stable_numerical_rank"])
        if rank_record is not None and rank_record.get("status") == "completed"
        else -1
    )
    return (-stable_rank, *result_priority(result))


def finite_reduction_certificate(
    result: PointSearchResult,
    rank_record: dict[str, Any],
    *,
    saturation_timeout: float,
    stack_bytes: int,
    certificate_prime_bound: int,
) -> dict[str, Any]:
    """Try to turn a stable numerical subset into an exact rank certificate."""

    points = tuple(rank_record["selected_exact_points"])
    coefficients = short_jacobian_coefficients(
        RANK21_CONSTRUCTION, result.candidate.parameter
    )
    signatures = find_mod2_reduction_certificate(
        coefficients, points, prime_bound=certificate_prime_bound
    )
    binary_rank = combined_mod2_rank(signatures, len(points))
    saturation: dict[str, Any] | None = None
    if binary_rank != len(points):
        saturated, saturation = saturate_exact_basis(
            coefficients,
            points,
            prime_bound=20,
            timeout=saturation_timeout,
            stack_bytes=stack_bytes,
        )
        points = saturated
        signatures = find_mod2_reduction_certificate(
            coefficients, points, prime_bound=certificate_prime_bound
        )
        binary_rank = combined_mod2_rank(signatures, len(points))
    certified = binary_rank == len(points)
    two_torsion_prime = (
        find_two_torsion_certificate_prime(coefficients, prime_bound=200)
        if certified
        else None
    )
    return {
        "status": "certified" if certified else "bounded-search-rank-deficient",
        "point_count": len(points),
        "point_sha256": point_digest(points),
        "small_prime_saturation": saturation,
        "certificate_prime_bound": certificate_prime_bound,
        "certificate_primes": [signature.prime for signature in signatures],
        "combined_exact_rank_over_F2": binary_rank,
        "two_torsion_certificate_prime": two_torsion_prime,
        "signatures": [
            {
                "prime": signature.prime,
                "group_order": signature.group_order,
                "doubled_subgroup_order": signature.doubled_subgroup_order,
                "quotient_dimension": signature.quotient_dimension,
                "rows": [list(row) for row in signature.rows],
            }
            for signature in signatures
        ],
        "certified_algebraic_rank_lower_bound": len(points) if certified else None,
    }


def local_profile_verification() -> list[dict[str, Any]]:
    records = []
    roots = {5: (1, 4), 7: (0, 3, 4), 11: (2, 9), 13: (3, 10)}
    expected = {5: 7, 7: 4, 11: 5, 13: 4}
    for prime, residues in roots.items():
        for residue in residues:
            forced = fixed_divisor_valuation(
                affine_variable_coefficients(
                    DISCRIMINANT_POLYNOMIAL, residue, prime
                ),
                prime,
            )
            if forced != expected[prime]:
                raise AssertionError("a pinned automatic local condition changed")
            condition = LocalCondition(
                f"auto-p{prime}-r{residue}", prime, 1, residue, forced, True
            )
            classification = classify_condition(condition)
            if (
                classification["reduction"] != "split multiplicative"
                or classification["conductor_exponent"] != 1
            ):
                raise AssertionError("a pinned split local classification changed")
            records.append(
                {
                    "prime": prime,
                    "residue": residue,
                    "forced_primitive_discriminant_valuation": forced,
                    "classification": classification,
                }
            )
    return records


def parse_positive_ints(value: str) -> tuple[int, ...]:
    try:
        values = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from error
    if not values or any(value < 1 for value in values):
        raise argparse.ArgumentTypeError("all values must be positive")
    return values


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-heights", type=parse_positive_ints, default=(5_000, 50_000, 250_000))
    parser.add_argument("--keep-counts", type=parse_positive_ints, default=(240, 24))
    parser.add_argument("--search-timeouts", type=parse_positive_ints, default=(3, 5, 20))
    parser.add_argument("--deep-height", type=int, default=1_000_000)
    parser.add_argument("--deep-keep", type=int, default=3)
    parser.add_argument("--deep-rank-threshold", type=int, default=15)
    parser.add_argument("--deep-timeout", type=float, default=75.0)
    parser.add_argument("--height-precisions", type=parse_positive_ints, default=(72, 120))
    parser.add_argument("--height-timeout", type=float, default=15.0)
    parser.add_argument("--conductor-timeout", type=float, default=20.0)
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=1_000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "artifacts/generated-results/elliptic_nagao_rank21_mutations.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if len(args.stage_heights) != 3 or len(args.keep_counts) != 2:
        raise SystemExit("this run requires three stage heights and two keep counts")
    if len(args.search_timeouts) != 3:
        raise SystemExit("provide one search timeout for each stage")
    if tuple(sorted(args.stage_heights)) != args.stage_heights:
        raise SystemExit("stage heights must be strictly increasing")
    if tuple(sorted(set(args.height_precisions))) != args.height_precisions:
        raise SystemExit("height precisions must be increasing and distinct")
    if min(
        *args.search_timeouts,
        args.deep_timeout,
        args.height_timeout,
        args.conductor_timeout,
        args.saturation_timeout,
    ) <= 0:
        raise SystemExit("all subprocess timeouts must be positive")
    if args.stack_bytes < 8_000_000 or args.certificate_prime_bound < 3:
        raise SystemExit("invalid PARI stack or certificate-prime bound")

    local_verification = local_profile_verification()
    if tuple(rational_residue(LEAD_PARAMETER, prime) for prime in CRT_PRIMES) != LEAD_PROFILE.residues:
        raise AssertionError("the lead parameter's exact local profile changed")
    if LEAD_PROFILE.crt_residue != 1901 or REQUESTED_PROFILE.crt_residue != 361:
        raise AssertionError("a pinned profile CRT residue changed")

    height_box_candidates = enumerate_mutations()
    constraint_candidates = enumerate_rank_friendly_lattice_mutations()
    combined_candidates = {
        candidate.parameter: candidate
        for candidate in (*height_box_candidates, *constraint_candidates)
    }
    # Prefer the explicitly rank-friendly provenance when the same rational
    # parameter appears in both independently generated populations.
    combined_candidates.update(
        {candidate.parameter: candidate for candidate in constraint_candidates}
    )
    candidates = tuple(
        sorted(
            combined_candidates.values(),
            key=lambda candidate: (
                candidate.radical_proxy["log_radical_upper_proxy"],
                max(abs(candidate.parameter.numerator), candidate.parameter.denominator),
                candidate.identifier,
            ),
        )
    )
    if len(candidates) < args.keep_counts[0]:
        raise SystemExit("the first keep count exceeds the mutation population")
    by_identifier = {candidate.identifier: candidate for candidate in candidates}

    stage_records: list[dict[str, Any]] = []
    survivors = list(candidates)
    latest_results: dict[str, PointSearchResult] = {}
    rank_records: dict[str, dict[str, Any]] = {}
    for stage_index, (height_bound, timeout) in enumerate(
        zip(args.stage_heights, args.search_timeouts), start=1
    ):
        results = [
            search_points(
                candidate,
                height_bound=height_bound,
                timeout=float(timeout),
                stack_bytes=args.stack_bytes,
            )
            for candidate in survivors
        ]
        results.sort(key=result_priority)
        for result in results:
            latest_results[result.candidate.identifier] = result

        # At the two expensive stages, replay heights for the exact leaders.
        stage_ranked_ids: list[str] = []
        if stage_index >= 2:
            rank_probe_count = min(32, len(results))
            for result in results[:rank_probe_count]:
                record = numerical_rank_record(
                    result,
                    precisions=args.height_precisions,
                    timeout=args.height_timeout,
                    stack_bytes=args.stack_bytes,
                )
                rank_records[result.candidate.identifier] = record
            results.sort(
                key=lambda result: rank_priority(
                    result, rank_records.get(result.candidate.identifier)
                )
            )
            stage_ranked_ids = [result.candidate.identifier for result in results[:rank_probe_count]]

        keep_count = (
            args.keep_counts[stage_index - 1]
            if stage_index <= len(args.keep_counts)
            else len(results)
        )
        retained = results[: min(keep_count, len(results))]
        stage_records.append(
            {
                "stage": stage_index,
                "quartic_naive_height_bound": height_bound,
                "search_timeout_seconds_per_candidate": timeout,
                "population_searched": len(results),
                "completed": sum(result.status == "completed" for result in results),
                "timeouts": sum(result.status == "timeout" for result in results),
                "errors": sum(result.status == "error" for result in results),
                "rank_replay_candidate_ids": stage_ranked_ids,
                "rank_replays": {
                    identifier: {
                        key: value
                        for key, value in rank_records[identifier].items()
                        if key != "selected_exact_points"
                    }
                    for identifier in stage_ranked_ids
                },
                "ranked_population": [
                    result_record(result, include_points=False) for result in results
                ],
                "retained_candidate_ids": [
                    result.candidate.identifier for result in retained
                ],
            }
        )
        survivors = [by_identifier[result.candidate.identifier] for result in retained]

    deepest_ranked = sorted(
        (latest_results[candidate.identifier] for candidate in survivors),
        key=lambda result: rank_priority(
            result, rank_records.get(result.candidate.identifier)
        ),
    )
    deep_eligible = [
        result
        for result in deepest_ranked
        if rank_records.get(result.candidate.identifier, {}).get("status") == "completed"
        and int(rank_records[result.candidate.identifier]["stable_numerical_rank"])
        >= args.deep_rank_threshold
    ][: args.deep_keep]
    # If no curve reaches the threshold, still deepen the best exact-point
    # leader once; this makes the negative bounded search informative.
    if not deep_eligible and deepest_ranked and args.deep_keep:
        deep_eligible = deepest_ranked[:1]

    deep_results: list[PointSearchResult] = []
    for prior in deep_eligible:
        result = search_points(
            prior.candidate,
            height_bound=args.deep_height,
            timeout=args.deep_timeout,
            stack_bytes=args.stack_bytes,
        )
        deep_results.append(result)
        rank_records[result.candidate.identifier] = numerical_rank_record(
            result,
            precisions=args.height_precisions,
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
        )

    final_pool = deep_results or deepest_ranked[: min(6, len(deepest_ranked))]
    exact_certificates: dict[str, dict[str, Any]] = {}
    conductor_replays: dict[str, dict[str, Any]] = {}
    certified_hits: list[dict[str, Any]] = []
    for result in final_pool:
        identifier = result.candidate.identifier
        rank_record = rank_records.get(identifier, {})
        if (
            rank_record.get("status") == "completed"
            and int(rank_record["stable_numerical_rank"]) >= 18
        ):
            try:
                exact_certificates[identifier] = finite_reduction_certificate(
                    result,
                    rank_record,
                    saturation_timeout=args.saturation_timeout,
                    stack_bytes=args.stack_bytes,
                    certificate_prime_bound=args.certificate_prime_bound,
                )
            except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
                exact_certificates[identifier] = {
                    "status": "error",
                    "error": str(error)[:500],
                }
        coefficients = short_jacobian_coefficients(
            RANK21_CONSTRUCTION, result.candidate.parameter
        )
        try:
            conductor = minimal_curve_data(
                coefficients,
                timeout=args.conductor_timeout,
                local_primes=CRT_PRIMES,
                stack_bytes=args.stack_bytes,
            )
            local_checks = {
                str(prime): discriminant_valuation_at_rational(
                    result.candidate.parameter.numerator,
                    result.candidate.parameter.denominator,
                    prime,
                )
                for prime in CRT_PRIMES
            }
            conductor_replays[identifier] = {
                "status": "completed",
                **conductor,
                "exact_primitive_discriminant_valuations": local_checks,
                "below_strict_log_conductor_target": (
                    Decimal(conductor["log_conductor"]) < TARGET_LOG_CONDUCTOR
                ),
            }
        except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
            conductor_replays[identifier] = {
                "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
                "error": str(error)[:500],
            }

        certificate = exact_certificates.get(identifier)
        conductor = conductor_replays.get(identifier)
        if (
            certificate is not None
            and certificate.get("status") == "certified"
            and int(certificate["certified_algebraic_rank_lower_bound"]) >= 21
            and conductor is not None
            and conductor.get("below_strict_log_conductor_target") is True
        ):
            certified_hits.append(
                {
                    "candidate_id": identifier,
                    "constructor_parameter": rational_to_string(result.candidate.parameter),
                    "certified_rank_lower_bound": certificate[
                        "certified_algebraic_rank_lower_bound"
                    ],
                    "conductor": conductor["conductor"],
                    "log_conductor": conductor["log_conductor"],
                }
            )

    lead_proxy = conductor_radical_proxy(LEAD_PARAMETER)
    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded exact CRT mutation and quartic-point experiment; numerical "
            "height rank is triage only and target claims require the stored "
            "finite-reduction and conductor certificates"
        ),
        "primary_source": PRIMARY_SOURCE,
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "certified_hits": certified_hits,
        },
        "profile_correction": {
            "lead_parameter": rational_to_string(LEAD_PARAMETER),
            "actual_residues_mod_5_7_11_13": list(LEAD_PROFILE.residues),
            "actual_crt_residue_mod_5005": LEAD_PROFILE.crt_residue,
            "separate_alternate_branch_residues": list(REQUESTED_PROFILE.residues),
            "separate_alternate_branch_crt_residue_mod_5005": REQUESTED_PROFILE.crt_residue,
        },
        "local_root_union": {
            "5": [1, 4],
            "7": [0, 3, 4],
            "11": [2, 9],
            "13": [3, 10],
            "exact_verification": local_verification,
        },
        "lead_calibration": {
            "constructor_parameter": rational_to_string(LEAD_PARAMETER),
            "radical_proxy": lead_proxy,
            "known_log_conductor_from_prior_exact_replay": "138.23366243820961694",
            "proxy_matches_known_log_conductor_to_displayed_precision": abs(
                Decimal(str(lead_proxy["log_radical_upper_proxy"]))
                - Decimal("138.23366243820961694")
            ) < Decimal("1e-12"),
        },
        "mutation_space": {
            "profiles": [
                {
                    "label": profile.label,
                    "residues_mod_5_7_11_13": list(profile.residues),
                    "crt_residue_mod_5005": profile.crt_residue,
                }
                for profile in DEFAULT_PROFILES
            ],
            "regions": [region.__dict__ for region in DEFAULT_REGIONS],
            "earlier_scan_exclusion": (
                "the denominator<=200 slice excludes |numerator|<=10000; the "
                "second slice uses previously unsearched denominators 201..1000"
            ),
            "conductor_proxy_limit": str(TARGET_LOG_CONDUCTOR),
            "candidate_count_after_exact_deduplication_and_proxy_filter": len(candidates),
            "height_box_candidate_count": len(height_box_candidates),
            "rank_friendly_constraint_candidate_count": len(constraint_candidates),
            "rank_friendly_constraint": {
                "good_primes": [31, 43],
                "coefficient_radius": 40,
                "selection_rule": (
                    "all nonsingular residues attaining the minimum a_p at "
                    "each good prime"
                ),
                "lookup_tables": {
                    str(prime): list(best_good_prime_residues(prime))
                    for prime in (31, 43)
                },
            },
        },
        "stages": stage_records,
        "deep_escalation": {
            "height_bound": args.deep_height,
            "rank_threshold": args.deep_rank_threshold,
            "keep_count": args.deep_keep,
            "candidate_ids": [result.candidate.identifier for result in deep_results],
            "records": [result_record(result, include_points=True) for result in deep_results],
            "rank_replays": {
                result.candidate.identifier: {
                    key: value
                    for key, value in rank_records[result.candidate.identifier].items()
                    if key != "selected_exact_points"
                }
                for result in deep_results
            },
        },
        "final_exact_certificates": exact_certificates,
        "final_conductor_replays": conductor_replays,
        "summary": {
            "mutation_candidates": len(candidates),
            "deepest_stable_numerical_rank": max(
                (
                    int(record["stable_numerical_rank"])
                    for record in rank_records.values()
                    if record.get("status") == "completed"
                ),
                default=0,
            ),
            "deepest_candidate_ids": [
                result.candidate.identifier for result in final_pool
            ],
            "certified_target_hit": bool(certified_hits),
        },
        "bounds": {
            "stage_quartic_naive_height_bounds": list(args.stage_heights),
            "stage_keep_counts": list(args.keep_counts),
            "stage_search_timeouts_seconds_per_candidate": list(args.search_timeouts),
            "deep_height_bound": args.deep_height,
            "deep_timeout_seconds_per_candidate": args.deep_timeout,
            "height_precisions": list(args.height_precisions),
            "height_timeout_seconds_per_replay": args.height_timeout,
            "conductor_timeout_seconds_per_candidate": args.conductor_timeout,
            "saturation_timeout_seconds": args.saturation_timeout,
            "certificate_prime_bound": args.certificate_prime_bound,
            "pari_stack_bytes": args.stack_bytes,
            "process_policy": (
                "all PARI jobs are synchronous foreground subprocesses; every "
                "call uses subprocess.run with a finite timeout"
            ),
        },
        "interpretation": (
            "A bounded quartic box that finds few points is not a rank upper "
            "bound.  Radical proxy values are selection features, not conductor "
            "certificates.  Only exact finite-reduction records certify rank."
        ),
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "invocation": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "reproducing_command": REPRODUCING_COMMAND,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True, default=str) + "\n")
    print(f"wrote {args.output}")
    print(json.dumps(artifact["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
