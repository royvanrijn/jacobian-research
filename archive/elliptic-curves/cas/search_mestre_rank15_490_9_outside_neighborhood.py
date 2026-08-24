#!/usr/bin/env python3
"""Conductor-first outside-box neighborhood of the rank-15 T=490/9 fibre.

The exact anchor has Mestre roots (0,7,121,128,183,194), constructor
parameter 490/9, certified algebraic rank at least 15, and log conductor
124.723916....  The earlier rational frontier exhaustively scored primitive
positive a/b with a<=30000 and b<=1000.  This standalone continuation rejects
that complete box and the anchor before doing any arithmetic.

The frozen outside-box population is the exact union of four constructions:

* nearest-lattice parameters for 552<=b<=4096 and offsets -4..4;
* the two determinant-one Farey rays adjacent to 490/9, m=61..2048;
* twelve discriminant-power congruence moduli on a fixed denominator grid;
* single and pair CRT classes for the four most favourable anchor traces in a
  fresh discovery-prime band, on a disjoint fixed denominator grid.

Every candidate gets fresh discovery/held local-trace scores and the exact
degree-20 homogeneous discriminant with all valuations through 997.  The
closed conductor panel consists of every rigorous radical proxy below 190
plus fixed trace/diversity quotas.  Conductors are computed before any point
search.  Completed curves with log conductor below 190 enter nested H=50000,
H=250000,H=1000000 tiers.  A stable numerical rank at least 16 immediately
triggers exact mod-3 finite-reduction certification.

No fixed-fibre alternate cover, skew chart, or hidden-direction search is
performed here.  Stable numerical rank is triage, not a theorem.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import gcd, log
from pathlib import Path
import platform
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

from mestre_root_tuples import SixRootMestreConstruction
from pari_bridge import pari_version
from search_mestre_0430313946_frontier import exact_point_stage
from search_mestre_rank14_pair_rational_frontier import (
    complete_radical,
    exact_local_trace_projective,
    family_coefficients,
    primes_up_to,
)
from search_mestre_root_tuple_scale import (
    CappedProcessTimeout,
    capped_minimal_curve_data,
    sha256_file,
)
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

ROOTS = (0, 7, 121, 128, 183, 194)
T0 = Q(490, 9)
CONSTRUCTION = SixRootMestreConstruction(tuple(Q(root) for root in ROOTS))
TARGET_LOG_CONDUCTOR = Decimal("182.72")
POINT_LOG_CONDUCTOR_GATE = Decimal("190")
PRIOR_NUMERATOR_BOUND = 30_000
PRIOR_DENOMINATOR_BOUND = 1_000

CERTIFICATE_SCRIPT = Path("elliptic-curves/cas/certify_mestre_rank15_490_9.py")
CERTIFICATE_ARTIFACT = Path(
    "artifacts/generated-results/elliptic_mestre_rank15_490_9.json"
)
EXPECTED_CERTIFICATE_SCRIPT_SHA256 = (
    "622f0d563f7b34d9a06635d79da992066b86b3797dff3f496c7b7959c8f7bd12"
)
EXPECTED_CERTIFICATE_ARTIFACT_SHA256 = (
    "50b2b9c8bd24bcb5533534446af6404f3a9a761b5f33e0e28e04dc572227f950"
)
FRONTIER_SCRIPT = Path(
    "elliptic-curves/cas/search_mestre_rank14_pair_rational_frontier.py"
)
FRONTIER_ARTIFACT = Path(
    "artifacts/generated-results/elliptic_mestre_rank14_pair_rational_frontier.json"
)
EXPECTED_FRONTIER_SCRIPT_SHA256 = (
    "2f6251c67e2eb3cee2eca37d7e866913e9d5de73d30e3bfcb253641454d40d5f"
)
EXPECTED_FRONTIER_ARTIFACT_SHA256 = (
    "87e2d278cc1ee0653d1a4f871c1e34ed3d03babe1c1cd2ffe6712b7608efaee7"
)

DISCOVERY_PRIMES = (
    379, 383, 389, 397, 401, 409, 419, 421, 431,
    433, 439, 443, 449, 457, 461, 463, 467,
)
HELD_PRIMES = (
    479, 487, 491, 499, 503, 509, 521, 523,
    541, 547, 557, 563, 569, 571, 577,
)
ANCHOR_TRACE_PRIMES = (419, 433, 401, 463)
TRACE_MUTATION_MODULI = ANCHOR_TRACE_PRIMES + tuple(
    ANCHOR_TRACE_PRIMES[left] * ANCHOR_TRACE_PRIMES[right]
    for left in range(len(ANCHOR_TRACE_PRIMES))
    for right in range(left + 1, len(ANCHOR_TRACE_PRIMES))
)
POWER_MUTATION_MODULI = (
    125, 2401, 83521, 22801, 1225, 1445,
    14161, 2567, 595, 5285, 17969, 89845,
)
TRIAL_PRIME_BOUND = 997
TRIAL_PRIMES = primes_up_to(TRIAL_PRIME_BOUND)
DISCRIMINANT_POLYNOMIAL = tuple(
    int(value) for value in CONSTRUCTION.primitive_discriminant_polynomial
)
if len(DISCRIMINANT_POLYNOMIAL) != 21:
    raise AssertionError("the degree-20 discriminant changed")

EXPECTED_GENERATOR_COUNT = 35_507
EXPECTED_GENERATOR_SHA256 = (
    "333ccb6ebf567296965613200e713af168f105fa9ee3e84dbb24f028fa74e368"
)
COMPLEMENTARY_ANNULUS_MANIFEST_SHA256 = (
    "6c2fd81b95bb0d7c0531bbccaeff6bcb7fecd1353420febab32492dd04aef4f2"
)
COMPLEMENTARY_ANNULUS_PARAMETER_COUNT = 168_097
PROXY_CONDUCTOR_GATE = 190.0
EXTRA_PROXY_GATE = 210.0
MAX_CONDUCTOR_PANEL = 48
H50000_KEEP = 16
H250000_KEEP = 6
H1000000_KEEP = 2
FINITE_CERTIFICATE_TRIGGER = 16
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/"
    "elliptic_mestre_rank15_490_9_outside_neighborhood.json"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_mestre_rank15_490_9_outside_neighborhood.py"
)


@dataclass(frozen=True)
class Candidate:
    parameter: Fraction
    sources: tuple[str, ...]
    discovery_units: int
    held_units: int
    feature: dict[str, Any]

    @property
    def combined_units(self) -> int:
        return self.discovery_units + self.held_units


def rational_string(value: Fraction) -> str:
    value = Q(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def file_record(path: Path) -> dict[str, str]:
    return {"path": str(path), "sha256": sha256_file(path)}


def validate_pinned_anchor(root: Path) -> dict[str, Any]:
    paths = (
        (CERTIFICATE_SCRIPT, EXPECTED_CERTIFICATE_SCRIPT_SHA256),
        (CERTIFICATE_ARTIFACT, EXPECTED_CERTIFICATE_ARTIFACT_SHA256),
        (FRONTIER_SCRIPT, EXPECTED_FRONTIER_SCRIPT_SHA256),
        (FRONTIER_ARTIFACT, EXPECTED_FRONTIER_ARTIFACT_SHA256),
    )
    for relative, expected in paths:
        if sha256_file(root / relative) != expected:
            raise AssertionError(f"the pinned input changed: {relative}")
    artifact = json.loads((root / CERTIFICATE_ARTIFACT).read_text())
    if (
        tuple(artifact["curve"]["roots"]) != ROOTS
        or Q(artifact["curve"]["parameter"]) != T0
        or artifact["claim"]["certified_algebraic_rank_lower_bound"] != 15
    ):
        raise AssertionError("the rank-15 anchor certificate changed")
    expected_coefficients = tuple(Q(value) for value in artifact["curve"]["weierstrass_coefficients"])
    if family_coefficients(1, T0) != expected_coefficients:
        raise AssertionError("the anchor short model changed")
    return artifact


def nearest_integer(numerator: int, denominator: int) -> int:
    if denominator <= 0:
        raise ValueError("nearest_integer expects a positive denominator")
    # floor(x+1/2), including for the signed CRT recentering quotient.
    return (2 * numerator + denominator) // (2 * denominator)


def belongs_to_prior_box(parameter: Fraction) -> bool:
    parameter = abs(Q(parameter))
    return (
        parameter.numerator <= PRIOR_NUMERATOR_BOUND
        and parameter.denominator <= PRIOR_DENOMINATOR_BOUND
    )


def generator_digest(population: dict[Fraction, set[str]]) -> str:
    payload = "".join(
        f"{rational_string(parameter)}|{','.join(sorted(population[parameter]))}\n"
        for parameter in sorted(population)
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def generate_population() -> dict[Fraction, set[str]]:
    population: dict[Fraction, set[str]] = {}

    def add(numerator: int, denominator: int, source: str) -> None:
        if numerator <= 0 or denominator <= 0:
            return
        parameter = Q(numerator, denominator)
        if parameter == T0 or belongs_to_prior_box(parameter):
            return
        population.setdefault(parameter, set()).add(source)

    for denominator in range(552, 4097):
        center = nearest_integer(490 * denominator, 9)
        for offset in range(-4, 5):
            add(center + offset, denominator, f"gauss-near:{offset:+d}")

    for multiplier in range(61, 2049):
        add(
            490 * multiplier + 381,
            9 * multiplier + 7,
            "farey-left:det+1",
        )
        add(
            490 * multiplier + 109,
            9 * multiplier + 2,
            "farey-right:det-1",
        )

    for modulus in POWER_MUTATION_MODULI:
        inverse_nine = pow(9, -1, modulus)
        for grid_index in range(256):
            denominator = 552 + 37 * grid_index
            if gcd(denominator, modulus) != 1:
                continue
            residue = 490 * inverse_nine * denominator % modulus
            center = nearest_integer(490 * denominator, 9)
            quotient = nearest_integer(center - residue, modulus)
            for shift in (-1, 0, 1):
                add(
                    residue + (quotient + shift) * modulus,
                    denominator,
                    f"discriminant-power:M={modulus}:shift={shift:+d}",
                )

    for modulus in TRACE_MUTATION_MODULI:
        inverse_nine = pow(9, -1, modulus)
        for grid_index in range(256):
            denominator = 553 + 41 * grid_index
            if gcd(denominator, modulus) != 1:
                continue
            residue = 490 * inverse_nine * denominator % modulus
            center = nearest_integer(490 * denominator, 9)
            quotient = nearest_integer(center - residue, modulus)
            for shift in (-1, 0, 1):
                add(
                    residue + (quotient + shift) * modulus,
                    denominator,
                    f"local-trace:M={modulus}:shift={shift:+d}",
                )

    if len(population) != EXPECTED_GENERATOR_COUNT:
        raise AssertionError("the frozen outside-box generator count changed")
    if generator_digest(population) != EXPECTED_GENERATOR_SHA256:
        raise AssertionError("the frozen outside-box generator digest changed")
    if any(parameter <= 0 or parameter == T0 for parameter in population):
        raise AssertionError("the generator leaked the anchor or a nonpositive T")
    if any(belongs_to_prior_box(parameter) for parameter in population):
        raise AssertionError("the generator leaked the prior rational box")
    return population


def llround(value: float) -> int:
    return int(value + 0.5) if value >= 0 else int(value - 0.5)


def local_contribution_units(trace: int | None, prime: int) -> int:
    if trace is None:
        return 0
    return llround(
        ((2 - trace) / (prime + 1 - trace)) * log(float(prime)) * 1.0e12
    )


def score_string(units: int) -> str:
    sign = "-" if units < 0 else ""
    absolute = abs(units)
    return (
        f"{sign}{absolute // 1_000_000_000_000}."
        f"{absolute % 1_000_000_000_000:012d}"
    )


def build_trace_tables(
    primes: Sequence[int],
) -> dict[int, tuple[int | None, ...]]:
    tables = {}
    for prime in primes:
        tables[prime] = tuple(
            exact_local_trace_projective(1, residue, 1, prime)
            for residue in range(prime)
        ) + (exact_local_trace_projective(1, 1, 0, prime),)
    return tables


def score_parameter(
    parameter: Fraction,
    primes: Sequence[int],
    tables: dict[int, tuple[int | None, ...]],
) -> int:
    parameter = Q(parameter)
    answer = 0
    for prime in primes:
        index = (
            prime
            if parameter.denominator % prime == 0
            else parameter.numerator * pow(parameter.denominator, -1, prime) % prime
        )
        answer += local_contribution_units(tables[prime][index], prime)
    return answer


def validate_anchor_trace_prime_selection(
    tables: dict[int, tuple[int | None, ...]]
) -> None:
    contributions = []
    for prime in DISCOVERY_PRIMES:
        index = T0.numerator * pow(T0.denominator, -1, prime) % prime
        contributions.append(
            (local_contribution_units(tables[prime][index], prime), prime)
        )
    selected = tuple(prime for _, prime in sorted(contributions, reverse=True)[:4])
    if selected != ANCHOR_TRACE_PRIMES:
        raise AssertionError("the pinned anchor-trace prime selection changed")


def homogeneous_discriminant(parameter: Fraction) -> int:
    parameter = abs(Q(parameter))
    numerator, denominator = parameter.numerator, parameter.denominator
    value = sum(
        coefficient
        * numerator**power
        * denominator ** (20 - power)
        for power, coefficient in enumerate(DISCRIMINANT_POLYNOMIAL)
    )
    if value == 0:
        raise ValueError("singular specialization")
    return value


def discriminant_feature(parameter: Fraction) -> dict[str, Any]:
    absolute = abs(homogeneous_discriminant(parameter))
    remaining = absolute
    known_radical = 1
    known_powerful = 1
    valuations = []
    for prime in TRIAL_PRIMES:
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        if exponent:
            valuations.append([prime, exponent])
            known_radical *= prime
            if exponent > 1:
                known_powerful *= prime ** (exponent - 1)
    denominator_radical = complete_radical(Q(parameter).denominator)
    radical_upper = known_radical * remaining * denominator_radical
    return {
        "absolute_homogeneous_discriminant": str(absolute),
        "trial_division_prime_bound": TRIAL_PRIME_BOUND,
        "small_prime_valuations": valuations,
        "residual_cofactor": str(remaining),
        "residual_cofactor_bit_length": remaining.bit_length(),
        "known_discriminant_radical": str(known_radical),
        "known_powerful_part": str(known_powerful),
        "denominator_radical": denominator_radical,
        "combined_radical_upper_bound": str(radical_upper),
        "log_combined_radical_upper_bound": log(radical_upper),
        "upper_bound_semantics": (
            "rad(residual)<=residual exactly; valuations p<=997 and rad(b) "
            "are exact"
        ),
    }


def analyze_population(
    population: dict[Fraction, set[str]]
) -> tuple[tuple[Candidate, ...], dict[str, Any]]:
    tables = build_trace_tables(DISCOVERY_PRIMES + HELD_PRIMES)
    validate_anchor_trace_prime_selection(tables)
    candidates = []
    digest = hashlib.sha256()
    source_counts: dict[str, int] = defaultdict(int)
    for index, parameter in enumerate(sorted(population), 1):
        sources = tuple(sorted(population[parameter]))
        for family in {source.split(":", 1)[0] for source in sources}:
            source_counts[family] += 1
        discovery = score_parameter(parameter, DISCOVERY_PRIMES, tables)
        held = score_parameter(parameter, HELD_PRIMES, tables)
        feature = discriminant_feature(parameter)
        candidate = Candidate(parameter, sources, discovery, held, feature)
        candidates.append(candidate)
        digest.update(
            (
                f"{rational_string(parameter)}|{','.join(sources)}|"
                f"{discovery}|{held}|"
                f"{feature['combined_radical_upper_bound']}|"
                f"{feature['known_powerful_part']}\n"
            ).encode()
        )
        if index % 5_000 == 0:
            print(f"feature {index}/{len(population)}", flush=True)
    anchor_discovery = score_parameter(T0, DISCOVERY_PRIMES, tables)
    anchor_held = score_parameter(T0, HELD_PRIMES, tables)
    return tuple(candidates), {
        "candidate_count": len(candidates),
        "generator_sha256": generator_digest(population),
        "feature_population_sha256": digest.hexdigest(),
        "source_family_counts_after_deduplication": dict(sorted(source_counts.items())),
        "anchor_discovery_score": score_string(anchor_discovery),
        "anchor_held_score": score_string(anchor_held),
        "anchor_combined_score": score_string(anchor_discovery + anchor_held),
        "selected_anchor_trace_primes": list(ANCHOR_TRACE_PRIMES),
        "fresh_discovery_primes": list(DISCOVERY_PRIMES),
        "fresh_held_primes": list(HELD_PRIMES),
    }


def source_families(candidate: Candidate) -> set[str]:
    return {source.split(":", 1)[0] for source in candidate.sources}


def select_conductor_panel(
    candidates: Sequence[Candidate],
) -> tuple[tuple[Candidate, ...], dict[str, Any]]:
    reasons: dict[Fraction, set[str]] = defaultdict(set)
    by_parameter = {candidate.parameter: candidate for candidate in candidates}
    proxy_gate = [
        candidate
        for candidate in candidates
        if candidate.feature["log_combined_radical_upper_bound"]
        < PROXY_CONDUCTOR_GATE
    ]
    for candidate in proxy_gate:
        reasons[candidate.parameter].add("all-rigorous-radical-proxy-below-190")

    extras = [
        candidate
        for candidate in candidates
        if candidate.feature["log_combined_radical_upper_bound"] < EXTRA_PROXY_GATE
    ]
    orders = {
        "highest-combined-fresh-trace": sorted(
            extras,
            key=lambda candidate: (
                -candidate.combined_units,
                candidate.feature["log_combined_radical_upper_bound"],
                candidate.parameter,
            ),
        )[:8],
        "highest-held-fresh-trace": sorted(
            extras,
            key=lambda candidate: (
                -candidate.held_units,
                -candidate.discovery_units,
                candidate.feature["log_combined_radical_upper_bound"],
            ),
        )[:8],
        "largest-known-powerful-part": sorted(
            extras,
            key=lambda candidate: (
                -int(candidate.feature["known_powerful_part"]),
                candidate.feature["log_combined_radical_upper_bound"],
            ),
        )[:6],
    }
    for label, order in orders.items():
        for candidate in order:
            reasons[candidate.parameter].add(label)
    for family in ("gauss-near", "farey-left", "farey-right", "discriminant-power", "local-trace"):
        order = sorted(
            (candidate for candidate in extras if family in source_families(candidate)),
            key=lambda candidate: (
                -candidate.combined_units,
                candidate.feature["log_combined_radical_upper_bound"],
            ),
        )[:2]
        for candidate in order:
            reasons[candidate.parameter].add(f"source-diversity:{family}")

    mandatory = {candidate.parameter for candidate in proxy_gate}
    optional = sorted(
        (parameter for parameter in reasons if parameter not in mandatory),
        key=lambda parameter: (
            -len(reasons[parameter]),
            -by_parameter[parameter].combined_units,
            by_parameter[parameter].feature["log_combined_radical_upper_bound"],
            parameter,
        ),
    )
    room = MAX_CONDUCTOR_PANEL - len(mandatory)
    if room < 0:
        raise AssertionError("the mandatory proxy population exceeded the panel cap")
    retained = mandatory | set(optional[:room])
    selected = tuple(
        sorted(
            (by_parameter[parameter] for parameter in retained),
            key=lambda candidate: (
                candidate.feature["log_combined_radical_upper_bound"],
                -candidate.combined_units,
                candidate.parameter,
            ),
        )
    )
    digest = hashlib.sha256()
    for candidate in selected:
        digest.update(
            (
                f"{rational_string(candidate.parameter)}|"
                f"{','.join(sorted(reasons[candidate.parameter]))}\n"
            ).encode()
        )
    return selected, {
        "all_proxy_below_190_retained": True,
        "proxy_below_190_count": len(proxy_gate),
        "extra_proxy_gate": EXTRA_PROXY_GATE,
        "maximum_panel_size": MAX_CONDUCTOR_PANEL,
        "selected_count": len(selected),
        "selected_sha256": digest.hexdigest(),
        "reasons": {
            rational_string(candidate.parameter): sorted(reasons[candidate.parameter])
            for candidate in selected
        },
    }


def conductor_worker(
    candidate: Candidate, *, timeout: float, stack_bytes: int
) -> tuple[Fraction, dict[str, Any]]:
    try:
        data = capped_minimal_curve_data(
            family_coefficients(1, candidate.parameter),
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        return candidate.parameter, {
            "status": "completed exact conductor",
            **data,
            "below_strict_182_72_target": (
                Decimal(data["log_conductor"]) < TARGET_LOG_CONDUCTOR
            ),
            "below_point_stage_190_gate": (
                Decimal(data["log_conductor"]) < POINT_LOG_CONDUCTOR_GATE
            ),
        }
    except CappedProcessTimeout:
        return candidate.parameter, {
            "status": "timeout",
            "timeout_seconds": timeout,
            "retried": False,
        }
    except Exception as error:
        return candidate.parameter, {
            "status": "error",
            "error": str(error)[:1000],
            "retried": False,
        }


def candidate_summary(candidate: Candidate) -> dict[str, Any]:
    return {
        "parameter": rational_string(candidate.parameter),
        "numerator": candidate.parameter.numerator,
        "denominator": candidate.parameter.denominator,
        "sources": list(candidate.sources),
        "discovery_score": score_string(candidate.discovery_units),
        "held_score": score_string(candidate.held_units),
        "combined_score": score_string(candidate.combined_units),
        "discriminant_feature": candidate.feature,
    }


def point_worker(
    candidate: Candidate,
    *,
    height_bound: int,
    point_timeout: float,
    height_timeout: float,
    ellrank_timeout: float,
    stack_bytes: int,
    mapping_cap: int,
) -> tuple[Fraction, dict[str, Any], tuple[tuple[Fraction, Fraction], ...]]:
    try:
        stage, subset = exact_point_stage(
            CONSTRUCTION,
            candidate.parameter,
            family_coefficients(1, candidate.parameter),
            height_bound=height_bound,
            point_timeout=point_timeout,
            height_timeout=height_timeout,
            ellrank_timeout=ellrank_timeout,
            stack_bytes=stack_bytes,
            mapping_cap=mapping_cap,
        )
        return candidate.parameter, stage, subset
    except CappedProcessTimeout:
        return candidate.parameter, {
            "status": "timeout",
            "timeout_seconds": point_timeout,
            "same_height_retry": False,
        }, ()
    except Exception as error:
        return candidate.parameter, {
            "status": "error",
            "error": str(error)[:1000],
            "same_height_retry": False,
        }, ()


def write_artifact(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--conductor-timeout", type=float, default=30.0)
    parser.add_argument("--h50000-timeout", type=float, default=25.0)
    parser.add_argument("--h250000-timeout", type=float, default=35.0)
    parser.add_argument("--h1000000-timeout", type=float, default=50.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--ellrank-timeout", type=float, default=12.0)
    parser.add_argument("--mapping-cap", type=int, default=512)
    parser.add_argument("--certificate-prime-bound", type=int, default=499)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--output", type=Path, default=root / DEFAULT_OUTPUT)
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if not 1 <= args.workers <= 8:
        raise SystemExit("workers must lie in [1,8]")
    timeouts = (
        args.conductor_timeout,
        args.h50000_timeout,
        args.h250000_timeout,
        args.h1000000_timeout,
        args.height_timeout,
        args.ellrank_timeout,
    )
    if min(timeouts) <= 0 or max(timeouts) > 60:
        raise SystemExit("all one-shot subprocess caps must lie in (0,60]")
    if args.mapping_cap < 32 or args.mapping_cap > 1024:
        raise SystemExit("mapping cap must lie in [32,1024]")
    if not 211 <= args.certificate_prime_bound <= 2000:
        raise SystemExit("certificate prime bound must lie in [211,2000]")


def main() -> None:
    args = build_parser().parse_args()
    validate_args(args)
    started = time.monotonic()
    script_path = Path(__file__).resolve()
    root = script_path.parents[2]
    anchor = validate_pinned_anchor(root)
    population = generate_population()
    candidates, population_record = analyze_population(population)
    by_parameter = {candidate.parameter: candidate for candidate in candidates}
    conductor_panel, selection = select_conductor_panel(candidates)
    print(
        f"population={len(population)} conductor_panel={len(conductor_panel)} "
        f"proxy<190={selection['proxy_below_190_count']}",
        flush=True,
    )

    conductor_records: dict[Fraction, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(
                conductor_worker,
                candidate,
                timeout=args.conductor_timeout,
                stack_bytes=args.stack_bytes,
            )
            for candidate in conductor_panel
        ]
        for completed, future in enumerate(as_completed(futures), 1):
            parameter, record = future.result()
            conductor_records[parameter] = record
            print(
                f"conductor {completed}/{len(futures)} "
                f"T={rational_string(parameter)} status={record['status']}",
                flush=True,
            )

    point_eligible = [
        by_parameter[parameter]
        for parameter, record in conductor_records.items()
        if record["status"] == "completed exact conductor"
        and Decimal(record["log_conductor"]) < POINT_LOG_CONDUCTOR_GATE
    ]
    point_eligible.sort(
        key=lambda candidate: (
            -candidate.combined_units,
            Decimal(conductor_records[candidate.parameter]["log_conductor"]),
            candidate.feature["log_combined_radical_upper_bound"],
            candidate.parameter,
        )
    )
    h50000_panel = tuple(point_eligible[:H50000_KEEP])

    checkpoint = {
        "schema_version": 1,
        "status": "conductor-first checkpoint before point search",
        "scope": {
            "roots": list(ROOTS),
            "anchor_parameter": rational_string(T0),
            "prior_box_exclusion": {
                "positive_primitive_numerator_maximum": PRIOR_NUMERATOR_BOUND,
                "positive_primitive_denominator_maximum": PRIOR_DENOMINATOR_BOUND,
                "all_reduced_parameters_in_box_excluded": True,
            },
            "generator_count": len(population),
            "generator_sha256": generator_digest(population),
        },
        "population": population_record,
        "conductor_selection": selection,
        "conductor_records": {
            rational_string(parameter): record
            for parameter, record in sorted(conductor_records.items())
        },
        "point_protocol": {
            "H50000_keep": H50000_KEEP,
            "H250000_keep": H250000_KEEP,
            "H1000000_keep": H1000000_KEEP,
            "finite_certificate_trigger": FINITE_CERTIFICATE_TRIGGER,
        },
        "checkpoint_written_before_any_point_search": True,
    }
    write_artifact(args.output, checkpoint)
    print(f"checkpoint written {args.output}", flush=True)

    stages: dict[str, dict[Fraction, dict[str, Any]]] = {}
    certified_signals: list[dict[str, Any]] = []

    def run_stage(
        name: str,
        panel: Sequence[Candidate],
        *,
        height_bound: int,
        point_timeout: float,
    ) -> None:
        records: dict[Fraction, dict[str, Any]] = {}
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = [
                executor.submit(
                    point_worker,
                    candidate,
                    height_bound=height_bound,
                    point_timeout=point_timeout,
                    height_timeout=args.height_timeout,
                    ellrank_timeout=args.ellrank_timeout,
                    stack_bytes=args.stack_bytes,
                    mapping_cap=args.mapping_cap,
                )
                for candidate in panel
            ]
            for completed, future in enumerate(as_completed(futures), 1):
                parameter, stage, subset = future.result()
                records[parameter] = stage
                rank = stage.get("stable_numerical_rank")
                print(
                    f"{name} {completed}/{len(futures)} "
                    f"T={rational_string(parameter)} status={stage['status']} rank={rank}",
                    flush=True,
                )
                if stage["status"] == "completed" and int(rank) >= FINITE_CERTIFICATE_TRIGGER:
                    print(
                        f"ALERT stable-rank-{rank} T={rational_string(parameter)} "
                        f"stage={name}",
                        flush=True,
                    )
                    certificate = mod3_independence_certificate(
                        family_coefficients(1, parameter),
                        subset,
                        prime_bound=args.certificate_prime_bound,
                    )
                    stage["immediate_finite_reduction_certificate"] = certificate
                    certified_signals.append(
                        {
                            "parameter": rational_string(parameter),
                            "stage": name,
                            "stable_numerical_rank": rank,
                            "certified_algebraic_rank_lower_bound": certificate[
                                "certified_algebraic_rank_lower_bound"
                            ],
                        }
                    )
                    print(
                        f"ALERT exact-rank-lower-bound="
                        f"{certificate['certified_algebraic_rank_lower_bound']} "
                        f"T={rational_string(parameter)}",
                        flush=True,
                    )
        stages[name] = records

    run_stage(
        "H50000",
        h50000_panel,
        height_bound=50_000,
        point_timeout=args.h50000_timeout,
    )
    h250000_panel = tuple(
        sorted(
            (
                candidate
                for candidate in h50000_panel
                if stages["H50000"].get(candidate.parameter, {}).get("status")
                == "completed"
            ),
            key=lambda candidate: (
                -int(stages["H50000"][candidate.parameter]["stable_numerical_rank"]),
                -candidate.combined_units,
                Decimal(conductor_records[candidate.parameter]["log_conductor"]),
            ),
        )[:H250000_KEEP]
    )
    run_stage(
        "H250000",
        h250000_panel,
        height_bound=250_000,
        point_timeout=args.h250000_timeout,
    )
    h1000000_panel = tuple(
        sorted(
            (
                candidate
                for candidate in h250000_panel
                if stages["H250000"].get(candidate.parameter, {}).get("status")
                == "completed"
            ),
            key=lambda candidate: (
                -int(stages["H250000"][candidate.parameter]["stable_numerical_rank"]),
                -candidate.combined_units,
                Decimal(conductor_records[candidate.parameter]["log_conductor"]),
            ),
        )[:H1000000_KEEP]
    )
    run_stage(
        "H1000000",
        h1000000_panel,
        height_bound=1_000_000,
        point_timeout=args.h1000000_timeout,
    )

    strict_subtarget = sorted(
        parameter
        for parameter, record in conductor_records.items()
        if record.get("below_strict_182_72_target", False)
    )
    maximum_stable_rank = max(
        (
            int(stage["stable_numerical_rank"])
            for records in stages.values()
            for stage in records.values()
            if stage.get("status") == "completed"
        ),
        default=0,
    )
    top_candidates = sorted(
        candidates,
        key=lambda candidate: (
            candidate.feature["log_combined_radical_upper_bound"],
            -candidate.combined_units,
            candidate.parameter,
        ),
    )[:64]
    artifact = {
        "schema_version": 1,
        "status": "completed outside-box conductor-first staged neighborhood",
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_below": "182.72",
            "alternative_rank_at_least": 30,
            "intermediate_exact_signal_trigger": 16,
        },
        "anchor": {
            "roots": list(ROOTS),
            "parameter": rational_string(T0),
            "certified_algebraic_rank_lower_bound": anchor["claim"][
                "certified_algebraic_rank_lower_bound"
            ],
            "log_conductor": anchor["curve"]["log_conductor"],
            "fixed_fibre_followup_excluded": True,
            "pinned_inputs": [
                file_record(root / CERTIFICATE_SCRIPT),
                file_record(root / CERTIFICATE_ARTIFACT),
                file_record(root / FRONTIER_SCRIPT),
                file_record(root / FRONTIER_ARTIFACT),
            ],
        },
        "scope": {
            "T_sign_quotient": "positive representative only; family is even in T",
            "prior_box_exclusion": {
                "numerator_maximum": PRIOR_NUMERATOR_BOUND,
                "denominator_maximum": PRIOR_DENOMINATOR_BOUND,
                "complete_box_excluded_after_exact_reduction": True,
            },
            "anchor_and_sign_mate_excluded": True,
            "generator_count": len(population),
            "generator_sha256": generator_digest(population),
            "gauss_near": {
                "raw_denominator_interval": [552, 4096],
                "nearest_numerator_offsets": list(range(-4, 5)),
            },
            "continued_fraction_farey_rays": {
                "multiplier_interval": [61, 2048],
                "left_vector": [381, 7],
                "right_vector": [109, 2],
                "determinants_with_anchor": [1, -1],
            },
            "discriminant_power_grid": {
                "raw_denominators": "552+37*j, 0<=j<=255",
                "moduli": list(POWER_MUTATION_MODULI),
                "nearest_residue_class_shifts": [-1, 0, 1],
            },
            "local_trace_grid": {
                "raw_denominators": "553+41*j, 0<=j<=255",
                "anchor_trace_primes": list(ANCHOR_TRACE_PRIMES),
                "single_and_pair_moduli": list(TRACE_MUTATION_MODULI),
                "nearest_residue_class_shifts": [-1, 0, 1],
            },
            "no_fixed_fibre_cover_or_skew_search": True,
            "coordinated_complementary_lane": {
                "raw_denominator_interval": [4097, 16000],
                "absolute_offset_interval": [5, 16],
                "evaluated_unique_primitive_parameter_count": (
                    COMPLEMENTARY_ANNULUS_PARAMETER_COUNT
                ),
                "canonical_manifest_sha256": (
                    COMPLEMENTARY_ANNULUS_MANIFEST_SHA256
                ),
                "overlap_exclusions": (
                    "both Farey rays and every divisor of a C/D raw-grid "
                    "denominator"
                ),
            },
        },
        "population": population_record,
        "top_64_candidates_by_rigorous_proxy": [
            candidate_summary(candidate) for candidate in top_candidates
        ],
        "conductor_selection": {
            **selection,
            "selected_candidates": [
                candidate_summary(candidate) for candidate in conductor_panel
            ],
        },
        "conductors": {
            "records": {
                rational_string(parameter): record
                for parameter, record in sorted(conductor_records.items())
            },
            "strict_subtarget_count": len(strict_subtarget),
            "strict_subtarget_parameters": [
                rational_string(parameter) for parameter in strict_subtarget
            ],
            "completed_below_190_count": len(point_eligible),
        },
        "point_protocol": {
            "conductor_first": True,
            "point_log_conductor_gate": str(POINT_LOG_CONDUCTOR_GATE),
            "H50000_panel": [
                rational_string(candidate.parameter) for candidate in h50000_panel
            ],
            "H250000_panel": [
                rational_string(candidate.parameter) for candidate in h250000_panel
            ],
            "H1000000_panel": [
                rational_string(candidate.parameter) for candidate in h1000000_panel
            ],
            "stage_records": {
                name: {
                    rational_string(parameter): record
                    for parameter, record in sorted(records.items())
                }
                for name, records in stages.items()
            },
            "maximum_stable_numerical_rank": maximum_stable_rank,
            "finite_certificate_trigger": FINITE_CERTIFICATE_TRIGGER,
            "certified_signals": certified_signals,
            "numerical_rank_scope_warning": "triage only unless a finite certificate is present",
        },
        "outcome": {
            "breakthrough_found": any(
                signal["certified_algebraic_rank_lower_bound"] >= 21
                and Decimal(
                    conductor_records[Q(signal["parameter"])]["log_conductor"]
                ) < TARGET_LOG_CONDUCTOR
                for signal in certified_signals
            ),
            "rank_signal_at_least_16": bool(certified_signals),
            "maximum_stable_numerical_rank": maximum_stable_rank,
            "strict_subtarget_conductor_count": len(strict_subtarget),
            "statement": (
                "The frozen outside-box neighborhood completed without an exact "
                "rank-21/subtarget or rank-30 breakthrough."
            ),
            "bounded_search_warning": (
                "This excludes only the declared rational population and nested "
                "point-height panels; it is not a rank upper bound."
            ),
        },
        "checkpoint": {
            "conductor_checkpoint_written_before_point_search": True,
            "final_artifact_replaced_checkpoint_after_all_stages": True,
        },
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "script_path": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "pari": pari_version(),
            "wall_seconds": time.monotonic() - started,
            "workers": args.workers,
            "all_processes_foreground_and_bounded": True,
            "one_attempt_per_declared_subprocess": True,
        },
    }
    write_artifact(args.output, artifact)
    print(
        f"wrote {args.output} strict={len(strict_subtarget)} "
        f"max_rank={maximum_stable_rank} exact_signals={len(certified_signals)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
