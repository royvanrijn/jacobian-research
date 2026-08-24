#!/usr/bin/env python3
"""Exact conductor-first Fermigier high-power CRT--Gauss checkpoint.

The parameter in this file is always the canonical adapter ``u=s/2``.  The
search is a finite reduced-basis ball, not an ordinary rational box.  It is
therefore disjoint in construction and by exact replay from the separately
frozen neighbourhood of ``u=28917/20``.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import gcd
from pathlib import Path
import platform
import subprocess
import time
from typing import Any, Iterable

from ecsearch.crt_lattice import (
    crt,
    enumerate_rational_representatives,
    hensel_lift_roots,
)
from ecsearch.fermigier import (
    FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS,
    fermigier_canonical_coefficients,
)
from ecsearch.fermigier_seed import homogenized_discriminant_factor
from search_mestre_root_tuple_scale import capped_minimal_curve_data


Q = Fraction
PRIMES = (89, 131, 137)
EXPONENT_PROFILES = (
    (2, 2, 2),
    (3, 2, 2),
    (2, 3, 2),
    (2, 2, 3),
    (3, 3, 2),
    (3, 2, 3),
    (2, 3, 3),
    (3, 3, 3),
    (4, 2, 2),
    (2, 4, 2),
    (2, 2, 4),
)
BASE_RADIUS = 64
HIGH_POWER_RADIUS = 16
TRIAL_PRIME_BOUND = 2_000
CONDUCTOR_KEEP = 8
TARGET_LOG_CONDUCTOR = Decimal("182.72")
IMPORTED_SEED = Q(673709, 29965)
AUDIT_RELATIVE = Path(
    "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_rank20_adapter_neighborhood_audit.json"
)
AUDIT_SHA256 = "0eef1ad22211d9b8f6b8cdcec3e1c8829322f2889195a2f1527b03465e799615"
EXPECTED_PRIOR_COUNT = 1_798
EXPECTED_PRIOR_DIGEST = "7c6884a10a434b4d0b3c463ccab97f1d24cf4dbeebfc33e4e250d1ced158bf91"
EXPECTED_RAW_COUNT = 31_524
EXPECTED_FRESH_COUNT = 31_520
EXPECTED_FRESH_DIGEST = "3185367a55b78b6b490eba416198e9f30e19bb38db5f69a73fbb3b4778724aff"
EXPECTED_PRIOR_INTERSECTION = (
    Q(951859, 1297744),
    Q(722299, 169126),
    Q(1174354, 111753),
    Q(673709, 29965),
)


@dataclass(frozen=True)
class Candidate:
    parameter: Fraction
    profiles: tuple[tuple[int, int, int], ...]
    radical: dict[str, Any]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational_string(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def rational_digest(values: Iterable[Fraction]) -> str:
    digest = hashlib.sha256()
    for value in sorted(set(Q(value) for value in values)):
        digest.update((rational_string(value) + "\n").encode())
    return digest.hexdigest()


def denominator_first_rational_digest(values: Iterable[Fraction]) -> str:
    digest = hashlib.sha256()
    for value in sorted(set(Q(value) for value in values), key=lambda item: (item.denominator, item.numerator)):
        digest.update((rational_string(value) + "\n").encode())
    return digest.hexdigest()


def primes_up_to(bound: int) -> tuple[int, ...]:
    sieve = bytearray(b"\x01") * (bound + 1)
    if bound >= 0:
        sieve[0] = 0
    if bound >= 1:
        sieve[1] = 0
    for prime in range(2, int(bound**0.5) + 1):
        if sieve[prime]:
            sieve[prime * prime : bound + 1 : prime] = b"\x00" * (
                (bound - prime * prime) // prime + 1
            )
    return tuple(index for index, flag in enumerate(sieve) if flag)


def generate_population() -> tuple[dict[Fraction, set[tuple[int, int, int]]], list[dict[str, Any]]]:
    population: dict[Fraction, set[tuple[int, int, int]]] = {}
    profile_records: list[dict[str, Any]] = []
    for profile in EXPONENT_PROFILES:
        radius = BASE_RADIUS if profile == (2, 2, 2) else HIGH_POWER_RADIUS
        lifted = tuple(
            tuple(hensel_lift_roots(FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS, prime, exponent))
            for prime, exponent in zip(PRIMES, profile, strict=True)
        )
        raw_count = 0
        for roots in product(*lifted):
            residue, modulus = crt(
                zip(roots, (prime**exponent for prime, exponent in zip(PRIMES, profile, strict=True)), strict=True)
            )
            representatives = enumerate_rational_representatives(
                residue, modulus, coefficient_radius=radius, weights=(1, 1)
            )
            raw_count += len(representatives)
            for representative in representatives:
                parameter = abs(Q(representative.numerator, representative.denominator))
                population.setdefault(parameter, set()).add(profile)
        profile_records.append(
            {
                "profile": list(profile),
                "coefficient_radius": radius,
                "root_counts": [len(values) for values in lifted],
                "root_combinations": len(tuple(product(*lifted))),
                "signed_representatives_before_sign_quotient": raw_count,
            }
        )
    return population, profile_records


def load_prior(root: Path) -> tuple[set[Fraction], dict[str, Any]]:
    path = root / AUDIT_RELATIVE
    observed = sha256_file(path)
    if observed != AUDIT_SHA256:
        raise AssertionError("the pinned canonical-parameter audit changed")
    audit = json.loads(path.read_text())
    prior = {
        Q(record["adapter_u"])
        for record in audit["prior_parameter_exclusion"]["parameters"]
    }
    if len(prior) != EXPECTED_PRIOR_COUNT or denominator_first_rational_digest(prior) != EXPECTED_PRIOR_DIGEST:
        raise AssertionError("the pinned prior canonical-u manifest changed")
    return prior, {
        "audit_path": str(path.relative_to(root)),
        "audit_sha256": observed,
        "canonical_adapter_parameter_count": len(prior),
        "canonical_adapter_parameter_sha256": denominator_first_rational_digest(prior),
    }


def in_broad_neighbourhood(parameter: Fraction) -> bool:
    """Replay both peer populations exactly in reduced positive coordinates."""

    numerator, denominator = parameter.numerator, parameter.denominator
    dense = denominator <= 1_200 and abs(20 * numerator - 28_917 * denominator) <= 40 * denominator
    nearest = (2 * 28_917 * denominator + 20) // 40
    deep = 1_201 <= denominator <= 20_000 and abs(numerator - nearest) <= 64
    return dense or deep


def radical_record(parameter: Fraction, trial_primes: tuple[int, ...]) -> dict[str, Any]:
    homogeneous = abs(homogenized_discriminant_factor(parameter.numerator, parameter.denominator))
    if homogeneous == 0:
        return {"singular": True}
    remaining = homogeneous
    factors: dict[str, int] = {}
    repeated_divisor = 1
    for prime in trial_primes:
        if remaining % prime:
            continue
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        factors[str(prime)] = exponent
        repeated_divisor *= prime ** (exponent - 1)
    radical_upper = homogeneous // repeated_divisor
    return {
        "singular": False,
        "homogenized_discriminant_factor": str(homogeneous),
        "trial_prime_bound": trial_primes[-1],
        "exact_trial_factorization": factors,
        "exact_repeated_prime_divisor": str(repeated_divisor),
        "exact_radical_upper_integer": str(radical_upper),
        "unfactored_cofactor": str(remaining),
        "unfactored_cofactor_bits": remaining.bit_length(),
        "explanation": (
            "rad(H) divides this integer because the unfactored cofactor is charged in full; "
            "this is an exact integer upper proxy, not a conductor bound"
        ),
    }


def conductor_probe(parameter: Fraction, *, timeout: float, stack_bytes: int) -> dict[str, Any]:
    started = time.monotonic()
    try:
        data = capped_minimal_curve_data(
            fermigier_canonical_coefficients(parameter), timeout=timeout, stack_bytes=stack_bytes
        )
        return {
            "status": "completed",
            "wall_seconds": time.monotonic() - started,
            **data,
            "below_strict_log_conductor_target": Decimal(data["log_conductor"]) < TARGET_LOG_CONDUCTOR,
        }
    except (subprocess.TimeoutExpired, TimeoutError) as error:
        return {
            "status": "timeout",
            "error_type": type(error).__name__,
            "timeout_seconds": timeout,
            "wall_seconds": time.monotonic() - started,
        }
    except Exception as error:  # a bounded external factorization can fail in several exact ways
        if "wall cap" in str(error) or "timeout" in str(error).lower():
            return {
                "status": "timeout",
                "error_type": type(error).__name__,
                "timeout_seconds": timeout,
                "wall_seconds": time.monotonic() - started,
            }
        return {
            "status": "error",
            "error_type": type(error).__name__,
            "error": str(error)[:1000],
            "wall_seconds": time.monotonic() - started,
        }


def result_digest(artifact: dict[str, Any]) -> str:
    stable = {
        "population": artifact["population"],
        "prior_exclusion": artifact["prior_exclusion"],
        "selected": [
            {
                "adapter_u": row["adapter_u"],
                "profiles": row["profiles"],
                "radical": row["radical"],
                "conductor": {
                    key: value
                    for key, value in row["conductor"].items()
                    if key not in {"wall_seconds", "timeout_seconds"}
                },
            }
            for row in artifact["selected"]
        ],
        "seed_calibration": {
            **artifact["seed_calibration"],
            "conductor": {
                key: value
                for key, value in artifact["seed_calibration"]["conductor"].items()
                if key not in {"wall_seconds", "timeout_seconds"}
            },
        },
        "outcome": artifact["outcome"],
    }
    return hashlib.sha256(json.dumps(stable, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def run(*, root: Path, conductor_timeout: float, stack_bytes: int, workers: int) -> dict[str, Any]:
    raw, profile_records = generate_population()
    if len(raw) != EXPECTED_RAW_COUNT or IMPORTED_SEED not in raw:
        raise AssertionError("the exact CRT--Gauss population changed")
    prior, prior_record = load_prior(root)
    prior_intersection = tuple(sorted(set(raw) & prior))
    if prior_intersection != EXPECTED_PRIOR_INTERSECTION:
        raise AssertionError("the exact prior intersection changed")
    broad_intersection = tuple(sorted(value for value in raw if in_broad_neighbourhood(value)))
    if broad_intersection:
        raise AssertionError("the CRT root ball overlaps the frozen broad neighbourhood")
    fresh = set(raw) - prior
    if len(fresh) != EXPECTED_FRESH_COUNT or rational_digest(fresh) != EXPECTED_FRESH_DIGEST:
        raise AssertionError("the fresh CRT--Gauss manifest changed")

    trial_primes = primes_up_to(TRIAL_PRIME_BOUND)
    candidates = [
        Candidate(value, tuple(sorted(raw[value])), radical_record(value, trial_primes))
        for value in fresh
    ]
    if any(candidate.radical["singular"] for candidate in candidates):
        raise AssertionError("a CRT representative specialized to a singular fibre")
    candidates.sort(
        key=lambda row: (
            int(row.radical["exact_radical_upper_integer"]),
            max(abs(row.parameter.numerator), row.parameter.denominator),
            row.parameter,
        )
    )
    selected = candidates[:CONDUCTOR_KEEP]
    conductor_records: dict[Fraction, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                conductor_probe,
                candidate.parameter,
                timeout=conductor_timeout,
                stack_bytes=stack_bytes,
            ): candidate.parameter
            for candidate in selected
        }
        futures[
            executor.submit(
                conductor_probe,
                IMPORTED_SEED,
                timeout=conductor_timeout,
                stack_bytes=stack_bytes,
            )
        ] = IMPORTED_SEED
        for future in as_completed(futures):
            conductor_records[futures[future]] = future.result()

    selected_records = [
        {
            "adapter_u": rational_string(candidate.parameter),
            "literal_shift_T": rational_string(2 * candidate.parameter),
            "profiles": [list(profile) for profile in candidate.profiles],
            "height": max(abs(candidate.parameter.numerator), candidate.parameter.denominator),
            "radical": candidate.radical,
            "conductor": conductor_records[candidate.parameter],
            "point_search": "not run; exact completed conductor below target required",
        }
        for candidate in selected
    ]
    seed_radical = radical_record(IMPORTED_SEED, trial_primes)
    feasible = [
        row for row in selected_records if row["conductor"].get("below_strict_log_conductor_target") is True
    ]
    seed_feasible = conductor_records[IMPORTED_SEED].get("below_strict_log_conductor_target") is True
    artifact: dict[str, Any] = {
        "schema_version": "elliptic-curves.fermigier-high-power-crt-gauss.v1",
        "status": "complete bounded exact population and conductor-first gate",
        "claim_level": "bounded computation; no rank claim",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "coordinate_normalization": {
            "search_coordinate": "canonical adapter u=s/2",
            "legacy_CAS_coordinate": "literal symmetric shift T=s=2u",
            "sign_quotient": "u and -u give the same canonical coefficients",
        },
        "population": {
            "primes": list(PRIMES),
            "exponent_profiles": [list(profile) for profile in EXPONENT_PROFILES],
            "profile_replay": profile_records,
            "raw_sign_quotiented_count": len(raw),
            "fresh_count": len(fresh),
            "fresh_parameter_sha256": rational_digest(fresh),
            "broad_neighbourhood_intersection_count": len(broad_intersection),
            "selection": f"first {CONDUCTOR_KEEP} by exact trial-radical upper integer",
            "trial_prime_bound": TRIAL_PRIME_BOUND,
        },
        "prior_exclusion": {
            **prior_record,
            "intersection_count": len(prior_intersection),
            "intersection": [rational_string(value) for value in prior_intersection],
            "imported_seed_role": "excluded from discovery; replayed only as calibration",
        },
        "selected": selected_records,
        "seed_calibration": {
            "adapter_u": rational_string(IMPORTED_SEED),
            "literal_shift_T": rational_string(2 * IMPORTED_SEED),
            "radical": seed_radical,
            "conductor": conductor_records[IMPORTED_SEED],
            "point_search": "not run unless the exact conductor gate completes below target",
        },
        "outcome": {
            "completed_conductors": sum(row["conductor"]["status"] == "completed" for row in selected_records),
            "conductor_timeouts_or_errors": sum(row["conductor"]["status"] != "completed" for row in selected_records),
            "strict_target_feasible_fresh_fibres": len(feasible),
            "seed_strict_target_feasible": seed_feasible,
            "point_search_calls": 0,
            "rank_calls": 0,
            "target_met": False,
            "interpretation": (
                "No completed exact conductor passed ln(N)<182.72; incomplete factorizations were not promoted."
                if not feasible and not seed_feasible
                else "At least one fibre requires a separately pinned point-search follow-up."
            ),
        },
        "parameters": {
            "conductor_keep": CONDUCTOR_KEEP,
            "conductor_timeout_seconds": conductor_timeout,
            "stack_bytes": stack_bytes,
            "workers": workers,
        },
        "software": {"python": platform.python_version()},
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves:elliptic-curves/cas .venv/bin/python "
            "elliptic-curves/cas/search_fermigier_high_power_crt_gauss.py"
        ),
    }
    artifact["result_sha256"] = result_digest(artifact)
    return artifact


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--conductor-timeout", type=float, default=25.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "archive/elliptic-curves/artifacts/generated-results/elliptic_fermigier_high_power_crt_gauss.json",
    )
    args = parser.parse_args()
    if not 1 <= args.conductor_timeout <= 60 or not 1 <= args.workers <= 4:
        raise SystemExit("bounded search requires timeout in [1,60] and workers in [1,4]")
    artifact = run(
        root=root,
        conductor_timeout=args.conductor_timeout,
        stack_bytes=args.stack_bytes,
        workers=args.workers,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps(artifact["outcome"], sort_keys=True))
    print(f"result_sha256={artifact['result_sha256']}")


if __name__ == "__main__":
    main()
