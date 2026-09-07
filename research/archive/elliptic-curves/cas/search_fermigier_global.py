#!/usr/bin/env python3
"""Leakage-free global box scan in the even Fermigier--Mestre family.

This lane is deliberately independent of the earlier record-residue, CRT,
and projective-height populations.  It exhausts every primitive ``T=a/b``
with ``0 <= a <= A`` and ``1 <= b <= B``; the sign quotient is exact because
the family is even in ``T``.  Every rational recoverable from every prior
``elliptic_fermigier*.json`` artifact is removed before post-scan selection.

A C++ helper retains a union of discovery-only frontiers.  One feature is the
Fermigier good-reduction score; the second records the exact congruence
``p^2 | H(a/b)`` for the degree-20 discriminant factor and hence one unit of
local discriminant-to-radical saving.  A disjoint held prime band is computed
but cannot affect the scanner heaps.  Exact integer radical proxies are then
used alongside the held features to select a declared conductor tranche.

PARI conductor/root-number calls and the H=50000, 250000, 1000000 quartic
point searches all have strict foreground time caps and no retry.  Numerical
height rank is triage only.  Any stable numerical signal at least 21 triggers
an exact saturation and finite-reduction independence attempt immediately.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import gcd, log
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any, Iterable, Sequence

from ek_k3 import primes_up_to, rational_to_string
from exhaustive_multiple_root_height import homogeneous_discriminant_factor
from fermigier_mestre import (
    DISCRIMINANT_FACTOR_COEFFICIENTS,
    FermigierMestreFamily,
    NORMALIZED_RECORD_PARAMETER,
)
from pari_bridge import pari_version
from search_fermigier_published_pair_fiber_products import (
    PARAMETER_KEYS,
    PARAMETER_LIST_KEYS,
)
from search_fermigier_rank22_accidental_slices import (
    conductor_probe,
    finite_reduction_attempt,
    specialized_quartic_screen,
)


Q = Fraction
TARGET_LOG_CONDUCTOR = Decimal("182.72")
SCORE_SCALE = 10**9
COMPOSITE_POWER_DIVISOR = 16
DISCOVERY_PRIME_MIN = 401
DISCOVERY_PRIME_MAX = 499
HELD_PRIME_MIN = 503
HELD_PRIME_MAX = 599
DEFAULT_A_MAX = 100_000
DEFAULT_B_MAX = 1_000
DEFAULT_FRONTIER_KEEP = 8_000
DEFAULT_PER_DENOMINATOR_KEEP = 3
DEFAULT_VALIDATION_KEEP = 384
DEFAULT_CONDUCTOR_KEEP = 48
DEFAULT_POINT_KEEP = 12
DEFAULT_STAGE_HEIGHTS = (50_000, 250_000, 1_000_000)
DEFAULT_STAGE_KEEPS = (4, 2)
DEFAULT_STAGE_TIMEOUTS = (8, 20, 60)
CALIBRATIONS = (
    ("published-rank22-E22", NORMALIZED_RECORD_PARAMETER),
    ("known-H1m-numerical-rank16", Q(1666, 9)),
    ("known-H1m-numerical-rank15", Q(3115, 3)),
)


@dataclass(frozen=True)
class ModularTable:
    prime: int
    rank_weights: tuple[int, ...]
    power_weights: tuple[int, ...]

    @property
    def power_modulus(self) -> int:
        return self.prime * self.prime


@dataclass(frozen=True)
class ScanCandidate:
    parameter: Fraction
    discovery_rank_scaled: int
    discovery_power_scaled: int
    held_rank_scaled: int
    held_power_scaled: int

    @property
    def identifier(self) -> str:
        return f"fermigier-global-{self.parameter.numerator}-{self.parameter.denominator}"

    @property
    def height(self) -> int:
        return max(self.parameter.numerator, self.parameter.denominator)

    @property
    def discovery_rank_score(self) -> float:
        return self.discovery_rank_scaled / SCORE_SCALE

    @property
    def discovery_power_score(self) -> float:
        return self.discovery_power_scaled / SCORE_SCALE

    @property
    def held_rank_score(self) -> float:
        return self.held_rank_scaled / SCORE_SCALE

    @property
    def held_power_score(self) -> float:
        return self.held_power_scaled / SCORE_SCALE

    @property
    def discovery_composite_scaled(self) -> int:
        return self.discovery_rank_scaled + self.discovery_power_scaled // COMPOSITE_POWER_DIVISOR

    @property
    def held_composite_scaled(self) -> int:
        return self.held_rank_scaled + self.held_power_scaled // COMPOSITE_POWER_DIVISOR


@dataclass(frozen=True)
class ProxyCandidate:
    scanned: ScanCandidate
    radical_proxy: dict[str, Any]

    @property
    def parameter(self) -> Fraction:
        return self.scanned.parameter

    @property
    def identifier(self) -> str:
        return self.scanned.identifier


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial_mod(value: int, modulus: int) -> int:
    answer = 0
    for coefficient in reversed(DISCRIMINANT_FACTOR_COEFFICIENTS):
        answer = (answer * value + coefficient) % modulus
    return answer


def build_modular_tables(lower: int, upper: int) -> tuple[ModularTable, ...]:
    """Build exact local-rank and p^2 discriminant-saving tables."""

    from search_record_residue_class import build_score_tables

    if lower > upper or lower < 5:
        raise ValueError("invalid modular prime band")
    local_tables = {
        table.prime: table
        for table in build_score_tables(upper, "fermigier-good")
        if lower <= table.prime <= upper
    }
    answer = []
    for prime in sorted(local_tables):
        source = local_tables[prime]
        rank_weights = tuple(
            round((0.0 if cell is None else cell.term) * SCORE_SCALE)
            for cell in source.cells
        ) + (0,)
        modulus = prime * prime
        power = [0] * (modulus + 1)
        saving = round(log(prime) * SCORE_SCALE)
        for residue in range(prime):
            if polynomial_mod(residue, prime) != 0:
                continue
            for lift_digit in range(prime):
                lifted = residue + lift_digit * prime
                if polynomial_mod(lifted, modulus) == 0:
                    power[lifted] = saving
        answer.append(ModularTable(prime, rank_weights, tuple(power)))
    if not answer:
        raise ValueError("the modular prime band is empty")
    return tuple(answer)


def scanner_input(
    discovery: Sequence[ModularTable], held: Sequence[ModularTable]
) -> str:
    if not discovery or not held or {row.prime for row in discovery} & {row.prime for row in held}:
        raise ValueError("scanner bands must be nonempty and disjoint")
    lines: list[str] = []
    for band in (discovery, held):
        lines.append(str(len(band)))
        for table in band:
            if len(table.rank_weights) != table.prime + 1:
                raise AssertionError("a rank table changed shape")
            if len(table.power_weights) != table.power_modulus + 1:
                raise AssertionError("a power table changed shape")
            lines.append(" ".join([str(table.prime), *(str(value) for value in table.rank_weights)]))
            lines.append(" ".join(str(value) for value in table.power_weights))
    return "\n".join(lines) + "\n"


def expected_primitive_pair_count(a_max: int, b_max: int) -> int:
    """Exact reference count including the unique primitive zero pair (0,1)."""

    return 1 + sum(
        gcd(numerator, denominator) == 1
        for denominator in range(1, b_max + 1)
        for numerator in range(1, a_max + 1)
    )


def parse_scanner_output(
    output: str, *, a_max: int, b_max: int
) -> tuple[tuple[ScanCandidate, ...], dict[str, int]]:
    lines = [line for line in output.splitlines() if line]
    if not lines or not lines[0].startswith("SUMMARY\t"):
        raise AssertionError("the scanner omitted its summary")
    parts = lines[0].split("\t")
    if len(parts) != 7:
        raise AssertionError("the scanner summary changed shape")
    summary = {
        "primitive_pairs_enumerated": int(parts[1]),
        "rank_frontier_count_before_union": int(parts[2]),
        "power_frontier_count_before_union": int(parts[3]),
        "composite_frontier_count_before_union": int(parts[4]),
        "per_denominator_frontier_count_before_union": int(parts[5]),
        "retained_union_count": int(parts[6]),
    }
    records: dict[Fraction, ScanCandidate] = {}
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != 7 or fields[0] != "ROW":
            raise AssertionError("the scanner emitted an invalid row")
        a, b, rank, power, held_rank, held_power = map(int, fields[1:])
        if not (0 <= a <= a_max and 1 <= b <= b_max):
            raise AssertionError("a scanner row escaped the declared box")
        if gcd(a, b) != 1:
            raise AssertionError("a scanner row is not primitive")
        parameter = Q(a, b)
        if parameter in records:
            raise AssertionError("the scanner emitted a duplicate rational")
        records[parameter] = ScanCandidate(parameter, rank, power, held_rank, held_power)
    if len(records) != summary["retained_union_count"]:
        raise AssertionError("scanner row count disagrees with summary")
    return tuple(records.values()), summary


def run_scanner(
    *,
    source: Path,
    a_max: int,
    b_max: int,
    frontier_keep: int,
    per_denominator_keep: int,
    discovery: Sequence[ModularTable],
    held: Sequence[ModularTable],
    compile_timeout: float,
    scan_timeout: float,
) -> tuple[tuple[ScanCandidate, ...], dict[str, Any]]:
    compiler = shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")
    if compiler is None:
        raise FileNotFoundError("a C++17 compiler was not found")
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="fermigier-global-") as directory:
        executable = Path(directory) / "scan"
        compiled = subprocess.run(
            [compiler, "-std=c++17", "-O3", "-DNDEBUG", str(source), "-o", str(executable)],
            text=True,
            capture_output=True,
            timeout=compile_timeout,
        )
        if compiled.returncode != 0:
            raise RuntimeError(f"global scanner compilation failed: {compiled.stderr[:1000]}")
        scan_started = time.monotonic()
        result = subprocess.run(
            [
                str(executable),
                str(a_max),
                str(b_max),
                str(frontier_keep),
                str(frontier_keep),
                str(frontier_keep),
                str(per_denominator_keep),
            ],
            input=scanner_input(discovery, held),
            text=True,
            capture_output=True,
            timeout=scan_timeout,
        )
        scan_seconds = time.monotonic() - scan_started
        if result.returncode != 0:
            raise RuntimeError(f"global scanner failed: {result.stderr[:1000]}")
        candidates, summary = parse_scanner_output(result.stdout, a_max=a_max, b_max=b_max)
    summary.update(
        {
            "scan_wall_seconds": scan_seconds,
            "compile_and_scan_wall_seconds": time.monotonic() - started,
            "scanner_source_sha256": sha256_file(source),
        }
    )
    return candidates, summary


def feature_scores(
    parameter: Fraction, tables: Sequence[ModularTable]
) -> tuple[int, int]:
    parameter = abs(Q(parameter))
    numerator, denominator = parameter.numerator, parameter.denominator
    rank = 0
    power = 0
    for table in tables:
        prime = table.prime
        if denominator % prime == 0:
            rank_index = prime
            power_index = table.power_modulus
        else:
            rank_index = numerator * pow(denominator, -1, prime) % prime
            power_index = numerator * pow(denominator, -1, table.power_modulus) % table.power_modulus
        rank += table.rank_weights[rank_index]
        power += table.power_weights[power_index]
    return rank, power


def rational_stream_sha256(values: Iterable[Fraction]) -> str:
    digest = hashlib.sha256()
    for value in sorted(set(map(abs, values))):
        digest.update((rational_to_string(value) + "\n").encode())
    return digest.hexdigest()


def extract_parameter_values(value: Any, *, key: str | None = None) -> set[Fraction]:
    """Conservatively recover rational parameter fields from old artifacts."""

    answer: set[Fraction] = set()
    normalized_key = key.lower() if key else ""
    scalar_key = (
        key in PARAMETER_KEYS
        or normalized_key == "t"
        or "parameter" in normalized_key
    )
    list_key = key in PARAMETER_LIST_KEYS or "parameters" in normalized_key
    if scalar_key and isinstance(value, (str, int)):
        try:
            answer.add(abs(Q(value)))
        except (ValueError, ZeroDivisionError):
            pass
    if list_key and isinstance(value, list):
        for item in value:
            if isinstance(item, (str, int)):
                try:
                    answer.add(abs(Q(item)))
                except (ValueError, ZeroDivisionError):
                    pass
    if isinstance(value, dict):
        for child_key, child in value.items():
            answer.update(extract_parameter_values(child, key=child_key))
    elif isinstance(value, list):
        for child in value:
            answer.update(extract_parameter_values(child))
    return answer


def prior_fermigier_parameters(
    artifact_directory: Path, output_path: Path
) -> tuple[set[Fraction], dict[str, Any]]:
    parameters: set[Fraction] = set()
    sources: dict[str, Any] = {}
    output_resolved = output_path.resolve()
    for path in sorted(artifact_directory.glob("elliptic_fermigier*.json")):
        if path.resolve() == output_resolved:
            continue
        raw = path.read_bytes()
        extracted = extract_parameter_values(json.loads(raw))
        parameters.update(extracted)
        sources[path.name] = {
            "sha256": hashlib.sha256(raw).hexdigest(),
            "extracted_parameter_count": len(extracted),
        }
    parameters.update(abs(parameter) for _, parameter in CALIBRATIONS)
    return parameters, {
        "canonicalization": "T -> abs(T), exact because the family is even",
        "artifact_glob": "elliptic_fermigier*.json",
        "artifact_source_count": len(sources),
        "artifact_sources": sources,
        "unique_prior_or_calibration_parameter_count": len(parameters),
        "parameter_stream_sha256": rational_stream_sha256(parameters),
    }


def radical_proxy(parameter: Fraction, *, trial_prime_bound: int) -> dict[str, Any]:
    parameter = abs(Q(parameter))
    value = abs(homogeneous_discriminant_factor(parameter.numerator, parameter.denominator))
    if value == 0:
        return {"singular": True}
    remaining = value
    radical_log = 0.0
    repeated_log_saving = 0.0
    factors: dict[str, int] = {}
    for prime in primes_up_to(trial_prime_bound):
        if remaining % prime:
            continue
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        factors[str(prime)] = exponent
        radical_log += log(prime)
        repeated_log_saving += (exponent - 1) * log(prime)
    if remaining > 1:
        radical_log += log(remaining)
    return {
        "singular": False,
        "homogeneous_discriminant_factor_decimal_digits": len(str(value)),
        "trial_prime_bound": trial_prime_bound,
        "trial_factors": factors,
        "unfactored_cofactor_decimal_digits": len(str(remaining)),
        "log_abs_homogeneous_discriminant_factor": log(value),
        "log_radical_upper_proxy": radical_log,
        "exact_trial_repeated_prime_log_saving": repeated_log_saving,
        "upper_proxy_explanation": (
            "the unfactored cofactor is charged in full, so this is an upper "
            "bound for log(rad(H_hom)); it is not an exact conductor bound"
        ),
    }


def attach_proxies(
    candidates: Sequence[ScanCandidate], *, trial_prime_bound: int
) -> tuple[ProxyCandidate, ...]:
    return tuple(
        ProxyCandidate(candidate, radical_proxy(candidate.parameter, trial_prime_bound=trial_prime_bound))
        for candidate in candidates
        if candidate.parameter != 0
    )


def _add_quota(
    selected: dict[Fraction, ProxyCandidate],
    ordered: Iterable[ProxyCandidate],
    count: int,
) -> None:
    added = 0
    for candidate in ordered:
        if candidate.parameter in selected:
            continue
        selected[candidate.parameter] = candidate
        added += 1
        if added >= count:
            break


def select_validation_population(
    candidates: Sequence[ProxyCandidate], *, keep: int
) -> tuple[ProxyCandidate, ...]:
    """Use held features and exact H arithmetic only after discovery closes."""

    if keep < 16:
        raise ValueError("the held population is too small for declared quotas")
    by_rank = sorted(candidates, key=lambda row: (-row.scanned.held_rank_scaled, row.identifier))
    by_power = sorted(candidates, key=lambda row: (-row.scanned.held_power_scaled, row.identifier))
    by_composite = sorted(candidates, key=lambda row: (-row.scanned.held_composite_scaled, row.identifier))
    by_proxy = sorted(candidates, key=lambda row: (row.radical_proxy["log_radical_upper_proxy"], row.identifier))
    selected: dict[Fraction, ProxyCandidate] = {}
    _add_quota(selected, by_rank, keep * 5 // 16)
    _add_quota(selected, by_power, keep * 4 // 16)
    _add_quota(selected, by_composite, keep * 4 // 16)
    _add_quota(selected, by_proxy, keep * 2 // 16)
    for lower, upper in ((1, 50), (51, 200), (201, 500), (501, 1000)):
        _add_quota(
            selected,
            (row for row in by_composite if lower <= row.parameter.denominator <= upper),
            max(1, keep // 64),
        )
    _add_quota(selected, by_composite, keep - len(selected))
    return tuple(selected.values())[: min(keep, len(candidates))]


def select_conductor_population(
    candidates: Sequence[ProxyCandidate], *, keep: int
) -> tuple[ProxyCandidate, ...]:
    by_rank = sorted(candidates, key=lambda row: (-row.scanned.held_rank_scaled, row.identifier))
    by_power = sorted(candidates, key=lambda row: (-row.scanned.held_power_scaled, row.identifier))
    by_composite = sorted(candidates, key=lambda row: (-row.scanned.held_composite_scaled, row.identifier))
    by_proxy = sorted(candidates, key=lambda row: (row.radical_proxy["log_radical_upper_proxy"], row.identifier))
    selected: dict[Fraction, ProxyCandidate] = {}
    _add_quota(selected, by_rank, keep // 4)
    _add_quota(selected, by_power, keep // 4)
    _add_quota(selected, by_composite, keep // 4)
    _add_quota(selected, by_proxy, keep // 4)
    _add_quota(selected, by_composite, keep - len(selected))
    return tuple(selected.values())[: min(keep, len(candidates))]


def candidate_record(candidate: ProxyCandidate) -> dict[str, Any]:
    row = candidate.scanned
    return {
        "candidate_id": row.identifier,
        "parameter_t": rational_to_string(row.parameter),
        "height": row.height,
        "discovery_rank_score_401_499": row.discovery_rank_score,
        "discovery_p2_saving_score_401_499": row.discovery_power_score,
        "discovery_composite_scaled": row.discovery_composite_scaled,
        "held_rank_score_503_599": row.held_rank_score,
        "held_p2_saving_score_503_599": row.held_power_score,
        "held_composite_scaled": row.held_composite_scaled,
        "radical_proxy": candidate.radical_proxy,
    }


def conductor_one(
    candidate: ProxyCandidate, *, timeout: float, stack_bytes: int
) -> tuple[str, dict[str, Any]]:
    data = conductor_probe(candidate.parameter, timeout=timeout, stack_bytes=stack_bytes)
    data["retried"] = False
    data["strict_per_call_timeout_seconds"] = timeout
    return candidate.identifier, data


def run_conductors(
    candidates: Sequence[ProxyCandidate],
    *,
    timeout: float,
    stack_bytes: int,
    workers: int,
) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                conductor_one, candidate, timeout=timeout, stack_bytes=stack_bytes
            ): candidate.identifier
            for candidate in candidates
        }
        for future in as_completed(futures):
            identifier, data = future.result()
            records[identifier] = data
    return records


def select_point_population(
    candidates: Sequence[ProxyCandidate],
    conductors: dict[str, dict[str, Any]],
    *,
    keep: int,
) -> tuple[ProxyCandidate, ...]:
    qualifying = [
        row
        for row in candidates
        if conductors[row.identifier].get("below_strict_log_conductor_target") is True
    ]
    qualifying.sort(
        key=lambda row: (
            int(conductors[row.identifier].get("root_number", 1)) != -1,
            -row.scanned.held_composite_scaled,
            Decimal(conductors[row.identifier]["log_conductor"]),
            row.identifier,
        )
    )
    unrestricted = sorted(
        candidates,
        key=lambda row: (-row.scanned.held_composite_scaled, row.identifier),
    )
    selected: dict[Fraction, ProxyCandidate] = {}
    _add_quota(selected, qualifying, keep * 3 // 4)
    _add_quota(selected, unrestricted, keep // 4)
    _add_quota(selected, qualifying, keep - len(selected))
    _add_quota(selected, unrestricted, keep - len(selected))
    return tuple(selected.values())[: min(keep, len(candidates))]


def point_one(
    candidate: ProxyCandidate,
    *,
    height: int,
    search_timeout: float,
    height_timeout: float,
    stack_bytes: int,
) -> tuple[str, dict[str, Any], tuple[tuple[Fraction, Fraction], ...]]:
    started = time.monotonic()
    try:
        screen, selected = specialized_quartic_screen(
            candidate.parameter,
            height_bound=height,
            search_timeout=search_timeout,
            height_timeout=height_timeout,
            precisions=(72, 120),
            stack_bytes=stack_bytes,
        )
        screen["total_wrapper_wall_seconds"] = time.monotonic() - started
        screen["retried"] = False
        return candidate.identifier, screen, selected
    except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
        return candidate.identifier, {
            "parameter_t": rational_to_string(candidate.parameter),
            "search": {
                "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
                "error": str(error)[:1000],
                "height_bound": height,
                "timeout_seconds": search_timeout,
                "retried": False,
            },
            "height_rank": None,
            "total_wrapper_wall_seconds": time.monotonic() - started,
            "retried": False,
        }, ()


def numerical_rank(screen: dict[str, Any]) -> int:
    rank = screen.get("height_rank")
    return int(rank["stable_numerical_rank"]) if isinstance(rank, dict) else -1


def point_priority(row: dict[str, Any]) -> tuple[Any, ...]:
    screen = row["screen"]
    return (
        -numerical_rank(screen),
        -int(screen.get("new_distinct_direct_images_modulo_sign", 0)),
        -row["candidate"].scanned.held_composite_scaled,
        row["candidate"].identifier,
    )


def run_point_stage(
    candidates: Sequence[ProxyCandidate],
    *,
    height: int,
    search_timeout: float,
    height_timeout: float,
    stack_bytes: int,
    workers: int,
    certificates: dict[str, dict[str, Any]],
    saturation_timeout: float,
    certificate_prime_bound: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                point_one,
                candidate,
                height=height,
                search_timeout=search_timeout,
                height_timeout=height_timeout,
                stack_bytes=stack_bytes,
            ): candidate
            for candidate in candidates
        }
        for future in as_completed(futures):
            candidate = futures[future]
            identifier, screen, selected = future.result()
            rank = numerical_rank(screen)
            if rank >= 21 and (
                identifier not in certificates
                or int(certificates[identifier].get("trigger_numerical_rank", -1)) < rank
            ):
                try:
                    certificate = finite_reduction_attempt(
                        FermigierMestreFamily.coefficients(candidate.parameter),
                        selected,
                        saturation_timeout=saturation_timeout,
                        stack_bytes=stack_bytes,
                        certificate_prime_bound=certificate_prime_bound,
                    )
                except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
                    certificate = {
                        "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
                        "error": str(error)[:1000],
                    }
                certificate["trigger_numerical_rank"] = rank
                certificate["trigger_height_bound"] = height
                certificate["retried"] = False
                certificates[identifier] = certificate
            rows.append({"candidate": candidate, "screen": screen})
    rows.sort(key=point_priority)
    return rows


def serializable_point_row(row: dict[str, Any]) -> dict[str, Any]:
    return {**candidate_record(row["candidate"]), "screen": row["screen"]}


def calibration_records(
    discovery: Sequence[ModularTable], held: Sequence[ModularTable], trial_bound: int
) -> list[dict[str, Any]]:
    rows = []
    for label, parameter in CALIBRATIONS:
        discovery_rank, discovery_power = feature_scores(parameter, discovery)
        held_rank, held_power = feature_scores(parameter, held)
        rows.append(
            {
                "label": label,
                "parameter_t": rational_to_string(parameter),
                "role": "separate calibration only; excluded before held selection",
                "discovery_rank_score": discovery_rank / SCORE_SCALE,
                "discovery_p2_saving_score": discovery_power / SCORE_SCALE,
                "held_rank_score": held_rank / SCORE_SCALE,
                "held_p2_saving_score": held_power / SCORE_SCALE,
                "radical_proxy": radical_proxy(parameter, trial_prime_bound=trial_bound),
            }
        )
    return rows


def parse_positive_ints(value: str) -> tuple[int, ...]:
    try:
        answer = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated positive integers") from error
    if not answer or any(item < 1 for item in answer):
        raise argparse.ArgumentTypeError("expected comma-separated positive integers")
    return answer


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a-max", type=int, default=DEFAULT_A_MAX)
    parser.add_argument("--b-max", type=int, default=DEFAULT_B_MAX)
    parser.add_argument("--frontier-keep", type=int, default=DEFAULT_FRONTIER_KEEP)
    parser.add_argument("--per-denominator-keep", type=int, default=DEFAULT_PER_DENOMINATOR_KEEP)
    parser.add_argument("--validation-keep", type=int, default=DEFAULT_VALIDATION_KEEP)
    parser.add_argument("--conductor-keep", type=int, default=DEFAULT_CONDUCTOR_KEEP)
    parser.add_argument("--point-keep", type=int, default=DEFAULT_POINT_KEEP)
    parser.add_argument("--trial-prime-bound", type=int, default=251)
    parser.add_argument("--compile-timeout", type=float, default=30)
    parser.add_argument("--scan-timeout", type=float, default=240)
    parser.add_argument("--conductor-timeout", type=float, default=20)
    parser.add_argument("--conductor-workers", type=int, default=4)
    parser.add_argument("--stage-heights", type=parse_positive_ints, default=DEFAULT_STAGE_HEIGHTS)
    parser.add_argument("--stage-keeps", type=parse_positive_ints, default=DEFAULT_STAGE_KEEPS)
    parser.add_argument("--stage-timeouts", type=parse_positive_ints, default=DEFAULT_STAGE_TIMEOUTS)
    parser.add_argument("--point-workers", type=int, default=2)
    parser.add_argument("--height-timeout", type=float, default=30)
    parser.add_argument("--saturation-timeout", type=float, default=40)
    parser.add_argument("--certificate-prime-bound", type=int, default=1000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "artifacts/generated-results/elliptic_fermigier_global.json",
    )
    return parser


def validate_arguments(args: argparse.Namespace) -> None:
    population = (
        args.a_max,
        args.b_max,
        args.frontier_keep,
        args.per_denominator_keep,
        args.validation_keep,
        args.conductor_keep,
        args.point_keep,
        args.trial_prime_bound,
        args.conductor_workers,
        args.point_workers,
    )
    if min(population) < 1:
        raise SystemExit("all population bounds must be positive")
    if not args.validation_keep >= args.conductor_keep >= args.point_keep:
        raise SystemExit("post-scan population caps must be nonincreasing")
    if len(args.stage_heights) != 3 or len(args.stage_keeps) != 2 or len(args.stage_timeouts) != 3:
        raise SystemExit("this lane requires three point stages and two retention caps")
    if tuple(sorted(set(args.stage_heights))) != args.stage_heights:
        raise SystemExit("point heights must be strictly increasing")
    if args.stage_keeps[1] > args.stage_keeps[0] or args.stage_keeps[0] > args.point_keep:
        raise SystemExit("point-stage caps must be nonincreasing")
    if min(
        args.compile_timeout,
        args.scan_timeout,
        args.conductor_timeout,
        args.height_timeout,
        args.saturation_timeout,
        *args.stage_timeouts,
    ) <= 0:
        raise SystemExit("all process time caps must be positive")
    if max(args.stage_timeouts) > 60:
        raise SystemExit("point search timeouts may not exceed the process-group cap of 60 seconds")


def main() -> None:
    args = build_parser().parse_args()
    validate_arguments(args)
    started = time.monotonic()
    root = Path(__file__).resolve().parents[2]
    output_path = args.output if args.output.is_absolute() else root / args.output
    scanner_source = root / "elliptic-curves/cas/scan_fermigier_global.cpp"

    discovery = build_modular_tables(DISCOVERY_PRIME_MIN, DISCOVERY_PRIME_MAX)
    held = build_modular_tables(HELD_PRIME_MIN, HELD_PRIME_MAX)
    scanned, scan_audit = run_scanner(
        source=scanner_source,
        a_max=args.a_max,
        b_max=args.b_max,
        frontier_keep=args.frontier_keep,
        per_denominator_keep=args.per_denominator_keep,
        discovery=discovery,
        held=held,
        compile_timeout=args.compile_timeout,
        scan_timeout=args.scan_timeout,
    )
    print(
        f"fermigier-global primitive={scan_audit['primitive_pairs_enumerated']} "
        f"discovery_union={len(scanned)}",
        flush=True,
    )

    prior, prior_audit = prior_fermigier_parameters(output_path.parent, output_path)
    prior_in_box = {
        parameter
        for parameter in prior
        if parameter.numerator <= args.a_max and parameter.denominator <= args.b_max
    }
    excluded_retained = {candidate.parameter for candidate in scanned if candidate.parameter in prior}
    novel_scanned = tuple(candidate for candidate in scanned if candidate.parameter not in prior)
    proxied = attach_proxies(novel_scanned, trial_prime_bound=args.trial_prime_bound)
    nonsingular = tuple(row for row in proxied if not row.radical_proxy["singular"])
    validation = select_validation_population(nonsingular, keep=args.validation_keep)
    conductor_population = select_conductor_population(validation, keep=args.conductor_keep)
    print(
        f"novel={len(nonsingular)} held={len(validation)} conductors={len(conductor_population)}",
        flush=True,
    )

    conductors = run_conductors(
        conductor_population,
        timeout=args.conductor_timeout,
        stack_bytes=args.stack_bytes,
        workers=args.conductor_workers,
    )
    point_population = select_point_population(
        conductor_population, conductors, keep=args.point_keep
    )
    certificates: dict[str, dict[str, Any]] = {}
    stages: list[dict[str, Any]] = []
    current = point_population
    for stage_index, (height, timeout) in enumerate(zip(args.stage_heights, args.stage_timeouts, strict=True)):
        stage_started = time.monotonic()
        rows = run_point_stage(
            current,
            height=height,
            search_timeout=timeout,
            height_timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
            workers=args.point_workers,
            certificates=certificates,
            saturation_timeout=args.saturation_timeout,
            certificate_prime_bound=args.certificate_prime_bound,
        )
        retain = args.stage_keeps[stage_index] if stage_index < len(args.stage_keeps) else 0
        stages.append(
            {
                "height_bound": height,
                "input_count": len(current),
                "strict_search_timeout_seconds_per_candidate": timeout,
                "height_replay_timeout_seconds_per_candidate": args.height_timeout,
                "retried_calls": 0,
                "elapsed_seconds": time.monotonic() - stage_started,
                "results": [serializable_point_row(row) for row in rows],
                "retained_for_next_stage": [row["candidate"].identifier for row in rows[:retain]],
            }
        )
        if stage_index < len(args.stage_keeps):
            current = tuple(row["candidate"] for row in rows[:retain])
        print(
            f"point-stage H={height} input={len(rows)} "
            f"best-rank={max((numerical_rank(row['screen']) for row in rows), default=-1)}",
            flush=True,
        )

    certified_hits = []
    candidate_by_id = {row.identifier: row for row in conductor_population}
    for identifier, certificate in certificates.items():
        exact_rank = certificate.get("certified_algebraic_rank_lower_bound")
        if exact_rank is None:
            continue
        conductor = conductors.get(identifier, {})
        qualifies_small = (
            int(exact_rank) >= 21
            and conductor.get("below_strict_log_conductor_target") is True
        )
        qualifies_rank30 = int(exact_rank) >= 30
        if qualifies_small or qualifies_rank30:
            certified_hits.append(
                {
                    "candidate_id": identifier,
                    "parameter_t": rational_to_string(candidate_by_id[identifier].parameter),
                    "certified_rank_lower_bound": exact_rank,
                    "conductor": conductor,
                    "route": "rank21-small-conductor" if qualifies_small else "rank30",
                }
            )

    completed_conductors = [data for data in conductors.values() if data.get("status") == "completed"]
    below_target = [data for data in completed_conductors if data.get("below_strict_log_conductor_target") is True]
    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded exhaustive global-box experiment; enumeration, exclusions, modular "
            "features, and point membership are exact; conductor/root numbers are PARI "
            "computations; height ranks are numerical triage only"
        ),
        "family": "normalized Fermigier--Mestre family, exactly even in T",
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "certified_hits": certified_hits,
        },
        "population": {
            "definition": "every primitive 0<=a<=A, 1<=b<=B; T=a/b; T and -T identified",
            "a_max": args.a_max,
            "b_max": args.b_max,
            **scan_audit,
            "expected_primitive_pair_count_reference": expected_primitive_pair_count(args.a_max, args.b_max),
            "zero_parameter_role": "enumerated exactly, then excluded as a degenerate generic-section parameter",
            "prior_or_calibration_parameters_in_full_box": len(prior_in_box),
            "prior_or_calibration_parameters_in_discovery_union": len(excluded_retained),
            "novel_nonsingular_discovery_union": len(nonsingular),
        },
        "prior_exclusion": prior_audit,
        "leakage_free_features": {
            "discovery_prime_band": [DISCOVERY_PRIME_MIN, DISCOVERY_PRIME_MAX],
            "discovery_primes": [row.prime for row in discovery],
            "held_prime_band": [HELD_PRIME_MIN, HELD_PRIME_MAX],
            "held_primes": [row.prime for row in held],
            "bands_disjoint": not ({row.prime for row in discovery} & {row.prime for row in held}),
            "rank_formula": "sum good ((2-a_p)/(p+1-a_p))*log(p); bad and denominator primes contribute zero",
            "power_formula": "sum log(p) for exact affine congruences p^2 | H(a/b); denominator primes contribute zero",
            "discovery_frontiers": ["rank", "p2-discriminant-saving", "rank + saving/16", "per-denominator composite"],
            "held_information_used_by_scanner_retention": False,
            "point_results_used_before_point-stage_escalation": False,
            "score_scale": SCORE_SCALE,
            "composite_power_divisor": COMPOSITE_POWER_DIVISOR,
        },
        "calibrations": calibration_records(discovery, held, args.trial_prime_bound),
        "selection": {
            "discovery_union_after_exclusion": len(nonsingular),
            "held_validation_count": len(validation),
            "held_validation_candidates": [candidate_record(row) for row in validation],
            "conductor_tranche_count": len(conductor_population),
            "conductor_tranche": [candidate_record(row) for row in conductor_population],
            "point_tranche_count": len(point_population),
            "point_tranche_ids": [row.identifier for row in point_population],
        },
        "conductor_replays": {
            identifier: {**candidate_record(candidate_by_id[identifier]), "conductor": conductors[identifier]}
            for identifier in sorted(conductors)
        },
        "conductor_summary": {
            "requested": len(conductor_population),
            "completed": len(completed_conductors),
            "below_strict_target": len(below_target),
            "timeouts": sum(data.get("status") == "timeout" for data in conductors.values()),
            "errors": sum(data.get("status") == "error" for data in conductors.values()),
            "best_log_conductor": (
                str(min(Decimal(data["log_conductor"]) for data in completed_conductors))
                if completed_conductors
                else None
            ),
        },
        "point_triage": {
            "stages": stages,
            "exact_certification_trigger": "stable numerical height rank >=21 at any stage",
            "exact_certificates": certificates,
            "all_calls_no_retry": True,
            "point_processes_started_in_new_process_groups": True,
        },
        "parameters": {
            "frontier_keep_per_global_feature": args.frontier_keep,
            "per_denominator_keep": args.per_denominator_keep,
            "validation_keep": args.validation_keep,
            "conductor_keep": args.conductor_keep,
            "point_keep": args.point_keep,
            "trial_prime_bound": args.trial_prime_bound,
            "stage_heights": list(args.stage_heights),
            "stage_keeps": list(args.stage_keeps),
            "stage_timeouts": list(args.stage_timeouts),
            "conductor_timeout_seconds": args.conductor_timeout,
            "height_timeout_seconds": args.height_timeout,
            "saturation_timeout_seconds": args.saturation_timeout,
            "stack_bytes": args.stack_bytes,
        },
        "frontier_summary": {
            "best_stable_numerical_rank": max(
                (
                    numerical_rank(row["screen"])
                    for stage in stages
                    for row in ({"screen": item["screen"]} for item in stage["results"])
                ),
                default=-1,
            ),
            "best_point_stage_rows": [
                stage["results"][0] for stage in stages if stage["results"]
            ],
            "scope_warning": "bounded searches and numerical height ranks do not upper-bound algebraic rank",
        },
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "reproducibility": {
            "reproducing_command": command,
            "script_sha256": sha256_file(script_path),
            "scanner_sha256": sha256_file(scanner_source),
            "elapsed_seconds": time.monotonic() - started,
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        },
    }
    if artifact["population"]["primitive_pairs_enumerated"] != artifact["population"]["expected_primitive_pair_count_reference"]:
        raise AssertionError("the exhaustive C++ population count failed its independent reference")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"wrote {output_path} completed_conductors={len(completed_conductors)} "
        f"below_target={len(below_target)} certified_hits={len(certified_hits)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
