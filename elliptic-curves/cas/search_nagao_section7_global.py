#!/usr/bin/env python3
"""Leakage-free global parameter search in Nagao's 1994 section-7 family.

The search exhausts every positive primitive ``T=a/b`` in one declared
rectangle.  A small C++ helper retains a bounded frontier using *only* a
training prime band.  A disjoint validation band is carried forward without
affecting that retention, and a third disjoint band is evaluated exactly by
PARI only after the rectangle has closed.  The known ``5081/47`` fiber and
every parameter in the earlier neighborhood artifact are excluded before any
new frontier selection.

Exact conductor and root-number computations gate a three-stage quartic point
search.  All twelve visible sections, all six linear companions, and all
three quadratic companions from ``nagao_1994_section7.py`` are predeclared and
removed from point yield.  Stable numerical height rank is triage only.  Every
signal of rank at least 17 immediately triggers saturation and an exact
finite-reduction independence attempt.  A target hit is recorded only after
that exact certificate and the required conductor replay agree.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import gcd, log
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any, Iterable, Sequence

from ek_k3 import primes_up_to, rational_to_string
from nagao_1994 import (
    PRIMARY_SOURCE,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    short_jacobian_coefficients,
)
from nagao_1994_section7 import (
    SECTION7_CONSTRUCTION,
    SECTION7_CONSTRUCTOR_PARAMETER,
    SECTION7_LINEAR_COMPANION_SECTIONS,
    SECTION7_QUADRATIC_COMPANION_SECTIONS,
    SECTION7_ROOTS,
    section7_primitive_quartic_coefficients,
)
from pari_bridge import minimal_curve_data, pari_version
from search_extra_points import parse_point_vector, run_gp, signless_quartic_points
from search_nagao_rank20_t5081_neighborhood import (
    MAX_PROXY_SURVIVORS as NEIGHBORHOOD_MAX_PROXY_SURVIVORS,
    PROXY_LIMIT as NEIGHBORHOOD_PROXY_LIMIT,
    SAVING_PRIMES,
    build_beam_strata,
    build_residue_tables,
    build_trace_beams,
    conductor_radical_proxy,
    generate_candidates as generate_neighborhood_candidates,
    homogenized_discriminant,
    learn_discriminant_root_balls,
    learn_local_trace_fingerprint,
    projective_index,
)
from search_nagao_rank21_unbiased import finite_reduction_certificate
from triage_nagao_rank13_finalists import (
    height_matrix_replay,
    point_digest,
    point_on_short_curve,
    quartic_gp_polynomial,
    stable_height_rank,
)


Q = Fraction
TARGET_LOG_CONDUCTOR = Decimal("182.72")
DEFAULT_A_MAX = 30_000
DEFAULT_B_MAX = 1_000
DEFAULT_GLOBAL_KEEP = 40_000
DEFAULT_PER_DENOMINATOR_KEEP = 20
DEFAULT_VALIDATION_KEEP = 800
DEFAULT_CONDUCTOR_KEEP = 96
DEFAULT_POINT_KEEP = 32
TRAINING_PRIME_MIN = 5
TRAINING_PRIME_MAX = 199
VALIDATION_PRIME_MIN = 211
VALIDATION_PRIME_MAX = 397
FINAL_PRIME_MIN = 401
FINAL_PRIME_MAX = 1_999
SCORE_SCALE = 10**12
PROXY_TRIAL_BOUND = 251
PROXY_FEASIBILITY_MARGIN = 196.0
CHECKPOINT_RANK = 17
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_section7_global.py"
)


@dataclass(frozen=True)
class ScanCandidate:
    parameter: Fraction
    training_scaled: int
    validation_scaled: int

    @property
    def identifier(self) -> str:
        return f"section7-global-{self.parameter.numerator}-{self.parameter.denominator}"

    @property
    def height(self) -> int:
        return max(self.parameter.numerator, self.parameter.denominator)

    @property
    def training_score(self) -> float:
        return self.training_scaled / SCORE_SCALE

    @property
    def validation_score(self) -> float:
        return self.validation_scaled / SCORE_SCALE


@dataclass(frozen=True)
class ProxyCandidate:
    scanned: ScanCandidate
    proxy: dict[str, Any]

    @property
    def parameter(self) -> Fraction:
        return self.scanned.parameter

    @property
    def identifier(self) -> str:
        return self.scanned.identifier

    @property
    def height(self) -> int:
        return self.scanned.height

    @property
    def training_score(self) -> float:
        return self.scanned.training_score

    @property
    def validation_score(self) -> float:
        return self.scanned.validation_score


@dataclass(frozen=True)
class ExactScoreCandidate:
    proxied: ProxyCandidate
    score_b2000: str
    good_primes: int
    bad_primes: int
    last_prime: int

    @property
    def parameter(self) -> Fraction:
        return self.proxied.parameter

    @property
    def identifier(self) -> str:
        return self.proxied.identifier

    @property
    def height(self) -> int:
        return self.proxied.height

    @property
    def training_score(self) -> float:
        return self.proxied.training_score

    @property
    def validation_score(self) -> float:
        return self.proxied.validation_score

    @property
    def proxy(self) -> dict[str, Any]:
        return self.proxied.proxy


@dataclass(frozen=True)
class ConductorReplay:
    candidate: ExactScoreCandidate
    status: str
    data: dict[str, Any]
    error: str | None = None


@dataclass(frozen=True)
class PointPool:
    candidate: ExactScoreCandidate
    height_bound: int
    status: str
    signed_points: int
    signless_points: int
    predeclared_abscissas_returned: int
    new_images: tuple[tuple[Fraction, Fraction], ...]
    seed_images: tuple[tuple[Fraction, Fraction], ...]
    coefficients: tuple[Fraction, ...]
    wall_seconds: float
    pari_milliseconds: int
    error: str | None = None

    @property
    def pool(self) -> tuple[tuple[Fraction, Fraction], ...]:
        return self.seed_images + self.new_images


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prime_band(
    tables: dict[int, Sequence[Any]], lower: int, upper: int
) -> dict[int, Sequence[Any]]:
    """Return one closed prime band from precomputed residue tables."""

    if lower > upper:
        raise ValueError("the prime-band lower bound exceeds its upper bound")
    return {
        prime: table
        for prime, table in tables.items()
        if lower <= prime <= upper
    }


def scaled_score_table_input(
    training: dict[int, Sequence[Any]], validation: dict[int, Sequence[Any]]
) -> str:
    """Serialize two disjoint score bands for the exhaustive C++ helper."""

    if not training or not validation or set(training) & set(validation):
        raise ValueError("scanner score bands must be nonempty and disjoint")
    lines: list[str] = []
    for band in (training, validation):
        lines.append(str(len(band)))
        for prime, table in sorted(band.items()):
            if len(table) != prime + 1:
                raise ValueError("a projective residue table has the wrong size")
            weights = [round(float(symbol.contribution) * SCORE_SCALE) for symbol in table]
            lines.append(" ".join([str(prime), *(str(value) for value in weights)]))
    return "\n".join(lines) + "\n"


def parse_scanner_output(
    output: str, *, a_max: int, b_max: int
) -> tuple[tuple[ScanCandidate, ...], dict[str, int]]:
    """Parse and independently validate the bounded C++ frontier."""

    lines = [line for line in output.splitlines() if line]
    if not lines or not lines[0].startswith("SUMMARY\t"):
        raise AssertionError("the global scanner omitted its summary")
    summary_parts = lines[0].split("\t")
    if len(summary_parts) != 5:
        raise AssertionError("the global scanner summary changed shape")
    summary = {
        "primitive_pairs_enumerated": int(summary_parts[1]),
        "global_frontier_count_before_union": int(summary_parts[2]),
        "per_denominator_frontier_count_before_union": int(summary_parts[3]),
        "retained_union_count": int(summary_parts[4]),
    }
    by_parameter: dict[Fraction, ScanCandidate] = {}
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != 5 or fields[0] != "ROW":
            raise AssertionError("the global scanner emitted an invalid row")
        numerator, denominator, training, validation = map(int, fields[1:])
        if not (1 <= numerator <= a_max and 1 <= denominator <= b_max):
            raise AssertionError("a scanner row escaped the declared rectangle")
        if gcd(numerator, denominator) != 1:
            raise AssertionError("a nonprimitive scanner row escaped")
        parameter = Q(numerator, denominator)
        candidate = ScanCandidate(parameter, training, validation)
        if parameter in by_parameter:
            raise AssertionError("the scanner frontier contains a duplicate")
        by_parameter[parameter] = candidate
    if len(by_parameter) != summary["retained_union_count"]:
        raise AssertionError("the scanner row count disagrees with its summary")
    return tuple(by_parameter.values()), summary


def run_global_scanner(
    *,
    source: Path,
    a_max: int,
    b_max: int,
    global_keep: int,
    per_denominator_keep: int,
    training: dict[int, Sequence[Any]],
    validation: dict[int, Sequence[Any]],
    compile_timeout: float,
    scan_timeout: float,
) -> tuple[tuple[ScanCandidate, ...], dict[str, Any]]:
    """Compile, synchronously run, and join the exhaustive scanner."""

    compiler = shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")
    if compiler is None:
        raise FileNotFoundError("a C++17 compiler was not found")
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="section7-global-") as directory:
        executable = Path(directory) / "scan"
        compile_result = subprocess.run(
            [
                compiler,
                "-std=c++17",
                "-O3",
                "-DNDEBUG",
                str(source),
                "-o",
                str(executable),
            ],
            text=True,
            capture_output=True,
            timeout=compile_timeout,
        )
        if compile_result.returncode != 0:
            raise RuntimeError(f"global scanner compilation failed: {compile_result.stderr[:1000]}")
        scan_started = time.monotonic()
        result = subprocess.run(
            [
                str(executable),
                str(a_max),
                str(b_max),
                str(global_keep),
                str(per_denominator_keep),
            ],
            input=scaled_score_table_input(training, validation),
            text=True,
            capture_output=True,
            timeout=scan_timeout,
        )
        scan_seconds = time.monotonic() - scan_started
        if result.returncode != 0:
            raise RuntimeError(f"global scanner failed: {result.stderr[:1000]}")
        candidates, summary = parse_scanner_output(
            result.stdout, a_max=a_max, b_max=b_max
        )
    summary.update(
        {
            "compile_and_scan_wall_seconds": time.monotonic() - started,
            "scan_wall_seconds": scan_seconds,
            "scanner_source_sha256": sha256_file(source),
            "score_scale": SCORE_SCALE,
        }
    )
    return candidates, summary


def expected_primitive_pair_count(a_max: int, b_max: int) -> int:
    """Slow exact reference count, used only by tests and small smoke scans."""

    return sum(
        gcd(numerator, denominator) == 1
        for denominator in range(1, b_max + 1)
        for numerator in range(1, a_max + 1)
    )


def prior_neighborhood_population(
    all_tables: dict[int, Sequence[Any]], artifact_path: Path
) -> tuple[set[Fraction], dict[str, Any]]:
    """Exactly reproduce the completed neighborhood population for exclusion."""

    fingerprint = learn_local_trace_fingerprint(all_tables)
    root_balls = learn_discriminant_root_balls()
    full, prefixes, _ = build_trace_beams(fingerprint, all_tables, width=4_000)
    strata = build_beam_strata(full, prefixes, root_balls, root_width=500)
    candidates, audit = generate_neighborhood_candidates(
        strata,
        proxy_limit=NEIGHBORHOOD_PROXY_LIMIT,
        max_survivors=NEIGHBORHOOD_MAX_PROXY_SURVIVORS,
    )
    stored = json.loads(artifact_path.read_text(encoding="utf-8"))
    expected_digest = stored["population"]["survivor_stream_sha256"]
    if audit["survivor_stream_sha256"] != expected_digest:
        raise AssertionError("the reproduced neighborhood population digest changed")
    return {candidate.parameter for candidate in candidates}, {
        "artifact": str(artifact_path),
        "artifact_sha256": sha256_file(artifact_path),
        "reproduced_count": len(candidates),
        "reproduced_stream_sha256": expected_digest,
    }


def attach_proxies(
    candidates: Sequence[ScanCandidate], *, trial_prime_bound: int
) -> tuple[ProxyCandidate, ...]:
    """Evaluate the exact homogeneous discriminant before later selection."""

    return tuple(
        ProxyCandidate(
            candidate,
            conductor_radical_proxy(
                candidate.parameter, trial_prime_bound=trial_prime_bound
            ),
        )
        for candidate in candidates
    )


def parameter_stream_sha256(candidates: Iterable[ProxyCandidate]) -> str:
    digest = hashlib.sha256()
    for candidate in sorted(
        candidates,
        key=lambda item: (item.parameter.denominator, item.parameter.numerator),
    ):
        digest.update(
            (
                f"{candidate.parameter}|{candidate.scanned.training_scaled}|"
                f"{candidate.scanned.validation_scaled}|"
                f"{candidate.proxy['log_radical_upper_proxy']!r}\n"
            ).encode()
        )
    return digest.hexdigest()


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
    """Advance the training frontier using held-out validation and proxy quotas."""

    if keep < 10:
        raise ValueError("the validation frontier is too small for declared quotas")
    by_validation = sorted(
        candidates,
        key=lambda item: (
            -item.validation_score,
            item.proxy["log_radical_upper_proxy"],
            item.height,
            item.identifier,
        ),
    )
    feasible = [
        item
        for item in by_validation
        if item.proxy["log_radical_upper_proxy"] < PROXY_FEASIBILITY_MARGIN
    ]
    by_proxy = sorted(
        candidates,
        key=lambda item: (
            item.proxy["log_radical_upper_proxy"],
            -item.validation_score,
            item.identifier,
        ),
    )
    by_training = sorted(
        candidates,
        key=lambda item: (
            -item.training_score,
            item.proxy["log_radical_upper_proxy"],
            item.identifier,
        ),
    )
    selected: dict[Fraction, ProxyCandidate] = {}
    _add_quota(selected, feasible, max(1, keep * 5 // 16))
    _add_quota(selected, by_validation, max(1, keep * 5 // 16))
    _add_quota(selected, by_proxy, max(1, keep // 8))
    _add_quota(selected, by_training, max(1, keep // 16))

    # Preserve denominator diversity without allowing it to dominate the tail.
    bins = ((1, 100), (101, 250), (251, 500), (501, 750), (751, 1_000_000_000))
    per_bin = max(1, keep // 80)
    for lower, upper in bins:
        _add_quota(
            selected,
            (
                item
                for item in by_validation
                if lower <= item.parameter.denominator <= upper
            ),
            per_bin,
        )
    _add_quota(selected, by_validation, keep - len(selected))
    answer = tuple(selected.values())[:keep]
    if len(answer) != min(keep, len(candidates)):
        raise AssertionError("the validation quota union did not fill")
    return answer


def gp_rational(value: Fraction) -> str:
    return f"({value.numerator}/{value.denominator})"


def exact_final_band_scores(
    candidates: Sequence[ProxyCandidate],
    *,
    lower: int,
    upper: int,
    batch_size: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[ExactScoreCandidate, ...]:
    """Evaluate a third, disjoint prime band on exact minimal models."""

    if batch_size < 1 or lower > upper or timeout <= 0:
        raise ValueError("invalid exact-score bounds")
    expected_primes = tuple(
        prime for prime in primes_up_to(upper) if prime >= lower
    )
    if not expected_primes:
        raise ValueError("the exact score band contains no primes")
    records: dict[str, tuple[str, int, int]] = {}
    for start in range(0, len(candidates), batch_size):
        batch = candidates[start : start + batch_size]
        commands = ["default(realprecision,80);"]
        for index, candidate in enumerate(batch):
            coefficients = short_jacobian_coefficients(
                SECTION7_CONSTRUCTION, candidate.parameter
            )
            vector = ",".join(gp_rational(Q(value)) for value in coefficients)
            commands.extend(
                (
                    f"E=ellminimalmodel(ellinit([{vector}]));",
                    "S=0;USED=0;BAD=0;",
                    (
                        f"forprime(p={lower},{upper},"
                        "if(valuation(E.disc,p)>0,BAD++,"
                        "A=ellap(E,p);S+=(2-A)/(p+1-A)*log(p);USED++));"
                    ),
                    f'print("ROW|{index}|",S,"|",USED,"|",BAD);',
                )
            )
        commands.append("quit")
        output, _ = run_gp(
            "\n".join(commands) + "\n",
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        observed = 0
        for line in output.splitlines():
            if not line.startswith("ROW|"):
                continue
            _, index_text, score, good, bad = line.split("|")
            candidate = batch[int(index_text)]
            records[candidate.identifier] = score, int(good), int(bad)
            observed += 1
        if observed != len(batch):
            raise RuntimeError("PARI omitted an exact final-band score")
    answer = tuple(
        ExactScoreCandidate(
            candidate,
            records[candidate.identifier][0],
            records[candidate.identifier][1],
            records[candidate.identifier][2],
            expected_primes[-1],
        )
        for candidate in candidates
    )
    return tuple(
        sorted(
            answer,
            key=lambda item: (
                -Decimal(item.score_b2000),
                -item.validation_score,
                item.proxy["log_radical_upper_proxy"],
                item.identifier,
            ),
        )
    )


def _add_exact_quota(
    selected: dict[Fraction, ExactScoreCandidate],
    ordered: Iterable[ExactScoreCandidate],
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


def select_conductor_population(
    candidates: Sequence[ExactScoreCandidate], *, keep: int
) -> tuple[ExactScoreCandidate, ...]:
    by_exact = sorted(candidates, key=lambda item: (-Decimal(item.score_b2000), item.identifier))
    by_validation = sorted(candidates, key=lambda item: (-item.validation_score, item.identifier))
    by_proxy = sorted(
        candidates,
        key=lambda item: (item.proxy["log_radical_upper_proxy"], item.identifier),
    )
    selected: dict[Fraction, ExactScoreCandidate] = {}
    _add_exact_quota(selected, by_exact, max(1, keep * 5 // 12))
    _add_exact_quota(selected, by_proxy, max(1, keep // 3))
    _add_exact_quota(selected, by_validation, max(1, keep // 6))
    for lower, upper in ((1, 100), (101, 250), (251, 500), (501, 750), (751, 10**9)):
        _add_exact_quota(
            selected,
            (item for item in by_exact if lower <= item.parameter.denominator <= upper),
            1,
        )
    _add_exact_quota(selected, by_exact, keep - len(selected))
    return tuple(selected.values())[:keep]


def replay_one_conductor(
    candidate: ExactScoreCandidate, *, timeout: float, stack_bytes: int
) -> ConductorReplay:
    try:
        data = minimal_curve_data(
            short_jacobian_coefficients(SECTION7_CONSTRUCTION, candidate.parameter),
            timeout=timeout,
            local_primes=SAVING_PRIMES,
            stack_bytes=stack_bytes,
        )
        data["below_strict_log_conductor_target"] = (
            Decimal(data["log_conductor"]) < TARGET_LOG_CONDUCTOR
        )
        return ConductorReplay(candidate, "completed", data)
    except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
        return ConductorReplay(
            candidate,
            "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
            {},
            str(error)[:500],
        )


def parallel_conductor_replay(
    candidates: Sequence[ExactScoreCandidate],
    *,
    timeout: float,
    stack_bytes: int,
    workers: int,
) -> tuple[ConductorReplay, ...]:
    records: dict[str, ConductorReplay] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                replay_one_conductor,
                candidate,
                timeout=timeout,
                stack_bytes=stack_bytes,
            ): candidate.identifier
            for candidate in candidates
        }
        for future in as_completed(futures):
            records[futures[future]] = future.result()
    return tuple(records[candidate.identifier] for candidate in candidates)


def select_point_population(
    replays: Sequence[ConductorReplay], *, keep: int
) -> tuple[ExactScoreCandidate, ...]:
    completed = [replay for replay in replays if replay.status == "completed"]
    root_minus = sorted(
        (
            replay.candidate
            for replay in completed
            if replay.data.get("below_strict_log_conductor_target") is True
            and int(replay.data["root_number"]) == -1
        ),
        key=lambda item: (-Decimal(item.score_b2000), item.identifier),
    )
    root_plus = sorted(
        (
            replay.candidate
            for replay in completed
            if replay.data.get("below_strict_log_conductor_target") is True
            and int(replay.data["root_number"]) == 1
        ),
        key=lambda item: (-Decimal(item.score_b2000), item.identifier),
    )
    by_exact = sorted(
        (replay.candidate for replay in replays),
        key=lambda item: (-Decimal(item.score_b2000), item.identifier),
    )
    by_proxy = sorted(
        (replay.candidate for replay in replays),
        key=lambda item: (item.proxy["log_radical_upper_proxy"], item.identifier),
    )
    selected: dict[Fraction, ExactScoreCandidate] = {}
    _add_exact_quota(selected, root_minus, max(1, keep // 2))
    _add_exact_quota(selected, root_plus, max(1, keep // 4))
    # The alternative rank-30 target has no conductor restriction.
    _add_exact_quota(selected, by_exact, max(1, keep // 4))
    _add_exact_quota(selected, by_proxy, max(1, keep // 8))
    _add_exact_quota(selected, by_exact, keep - len(selected))
    return tuple(selected.values())[:keep]


def exact_predeclared_seeds(
    parameter: Fraction,
) -> tuple[
    tuple[tuple[Fraction, Fraction], ...],
    tuple[tuple[Fraction, Fraction], ...],
    tuple[Fraction, ...],
]:
    """Specialize and deduplicate all 12+6+3 known generic sections."""

    parameter = Q(parameter)
    quartic_points = (
        primitive_visible_points(SECTION7_CONSTRUCTION, parameter)
        + tuple(section.point(parameter) for section in SECTION7_LINEAR_COMPANION_SECTIONS)
        + tuple(section.point(parameter) for section in SECTION7_QUADRATIC_COMPANION_SECTIONS)
    )
    if len(quartic_points) != 21:
        raise AssertionError("the section-7 predeclared section count changed")
    coefficients = short_jacobian_coefficients(SECTION7_CONSTRUCTION, parameter)
    images = tuple(
        quartic_point_to_short_jacobian(SECTION7_CONSTRUCTION, parameter, point)
        for point in quartic_points
    )
    if any(not point_on_short_curve(coefficients, point) for point in images):
        raise AssertionError("a predeclared section missed the exact Jacobian")
    quartic_by_x: dict[Fraction, tuple[Fraction, Fraction]] = {}
    jacobian_by_x: dict[Fraction, tuple[Fraction, Fraction]] = {}
    for quartic_point, image in zip(quartic_points, images):
        quartic_by_x.setdefault(quartic_point[0], quartic_point)
        jacobian_by_x.setdefault(image[0], image)
    return tuple(quartic_by_x.values()), tuple(jacobian_by_x.values()), coefficients


def bounded_quartic_points(
    parameter: Fraction,
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], float, int]:
    quartic = section7_primitive_quartic_coefficients(parameter)
    program = "\n".join(
        (
            f"Q={quartic_gp_polynomial(quartic)};",
            "gettime();",
            f"R=hyperellratpoints(Q,{height_bound});",
            'print("PARI_MILLISECONDS ",gettime());',
            'print("POINTS ",R);',
            "quit",
        )
    ) + "\n"
    output, wall_seconds = run_gp(program, timeout=timeout, stack_bytes=stack_bytes)
    match = re.search(r"PARI_MILLISECONDS (\d+)", output)
    if match is None or "POINTS " not in output:
        raise AssertionError("PARI omitted bounded point output")
    return (
        parse_point_vector(output.split("POINTS ", 1)[1]),
        wall_seconds,
        int(match.group(1)),
    )


def map_new_points(
    parameter: Fraction,
    raw_points: Sequence[tuple[Fraction, Fraction]],
    seed_quartic: Sequence[tuple[Fraction, Fraction]],
    seed_images: Sequence[tuple[Fraction, Fraction]],
) -> tuple[tuple[tuple[Fraction, Fraction], ...], int]:
    coefficients = short_jacobian_coefficients(SECTION7_CONSTRUCTION, parameter)
    seed_quartic_x = {point[0] for point in seed_quartic}
    seen_image_x = {point[0] for point in seed_images}
    new_images = []
    returned = 0
    for quartic_point in signless_quartic_points(tuple(raw_points)):
        if quartic_point[1] == 0:
            continue
        if quartic_point[0] in seed_quartic_x:
            returned += 1
        image = quartic_point_to_short_jacobian(
            SECTION7_CONSTRUCTION, parameter, quartic_point
        )
        if not point_on_short_curve(coefficients, image):
            raise AssertionError("a searched point missed the exact Jacobian")
        if image[0] in seen_image_x:
            continue
        seen_image_x.add(image[0])
        new_images.append(image)
    return tuple(new_images), returned


def search_one_pool(
    candidate: ExactScoreCandidate,
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> PointPool:
    started = time.monotonic()
    try:
        seed_quartic, seed_images, coefficients = exact_predeclared_seeds(
            candidate.parameter
        )
        raw, wall_seconds, milliseconds = bounded_quartic_points(
            candidate.parameter,
            height_bound=height_bound,
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        new_images, returned = map_new_points(
            candidate.parameter, raw, seed_quartic, seed_images
        )
        return PointPool(
            candidate,
            height_bound,
            "completed",
            len(raw),
            len({point[0] for point in raw}),
            returned,
            new_images,
            seed_images,
            coefficients,
            wall_seconds,
            milliseconds,
        )
    except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
        return PointPool(
            candidate,
            height_bound,
            "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
            0,
            0,
            0,
            (),
            (),
            (),
            time.monotonic() - started,
            0,
            str(error)[:500],
        )


def parallel_point_search(
    candidates: Sequence[ExactScoreCandidate],
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
    workers: int,
) -> tuple[PointPool, ...]:
    records: dict[str, PointPool] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                search_one_pool,
                candidate,
                height_bound=height_bound,
                timeout=timeout,
                stack_bytes=stack_bytes,
            ): candidate.identifier
            for candidate in candidates
        }
        for future in as_completed(futures):
            records[futures[future]] = future.result()
    return tuple(records[candidate.identifier] for candidate in candidates)


def rank_one_pool(
    pool: PointPool,
    *,
    precisions: tuple[int, ...],
    timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    if pool.status != "completed":
        return {"status": "not-run", "reason": "point search did not complete"}
    try:
        runs = height_matrix_replay(
            pool.coefficients,
            pool.pool,
            precisions=precisions,
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        rank = stable_height_rank(runs)
        indices = tuple(runs[-1]["subset_indices_one_based"])
        selected = tuple(pool.pool[index - 1] for index in indices)
        return {
            "status": "completed",
            "stable_numerical_rank": rank,
            "precision_runs": list(runs),
            "selected_subset_indices_one_based": list(indices),
            "selected_points": selected,
            "selected_point_sha256": point_digest(selected),
        }
    except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
        return {
            "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
            "error": str(error)[:500],
        }


def parallel_rank_replay(
    pools: Sequence[PointPool],
    *,
    precisions: tuple[int, ...],
    timeout: float,
    stack_bytes: int,
    workers: int,
) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                rank_one_pool,
                pool,
                precisions=precisions,
                timeout=timeout,
                stack_bytes=stack_bytes,
            ): pool.candidate.identifier
            for pool in pools
        }
        for future in as_completed(futures):
            records[futures[future]] = future.result()
    return records


def point_priority(
    pool: PointPool, rank: dict[str, Any], conductor_by_id: dict[str, ConductorReplay]
) -> tuple[Any, ...]:
    replay = conductor_by_id.get(pool.candidate.identifier)
    root_minus_feasible = bool(
        replay is not None
        and replay.status == "completed"
        and replay.data.get("below_strict_log_conductor_target") is True
        and int(replay.data["root_number"]) == -1
    )
    return (
        rank.get("status") != "completed",
        -int(rank.get("stable_numerical_rank", -1)),
        not root_minus_feasible,
        -len(pool.new_images),
        -Decimal(pool.candidate.score_b2000),
        pool.candidate.identifier,
    )


def score_candidate_record(candidate: ExactScoreCandidate) -> dict[str, Any]:
    return {
        "candidate_id": candidate.identifier,
        "constructor_parameter_T": rational_to_string(candidate.parameter),
        "height": candidate.height,
        "training_score_5_199": candidate.training_score,
        "heldout_validation_score_211_397": candidate.validation_score,
        "exact_final_score_401_1999": candidate.score_b2000,
        "exact_final_good_primes": candidate.good_primes,
        "exact_final_bad_primes_skipped": candidate.bad_primes,
        "radical_proxy": candidate.proxy,
    }


def point_pool_record(pool: PointPool, rank: dict[str, Any]) -> dict[str, Any]:
    return {
        **score_candidate_record(pool.candidate),
        "point_search": {
            "height_bound": pool.height_bound,
            "status": pool.status,
            "signed_points": pool.signed_points,
            "distinct_quartic_abscissas": pool.signless_points,
            "predeclared_abscissas_returned": pool.predeclared_abscissas_returned,
            "new_distinct_jacobian_sign_pairs_beyond_21_predeclared": len(pool.new_images),
            "new_image_sha256": point_digest(pool.new_images),
            "wall_seconds": pool.wall_seconds,
            "pari_milliseconds": pool.pari_milliseconds,
            **({"error": pool.error} if pool.error else {}),
        },
        "height_rank": {
            key: value for key, value in rank.items() if key != "selected_points"
        },
    }


def parse_positive_ints(value: str) -> tuple[int, ...]:
    try:
        answer = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from error
    if not answer or any(item < 1 for item in answer):
        raise argparse.ArgumentTypeError("integers must be positive")
    return answer


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts/generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a-max", type=int, default=DEFAULT_A_MAX)
    parser.add_argument("--b-max", type=int, default=DEFAULT_B_MAX)
    parser.add_argument("--global-keep", type=int, default=DEFAULT_GLOBAL_KEEP)
    parser.add_argument(
        "--per-denominator-keep", type=int, default=DEFAULT_PER_DENOMINATOR_KEEP
    )
    parser.add_argument("--validation-keep", type=int, default=DEFAULT_VALIDATION_KEEP)
    parser.add_argument("--conductor-keep", type=int, default=DEFAULT_CONDUCTOR_KEEP)
    parser.add_argument("--point-keep", type=int, default=DEFAULT_POINT_KEEP)
    parser.add_argument("--compile-timeout", type=float, default=30.0)
    parser.add_argument("--scan-timeout", type=float, default=180.0)
    parser.add_argument("--exact-score-batch", type=int, default=80)
    parser.add_argument("--exact-score-timeout", type=float, default=120.0)
    parser.add_argument("--conductor-timeout", type=float, default=25.0)
    parser.add_argument("--conductor-workers", type=int, default=4)
    parser.add_argument(
        "--stage-heights", type=parse_positive_ints, default=(50_000, 250_000, 1_000_000)
    )
    parser.add_argument("--stage-keeps", type=parse_positive_ints, default=(10, 3))
    parser.add_argument("--stage-timeouts", type=parse_positive_ints, default=(6, 25, 100))
    parser.add_argument("--stage-workers", type=parse_positive_ints, default=(4, 4, 2))
    parser.add_argument("--height-precisions", type=parse_positive_ints, default=(72, 120))
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--saturation-timeout", type=float, default=40.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=1_000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--prior-neighborhood",
        type=Path,
        default=generated / "elliptic_nagao_rank20_t5081_neighborhood.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_nagao_section7_global.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    positive = (
        args.a_max,
        args.b_max,
        args.global_keep,
        args.per_denominator_keep,
        args.validation_keep,
        args.conductor_keep,
        args.point_keep,
        args.exact_score_batch,
        args.conductor_workers,
    )
    if min(positive) < 1:
        raise SystemExit("all population bounds must be positive")
    if not (
        args.validation_keep >= args.conductor_keep >= args.point_keep
    ):
        raise SystemExit("post-scan population caps must be nonincreasing")
    if len(args.stage_heights) != 3 or len(args.stage_keeps) != 2:
        raise SystemExit("this lane requires three point stages and two retention counts")
    if len(args.stage_timeouts) != 3 or len(args.stage_workers) != 3:
        raise SystemExit("provide one timeout and worker count per point stage")
    if tuple(sorted(set(args.stage_heights))) != args.stage_heights:
        raise SystemExit("point heights must be strictly increasing")
    if tuple(sorted(set(args.height_precisions))) != args.height_precisions:
        raise SystemExit("height precisions must be strictly increasing")
    if args.stage_keeps[1] > args.stage_keeps[0] or args.stage_keeps[0] > args.point_keep:
        raise SystemExit("point-stage populations must be nonincreasing")
    if min(
        args.compile_timeout,
        args.scan_timeout,
        args.exact_score_timeout,
        args.conductor_timeout,
        *args.stage_timeouts,
        args.height_timeout,
        args.saturation_timeout,
    ) <= 0:
        raise SystemExit("all subprocess timeouts must be positive")

    started = time.monotonic()
    score_tables = build_residue_tables(VALIDATION_PRIME_MAX)
    training_tables = prime_band(
        score_tables, TRAINING_PRIME_MIN, TRAINING_PRIME_MAX
    )
    validation_tables = prime_band(
        score_tables, VALIDATION_PRIME_MIN, VALIDATION_PRIME_MAX
    )
    if set(training_tables) & set(validation_tables):
        raise AssertionError("training and validation prime bands overlap")

    root = Path(__file__).resolve().parents[2]
    scanner_source = root / "elliptic-curves/cas/scan_nagao_section7_global.cpp"
    scanned, scan_audit = run_global_scanner(
        source=scanner_source,
        a_max=args.a_max,
        b_max=args.b_max,
        global_keep=args.global_keep,
        per_denominator_keep=args.per_denominator_keep,
        training=training_tables,
        validation=validation_tables,
        compile_timeout=args.compile_timeout,
        scan_timeout=args.scan_timeout,
    )
    print(
        f"section7-global primitive={scan_audit['primitive_pairs_enumerated']} "
        f"training_frontier={len(scanned)}",
        flush=True,
    )

    neighborhood_tables = prime_band(score_tables, 5, 200)
    prior, prior_audit = prior_neighborhood_population(
        neighborhood_tables, args.prior_neighborhood
    )
    calibration = SECTION7_CONSTRUCTOR_PARAMETER
    excluded_prior = sum(candidate.parameter in prior for candidate in scanned)
    excluded_calibration = sum(candidate.parameter == calibration for candidate in scanned)
    novel_scanned = tuple(
        candidate
        for candidate in scanned
        if candidate.parameter not in prior and candidate.parameter != calibration
    )
    if not novel_scanned:
        raise SystemExit("the prior-population exclusion emptied the global frontier")

    proxy_started = time.monotonic()
    proxied = attach_proxies(novel_scanned, trial_prime_bound=PROXY_TRIAL_BOUND)
    proxy_seconds = time.monotonic() - proxy_started
    validation_population = select_validation_population(
        proxied, keep=min(args.validation_keep, len(proxied))
    )
    print(
        f"section7-global novel={len(proxied)} validation={len(validation_population)} "
        f"proxy_seconds={proxy_seconds:.3f}",
        flush=True,
    )

    exact_scores = exact_final_band_scores(
        validation_population,
        lower=FINAL_PRIME_MIN,
        upper=FINAL_PRIME_MAX,
        batch_size=args.exact_score_batch,
        timeout=args.exact_score_timeout,
        stack_bytes=args.stack_bytes,
    )
    conductor_population = select_conductor_population(
        exact_scores, keep=min(args.conductor_keep, len(exact_scores))
    )
    conductor_replays = parallel_conductor_replay(
        conductor_population,
        timeout=args.conductor_timeout,
        stack_bytes=args.stack_bytes,
        workers=args.conductor_workers,
    )
    conductor_by_id = {
        replay.candidate.identifier: replay for replay in conductor_replays
    }
    completed_conductors = [
        replay for replay in conductor_replays if replay.status == "completed"
    ]
    subtarget = [
        replay
        for replay in completed_conductors
        if replay.data.get("below_strict_log_conductor_target") is True
    ]
    print(
        f"section7-global exact_scores={len(exact_scores)} "
        f"conductors={len(completed_conductors)}/{len(conductor_replays)} "
        f"subtarget={len(subtarget)} root_minus="
        f"{sum(int(replay.data['root_number']) == -1 for replay in subtarget)}",
        flush=True,
    )

    point_population = select_point_population(
        conductor_replays, keep=min(args.point_keep, len(conductor_replays))
    )
    if not point_population:
        raise SystemExit("the exact conductor stage produced no point population")

    stages: list[dict[str, Any]] = []
    checkpoints: dict[str, dict[str, Any]] = {}
    checkpoint_pools: dict[str, PointPool] = {}
    current = tuple(point_population)
    final_pools: tuple[PointPool, ...] = ()
    final_ranks: dict[str, dict[str, Any]] = {}
    for stage_index, (height, timeout, workers) in enumerate(
        zip(args.stage_heights, args.stage_timeouts, args.stage_workers), start=1
    ):
        pools = parallel_point_search(
            current,
            height_bound=height,
            timeout=float(timeout),
            stack_bytes=args.stack_bytes,
            workers=workers,
        )
        ranks = parallel_rank_replay(
            pools,
            precisions=args.height_precisions,
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
            workers=workers,
        )
        ordered = tuple(
            sorted(
                pools,
                key=lambda pool: point_priority(
                    pool, ranks[pool.candidate.identifier], conductor_by_id
                ),
            )
        )
        for pool in ordered:
            rank = ranks[pool.candidate.identifier]
            if (
                rank.get("status") != "completed"
                or int(rank["stable_numerical_rank"]) < CHECKPOINT_RANK
                or pool.candidate.identifier in checkpoints
            ):
                continue
            try:
                certificate = finite_reduction_certificate(
                    pool,
                    rank,
                    saturation_timeout=args.saturation_timeout,
                    certificate_prime_bound=args.certificate_prime_bound,
                    stack_bytes=args.stack_bytes,
                )
            except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
                certificate = {
                    "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
                    "error": str(error)[:500],
                }
            certificate["trigger_height"] = height
            certificate["trigger_stable_numerical_rank"] = rank["stable_numerical_rank"]
            certificate["exact_selected_points"] = [
                [rational_to_string(point[0]), rational_to_string(point[1])]
                for point in rank.get("selected_points", ())
            ]
            checkpoints[pool.candidate.identifier] = certificate
            checkpoint_pools[pool.candidate.identifier] = pool

        keep_count = (
            args.stage_keeps[stage_index - 1]
            if stage_index <= len(args.stage_keeps)
            else len(ordered)
        )
        retained = ordered[: min(keep_count, len(ordered))]
        stages.append(
            {
                "stage": stage_index,
                "quartic_naive_height_bound": height,
                "population_searched": len(pools),
                "completed_point_searches": sum(pool.status == "completed" for pool in pools),
                "point_search_timeouts": sum(pool.status == "timeout" for pool in pools),
                "ranked_population": [
                    point_pool_record(pool, ranks[pool.candidate.identifier])
                    for pool in ordered
                ],
                "retained_candidate_ids": [pool.candidate.identifier for pool in retained],
            }
        )
        current = tuple(pool.candidate for pool in retained)
        final_pools, final_ranks = ordered, ranks
        print(
            f"section7-global H={height} n={len(pools)} best_rank="
            f"{max((int(record.get('stable_numerical_rank', 0)) for record in ranks.values()), default=0)}",
            flush=True,
        )

    certified_hits = []
    for identifier, certificate in checkpoints.items():
        replay = conductor_by_id.get(identifier)
        rank = certificate.get("certified_algebraic_rank_lower_bound")
        if certificate.get("status") != "certified" or rank is None:
            continue
        if int(rank) >= 30 or (
            int(rank) >= 21
            and replay is not None
            and replay.status == "completed"
            and replay.data.get("below_strict_log_conductor_target") is True
        ):
            pool = checkpoint_pools[identifier]
            certified_hits.append(
                {
                    "candidate_id": identifier,
                    "constructor_parameter_T": rational_to_string(pool.candidate.parameter),
                    "certified_rank_lower_bound": int(rank),
                    "conductor": replay.data.get("conductor") if replay else None,
                    "log_conductor": replay.data.get("log_conductor") if replay else None,
                    "root_number": replay.data.get("root_number") if replay else None,
                }
            )

    # The calibration is measured only after all population decisions close.
    calibration_scanned = ScanCandidate(
        calibration,
        round(
            sum(
                training_tables[prime][
                    projective_index(calibration.numerator, calibration.denominator, prime)
                ].contribution
                for prime in training_tables
            )
            * SCORE_SCALE
        ),
        round(
            sum(
                validation_tables[prime][
                    projective_index(calibration.numerator, calibration.denominator, prime)
                ].contribution
                for prime in validation_tables
            )
            * SCORE_SCALE
        ),
    )
    calibration_proxy = ProxyCandidate(
        calibration_scanned,
        conductor_radical_proxy(calibration, trial_prime_bound=PROXY_TRIAL_BOUND),
    )
    calibration_exact = exact_final_band_scores(
        (calibration_proxy,),
        lower=FINAL_PRIME_MIN,
        upper=FINAL_PRIME_MAX,
        batch_size=1,
        timeout=args.exact_score_timeout,
        stack_bytes=args.stack_bytes,
    )[0]
    _, calibration_seed_images, calibration_coefficients = exact_predeclared_seeds(calibration)
    calibration_height_runs = height_matrix_replay(
        calibration_coefficients,
        calibration_seed_images,
        precisions=args.height_precisions,
        timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    if stable_height_rank(calibration_height_runs) != 12:
        raise AssertionError("the 21-section generic calibration baseline changed")

    script_path = Path(__file__).resolve()
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": (
            "bounded exhaustive global section-7 parameter search; exact certificates "
            "only, numerical height ranks are triage"
        ),
        "primary_source": PRIMARY_SOURCE,
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "certified_hits": certified_hits,
        },
        "family": {
            "roots_in_source_order": list(SECTION7_ROOTS),
            "constructor_parameter_convention": "T=2t_paper",
            "predeclared_generic_sections": {
                "visible": 12,
                "linear_companions": len(SECTION7_LINEAR_COMPANION_SECTIONS),
                "quadratic_companions": len(SECTION7_QUADRATIC_COMPANION_SECTIONS),
                "all_used_for_exact_decontamination": True,
                "generic_span_rank": 12,
            },
        },
        "population": {
            "definition": (
                "every positive primitive (a,b) with 1<=a<=A_MAX and "
                "1<=b<=B_MAX, before bounded training-only frontier retention"
            ),
            "a_max": args.a_max,
            "b_max": args.b_max,
            "all_pairs_rectangle_size": args.a_max * args.b_max,
            **scan_audit,
            "frontier_retention_uses_training_band_only": True,
            "global_keep": args.global_keep,
            "per_denominator_keep": args.per_denominator_keep,
            "prior_neighborhood_exclusion": prior_audit,
            "retained_rows_overlapping_prior_neighborhood": excluded_prior,
            "retained_rows_equal_to_calibration": excluded_calibration,
            "novel_training_frontier_count": len(proxied),
            "novel_training_frontier_sha256": parameter_stream_sha256(proxied),
            "proxy_trial_prime_bound": PROXY_TRIAL_BOUND,
            "proxy_wall_seconds": proxy_seconds,
        },
        "leakage_free_scoring": {
            "training_prime_band": [TRAINING_PRIME_MIN, TRAINING_PRIME_MAX],
            "training_primes": list(training_tables),
            "heldout_validation_prime_band": [VALIDATION_PRIME_MIN, VALIDATION_PRIME_MAX],
            "heldout_validation_primes": list(validation_tables),
            "exact_final_prime_band": [FINAL_PRIME_MIN, FINAL_PRIME_MAX],
            "bands_pairwise_disjoint": True,
            "validation_population_count": len(validation_population),
            "validation_selection_union": (
                "heldout score leaders, heldout score leaders under proxy 196, "
                "proxy leaders, training leaders, and denominator-bin leaders"
            ),
            "exact_final_ranked_population": [
                score_candidate_record(candidate) for candidate in exact_scores
            ],
        },
        "calibration_only_after_selection": {
            **score_candidate_record(calibration_exact),
            "excluded_before_every_population_selection": True,
            "known_exact_rank_lower_bound": 20,
            "known_exact_log_conductor": "174.249816228548038353904973690348230789531837166898175962879",
            "predeclared_21_section_stable_height_rank": 12,
            "height_precision_runs": list(calibration_height_runs),
        },
        "conductor_and_parity": {
            "population_selected": len(conductor_population),
            "completed": len(completed_conductors),
            "below_strict_target": len(subtarget),
            "below_target_root_minus": sum(
                int(replay.data["root_number"]) == -1 for replay in subtarget
            ),
            "records": [
                {
                    "candidate_id": replay.candidate.identifier,
                    "constructor_parameter_T": rational_to_string(replay.candidate.parameter),
                    "status": replay.status,
                    **replay.data,
                    **({"error": replay.error} if replay.error else {}),
                }
                for replay in conductor_replays
            ],
        },
        "point_triage": {
            "initial_population": [candidate.identifier for candidate in point_population],
            "selection_priority": (
                "subtarget root -1, exceptional subtarget root +1, unrestricted "
                "final-score leaders for rank 30, and proxy leaders"
            ),
            "all_21_generic_sections_decontaminated": True,
            "height_precisions": list(args.height_precisions),
            "stages": stages,
        },
        "exact_checkpoints_stable_numerical_rank_at_least_17": checkpoints,
        "final_frontier": [
            {
                "candidate_id": pool.candidate.identifier,
                "constructor_parameter_T": rational_to_string(pool.candidate.parameter),
                "stable_numerical_rank": final_ranks[pool.candidate.identifier].get(
                    "stable_numerical_rank"
                ),
                "new_points_beyond_21_predeclared": len(pool.new_images),
                "exact_final_score_401_1999": pool.candidate.score_b2000,
                "conductor": (
                    conductor_by_id[pool.candidate.identifier].data
                    if pool.candidate.identifier in conductor_by_id
                    else None
                ),
            }
            for pool in final_pools
        ],
        "bounds_and_caveats": {
            "rectangle_exhaustive_before_frontier_retention": True,
            "frontiers_and_point_searches_bounded": True,
            "conductor_proxy_is_not_a_conductor_bound": True,
            "numerical_height_rank_is_not_a_rank_certificate": True,
            "target_hit_requires_exact_finite_reduction_and_conductor_replay": True,
            "all_subprocesses_synchronous_joined_and_finitely_timed": True,
        },
        "reproducibility": {
            "command": REPRODUCING_COMMAND,
            "argv": [shlex.join(sys.argv)],
            "python": platform.python_version(),
            "pari": pari_version(),
            "script_sha256": sha256_file(script_path),
            "scanner_sha256": sha256_file(scanner_source),
            "total_wall_seconds": time.monotonic() - started,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        f"wrote {args.output} certified_hits={len(certified_hits)} "
        f"checkpoints={len(checkpoints)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
