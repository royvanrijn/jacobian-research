#!/usr/bin/env python3
"""Exact max-root-300 Mestre census and bounded rank-aware fiber tranche.

Frozen max-root-200 code and artifacts are read-only inputs.  A new range
enumerator first reproduces the complete old obstruction stream record for
record, then exhausts diameters 201--300 in four balanced, disjoint ranges.
Every emitted record is replayed with an independent Python obstruction and
normalization check.  Generic nonsingularity is classified exactly with 21
specializations of the degree-at-most-20 discriminant.

The screen stage evaluates T=1,...,8 for every genuinely new nonsingular
nonreflection family.  All admissible visible point sets receive exact mod-3
finite-reduction signatures.  A fixed 64-family rank-aware/diversity tranche
then receives conductor-first computation and one H=5000 search.  Stable rank
at least 15 triggers an immediate exact independence certificate attempt.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import gcd
import os
from pathlib import Path
import platform
import shlex
import shutil
import sys
import tempfile
import time
from typing import Any, Iterable, Sequence

from mestre_root_tuples import normalize_integer_root_tuple
from search_mestre_root_tuple_scale import (
    CappedProcessTimeout,
    TARGET_LOG_CONDUCTOR,
    capped_minimal_curve_data,
    finite_reduction_attempt,
    point_record,
    run_capped_process,
    sha256_file,
    tuple_digest,
)
from search_mestre_root_tuple_scale_max100 import search_h5000, stable_json_digest
from search_mestre_root_tuple_scale_max200 import (
    fiber_rank_key,
    generic_classification,
    independent_normalized_count,
    mestre_obstruction_integer,
    mod3_independence_certificate,
    reflection_symmetric,
    screen_family,
    visible_points_and_coefficients,
)


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

MAX_ROOT = 300
PRIOR_MAX_ROOT = 200
EXTENSION_RANGES = ((201, 243), (244, 268), (269, 286), (287, 300))
PARAMETERS = tuple(range(1, 9))
SELECTED_FAMILY_COUNT = 64
GLOBAL_RANK_AWARE_FAMILY_KEEP = 34
DIVERSITY_KEEP_PER_DECILE = 3
STRONG_GAIN_TRIGGER = 15
EXACT_GAIN_PRIME_BOUND = 499
STACK_BYTES = 256_000_000
SCREEN_FAMILY_CAP = 5_000

EXPECTED_MAX200_CENSUS_SHA256 = (
    "7270769007f9c130fce8b1813164373de9c6a5eb1c6d86cfe71b8c96fada161b"
)
EXPECTED_MAX200_CENSUS_RESULT_SHA256 = (
    "352afa99131ba9af68a7e29309387cfb06b38e10762cef7e1b3866b6c1b38a90"
)
EXPECTED_MAX200_CPP_SHA256 = (
    "56f5111765315fefea45628066a2971894bb963469d48823f08b717ba91c0c3a"
)
EXPECTED_MAX200_DRIVER_SHA256 = (
    "405a2b9f7653c89af0e3e6caf2e77765cb4bfc88fccf88edffa67d3435aebf24"
)
EXPECTED_MAX200_SCREEN_SHA256 = (
    "5e1b53e187520735efba46fc8fd9cbdd4dfd4284545a815f6416baf3be84f342"
)


@dataclass(frozen=True)
class RangeResult:
    first: int
    last: int
    normalized_count: int
    obstruction_count: int
    reflection_count: int
    nonreflection_count: int
    nonreflection_roots: tuple[tuple[int, ...], ...]
    stdout_sha256: str
    stdout: str


def parse_range_output(stdout: str, first: int, last: int) -> RangeResult:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not lines or lines[0] != "MESTRE_ROOT_TUPLES_V3_RANGE":
        raise AssertionError("the max300 range enumerator omitted its V3 header")
    summary = lines[-1].split()
    if summary[0] != "S" or len(summary) != 7:
        raise AssertionError("the max300 range summary changed")
    if (int(summary[1]), int(summary[2])) != (first, last):
        raise AssertionError("the max300 range bounds changed")
    nonreflection = []
    reflection_count = 0
    for line in lines[1:-1]:
        fields = line.split()
        if fields[0] != "R" or len(fields) != 8:
            raise AssertionError("a malformed max300 tuple record was emitted")
        roots = tuple(map(int, fields[1:7]))
        flag = int(fields[7])
        if flag not in (0, 1):
            raise AssertionError("a malformed max300 reflection flag was emitted")
        if flag:
            reflection_count += 1
        else:
            nonreflection.append(roots)
    normalized, obstruction, reflected, nonreflected = map(int, summary[3:])
    if (
        obstruction != len(lines) - 2
        or reflected != reflection_count
        or nonreflected != len(nonreflection)
        or reflected + nonreflected != obstruction
    ):
        raise AssertionError("the max300 stream disagrees with its summary")
    return RangeResult(
        first, last, normalized, obstruction, reflected, nonreflected,
        tuple(nonreflection), hashlib.sha256(stdout.encode()).hexdigest(), stdout,
    )


def old_record_lines(stdout: str) -> tuple[str, ...]:
    lines = tuple(line.strip() for line in stdout.splitlines() if line.strip())
    if not lines or lines[0] != "MESTRE_ROOT_TUPLES_V2" or not lines[-1].startswith("S "):
        raise AssertionError("the frozen max200 enumerator output changed")
    return lines[1:-1]


def new_record_lines(result: RangeResult) -> tuple[str, ...]:
    lines = tuple(line.strip() for line in result.stdout.splitlines() if line.strip())
    return lines[1:-1]


def record_root_and_flag(line: str) -> tuple[tuple[int, ...], int]:
    fields = line.split()
    if fields[0] != "R" or len(fields) != 8:
        raise AssertionError("a malformed enumerator record escaped parsing")
    return tuple(map(int, fields[1:7])), int(fields[7])


def digest_record_roots(results: Sequence[RangeResult], flag: int | None = None) -> str:
    digest = hashlib.sha256()
    first = True
    for result in results:
        for line in new_record_lines(result):
            roots, record_flag = record_root_and_flag(line)
            if flag is not None and record_flag != flag:
                continue
            if not first:
                digest.update(b"\n")
            digest.update(",".join(map(str, roots)).encode())
            first = False
    return digest.hexdigest()


def verify_all_records(results: Sequence[RangeResult]) -> None:
    previous: tuple[int, ...] | None = None
    for result in results:
        for line in new_record_lines(result):
            roots, flag = record_root_and_flag(line)
            if (
                len(roots) != 6
                or roots[0] != 0
                or any(roots[index] >= roots[index + 1] for index in range(5))
                or not result.first <= roots[-1] <= result.last
                or gcd(*roots[1:]) != 1
                or normalize_integer_root_tuple(roots) != roots
            ):
                raise AssertionError("the max300 enumerator normalization replay failed")
            if mestre_obstruction_integer(roots) != 0:
                raise AssertionError("the max300 Python obstruction replay failed")
            if reflection_symmetric(roots) != bool(flag):
                raise AssertionError("the max300 Python reflection replay failed")
            if previous is not None and (roots[-1], roots) <= (previous[-1], previous):
                raise AssertionError("the concatenated max300 stream order changed")
            previous = roots


def compile_sources(
    new_source: Path,
    old_source: Path,
    directory: Path,
    *,
    timeout: float,
) -> tuple[Path, Path, dict[str, float]]:
    compiler = shutil.which("c++")
    if compiler is None:
        raise FileNotFoundError("a C++17 compiler is required")
    new_binary = directory / "max300"
    old_binary = directory / "max200"
    started = time.monotonic()
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [
            executor.submit(
                run_capped_process,
                (
                    compiler, "-std=c++17", "-O3", "-DNDEBUG",
                    str(source), "-o", str(binary),
                ),
                timeout=timeout,
            )
            for source, binary in ((new_source, new_binary), (old_source, old_binary))
        ]
        for future in futures:
            future.result()
    return new_binary, old_binary, {"parallel_compile_wall_seconds": time.monotonic() - started}


def enumerate_extension_and_prefix(
    new_binary: Path,
    old_binary: Path,
    *,
    timeout: float,
    workers: int,
) -> tuple[RangeResult, tuple[RangeResult, ...], str, dict[str, Any]]:
    started = time.monotonic()
    with ThreadPoolExecutor(max_workers=workers) as executor:
        extension_futures = [
            executor.submit(
                run_capped_process,
                (str(new_binary), str(first), str(last)),
                timeout=timeout,
            )
            for first, last in EXTENSION_RANGES
        ]
        extension_outputs = [future.result()[0] for future in extension_futures]
    extension_wall = time.monotonic() - started
    started = time.monotonic()
    with ThreadPoolExecutor(max_workers=2) as executor:
        new_future = executor.submit(
            run_capped_process,
            (str(new_binary), "5", str(PRIOR_MAX_ROOT)),
            timeout=timeout,
        )
        old_future = executor.submit(
            run_capped_process,
            (str(old_binary), str(PRIOR_MAX_ROOT)),
            timeout=timeout,
        )
        new_prefix_stdout = new_future.result()[0]
        old_prefix_stdout = old_future.result()[0]
    prefix_wall = time.monotonic() - started
    prefix = parse_range_output(new_prefix_stdout, 5, PRIOR_MAX_ROOT)
    extension = tuple(
        parse_range_output(stdout, first, last)
        for stdout, (first, last) in zip(extension_outputs, EXTENSION_RANGES)
    )
    if new_record_lines(prefix) != old_record_lines(old_prefix_stdout):
        raise AssertionError("the new enumerator's max200 prefix changed a record")
    return prefix, extension, hashlib.sha256(old_prefix_stdout.encode()).hexdigest(), {
        "extension_parallel_wall_seconds": extension_wall,
        "prefix_pair_parallel_wall_seconds": prefix_wall,
        "extension_ranges": [list(bounds) for bounds in EXTENSION_RANGES],
        "extension_range_stdout_sha256": [result.stdout_sha256 for result in extension],
    }


def classify_chunk(
    roots: tuple[tuple[int, ...], ...]
) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...], dict[tuple[int, ...], int]]:
    return generic_classification(roots)


def parallel_generic_classification(
    roots: tuple[tuple[int, ...], ...], workers: int
) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...], dict[tuple[int, ...], int]]:
    if not roots:
        return (), (), {}
    chunks = tuple(
        tuple(roots[index::workers]) for index in range(workers) if roots[index::workers]
    )
    with ProcessPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(classify_chunk, chunks))
    nonsingular_set = {roots for result in results for roots in result[0]}
    singular_set = {roots for result in results for roots in result[1]}
    witnesses = {roots: witness for result in results for roots, witness in result[2].items()}
    nonsingular = tuple(roots for roots in roots if roots in nonsingular_set)
    singular = tuple(roots for roots in roots if roots in singular_set)
    if len(nonsingular) + len(singular) != len(roots) or set(nonsingular) & set(singular):
        raise AssertionError("parallel generic classification lost a family")
    return nonsingular, singular, witnesses


def exclusive_json_write(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")


def frozen_paths(root: Path) -> dict[str, Path]:
    return {
        "census": root / "artifacts/generated-results/elliptic_mestre_root_tuple_scale_max200_census.json",
        "screen": root / "artifacts/generated-results/elliptic_mestre_root_tuple_scale_max200.json",
        "source": root / "elliptic-curves/cas/enumerate_mestre_root_tuples_scale_max200.cpp",
        "driver": root / "elliptic-curves/cas/search_mestre_root_tuple_scale_max200.py",
    }


def load_frozen_max200(root: Path) -> tuple[dict[str, Any], dict[str, str]]:
    paths = frozen_paths(root)
    expected = {
        "census": EXPECTED_MAX200_CENSUS_SHA256,
        "screen": EXPECTED_MAX200_SCREEN_SHA256,
        "source": EXPECTED_MAX200_CPP_SHA256,
        "driver": EXPECTED_MAX200_DRIVER_SHA256,
    }
    actual = {name: sha256_file(path) for name, path in paths.items()}
    if actual != expected:
        raise SystemExit("a frozen max200 input changed")
    census = json.loads(paths["census"].read_text())
    if census["result_sha256"] != EXPECTED_MAX200_CENSUS_RESULT_SHA256:
        raise AssertionError("the frozen max200 census result digest changed")
    return census, actual


def build_census(args: argparse.Namespace, root: Path) -> dict[str, Any]:
    if args.census_output.exists():
        raise SystemExit("refusing to overwrite the max300 census artifact")
    started = time.monotonic()
    frozen, frozen_hashes = load_frozen_max200(root)
    script_path = Path(__file__).resolve()
    source = script_path.with_name("enumerate_mestre_root_tuples_scale_max300.cpp")
    old_source = frozen_paths(root)["source"]
    timings: dict[str, Any] = {}
    try:
        with tempfile.TemporaryDirectory(prefix="mestre-root-300-") as directory:
            new_binary, old_binary, compile_timings = compile_sources(
                source, old_source, Path(directory), timeout=args.compile_timeout
            )
            timings.update(compile_timings)
            prefix, extension, old_stream_sha256, enumeration_timings = (
                enumerate_extension_and_prefix(
                    new_binary,
                    old_binary,
                    timeout=args.enumeration_timeout,
                    workers=args.enumeration_workers,
                )
            )
            timings.update(enumeration_timings)
    except CappedProcessTimeout as error:
        artifact = {
            "schema_version": 1,
            "status": "max300 census capped cleanly before completeness",
            "scope": {"max_root": MAX_ROOT, "complete": False},
            "cap": {
                "compile_timeout_seconds": args.compile_timeout,
                "enumeration_timeout_seconds_per_disjoint_range": args.enumeration_timeout,
                "error": str(error),
                "no_retry": True,
            },
            "frozen_max200_inputs": frozen_hashes,
            "provenance": {
                "script_sha256": sha256_file(script_path),
                "compiled_source_sha256": sha256_file(source),
                "owned_processes_remaining": 0,
            },
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        artifact["result_sha256"] = stable_json_digest(artifact)
        exclusive_json_write(args.census_output, artifact)
        return artifact

    results = (prefix, *extension)
    replay_started = time.monotonic()
    verify_all_records(results)
    timings["independent_record_replay_wall_seconds"] = time.monotonic() - replay_started
    frozen_nonreflection = tuple(
        tuple(map(int, roots)) for roots in frozen["tuple_populations"]["nonreflection_roots"]
    )
    if prefix.nonreflection_roots != frozen_nonreflection:
        raise AssertionError("the max200 nonreflection prefix changed record for record")
    frozen_census = frozen["census"]
    if (
        prefix.normalized_count
        != frozen_census["affine_normalized_primitive_reflection_quotient_count"]
        or prefix.obstruction_count != frozen_census["degree_five_obstruction_zero_count"]
        or prefix.reflection_count != frozen_census["reflection_obstruction_zero_count"]
        or prefix.nonreflection_count != frozen_census["nonreflection_obstruction_zero_count"]
        or digest_record_roots((prefix,)) != frozen_census["obstruction_tuple_sha256"]
        or digest_record_roots((prefix,), 1) != frozen_census["reflection_tuple_sha256"]
        or digest_record_roots((prefix,), 0) != frozen_census["nonreflection_tuple_sha256"]
    ):
        raise AssertionError("the max200 prefix counts or tuple digests changed")

    extension_nonreflection = tuple(
        roots for result in extension for roots in result.nonreflection_roots
    )
    classification_started = time.monotonic()
    extension_nonsingular, extension_singular, witnesses = parallel_generic_classification(
        extension_nonreflection, args.workers
    )
    timings["generic_classification_wall_seconds"] = time.monotonic() - classification_started
    if set(witnesses.values()) != ({1} if witnesses else set()):
        raise AssertionError("a new nonsingular family lost the T=1 witness pattern")
    old_nonsingular = tuple(
        tuple(map(int, roots))
        for roots in frozen["tuple_populations"]["generically_nonsingular_nonreflection_roots"]
    )
    old_singular = tuple(
        tuple(map(int, roots))
        for roots in frozen["tuple_populations"]["generically_singular_nonreflection_roots"]
    )
    full_nonsingular = old_nonsingular + extension_nonsingular
    full_singular = old_singular + extension_singular
    independent = independent_normalized_count(MAX_ROOT)
    total_normalized = sum(result.normalized_count for result in results)
    total_obstruction = sum(result.obstruction_count for result in results)
    total_reflection = sum(result.reflection_count for result in results)
    total_nonreflection = sum(result.nonreflection_count for result in results)
    if total_normalized != independent["primitive_reflection_orbit_count"]:
        raise AssertionError("the max300 enumerator disagrees with Burnside/Mobius")
    if total_reflection != independent["primitive_reflection_fixed_count"]:
        raise AssertionError("the max300 reflection family count changed")
    if total_nonreflection != len(full_nonsingular) + len(full_singular):
        raise AssertionError("the max300 generic classification count changed")
    if any(roots[-1] <= PRIOR_MAX_ROOT for roots in extension_nonsingular):
        raise AssertionError("a max200 family escaped prefix exclusion")

    timings["total_wall_seconds"] = time.monotonic() - started
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "complete exact max-root-300 census; no specialization arithmetic",
        "scope": {
            "max_root": MAX_ROOT,
            "complete_diameter_prefix": [5, MAX_ROOT],
            "open_diameter_remainder": [],
            "specialization_arithmetic_run": False,
        },
        "frozen_max200_inputs": {
            **frozen_hashes,
            "all_frozen_files_read_only": True,
        },
        "exact_max200_prefix_recovery": {
            "record_for_record_equal_to_frozen_enumerator": True,
            "new_prefix_stream_sha256": prefix.stdout_sha256,
            "frozen_rerun_stream_sha256": old_stream_sha256,
            "counts": [
                prefix.normalized_count, prefix.obstruction_count,
                prefix.reflection_count, prefix.nonreflection_count,
            ],
            "obstruction_tuple_sha256": digest_record_roots((prefix,)),
            "reflection_tuple_sha256": digest_record_roots((prefix,), 1),
            "nonreflection_tuple_sha256": digest_record_roots((prefix,), 0),
            "nonreflection_records_equal_frozen_census": True,
        },
        "census": {
            "affine_normalized_primitive_reflection_quotient_count": total_normalized,
            "degree_five_obstruction_zero_count": total_obstruction,
            "reflection_obstruction_zero_count": total_reflection,
            "nonreflection_obstruction_zero_count": total_nonreflection,
            "nonreflection_generically_nonsingular_count": len(full_nonsingular),
            "nonreflection_generically_singular_count": len(full_singular),
            "genuinely_new_diameter_201_to_300_family_count": len(extension_nonsingular),
            "genuinely_new_diameter_201_to_300_family_sha256": tuple_digest(extension_nonsingular),
            "obstruction_tuple_sha256": digest_record_roots(results),
            "reflection_tuple_sha256": digest_record_roots(results, 1),
            "nonreflection_tuple_sha256": digest_record_roots(results, 0),
            "nonsingular_nonreflection_tuple_sha256": tuple_digest(full_nonsingular),
            "singular_nonreflection_tuple_sha256": tuple_digest(full_singular),
            "extension_nonreflection_tuple_sha256": tuple_digest(extension_nonreflection),
            "all_normalization_reflection_and_obstruction_gates_replayed_in_Python": True,
            "generic_singularity_test_parameter_count": 21,
            "generic_discriminant_degree_upper_bound": 20,
            "all_new_nonsingularity_witness_parameters": sorted(set(witnesses.values())),
            "independent_burnside_mobius_count": independent,
            "absolute_obstruction_expression_bound": str(47_520 * MAX_ROOT**5),
            "bound_is_strictly_below_signed_128_max": 47_520 * MAX_ROOT**5 < 2**127,
        },
        "tuple_populations": {
            "genuinely_new_nonsingular_roots": [list(roots) for roots in extension_nonsingular],
            "genuinely_new_singular_roots": [list(roots) for roots in extension_singular],
            "genuinely_new_nonreflection_roots": [list(roots) for roots in extension_nonreflection],
        },
        "parameters": {
            "compile_timeout_seconds": args.compile_timeout,
            "enumeration_timeout_seconds_per_disjoint_range": args.enumeration_timeout,
            "enumeration_workers": args.enumeration_workers,
            "classification_workers": args.workers,
            "screen_family_cap": SCREEN_FAMILY_CAP,
        },
        "timings": timings,
        "provenance": {
            "script": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "compiled_source": str(source.relative_to(root)),
            "compiled_source_sha256": sha256_file(source),
            "reproducing_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "all_range_processes_use_disjoint_diameter_intervals": True,
            "subprocesses_run_in_foreground_process_groups": True,
            "temporary_enumerator_binaries_removed": True,
            "owned_processes_remaining": 0,
        },
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    artifact["result_sha256"] = stable_json_digest(
        {
            "prefix": artifact["exact_max200_prefix_recovery"],
            "census": artifact["census"],
            "new_roots": artifact["tuple_populations"]["genuinely_new_nonsingular_roots"],
        }
    )
    exclusive_json_write(args.census_output, artifact)
    return artifact


def select_rank_aware_diversity_leaders(
    family_records: Sequence[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    best_by_family = []
    for family in family_records:
        admissible = [record for record in family["fibers"] if record["admissible"]]
        if not admissible:
            continue
        best = dict(sorted(admissible, key=fiber_rank_key)[0])
        best.update(
            {
                "roots": family["roots"],
                "diameter": family["diameter"],
                "diameter_decile": family["diameter_decile"],
            }
        )
        best_by_family.append(best)
    if len(best_by_family) < SELECTED_FAMILY_COUNT:
        raise AssertionError("fewer than 64 new max300 families have an admissible fiber")
    selected: dict[tuple[int, ...], dict[str, Any]] = {}
    for record in sorted(best_by_family, key=fiber_rank_key)[:GLOBAL_RANK_AWARE_FAMILY_KEEP]:
        item = dict(record)
        item["selection_stratum"] = "top-34 exact mod-3 rank/local-score leader"
        selected[tuple(item["roots"])] = item
    diversity_counts = {}
    for lower in range(201, 301, 10):
        upper = lower + 9
        pool = [
            record for record in best_by_family
            if lower <= record["diameter"] <= upper
            and tuple(record["roots"]) not in selected
        ]
        chosen = sorted(pool, key=fiber_rank_key)[:DIVERSITY_KEEP_PER_DECILE]
        if len(chosen) != DIVERSITY_KEEP_PER_DECILE:
            raise AssertionError(f"diameter decile {lower}-{upper} lacks leaders")
        diversity_counts[f"{lower}-{upper}"] = len(chosen)
        for record in chosen:
            item = dict(record)
            item["selection_stratum"] = f"rank-aware diversity-{lower}-{upper}"
            selected[tuple(item["roots"])] = item
    result = sorted(selected.values(), key=lambda record: record["identifier"])
    if len(result) != SELECTED_FAMILY_COUNT:
        raise AssertionError("the max300 rank-aware tranche must contain 64 families")
    return result, {
        "one_fiber_per_family": True,
        "global_rank_aware_family_keep": GLOBAL_RANK_AWARE_FAMILY_KEEP,
        "diversity_keep_per_diameter_decile": DIVERSITY_KEEP_PER_DECILE,
        "diversity_decile_counts": diversity_counts,
        "selected_family_count": len(result),
        "selected_identifier_sha256": hashlib.sha256(
            "\n".join(record["identifier"] for record in result).encode()
        ).hexdigest(),
    }


def screen_result_digest(artifact: dict[str, Any]) -> str:
    return stable_json_digest(
        {
            "checkpoint": artifact["input"]["census_checkpoint_sha256"],
            "population": artifact["complete_panel_screen"]["population"],
            "screen_digest": artifact["complete_panel_screen"]["family_records_sha256"],
            "selection": artifact["rank_aware_diversity_selection"]["selected_identifier_sha256"],
            "followup": [
                [
                    record["identifier"],
                    record["conductor_phase"]["status"],
                    record["conductor_phase"].get("conductor"),
                    record["point_triage"]["status"],
                    record["point_triage"].get("stable_numerical_rank"),
                    record["point_triage"].get("pool_point_sha256"),
                ]
                for record in artifact["leader_followup"]["records"]
            ],
            "target": artifact["target"],
        }
    )


def build_screen(args: argparse.Namespace, root: Path) -> dict[str, Any]:
    if args.output.exists():
        raise SystemExit("refusing to overwrite the max300 screen artifact")
    frozen, frozen_hashes = load_frozen_max200(root)
    checkpoint = json.loads(args.census_output.read_text())
    checkpoint_sha256 = sha256_file(args.census_output)
    source = Path(__file__).with_name("enumerate_mestre_root_tuples_scale_max300.cpp")
    if (
        checkpoint.get("status")
        != "complete exact max-root-300 census; no specialization arithmetic"
        or checkpoint["scope"]["max_root"] != MAX_ROOT
        or checkpoint["provenance"]["compiled_source_sha256"] != sha256_file(source)
        or checkpoint["frozen_max200_inputs"] != {
            **frozen_hashes, "all_frozen_files_read_only": True
        }
    ):
        raise AssertionError("the max300 census checkpoint is incompatible")
    new_roots = tuple(
        tuple(map(int, roots))
        for roots in checkpoint["tuple_populations"]["genuinely_new_nonsingular_roots"]
    )
    census = checkpoint["census"]
    if (
        len(new_roots) != census["genuinely_new_diameter_201_to_300_family_count"]
        or tuple_digest(new_roots)
        != census["genuinely_new_diameter_201_to_300_family_sha256"]
        or any(not PRIOR_MAX_ROOT < roots[-1] <= MAX_ROOT for roots in new_roots)
    ):
        raise AssertionError("the max300 new-family checkpoint changed")
    if len(new_roots) > SCREEN_FAMILY_CAP:
        artifact = {
            "schema_version": 1,
            "status": "max300 panel screen capped before launch: family census excessive",
            "input": {
                "census_checkpoint_sha256": checkpoint_sha256,
                "new_family_count": len(new_roots),
            },
            "cap": {"screen_family_cap": SCREEN_FAMILY_CAP, "calls_launched": 0},
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        }
        artifact["result_sha256"] = stable_json_digest(artifact)
        exclusive_json_write(args.output, artifact)
        return artifact

    started = time.monotonic()
    if args.workers == 1:
        family_records = [screen_family(roots) for roots in new_roots]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            family_records = list(executor.map(screen_family, new_roots, chunksize=1))
    screen_wall = time.monotonic() - started
    if [tuple(record["roots"]) for record in family_records] != list(new_roots):
        raise AssertionError("parallel max300 panel screening changed order")
    selected, selection = select_rank_aware_diversity_leaders(family_records)
    selected_by_id = {record["identifier"]: record for record in selected}
    followup = []
    runtime: dict[str, tuple[tuple[int, ...], int]] = {}
    conductor_started = time.monotonic()
    for identifier in sorted(selected_by_id):
        selection_record = selected_by_id[identifier]
        roots = tuple(selection_record["roots"])
        parameter = int(selection_record["parameter"])
        record: dict[str, Any] = {
            "identifier": identifier,
            "roots": list(roots),
            "parameter": parameter,
            "selection_stratum": selection_record["selection_stratum"],
            "panel_visible_rank_lower_bound": selection_record[
                "mod3_finite_reduction_certificate"
            ]["certified_algebraic_rank_lower_bound"],
        }
        try:
            _, coefficients, visible_points = visible_points_and_coefficients(roots, parameter)
            record["exact_visible_points"] = [point_record(point) for point in visible_points]
            conductor = capped_minimal_curve_data(
                coefficients, timeout=args.conductor_timeout, stack_bytes=STACK_BYTES
            )
            record["conductor_phase"] = {
                "status": "completed exact PARI minimal-model/conductor computation",
                **conductor,
                "below_strict_log_conductor_target_numerically": (
                    Decimal(conductor["log_conductor"]) < TARGET_LOG_CONDUCTOR
                ),
            }
            runtime[identifier] = roots, parameter
        except CappedProcessTimeout:
            record["conductor_phase"] = {
                "status": "timeout-no-retry", "timeout_seconds": args.conductor_timeout
            }
        except Exception as error:
            record["conductor_phase"] = {
                "status": "error-no-retry", "error": str(error)[:1000]
            }
        followup.append(record)
    conductor_wall = time.monotonic() - conductor_started

    point_started = time.monotonic()
    target_hits = []
    for position, record in enumerate(followup, start=1):
        identifier = record["identifier"]
        if identifier not in runtime:
            record["point_triage"] = {"status": "not attempted after incomplete conductor"}
            record["immediate_exact_gain_attempt"] = {"status": "not attempted"}
            continue
        roots, parameter = runtime[identifier]
        try:
            triage, subset = search_h5000(
                roots,
                parameter,
                point_timeout=args.point_timeout,
                height_timeout=args.height_timeout,
            )
            record["point_triage"] = triage
            stable_rank = int(triage["stable_numerical_rank"])
            if stable_rank >= STRONG_GAIN_TRIGGER and subset is not None:
                print(
                    f"EARLY_SIGNAL {identifier} stable_rank={stable_rank}", flush=True
                )
                _, coefficients, _ = visible_points_and_coefficients(roots, parameter)
                mod3 = mod3_independence_certificate(
                    coefficients, subset, prime_bound=EXACT_GAIN_PRIME_BOUND
                )
                record["immediate_exact_gain_attempt"] = {"mod3": mod3}
                certified = mod3["certified_algebraic_rank_lower_bound"]
                if certified < len(subset):
                    mod2 = finite_reduction_attempt(
                        coefficients, subset, prime_bound=EXACT_GAIN_PRIME_BOUND
                    )
                    record["immediate_exact_gain_attempt"]["mod2"] = mod2
                    certified = max(
                        certified, mod2["certified_algebraic_rank_lower_bound"] or 0
                    )
                record["immediate_exact_gain_attempt"][
                    "best_certified_algebraic_rank_lower_bound"
                ] = certified
                print(
                    f"EXACT_SIGNAL {identifier} certified_rank={certified}", flush=True
                )
                conductor = record["conductor_phase"]
                if certified >= 30 or (
                    certified >= 21
                    and conductor["below_strict_log_conductor_target_numerically"]
                ):
                    target_hits.append(
                        {
                            "identifier": identifier,
                            "certified_algebraic_rank_lower_bound": certified,
                            "conductor": conductor["conductor"],
                            "log_conductor": conductor["log_conductor"],
                        }
                    )
            else:
                record["immediate_exact_gain_attempt"] = {
                    "status": "not triggered",
                    "trigger_stable_numerical_rank": STRONG_GAIN_TRIGGER,
                }
        except CappedProcessTimeout:
            record["point_triage"] = {
                "status": "timeout-no-retry",
                "point_timeout_seconds": args.point_timeout,
                "height_timeout_seconds": args.height_timeout,
            }
            record["immediate_exact_gain_attempt"] = {"status": "not attempted"}
        except Exception as error:
            record["point_triage"] = {
                "status": "error-no-retry", "error": str(error)[:1000]
            }
            record["immediate_exact_gain_attempt"] = {"status": "not attempted"}
        if position % 16 == 0:
            print(f"H5000 {position}/{len(followup)}", flush=True)
    point_wall = time.monotonic() - point_started

    all_fibers = [fiber for family in family_records for fiber in family["fibers"]]
    admissible = [fiber for fiber in all_fibers if fiber["admissible"]]
    completed_points = [
        record for record in followup
        if record["point_triage"]["status"].startswith("completed")
    ]
    rank_histogram = Counter(
        str(fiber["mod3_finite_reduction_certificate"]["certified_algebraic_rank_lower_bound"])
        for fiber in admissible
    )
    script_path = Path(__file__).resolve()
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "complete max-root-300 exact panel screen and bounded leader followup",
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "alternative_rank_at_least": 30,
            "hits": target_hits,
        },
        "input": {
            "census_checkpoint": str(args.census_output.relative_to(root)),
            "census_checkpoint_sha256": checkpoint_sha256,
            "census_result_sha256": checkpoint["result_sha256"],
            "frozen_max200_inputs": frozen_hashes,
        },
        "complete_panel_screen": {
            "protocol": {
                "integer_parameters": [1, 8],
                "every_genuinely_new_nonsingular_nonreflection_family_screened": True,
                "every_admissible_panel_fiber_maps_all_visible_points_exactly": True,
                "every_admissible_panel_fiber_receives_exact_mod3_finite_reduction_signature": True,
                "finite_reduction_prime_bound": 251,
                "finite_reduction_descent_modulus": 3,
                "selection_uses_exact_finite_reduction_rank": True,
                "selection_uses_visible_count_alone": False,
                "selection_uses_conductor": False,
            },
            "population": {
                "new_family_count": len(family_records),
                "proposed_panel_fiber_count": len(all_fibers),
                "admissible_panel_fiber_count": len(admissible),
                "inadmissible_panel_fiber_count": len(all_fibers) - len(admissible),
                "exact_mod3_certificate_count": len(admissible),
                "visible_rank_lower_bound_histogram": dict(sorted(rank_histogram.items())),
                "maximum_visible_certified_rank_lower_bound": max(
                    fiber["mod3_finite_reduction_certificate"]["certified_algebraic_rank_lower_bound"]
                    for fiber in admissible
                ),
            },
            "family_records_sha256": stable_json_digest(family_records),
            "family_records": family_records,
        },
        "rank_aware_diversity_selection": {
            "population_closed_before_conductor_calls": True,
            "selection_key": (
                "descending exact mod-3 visible rank lower bound, descending fixed-prime "
                "local score/coverage, ascending coefficient height/parameter/id"
            ),
            **selection,
            "selected_records": selected,
        },
        "leader_followup": {
            "protocol": {
                "all_selected_leaders_receive_conductor_first": True,
                "conductor_population_closed_before_any_point_or_height_call": True,
                "bounded_point_height": 5_000,
                "no_retries": True,
                "immediate_exact_gain_trigger_stable_numerical_rank": STRONG_GAIN_TRIGGER,
            },
            "population": {
                "selected_leaders": len(followup),
                "conductor_completed": sum(
                    record["conductor_phase"]["status"].startswith("completed")
                    for record in followup
                ),
                "conductor_timeouts": sum(
                    record["conductor_phase"]["status"] == "timeout-no-retry"
                    for record in followup
                ),
                "conductor_errors": sum(
                    record["conductor_phase"]["status"] == "error-no-retry"
                    for record in followup
                ),
                "subtarget_conductors": sum(
                    record["conductor_phase"].get(
                        "below_strict_log_conductor_target_numerically"
                    ) is True
                    for record in followup
                ),
                "point_search_completed": len(completed_points),
                "point_search_timeouts": sum(
                    record["point_triage"]["status"] == "timeout-no-retry"
                    for record in followup
                ),
                "point_search_errors": sum(
                    record["point_triage"]["status"] == "error-no-retry"
                    for record in followup
                ),
                "maximum_stable_numerical_rank": max(
                    (record["point_triage"]["stable_numerical_rank"] for record in completed_points),
                    default=None,
                ),
                "stable_numerical_rank_histogram": dict(
                    sorted(Counter(
                        str(record["point_triage"]["stable_numerical_rank"])
                        for record in completed_points
                    ).items())
                ),
                "immediate_exact_gain_attempts": sum(
                    record["immediate_exact_gain_attempt"].get("status")
                    not in {"not triggered", "not attempted"}
                    for record in followup
                ),
            },
            "records": followup,
        },
        "parameters": {
            "workers": args.workers,
            "conductor_timeout_seconds": args.conductor_timeout,
            "point_timeout_seconds": args.point_timeout,
            "height_timeout_seconds": args.height_timeout,
            "stack_bytes": STACK_BYTES,
            "exact_gain_prime_bound": EXACT_GAIN_PRIME_BOUND,
        },
        "timings": {
            "complete_exact_panel_screen_wall_seconds": screen_wall,
            "conductor_phase_wall_seconds": conductor_wall,
            "point_phase_wall_seconds": point_wall,
            "total_wall_seconds": time.monotonic() - started,
        },
        "provenance": {
            "script": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "reproducing_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "worker_pool_joined_before_write": True,
            "external_subprocesses_run_in_foreground_process_groups": True,
            "owned_processes_remaining": 0,
        },
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    artifact["result_sha256"] = screen_result_digest(artifact)
    exclusive_json_write(args.output, artifact)
    return artifact


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    generated = root / "artifacts/generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("census", "screen"), required=True)
    parser.add_argument("--compile-timeout", type=float, default=30.0)
    parser.add_argument("--enumeration-timeout", type=float, default=90.0)
    parser.add_argument("--enumeration-workers", type=int, default=4)
    parser.add_argument("--conductor-timeout", type=float, default=8.0)
    parser.add_argument("--point-timeout", type=float, default=12.0)
    parser.add_argument("--height-timeout", type=float, default=12.0)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--census-output",
        type=Path,
        default=generated / "elliptic_mestre_root_tuple_scale_max300_census.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_mestre_root_tuple_scale_max300.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 1 <= args.workers <= 8 or not 1 <= args.enumeration_workers <= 4:
        raise SystemExit("worker counts are outside their declared caps")
    if not 0 < args.compile_timeout <= 60 or not 0 < args.enumeration_timeout <= 120:
        raise SystemExit("compile/enumeration timeout is outside its declared cap")
    if any(
        timeout <= 0 or timeout > 30
        for timeout in (args.conductor_timeout, args.point_timeout, args.height_timeout)
    ):
        raise SystemExit("PARI subprocess caps must lie in (0,30]")
    root = Path(__file__).resolve().parents[2]
    result = build_census(args, root) if args.stage == "census" else build_screen(args, root)
    if args.stage == "census":
        print(
            json.dumps(
                {
                    "status": result["status"],
                    "census": result.get("census"),
                    "result_sha256": result["result_sha256"],
                },
                sort_keys=True,
            ),
            flush=True,
        )
    else:
        print(
            json.dumps(
                {
                    "status": result["status"],
                    "panel_population": result.get("complete_panel_screen", {}).get("population"),
                    "leader_population": result.get("leader_followup", {}).get("population"),
                    "target_hits": result.get("target", {}).get("hits"),
                    "result_sha256": result["result_sha256"],
                },
                sort_keys=True,
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
