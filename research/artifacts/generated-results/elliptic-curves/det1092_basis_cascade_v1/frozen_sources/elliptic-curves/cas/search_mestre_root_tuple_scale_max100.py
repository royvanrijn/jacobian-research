#!/usr/bin/env python3
"""Exact max-root-100 Mestre census and a rank-blind bounded fiber tranche.

The max-root-50 C++ enumerator, driver, test, and artifact are frozen inputs.
This standalone continuation compiles that unchanged enumerator and exhausts
the full bound 100.  Every emitted obstruction-zero tuple is independently
replayed in Python; nonreflection tuples then receive the exact 21-value
generic nonsingularity test from the frozen driver.  Families of diameter at
most 50 are excluded exactly.

The specialization tranche is chosen without conductor, point, or rank data.
Twenty families maximize a fixed small-prime local score; four additional
families in each diameter decile minimize a deterministic geometry key.  All
T=1,...,8 fibers of those forty families pass through exact admissibility and
then one conductor computation.  Before conductors are run, sixty-four fibers
are fixed for one H=5000 point search: the strongest local fiber in every
family plus the strongest second fiber in twenty-four distinct families.
Numerical height ranks are triage only.  Stable numerical rank at least 21
triggers an exact finite-reduction independence attempt immediately.
"""

from __future__ import annotations

import argparse
from collections import Counter
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

from ek_k3 import fraction_mod, legendre_symbol, rational_to_string
from mestre_root_tuples import SixRootMestreConstruction
from search_mestre_root_tuple_scale import (
    CappedProcessTimeout,
    EnumerationResult,
    TARGET_LOG_CONDUCTOR,
    bounded_quartic_points,
    canonical_signless_points,
    capped_minimal_curve_data,
    classify_nonreflection,
    finite_reduction_attempt,
    height_matrix_replay,
    numerical_subset,
    point_digest,
    point_record,
    primitive_visible_points,
    quartic_point_to_jacobian,
    quartic_value,
    run_capped_process,
    sha256_file,
    tuple_digest,
    verify_enumerator_records,
)


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

MAX_ROOT = 100
PRIOR_MAX_ROOT = 50
PARAMETERS = tuple(range(1, 9))
LOCAL_PRIMES = (11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
MINIMUM_LOCAL_COVERAGE = 8
LOCAL_FAMILY_KEEP = 20
GEOMETRY_KEEP_PER_DECILE = 4
H5000_EXTRA_SECOND_FAMILY_KEEP = 24
H5000_HEIGHT = 5_000
H5000_MAPPING_CAP = 128
STACK_BYTES = 256_000_000
CERTIFICATE_PRIME_BOUND = 500

FROZEN_CPP_SHA256 = (
    "31650333800698201819eddc91bf228089824bca026c629c9360683324a69eb5"
)
FROZEN_DRIVER_SHA256 = (
    "5e7228b95ae995019fbc50b9f7667de41e06a86b4490f0feacff5702bb5cc174"
)
FROZEN_TEST_SHA256 = (
    "a3930892e7e574161c0713c6c9b7c9f5aee0aa74e8e7acb89c250e2f9975d7c3"
)
FROZEN_ARTIFACT_SHA256 = (
    "fd2dccb1fd08aad70857df7ca19df77bd521e2be017b98f5579a748fd26cfc14"
)
EXPECTED_MAX100_COUNTS = (36_475_792, 33_945, 33_168, 777)
EXPECTED_MAX100_OBSTRUCTION_SHA256 = (
    "52f938de75951011526c59355036f4ed377c70617fc8826e9a1a29a178c5152d"
)
EXPECTED_MAX100_NONREFLECTION_SHA256 = (
    "87d7aa0d5a8fc2160d2dd6b8d8dc0eec68a54b73f83b09f4d4f0817ef42d775a"
)
EXPECTED_MAX100_NONSINGULAR_SHA256 = (
    "e92e9cd0be8fc8006275797df2752b714df0237ae27ce2b3ba4829c988681973"
)
EXPECTED_MAX50_NONSINGULAR_SHA256 = (
    "a892b1824af10c5e2dd428478778e8cc8db69fa057cfde03dbc10f1ab273433e"
)
EXPECTED_NEW_FAMILY_SHA256 = (
    "2a631570f830ac3da46dc0b6f414ea93614805d68de9705b2c2c5d8cb170971b"
)


def stable_json_digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def fraction_record(value: Fraction) -> dict[str, int | str]:
    value = Q(value)
    return {
        "value": rational_to_string(value),
        "numerator": value.numerator,
        "denominator": value.denominator,
    }


def parse_enumerator_output(stdout: str, max_root: int) -> EnumerationResult:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not lines or lines[0] != "MESTRE_ROOT_TUPLES_V1":
        raise AssertionError("the frozen enumerator omitted its format header")
    roots: list[tuple[int, ...]] = []
    reflection: list[tuple[int, ...]] = []
    nonreflection: list[tuple[int, ...]] = []
    summary: tuple[int, ...] | None = None
    for line in lines[1:]:
        fields = line.split()
        if fields[0] == "R":
            if len(fields) != 8:
                raise AssertionError("malformed frozen-enumerator root record")
            item = tuple(int(value) for value in fields[1:7])
            flag = int(fields[7])
            if flag not in (0, 1):
                raise AssertionError("malformed frozen-enumerator reflection flag")
            roots.append(item)
            (reflection if flag else nonreflection).append(item)
        elif fields[0] == "S":
            if len(fields) != 6 or summary is not None:
                raise AssertionError("malformed frozen-enumerator summary")
            summary = tuple(int(value) for value in fields[1:])
        else:
            raise AssertionError("unknown frozen-enumerator output record")
    if summary is None:
        raise AssertionError("the frozen enumerator omitted its summary")
    declared, normalized, obstruction, reflected, nonreflected = summary
    if (
        declared != max_root
        or obstruction != len(roots)
        or reflected != len(reflection)
        or nonreflected != len(nonreflection)
        or reflected + nonreflected != obstruction
    ):
        raise AssertionError("the frozen-enumerator summary disagrees with its stream")
    if roots != sorted(roots, key=lambda item: (item[-1], item)):
        raise AssertionError("the frozen-enumerator order changed")
    return EnumerationResult(
        max_root=max_root,
        normalized_count=normalized,
        obstruction_count=obstruction,
        reflection_count=reflected,
        nonreflection_count=nonreflected,
        obstruction_roots=tuple(roots),
        reflection_roots=tuple(reflection),
        nonreflection_roots=tuple(nonreflection),
    )


def compiled_enumeration_max100(
    source: Path,
    *,
    compile_timeout: float,
    enumeration_timeout: float,
) -> tuple[EnumerationResult, dict[str, float]]:
    compiler = shutil.which("c++")
    if compiler is None:
        raise FileNotFoundError("a C++17 compiler is required")
    timings: dict[str, float] = {}
    with tempfile.TemporaryDirectory(prefix="mestre-root-100-") as directory:
        binary = Path(directory) / "enumerator"
        started = time.monotonic()
        run_capped_process(
            (
                compiler,
                "-std=c++17",
                "-O3",
                "-DNDEBUG",
                str(source),
                "-o",
                str(binary),
            ),
            timeout=compile_timeout,
        )
        timings["compile_wall_seconds"] = time.monotonic() - started
        started = time.monotonic()
        stdout, _ = run_capped_process(
            (str(binary), str(MAX_ROOT)), timeout=enumeration_timeout
        )
        timings["enumeration_wall_seconds"] = time.monotonic() - started
    return parse_enumerator_output(stdout, MAX_ROOT), timings


def curve_discriminant(coefficients: Sequence[Fraction]) -> Fraction:
    a1, a2, a3, a4, a6 = map(Q, coefficients)
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
    return -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6


def local_trace(coefficients: Sequence[Fraction], prime: int) -> int | None:
    discriminant = curve_discriminant(coefficients)
    try:
        if fraction_mod(discriminant, prime) == 0:
            return None
        reduced = [fraction_mod(Q(value), prime) for value in coefficients]
    except ValueError:
        return None
    a1, a2, a3, a4, a6 = reduced
    character_sum = 0
    for x_value in range(prime):
        linear_y = (a1 * x_value + a3) % prime
        rhs = (
            x_value**3 + a2 * x_value**2 + a4 * x_value + a6
        ) % prime
        character_sum += legendre_symbol(linear_y**2 + 4 * rhs, prime)
    return -character_sum


def integer_radical(value: int) -> int:
    value = abs(value)
    radical = 1
    prime = 2
    while prime * prime <= value:
        if value % prime == 0:
            radical *= prime
            while value % prime == 0:
                value //= prime
        prime += 1 if prime == 2 else 2
    if value > 1:
        radical *= value
    return radical


def geometry_features(roots: tuple[int, ...]) -> dict[str, int]:
    vandermonde = 1
    for right in range(1, len(roots)):
        for left in range(right):
            vandermonde *= roots[right] - roots[left]
    radical = integer_radical(vandermonde)
    return {
        "root_difference_vandermonde": vandermonde,
        "root_difference_radical": radical,
        "root_difference_powerful_savings": vandermonde // radical,
    }


def coefficient_height(coefficients: Sequence[Fraction]) -> int:
    return max(
        max(abs(Q(value).numerator), Q(value).denominator)
        for value in coefficients
    )


def family_feature_record(roots: tuple[int, ...]) -> dict[str, Any]:
    construction = SixRootMestreConstruction(tuple(Q(root) for root in roots))
    parameters: list[dict[str, Any]] = []
    for integer_parameter in PARAMETERS:
        parameter = Q(integer_parameter)
        discriminant = construction.quartic_discriminant(parameter)
        degeneracy = construction.visible_point_degeneracy(parameter)
        admissible = (
            discriminant != 0
            and degeneracy.collision_loss == 0
            and degeneracy.zero_ordinates == 0
        )
        record: dict[str, Any] = {
            "parameter": integer_parameter,
            "admissible": admissible,
            "collision_loss": degeneracy.collision_loss,
            "zero_ordinates": degeneracy.zero_ordinates,
        }
        if admissible:
            coefficients = construction.primitive_jacobian_coefficients(parameter)
            traces = []
            score = Q(0)
            for prime in LOCAL_PRIMES:
                trace = local_trace(coefficients, prime)
                if trace is None:
                    continue
                traces.append((prime, trace))
                score += Q(2 - trace, prime + 1 - trace)
            record.update(
                {
                    "coefficient_height": coefficient_height(coefficients),
                    "good_prime_coverage": len(traces),
                    "local_score": fraction_record(score),
                    "local_traces": [
                        {"prime": prime, "trace": trace}
                        for prime, trace in traces
                    ],
                }
            )
        parameters.append(record)
    scored = [
        record
        for record in parameters
        if record["admissible"]
        and record["good_prime_coverage"] >= MINIMUM_LOCAL_COVERAGE
    ]
    scored.sort(
        key=lambda record: (
            -Q(record["local_score"]["value"]),
            -record["good_prime_coverage"],
            record["coefficient_height"],
            record["parameter"],
        )
    )
    geometry = geometry_features(roots)
    result = {
        "roots": list(roots),
        "diameter": roots[-1],
        "diameter_decile": f"{10 * ((roots[-1] - 1) // 10) + 1}-{10 * ((roots[-1] - 1) // 10) + 10}",
        **geometry,
        "admissible_parameter_count": sum(
            record["admissible"] for record in parameters
        ),
        "scorable_parameter_count": len(scored),
        "minimum_coefficient_height": min(
            (
                record["coefficient_height"]
                for record in parameters
                if record["admissible"]
            ),
            default=None,
        ),
        "parameter_features": parameters,
    }
    if scored:
        result.update(
            {
                "best_local_parameter": scored[0]["parameter"],
                "best_local_score": scored[0]["local_score"],
                "second_local_score": scored[min(1, len(scored) - 1)][
                    "local_score"
                ],
                "minimum_good_prime_coverage": min(
                    record["good_prime_coverage"] for record in scored
                ),
            }
        )
    else:
        result.update(
            {
                "best_local_parameter": None,
                "best_local_score": None,
                "second_local_score": None,
                "minimum_good_prime_coverage": None,
            }
        )
    return result


def local_family_key(record: dict[str, Any]) -> tuple[Any, ...]:
    return (
        -Q(record["best_local_score"]["value"]),
        -Q(record["second_local_score"]["value"]),
        -record["minimum_good_prime_coverage"],
        record["minimum_coefficient_height"],
        tuple(record["roots"]),
    )


def geometry_family_key(record: dict[str, Any]) -> tuple[Any, ...]:
    return (
        -record["root_difference_powerful_savings"],
        record["root_difference_radical"],
        record["minimum_coefficient_height"],
        tuple(record["roots"]),
    )


def select_family_tranche(
    records: Sequence[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    eligible = [record for record in records if record["scorable_parameter_count"]]
    if len(eligible) < 40:
        raise AssertionError("fewer than forty families pass the fixed local coverage gate")
    selected: dict[tuple[int, ...], dict[str, Any]] = {}
    for record in sorted(eligible, key=local_family_key)[:LOCAL_FAMILY_KEEP]:
        copy = dict(record)
        copy["selection_stratum"] = "top-20 fixed-panel local score"
        selected[tuple(record["roots"])] = copy
    geometry_counts: dict[str, int] = {}
    for lower in (51, 61, 71, 81, 91):
        upper = lower + 9
        pool = [
            record
            for record in eligible
            if lower <= record["diameter"] <= upper
            and tuple(record["roots"]) not in selected
        ]
        chosen = sorted(pool, key=geometry_family_key)[:GEOMETRY_KEEP_PER_DECILE]
        if len(chosen) != GEOMETRY_KEEP_PER_DECILE:
            raise AssertionError("a diameter decile lacked four geometry candidates")
        label = f"{lower}-{upper}"
        geometry_counts[label] = len(chosen)
        for record in chosen:
            copy = dict(record)
            copy["selection_stratum"] = f"geometry-{label}"
            selected[tuple(record["roots"])] = copy
    result = sorted(selected.values(), key=lambda record: tuple(record["roots"]))
    if len(result) != 40:
        raise AssertionError("the rank-blind family tranche must contain forty families")
    return result, {
        "selection_eligible_family_count": len(eligible),
        "local_coverage_ineligible_family_count": len(records) - len(eligible),
        "local_keep": LOCAL_FAMILY_KEEP,
        "geometry_keep_per_decile": GEOMETRY_KEEP_PER_DECILE,
        "geometry_decile_counts": geometry_counts,
        "selected_family_count": len(result),
        "selected_family_sha256": tuple_digest(
            tuple(tuple(record["roots"]) for record in result)
        ),
    }


def parameter_local_key(record: dict[str, Any]) -> tuple[Any, ...]:
    return (
        -Q(record["local_score"]["value"]),
        -record["good_prime_coverage"],
        record["coefficient_height"],
        record["parameter"],
    )


def predeclare_h5000_fibers(
    selected: Sequence[dict[str, Any]],
) -> tuple[list[str], dict[str, Any]]:
    primary: list[tuple[tuple[int, ...], dict[str, Any]]] = []
    second: list[tuple[tuple[int, ...], dict[str, Any]]] = []
    for family in selected:
        roots = tuple(family["roots"])
        scorable = [
            record
            for record in family["parameter_features"]
            if record["admissible"]
            and record["good_prime_coverage"] >= MINIMUM_LOCAL_COVERAGE
        ]
        scorable.sort(key=parameter_local_key)
        primary.append((roots, scorable[0]))
        if len(scorable) > 1:
            second.append((roots, scorable[1]))
    second.sort(key=lambda item: (*parameter_local_key(item[1]), item[0]))
    chosen = primary + second[:H5000_EXTRA_SECOND_FAMILY_KEEP]
    identifiers = sorted(
        "r" + "_".join(map(str, roots)) + f"_t{record['parameter']}"
        for roots, record in chosen
    )
    if len(identifiers) != 64 or len(set(identifiers)) != 64:
        raise AssertionError("the rank-blind H5000 tranche must contain 64 fibers")
    return identifiers, {
        "one_best_local_fiber_per_selected_family": 40,
        "best_second_fiber_from_distinct_families": 24,
        "predeclared_fiber_count": len(identifiers),
        "predeclared_fiber_sha256": hashlib.sha256(
            "\n".join(identifiers).encode()
        ).hexdigest(),
    }


def search_h5000(
    roots: tuple[int, ...],
    parameter: int,
    *,
    point_timeout: float,
    height_timeout: float,
) -> tuple[dict[str, Any], tuple[tuple[Fraction, Fraction], ...] | None]:
    construction = SixRootMestreConstruction(tuple(Q(root) for root in roots))
    parameter_q = Q(parameter)
    coefficients = construction.primitive_jacobian_coefficients(parameter_q)
    visible_quartic = primitive_visible_points(construction, parameter_q)
    visible_jacobian = tuple(
        quartic_point_to_jacobian(construction, parameter_q, point)
        for point in visible_quartic
    )
    raw = bounded_quartic_points(
        construction.primitive_quartic_coefficients(parameter_q),
        height_bound=H5000_HEIGHT,
        timeout=point_timeout,
        stack_bytes=STACK_BYTES,
    )
    signless = canonical_signless_points(raw)
    retained = signless[:H5000_MAPPING_CAP]
    quartic_coefficients = construction.primitive_quartic_coefficients(parameter_q)
    if any(
        point[1] ** 2 != quartic_value(quartic_coefficients, point[0])
        for point in retained
    ):
        raise AssertionError("H5000 returned a point off the exact quartic")
    searched_jacobian = tuple(
        quartic_point_to_jacobian(construction, parameter_q, point)
        for point in retained
    )
    pool_by_x = {point[0]: point for point in visible_jacobian}
    for point in searched_jacobian:
        pool_by_x.setdefault(point[0], point)
    pool = tuple(pool_by_x.values())
    height = height_matrix_replay(
        coefficients,
        pool,
        precisions=(72, 120),
        timeout=height_timeout,
        stack_bytes=STACK_BYTES,
    )
    stable_rank = int(height[-1]["numerical_rank"])
    subset = numerical_subset(pool, height)
    return (
        {
            "status": "completed exact H5000 point checks and numerical height triage",
            "height_bound": H5000_HEIGHT,
            "signed_points_returned": len(raw),
            "distinct_nonzero_ordinate_abscissas": len(signless),
            "abscissas_retained_for_mapping": len(retained),
            "mapping_cap": H5000_MAPPING_CAP,
            "mapping_truncated": len(signless) > len(retained),
            "visible_quartic_point_count": len(visible_quartic),
            "visible_jacobian_point_count": len(visible_jacobian),
            "pool_point_count_modulo_inverse": len(pool),
            "pool_point_sha256": point_digest(pool),
            "height_matrix_runs": list(height),
            "stable_numerical_rank": stable_rank,
            "numerical_subset": [point_record(point) for point in subset],
            "numerical_rank_is_not_an_independence_certificate": True,
        },
        subset,
    )


def result_digest(artifact: dict[str, Any]) -> str:
    conductors = []
    for record in artifact["specialization_screen"]["conductor_records"]:
        phase = record["conductor_phase"]
        conductors.append(
            [
                record["identifier"],
                phase["status"],
                phase.get("conductor"),
                phase.get("log_conductor"),
                phase.get("root_number"),
            ]
        )
    triage = []
    for record in artifact["specialization_screen"]["h5000_records"]:
        phase = record["point_triage"]
        triage.append(
            [
                record["identifier"],
                phase["status"],
                phase.get("stable_numerical_rank"),
                phase.get("pool_point_sha256"),
            ]
        )
    return stable_json_digest(
        {
            "census": artifact["census"],
            "selection_digest": artifact["rank_blind_selection"][
                "selected_family_sha256"
            ],
            "h5000_population_digest": artifact["rank_blind_selection"][
                "predeclared_h5000_fibers"
            ]["predeclared_fiber_sha256"],
            "conductors": conductors,
            "triage": triage,
            "target_hits": artifact["target"]["hits"],
        }
    )


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compile-timeout", type=float, default=30.0)
    parser.add_argument("--enumeration-timeout", type=float, default=30.0)
    parser.add_argument("--conductor-timeout", type=float, default=8.0)
    parser.add_argument("--point-timeout", type=float, default=12.0)
    parser.add_argument("--height-timeout", type=float, default=12.0)
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "artifacts"
        / "generated-results"
        / "elliptic_mestre_root_tuple_scale_max100.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    caps = (
        args.compile_timeout,
        args.enumeration_timeout,
        args.conductor_timeout,
        args.point_timeout,
        args.height_timeout,
    )
    if min(caps) <= 0 or max(caps) > 30:
        raise SystemExit("all subprocess caps must lie in (0,30]")
    if args.output.exists():
        raise SystemExit("refusing to overwrite the max-root-100 artifact")

    root = Path(__file__).resolve().parents[2]
    cas = root / "elliptic-curves" / "cas"
    tests = root / "elliptic-curves" / "tests"
    generated = root / "artifacts" / "generated-results"
    frozen_cpp = cas / "enumerate_mestre_root_tuples_scale.cpp"
    frozen_driver = cas / "search_mestre_root_tuple_scale.py"
    frozen_test = tests / "test_search_mestre_root_tuple_scale.py"
    frozen_artifact = generated / "elliptic_mestre_root_tuple_scale.json"
    observed_frozen = {
        "compiled_source_sha256": sha256_file(frozen_cpp),
        "driver_sha256": sha256_file(frozen_driver),
        "test_sha256": sha256_file(frozen_test),
        "artifact_sha256": sha256_file(frozen_artifact),
    }
    if observed_frozen != {
        "compiled_source_sha256": FROZEN_CPP_SHA256,
        "driver_sha256": FROZEN_DRIVER_SHA256,
        "test_sha256": FROZEN_TEST_SHA256,
        "artifact_sha256": FROZEN_ARTIFACT_SHA256,
    }:
        raise AssertionError("a frozen max-root-50 input changed")

    started = time.monotonic()
    enumeration, enumeration_timings = compiled_enumeration_max100(
        frozen_cpp,
        compile_timeout=args.compile_timeout,
        enumeration_timeout=args.enumeration_timeout,
    )
    expected_counts = (
        enumeration.normalized_count,
        enumeration.obstruction_count,
        enumeration.reflection_count,
        enumeration.nonreflection_count,
    )
    if expected_counts != EXPECTED_MAX100_COUNTS:
        raise AssertionError("the exact max-root-100 census count changed")
    verify_started = time.monotonic()
    verify_enumerator_records(enumeration)
    verification_wall_seconds = time.monotonic() - verify_started
    classify_started = time.monotonic()
    nonsingular, singular, witnesses = classify_nonreflection(enumeration)
    classification_wall_seconds = time.monotonic() - classify_started
    if (
        tuple_digest(enumeration.obstruction_roots)
        != EXPECTED_MAX100_OBSTRUCTION_SHA256
        or tuple_digest(enumeration.nonreflection_roots)
        != EXPECTED_MAX100_NONREFLECTION_SHA256
        or tuple_digest(nonsingular) != EXPECTED_MAX100_NONSINGULAR_SHA256
    ):
        raise AssertionError("the exact max-root-100 tuple digest changed")
    old_families = tuple(roots for roots in nonsingular if roots[-1] <= PRIOR_MAX_ROOT)
    new_families = tuple(roots for roots in nonsingular if roots[-1] > PRIOR_MAX_ROOT)
    if (
        len(old_families) != 44
        or tuple_digest(old_families) != EXPECTED_MAX50_NONSINGULAR_SHA256
        or len(new_families) != 191
        or tuple_digest(new_families) != EXPECTED_NEW_FAMILY_SHA256
        or set(witnesses.values()) != {1}
    ):
        raise AssertionError("the exact max-root-50 exclusion boundary changed")

    feature_started = time.monotonic()
    feature_records = [family_feature_record(roots) for roots in new_families]
    full_feature_digest = stable_json_digest(feature_records)
    selected, selection = select_family_tranche(feature_records)
    h5000_ids, h5000_selection = predeclare_h5000_fibers(selected)
    feature_wall_seconds = time.monotonic() - feature_started
    selected_roots = {tuple(record["roots"]) for record in selected}
    selected_feature = {
        tuple(record["roots"]): record for record in selected
    }

    conductor_records: list[dict[str, Any]] = []
    inadmissible: list[dict[str, Any]] = []
    runtime: dict[str, tuple[tuple[int, ...], int]] = {}
    conductor_started = time.monotonic()
    for roots in sorted(selected_roots):
        construction = SixRootMestreConstruction(tuple(Q(root) for root in roots))
        features_by_parameter = {
            record["parameter"]: record
            for record in selected_feature[roots]["parameter_features"]
        }
        for parameter in PARAMETERS:
            identifier = "r" + "_".join(map(str, roots)) + f"_t{parameter}"
            feature = features_by_parameter[parameter]
            if not feature["admissible"]:
                inadmissible.append(
                    {
                        "identifier": identifier,
                        "roots": list(roots),
                        "parameter": parameter,
                        "reason": (
                            "singular quartic or visible-point collision/zero ordinate"
                        ),
                        "collision_loss": feature["collision_loss"],
                        "zero_ordinates": feature["zero_ordinates"],
                    }
                )
                continue
            record: dict[str, Any] = {
                "identifier": identifier,
                "roots": list(roots),
                "parameter": parameter,
                "family_selection_stratum": selected_feature[roots][
                    "selection_stratum"
                ],
                "predeclared_for_h5000": identifier in h5000_ids,
                "local_score": feature["local_score"],
                "good_prime_coverage": feature["good_prime_coverage"],
            }
            try:
                coefficients = construction.primitive_jacobian_coefficients(Q(parameter))
                conductor = capped_minimal_curve_data(
                    coefficients,
                    timeout=args.conductor_timeout,
                    stack_bytes=STACK_BYTES,
                )
                record["conductor_phase"] = {
                    "status": "completed exact PARI minimal-model/conductor computation",
                    **conductor,
                    "below_strict_log_conductor_target_numerically": (
                        Decimal(conductor["log_conductor"]) < TARGET_LOG_CONDUCTOR
                    ),
                }
                runtime[identifier] = (roots, parameter)
            except CappedProcessTimeout:
                record["conductor_phase"] = {
                    "status": "timeout-no-retry",
                    "timeout_seconds": args.conductor_timeout,
                }
            except Exception as error:
                record["conductor_phase"] = {
                    "status": "error-no-retry",
                    "error": str(error)[:1000],
                }
            conductor_records.append(record)
    conductor_wall_seconds = time.monotonic() - conductor_started

    # No point or rank call occurs before the complete conductor phase above.
    h5000_records: list[dict[str, Any]] = []
    target_hits: list[dict[str, Any]] = []
    point_started = time.monotonic()
    records_by_id = {record["identifier"]: record for record in conductor_records}
    for identifier in h5000_ids:
        if identifier not in records_by_id:
            h5000_records.append(
                {
                    "identifier": identifier,
                    "point_triage": {
                        "status": "not attempted: exact admissibility gate excluded fiber"
                    },
                }
            )
            continue
        conductor_record = records_by_id[identifier]
        if identifier not in runtime:
            h5000_records.append(
                {
                    "identifier": identifier,
                    "roots": conductor_record["roots"],
                    "parameter": conductor_record["parameter"],
                    "conductor_phase": conductor_record["conductor_phase"],
                    "point_triage": {
                        "status": "not attempted after incomplete conductor"
                    },
                }
            )
            continue
        roots, parameter = runtime[identifier]
        triage_record: dict[str, Any] = {
            "identifier": identifier,
            "roots": list(roots),
            "parameter": parameter,
            "conductor_phase": conductor_record["conductor_phase"],
        }
        try:
            triage, subset = search_h5000(
                roots,
                parameter,
                point_timeout=args.point_timeout,
                height_timeout=args.height_timeout,
            )
            triage_record["point_triage"] = triage
            stable_rank = triage["stable_numerical_rank"]
            if stable_rank >= 21 and subset is not None:
                construction = SixRootMestreConstruction(
                    tuple(Q(root) for root in roots)
                )
                coefficients = construction.primitive_jacobian_coefficients(Q(parameter))
                certificate = finite_reduction_attempt(
                    coefficients, subset, prime_bound=CERTIFICATE_PRIME_BOUND
                )
                triage_record["finite_reduction_attempt"] = certificate
                certified = certificate["certified_algebraic_rank_lower_bound"]
                below_target = conductor_record["conductor_phase"][
                    "below_strict_log_conductor_target_numerically"
                ]
                if certified is not None and (
                    certified >= 30 or (certified >= 21 and below_target)
                ):
                    target_hits.append(
                        {
                            "identifier": identifier,
                            "certified_algebraic_rank_lower_bound": certified,
                            "conductor": conductor_record["conductor_phase"]["conductor"],
                            "log_conductor": conductor_record["conductor_phase"][
                                "log_conductor"
                            ],
                        }
                    )
            else:
                triage_record["finite_reduction_attempt"] = {
                    "status": "not triggered",
                    "trigger_stable_numerical_rank": 21,
                }
        except CappedProcessTimeout:
            triage_record["point_triage"] = {
                "status": "timeout-no-retry",
                "point_timeout_seconds": args.point_timeout,
                "height_timeout_seconds": args.height_timeout,
            }
        except Exception as error:
            triage_record["point_triage"] = {
                "status": "error-no-retry",
                "error": str(error)[:1000],
            }
        h5000_records.append(triage_record)
    point_wall_seconds = time.monotonic() - point_started

    completed_triage = [
        record
        for record in h5000_records
        if record["point_triage"]["status"].startswith("completed")
    ]
    script_path = Path(__file__).resolve()
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": (
            "complete exact max-root-100 census and bounded rank-blind "
            "conductor-first specialization tranche"
        ),
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": "182.72",
            "alternative_rank_at_least": 30,
            "hits": target_hits,
        },
        "frozen_max50_inputs": {
            **observed_frozen,
            "all_frozen_files_read_only": True,
        },
        "census": {
            "max_root": MAX_ROOT,
            "complete_diameter_prefix": [5, MAX_ROOT],
            "open_diameter_remainder": [],
            "affine_normalized_primitive_reflection_quotient_count": (
                enumeration.normalized_count
            ),
            "degree_five_obstruction_zero_count": enumeration.obstruction_count,
            "reflection_obstruction_zero_count": enumeration.reflection_count,
            "nonreflection_obstruction_zero_count": enumeration.nonreflection_count,
            "nonreflection_generically_nonsingular_count": len(nonsingular),
            "nonreflection_generically_singular_count": len(singular),
            "obstruction_tuple_sha256": tuple_digest(enumeration.obstruction_roots),
            "nonreflection_tuple_sha256": tuple_digest(
                enumeration.nonreflection_roots
            ),
            "nonsingular_nonreflection_tuple_sha256": tuple_digest(nonsingular),
            "max50_nonsingular_family_count": len(old_families),
            "max50_nonsingular_family_sha256": tuple_digest(old_families),
            "genuinely_new_diameter_51_to_100_family_count": len(new_families),
            "genuinely_new_diameter_51_to_100_family_sha256": tuple_digest(
                new_families
            ),
            "all_generic_nonsingularity_witness_parameters": sorted(
                set(witnesses.values())
            ),
            "exact_obstruction_replayed_in_python": True,
            "exact_normalization_and_reflection_gates_replayed_in_python": True,
            "generic_nonsingularity_test_parameter_count": 21,
        },
        "rank_blind_selection": {
            "selection_uses_conductor": False,
            "selection_uses_point_search": False,
            "selection_uses_numerical_or_algebraic_rank": False,
            "fixed_local_primes": list(LOCAL_PRIMES),
            "minimum_local_prime_coverage": MINIMUM_LOCAL_COVERAGE,
            "local_score_formula": (
                "sum over good fixed-panel primes of (2-a_p)/(p+1-a_p)"
            ),
            "geometry_proxy": (
                "descending Vandermonde/radical(Vandermonde), then ascending "
                "radical and primitive-Jacobian coefficient height"
            ),
            "full_new_family_feature_sha256": full_feature_digest,
            **selection,
            "selected_families": selected,
            "predeclared_h5000_fibers": h5000_selection,
            "predeclared_h5000_identifiers": h5000_ids,
        },
        "specialization_screen": {
            "protocol": {
                "integer_parameters": [1, 8],
                "conductor_population_closed_before_any_point_or_rank_call": True,
                "all_admissible_selected_family_fibers_receive_conductor": True,
                "h5000_population_predeclared_before_conductor_phase": True,
                "h5000_searches_run_regardless_of_conductor_target_when_conductor_completed": True,
                "no_retries": True,
                "finite_reduction_trigger_stable_numerical_rank": 21,
            },
            "population": {
                "new_family_universe": len(new_families),
                "selected_family_count": len(selected),
                "proposed_integer_fibers": len(selected) * len(PARAMETERS),
                "inadmissible_fibers": len(inadmissible),
                "admissible_conductor_records": len(conductor_records),
                "conductor_completed": sum(
                    record["conductor_phase"]["status"].startswith("completed")
                    for record in conductor_records
                ),
                "conductor_timeouts": sum(
                    record["conductor_phase"]["status"] == "timeout-no-retry"
                    for record in conductor_records
                ),
                "conductor_errors": sum(
                    record["conductor_phase"]["status"] == "error-no-retry"
                    for record in conductor_records
                ),
                "subtarget_conductors": sum(
                    record["conductor_phase"].get(
                        "below_strict_log_conductor_target_numerically"
                    )
                    is True
                    for record in conductor_records
                ),
                "predeclared_h5000_fibers": len(h5000_ids),
                "h5000_completed": len(completed_triage),
                "h5000_timeouts": sum(
                    record["point_triage"]["status"] == "timeout-no-retry"
                    for record in h5000_records
                ),
                "h5000_errors": sum(
                    record["point_triage"]["status"] == "error-no-retry"
                    for record in h5000_records
                ),
                "maximum_stable_numerical_rank": max(
                    (
                        record["point_triage"]["stable_numerical_rank"]
                        for record in completed_triage
                    ),
                    default=None,
                ),
                "stable_numerical_rank_histogram": dict(
                    sorted(
                        Counter(
                            str(record["point_triage"]["stable_numerical_rank"])
                            for record in completed_triage
                        ).items()
                    )
                ),
            },
            "inadmissible_fibers": inadmissible,
            "conductor_records": conductor_records,
            "h5000_records": h5000_records,
        },
        "parameters": {
            "compile_timeout_seconds": args.compile_timeout,
            "enumeration_timeout_seconds": args.enumeration_timeout,
            "conductor_timeout_seconds": args.conductor_timeout,
            "point_timeout_seconds": args.point_timeout,
            "height_timeout_seconds": args.height_timeout,
            "stack_bytes": STACK_BYTES,
            "H5000_height_bound": H5000_HEIGHT,
            "H5000_mapping_cap": H5000_MAPPING_CAP,
            "certificate_prime_bound": CERTIFICATE_PRIME_BOUND,
        },
        "timings": {
            **enumeration_timings,
            "independent_enumerator_replay_wall_seconds": verification_wall_seconds,
            "generic_nonsingularity_classification_wall_seconds": (
                classification_wall_seconds
            ),
            "rank_blind_feature_and_selection_wall_seconds": feature_wall_seconds,
            "conductor_phase_wall_seconds": conductor_wall_seconds,
            "H5000_phase_wall_seconds": point_wall_seconds,
            "total_wall_seconds": time.monotonic() - started,
        },
        "provenance": {
            "script": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "reproducing_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "subprocesses_run_in_foreground_process_groups": True,
            "whole_process_group_killed_and_reaped_on_timeout": True,
            "temporary_enumerator_binary_removed": True,
            "owned_processes_remaining": 0,
        },
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    artifact["result_sha256"] = result_digest(artifact)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    descriptor = os.open(args.output, flags, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(
        json.dumps(
            {
                "new_families": len(new_families),
                "selected_families": len(selected),
                "conductors": len(conductor_records),
                "h5000_completed": len(completed_triage),
                "maximum_stable_numerical_rank": artifact[
                    "specialization_screen"
                ]["population"]["maximum_stable_numerical_rank"],
                "target_hits": target_hits,
                "result_sha256": artifact["result_sha256"],
            },
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
