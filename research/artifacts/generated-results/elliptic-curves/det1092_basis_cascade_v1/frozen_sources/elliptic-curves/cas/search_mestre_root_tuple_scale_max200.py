#!/usr/bin/env python3
"""Exact max-root-200 Mestre census and rank-aware bounded fiber tranche.

The max-root-50/max-root-100 programs and artifacts are frozen inputs.  A new
compiled enumerator exhausts primitive affine-normalized six-root tuples
through diameter 200.  Its complete max-root-100 prefix is compared record for
record with the frozen enumerator.  Every emitted obstruction-zero tuple is
then replayed by a separate Python integer formula, and all nonreflection
families receive the exact 21-value generic-singularity classification.

The calculation has two explicit stages.  ``--stage census`` writes the exact
census checkpoint before any specialization arithmetic.  ``--stage screen``
loads that immutable checkpoint and screens T=1,...,8 in every genuinely new
nonsingular family.  The screen maps all visible points exactly and computes
exact mod-3 finite-reduction independence signatures; point count, root
geometry, or radical size alone never selects a leader.  Sixty-four distinct
rank-aware/diversity family leaders then receive conductors and one bounded
H=5000 point search.  A stable numerical rank of at least fourteen triggers an
immediate exact finite-reduction certificate attempt.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import comb, gcd, isqrt
import os
from pathlib import Path
import platform
import shlex
import shutil
import sys
import tempfile
import time
from typing import Any, Iterable, Sequence

from ek_k3 import rational_to_string
from mestre_root_tuples import (
    SixRootMestreConstruction,
    normalize_integer_root_tuple,
)
from mod2_reduction_independence import (
    _primes_up_to,
    _reduce_rational,
    combined_mod2_rank,
    finite_add,
    finite_curve_points,
    finite_multiply,
    finite_subtract,
)
from search_mestre_root_tuple_scale import (
    CappedProcessTimeout,
    EnumerationResult,
    TARGET_LOG_CONDUCTOR,
    capped_minimal_curve_data,
    finite_reduction_attempt,
    point_digest,
    point_record,
    run_capped_process,
    sha256_file,
    tuple_digest,
)
from search_mestre_root_tuple_scale_max100 import (
    EXPECTED_MAX100_COUNTS,
    EXPECTED_MAX100_NONREFLECTION_SHA256,
    EXPECTED_MAX100_NONSINGULAR_SHA256,
    EXPECTED_MAX100_OBSTRUCTION_SHA256,
    coefficient_height,
    compiled_enumeration_max100,
    curve_discriminant,
    fraction_record,
    local_trace,
    search_h5000,
    stable_json_digest,
)


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

MAX_ROOT = 200
PRIOR_MAX_ROOT = 100
PARAMETERS = tuple(range(1, 9))
MODULUS = 3
MOD3_PRIME_BOUND = 251
LOCAL_PRIMES = (11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
GLOBAL_RANK_AWARE_FAMILY_KEEP = 34
DIVERSITY_KEEP_PER_DECILE = 3
SELECTED_FAMILY_COUNT = 64
H5000_HEIGHT = 5_000
STACK_BYTES = 256_000_000
STRONG_GAIN_TRIGGER = 14
EXACT_GAIN_PRIME_BOUND = 499

FROZEN_MAX100_CPP_SHA256 = (
    "31650333800698201819eddc91bf228089824bca026c629c9360683324a69eb5"
)
FROZEN_MAX100_DRIVER_SHA256 = (
    "34677c38be30aa15e99b3239a6d487a51c158fa33326826d37ceead310555600"
)
FROZEN_MAX100_TEST_SHA256 = (
    "e29826713ab9ba3045701891c1ae5bf07f9c1c69a86c171b11c86bc3c576f965"
)
FROZEN_MAX100_ARTIFACT_SHA256 = (
    "63dcd39555ad8b39c7b584a16663164bf73e6c6c59906b6a230bfa9b9f65a3bb"
)
EXPECTED_MAX100_NONSINGULAR_COUNT = 235


@dataclass(frozen=True)
class ModLReductionSignature:
    prime: int
    group_order: int
    multiple_subgroup_order: int
    quotient_dimension: int
    rows: tuple[tuple[int, ...], ...]


def parse_enumerator_output(stdout: str, max_root: int) -> EnumerationResult:
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not lines or lines[0] != "MESTRE_ROOT_TUPLES_V2":
        raise AssertionError("the max-root-200 enumerator omitted its V2 header")
    roots: list[tuple[int, ...]] = []
    reflection: list[tuple[int, ...]] = []
    nonreflection: list[tuple[int, ...]] = []
    summary: tuple[int, ...] | None = None
    for line in lines[1:]:
        fields = line.split()
        if fields[0] == "R":
            if len(fields) != 8:
                raise AssertionError("malformed max-root-200 tuple record")
            item = tuple(int(value) for value in fields[1:7])
            flag = int(fields[7])
            if flag not in (0, 1):
                raise AssertionError("malformed reflection flag")
            roots.append(item)
            (reflection if flag else nonreflection).append(item)
        elif fields[0] == "S":
            if len(fields) != 6 or summary is not None:
                raise AssertionError("malformed max-root-200 summary")
            summary = tuple(int(value) for value in fields[1:])
        else:
            raise AssertionError("unknown max-root-200 output record")
    if summary is None:
        raise AssertionError("the max-root-200 enumerator omitted its summary")
    declared, normalized, obstruction, reflected, nonreflected = summary
    if (
        declared != max_root
        or obstruction != len(roots)
        or reflected != len(reflection)
        or nonreflected != len(nonreflection)
        or reflected + nonreflected != obstruction
    ):
        raise AssertionError("the max-root-200 stream disagrees with its summary")
    if roots != sorted(roots, key=lambda item: (item[-1], item)):
        raise AssertionError("the max-root-200 tuple order changed")
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


def compiled_enumeration_max200(
    source: Path,
    max_root: int,
    *,
    compile_timeout: float,
    enumeration_timeout: float,
) -> tuple[EnumerationResult, str, dict[str, float]]:
    if max_root < 5 or max_root > MAX_ROOT:
        raise ValueError("max_root must lie in [5,200]")
    compiler = shutil.which("c++")
    if compiler is None:
        raise FileNotFoundError("a C++17 compiler is required")
    timings: dict[str, float] = {}
    with tempfile.TemporaryDirectory(prefix="mestre-root-200-") as directory:
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
            (str(binary), str(max_root)), timeout=enumeration_timeout
        )
        timings["enumeration_wall_seconds"] = time.monotonic() - started
    return (
        parse_enumerator_output(stdout, max_root),
        hashlib.sha256(stdout.encode()).hexdigest(),
        timings,
    )


def mestre_obstruction_integer(roots: Sequence[int]) -> int:
    symmetric = [1, 0, 0, 0, 0, 0]
    for root in roots:
        for degree in range(5, 0, -1):
            symmetric[degree] += root * symmetric[degree - 1]
    s1, s2, s3, s4, s5 = symmetric[1:]
    return (
        -s1**5
        + 6 * s1**3 * s2
        - 7 * s1**2 * s3
        - 8 * s1 * s2**2
        + 8 * s1 * s4
        + 12 * s2 * s3
        - 24 * s5
    )


def reflection_symmetric(roots: Sequence[int]) -> bool:
    return all(roots[index] + roots[5 - index] == roots[-1] for index in range(3))


def verify_enumerator_records_fast(enumeration: EnumerationResult) -> None:
    """Independently replay every emitted V2 tuple in expected linear time."""

    reflected = set(enumeration.reflection_roots)
    if len(reflected) != enumeration.reflection_count:
        raise AssertionError("the reflection stream contains a duplicate")
    seen: set[tuple[int, ...]] = set()
    for roots in enumeration.obstruction_roots:
        if roots in seen:
            raise AssertionError("the V2 enumerator emitted a duplicate")
        seen.add(roots)
        if (
            len(roots) != 6
            or roots[0] != 0
            or any(roots[index] >= roots[index + 1] for index in range(5))
            or roots[-1] > enumeration.max_root
        ):
            raise AssertionError("the V2 enumerator emitted invalid roots")
        if gcd(*roots[1:]) != 1 or normalize_integer_root_tuple(roots) != roots:
            raise AssertionError("the V2 normalization gate changed")
        if mestre_obstruction_integer(roots) != 0:
            raise AssertionError("the V2 obstruction replay failed")
        if reflection_symmetric(roots) != (roots in reflected):
            raise AssertionError("the V2 reflection gate changed")


def mobius(value: int) -> int:
    result = 1
    remaining = value
    prime = 2
    while prime * prime <= remaining:
        if remaining % prime == 0:
            remaining //= prime
            result = -result
            if remaining % prime == 0:
                return 0
            while remaining % prime == 0:
                remaining //= prime
        prime += 1 if prime == 2 else 2
    if remaining > 1:
        result = -result
    return result


def independent_normalized_count(max_root: int) -> dict[str, int]:
    """Count primitive reflection orbits without traversing the five-loop box."""

    primitive = sum(
        mobius(divisor) * comb(max_root // divisor, 5)
        for divisor in range(1, max_root + 1)
        if max_root // divisor >= 5
    )
    symmetric = 0
    for diameter in range(5, max_root + 1):
        for left in range(1, (diameter - 1) // 2 + 1):
            for middle in range(left + 1, (diameter - 1) // 2 + 1):
                if gcd(left, middle, diameter) == 1:
                    symmetric += 1
    if (primitive + symmetric) % 2:
        raise AssertionError("Burnside orbit count is not integral")
    return {
        "primitive_unquotiented_count": primitive,
        "primitive_reflection_fixed_count": symmetric,
        "primitive_reflection_orbit_count": (primitive + symmetric) // 2,
    }


def generic_classification(
    roots_population: Sequence[tuple[int, ...]],
) -> tuple[
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, ...], ...],
    dict[tuple[int, ...], int],
]:
    nonsingular: list[tuple[int, ...]] = []
    singular: list[tuple[int, ...]] = []
    witnesses: dict[tuple[int, ...], int] = {}
    for roots in roots_population:
        construction = SixRootMestreConstruction(tuple(Q(root) for root in roots))
        witness = next(
            (
                parameter
                for parameter in range(1, 22)
                if construction.quartic_discriminant(Q(parameter)) != 0
            ),
            None,
        )
        if witness is None:
            singular.append(roots)
        else:
            nonsingular.append(roots)
            witnesses[roots] = witness
    return tuple(nonsingular), tuple(singular), witnesses


def gf_l_rank_and_pivots(
    rows: Iterable[Sequence[int]], column_count: int, modulus: int
) -> tuple[int, tuple[int, ...]]:
    if modulus < 2 or any(modulus % divisor == 0 for divisor in range(2, isqrt(modulus) + 1)):
        raise ValueError("the row modulus must be prime")
    matrix = []
    for row in rows:
        item = [int(value) % modulus for value in row]
        if len(item) != column_count:
            raise ValueError("a finite-reduction row has the wrong width")
        if any(item):
            matrix.append(item)
    rank = 0
    pivots: list[int] = []
    for column in range(column_count):
        pivot = next(
            (
                index
                for index in range(rank, len(matrix))
                if matrix[index][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, modulus)
        matrix[rank] = [(value * inverse) % modulus for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank or matrix[index][column] == 0:
                continue
            multiplier = matrix[index][column]
            matrix[index] = [
                (left - multiplier * right) % modulus
                for left, right in zip(matrix[index], matrix[rank])
            ]
        pivots.append(column)
        rank += 1
        if rank == len(matrix):
            break
    return rank, tuple(pivots)


def mod_l_reduction_signature(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    prime: int,
    *,
    modulus: int = MODULUS,
) -> ModLReductionSignature:
    """Return exact rows for E(F_p)/modulus*E(F_p)."""

    if prime == modulus:
        raise ValueError("the reduction prime must differ from the quotient modulus")
    if len(coefficients) != 5 or any(Q(value) for value in coefficients[:3]):
        raise ValueError("the certificate requires a short Weierstrass model")
    coefficient_a = _reduce_rational(Q(coefficients[3]), prime)
    coefficient_b = _reduce_rational(Q(coefficients[4]), prime)
    discriminant = -16 * (4 * coefficient_a**3 + 27 * coefficient_b**2)
    if discriminant % prime == 0:
        raise ValueError("the curve has bad reduction")
    finite_points = finite_curve_points(coefficient_a, coefficient_b, prime)
    multiples = {
        finite_multiply(point, modulus, coefficient_a, prime)
        for point in finite_points
    }
    span = [None]
    coordinates: list[tuple[int, ...]] = [()]
    basis = []
    for point in finite_points:
        if any(
            finite_subtract(point, representative, coefficient_a, prime) in multiples
            for representative in span
        ):
            continue
        basis.append(point)
        old_span = tuple(span)
        old_coordinates = tuple(coordinates)
        span = []
        coordinates = []
        for scalar in range(modulus):
            multiple = finite_multiply(point, scalar, coefficient_a, prime)
            span.extend(
                finite_add(representative, multiple, coefficient_a, prime)
                for representative in old_span
            )
            coordinates.extend(
                coordinate + (scalar,) for coordinate in old_coordinates
            )
    if len(span) * len(multiples) != len(finite_points):
        raise AssertionError("finite quotient representatives do not cover the group")
    quotient_order = len(span)
    quotient_dimension = 0
    while quotient_order > 1 and quotient_order % modulus == 0:
        quotient_order //= modulus
        quotient_dimension += 1
    if quotient_order != 1 or quotient_dimension != len(basis):
        raise AssertionError("the finite quotient is not the expected vector space")
    rows = [[0] * len(points) for _ in basis]
    for point_index, point in enumerate(points):
        reduced = (
            _reduce_rational(Q(point[0]), prime),
            _reduce_rational(Q(point[1]), prime),
        )
        coordinate_index = next(
            (
                index
                for index, representative in enumerate(span)
                if finite_subtract(reduced, representative, coefficient_a, prime)
                in multiples
            ),
            None,
        )
        if coordinate_index is None:
            raise AssertionError("a rational point missed every finite quotient coset")
        for basis_index, value in enumerate(coordinates[coordinate_index]):
            rows[basis_index][point_index] = value
    return ModLReductionSignature(
        prime=prime,
        group_order=len(finite_points),
        multiple_subgroup_order=len(multiples),
        quotient_dimension=quotient_dimension,
        rows=tuple(tuple(row) for row in rows),
    )


def mod3_independence_certificate(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    prime_bound: int,
) -> dict[str, Any]:
    if not points:
        raise ValueError("the mod-3 certificate requires at least one point")
    selected: list[ModLReductionSignature] = []
    rows: list[tuple[int, ...]] = []
    current_rank = 0
    torsion_certificate: dict[str, int] | None = None
    for prime in _primes_up_to(prime_bound):
        if prime in (2, MODULUS):
            continue
        try:
            signature = mod_l_reduction_signature(
                coefficients, points, prime, modulus=MODULUS
            )
        except ValueError:
            continue
        if torsion_certificate is None and signature.group_order % MODULUS:
            torsion_certificate = {
                "prime": prime,
                "group_order": signature.group_order,
            }
        candidate_rows = [*rows, *signature.rows]
        candidate_rank, _ = gf_l_rank_and_pivots(
            candidate_rows, len(points), MODULUS
        )
        if candidate_rank > current_rank:
            selected.append(signature)
            rows.extend(signature.rows)
            current_rank = candidate_rank
        if current_rank == len(points) and torsion_certificate is not None:
            break
    if torsion_certificate is None:
        raise AssertionError("no good reduction excluded rational 3-torsion")
    rank, pivots = gf_l_rank_and_pivots(rows, len(points), MODULUS)
    if rank != current_rank or len(pivots) != rank:
        raise AssertionError("the combined mod-3 rank replay changed")
    subset = tuple(points[index] for index in pivots)
    signatures = [asdict(signature) for signature in selected]
    return {
        "status": "certified exact algebraic rank lower bound",
        "descent_modulus": MODULUS,
        "certificate_prime_bound": prime_bound,
        "point_count": len(points),
        "point_sha256": point_digest(points),
        "certificate_primes": [signature.prime for signature in selected],
        "combined_exact_rank_over_F3": rank,
        "independent_subset_indices_one_based": [index + 1 for index in pivots],
        "independent_subset_sha256": point_digest(subset),
        "rational_3_torsion_exclusion": {
            **torsion_certificate,
            "reason": "rational prime-to-p torsion injects at good reduction",
        },
        "signatures": signatures,
        "certified_algebraic_rank_lower_bound": rank,
    }


def visible_points_and_coefficients(
    roots: tuple[int, ...], parameter: int
) -> tuple[
    SixRootMestreConstruction,
    tuple[Fraction, ...],
    tuple[tuple[Fraction, Fraction], ...],
]:
    from search_mestre_root_tuple_scale import (
        primitive_visible_points,
        quartic_point_to_jacobian,
    )

    construction = SixRootMestreConstruction(tuple(Q(root) for root in roots))
    parameter_q = Q(parameter)
    coefficients = construction.primitive_jacobian_coefficients(parameter_q)
    quartic_points = primitive_visible_points(construction, parameter_q)
    jacobian_points = tuple(
        quartic_point_to_jacobian(construction, parameter_q, point)
        for point in quartic_points
    )
    return construction, coefficients, jacobian_points


def screen_family(roots: tuple[int, ...]) -> dict[str, Any]:
    construction = SixRootMestreConstruction(tuple(Q(root) for root in roots))
    fibers = []
    for parameter in PARAMETERS:
        parameter_q = Q(parameter)
        degeneracy = construction.visible_point_degeneracy(parameter_q)
        admissible = (
            construction.quartic_discriminant(parameter_q) != 0
            and degeneracy.collision_loss == 0
            and degeneracy.zero_ordinates == 0
        )
        identifier = "r" + "_".join(map(str, roots)) + f"_t{parameter}"
        record: dict[str, Any] = {
            "identifier": identifier,
            "parameter": parameter,
            "admissible": admissible,
            "collision_loss": degeneracy.collision_loss,
            "zero_ordinates": degeneracy.zero_ordinates,
        }
        if admissible:
            _, coefficients, visible_points = visible_points_and_coefficients(
                roots, parameter
            )
            certificate = mod3_independence_certificate(
                coefficients, visible_points, prime_bound=MOD3_PRIME_BOUND
            )
            traces = []
            local_score = Q(0)
            for prime in LOCAL_PRIMES:
                trace = local_trace(coefficients, prime)
                if trace is None:
                    continue
                traces.append((prime, trace))
                local_score += Q(2 - trace, prime + 1 - trace)
            record.update(
                {
                    "coefficient_height": coefficient_height(coefficients),
                    "exact_visible_jacobian_point_count": len(visible_points),
                    "exact_visible_jacobian_point_sha256": point_digest(visible_points),
                    "mod3_finite_reduction_certificate": certificate,
                    "good_local_prime_coverage": len(traces),
                    "local_score": fraction_record(local_score),
                    "local_traces": [
                        {"prime": prime, "trace": trace}
                        for prime, trace in traces
                    ],
                }
            )
        fibers.append(record)
    return {
        "roots": list(roots),
        "diameter": roots[-1],
        "diameter_decile": (
            f"{10 * ((roots[-1] - 1) // 10) + 1}-"
            f"{10 * ((roots[-1] - 1) // 10) + 10}"
        ),
        "fibers": fibers,
    }


def fiber_rank_key(record: dict[str, Any]) -> tuple[Any, ...]:
    certificate = record["mod3_finite_reduction_certificate"]
    return (
        -certificate["certified_algebraic_rank_lower_bound"],
        -Q(record["local_score"]["value"]),
        -record["good_local_prime_coverage"],
        record["coefficient_height"],
        record["parameter"],
        record["identifier"],
    )


def select_rank_aware_diversity_leaders(
    family_records: Sequence[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    best_by_family = []
    for family in family_records:
        admissible = [record for record in family["fibers"] if record["admissible"]]
        if not admissible:
            continue
        best = dict(sorted(admissible, key=fiber_rank_key)[0])
        best["roots"] = family["roots"]
        best["diameter"] = family["diameter"]
        best["diameter_decile"] = family["diameter_decile"]
        best_by_family.append(best)
    if len(best_by_family) < SELECTED_FAMILY_COUNT:
        raise AssertionError("fewer than 64 new families have an admissible panel fiber")
    selected: dict[tuple[int, ...], dict[str, Any]] = {}
    for record in sorted(best_by_family, key=fiber_rank_key)[
        :GLOBAL_RANK_AWARE_FAMILY_KEEP
    ]:
        item = dict(record)
        item["selection_stratum"] = "top-34 exact mod-3 rank/local-score leader"
        selected[tuple(item["roots"])] = item
    diversity_counts = {}
    for lower in range(101, 201, 10):
        upper = lower + 9
        pool = [
            record
            for record in best_by_family
            if lower <= record["diameter"] <= upper
            and tuple(record["roots"]) not in selected
        ]
        chosen = sorted(pool, key=fiber_rank_key)[:DIVERSITY_KEEP_PER_DECILE]
        if len(chosen) != DIVERSITY_KEEP_PER_DECILE:
            raise AssertionError(f"diameter decile {lower}-{upper} lacks diversity leaders")
        diversity_counts[f"{lower}-{upper}"] = len(chosen)
        for record in chosen:
            item = dict(record)
            item["selection_stratum"] = f"rank-aware diversity-{lower}-{upper}"
            selected[tuple(item["roots"])] = item
    result = sorted(selected.values(), key=lambda record: record["identifier"])
    if len(result) != SELECTED_FAMILY_COUNT:
        raise AssertionError("the rank-aware/diversity tranche must contain 64 families")
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


def frozen_max100_inputs(root: Path) -> tuple[dict[str, str], dict[str, Path]]:
    cas = root / "elliptic-curves" / "cas"
    tests = root / "elliptic-curves" / "tests"
    generated = root / "artifacts" / "generated-results"
    paths = {
        "compiled_source": cas / "enumerate_mestre_root_tuples_scale.cpp",
        "driver": cas / "search_mestre_root_tuple_scale_max100.py",
        "test": tests / "test_search_mestre_root_tuple_scale_max100.py",
        "artifact": generated / "elliptic_mestre_root_tuple_scale_max100.json",
    }
    observed = {f"{name}_sha256": sha256_file(path) for name, path in paths.items()}
    expected = {
        "compiled_source_sha256": FROZEN_MAX100_CPP_SHA256,
        "driver_sha256": FROZEN_MAX100_DRIVER_SHA256,
        "test_sha256": FROZEN_MAX100_TEST_SHA256,
        "artifact_sha256": FROZEN_MAX100_ARTIFACT_SHA256,
    }
    if observed != expected:
        raise AssertionError("a frozen max-root-100 input changed")
    return observed, paths


def exclusive_json_write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def build_census(args: argparse.Namespace, root: Path) -> dict[str, Any]:
    if args.census_output.exists():
        raise SystemExit("refusing to overwrite the max-root-200 census checkpoint")
    frozen, frozen_paths = frozen_max100_inputs(root)
    source = root / "elliptic-curves" / "cas" / "enumerate_mestre_root_tuples_scale_max200.cpp"
    started = time.monotonic()
    prefix, prefix_stream_sha256, prefix_timings = compiled_enumeration_max200(
        source,
        PRIOR_MAX_ROOT,
        compile_timeout=args.compile_timeout,
        enumeration_timeout=args.enumeration_timeout,
    )
    frozen_prefix, _ = compiled_enumeration_max100(
        frozen_paths["compiled_source"],
        compile_timeout=args.compile_timeout,
        enumeration_timeout=args.enumeration_timeout,
    )
    if prefix != frozen_prefix:
        raise AssertionError("the V2 enumerator failed exact max-root-100 recovery")
    full, full_stream_sha256, full_timings = compiled_enumeration_max200(
        source,
        MAX_ROOT,
        compile_timeout=args.compile_timeout,
        enumeration_timeout=args.enumeration_timeout,
    )
    verify_started = time.monotonic()
    verify_enumerator_records_fast(full)
    independent_count = independent_normalized_count(MAX_ROOT)
    if independent_count["primitive_reflection_orbit_count"] != full.normalized_count:
        raise AssertionError("the independent Burnside/Mobius count disagrees")
    verification_wall_seconds = time.monotonic() - verify_started
    classification_started = time.monotonic()
    nonsingular, singular, witnesses = generic_classification(full.nonreflection_roots)
    classification_wall_seconds = time.monotonic() - classification_started
    old_nonsingular = tuple(roots for roots in nonsingular if roots[-1] <= PRIOR_MAX_ROOT)
    new_nonsingular = tuple(roots for roots in nonsingular if roots[-1] > PRIOR_MAX_ROOT)
    if (
        len(old_nonsingular) != EXPECTED_MAX100_NONSINGULAR_COUNT
        or tuple_digest(old_nonsingular) != EXPECTED_MAX100_NONSINGULAR_SHA256
    ):
        raise AssertionError("the exact max-root-100 nonsingular prefix changed")
    if set(witnesses.values()) != {1}:
        raise AssertionError("the nonsingularity witness population changed")
    script_path = Path(__file__).resolve()
    checkpoint: dict[str, Any] = {
        "schema_version": 1,
        "status": "complete exact max-root-200 census checkpoint; no specialization calls",
        "scope": {
            "max_root": MAX_ROOT,
            "complete_diameter_prefix": [5, MAX_ROOT],
            "open_diameter_remainder": [],
            "specialization_arithmetic_run": False,
        },
        "frozen_max100_inputs": {**frozen, "all_frozen_files_read_only": True},
        "exact_max100_prefix_recovery": {
            "record_for_record_equal_to_frozen_enumerator": True,
            "counts": [
                prefix.normalized_count,
                prefix.obstruction_count,
                prefix.reflection_count,
                prefix.nonreflection_count,
            ],
            "expected_counts": list(EXPECTED_MAX100_COUNTS),
            "obstruction_tuple_sha256": tuple_digest(prefix.obstruction_roots),
            "expected_obstruction_tuple_sha256": EXPECTED_MAX100_OBSTRUCTION_SHA256,
            "nonreflection_tuple_sha256": tuple_digest(prefix.nonreflection_roots),
            "expected_nonreflection_tuple_sha256": (
                EXPECTED_MAX100_NONREFLECTION_SHA256
            ),
            "new_enumerator_stream_sha256": prefix_stream_sha256,
        },
        "census": {
            "affine_normalized_primitive_reflection_quotient_count": full.normalized_count,
            "independent_burnside_mobius_count": independent_count,
            "degree_five_obstruction_zero_count": full.obstruction_count,
            "reflection_obstruction_zero_count": full.reflection_count,
            "nonreflection_obstruction_zero_count": full.nonreflection_count,
            "nonreflection_generically_nonsingular_count": len(nonsingular),
            "nonreflection_generically_singular_count": len(singular),
            "obstruction_tuple_sha256": tuple_digest(full.obstruction_roots),
            "reflection_tuple_sha256": tuple_digest(full.reflection_roots),
            "nonreflection_tuple_sha256": tuple_digest(full.nonreflection_roots),
            "nonsingular_nonreflection_tuple_sha256": tuple_digest(nonsingular),
            "singular_nonreflection_tuple_sha256": tuple_digest(singular),
            "max100_nonsingular_family_count": len(old_nonsingular),
            "max100_nonsingular_family_sha256": tuple_digest(old_nonsingular),
            "genuinely_new_diameter_101_to_200_family_count": len(new_nonsingular),
            "genuinely_new_diameter_101_to_200_family_sha256": tuple_digest(
                new_nonsingular
            ),
            "all_obstruction_records_replayed_by_independent_python_integer_formula": True,
            "all_normalization_and_reflection_gates_replayed_in_python": True,
            "generic_singularity_test_parameter_count": 21,
            "generic_discriminant_degree_upper_bound": 20,
            "all_nonsingularity_witness_parameters": sorted(set(witnesses.values())),
            "raw_max200_stream_sha256": full_stream_sha256,
        },
        "tuple_populations": {
            "nonreflection_roots": [list(roots) for roots in full.nonreflection_roots],
            "generically_nonsingular_nonreflection_roots": [
                list(roots) for roots in nonsingular
            ],
            "generically_singular_nonreflection_roots": [
                list(roots) for roots in singular
            ],
            "genuinely_new_nonsingular_roots": [
                list(roots) for roots in new_nonsingular
            ],
        },
        "arithmetic_bounds": {
            "signed_integer_type_bits": 128,
            "absolute_obstruction_expression_bound": str(47_520 * MAX_ROOT**5),
            "signed_128_max": str(2**127 - 1),
            "bound_is_strictly_below_signed_128_max": 47_520 * MAX_ROOT**5 < 2**127,
        },
        "parameters": {
            "compile_timeout_seconds": args.compile_timeout,
            "enumeration_timeout_seconds": args.enumeration_timeout,
        },
        "timings": {
            "max100_prefix": prefix_timings,
            "max200_full": full_timings,
            "independent_integer_replay_wall_seconds": verification_wall_seconds,
            "generic_classification_wall_seconds": classification_wall_seconds,
            "total_wall_seconds": time.monotonic() - started,
        },
        "provenance": {
            "script": str(script_path.relative_to(root)),
            "script_sha256": sha256_file(script_path),
            "compiled_source": str(source.relative_to(root)),
            "compiled_source_sha256": sha256_file(source),
            "reproducing_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "subprocesses_run_in_foreground_process_groups": True,
            "whole_process_group_killed_and_reaped_on_timeout": True,
            "temporary_enumerator_binaries_removed": True,
            "owned_processes_remaining": 0,
        },
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    checkpoint["result_sha256"] = stable_json_digest(
        {
            "prefix": checkpoint["exact_max100_prefix_recovery"],
            "census": checkpoint["census"],
            "new_roots": checkpoint["tuple_populations"][
                "genuinely_new_nonsingular_roots"
            ],
        }
    )
    exclusive_json_write(args.census_output, checkpoint)
    return checkpoint


def screen_result_digest(artifact: dict[str, Any]) -> str:
    conductors = [
        [
            record["identifier"],
            record["conductor_phase"]["status"],
            record["conductor_phase"].get("conductor"),
            record["conductor_phase"].get("log_conductor"),
            record["conductor_phase"].get("root_number"),
        ]
        for record in artifact["leader_followup"]["records"]
    ]
    points = [
        [
            record["identifier"],
            record["point_triage"]["status"],
            record["point_triage"].get("stable_numerical_rank"),
            record["point_triage"].get("pool_point_sha256"),
        ]
        for record in artifact["leader_followup"]["records"]
    ]
    return stable_json_digest(
        {
            "checkpoint_sha256": artifact["input"]["census_checkpoint_sha256"],
            "screen_population": artifact["complete_panel_screen"]["population"],
            "screen_digest": artifact["complete_panel_screen"][
                "family_records_sha256"
            ],
            "selection": artifact["rank_aware_diversity_selection"][
                "selected_identifier_sha256"
            ],
            "conductors": conductors,
            "points": points,
            "target_hits": artifact["target"]["hits"],
        }
    )


def build_screen(args: argparse.Namespace, root: Path) -> dict[str, Any]:
    if args.output.exists():
        raise SystemExit("refusing to overwrite the max-root-200 screen artifact")
    frozen, _ = frozen_max100_inputs(root)
    checkpoint = json.loads(args.census_output.read_text())
    checkpoint_sha256 = sha256_file(args.census_output)
    if (
        checkpoint["schema_version"] != 1
        or checkpoint["scope"]["max_root"] != MAX_ROOT
        or checkpoint["frozen_max100_inputs"] != {
            **frozen,
            "all_frozen_files_read_only": True,
        }
        or checkpoint["provenance"]["compiled_source_sha256"]
        != sha256_file(
            root
            / "elliptic-curves"
            / "cas"
            / "enumerate_mestre_root_tuples_scale_max200.cpp"
        )
    ):
        raise AssertionError("the max-root-200 census checkpoint is incompatible")
    new_roots = tuple(
        tuple(int(value) for value in roots)
        for roots in checkpoint["tuple_populations"]["genuinely_new_nonsingular_roots"]
    )
    census = checkpoint["census"]
    if (
        len(new_roots) != census["genuinely_new_diameter_101_to_200_family_count"]
        or tuple_digest(new_roots)
        != census["genuinely_new_diameter_101_to_200_family_sha256"]
        or any(roots[-1] <= PRIOR_MAX_ROOT or roots[-1] > MAX_ROOT for roots in new_roots)
    ):
        raise AssertionError("the checkpoint's new-family population changed")
    started = time.monotonic()
    if args.workers == 1:
        family_records = [screen_family(roots) for roots in new_roots]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            family_records = list(executor.map(screen_family, new_roots, chunksize=1))
    screen_wall_seconds = time.monotonic() - started
    if [tuple(record["roots"]) for record in family_records] != list(new_roots):
        raise AssertionError("parallel family screening changed deterministic order")
    selected, selection = select_rank_aware_diversity_leaders(family_records)
    selected_by_id = {record["identifier"]: record for record in selected}
    conductor_started = time.monotonic()
    followup_records = []
    runtime: dict[str, tuple[tuple[int, ...], int]] = {}
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
            _, coefficients, visible_points = visible_points_and_coefficients(
                roots, parameter
            )
            record["exact_visible_points"] = [point_record(point) for point in visible_points]
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
        followup_records.append(record)
    conductor_wall_seconds = time.monotonic() - conductor_started

    # The complete selected conductor population closes before any bounded
    # point search or numerical height call occurs.
    point_started = time.monotonic()
    target_hits = []
    for record in followup_records:
        identifier = record["identifier"]
        if identifier not in runtime:
            record["point_triage"] = {
                "status": "not attempted after incomplete conductor"
            }
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
            stable_rank = triage["stable_numerical_rank"]
            if stable_rank >= STRONG_GAIN_TRIGGER and subset is not None:
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
                        certified,
                        mod2["certified_algebraic_rank_lower_bound"] or 0,
                    )
                record["immediate_exact_gain_attempt"][
                    "best_certified_algebraic_rank_lower_bound"
                ] = certified
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
        except Exception as error:
            record["point_triage"] = {
                "status": "error-no-retry",
                "error": str(error)[:1000],
            }
    point_wall_seconds = time.monotonic() - point_started
    all_fibers = [fiber for family in family_records for fiber in family["fibers"]]
    admissible = [record for record in all_fibers if record["admissible"]]
    completed_points = [
        record
        for record in followup_records
        if record["point_triage"]["status"].startswith("completed")
    ]
    rank_histogram = Counter(
        str(record["mod3_finite_reduction_certificate"][
            "certified_algebraic_rank_lower_bound"
        ])
        for record in admissible
    )
    script_path = Path(__file__).resolve()
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "status": "complete max-root-200 exact panel screen and bounded leader followup",
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
            "frozen_max100_inputs": frozen,
        },
        "complete_panel_screen": {
            "protocol": {
                "integer_parameters": [1, 8],
                "every_genuinely_new_nonsingular_nonreflection_family_screened": True,
                "every_admissible_panel_fiber_maps_all_visible_points_exactly": True,
                "every_admissible_panel_fiber_receives_exact_mod3_finite_reduction_signature": True,
                "finite_reduction_prime_bound": MOD3_PRIME_BOUND,
                "finite_reduction_descent_modulus": MODULUS,
                "selection_uses_exact_finite_reduction_rank": True,
                "selection_uses_visible_count_alone": False,
                "selection_uses_root_radical_or_conductor": False,
            },
            "population": {
                "new_family_count": len(family_records),
                "proposed_panel_fiber_count": len(all_fibers),
                "admissible_panel_fiber_count": len(admissible),
                "inadmissible_panel_fiber_count": len(all_fibers) - len(admissible),
                "exact_mod3_certificate_count": len(admissible),
                "visible_rank_lower_bound_histogram": dict(sorted(rank_histogram.items())),
                "maximum_visible_certified_rank_lower_bound": max(
                    record["mod3_finite_reduction_certificate"][
                        "certified_algebraic_rank_lower_bound"
                    ]
                    for record in admissible
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
                "bounded_point_height": H5000_HEIGHT,
                "no_retries": True,
                "immediate_exact_gain_trigger_stable_numerical_rank": (
                    STRONG_GAIN_TRIGGER
                ),
            },
            "population": {
                "selected_leaders": len(followup_records),
                "conductor_completed": sum(
                    record["conductor_phase"]["status"].startswith("completed")
                    for record in followup_records
                ),
                "conductor_timeouts": sum(
                    record["conductor_phase"]["status"] == "timeout-no-retry"
                    for record in followup_records
                ),
                "conductor_errors": sum(
                    record["conductor_phase"]["status"] == "error-no-retry"
                    for record in followup_records
                ),
                "subtarget_conductors": sum(
                    record["conductor_phase"].get(
                        "below_strict_log_conductor_target_numerically"
                    )
                    is True
                    for record in followup_records
                ),
                "point_search_completed": len(completed_points),
                "point_search_timeouts": sum(
                    record["point_triage"]["status"] == "timeout-no-retry"
                    for record in followup_records
                ),
                "point_search_errors": sum(
                    record["point_triage"]["status"] == "error-no-retry"
                    for record in followup_records
                ),
                "maximum_stable_numerical_rank": max(
                    (
                        record["point_triage"]["stable_numerical_rank"]
                        for record in completed_points
                    ),
                    default=None,
                ),
                "stable_numerical_rank_histogram": dict(
                    sorted(
                        Counter(
                            str(record["point_triage"]["stable_numerical_rank"])
                            for record in completed_points
                        ).items()
                    )
                ),
                "immediate_exact_gain_attempts": sum(
                    record["immediate_exact_gain_attempt"].get("status")
                    != "not triggered"
                    for record in completed_points
                ),
            },
            "records": followup_records,
        },
        "parameters": {
            "workers": args.workers,
            "conductor_timeout_seconds": args.conductor_timeout,
            "point_timeout_seconds": args.point_timeout,
            "height_timeout_seconds": args.height_timeout,
            "stack_bytes": STACK_BYTES,
            "mod3_prime_bound": MOD3_PRIME_BOUND,
            "exact_gain_prime_bound": EXACT_GAIN_PRIME_BOUND,
        },
        "timings": {
            "complete_exact_panel_screen_wall_seconds": screen_wall_seconds,
            "conductor_phase_wall_seconds": conductor_wall_seconds,
            "point_phase_wall_seconds": point_wall_seconds,
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
            "whole_process_group_killed_and_reaped_on_timeout": True,
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
    generated = root / "artifacts" / "generated-results"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("census", "screen"), required=True)
    parser.add_argument("--compile-timeout", type=float, default=30.0)
    parser.add_argument("--enumeration-timeout", type=float, default=60.0)
    parser.add_argument("--conductor-timeout", type=float, default=8.0)
    parser.add_argument("--point-timeout", type=float, default=12.0)
    parser.add_argument("--height-timeout", type=float, default=12.0)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--census-output",
        type=Path,
        default=generated / "elliptic_mestre_root_tuple_scale_max200_census.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=generated / "elliptic_mestre_root_tuple_scale_max200.json",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 1 <= args.workers <= 16:
        raise SystemExit("workers must lie in [1,16]")
    if not 0 < args.compile_timeout <= 60 or not 0 < args.enumeration_timeout <= 120:
        raise SystemExit("compile/enumeration timeout is outside its declared cap")
    if any(
        timeout <= 0 or timeout > 30
        for timeout in (
            args.conductor_timeout,
            args.point_timeout,
            args.height_timeout,
        )
    ):
        raise SystemExit("PARI subprocess caps must lie in (0,30]")
    root = Path(__file__).resolve().parents[2]
    if args.stage == "census":
        result = build_census(args, root)
        print(
            json.dumps(
                {
                    "census": result["census"],
                    "result_sha256": result["result_sha256"],
                },
                sort_keys=True,
            ),
            flush=True,
        )
    else:
        result = build_screen(args, root)
        print(
            json.dumps(
                {
                    "panel_population": result["complete_panel_screen"]["population"],
                    "leader_population": result["leader_followup"]["population"],
                    "target_hits": result["target"]["hits"],
                    "result_sha256": result["result_sha256"],
                },
                sort_keys=True,
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
