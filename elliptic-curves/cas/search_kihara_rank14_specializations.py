#!/usr/bin/env python3
"""Bounded conductor-first specialization search in Kihara's rank-14 family.

The deterministic population has three disjoint sources: every positive
reduced rational of a fixed height, an extended integer interval, and exact
Gauss-reduced rational representatives of pairs of Hensel-lifted roots of the
degree-398 discriminant-frontier factor.  Selection for expensive conductor
work uses only exact discriminant geometry.  Point yield is first observed in
a later stage and is therefore unavailable to the conductor prefilter.

For retained specializations, PARI computes a minimal model and conductor
under a strict process-group timeout.  Sub-threshold curves (plus a declared
small conductor fallback) are searched in the global quartic chart and in
charts centered at all fifteen known sections.  All returned points and their
covariant images are checked over ``QQ`` before two-precision height ranks are
recorded.  Numerical height rank is explicitly not an exact rank certificate.
If rank 18 or more appears, the script immediately attempts exact mod-2
finite-reduction certification, including a capped small-prime saturation.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
from itertools import combinations, product
import json
from math import comb, gcd, isqrt, lcm, log
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import time
from typing import Any, Iterable, Sequence

from sympy import primerange

from alternate_quartic_covers import point_on_short_curve, short_add
from crt_lattice import crt_pair, short_rational_representatives
from kihara_discriminant_geometry import (
    derive_discriminant_geometry,
    hensel_lift_simple_t_root,
    homogeneous_frontier_value,
    polynomial_derivative_value,
    roots_mod_prime_t,
)
from kihara_rank14 import (
    binary_invariants,
    kihara_specialization,
    known_quartic_points,
)
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from nagao_1994 import quartic_covariants_at


Q = Fraction
TARGET_LOG_CONDUCTOR = Decimal("182.72")
CRT_PRIMES = (11, 17, 19, 41)
SMALL_RADICAL_PRIME_BOUND = 97
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_kihara_rank14_specializations.py "
    "--output artifacts/generated-results/elliptic_kihara_rank14_specializations.json"
)


@dataclass(frozen=True)
class PopulationCandidate:
    parameter_t: Fraction
    origins: tuple[str, ...]
    local_conditions: tuple[tuple[int, int, int], ...] = ()


@dataclass(frozen=True)
class NormalizedSpecialization:
    parameter_t: Fraction
    quartic_coefficients: tuple[Fraction, ...]
    ordinate_scale: Fraction
    short_coefficients: tuple[Fraction, ...]
    short_discriminant: Fraction


def q_string(value: Fraction) -> str:
    value = Q(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def point_record(point: tuple[Fraction, Fraction]) -> list[str]:
    return [q_string(point[0]), q_string(point[1])]


def _largest_square_divisor(
    value: int, *, allowed_squarefree_kernels: Sequence[int]
) -> int:
    """Return the square part after an exact small-kernel recognition.

    In Kihara's printed normalization the common denominator is a square and
    the common numerator has squarefree kernel in ``{1,2,5,10}``.  Recognizing
    this by integer square roots avoids asking a general-purpose factorizer to
    split hundreds-digit perfect squares.  Failure is explicit rather than an
    unverified normalization.
    """

    if value == 0:
        raise ValueError("the square divisor of zero is undefined")
    value = abs(value)
    matches = []
    for kernel in allowed_squarefree_kernels:
        if kernel <= 0 or value % kernel:
            continue
        root = isqrt(value // kernel)
        if root * root == value // kernel:
            matches.append(root * root)
    if len(matches) != 1:
        raise ValueError(
            "the Kihara square-content kernel left the certified small set"
        )
    return matches[0]


def normalized_specialization(parameter_t: Fraction) -> NormalizedSpecialization:
    """Remove the maximal rational square common to the quartic coefficients."""

    parameter_t = abs(Q(parameter_t))
    if not parameter_t:
        raise ValueError("Kihara's parameter t must be nonzero")
    raw = kihara_specialization(parameter_t).quartic_coefficients
    common_denominator = 1
    for coefficient in raw:
        common_denominator = lcm(common_denominator, coefficient.denominator)
    integer_coefficients = tuple(
        coefficient.numerator * (common_denominator // coefficient.denominator)
        for coefficient in raw
    )
    common_numerator = 0
    for coefficient in integer_coefficients:
        common_numerator = gcd(common_numerator, abs(coefficient))
    square_scale_squared = Q(
        _largest_square_divisor(
            common_numerator, allowed_squarefree_kernels=(1, 2, 5, 10)
        ),
        _largest_square_divisor(
            common_denominator, allowed_squarefree_kernels=(1,)
        ),
    )
    numerator_root = isqrt(square_scale_squared.numerator)
    denominator_root = isqrt(square_scale_squared.denominator)
    if (
        numerator_root**2 != square_scale_squared.numerator
        or denominator_root**2 != square_scale_squared.denominator
    ):
        raise AssertionError("the extracted square content was not a rational square")
    ordinate_scale = Q(numerator_root, denominator_root)
    quartic = tuple(coefficient / square_scale_squared for coefficient in raw)
    invariant_i, invariant_j = binary_invariants(quartic)
    short = (Q(0), Q(0), Q(0), -27 * invariant_i, -27 * invariant_j)
    discriminant = -16 * (4 * short[3] ** 3 + 27 * short[4] ** 2)
    if not discriminant:
        raise ValueError("the specialization is singular")
    return NormalizedSpecialization(
        parameter_t, quartic, ordinate_scale, short, discriminant
    )


def quartic_value(coefficients: Sequence[Fraction], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * Q(value) + Q(coefficient)
    return answer


def shifted_polynomial_coefficients(
    coefficients: Sequence[Fraction], center: Fraction
) -> tuple[Fraction, ...]:
    center = Q(center)
    return tuple(
        sum(
            Q(coefficients[source_degree])
            * comb(source_degree, target_degree)
            * center ** (source_degree - target_degree)
            for source_degree in range(target_degree, len(coefficients))
        )
        for target_degree in range(len(coefficients))
    )


def _merge_candidate(
    candidates: dict[Fraction, PopulationCandidate], candidate: PopulationCandidate
) -> None:
    parameter = abs(candidate.parameter_t)
    if not parameter:
        return
    previous = candidates.get(parameter)
    if previous is None:
        candidates[parameter] = PopulationCandidate(
            parameter,
            tuple(sorted(set(candidate.origins))),
            candidate.local_conditions,
        )
        return
    conditions = set(previous.local_conditions)
    conditions.update(candidate.local_conditions)
    candidates[parameter] = PopulationCandidate(
        parameter,
        tuple(sorted(set(previous.origins) | set(candidate.origins))),
        tuple(sorted(conditions)),
    )


def build_population(
    frontier_coefficients_z: Sequence[int],
    *,
    rational_height: int,
    integer_bound: int,
    crt_exponent: int,
    crt_coefficient_radius: int,
    crt_representatives_per_class: int,
) -> tuple[tuple[PopulationCandidate, ...], dict[str, Any]]:
    if rational_height < 1 or integer_bound < rational_height:
        raise ValueError("population bounds are inconsistent")
    candidates: dict[Fraction, PopulationCandidate] = {}
    for numerator in range(1, rational_height + 1):
        for denominator in range(1, rational_height + 1):
            if gcd(numerator, denominator) == 1:
                _merge_candidate(
                    candidates,
                    PopulationCandidate(Q(numerator, denominator), ("low-height",)),
                )
    for parameter in range(rational_height + 1, integer_bound + 1):
        _merge_candidate(
            candidates,
            PopulationCandidate(Q(parameter), ("extended-integer",)),
        )

    lifted_roots: dict[int, tuple[int, ...]] = {}
    for prime in CRT_PRIMES:
        simple_roots = tuple(
            root
            for root in roots_mod_prime_t(frontier_coefficients_z, prime)
            if (
                2
                * root
                * polynomial_derivative_value(frontier_coefficients_z, root * root)
            )
            % prime
        )
        lifted_roots[prime] = tuple(
            hensel_lift_simple_t_root(
                frontier_coefficients_z, root, prime, crt_exponent
            )
            for root in simple_roots
        )

    crt_class_count = 0
    crt_generated_count = 0
    for left_prime, right_prime in combinations(CRT_PRIMES, 2):
        left_modulus = left_prime**crt_exponent
        right_modulus = right_prime**crt_exponent
        origin = f"crt-{left_prime}-{right_prime}-power-{crt_exponent}"
        for left_root, right_root in product(
            lifted_roots[left_prime], lifted_roots[right_prime]
        ):
            crt_class_count += 1
            residue, modulus = crt_pair(
                left_root, left_modulus, right_root, right_modulus
            )
            representatives = short_rational_representatives(
                residue,
                modulus,
                coefficient_radius=crt_coefficient_radius,
                limit=crt_representatives_per_class,
            )
            for representative in representatives:
                parameter = abs(Q(representative.numerator, representative.denominator))
                if not parameter:
                    continue
                actual_conditions = tuple(
                    (
                        prime,
                        crt_exponent,
                        parameter.numerator
                        * pow(parameter.denominator, -1, prime**crt_exponent)
                        % (prime**crt_exponent),
                    )
                    for prime in (left_prime, right_prime)
                )
                for prime, exponent, _ in actual_conditions:
                    if (
                        homogeneous_frontier_value(
                            frontier_coefficients_z,
                            parameter.numerator,
                            parameter.denominator,
                        )
                        % (prime**exponent)
                    ):
                        raise AssertionError("a CRT representative lost its forced divisibility")
                _merge_candidate(
                    candidates,
                    PopulationCandidate(parameter, (origin,), actual_conditions),
                )
                crt_generated_count += 1

    ordered = tuple(
        sorted(
            candidates.values(),
            key=lambda candidate: (
                max(candidate.parameter_t.numerator, candidate.parameter_t.denominator),
                candidate.parameter_t,
            ),
        )
    )
    metadata = {
        "positive_t_only_reason": "the printed family and discriminant are invariant under t -> -t",
        "low_height_bound": rational_height,
        "extended_integer_bound": integer_bound,
        "crt_primes": list(CRT_PRIMES),
        "crt_hensel_exponent": crt_exponent,
        "crt_coefficient_radius": crt_coefficient_radius,
        "crt_representatives_per_class": crt_representatives_per_class,
        "lifted_roots": {
            str(prime): list(roots) for prime, roots in lifted_roots.items()
        },
        "crt_class_count": crt_class_count,
        "crt_representatives_before_deduplication": crt_generated_count,
        "distinct_candidate_count": len(ordered),
    }
    return ordered, metadata


def _integer_log(value: int) -> float:
    return log(abs(value)) if abs(value) > 1 else 0.0


def geometry_record(
    candidate: PopulationCandidate,
    frontier_coefficients_z: Sequence[int],
    *,
    radical_prime_bound: int,
) -> dict[str, Any]:
    specialization = normalized_specialization(candidate.parameter_t)
    combined = abs(
        specialization.short_discriminant.numerator
        * specialization.short_discriminant.denominator
    )
    cofactor = combined
    valuations: list[list[int]] = []
    radical_log = 0.0
    for prime in primerange(2, radical_prime_bound + 1):
        exponent = 0
        while cofactor % prime == 0:
            exponent += 1
            cofactor //= prime
        if exponent:
            valuations.append([int(prime), exponent])
            radical_log += log(int(prime))
    radical_log += _integer_log(cofactor)
    discriminant_log = _integer_log(combined)
    frontier_value = homogeneous_frontier_value(
        frontier_coefficients_z,
        candidate.parameter_t.numerator,
        candidate.parameter_t.denominator,
    )
    forced_valuations = []
    for prime, requested_exponent, residue in candidate.local_conditions:
        value = abs(frontier_value)
        actual_exponent = 0
        while value and value % prime == 0:
            actual_exponent += 1
            value //= prime
        forced_valuations.append(
            {
                "prime": prime,
                "requested_exponent": requested_exponent,
                "actual_exponent": actual_exponent,
                "parameter_residue": residue,
            }
        )
    return {
        "parameter_t": q_string(candidate.parameter_t),
        "parameter_height": max(
            candidate.parameter_t.numerator, candidate.parameter_t.denominator
        ),
        "origins": list(candidate.origins),
        "normalized_short_discriminant_log_abs_num_times_den": f"{discriminant_log:.12f}",
        "small_prime_radical_upper_proxy_log": f"{radical_log:.12f}",
        "powerful_gain_log": f"{discriminant_log - radical_log:.12f}",
        "small_prime_valuations": valuations,
        "unfactored_cofactor_digits": len(str(cofactor)),
        "frontier_homogeneous_value_bits": abs(frontier_value).bit_length(),
        "forced_frontier_valuations": forced_valuations,
    }


def select_conductor_candidates(
    records: Sequence[dict[str, Any]], *, geometry_keep: int, crt_origin_keep: int
) -> tuple[str, ...]:
    ordered = sorted(
        records,
        key=lambda record: (
            Decimal(record["small_prime_radical_upper_proxy_log"]),
            Decimal(record["normalized_short_discriminant_log_abs_num_times_den"]),
            int(record["parameter_height"]),
            Q(record["parameter_t"]),
        ),
    )
    selected = {record["parameter_t"] for record in ordered[:geometry_keep]}
    selected.update(("1", "2", "3"))
    crt_records = [
        record for record in ordered if any(origin.startswith("crt-") for origin in record["origins"])
    ]
    selected.update(record["parameter_t"] for record in crt_records[:crt_origin_keep])
    return tuple(
        record["parameter_t"]
        for record in ordered
        if record["parameter_t"] in selected
    )


def _terminate_process_group(process: subprocess.Popen[str]) -> None:
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
        process.wait(timeout=2)
    except ProcessLookupError:
        return
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait(timeout=2)


def run_gp_capped(
    program: str, *, timeout: float, stack_bytes: int
) -> tuple[str | None, dict[str, Any]]:
    if timeout <= 0 or timeout > 60:
        raise ValueError("every GP timeout must lie in (0,60]")
    executable = shutil.which("gp")
    if executable is None:
        return None, {"status": "unavailable", "wall_seconds": 0.0}
    process = subprocess.Popen(
        [executable, "-q", "-s", str(stack_bytes)],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    started = time.monotonic()
    try:
        stdout, stderr = process.communicate(program, timeout=timeout)
    except subprocess.TimeoutExpired:
        _terminate_process_group(process)
        return None, {
            "status": "timeout",
            "timeout_seconds": timeout,
            "wall_seconds": time.monotonic() - started,
        }
    except BaseException:
        _terminate_process_group(process)
        raise
    elapsed = time.monotonic() - started
    if process.returncode != 0 or "***" in stderr:
        return None, {
            "status": "pari_error",
            "wall_seconds": elapsed,
            "error": " ".join(stderr.split())[:1000],
        }
    return stdout, {"status": "completed", "wall_seconds": elapsed}


def gp_rational(value: Fraction) -> str:
    return f"({q_string(Q(value))})"


def gp_curve(coefficients: Sequence[Fraction]) -> str:
    return ",".join(gp_rational(Q(value)) for value in coefficients)


def gp_point(point: tuple[Fraction, Fraction]) -> str:
    return f"[{gp_rational(point[0])},{gp_rational(point[1])}]"


def gp_quartic(coefficients: Sequence[Fraction]) -> str:
    return "+".join(
        f"{gp_rational(Q(coefficient))}*x^{degree}"
        for degree, coefficient in enumerate(coefficients)
    )


def _block(lines: Sequence[str], name: str) -> list[str]:
    start = lines.index(f"{name}_BEGIN") + 1
    end = lines.index(f"{name}_END")
    return list(lines[start:end])


def conductor_probe(
    specialization: NormalizedSpecialization,
    *,
    timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    program = "\n".join(
        (
            "default(realprecision,80);",
            f"E=ellinit([{gp_curve(specialization.short_coefficients)}]);",
            "Em=ellminimalmodel(E);",
            "G=ellglobalred(Em);",
            'print("MODEL_BEGIN");',
            "print(Em.a1);print(Em.a2);print(Em.a3);print(Em.a4);print(Em.a6);",
            'print("MODEL_END");',
            'print("CONDUCTOR_BEGIN");print(G[1]);print("CONDUCTOR_END");',
            'print("LOG_BEGIN");print(log(G[1]));print("LOG_END");',
            'print("DISC_BEGIN");print(Em.disc);print("DISC_END");',
            'print("ROOT_BEGIN");print(ellrootno(Em));print("ROOT_END");',
            "quit",
        )
    ) + "\n"
    output, process_record = run_gp_capped(
        program, timeout=timeout, stack_bytes=stack_bytes
    )
    answer = dict(process_record)
    answer.update(
        {
            "timeout_seconds": timeout,
            "pari_stack_bytes": stack_bytes,
            "normalized_quartic_coefficients": [
                q_string(value) for value in specialization.quartic_coefficients
            ],
            "normalized_short_coefficients": [
                q_string(value) for value in specialization.short_coefficients
            ],
        }
    )
    if output is None:
        return answer
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    model = tuple(int(value) for value in _block(lines, "MODEL"))
    conductor = int(_block(lines, "CONDUCTOR")[0])
    log_conductor = _block(lines, "LOG")[0]
    answer.update(
        {
            "minimal_model": list(model),
            "conductor": str(conductor),
            "log_conductor": log_conductor,
            "minimal_discriminant": _block(lines, "DISC")[0],
            "root_number": int(_block(lines, "ROOT")[0]),
            "below_target": Decimal(log_conductor) < TARGET_LOG_CONDUCTOR,
        }
    )
    return answer


POINT_PATTERN = re.compile(r"\[(-?\d+(?:/\d+)?),\s*(-?\d+(?:/\d+)?)\]")


def parse_points(text: str) -> tuple[tuple[Fraction, Fraction], ...]:
    return tuple((Q(x_value), Q(y_value)) for x_value, y_value in POINT_PATTERN.findall(text))


def bounded_centered_quartic_search(
    specialization: NormalizedSpecialization,
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    known = tuple(
        (x_value, y_value / specialization.ordinate_scale)
        for x_value, y_value in known_quartic_points(specialization.parameter_t)
    )
    centers = (Q(0),) + tuple(point[0] for point in known)
    commands = ["gettime();"]
    for index, center in enumerate(centers):
        shifted = shifted_polynomial_coefficients(
            specialization.quartic_coefficients, center
        )
        commands.extend(
            (
                f"Q={gp_quartic(shifted)};",
                f"R=hyperellratpoints(Q,{height_bound});",
                f'print("ROW|{index}|",R);',
            )
        )
    commands.extend(('print("PARI_MILLISECONDS ",gettime());', "quit"))
    output, process_record = run_gp_capped(
        "\n".join(commands) + "\n",
        timeout=timeout,
        stack_bytes=stack_bytes,
    )
    record = dict(process_record)
    record.update(
        {
            "height_bound": height_bound,
            "timeout_seconds": timeout,
            "pari_stack_bytes": stack_bytes,
            "chart_count": len(centers),
            "chart_centers": "0 plus all fifteen exact printed abscissae",
        }
    )
    if output is None:
        return (), record
    found: list[tuple[Fraction, Fraction]] = []
    row_count = 0
    for line in output.splitlines():
        if not line.startswith("ROW|"):
            continue
        _, index_text, vector_text = line.split("|", 2)
        center = centers[int(index_text)]
        row_count += 1
        found.extend((offset + center, ordinate) for offset, ordinate in parse_points(vector_text))
    if row_count != len(centers):
        raise AssertionError("PARI omitted one or more centered-search charts")
    distinct = tuple(dict.fromkeys(found))
    if any(
        y_value**2
        != quartic_value(specialization.quartic_coefficients, x_value)
        for x_value, y_value in distinct
    ):
        raise AssertionError("PARI returned a point off the exact quartic")
    milliseconds = re.search(r"^PARI_MILLISECONDS (\d+)$", output, re.MULTILINE)
    if milliseconds is None:
        raise AssertionError("PARI omitted its point-search timing")
    record.update(
        {
            "pari_milliseconds": int(milliseconds.group(1)),
            "signed_points_before_deduplication": len(found),
            "distinct_signed_points": len(distinct),
            "distinct_abscissae": len({point[0] for point in distinct}),
        }
    )
    return distinct, record


def _negate(point: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return point[0], -point[1]


def covariant_difference_pool(
    specialization: NormalizedSpecialization,
    searched_points: Iterable[tuple[Fraction, Fraction]],
) -> tuple[
    tuple[tuple[Fraction, Fraction], ...],
    tuple[dict[str, Any], ...],
]:
    known = tuple(
        (x_value, y_value / specialization.ordinate_scale)
        for x_value, y_value in known_quartic_points(specialization.parameter_t)
    )

    def image(point: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
        x_value, y_value = point
        if not y_value or y_value**2 != quartic_value(
            specialization.quartic_coefficients, x_value
        ):
            raise ValueError("the covariant map needs a checked nonzero ordinate")
        g_value, h_value = quartic_covariants_at(
            specialization.quartic_coefficients, x_value
        )
        result = (36 * g_value / y_value**2, 108 * h_value / y_value**3)
        if not point_on_short_curve(specialization.short_coefficients, result):
            raise AssertionError("a covariant image missed the exact short curve")
        return result

    images = tuple(image(point) for point in known)
    origin_negative = _negate(images[14])
    baseline = tuple(
        short_add(specialization.short_coefficients, value, origin_negative)
        for value in images[:14]
    )
    if any(value is None for value in baseline):
        raise AssertionError("a known section collided with the selected origin")
    known_x = {point[0] for point in known}
    signless: dict[Fraction, tuple[Fraction, Fraction]] = {}
    for point in searched_points:
        signless.setdefault(point[0], point)
    extras: list[tuple[Fraction, Fraction]] = []
    records: list[dict[str, Any]] = []
    for x_value, y_value in signless.values():
        if x_value in known_x or not y_value:
            continue
        jacobian_image = image((x_value, y_value))
        difference = short_add(
            specialization.short_coefficients, jacobian_image, origin_negative
        )
        if difference is None:
            continue
        if difference in extras or _negate(difference) in extras:
            continue
        extras.append(difference)
        records.append(
            {
                "quartic_point": point_record((x_value, y_value)),
                "jacobian_covariant_image": point_record(jacobian_image),
                "jacobian_difference_from_P15_image": point_record(difference),
                "exact_quartic_membership": True,
                "exact_jacobian_membership": True,
            }
        )
    return tuple(baseline) + tuple(extras), tuple(records)  # type: ignore[arg-type]


def _parse_integer_vector(text: str) -> list[int]:
    match = re.search(r"\[(.*?)\]", text)
    if match is None or not match.group(1).strip():
        return []
    return [int(value.strip()) for value in match.group(1).split(",")]


def height_replay(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    precisions: Sequence[int],
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[dict[str, Any], ...], dict[str, Any]]:
    if any(not point_on_short_curve(coefficients, point) for point in points):
        raise AssertionError("a height point missed the exact curve")
    commands = [
        f"E=ellinit([{gp_curve(coefficients)}]);",
        f"P=[{','.join(gp_point(point) for point in points)}];",
        "gettime();",
    ]
    for precision in precisions:
        commands.extend(
            (
                f"default(realprecision,{precision});",
                "H=ellheightmatrix(E,P);IX=matindexrank(H);K=vecextract(P,IX[2]);HK=ellheightmatrix(E,K);",
                f'print("HEIGHT_{precision}_BEGIN");',
                "print(matrank(H));print(Vec(IX[2]));print(matdet(HK));",
                f'print("HEIGHT_{precision}_END");',
            )
        )
    commands.extend(('print("PARI_MILLISECONDS ",gettime());', "quit"))
    output, process_record = run_gp_capped(
        "\n".join(commands) + "\n", timeout=timeout, stack_bytes=stack_bytes
    )
    if output is None:
        return (), process_record
    lines = [line.strip() for line in output.splitlines() if line.strip()]
    runs = []
    for precision in precisions:
        values = _block(lines, f"HEIGHT_{precision}")
        runs.append(
            {
                "decimal_precision": precision,
                "numerical_rank": int(values[0]),
                "subset_indices_one_based": _parse_integer_vector(values[1]),
                "subset_height_determinant": values[2],
            }
        )
    ranks = {run["numerical_rank"] for run in runs}
    subsets = {tuple(run["subset_indices_one_based"]) for run in runs}
    stability = dict(process_record)
    stability.update(
        {
            "stable_across_precisions": len(ranks) == 1 and len(subsets) == 1,
            "stable_numerical_rank": next(iter(ranks)) if len(ranks) == 1 and len(subsets) == 1 else None,
        }
    )
    milliseconds = re.search(r"^PARI_MILLISECONDS (\d+)$", output, re.MULTILINE)
    if milliseconds:
        stability["pari_milliseconds"] = int(milliseconds.group(1))
    return tuple(runs), stability


def _signature_records(signatures: Sequence[Any]) -> list[dict[str, Any]]:
    return [
        {
            "prime": signature.prime,
            "group_order": signature.group_order,
            "doubled_subgroup_order": signature.doubled_subgroup_order,
            "quotient_dimension": signature.quotient_dimension,
            "rows": [list(row) for row in signature.rows],
        }
        for signature in signatures
    ]


def _parse_saturation_points(output: str) -> tuple[tuple[Fraction, Fraction], ...]:
    marker = re.search(r"SATURATION_BEGIN\n(.*?)\nSATURATION_END", output, re.DOTALL)
    if marker is None:
        raise AssertionError("PARI omitted saturation points")
    return parse_points(marker.group(1))


def exact_rank_certificate_attempt(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    subset_indices_one_based: Sequence[int],
    prime_bound: int,
    saturation_timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    basis = tuple(points[index - 1] for index in subset_indices_one_based)
    if len(basis) < 18:
        return {"status": "not_attempted", "reason": "stable numerical rank below 18"}

    def reductions(candidate_basis: Sequence[tuple[Fraction, Fraction]]) -> dict[str, Any]:
        signatures = find_mod2_reduction_certificate(
            coefficients, candidate_basis, prime_bound=prime_bound
        )
        rank = combined_mod2_rank(signatures, len(candidate_basis))
        result: dict[str, Any] = {
            "point_count": len(candidate_basis),
            "combined_mod2_rank": rank,
            "signatures": _signature_records(signatures),
        }
        if rank == len(candidate_basis):
            result["no_rational_2_torsion_prime"] = find_two_torsion_certificate_prime(
                coefficients, prime_bound=200
            )
            result["certified_rank_lower_bound"] = len(candidate_basis)
        return result

    direct = reductions(basis)
    if direct["combined_mod2_rank"] == len(basis):
        return {"status": "certified_directly", "direct": direct}
    program = "\n".join(
        (
            f"E=ellinit([{gp_curve(coefficients)}]);",
            f"P=[{','.join(gp_point(point) for point in basis)}];",
            "S=ellsaturation(E,P,20);",
            'print("SATURATION_BEGIN");print(S);print("SATURATION_END");',
            "quit",
        )
    ) + "\n"
    output, process_record = run_gp_capped(
        program, timeout=saturation_timeout, stack_bytes=stack_bytes
    )
    if output is None:
        return {
            "status": "direct_rank_deficient_saturation_failed",
            "direct": direct,
            "saturation_process": process_record,
        }
    saturated = _parse_saturation_points(output)
    if len(saturated) != len(basis) or any(
        not point_on_short_curve(coefficients, point) for point in saturated
    ):
        raise AssertionError("PARI saturation returned an invalid exact basis")
    replay = reductions(saturated)
    return {
        "status": (
            "certified_after_saturation"
            if replay["combined_mod2_rank"] == len(saturated)
            else "bounded_certificate_search_rank_deficient"
        ),
        "direct": direct,
        "saturation_process": process_record,
        "saturated_basis": [point_record(point) for point in saturated],
        "saturated_reduction": replay,
    }


def _sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def build_search(args: argparse.Namespace) -> dict[str, Any]:
    geometry = derive_discriminant_geometry()
    frontier = geometry.frontier_coefficients_z
    population, population_metadata = build_population(
        frontier,
        rational_height=args.rational_height,
        integer_bound=args.integer_bound,
        crt_exponent=args.crt_exponent,
        crt_coefficient_radius=args.crt_coefficient_radius,
        crt_representatives_per_class=args.crt_representatives_per_class,
    )
    geometry_records = tuple(
        geometry_record(
            candidate, frontier, radical_prime_bound=args.radical_prime_bound
        )
        for candidate in population
    )
    selected_parameters = select_conductor_candidates(
        geometry_records,
        geometry_keep=args.geometry_keep,
        crt_origin_keep=args.crt_origin_keep,
    )
    geometry_by_t = {record["parameter_t"]: record for record in geometry_records}

    conductor_records = []
    specializations: dict[str, NormalizedSpecialization] = {}
    for parameter_text in selected_parameters:
        specialization = normalized_specialization(Q(parameter_text))
        specializations[parameter_text] = specialization
        probe = conductor_probe(
            specialization,
            timeout=args.conductor_timeout,
            stack_bytes=args.stack_bytes,
        )
        conductor_records.append(
            {
                "parameter_t": parameter_text,
                "geometry_selection_objective": geometry_by_t[parameter_text][
                    "small_prime_radical_upper_proxy_log"
                ],
                "pari": probe,
            }
        )

    completed = [
        record for record in conductor_records if record["pari"]["status"] == "completed"
    ]
    subthreshold = [record for record in completed if record["pari"]["below_target"]]
    completed_order = sorted(
        completed,
        key=lambda record: (Decimal(record["pari"]["log_conductor"]), Q(record["parameter_t"])),
    )
    point_parameters = {record["parameter_t"] for record in subthreshold}
    point_parameters.update(
        record["parameter_t"]
        for record in completed_order[: args.point_fallback_keep]
    )
    point_parameters = set(
        sorted(
            point_parameters,
            key=lambda text: (
                Decimal(next(record for record in completed if record["parameter_t"] == text)["pari"]["log_conductor"]),
                Q(text),
            ),
        )[: args.point_candidate_cap]
    )

    initial_searches: dict[str, tuple[tuple[Fraction, Fraction], ...]] = {}
    point_records: dict[str, dict[str, Any]] = {}
    for parameter_text in sorted(point_parameters, key=Q):
        specialization = specializations[parameter_text]
        found, search_record = bounded_centered_quartic_search(
            specialization,
            height_bound=args.point_height,
            timeout=args.point_timeout,
            stack_bytes=args.stack_bytes,
        )
        initial_searches[parameter_text] = found
        known_x = {point[0] for point in known_quartic_points(Q(parameter_text))}
        point_records[parameter_text] = {
            "parameter_t": parameter_text,
            "initial_search": search_record,
            "initial_unexpected_abscissa_count": len(
                {point[0] for point in found if point[0] not in known_x}
            ),
        }

    escalation_order = sorted(
        point_records,
        key=lambda text: (
            -point_records[text]["initial_unexpected_abscissa_count"],
            Decimal(next(record for record in completed if record["parameter_t"] == text)["pari"]["log_conductor"]),
            Q(text),
        ),
    )[: args.escalation_keep]
    for parameter_text in escalation_order:
        specialization = specializations[parameter_text]
        found, search_record = bounded_centered_quartic_search(
            specialization,
            height_bound=args.escalated_height,
            timeout=args.escalated_timeout,
            stack_bytes=args.stack_bytes,
        )
        initial_searches[parameter_text] = tuple(
            dict.fromkeys((*initial_searches[parameter_text], *found))
        )
        point_records[parameter_text]["escalated_search"] = search_record

    exact_certificates: dict[str, dict[str, Any]] = {}
    for parameter_text in sorted(point_records, key=Q):
        specialization = specializations[parameter_text]
        pool, extras = covariant_difference_pool(
            specialization, initial_searches[parameter_text]
        )
        runs, stability = height_replay(
            specialization.short_coefficients,
            pool,
            precisions=args.precisions,
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
        )
        record = point_records[parameter_text]
        record.update(
            {
                "exact_known_difference_count": 14,
                "unexpected_exact_points": list(extras),
                "unexpected_independent_direction_upper_count": len(extras),
                "height_pool_count": len(pool),
                "height_runs": list(runs),
                "height_stability": stability,
            }
        )
        stable_rank = stability.get("stable_numerical_rank")
        if stable_rank is not None:
            record["stable_numerical_rank_gain_over_generic_14"] = int(stable_rank) - 14
        if stable_rank is not None and int(stable_rank) >= 18:
            subset = runs[-1]["subset_indices_one_based"]
            exact_certificates[parameter_text] = exact_rank_certificate_attempt(
                specialization.short_coefficients,
                pool,
                subset_indices_one_based=subset,
                prime_bound=args.certificate_prime_bound,
                saturation_timeout=args.saturation_timeout,
                stack_bytes=args.stack_bytes,
            )

    stable_rows = [
        record
        for record in point_records.values()
        if record["height_stability"].get("stable_numerical_rank") is not None
    ]
    strongest = max(
        stable_rows,
        key=lambda record: (
            int(record["height_stability"]["stable_numerical_rank"]),
            -Decimal(next(row for row in completed if row["parameter_t"] == record["parameter_t"])["pari"]["log_conductor"]),
        ),
        default=None,
    )
    basis_artifact = Path(args.basis_artifact)
    basis_digest = hashlib.sha256(basis_artifact.read_bytes()).hexdigest()
    frontier_digest = _sha256_lines(str(value) for value in frontier)
    reached_target = any(
        certificate.get("status") in ("certified_directly", "certified_after_saturation")
        and (
            certificate.get("direct", {}).get("certified_rank_lower_bound", 0) >= 21
            or certificate.get("saturated_reduction", {}).get("certified_rank_lower_bound", 0) >= 21
        )
        and next(record for record in completed if record["parameter_t"] == parameter)["pari"]["below_target"]
        for parameter, certificate in exact_certificates.items()
    )
    return {
        "schema": "elliptic-kihara-rank14-specialization-search-v1",
        "reproduce_command": REPRODUCING_COMMAND,
        "status": (
            "target curve exactly certified"
            if reached_target
            else "negative bounded specialization search; no target curve certified"
        ),
        "target": {
            "algebraic_rank_at_least": 21,
            "strict_log_conductor_less_than": str(TARGET_LOG_CONDUCTOR),
        },
        "generic_rank14_basis_input": {
            "path": str(basis_artifact),
            "sha256": basis_digest,
            "exact_specialized_rank_lower_bound": 14,
        },
        "discriminant_geometry": {
            "scaled_discriminant_degree": geometry.discriminant_degree,
            "factor_signature_degree_exponent": [list(pair) for pair in geometry.factor_signature],
            "frontier_description": "irreducible even degree 398, represented by primitive degree-199 f(z) with z=t^2",
            "frontier_degree_z": geometry.frontier_degree_z,
            "frontier_coefficients_sha256": frontier_digest,
            "frontier_first_three_coefficients": list(frontier[:3]),
            "frontier_last_three_coefficients": list(frontier[-3:]),
            "conductor_geometry_interpretation": (
                "all repeated low-degree factors are conductor-friendly; the squarefree degree-398 factor is the asymptotic conductor bottleneck"
            ),
        },
        "population": population_metadata,
        "geometry_prefilter": {
            "uses_point_or_rank_data": False,
            "small_prime_bound": args.radical_prime_bound,
            "objective": (
                "log radical upper proxy: exact distinct small-prime support plus the unfactored cofactor treated as squarefree"
            ),
            "records": list(geometry_records),
            "selected_parameters": list(selected_parameters),
        },
        "conductor_stage": {
            "candidate_count": len(conductor_records),
            "completed_count": len(completed),
            "subthreshold_count": len(subthreshold),
            "records": conductor_records,
        },
        "point_stage_selection": {
            "rule": (
                "all completed sub-threshold conductors, then the declared lowest-conductor completed fallbacks, capped before point search"
            ),
            "fallback_keep": args.point_fallback_keep,
            "candidate_cap": args.point_candidate_cap,
            "parameters": sorted(point_parameters, key=Q),
            "escalation_rule": "unexpected exact abscissa yield, then lower conductor",
            "escalated_parameters": escalation_order,
        },
        "point_stage": {
            "records": [point_records[text] for text in sorted(point_records, key=Q)],
            "strongest_stable_numerical_result": strongest,
            "status_warning": (
                "exact point membership and bounded enumeration are certified; height-matrix ranks are numerical evidence only"
            ),
        },
        "exact_rank_certificate_attempts": exact_certificates,
        "target_reached": reached_target,
        "bounded_scope": {
            "rational_height": args.rational_height,
            "integer_bound": args.integer_bound,
            "point_height": args.point_height,
            "escalated_height": args.escalated_height,
            "per_conductor_timeout_seconds": args.conductor_timeout,
            "per_point_search_timeout_seconds": args.point_timeout,
            "per_escalated_search_timeout_seconds": args.escalated_timeout,
            "height_timeout_seconds": args.height_timeout,
            "all_external_processes_foreground": True,
            "all_external_processes_have_separate_process_groups": True,
            "timeout_cleanup": "SIGTERM followed by SIGKILL after two seconds if necessary",
        },
        "software": {
            "python": platform.python_version(),
            "pari_gp": shutil.which("gp"),
        },
    }


def parse_precisions(value: str) -> tuple[int, ...]:
    try:
        result = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("precisions must be integers") from error
    if len(result) < 2 or result != tuple(sorted(set(result))) or result[0] < 32:
        raise argparse.ArgumentTypeError("provide at least two increasing precisions >=32")
    return result


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rational-height", type=int, default=20)
    parser.add_argument("--integer-bound", type=int, default=128)
    parser.add_argument("--crt-exponent", type=int, default=2)
    parser.add_argument("--crt-coefficient-radius", type=int, default=12)
    parser.add_argument("--crt-representatives-per-class", type=int, default=3)
    parser.add_argument("--radical-prime-bound", type=int, default=SMALL_RADICAL_PRIME_BOUND)
    parser.add_argument("--geometry-keep", type=int, default=12)
    parser.add_argument("--crt-origin-keep", type=int, default=2)
    parser.add_argument("--conductor-timeout", type=float, default=4.0)
    parser.add_argument("--point-fallback-keep", type=int, default=2)
    parser.add_argument("--point-candidate-cap", type=int, default=4)
    parser.add_argument("--point-height", type=int, default=5_000)
    parser.add_argument("--point-timeout", type=float, default=12.0)
    parser.add_argument("--escalation-keep", type=int, default=1)
    parser.add_argument("--escalated-height", type=int, default=100_000)
    parser.add_argument("--escalated-timeout", type=float, default=20.0)
    parser.add_argument("--precisions", type=parse_precisions, default=(72, 120))
    parser.add_argument("--height-timeout", type=float, default=15.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=500)
    parser.add_argument("--saturation-timeout", type=float, default=20.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--basis-artifact",
        type=Path,
        default=root / "artifacts" / "generated-results" / "elliptic_kihara_rank14_basis.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root / "artifacts" / "generated-results" / "elliptic_kihara_rank14_specializations.json",
    )
    return parser


def validate_args(args: argparse.Namespace) -> None:
    positive_integers = (
        args.rational_height,
        args.integer_bound,
        args.crt_exponent,
        args.crt_coefficient_radius,
        args.crt_representatives_per_class,
        args.radical_prime_bound,
        args.geometry_keep,
        args.point_candidate_cap,
        args.point_height,
        args.escalated_height,
        args.certificate_prime_bound,
        args.stack_bytes,
    )
    if any(value <= 0 for value in positive_integers):
        raise SystemExit("all search bounds must be positive")
    if args.integer_bound < args.rational_height:
        raise SystemExit("integer bound must be at least rational height")
    if args.crt_exponent < 2:
        raise SystemExit("CRT Hensel exponent must be at least two")
    if args.point_height >= args.escalated_height:
        raise SystemExit("escalated point height must exceed initial height")
    if args.crt_origin_keep < 0 or args.point_fallback_keep < 0 or args.escalation_keep < 0:
        raise SystemExit("keep counts cannot be negative")
    for timeout in (
        args.conductor_timeout,
        args.point_timeout,
        args.escalated_timeout,
        args.height_timeout,
        args.saturation_timeout,
    ):
        if timeout <= 0 or timeout > 60:
            raise SystemExit("every external-process timeout must lie in (0,60]")


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    validate_args(args)
    result = build_search(args)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded, encoding="utf-8")
    strongest = result["point_stage"]["strongest_stable_numerical_result"]
    print(f"population={result['population']['distinct_candidate_count']}")
    print(f"conductor_completed={result['conductor_stage']['completed_count']}")
    print(f"subthreshold={result['conductor_stage']['subthreshold_count']}")
    if strongest:
        print(
            "strongest_stable_numerical_rank="
            f"{strongest['height_stability']['stable_numerical_rank']} "
            f"at t={strongest['parameter_t']}"
        )
    print(f"target_reached={result['target_reached']}")


if __name__ == "__main__":
    main()
