#!/usr/bin/env python3
"""Exact max-root-50 Mestre tuple enumeration and bounded fiber screen.

This is a standalone higher-dimensional continuation of the small root-tuple
survey.  A compiled 128-bit enumerator exhausts primitive integer root tuples
modulo translation, integral scaling, and reflection.  Python independently
replays Mestre's degree-five obstruction and then rejects every nonreflection
tuple whose quartic discriminant vanishes at 21 consecutive parameters.  The
discriminant has degree at most 20, so this is an exact generic-singularity
test rather than a sampling heuristic.

The specialization stage is deliberately leakage-free.  First, every
admissible integer parameter in the declared box receives an exact PARI
conductor computation.  Only after that phase has closed are the twelve
displayed quartic points mapped exactly to the Jacobian and ranked by
two-precision numerical height matrices.  A fixed number of the strongest
signals then receives one bounded ``hyperellratpoints`` search and an
effort-zero ``ellrank`` probe.  Numerical ranks are triage evidence only.  An
exact finite-reduction independence attempt is triggered immediately, and
only, if a stable numerical rank reaches 21.
"""

from __future__ import annotations

import argparse
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
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
from typing import Any, Iterable, Sequence

from ek_k3 import rational_to_string
from mestre_root_tuples import (
    SixRootMestreConstruction,
    mestre_quartic_condition,
    normalize_integer_root_tuple,
)
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)


Q = Fraction
TARGET_LOG_CONDUCTOR = Decimal("182.72")
NAGAO_NORMALIZED_ROOTS = (0, 1, 27, 28, 31, 34)
DEFAULT_OUTPUT = Path(
    "archive/elliptic-curves/artifacts/generated-results/elliptic_mestre_root_tuple_scale.json"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_mestre_root_tuple_scale.py"
)
MAX_SCRIPT_ROOT = 60
MAX_PARAMETER_BOUND = 12
MAX_POINT_KEEP = 30
MAX_HEIGHT_BOUND = 50_000
MAX_SUBPROCESS_TIMEOUT = 30.0


class CappedProcessTimeout(RuntimeError):
    """A foreground subprocess exceeded its declared wall cap."""


@dataclass(frozen=True)
class EnumerationResult:
    max_root: int
    normalized_count: int
    obstruction_count: int
    reflection_count: int
    nonreflection_count: int
    obstruction_roots: tuple[tuple[int, ...], ...]
    reflection_roots: tuple[tuple[int, ...], ...]
    nonreflection_roots: tuple[tuple[int, ...], ...]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def tuple_digest(roots: Iterable[Sequence[int]]) -> str:
    text = "\n".join(",".join(str(value) for value in item) for item in roots)
    return hashlib.sha256(text.encode()).hexdigest()


def point_digest(points: Iterable[tuple[Fraction, Fraction]]) -> str:
    text = "\n".join(
        f"{rational_to_string(x_value)},{rational_to_string(y_value)}"
        for x_value, y_value in points
    )
    return hashlib.sha256(text.encode()).hexdigest()


def run_capped_process(
    command: Sequence[str],
    *,
    timeout: float,
    input_text: str | None = None,
) -> tuple[str, str]:
    """Run one foreground process group and kill the whole group on timeout."""

    process = subprocess.Popen(
        tuple(command),
        stdin=subprocess.PIPE if input_text is not None else subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    try:
        stdout, stderr = process.communicate(input=input_text, timeout=timeout)
    except subprocess.TimeoutExpired as error:
        os.killpg(process.pid, signal.SIGKILL)
        process.communicate()
        raise CappedProcessTimeout(
            f"subprocess exceeded its {timeout:g}-second wall cap"
        ) from error
    if process.returncode != 0:
        raise RuntimeError(
            f"subprocess exited {process.returncode}: {' '.join(stderr.split())[:1000]}"
        )
    return stdout, stderr


def compiled_enumeration(
    max_root: int,
    *,
    compiler: str = "c++",
    compile_timeout: float = 30.0,
    enumeration_timeout: float = 30.0,
) -> EnumerationResult:
    """Compile in a temporary directory and parse the exhaustive tuple stream."""

    if max_root < 5 or max_root > MAX_SCRIPT_ROOT:
        raise ValueError(f"max_root must lie in [5,{MAX_SCRIPT_ROOT}]")
    executable = shutil.which(compiler)
    if executable is None:
        raise FileNotFoundError(f"C++ compiler {compiler!r} was not found")
    source = Path(__file__).with_name("enumerate_mestre_root_tuples_scale.cpp")
    with tempfile.TemporaryDirectory(prefix="mestre-root-tuples-") as directory:
        binary = Path(directory) / "enumerator"
        run_capped_process(
            (
                executable,
                "-std=c++17",
                "-O3",
                "-DNDEBUG",
                str(source),
                "-o",
                str(binary),
            ),
            timeout=compile_timeout,
        )
        stdout, _ = run_capped_process(
            (str(binary), str(max_root)), timeout=enumeration_timeout
        )

    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not lines or lines[0] != "MESTRE_ROOT_TUPLES_V1":
        raise AssertionError("the compiled enumerator omitted its format header")
    roots: list[tuple[int, ...]] = []
    reflection: list[tuple[int, ...]] = []
    nonreflection: list[tuple[int, ...]] = []
    summary: tuple[int, ...] | None = None
    for line in lines[1:]:
        fields = line.split()
        if fields[0] == "R":
            if len(fields) != 8:
                raise AssertionError("malformed root record from compiled enumerator")
            item = tuple(int(value) for value in fields[1:7])
            flag = int(fields[7])
            if flag not in (0, 1):
                raise AssertionError("malformed reflection flag")
            roots.append(item)
            (reflection if flag else nonreflection).append(item)
        elif fields[0] == "S":
            if len(fields) != 6 or summary is not None:
                raise AssertionError("malformed enumerator summary")
            summary = tuple(int(value) for value in fields[1:])
        else:
            raise AssertionError("unknown compiled-enumerator record")
    if summary is None:
        raise AssertionError("the compiled enumerator omitted its summary")
    declared_root, normalized, obstruction, reflected, nonreflected = summary
    if declared_root != max_root:
        raise AssertionError("the compiled enumerator changed the root bound")
    if obstruction != len(roots) or reflected != len(reflection):
        raise AssertionError("compiled-enumerator counts disagree with its stream")
    if nonreflected != len(nonreflection) or reflected + nonreflected != obstruction:
        raise AssertionError("compiled-enumerator reflection counts disagree")
    if roots != sorted(roots, key=lambda item: (item[-1], item)):
        raise AssertionError("the compiled-enumerator order changed")
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


def reflection_symmetric(roots: Sequence[int]) -> bool:
    total = roots[0] + roots[-1]
    return all(roots[index] + roots[-1 - index] == total for index in range(3))


def verify_enumerator_records(enumeration: EnumerationResult) -> None:
    """Independently replay every emitted normalization and obstruction."""

    seen: set[tuple[int, ...]] = set()
    for roots in enumeration.obstruction_roots:
        if roots in seen:
            raise AssertionError("the compiled enumerator emitted a duplicate")
        seen.add(roots)
        if len(roots) != 6 or roots != tuple(sorted(roots)) or len(set(roots)) != 6:
            raise AssertionError("the compiled enumerator emitted invalid roots")
        if roots[0] != 0 or roots[-1] > enumeration.max_root:
            raise AssertionError("the compiled enumerator escaped its bound")
        if gcd(*roots[1:]) != 1 or normalize_integer_root_tuple(roots) != roots:
            raise AssertionError("the compiled enumerator missed normalization")
        construction = SixRootMestreConstruction(tuple(Q(root) for root in roots))
        if mestre_quartic_condition(construction.polynomial) != 0:
            raise AssertionError("the compiled enumerator missed Mestre's obstruction")
        expected_reflection = reflection_symmetric(roots)
        if expected_reflection != (roots in enumeration.reflection_roots):
            raise AssertionError("the compiled reflection gate changed")


def generic_nonsingularity_witness(
    construction: SixRootMestreConstruction,
) -> int | None:
    """Return a nonzero exact discriminant value, or prove it identically zero.

    The binary-quartic discriminant has degree at most 20 in the parameter.
    Consequently vanishing at the 21 distinct integers 1,...,21 is equivalent
    to generic singularity.
    """

    for parameter in range(1, 22):
        if construction.quartic_discriminant(Q(parameter)) != 0:
            return parameter
    return None


def classify_nonreflection(
    enumeration: EnumerationResult,
) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...], dict[tuple[int, ...], int]]:
    nonsingular: list[tuple[int, ...]] = []
    singular: list[tuple[int, ...]] = []
    witnesses: dict[tuple[int, ...], int] = {}
    for roots in enumeration.nonreflection_roots:
        construction = SixRootMestreConstruction(tuple(Q(root) for root in roots))
        witness = generic_nonsingularity_witness(construction)
        if witness is None:
            singular.append(roots)
        else:
            nonsingular.append(roots)
            witnesses[roots] = witness
    return tuple(nonsingular), tuple(singular), witnesses


def max14_calibration(
    enumeration: EnumerationResult,
    bounded_enumeration: EnumerationResult | None = None,
) -> dict[str, Any]:
    roots = tuple(item for item in enumeration.obstruction_roots if item[-1] <= 14)
    normalized_count = 1023
    if bounded_enumeration is not None:
        if bounded_enumeration.max_root != 14:
            raise ValueError("the bounded calibration enumeration must have max_root=14")
        if bounded_enumeration.obstruction_roots != roots:
            raise AssertionError("the max-root-14 obstruction population changed")
        normalized_count = bounded_enumeration.normalized_count
    nonsingular = []
    for item in roots:
        construction = SixRootMestreConstruction(tuple(Q(root) for root in item))
        if generic_nonsingularity_witness(construction) is not None:
            nonsingular.append(item)
    nonreflection = [item for item in nonsingular if not reflection_symmetric(item)]
    return {
        "affine_normalized_root_tuples": normalized_count,
        "obstruction_zero": len(roots),
        "generically_nonsingular": len(nonsingular),
        "generically_nonsingular_nonreflection": len(nonreflection),
        "nonsingular_nonreflection_tuples": [list(item) for item in nonreflection],
    }


def gp_rational(value: Fraction) -> str:
    return f"({rational_to_string(Q(value))})"


def gp_vector(point: tuple[Fraction, Fraction]) -> str:
    return f"[{gp_rational(point[0])},{gp_rational(point[1])}]"


def gp_block(lines: Sequence[str], name: str) -> list[str]:
    start = lines.index(f"{name}_BEGIN") + 1
    end = lines.index(f"{name}_END")
    return list(lines[start:end])


def capped_gp(
    program: str, *, timeout: float, stack_bytes: int
) -> tuple[str, str]:
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    stdout, stderr = run_capped_process(
        (executable, "-q", "-s", str(stack_bytes)),
        timeout=timeout,
        input_text=program,
    )
    if "***" in stderr:
        raise RuntimeError(f"PARI/GP failed: {' '.join(stderr.split())[:1000]}")
    return stdout, stderr


def capped_minimal_curve_data(
    coefficients: Sequence[Fraction], *, timeout: float, stack_bytes: int
) -> dict[str, Any]:
    vector = ",".join(gp_rational(Q(value)) for value in coefficients)
    program = "\n".join(
        (
            "default(realprecision,60);",
            f"E=ellinit([{vector}]);",
            "Em=ellminimalmodel(E);",
            "G=ellglobalred(Em);",
            'print("MODEL_BEGIN");',
            "print(Em.a1);print(Em.a2);print(Em.a3);print(Em.a4);print(Em.a6);",
            'print("MODEL_END");',
            'print("CONDUCTOR_BEGIN");print(G[1]);print("CONDUCTOR_END");',
            'print("LOG_CONDUCTOR_BEGIN");print(log(G[1]));print("LOG_CONDUCTOR_END");',
            'print("DISCRIMINANT_BEGIN");print(Em.disc);print("DISCRIMINANT_END");',
            'print("ROOT_NUMBER_BEGIN");print(ellrootno(Em));print("ROOT_NUMBER_END");',
            "quit",
        )
    ) + "\n"
    stdout, _ = capped_gp(program, timeout=timeout, stack_bytes=stack_bytes)
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    return {
        "minimal_model": [int(value) for value in gp_block(lines, "MODEL")],
        "conductor": str(int(gp_block(lines, "CONDUCTOR")[0])),
        "log_conductor": gp_block(lines, "LOG_CONDUCTOR")[0],
        "minimal_discriminant": str(int(gp_block(lines, "DISCRIMINANT")[0])),
        "root_number": int(gp_block(lines, "ROOT_NUMBER")[0]),
    }


def quartic_value(coefficients: Sequence[Fraction], x_value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * Q(x_value) + Q(coefficient)
    return answer


def primitive_visible_points(
    construction: SixRootMestreConstruction, parameter: Fraction
) -> tuple[tuple[Fraction, Fraction], ...]:
    scale = construction.quartic_square_scale
    coefficients = construction.primitive_quartic_coefficients(parameter)
    points = tuple(
        (x_value, y_value / scale)
        for x_value, y_value in construction.visible_points(parameter)
    )
    if any(y_value**2 != quartic_value(coefficients, x_value) for x_value, y_value in points):
        raise AssertionError("a primitive visible point missed the exact quartic")
    return points


def quartic_point_to_jacobian(
    construction: SixRootMestreConstruction,
    parameter: Fraction,
    point: tuple[Fraction, Fraction],
) -> tuple[Fraction, Fraction]:
    coefficients = construction.primitive_quartic_coefficients(parameter)
    e, d, c, b, a = coefficients
    x_value, z_value = (Q(value) for value in point)
    if z_value == 0 or z_value**2 != quartic_value(coefficients, x_value):
        raise ValueError("the affine quartic point is invalid for the covariant map")
    g0 = b**2 / Q(16) - a * c / Q(6)
    g1 = b * c / Q(12) - a * d / Q(2)
    g2 = c**2 / Q(12) - b * d / Q(8) - a * e
    g3 = c * d / Q(12) - b * e / Q(2)
    g4 = d**2 / Q(16) - c * e / Q(6)
    u_x = 4 * a * x_value**3 + 3 * b * x_value**2 + 2 * c * x_value + d
    u_y = b * x_value**3 + 2 * c * x_value**2 + 3 * d * x_value + 4 * e
    g_value = (
        g0 * x_value**4
        + g1 * x_value**3
        + g2 * x_value**2
        + g3 * x_value
        + g4
    )
    g_x = 4 * g0 * x_value**3 + 3 * g1 * x_value**2 + 2 * g2 * x_value + g3
    g_y = g1 * x_value**3 + 2 * g2 * x_value**2 + 3 * g3 * x_value + 4 * g4
    h_value = (u_x * g_y - u_y * g_x) / Q(8)
    answer = (36 * g_value / z_value**2, 108 * h_value / z_value**3)
    _, _, _, coefficient_a, coefficient_b = construction.primitive_jacobian_coefficients(
        parameter
    )
    if answer[1] ** 2 != answer[0] ** 3 + coefficient_a * answer[0] + coefficient_b:
        raise AssertionError("the exact binary-quartic covariant identity failed")
    return answer


def point_on_short_curve(
    coefficients: Sequence[Fraction], point: tuple[Fraction, Fraction]
) -> bool:
    if len(coefficients) != 5 or any(Q(value) for value in coefficients[:3]):
        return False
    x_value, y_value = (Q(value) for value in point)
    return y_value**2 == x_value**3 + Q(coefficients[3]) * x_value + Q(coefficients[4])


def parse_vecsmall(text: str) -> list[int]:
    match = re.search(r"Vecsmall\(\[(.*?)\]\)", text)
    if match is None:
        raise AssertionError("PARI omitted the height-subset index vector")
    return [int(value) for value in match.group(1).split(",") if value]


def height_matrix_replay(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    precisions: tuple[int, ...],
    timeout: float,
    stack_bytes: int,
) -> tuple[dict[str, Any], ...]:
    if not points or any(not point_on_short_curve(coefficients, point) for point in points):
        raise ValueError("height replay requires exact short-Weierstrass points")
    curve = ",".join(gp_rational(Q(value)) for value in coefficients)
    point_vector = ",".join(gp_vector(point) for point in points)
    commands = [f"E=ellinit([{curve}]);", f"P=[{point_vector}];"]
    for precision in precisions:
        commands.extend(
            (
                f"default(realprecision,{precision});",
                "H=ellheightmatrix(E,P);IX=matindexrank(H);K=vecextract(P,IX[2]);",
                "HK=ellheightmatrix(E,K);EV=mateigen(HK,1)[1];",
                f'print("HEIGHT_{precision}_BEGIN");',
                "print(matrank(H));print(IX[2]);print(matdet(HK));",
                "print(vecmin(EV));print(vecmax(EV));",
                f'print("HEIGHT_{precision}_END");',
            )
        )
    commands.append("quit")
    stdout, _ = capped_gp(
        "\n".join(commands) + "\n", timeout=timeout, stack_bytes=stack_bytes
    )
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    records = []
    for precision in precisions:
        values = gp_block(lines, f"HEIGHT_{precision}")
        records.append(
            {
                "decimal_precision": precision,
                "numerical_rank": int(values[0]),
                "subset_indices_one_based": parse_vecsmall(values[1]),
                "subset_height_determinant": values[2],
                "subset_smallest_eigenvalue": values[3],
                "subset_largest_eigenvalue": values[4],
            }
        )
    ranks = {record["numerical_rank"] for record in records}
    subsets = {tuple(record["subset_indices_one_based"]) for record in records}
    if len(ranks) != 1 or len(subsets) != 1:
        raise AssertionError("the numerical height signal changed with precision")
    return tuple(records)


def bounded_quartic_points(
    coefficients: Sequence[Fraction],
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[Fraction, Fraction], ...]:
    polynomial = "+".join(
        f"{gp_rational(Q(coefficient))}*x^{power}"
        for power, coefficient in enumerate(coefficients)
    )
    program = (
        f"Q={polynomial};R=hyperellratpoints(Q,{height_bound});"
        'print("POINTS_BEGIN");print(R);print("POINTS_END");quit\n'
    )
    stdout, _ = capped_gp(program, timeout=timeout, stack_bytes=stack_bytes)
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    text = " ".join(gp_block(lines, "POINTS"))
    return tuple(
        (Q(x_value), Q(y_value))
        for x_value, y_value in re.findall(
            r"\[(-?\d+(?:/\d+)?),\s*(-?\d+(?:/\d+)?)\]", text
        )
    )


def ellrank_probe(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    curve = ",".join(gp_rational(Q(value)) for value in coefficients)
    point_vector = ",".join(gp_vector(point) for point in points)
    program = "\n".join(
        (
            f"E=ellinit([{curve}]);P=[{point_vector}];R=ellrank(E,0,P);",
            'print("RANK_BEGIN");print(R[1]);print(R[2]);print(#R[4]);print("RANK_END");',
            "quit",
        )
    ) + "\n"
    stdout, _ = capped_gp(program, timeout=timeout, stack_bytes=stack_bytes)
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    values = gp_block(lines, "RANK")
    return {
        "status": "completed PARI effort-zero computation; not a portable rank proof",
        "lower_bound": int(values[0]),
        "upper_bound": int(values[1]),
        "returned_independent_points": int(values[2]),
        "effort": 0,
    }


def point_record(point: tuple[Fraction, Fraction]) -> dict[str, str]:
    return {
        "x": rational_to_string(point[0]),
        "y": rational_to_string(point[1]),
    }


def finite_reduction_attempt(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    prime_bound: int,
) -> dict[str, Any]:
    signatures = find_mod2_reduction_certificate(
        coefficients, points, prime_bound=prime_bound
    )
    exact_rank = combined_mod2_rank(signatures, len(points))
    certified = exact_rank == len(points)
    return {
        "status": "certified" if certified else "bounded-search-rank-deficient",
        "point_count": len(points),
        "point_sha256": point_digest(points),
        "certificate_prime_bound": prime_bound,
        "certificate_primes": [signature.prime for signature in signatures],
        "combined_exact_rank_over_F2": exact_rank,
        "two_torsion_certificate_prime": (
            find_two_torsion_certificate_prime(coefficients, prime_bound=200)
            if certified
            else None
        ),
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


def numerical_subset(
    points: Sequence[tuple[Fraction, Fraction]], records: Sequence[dict[str, Any]]
) -> tuple[tuple[Fraction, Fraction], ...]:
    return tuple(
        points[index - 1] for index in records[-1]["subset_indices_one_based"]
    )


def canonical_signless_points(
    points: Iterable[tuple[Fraction, Fraction]],
) -> tuple[tuple[Fraction, Fraction], ...]:
    by_x: dict[Fraction, Fraction] = {}
    for x_value, y_value in points:
        if y_value != 0:
            by_x[Q(x_value)] = abs(Q(y_value))
    return tuple(
        (x_value, by_x[x_value])
        for x_value in sorted(
            by_x,
            key=lambda value: (
                max(abs(value.numerator), value.denominator),
                value.denominator,
                value.numerator,
            ),
        )
    )


def pari_version_capped() -> str:
    stdout, _ = capped_gp("print(version());quit\n", timeout=5, stack_bytes=8_000_000)
    return stdout.strip()


def run_specialization_screen(
    families: Sequence[tuple[int, ...]],
    *,
    parameter_bound: int,
    point_keep: int,
    height_bound: int,
    escalation_height_bound: int,
    max_search_abscissas: int,
    conductor_timeout: float,
    point_timeout: float,
    escalation_timeout: float,
    height_timeout: float,
    ellrank_timeout: float,
    stack_bytes: int,
    certificate_prime_bound: int,
) -> dict[str, Any]:
    """Run the conductor-closed phase, then exact-point/rank triage."""

    constructions = {
        roots: SixRootMestreConstruction(tuple(Q(root) for root in roots))
        for roots in families
    }
    records: list[dict[str, Any]] = []
    excluded: list[dict[str, Any]] = []
    runtime: dict[str, dict[str, Any]] = {}

    # Phase one is conductor-only.  No point, height, or rank information is
    # computed until this complete predeclared population has closed.
    for roots in families:
        construction = constructions[roots]
        for integer_parameter in range(1, parameter_bound + 1):
            parameter = Q(integer_parameter)
            discriminant = construction.quartic_discriminant(parameter)
            if discriminant == 0:
                excluded.append(
                    {
                        "roots": list(roots),
                        "parameter": integer_parameter,
                        "reason": "singular specialized quartic",
                    }
                )
                continue
            degeneracy = construction.visible_point_degeneracy(parameter)
            if degeneracy.collision_loss or degeneracy.zero_ordinates:
                excluded.append(
                    {
                        "roots": list(roots),
                        "parameter": integer_parameter,
                        "reason": "visible-point collision or zero ordinate",
                        "collision_loss": degeneracy.collision_loss,
                        "zero_ordinates": degeneracy.zero_ordinates,
                    }
                )
                continue
            identifier = "r" + "_".join(str(value) for value in roots) + f"_t{integer_parameter}"
            coefficients = construction.primitive_jacobian_coefficients(parameter)
            record: dict[str, Any] = {
                "identifier": identifier,
                "roots": list(roots),
                "parameter": integer_parameter,
                "known_nagao_calibration_tuple": roots == NAGAO_NORMALIZED_ROOTS,
                "admissibility": {
                    "exact_nonzero_quartic_discriminant": True,
                    "visible_collision_loss": 0,
                    "visible_zero_ordinates": 0,
                },
            }
            try:
                conductor = capped_minimal_curve_data(
                    coefficients,
                    timeout=conductor_timeout,
                    stack_bytes=stack_bytes,
                )
                record["conductor_phase"] = {
                    "status": "completed exact PARI minimal-model/conductor computation",
                    **conductor,
                    "below_strict_log_conductor_target_numerically": (
                        Decimal(conductor["log_conductor"]) < TARGET_LOG_CONDUCTOR
                    ),
                }
                runtime[identifier] = {
                    "roots": roots,
                    "construction": construction,
                    "parameter": parameter,
                    "coefficients": coefficients,
                }
            except CappedProcessTimeout:
                record["conductor_phase"] = {
                    "status": "timeout",
                    "timeout_seconds": conductor_timeout,
                }
            except Exception as error:
                record["conductor_phase"] = {
                    "status": "error",
                    "error": str(error)[:1000],
                }
            records.append(record)

    conductor_population_closed = True

    # Phase two maps the exact displayed points and obtains numerical signals
    # for every successfully conducted member, without retroactive filtering.
    for record in records:
        identifier = record["identifier"]
        if identifier not in runtime:
            continue
        item = runtime[identifier]
        construction = item["construction"]
        parameter = item["parameter"]
        coefficients = item["coefficients"]
        try:
            quartic_points = primitive_visible_points(construction, parameter)
            jacobian_points = tuple(
                quartic_point_to_jacobian(construction, parameter, point)
                for point in quartic_points
            )
            if len(quartic_points) != 12 or len({point[0] for point in quartic_points}) != 12:
                raise AssertionError("an admissible fiber lost a visible abscissa")
            height = height_matrix_replay(
                coefficients,
                jacobian_points,
                precisions=(72, 120),
                timeout=height_timeout,
                stack_bytes=stack_bytes,
            )
            stable_rank = int(height[-1]["numerical_rank"])
            subset = numerical_subset(jacobian_points, height)
            record["visible_point_triage"] = {
                "status": "completed exact point checks and numerical height triage",
                "quartic_point_count": len(quartic_points),
                "distinct_quartic_abscissas": len({point[0] for point in quartic_points}),
                "quartic_point_sha256": point_digest(quartic_points),
                "jacobian_point_count": len(jacobian_points),
                "jacobian_point_sha256": point_digest(jacobian_points),
                "exact_quartic_and_jacobian_membership_checked": True,
                "height_matrix_runs": list(height),
                "stable_numerical_rank": stable_rank,
                "numerical_subset": [point_record(point) for point in subset],
                "numerical_rank_is_not_an_independence_certificate": True,
            }
            item.update(
                {
                    "quartic_points": quartic_points,
                    "jacobian_points": jacobian_points,
                    "visible_height": height,
                }
            )
        except CappedProcessTimeout:
            record["visible_point_triage"] = {
                "status": "timeout",
                "timeout_seconds": height_timeout,
            }
        except Exception as error:
            record["visible_point_triage"] = {
                "status": "error",
                "error": str(error)[:1000],
            }

    ranked = [
        record
        for record in records
        if record.get("visible_point_triage", {}).get("status", "").startswith("completed")
    ]
    ranked.sort(
        key=lambda record: (
            -int(record["visible_point_triage"]["stable_numerical_rank"]),
            Decimal(record["conductor_phase"]["log_conductor"]),
            tuple(record["roots"]),
            int(record["parameter"]),
        )
    )
    finalists = ranked[:point_keep]
    finalist_records: list[dict[str, Any]] = []
    target_hits: list[dict[str, Any]] = []
    for selection_position, selected in enumerate(finalists, 1):
        identifier = selected["identifier"]
        item = runtime[identifier]
        construction = item["construction"]
        parameter = item["parameter"]
        coefficients = item["coefficients"]
        visible_quartic = item["quartic_points"]
        visible_jacobian = item["jacobian_points"]
        finalist: dict[str, Any] = {
            "selection_position": selection_position,
            "identifier": identifier,
            "roots": selected["roots"],
            "parameter": selected["parameter"],
            "selection_rule": (
                "descending stable visible-point numerical rank, then ascending "
                "exactly computed log conductor, roots, and parameter"
            ),
            "input_visible_stable_numerical_rank": selected["visible_point_triage"][
                "stable_numerical_rank"
            ],
            "conductor": selected["conductor_phase"],
        }
        try:
            raw_points = bounded_quartic_points(
                construction.primitive_quartic_coefficients(parameter),
                height_bound=height_bound,
                timeout=point_timeout,
                stack_bytes=stack_bytes,
            )
            signless = canonical_signless_points(raw_points)
            if any(
                point[1] ** 2
                != quartic_value(
                    construction.primitive_quartic_coefficients(parameter), point[0]
                )
                for point in signless
            ):
                raise AssertionError("PARI returned a point off the exact quartic")
            selected_signless = signless[:max_search_abscissas]
            searched_jacobian = tuple(
                quartic_point_to_jacobian(construction, parameter, point)
                for point in selected_signless
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
                stack_bytes=stack_bytes,
            )
            stable_rank = int(height[-1]["numerical_rank"])
            subset = numerical_subset(pool, height)
            try:
                ellrank = ellrank_probe(
                    coefficients,
                    subset,
                    timeout=ellrank_timeout,
                    stack_bytes=stack_bytes,
                )
            except CappedProcessTimeout:
                ellrank = {
                    "status": "timeout",
                    "timeout_seconds": ellrank_timeout,
                }
            except Exception as error:
                ellrank = {"status": "error", "error": str(error)[:1000]}
            finalist.update(
                {
                    "point_search": {
                        "status": "complete bounded PARI hyperellratpoints enumeration",
                        "height_bound": height_bound,
                        "signed_points_returned": len(raw_points),
                        "distinct_nonzero_ordinate_abscissas": len(signless),
                        "abscissas_retained_for_mapping": len(selected_signless),
                        "mapping_cap": max_search_abscissas,
                        "mapping_truncated": len(signless) > len(selected_signless),
                        "visible_abscissas_returned": sum(
                            point[0] in {visible[0] for visible in visible_quartic}
                            for point in signless
                        ),
                        "all_retained_points_checked_exactly": True,
                    },
                    "augmented_rank_triage": {
                        "pool_point_count_modulo_inverse": len(pool),
                        "pool_point_sha256": point_digest(pool),
                        "height_matrix_runs": list(height),
                        "stable_numerical_rank": stable_rank,
                        "numerical_subset": [point_record(point) for point in subset],
                        "effort_zero_ellrank": ellrank,
                        "numerical_rank_is_not_an_independence_certificate": True,
                    },
                }
            )
            if stable_rank >= 21:
                certificate = finite_reduction_attempt(
                    coefficients, subset, prime_bound=certificate_prime_bound
                )
                finalist["finite_reduction_attempt"] = certificate
                certified_rank = certificate["certified_algebraic_rank_lower_bound"]
                below_target = Decimal(selected["conductor_phase"]["log_conductor"]) < TARGET_LOG_CONDUCTOR
                if certified_rank is not None and (
                    certified_rank >= 30 or (certified_rank >= 21 and below_target)
                ):
                    target_hits.append(
                        {
                            "identifier": identifier,
                            "certified_rank_lower_bound": certified_rank,
                            "conductor": selected["conductor_phase"]["conductor"],
                            "log_conductor": selected["conductor_phase"]["log_conductor"],
                        }
                    )
            else:
                finalist["finite_reduction_attempt"] = {
                    "status": "not triggered",
                    "trigger_stable_numerical_rank": 21,
                }
        except CappedProcessTimeout:
            finalist["point_search"] = {
                "status": "timeout",
                "timeout_seconds": point_timeout,
            }
        except Exception as error:
            finalist["point_search"] = {
                "status": "error",
                "error": str(error)[:1000],
            }
        finalist_records.append(finalist)

    # One deterministic deeper replay follows the strongest augmented signal.
    # This is selected only after every shallow finalist has closed, and hence
    # cannot alter the conductor or shallow-search populations.
    escalation_candidates = [
        record for record in finalist_records if "augmented_rank_triage" in record
    ]
    escalation_candidates.sort(
        key=lambda record: (
            -int(record["augmented_rank_triage"]["stable_numerical_rank"]),
            Decimal(record["conductor"]["log_conductor"]),
            tuple(record["roots"]),
            int(record["parameter"]),
        )
    )
    if escalation_height_bound > height_bound and escalation_candidates:
        escalated = escalation_candidates[0]
        item = runtime[escalated["identifier"]]
        construction = item["construction"]
        parameter = item["parameter"]
        coefficients = item["coefficients"]
        visible_jacobian = item["jacobian_points"]
        try:
            raw_points = bounded_quartic_points(
                construction.primitive_quartic_coefficients(parameter),
                height_bound=escalation_height_bound,
                timeout=escalation_timeout,
                stack_bytes=stack_bytes,
            )
            signless = canonical_signless_points(raw_points)
            selected_signless = signless[:max_search_abscissas]
            searched_jacobian = tuple(
                quartic_point_to_jacobian(construction, parameter, point)
                for point in selected_signless
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
                stack_bytes=stack_bytes,
            )
            stable_rank = int(height[-1]["numerical_rank"])
            subset = numerical_subset(pool, height)
            try:
                ellrank = ellrank_probe(
                    coefficients,
                    subset,
                    timeout=ellrank_timeout,
                    stack_bytes=stack_bytes,
                )
            except CappedProcessTimeout:
                ellrank = {
                    "status": "timeout",
                    "timeout_seconds": ellrank_timeout,
                }
            except Exception as error:
                ellrank = {"status": "error", "error": str(error)[:1000]}
            escalated["single_strongest_signal_escalation"] = {
                "selection_rule": (
                    "maximum shallow augmented stable numerical rank, then "
                    "ascending precomputed conductor, roots, and parameter"
                ),
                "status": "complete bounded PARI hyperellratpoints enumeration",
                "height_bound": escalation_height_bound,
                "signed_points_returned": len(raw_points),
                "distinct_nonzero_ordinate_abscissas": len(signless),
                "abscissas_retained_for_mapping": len(selected_signless),
                "mapping_cap": max_search_abscissas,
                "mapping_truncated": len(signless) > len(selected_signless),
                "pool_point_count_modulo_inverse": len(pool),
                "pool_point_sha256": point_digest(pool),
                "height_matrix_runs": list(height),
                "stable_numerical_rank": stable_rank,
                "numerical_subset": [point_record(point) for point in subset],
                "effort_zero_ellrank": ellrank,
                "numerical_rank_is_not_an_independence_certificate": True,
            }
            if stable_rank >= 21:
                certificate = finite_reduction_attempt(
                    coefficients, subset, prime_bound=certificate_prime_bound
                )
                escalated["single_strongest_signal_escalation"][
                    "finite_reduction_attempt"
                ] = certificate
                certified_rank = certificate["certified_algebraic_rank_lower_bound"]
                below_target = Decimal(escalated["conductor"]["log_conductor"]) < TARGET_LOG_CONDUCTOR
                if certified_rank is not None and (
                    certified_rank >= 30 or (certified_rank >= 21 and below_target)
                ):
                    target_hits.append(
                        {
                            "identifier": escalated["identifier"],
                            "certified_rank_lower_bound": certified_rank,
                            "conductor": escalated["conductor"]["conductor"],
                            "log_conductor": escalated["conductor"]["log_conductor"],
                            "source": "single strongest-signal escalation",
                        }
                    )
            else:
                escalated["single_strongest_signal_escalation"][
                    "finite_reduction_attempt"
                ] = {
                    "status": "not triggered",
                    "trigger_stable_numerical_rank": 21,
                }
        except CappedProcessTimeout:
            escalated["single_strongest_signal_escalation"] = {
                "status": "timeout",
                "timeout_seconds": escalation_timeout,
            }
        except Exception as error:
            escalated["single_strongest_signal_escalation"] = {
                "status": "error",
                "error": str(error)[:1000],
            }

    rank_histogram: dict[str, int] = {}
    for record in ranked:
        key = str(record["visible_point_triage"]["stable_numerical_rank"])
        rank_histogram[key] = rank_histogram.get(key, 0) + 1
    completed_augmented = [
        record
        for record in finalist_records
        if "augmented_rank_triage" in record
    ]
    return {
        "protocol": {
            "conductor_population_closed_before_point_or_rank_triage": conductor_population_closed,
            "integer_parameters": [1, parameter_bound],
            "conductor_timeout_seconds_per_fiber": conductor_timeout,
            "height_precisions": [72, 120],
            "height_timeout_seconds_per_fiber": height_timeout,
            "point_search_height_bound": height_bound,
            "point_search_timeout_seconds_per_finalist": point_timeout,
            "single_strongest_signal_escalation_height_bound": escalation_height_bound,
            "single_strongest_signal_escalation_timeout_seconds": escalation_timeout,
            "point_mapping_abscissa_cap": max_search_abscissas,
            "effort_zero_ellrank_timeout_seconds_per_finalist": ellrank_timeout,
            "finite_reduction_trigger_rank": 21,
            "finite_reduction_prime_bound": certificate_prime_bound,
        },
        "population": {
            "families": len(families),
            "proposed_integer_fibers": len(families) * parameter_bound,
            "inadmissible_fibers": len(excluded),
            "admissible_fibers": len(records),
            "conductor_completed": sum(
                record["conductor_phase"]["status"].startswith("completed")
                for record in records
            ),
            "conductor_timeouts": sum(
                record["conductor_phase"]["status"] == "timeout" for record in records
            ),
            "conductor_errors": sum(
                record["conductor_phase"]["status"] == "error" for record in records
            ),
            "visible_triage_completed": len(ranked),
            "visible_rank_histogram": rank_histogram,
            "maximum_visible_stable_numerical_rank": max(
                (
                    int(record["visible_point_triage"]["stable_numerical_rank"])
                    for record in ranked
                ),
                default=None,
            ),
        },
        "inadmissible_fibers": excluded,
        "conductor_first_fiber_records": records,
        "point_search_finalists": finalist_records,
        "maximum_augmented_stable_numerical_rank": max(
            (
                int(record["augmented_rank_triage"]["stable_numerical_rank"])
                for record in completed_augmented
            ),
            default=None,
        ),
        "maximum_escalated_stable_numerical_rank": max(
            (
                int(record["single_strongest_signal_escalation"]["stable_numerical_rank"])
                for record in finalist_records
                if "stable_numerical_rank"
                in record.get("single_strongest_signal_escalation", {})
            ),
            default=None,
        ),
        "target_hits": target_hits,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-root", type=int, default=50)
    parser.add_argument("--prior-max-root", type=int, default=14)
    parser.add_argument("--parameter-bound", type=int, default=8)
    parser.add_argument("--point-keep", type=int, default=25)
    parser.add_argument("--height-bound", type=int, default=5_000)
    parser.add_argument("--escalation-height-bound", type=int, default=50_000)
    parser.add_argument("--max-search-abscissas", type=int, default=256)
    parser.add_argument("--compiler", default="c++")
    parser.add_argument("--compile-timeout", type=float, default=30.0)
    parser.add_argument("--enumeration-timeout", type=float, default=30.0)
    parser.add_argument("--conductor-timeout", type=float, default=5.0)
    parser.add_argument("--point-timeout", type=float, default=10.0)
    parser.add_argument("--escalation-timeout", type=float, default=20.0)
    parser.add_argument("--height-timeout", type=float, default=10.0)
    parser.add_argument("--ellrank-timeout", type=float, default=10.0)
    parser.add_argument("--stack-bytes", type=int, default=256_000_000)
    parser.add_argument("--certificate-prime-bound", type=int, default=500)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.max_root < 34 or args.max_root > MAX_SCRIPT_ROOT:
        raise SystemExit(f"--max-root must lie in [34,{MAX_SCRIPT_ROOT}]")
    if args.prior_max_root < 5 or args.prior_max_root >= args.max_root:
        raise SystemExit("--prior-max-root must lie below --max-root")
    if args.parameter_bound < 1 or args.parameter_bound > MAX_PARAMETER_BOUND:
        raise SystemExit(f"--parameter-bound must lie in [1,{MAX_PARAMETER_BOUND}]")
    if args.point_keep < 1 or args.point_keep > MAX_POINT_KEEP:
        raise SystemExit(f"--point-keep must lie in [1,{MAX_POINT_KEEP}]")
    if args.height_bound < 1 or args.height_bound > MAX_HEIGHT_BOUND:
        raise SystemExit(f"--height-bound must lie in [1,{MAX_HEIGHT_BOUND}]")
    if (
        args.escalation_height_bound < args.height_bound
        or args.escalation_height_bound > MAX_HEIGHT_BOUND
    ):
        raise SystemExit(
            "--escalation-height-bound must be at least --height-bound and at "
            f"most {MAX_HEIGHT_BOUND}"
        )
    if args.max_search_abscissas < 12 or args.max_search_abscissas > 512:
        raise SystemExit("--max-search-abscissas must lie in [12,512]")
    timeouts = (
        args.compile_timeout,
        args.enumeration_timeout,
        args.conductor_timeout,
        args.point_timeout,
        args.escalation_timeout,
        args.height_timeout,
        args.ellrank_timeout,
    )
    if min(timeouts) <= 0 or max(timeouts) > MAX_SUBPROCESS_TIMEOUT:
        raise SystemExit(
            f"all subprocess timeouts must lie in (0,{MAX_SUBPROCESS_TIMEOUT:g}]"
        )
    if args.stack_bytes < 8_000_000 or args.stack_bytes > 1_000_000_000:
        raise SystemExit("--stack-bytes must lie in [8000000,1000000000]")
    if args.certificate_prime_bound < 3 or args.certificate_prime_bound > 2000:
        raise SystemExit("--certificate-prime-bound must lie in [3,2000]")


def main() -> None:
    args = build_parser().parse_args()
    validate_args(args)
    script_path = Path(__file__).resolve()
    cpp_path = script_path.with_name("enumerate_mestre_root_tuples_scale.cpp")
    repo_root = script_path.parents[2]
    output = args.output if args.output.is_absolute() else repo_root / args.output

    enumeration = compiled_enumeration(
        args.max_root,
        compiler=args.compiler,
        compile_timeout=args.compile_timeout,
        enumeration_timeout=args.enumeration_timeout,
    )
    verify_enumerator_records(enumeration)
    calibration_enumeration = compiled_enumeration(
        14,
        compiler=args.compiler,
        compile_timeout=args.compile_timeout,
        enumeration_timeout=args.enumeration_timeout,
    )
    verify_enumerator_records(calibration_enumeration)
    nonsingular, singular, witnesses = classify_nonreflection(enumeration)
    calibration = max14_calibration(enumeration, calibration_enumeration)
    expected_calibration = {
        "affine_normalized_root_tuples": 1023,
        "obstruction_zero": 68,
        "generically_nonsingular": 59,
        "generically_nonsingular_nonreflection": 2,
    }
    for key, value in expected_calibration.items():
        if calibration[key] != value:
            raise AssertionError(f"max-root-14 calibration changed at {key}")
    if NAGAO_NORMALIZED_ROOTS not in nonsingular:
        raise AssertionError("the scale-up missed Nagao's normalized root tuple")

    screen_families = tuple(
        roots for roots in nonsingular if roots[-1] > args.prior_max_root
    )
    specialization = run_specialization_screen(
        screen_families,
        parameter_bound=args.parameter_bound,
        point_keep=args.point_keep,
        height_bound=args.height_bound,
        escalation_height_bound=args.escalation_height_bound,
        max_search_abscissas=args.max_search_abscissas,
        conductor_timeout=args.conductor_timeout,
        point_timeout=args.point_timeout,
        escalation_timeout=args.escalation_timeout,
        height_timeout=args.height_timeout,
        ellrank_timeout=args.ellrank_timeout,
        stack_bytes=args.stack_bytes,
        certificate_prime_bound=args.certificate_prime_bound,
    )

    target_hits = specialization["target_hits"]
    artifact: dict[str, Any] = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": (
            "complete bounded exact max-root enumeration and conductor-first "
            "specialization screen; numerical height ranks are triage only"
        ),
        "target": {
            "rank_at_least": 21,
            "log_conductor_strict_upper_bound": "182.72",
            "alternative_rank_at_least": 30,
            "hits": target_hits,
            "explanation": (
                "no exact finite-reduction rank certificate reached either target"
                if not target_hits
                else "at least one exact finite-reduction target certificate was produced"
            ),
        },
        "method": {
            "enumeration": (
                "compiled exhaustive nested enumeration with signed 128-bit exact "
                "arithmetic; primitive translation/scale normalization and a "
                "lexicographic reflection quotient"
            ),
            "quartic_filter": (
                "exact degree-five Mestre obstruction, independently replayed in Python"
            ),
            "nonsingularity_filter": (
                "exact quartic discriminant at T=1,...,21; its degree is at most "
                "20, so 21 zeros prove generic singularity"
            ),
            "rank_viability_gate": (
                "reflection-symmetric tuples are excluded before specialization; "
                "passing is necessary but does not prove section independence"
            ),
            "selection_leakage_control": (
                "the full declared admissible fiber population receives conductor "
                "computations before any exact-point or numerical-rank triage begins"
            ),
            "rank_semantics": (
                "exact point membership plus precision-stable numerical height rank; "
                "only finite-reduction independence counts as an algebraic rank proof"
            ),
        },
        "parameters": {
            "max_root": args.max_root,
            "prior_max_root": args.prior_max_root,
            "parameter_bound": args.parameter_bound,
            "point_keep": args.point_keep,
            "height_bound": args.height_bound,
            "escalation_height_bound": args.escalation_height_bound,
            "max_search_abscissas": args.max_search_abscissas,
            "compile_timeout": args.compile_timeout,
            "enumeration_timeout": args.enumeration_timeout,
            "conductor_timeout": args.conductor_timeout,
            "point_timeout": args.point_timeout,
            "escalation_timeout": args.escalation_timeout,
            "height_timeout": args.height_timeout,
            "ellrank_timeout": args.ellrank_timeout,
            "stack_bytes": args.stack_bytes,
            "certificate_prime_bound": args.certificate_prime_bound,
            "output": str(args.output),
        },
        "enumeration": {
            "affine_normalized_primitive_reflection_quotient_count": enumeration.normalized_count,
            "degree_five_obstruction_zero_count": enumeration.obstruction_count,
            "reflection_obstruction_zero_count": enumeration.reflection_count,
            "nonreflection_obstruction_zero_count": enumeration.nonreflection_count,
            "nonreflection_generically_nonsingular_count": len(nonsingular),
            "nonreflection_generically_singular_count": len(singular),
            "obstruction_tuple_sha256": tuple_digest(enumeration.obstruction_roots),
            "nonreflection_tuple_sha256": tuple_digest(enumeration.nonreflection_roots),
            "nonsingular_nonreflection_tuple_sha256": tuple_digest(nonsingular),
            "generically_nonsingular_nonreflection_tuples": [
                {
                    "roots": list(roots),
                    "discriminant_witness_parameter": witnesses[roots],
                    "known_nagao_calibration_tuple": roots == NAGAO_NORMALIZED_ROOTS,
                }
                for roots in nonsingular
            ],
            "generically_singular_nonreflection_tuples": [
                list(roots) for roots in singular
            ],
            "max_root_14_calibration": calibration,
            "nagao_calibration": {
                "source_tuple": [-17, -16, 10, 11, 14, 17],
                "normalized_tuple": list(NAGAO_NORMALIZED_ROOTS),
                "recovered": True,
                "diameter": 34,
            },
            "families_beyond_prior_bound": len(screen_families),
            "new_families_excluding_known_nagao_calibration": (
                len(screen_families) - int(NAGAO_NORMALIZED_ROOTS in screen_families)
            ),
        },
        "specialization_screen": specialization,
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "pari_gp": pari_version_capped(),
            "compiler": shutil.which(args.compiler),
        },
        "provenance": {
            "script": str(script_path.relative_to(repo_root)),
            "script_sha256": sha256_file(script_path),
            "compiled_source": str(cpp_path.relative_to(repo_root)),
            "compiled_source_sha256": sha256_file(cpp_path),
            "reproducing_command": REPRODUCING_COMMAND,
            "temporary_binary_removed_after_each_enumeration": True,
            "subprocesses_run_in_foreground_process_groups": True,
            "whole_process_group_killed_and_reaped_on_timeout": True,
        },
    }
    digest_payload = {
        "enumeration": artifact["enumeration"],
        "specialization_screen": artifact["specialization_screen"],
        "target": artifact["target"],
    }
    artifact["result_sha256"] = hashlib.sha256(
        json.dumps(digest_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        f"normalized={enumeration.normalized_count} obstruction={enumeration.obstruction_count} "
        f"nonreflection={enumeration.nonreflection_count} nonsingular={len(nonsingular)}"
    )
    print(
        f"fibers={specialization['population']['admissible_fibers']} "
        f"max_visible_rank={specialization['population']['maximum_visible_stable_numerical_rank']} "
        f"max_augmented_rank={specialization['maximum_augmented_stable_numerical_rank']} "
        f"max_escalated_rank={specialization['maximum_escalated_stable_numerical_rank']} "
        f"hits={len(target_hits)}"
    )
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
