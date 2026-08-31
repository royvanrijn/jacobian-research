#!/usr/bin/env python3
"""Exact specialization and baseline certification for q12o5867.

The rootless family is stored as a short Weierstrass equation with weights
``(A,B,x,y)=(8,12,4,6)`` on the projective parameter line.  This module keeps
those weights explicit, obtains PARI's exact global-minimal change, transports
all seventeen polynomial sections, and delegates independence to the common
finite-quotient certificate implementation.

This is a specialization adapter, not a point search or a rank upper bound.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, lcm
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Iterable, Sequence


ELLIPTIC_ROOT = Path(__file__).resolve().parents[1]
CAS = ELLIPTIC_ROOT / "cas"
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))

from elliptic_candidate_record import (  # noqa: E402
    WeierstrassChange,
    build_finite_quotient_certificate,
    change_weierstrass_model,
    is_on_weierstrass_curve,
    source_point_to_target,
    verify_finite_quotient_certificate,
    weierstrass_invariants,
)
from pari_bridge import pari_version  # noqa: E402


Q = Fraction
RationalPoint = tuple[Fraction, Fraction]
WeierstrassModel = tuple[Fraction, Fraction, Fraction, Fraction, Fraction]

MODEL_STATUS = "PASS_EXACT_QQ_Q12O5867_SMOOTH_RR_ROOTLESS_JACOBIAN"
SECTION_STATUS = "PASS_EXACT_QQ_Q12O5867_ROOTLESS_17_SELECTED_SECTIONS"


def fraction_text(value: Fraction | int) -> str:
    value_q = Q(value)
    if value_q.denominator == 1:
        return str(value_q.numerator)
    return f"{value_q.numerator}/{value_q.denominator}"


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_projective_parameter(a: int, b: int) -> tuple[int, int]:
    """Return the unique primitive representative with ``b >= 0``."""

    if a == 0 and b == 0:
        raise ValueError("(0:0) is not a projective parameter")
    common = gcd(abs(a), abs(b))
    a //= common
    b //= common
    if b < 0 or (b == 0 and a < 0):
        a, b = -a, -b
    return a, b


def homogeneous_value(
    coefficients_low_to_high: Sequence[Fraction | int | str],
    a: int,
    b: int,
    weight: int,
) -> Fraction:
    """Evaluate ``sum c_i*a^i*b^(weight-i)`` exactly."""

    if len(coefficients_low_to_high) > weight + 1:
        raise ValueError("a polynomial exceeds its declared homogeneous weight")
    return sum(
        (
            Q(coefficient)
            * a**index
            * b ** (weight - index)
        )
        for index, coefficient in enumerate(coefficients_low_to_high)
    )


@dataclass(frozen=True)
class Q12O5867Data:
    model_path: Path
    sections_path: Path
    a_coefficients: tuple[Fraction, ...]
    b_coefficients: tuple[Fraction, ...]
    section_coefficients: tuple[tuple[tuple[Fraction, ...], tuple[Fraction, ...]], ...]
    model_sha256: str
    sections_sha256: str


def load_q12o5867_data(model_path: Path, sections_path: Path) -> Q12O5867Data:
    model_record = json.loads(model_path.read_text())
    section_record = json.loads(sections_path.read_text())
    if model_record.get("status") != MODEL_STATUS:
        raise ValueError("the q12o5867 model does not have its expected exact status")
    if section_record.get("status") != SECTION_STATUS:
        raise ValueError("the q12o5867 sections do not have their expected exact status")
    child = model_record["child"]
    a_coefficients = tuple(
        Q(value) for value in child["minimal_A_coefficients_low_to_high"]
    )
    b_coefficients = tuple(
        Q(value) for value in child["minimal_B_coefficients_low_to_high"]
    )
    if len(a_coefficients) != 9 or len(b_coefficients) != 13:
        raise ValueError("the q12o5867 coefficient degrees changed")
    sections = []
    for expected_index, record in enumerate(section_record["sections"]):
        if int(record["basis_index"]) != expected_index:
            raise ValueError("the ordered q12o5867 basis indices changed")
        section = record["section"]
        if not section.get("exact_weierstrass_identity"):
            raise ValueError("a selected section lacks its exact identity flag")
        x_coefficients = tuple(Q(value) for value in section["x_coefficients_low_to_high"])
        y_coefficients = tuple(Q(value) for value in section["y_coefficients_low_to_high"])
        if len(x_coefficients) > 5 or len(y_coefficients) > 7:
            raise ValueError("a selected section exceeds weights (4,6)")
        sections.append((x_coefficients, y_coefficients))
    if len(sections) != 17:
        raise ValueError("the selected q12o5867 basis no longer has 17 sections")
    return Q12O5867Data(
        model_path=model_path,
        sections_path=sections_path,
        a_coefficients=a_coefficients,
        b_coefficients=b_coefficients,
        section_coefficients=tuple(sections),
        model_sha256=sha256_file(model_path),
        sections_sha256=sha256_file(sections_path),
    )


@dataclass(frozen=True)
class ProjectiveSpecialization:
    a: int
    b: int
    coefficient_a: Fraction
    coefficient_b: Fraction
    points: tuple[RationalPoint, ...]

    @property
    def model(self) -> WeierstrassModel:
        return Q(0), Q(0), Q(0), self.coefficient_a, self.coefficient_b


def evaluate_projective_specialization(
    data: Q12O5867Data, a: int, b: int
) -> ProjectiveSpecialization:
    a, b = normalize_projective_parameter(a, b)
    coefficient_a = homogeneous_value(data.a_coefficients, a, b, 8)
    coefficient_b = homogeneous_value(data.b_coefficients, a, b, 12)
    model = (Q(0), Q(0), Q(0), coefficient_a, coefficient_b)
    if weierstrass_invariants(model)["discriminant"] == 0:
        raise ValueError(f"q12o5867 has a singular fibre at ({a}:{b})")
    points = tuple(
        (
            homogeneous_value(x_coefficients, a, b, 4),
            homogeneous_value(y_coefficients, a, b, 6),
        )
        for x_coefficients, y_coefficients in data.section_coefficients
    )
    if any(not is_on_weierstrass_curve(model, point) for point in points):
        raise AssertionError("a homogenized q12o5867 section missed the fibre")
    return ProjectiveSpecialization(a, b, coefficient_a, coefficient_b, points)


def _run_gp(program: str, *, timeout: float, stack_bytes: int) -> str:
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    try:
        completed = subprocess.run(
            [executable, "-q", "-f", "-s", str(stack_bytes)],
            input=program,
            text=True,
            capture_output=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(f"PARI/GP exceeded the exact timeout of {timeout}s") from error
    fatal = [
        line
        for line in completed.stderr.splitlines()
        if "***" in line and "Warning:" not in line
    ]
    if completed.returncode != 0 or fatal:
        detail = "\n".join(fatal) or completed.stderr.strip()
        raise RuntimeError(f"PARI/GP failed: {detail}")
    return completed.stdout


def factor_integer_with_pari(
    value: int, *, timeout: float = 120.0, stack_bytes: int = 128_000_000
) -> tuple[tuple[int, int], ...]:
    """Return and independently multiply-check PARI's complete factorization."""

    if value < 1:
        raise ValueError("factorization input must be positive")
    if value == 1:
        return ()
    program = "\n".join(
        (
            f"N={value};F=factor(N);",
            'for(i=1,matsize(F)[1],print("FACTOR|",F[i,1],"|",F[i,2]));',
            'print("FACTOR_END");',
            "quit",
        )
    ) + "\n"
    output = _run_gp(program, timeout=timeout, stack_bytes=stack_bytes)
    factors = tuple(
        (int(fields[1]), int(fields[2]))
        for line in output.splitlines()
        if (fields := line.strip().split("|"))[0] == "FACTOR"
    )
    if "FACTOR_END" not in output.splitlines():
        raise RuntimeError("PARI factorization output was truncated")
    product = 1
    for prime, exponent in factors:
        if prime < 2 or exponent < 1:
            raise AssertionError("PARI emitted an invalid factor record")
        product *= prime**exponent
    if product != value:
        raise AssertionError("PARI factorization did not multiply back exactly")
    return factors


def _valuation(value: int, prime: int) -> int:
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def integral_scaling_factor(
    coefficient_a: Fraction,
    coefficient_b: Fraction,
    *,
    timeout: float = 120.0,
    stack_bytes: int = 128_000_000,
) -> tuple[int, dict[str, Any]]:
    """Choose the least positive integral scale from denominator valuations.

    If complete factorization is unavailable, the denominator lcm is a larger
    but unconditional scale and still preserves an exact pre-minimal result.
    """

    denominator_a = coefficient_a.denominator
    denominator_b = coefficient_b.denominator
    denominator_lcm = lcm(denominator_a, denominator_b)
    try:
        factors = factor_integer_with_pari(
            denominator_lcm, timeout=timeout, stack_bytes=stack_bytes
        )
        scale = 1
        factor_records = []
        for prime, lcm_exponent in factors:
            exponent_a = _valuation(denominator_a, prime)
            exponent_b = _valuation(denominator_b, prime)
            scale_exponent = max((exponent_a + 3) // 4, (exponent_b + 5) // 6)
            scale *= prime**scale_exponent
            factor_records.append(
                {
                    "prime": str(prime),
                    "lcm_exponent": lcm_exponent,
                    "A_denominator_exponent": exponent_a,
                    "B_denominator_exponent": exponent_b,
                    "scale_exponent": scale_exponent,
                }
            )
        method: dict[str, Any] = {
            "method": "valuation-optimal-from-complete-PARI-factorization",
            "denominator_lcm_factors": factor_records,
        }
    except (FileNotFoundError, RuntimeError) as error:
        scale = denominator_lcm
        method = {
            "method": "unconditional-denominator-lcm-fallback",
            "factorization_blocker": f"{type(error).__name__}: {error}",
        }
    if (coefficient_a * scale**4).denominator != 1:
        raise AssertionError("the chosen scale did not integralize A")
    if (coefficient_b * scale**6).denominator != 1:
        raise AssertionError("the chosen scale did not integralize B")
    method.update(
        {
            "A_denominator": str(denominator_a),
            "B_denominator": str(denominator_b),
            "scale": str(scale),
            "scale_bits": scale.bit_length(),
        }
    )
    return scale, method


def integralize_specialization(
    specialization: ProjectiveSpecialization,
    scale: int,
) -> tuple[WeierstrassModel, tuple[RationalPoint, ...]]:
    model = (
        Q(0),
        Q(0),
        Q(0),
        specialization.coefficient_a * scale**4,
        specialization.coefficient_b * scale**6,
    )
    if any(value.denominator != 1 for value in model):
        raise AssertionError("the pre-minimal short model is not integral")
    points = tuple(
        (x_coordinate * scale**2, y_coordinate * scale**3)
        for x_coordinate, y_coordinate in specialization.points
    )
    if any(not is_on_weierstrass_curve(model, point) for point in points):
        raise AssertionError("an integralized q12o5867 point missed the curve")
    return model, points


def global_minimal_model_with_change(
    source_model: Sequence[Fraction | int],
    *,
    timeout: float = 300.0,
    stack_bytes: int = 512_000_000,
) -> tuple[WeierstrassModel, WeierstrassChange, dict[str, Any]]:
    """Run PARI ``ellminimalmodel`` and exactly replay its returned change."""

    if len(source_model) != 5 or any(Q(value).denominator != 1 for value in source_model):
        raise ValueError("PARI minimalization input must be an integral model")
    vector = ",".join(str(Q(value).numerator) for value in source_model)
    program = "\n".join(
        (
            f"E=ellinit([{vector}]);",
            "v=0;M=ellminimalmodel(E,&v);",
            'print("MIN_MODEL|",M.a1,"|",M.a2,"|",M.a3,"|",M.a4,"|",M.a6);',
            'print("MIN_CHANGE|",v[1],"|",v[2],"|",v[3],"|",v[4]);',
            "quit",
        )
    ) + "\n"
    output = _run_gp(program, timeout=timeout, stack_bytes=stack_bytes)
    records: dict[str, list[str]] = {}
    for line in output.splitlines():
        fields = line.strip().split("|")
        if fields[0] in {"MIN_MODEL", "MIN_CHANGE"}:
            if fields[0] in records:
                raise RuntimeError(f"PARI duplicated {fields[0]}")
            records[fields[0]] = fields[1:]
    if len(records.get("MIN_MODEL", ())) != 5 or len(records.get("MIN_CHANGE", ())) != 4:
        raise RuntimeError("PARI minimal-model output was incomplete")
    minimal_model = tuple(Q(value) for value in records["MIN_MODEL"])
    change = WeierstrassChange.from_values(records["MIN_CHANGE"])
    if any(value.denominator != 1 for value in minimal_model):
        raise AssertionError("PARI's global minimal model is not integral")
    if change_weierstrass_model(source_model, change) != minimal_model:
        raise AssertionError("PARI's minimal change failed exact model replay")
    source_discriminant = weierstrass_invariants(source_model)["discriminant"]
    target_discriminant = weierstrass_invariants(minimal_model)["discriminant"]
    if source_discriminant != change.u**12 * target_discriminant:
        raise AssertionError("the minimal change failed discriminant scaling")
    metadata = {
        "engine": "PARI/GP ellminimalmodel(E,&v)",
        "pari_version": pari_version(),
        "timeout_seconds": timeout,
        "stack_bytes": stack_bytes,
        "program_sha256": sha256(program.encode()).hexdigest(),
        "exact_model_change_replay": True,
        "exact_discriminant_scaling": True,
    }
    return minimal_model, change, metadata


def short_certificate_model(
    minimal_model: Sequence[Fraction | int],
) -> tuple[WeierstrassModel, WeierstrassChange]:
    """Return the integral short model ``[-27*c4,-54*c6]`` of a minimal model."""

    a1, _a2, a3, _a4, _a6 = (Q(value) for value in minimal_model)
    invariants = weierstrass_invariants(minimal_model)
    b2 = invariants["b2"]
    change = WeierstrassChange(
        Q(1, 6),
        -b2 / 12,
        -a1 / 2,
        a1 * b2 / 24 - a3 / 2,
    )
    expected = (
        Q(0),
        Q(0),
        Q(0),
        -27 * invariants["c4"],
        -54 * invariants["c6"],
    )
    if change_weierstrass_model(minimal_model, change) != expected:
        raise AssertionError("the canonical short-model change failed exact replay")
    if any(value.denominator != 1 for value in expected):
        raise AssertionError("the canonical certificate model is not integral")
    return expected, change


def _model_record(model: Iterable[Fraction | int]) -> list[str]:
    return [fraction_text(value) for value in model]


def _points_record(points: Iterable[RationalPoint]) -> list[list[str]]:
    return [[fraction_text(x), fraction_text(y)] for x, y in points]


def build_specialization_record(
    data: Q12O5867Data,
    a: int,
    b: int,
    *,
    relation_primes: Sequence[int] = (2, 3, 5),
    reduction_prime_bound: int = 500,
    gp_timeout: float = 300.0,
    gp_stack_bytes: int = 512_000_000,
) -> dict[str, Any]:
    """Build the strongest exact specialization record available."""

    specialization = evaluate_projective_specialization(data, a, b)
    parameter = {
        "normalized_projective": [specialization.a, specialization.b],
        "chart": "infinity" if specialization.b == 0 else "finite",
        "affine_value": (
            None
            if specialization.b == 0
            else fraction_text(Q(specialization.a, specialization.b))
        ),
        "homogeneous_weights": {"A": 8, "B": 12, "x": 4, "y": 6},
    }
    record: dict[str, Any] = {
        "schema": "elliptic-curves.q12o5867-specialization.v1",
        "family": "q12o5867-rootless-rank17",
        "parameter": parameter,
        "inputs": {
            "model": {"path": str(data.model_path), "sha256": data.model_sha256},
            "sections": {
                "path": str(data.sections_path),
                "sha256": data.sections_sha256,
            },
        },
        "projective_specialization": {
            "model": _model_record(specialization.model),
            "nonsingular": True,
            "section_count": len(specialization.points),
            "all_sections_on_curve": True,
        },
        "claim_boundary": [
            "The finite-quotient result is a rank lower bound only.",
            "No saturation, Selmer upper bound, exact rank, or conductor is claimed.",
            "PARI ellminimalmodel supplies the computational global-minimal-model step.",
        ],
    }
    scale, scale_record = integral_scaling_factor(
        specialization.coefficient_a,
        specialization.coefficient_b,
        timeout=min(gp_timeout, 120.0),
        stack_bytes=min(gp_stack_bytes, 128_000_000),
    )
    integral_model, integral_points = integralize_specialization(specialization, scale)
    record["integral_short_specialization"] = {
        "scaling": scale_record,
        "model": _model_record(integral_model),
        "points": _points_record(integral_points),
        "section_count": 17,
        "all_sections_on_curve": True,
    }
    try:
        minimal_model, minimal_change, minimal_metadata = global_minimal_model_with_change(
            integral_model, timeout=gp_timeout, stack_bytes=gp_stack_bytes
        )
    except (FileNotFoundError, RuntimeError, AssertionError, ValueError) as error:
        record["status"] = "BLOCKED_AFTER_EXACT_INTEGRAL_Q12O5867_SPECIALIZATION"
        record["blocker"] = {
            "stage": "global-minimal-model-with-change",
            "exception": type(error).__name__,
            "detail": str(error),
        }
        return record

    minimal_points = tuple(
        source_point_to_target(point, minimal_change) for point in integral_points
    )
    if any(not is_on_weierstrass_curve(minimal_model, point) for point in minimal_points):
        raise AssertionError("a transported section missed the global minimal model")
    record["global_minimal_specialization"] = {
        "model": _model_record(minimal_model),
        "integral_short_to_minimal_change_u_r_s_t": minimal_change.to_record(),
        "points": _points_record(minimal_points),
        "section_count": 17,
        "all_sections_on_curve": True,
        **minimal_metadata,
    }

    try:
        certificate_model, minimal_to_short = short_certificate_model(minimal_model)
        certificate_points = tuple(
            source_point_to_target(point, minimal_to_short) for point in minimal_points
        )
        if any(
            not is_on_weierstrass_curve(certificate_model, point)
            for point in certificate_points
        ):
            raise AssertionError(
                "a transported section missed the certificate short model"
            )
        certificate_attempts = []
        successful_certificate = None
        for relation_prime in relation_primes:
            certificate = build_finite_quotient_certificate(
                certificate_model,
                certificate_points,
                relation_prime=int(relation_prime),
                prime_bound=reduction_prime_bound,
            )
            verify_finite_quotient_certificate(
                certificate_model, certificate_points, certificate
            )
            certificate_attempts.append(certificate)
            if certificate["certified_independent"]:
                successful_certificate = certificate
                break
    except (RuntimeError, AssertionError, ValueError) as error:
        record["status"] = "BLOCKED_AFTER_EXACT_GLOBAL_MINIMAL_Q12O5867_SPECIALIZATION"
        record["blocker"] = {
            "stage": "finite-quotient-independence-certificate",
            "exception": type(error).__name__,
            "detail": str(error),
        }
        return record
    record["finite_quotient_independence"] = {
        "certificate_short_model": _model_record(certificate_model),
        "minimal_to_certificate_short_change_u_r_s_t": minimal_to_short.to_record(),
        "points": _points_record(certificate_points),
        "all_sections_on_curve": True,
        "attempts": certificate_attempts,
        "successful_relation_prime": (
            None
            if successful_certificate is None
            else successful_certificate["relation_prime"]
        ),
        "certified_independent": successful_certificate is not None,
        "certified_rank_lower_bound": 17 if successful_certificate is not None else None,
    }
    if successful_certificate is None:
        record["status"] = "PASS_EXACT_SPECIALIZATION_WITH_BOUNDED_CERTIFICATE_FAILURE"
    else:
        record["status"] = "PASS_EXACT_Q12O5867_SPECIALIZED_GENERIC_RANK17_LOWER_BOUND"
    return record
