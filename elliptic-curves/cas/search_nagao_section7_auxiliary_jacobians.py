#!/usr/bin/env python3
"""Constructive section-7 search through auxiliary elliptic slice Jacobians.

At Nagao's certified rank-20 specialization ``T0=5081/47``, two accidental
quartic points define the four declared slices

``x = T + n_plus`` and ``x = -T + n_minus``.

Substitution into the section-7 quartic leaves a genus-one quartic in ``T``.
The known point above ``T0`` makes that quartic an elliptic curve.  This script
constructs the exact pointed-quartic birational map, heuristically selects a
rank-seven auxiliary subgroup from the intersections with the known generic
sections, and transports bounded-coefficient group-law combinations back to
rational section-7 parameters.

This is not another ``hyperellratpoints`` box: every emitted parameter has
projective height strictly greater than 200,000.  All intersections with the
12 visible and 9 companion section abscissae, plus the sibling slice pilot's
11 subtarget parameters, are explicitly excluded.  Exact discriminant-radical
proxies precede a strictly capped exact-conductor population.  A numerical
auxiliary height rank is triage, not a rank theorem for either curve.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import itertools
import json
from math import comb, isqrt
from pathlib import Path
import platform
import re
import shlex
import subprocess
import sys
from typing import Any, Iterable, Sequence

from ek_k3 import rational_to_string
from nagao_1994 import quartic_value, short_jacobian_coefficients
from nagao_1994_section7 import (
    SECTION7_CONSTRUCTION,
    SECTION7_CONSTRUCTOR_PARAMETER,
    SECTION7_LINEAR_COMPANION_SECTIONS,
    SECTION7_QUADRATIC_COMPANION_SECTIONS,
    SECTION7_ROOTS,
    section7_primitive_quartic_coefficients,
)
from pari_bridge import minimal_curve_data, pari_version
from search_extra_points import gp_rational, gp_vector, run_gp
from search_nagao_rank20_t5081_neighborhood import (
    conductor_radical_proxy,
    homogenized_discriminant,
)
from triage_nagao_rank13_finalists import (
    height_matrix_replay,
    point_on_short_curve,
    stable_height_rank,
)


Q = Fraction
T0 = SECTION7_CONSTRUCTOR_PARAMETER
NAIVE_HEIGHT_EXCLUSION = 200_000
MAX_ABSOLUTE_COEFFICIENT = 3
MAX_L1_NORM = 5
EXACT_CONDUCTOR_KEEP = 64
PROXY_TRIAL_BOUND = 2_000
TARGET_LOG_CONDUCTOR = Decimal("182.72")
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_section7_auxiliary_jacobians.py "
    "--skip-ellrank"
)

# These are exact outputs of the independent bounded slice pilot.  They are
# exclusions here, never training examples or discoveries of this lane.
SIBLING_BOUNDED_SUBTARGET_PARAMETERS = tuple(
    Q(value)
    for value in (
        "-1529/17",
        "-3795/94",
        "6115/6",
        "5605/94",
        "2920/47",
        "7649/47",
        "4941/47",
        "241/47",
        "13108/91",
        "-3190/47",
        "-3754/47",
    )
)


@dataclass(frozen=True)
class SliceSpecification:
    label: str
    slope: Fraction
    intercept: Fraction
    base_x: Fraction
    base_y: Fraction
    parent_accidental_label: str


PRIMARY_X0 = Q(11461, 47)
PRIMARY_Y0 = Q(214287355335, 103823)
SECONDARY_X0 = Q(-145339, 4277)
SECONDARY_Y0 = Q(731165923457085, 18292729)
SLICE_SPECIFICATIONS = (
    SliceSpecification(
        "primary-plus",
        Q(1),
        Q(6380, 47),
        PRIMARY_X0,
        PRIMARY_Y0,
        "T0 accidental x=11461/47",
    ),
    SliceSpecification(
        "primary-minus",
        Q(-1),
        PRIMARY_X0 + T0,
        PRIMARY_X0,
        PRIMARY_Y0,
        "T0 accidental x=11461/47",
    ),
    SliceSpecification(
        "secondary-plus",
        Q(1),
        Q(-12930, 91),
        SECONDARY_X0,
        SECONDARY_Y0,
        "T0 accidental x=-145339/4277",
    ),
    SliceSpecification(
        "secondary-minus",
        Q(-1),
        SECONDARY_X0 + T0,
        SECONDARY_X0,
        SECONDARY_Y0,
        "T0 accidental x=-145339/4277",
    ),
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def polynomial_add(*values: Sequence[Fraction]) -> tuple[Fraction, ...]:
    length = max(len(value) for value in values)
    answer = [Q(0)] * length
    for value in values:
        for index, coefficient in enumerate(value):
            answer[index] += Q(coefficient)
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def polynomial_multiply(
    left: Sequence[Fraction], right: Sequence[Fraction]
) -> tuple[Fraction, ...]:
    answer = [Q(0)] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            answer[left_index + right_index] += Q(left_value) * Q(right_value)
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def polynomial_power(
    value: Sequence[Fraction], exponent: int
) -> tuple[Fraction, ...]:
    if exponent < 0:
        raise ValueError("a polynomial exponent must be nonnegative")
    answer = (Q(1),)
    for _ in range(exponent):
        answer = polynomial_multiply(answer, value)
    return answer


def polynomial_value(
    coefficients: Sequence[Fraction], value: Fraction
) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * Q(value) + Q(coefficient)
    return answer


# Coefficients in ascending T-order for the ascending X-coefficients e,d,c,b,a.
SECTION7_QUARTIC_COEFFICIENT_POLYNOMIALS = (
    (Q(557726319412900), Q(0), Q(23718659440), Q(0), Q(-910748), Q(0), Q(9)),
    tuple(
        Q(18) * value
        for value in (Q(-884640359570), Q(0), Q(-17901331), Q(0), Q(354))
    ),
    tuple(
        Q(-3) * value
        for value in (Q(-52052853547), Q(0), Q(-668642), Q(0), Q(6))
    ),
    tuple(Q(-54) * value for value in (Q(11538729), Q(0), Q(118))),
    tuple(Q(9) * value for value in (Q(96714), Q(0), Q(1))),
)


def slice_quartic_polynomial(specification: SliceSpecification) -> tuple[Fraction, ...]:
    linear = specification.intercept, specification.slope
    answer = polynomial_add(
        *(
            polynomial_multiply(
                coefficient_polynomial, polynomial_power(linear, x_power)
            )
            for x_power, coefficient_polynomial in enumerate(
                SECTION7_QUARTIC_COEFFICIENT_POLYNOMIALS
            )
        )
    )
    if len(answer) != 5:
        raise AssertionError(
            f"slice {specification.label} has degree {len(answer) - 1}, not four"
        )
    for parameter in (Q(1), Q(2), T0):
        expected = quartic_value(
            section7_primitive_quartic_coefficients(parameter),
            specification.slope * parameter + specification.intercept,
        )
        if polynomial_value(answer, parameter) != expected:
            raise AssertionError("the expanded slice polynomial failed exactly")
    return answer


def translate_polynomial(
    coefficients: Sequence[Fraction], translation: Fraction
) -> tuple[Fraction, ...]:
    """Return ascending coefficients of f(translation + u)."""

    answer = [Q(0)] * len(coefficients)
    for power, coefficient in enumerate(coefficients):
        for shifted_power in range(power + 1):
            answer[shifted_power] += (
                Q(coefficient)
                * comb(power, shifted_power)
                * Q(translation) ** (power - shifted_power)
            )
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


@dataclass(frozen=True)
class AuxiliarySlice:
    specification: SliceSpecification
    quartic_coefficients: tuple[Fraction, ...]
    shifted_coefficients: tuple[Fraction, ...]
    weierstrass_coefficients: tuple[Fraction, ...]

    @property
    def base_y(self) -> Fraction:
        return self.specification.base_y

    def quartic_value(self, parameter: Fraction) -> Fraction:
        return polynomial_value(self.quartic_coefficients, Q(parameter))

    def forward(
        self, point: tuple[Fraction, Fraction]
    ) -> tuple[Fraction, Fraction]:
        """Map the pointed quartic to its generalized Weierstrass model."""

        parameter, ordinate = (Q(value) for value in point)
        if ordinate**2 != self.quartic_value(parameter):
            raise ValueError("the point missed the auxiliary quartic")
        u_value = parameter - T0
        if u_value == 0:
            raise ValueError("the chosen base point maps to infinity")
        _, d_value, c_value, _, _ = self.shifted_coefficients
        q_value = self.base_y
        x_value = (
            2 * q_value * (ordinate + q_value) + d_value * u_value
        ) / u_value**2
        y_value = (
            4 * q_value**2 * (ordinate + q_value)
            + 2
            * q_value
            * (d_value * u_value + c_value * u_value**2)
            - d_value**2 * u_value**2 / (2 * q_value)
        ) / u_value**3
        image = x_value, y_value
        if not point_on_short_curve(self.weierstrass_coefficients, image):
            raise AssertionError("the pointed-quartic map missed its auxiliary curve")
        return image

    def inverse(
        self, point: tuple[Fraction, Fraction] | None
    ) -> tuple[Fraction, Fraction] | None:
        """Invert a nonexceptional auxiliary point back to (T,y)."""

        if point is None:
            return T0, self.base_y
        if not point_on_short_curve(self.weierstrass_coefficients, point):
            raise ValueError("the point missed the auxiliary Weierstrass model")
        x_value, y_value = point
        if y_value == 0:
            return None
        _, d_value, c_value, _, _ = self.shifted_coefficients
        q_value = self.base_y
        u_value = (
            4 * q_value**2 * (x_value + c_value) - d_value**2
        ) / (2 * q_value * y_value)
        if u_value == 0:
            return None
        ordinate = (
            (x_value * u_value**2 - d_value * u_value) / (2 * q_value)
            - q_value
        )
        answer = T0 + u_value, ordinate
        if ordinate**2 != self.quartic_value(answer[0]):
            raise AssertionError("the inverse auxiliary map failed exactly")
        return answer


def make_auxiliary_slice(specification: SliceSpecification) -> AuxiliarySlice:
    quartic = slice_quartic_polynomial(specification)
    shifted = translate_polynomial(quartic, T0)
    if len(shifted) != 5 or shifted[0] != specification.base_y**2:
        raise AssertionError("the declared accidental base point missed its slice")
    _, d_value, c_value, b_value, a_value = shifted
    q_value = specification.base_y
    a1 = d_value / q_value
    a2 = c_value - d_value**2 / (4 * q_value**2)
    a3 = 2 * q_value * b_value
    a4 = -4 * q_value**2 * a_value
    a6 = a_value * (d_value**2 - 4 * q_value**2 * c_value)
    return AuxiliarySlice(specification, quartic, shifted, (a1, a2, a3, a4, a6))


def rational_square_root(value: Fraction) -> Fraction | None:
    value = Q(value)
    if value < 0:
        return None
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator**2 != value.numerator or denominator**2 != value.denominator:
        return None
    return Q(numerator, denominator)


def generic_intersection_parameters(
    auxiliary: AuxiliarySlice,
) -> tuple[tuple[Fraction, tuple[str, ...]], ...]:
    """Solve all 12 visible and 9 companion abscissa intersections exactly."""

    specification = auxiliary.specification
    found: dict[Fraction, set[str]] = {}
    for root in SECTION7_ROOTS:
        for slope in (Q(1), Q(-1)):
            if slope == specification.slope:
                if Q(root) == specification.intercept:
                    raise AssertionError("a slice unexpectedly equals a visible section")
                continue
            parameter = Q(
                Q(root) - specification.intercept,
                specification.slope - slope,
            )
            found.setdefault(parameter, set()).add(
                f"visible-x={root}{'+' if slope == 1 else '-'}T"
            )
    for section in SECTION7_LINEAR_COMPANION_SECTIONS:
        if section.slope == specification.slope:
            if section.intercept == specification.intercept:
                raise AssertionError("a slice unexpectedly equals a linear companion")
            continue
        parameter = Q(
            section.intercept - specification.intercept,
            specification.slope - section.slope,
        )
        found.setdefault(parameter, set()).add(section.label)
    for section in SECTION7_QUADRATIC_COMPANION_SECTIONS:
        coefficient_a = section.quadratic_coefficient
        coefficient_b = -specification.slope
        coefficient_c = section.constant_coefficient - specification.intercept
        square_root = rational_square_root(
            coefficient_b**2 - 4 * coefficient_a * coefficient_c
        )
        if square_root is None:
            continue
        for sign in (1, -1):
            parameter = Q(
                -coefficient_b + sign * square_root,
                2 * coefficient_a,
            )
            found.setdefault(parameter, set()).add(section.label)

    for parameter in found:
        if rational_square_root(auxiliary.quartic_value(parameter)) is None:
            raise AssertionError("a declared generic intersection was not rational")
    return tuple(
        (parameter, tuple(sorted(labels)))
        for parameter, labels in sorted(
            found.items(),
            key=lambda item: (
                max(abs(item[0].numerator), item[0].denominator),
                item[0],
            ),
        )
    )


WeierstrassPoint = tuple[Fraction, Fraction] | None


def weierstrass_negate(
    coefficients: Sequence[Fraction], point: WeierstrassPoint
) -> WeierstrassPoint:
    if point is None:
        return None
    a1, _, a3, _, _ = (Q(value) for value in coefficients)
    return point[0], -point[1] - a1 * point[0] - a3


def weierstrass_add(
    coefficients: Sequence[Fraction],
    left: WeierstrassPoint,
    right: WeierstrassPoint,
) -> WeierstrassPoint:
    if left is None:
        return right
    if right is None:
        return left
    a1, a2, a3, a4, a6 = (Q(value) for value in coefficients)
    x1, y1 = left
    x2, y2 = right
    if x1 == x2 and y1 + y2 + a1 * x1 + a3 == 0:
        return None
    if x1 == x2:
        denominator = 2 * y1 + a1 * x1 + a3
        if denominator == 0:
            return None
        slope = (3 * x1**2 + 2 * a2 * x1 + a4 - a1 * y1) / denominator
        intercept = (
            -x1**3 + a4 * x1 + 2 * a6 - a3 * y1
        ) / denominator
    else:
        slope = (y2 - y1) / (x2 - x1)
        intercept = (y1 * x2 - y2 * x1) / (x2 - x1)
    x3 = slope**2 + a1 * slope - a2 - x1 - x2
    y3 = -(slope + a1) * x3 - intercept - a3
    answer = x3, y3
    if not point_on_short_curve(coefficients, answer):
        raise AssertionError("the exact generalized group law failed")
    return answer


def weierstrass_multiply(
    coefficients: Sequence[Fraction], point: WeierstrassPoint, scalar: int
) -> WeierstrassPoint:
    if scalar < 0:
        return weierstrass_multiply(
            coefficients, weierstrass_negate(coefficients, point), -scalar
        )
    answer: WeierstrassPoint = None
    addend = point
    while scalar:
        if scalar & 1:
            answer = weierstrass_add(coefficients, answer, addend)
        addend = weierstrass_add(coefficients, addend, addend)
        scalar >>= 1
    return answer


def exact_seed_points(
    auxiliary: AuxiliarySlice,
) -> tuple[
    tuple[tuple[Fraction, Fraction], ...],
    tuple[dict[str, Any], ...],
]:
    """Map both signs above every rational generic intersection."""

    points: list[tuple[Fraction, Fraction]] = []
    records = []
    for parameter, labels in generic_intersection_parameters(auxiliary):
        square_root = rational_square_root(auxiliary.quartic_value(parameter))
        if square_root is None:
            raise AssertionError("the generic ordinate disappeared")
        images = []
        if parameter != T0:
            for ordinate in (square_root, -square_root):
                image = auxiliary.forward((parameter, ordinate))
                if image not in points:
                    points.append(image)
                images.append(image)
        records.append(
            {
                "constructor_parameter": rational_to_string(parameter),
                "section_labels": list(labels),
                "ordinate_signs_mapped": len(images),
            }
        )
    return tuple(points), tuple(records)


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
            f"E=ellinit([{curve}]);",
            f"P=[{point_vector}];",
            "R=ellrank(E,0,P);",
            'print("RANK_BOUNDS ",R[1],"|",R[2],"|",#R[4]);',
            "quit",
        )
    ) + "\n"
    try:
        output, wall_seconds = run_gp(
            program, timeout=timeout, stack_bytes=stack_bytes
        )
    except subprocess.TimeoutExpired:
        return {
            "status": "timeout",
            "timeout_seconds": timeout,
            "one_attempt_no_retry": True,
        }
    except RuntimeError as error:
        return {
            "status": "error",
            "timeout_seconds": timeout,
            "one_attempt_no_retry": True,
            "error": str(error)[:500],
        }
    match = re.search(r"RANK_BOUNDS (\d+)\|(\d+)\|(\d+)", output)
    if match is None:
        raise AssertionError("PARI ellrank omitted its bounded result")
    return {
        "status": "completed",
        "lower_bound": int(match.group(1)),
        "upper_bound": int(match.group(2)),
        "returned_generator_count": int(match.group(3)),
        "effort": 0,
        "wall_seconds": wall_seconds,
        "scope": "PARI ellrank triage on the auxiliary curve",
    }


def coefficient_vectors(
    dimension: int,
    *,
    maximum_absolute_coefficient: int = MAX_ABSOLUTE_COEFFICIENT,
    maximum_l1_norm: int = MAX_L1_NORM,
) -> tuple[tuple[int, ...], ...]:
    if dimension < 1 or maximum_absolute_coefficient < 1 or maximum_l1_norm < 1:
        raise ValueError("invalid coefficient-vector bounds")
    return tuple(
        vector
        for vector in itertools.product(
            range(-maximum_absolute_coefficient, maximum_absolute_coefficient + 1),
            repeat=dimension,
        )
        if 0 < sum(abs(value) for value in vector) <= maximum_l1_norm
    )


def linear_combination(
    coefficients: Sequence[Fraction],
    basis: Sequence[tuple[Fraction, Fraction]],
    vector: Sequence[int],
) -> WeierstrassPoint:
    if len(basis) != len(vector):
        raise ValueError("the coefficient vector has the wrong dimension")
    answer: WeierstrassPoint = None
    for point, scalar in zip(basis, vector):
        if scalar:
            answer = weierstrass_add(
                coefficients,
                answer,
                weierstrass_multiply(coefficients, point, scalar),
            )
    return answer


def canonical_parameter(parameter: Fraction) -> Fraction:
    """The section-7 Jacobian is even in T; retain its nonnegative representative."""

    return abs(Q(parameter))


def projective_height(parameter: Fraction) -> int:
    parameter = Q(parameter)
    return max(abs(parameter.numerator), parameter.denominator)


def rational_record(value: Fraction) -> str:
    return rational_to_string(Q(value))


def signature_records(signatures: Iterable[tuple[str, tuple[int, ...]]]) -> list[dict[str, Any]]:
    return [
        {"slice": label, "coefficient_vector": list(vector)}
        for label, vector in signatures
    ]


def generate_candidates(
    slices: Sequence[AuxiliarySlice],
    selected_bases: dict[str, tuple[tuple[Fraction, Fraction], ...]],
    vectors: Sequence[tuple[int, ...]],
) -> tuple[dict[Fraction, set[tuple[str, tuple[int, ...]]]], dict[str, Any]]:
    generic_parameters = {
        canonical_parameter(parameter)
        for auxiliary in slices
        for parameter, _ in generic_intersection_parameters(auxiliary)
    }
    sibling_parameters = {
        canonical_parameter(parameter)
        for parameter in SIBLING_BOUNDED_SUBTARGET_PARAMETERS
    }
    candidates: dict[Fraction, set[tuple[str, tuple[int, ...]]]] = {}
    counts = {
        "group_elements_attempted": 0,
        "exceptional_inverse_images": 0,
        "zero_parameters": 0,
        "naive_height_at_most_200000": 0,
        "generic_intersections": 0,
        "sibling_bounded_subtarget_parameters": 0,
        "singular_section7_parameters": 0,
    }
    for auxiliary in slices:
        basis = selected_bases[auxiliary.specification.label]
        for vector in vectors:
            counts["group_elements_attempted"] += 1
            point = linear_combination(
                auxiliary.weierstrass_coefficients, basis, vector
            )
            inverse = auxiliary.inverse(point)
            if inverse is None:
                counts["exceptional_inverse_images"] += 1
                continue
            signed_parameter, ordinate = inverse
            if ordinate**2 != auxiliary.quartic_value(signed_parameter):
                raise AssertionError("a generated inverse missed the auxiliary quartic")
            parameter = canonical_parameter(signed_parameter)
            if parameter == 0:
                counts["zero_parameters"] += 1
                continue
            if projective_height(parameter) <= NAIVE_HEIGHT_EXCLUSION:
                counts["naive_height_at_most_200000"] += 1
                continue
            if parameter in generic_parameters:
                counts["generic_intersections"] += 1
                continue
            if parameter in sibling_parameters:
                counts["sibling_bounded_subtarget_parameters"] += 1
                continue
            try:
                homogenized_discriminant(parameter)
            except ValueError:
                counts["singular_section7_parameters"] += 1
                continue
            candidates.setdefault(parameter, set()).add(
                (auxiliary.specification.label, tuple(vector))
            )
    counts["exact_unique_parameters"] = len(candidates)
    counts["generic_parameter_exclusion_count"] = len(generic_parameters)
    counts["sibling_parameter_exclusion_count"] = len(sibling_parameters)
    return candidates, counts


def candidate_stream_digest(parameters: Iterable[Fraction]) -> str:
    text = "\n".join(
        rational_to_string(parameter)
        for parameter in sorted(
            parameters,
            key=lambda value: (
                projective_height(value),
                value.numerator,
                value.denominator,
            ),
        )
    )
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def proxy_population(
    candidates: dict[Fraction, set[tuple[str, tuple[int, ...]]]],
) -> tuple[dict[str, Any], ...]:
    records = []
    for parameter, sources in candidates.items():
        proxy = conductor_radical_proxy(
            parameter, trial_prime_bound=PROXY_TRIAL_BOUND
        )
        records.append(
            {
                "constructor_parameter": rational_to_string(parameter),
                "projective_height": projective_height(parameter),
                "projective_height_bits": projective_height(parameter).bit_length(),
                "radical_proxy": proxy,
                "source_count": len(sources),
                "first_sources": signature_records(sorted(sources)[:4]),
            }
        )
    return tuple(
        sorted(
            records,
            key=lambda record: (
                record["radical_proxy"]["log_radical_upper_proxy"],
                record["projective_height"],
                Q(record["constructor_parameter"]),
            ),
        )
    )


def exact_conductor_population(
    proxy_records: Sequence[dict[str, Any]],
    *,
    keep_count: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[dict[str, Any], ...], tuple[dict[str, Any], ...]]:
    completed = []
    errors = []
    for index, record in enumerate(proxy_records[:keep_count], start=1):
        parameter = Q(record["constructor_parameter"])
        coefficients = short_jacobian_coefficients(
            SECTION7_CONSTRUCTION, parameter
        )
        try:
            conductor = minimal_curve_data(
                coefficients,
                timeout=timeout,
                stack_bytes=stack_bytes,
            )
        except (subprocess.TimeoutExpired, RuntimeError) as error:
            errors.append(
                {
                    "constructor_parameter": rational_to_string(parameter),
                    "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
                    "one_attempt_no_retry": True,
                    "error": str(error)[:500],
                }
            )
            print(
                f"conductor {index}/{min(keep_count, len(proxy_records))} "
                f"T={parameter} status={errors[-1]['status']}",
                flush=True,
            )
            continue
        exact = {
            **record,
            "conductor": str(conductor["conductor"]),
            "log_conductor": conductor["log_conductor"],
            "root_number": conductor["root_number"],
            "minimal_discriminant": str(conductor["minimal_discriminant"]),
            "minimal_model": [str(value) for value in conductor["minimal_model"]],
            "below_strict_log_conductor_target": (
                Decimal(conductor["log_conductor"]) < TARGET_LOG_CONDUCTOR
            ),
            "status": "completed",
        }
        completed.append(exact)
        print(
            f"conductor {index}/{min(keep_count, len(proxy_records))} "
            f"T={parameter} lnN={conductor['log_conductor']} "
            f"subtarget={exact['below_strict_log_conductor_target']}",
            flush=True,
        )
    return tuple(completed), tuple(errors)


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--ellrank-timeout", type=float, default=20.0)
    parser.add_argument(
        "--skip-ellrank",
        action="store_true",
        help="use the two-precision height replay alone for auxiliary-rank triage",
    )
    parser.add_argument("--conductor-timeout", type=float, default=15.0)
    parser.add_argument("--stack-bytes", type=int, default=256_000_000)
    parser.add_argument("--conductor-keep", type=int, default=EXACT_CONDUCTOR_KEEP)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts"
            / "generated-results"
            / "elliptic_nagao_section7_auxiliary_jacobians.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 0 < args.height_timeout <= 60:
        raise SystemExit("--height-timeout must be in (0,60]")
    if not 0 < args.ellrank_timeout <= 60:
        raise SystemExit("--ellrank-timeout must be in (0,60]")
    if not 0 < args.conductor_timeout <= 60:
        raise SystemExit("--conductor-timeout must be in (0,60]")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes is too small")
    if not 1 <= args.conductor_keep <= 256:
        raise SystemExit("--conductor-keep must be in [1,256]")

    slices = tuple(make_auxiliary_slice(specification) for specification in SLICE_SPECIFICATIONS)
    slice_records = []
    selected_bases: dict[str, tuple[tuple[Fraction, Fraction], ...]] = {}
    for auxiliary in slices:
        seed_points, intersections = exact_seed_points(auxiliary)
        height_runs = height_matrix_replay(
            auxiliary.weierstrass_coefficients,
            seed_points,
            precisions=(72, 120),
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
        )
        auxiliary_rank = stable_height_rank(height_runs)
        selected_indices = tuple(height_runs[-1]["subset_indices_one_based"])
        selected = tuple(seed_points[index - 1] for index in selected_indices)
        if auxiliary_rank != len(selected):
            raise AssertionError("the auxiliary numerical subset has the wrong size")
        selected_bases[auxiliary.specification.label] = selected
        ellrank = (
            {
                "status": "not_run",
                "reason": (
                    "the declared two-precision height replay is sufficient for "
                    "heuristic auxiliary-rank triage"
                ),
            }
            if args.skip_ellrank
            else ellrank_probe(
                auxiliary.weierstrass_coefficients,
                selected,
                timeout=args.ellrank_timeout,
                stack_bytes=args.stack_bytes,
            )
        )
        slice_records.append(
            {
                "label": auxiliary.specification.label,
                "slice_abscissa": (
                    f"{rational_to_string(auxiliary.specification.slope)}*T+"
                    f"{rational_to_string(auxiliary.specification.intercept)}"
                ),
                "slope": rational_to_string(auxiliary.specification.slope),
                "intercept": rational_to_string(auxiliary.specification.intercept),
                "base_point": {
                    "T": rational_to_string(T0),
                    "x": rational_to_string(auxiliary.specification.base_x),
                    "y": rational_to_string(auxiliary.specification.base_y),
                    "source": auxiliary.specification.parent_accidental_label,
                    "maps_to_auxiliary_infinity": True,
                },
                "quartic_coefficients_ascending_in_T": [
                    rational_to_string(value)
                    for value in auxiliary.quartic_coefficients
                ],
                "shifted_coefficients_ascending_in_u_T_minus_T0": [
                    rational_to_string(value)
                    for value in auxiliary.shifted_coefficients
                ],
                "auxiliary_weierstrass_coefficients_a1_a2_a3_a4_a6": [
                    rational_to_string(value)
                    for value in auxiliary.weierstrass_coefficients
                ],
                "generic_intersections": list(intersections),
                "generic_intersection_parameter_count": len(intersections),
                "mapped_signed_seed_count": len(seed_points),
                "height_replay": list(height_runs),
                "stable_numerical_auxiliary_rank": auxiliary_rank,
                "selected_seed_indices_one_based": list(selected_indices),
                "pari_ellrank_effort_zero": ellrank,
                "rank_status": (
                    "heuristic_numerical_lower_direction_count; not an exact "
                    "auxiliary rank certificate"
                ),
            }
        )
        print(
            f"slice {auxiliary.specification.label} generic={len(intersections)} "
            f"seed={len(seed_points)} stable_aux_rank={auxiliary_rank} "
            f"ellrank={ellrank['status']}",
            flush=True,
        )

    dimension = len(next(iter(selected_bases.values())))
    vectors = coefficient_vectors(dimension)
    candidates, generation_counts = generate_candidates(
        slices, selected_bases, vectors
    )
    proxy_records = proxy_population(candidates)
    checkpoint = {
        "vector_count_per_slice": len(vectors),
        "slice_count": len(slices),
        **generation_counts,
        "candidate_stream_sha256": candidate_stream_digest(candidates),
        "proxy_population_count": len(proxy_records),
        "minimum_log_radical_upper_proxy": (
            proxy_records[0]["radical_proxy"]["log_radical_upper_proxy"]
            if proxy_records
            else None
        ),
    }
    print("GENERATION_CHECKPOINT " + json.dumps(checkpoint, sort_keys=True), flush=True)

    conductor_records, conductor_errors = exact_conductor_population(
        proxy_records,
        keep_count=min(args.conductor_keep, len(proxy_records)),
        timeout=args.conductor_timeout,
        stack_bytes=args.stack_bytes,
    )
    subtarget = tuple(
        record
        for record in conductor_records
        if record["below_strict_log_conductor_target"]
    )
    artifact = {
        "schema_version": 1,
        "status": "bounded auxiliary-Jacobian constructive search complete",
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
        },
        "method": {
            "constructor_parameter_T0": rational_to_string(T0),
            "slice_count": len(slices),
            "no_naive_hyperellratpoints_calls": True,
            "naive_projective_height_exclusion": NAIVE_HEIGHT_EXCLUSION,
            "group_coefficient_absolute_bound": MAX_ABSOLUTE_COEFFICIENT,
            "group_coefficient_l1_bound": MAX_L1_NORM,
            "pointed_quartic_map": (
                "exact birational transformation to generalized Weierstrass "
                "form with the accidental T0 point as infinity"
            ),
            "canonicalization": (
                "T and -T have the same section-7 Jacobian coefficients; exact "
                "conductor work keeps |T|"
            ),
        },
        "sibling_bounded_exclusions": {
            "source_role": "independent pilot outputs excluded before selection",
            "constructor_parameters": [
                rational_to_string(value)
                for value in SIBLING_BOUNDED_SUBTARGET_PARAMETERS
            ],
        },
        "slices": slice_records,
        "generation": checkpoint,
        "proxy_selection": {
            "trial_prime_bound": PROXY_TRIAL_BOUND,
            "exact_proxy_population_count": len(proxy_records),
            "exact_conductor_keep": min(args.conductor_keep, len(proxy_records)),
            "top_proxy_records": list(proxy_records[: min(128, len(proxy_records))]),
        },
        "exact_conductors": {
            "attempted": min(args.conductor_keep, len(proxy_records)),
            "completed": len(conductor_records),
            "errors_or_timeouts": list(conductor_errors),
            "records": list(conductor_records),
            "subtarget_count": len(subtarget),
            "subtarget_records": list(subtarget),
        },
        "outcome": {
            "subtarget_parameters_for_rank_triage": [
                record["constructor_parameter"] for record in subtarget
            ],
            "rank21_certified_in_this_lane": False,
            "interpretation": (
                "the auxiliary group law constructs exact section-7 points; "
                "a subtarget conductor is only a candidate for later rank triage"
            ),
        },
        "bounded_scope": {
            "every_subprocess_timeout_at_most_60_seconds": True,
            "one_attempt_no_retry": True,
            "no_detached_processes": True,
            "negative_result_is_not_a_rank_upper_bound": True,
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "actual_command": " ".join(
            shlex.quote(part) for part in [sys.executable, *sys.argv]
        ),
        "script_sha256": file_sha256(Path(__file__).resolve()),
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "platform": platform.platform(),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"wrote {args.output}: candidates={len(candidates)} "
        f"conductors={len(conductor_records)} subtarget={len(subtarget)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
