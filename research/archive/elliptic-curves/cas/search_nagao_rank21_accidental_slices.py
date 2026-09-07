#!/usr/bin/env python3
"""Construct specialization slices through accidental points at ``T=5081/47``.

The section-7 rank-20 fiber contains 34 exact affine quartic abscissas in the
pinned certificate search.  This script reconstructs that set independently,
removes the twelve visible, six linear-companion, and three
quadratic-companion generic sections, and uses each remaining point
``(T0,x0,y0)`` to define the exact slices

``x = m*T + (x0-m*T0)``,  ``-8 <= m <= 8``.

Every resulting equation ``y^2=Q_T(mT+n)`` is square-normalized over ``Q``
and classified by degree and genus.  The leading terms cancel for ``m=+/-1``:
these 32 priority slices are genus-one quartics.  They receive one bounded
height search and two deterministic centered Mobius charts.  All subprocesses
are one-shot foreground process groups with strict caps and no retries.

Returned points are mapped back exactly.  A new parameter is called forced
only after its slice abscissa is checked against *all 21* generic-section
abscissas at that new parameter.  Conductor radical proxies rank the surviving
parameters; exact conductor/rank triage is deliberately left to the separate
candidate triage pass.

This is a bounded construction experiment, not a rank certificate.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
import hashlib
from math import comb, gcd, isqrt, lcm
import json
from pathlib import Path
import platform
import shlex
import sys
import time
from typing import Any, Iterable, Sequence

import sympy as sp

from certify_nagao_rank20_t5081 import (
    CHART_COUNT,
    CHART_HEIGHT,
    CONSTRUCTION,
    EXPECTED_TOTAL_ABSCISSAS,
    EXPECTED_UNIFORM_ABSCISSAS,
    EXPECTED_UNIFORM_SIGNED_POINTS,
    PARAMETER_T,
    UNIFORM_HEIGHT,
    exact_curve_data,
)
from ek_k3 import rational_to_string
from nagao_1994 import (
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    quartic_value,
    short_jacobian_coefficients,
)
from nagao_1994_section7 import (
    SECTION7_LINEAR_COMPANION_SECTIONS,
    SECTION7_QUADRATIC_COMPANION_SECTIONS,
    SECTION7_ROOTS,
    section7_primitive_quartic_coefficients,
)
from pari_bridge import pari_version
from search_elkies_klagsbrun_rank30 import (
    poly_add,
    poly_evaluate,
    poly_multiply,
)
from search_extra_points import signless_quartic_points
from search_nagao_rank20_t5081_neighborhood import conductor_radical_proxy
from search_nagao_rank21_t6793_skew import (
    map_run_points,
    optimized_cross_ratio_charts,
)
from search_nagao_rank21_t956_skew import search_original_quartic
from search_nagao_u42_skew_height import centered_unimodular_matrix, transform_binary_quartic
from triage_nagao_rank13_finalists import point_on_short_curve


Q = Fraction
T0 = PARAMETER_T
SLOPES = tuple(range(-8, 9))
PRIORITY_SLOPES = (-1, 1)
EXPECTED_RECONSTRUCTION_SHA256 = (
    "2079000c17c3a45043eff926cc594935355d897f85d633ba937e9eef7fbb9306"
)
# Exact signless output of the all-16 minimum-intercept m=+/-1 identity pilot
# at H=200000.  Its temporary stdout was not retained as an artifact; this
# script pins the complete reported list and independently accepts a parameter
# only when the normalized auxiliary value is an exact rational square.
PINNED_IDENTITY_HEIGHT_200000_PARAMETERS = tuple(
    Q(value)
    for value in (
        "-3251/1274", "16255/1456", "6465/91", "17935/182", "13108/91",
        "36135/182", "18295/91", "22208/91", "97883/182", "6115/6",
        "38525/821", "122205/1642", "98458/821", "145255/821", "180558/821",
        "4471/678", "11558/339", "53965/678", "45458/339", "92611/678",
        "121765/678", "25841/113", "131489/226", "6375/8186", "115745/4093",
        "-34484/307", "-19477/491", "1636/491", "5727/982", "68722/8347",
        "34361/2455", "29623/491", "103927/982", "65466/491", "-30515/1739",
        "34615/3478", "96432/1739", "195555/1739", "-14527/5009", "25/894",
        "38467/894", "20351/447", "89425/894", "65051/447", "154687/894",
        "115669/298", "-3315/2657", "139505/5314", "190646/2657", "24725/4106",
        "105774/2053", "-72193/658", "-1529/17", "-3754/47", "-3190/47",
        "-128873/2162", "-3795/94", "-44617/16685", "241/47", "42227/1034",
        "5605/94", "2920/47", "76490/799", "43765/446", "4941/47",
        "115945/1026", "7649/47", "42227/235", "154055/752", "-195713/39309",
        "-26162/5361", "-89650/1739", "-21051/3478", "84250/1739",
        "177195/3478", "-131490/1081", "-164609/2162", "-23390/1081",
        "-41375/2162", "51591/2162", "-181690/57017",
    )
)
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/elliptic_nagao_rank21_accidental_slices.json"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank21_accidental_slices.py"
)


# The coefficients of Q_T(x), first by ascending x-power and then by
# ascending T-power.  They are the exact explicit normalization in
# nagao_1994_section7.py, represented here as bivariate polynomial data so
# substitutions remain transparent and independently testable.
SECTION7_BIVARIATE_COEFFICIENTS: tuple[tuple[Fraction, ...], ...] = (
    (Q(557726319412900), Q(0), Q(23718659440), Q(0), Q(-910748), Q(0), Q(9)),
    (Q(-15923526472260), Q(0), Q(-322223958), Q(0), Q(6372)),
    (Q(156158560641), Q(0), Q(2005926), Q(0), Q(-18)),
    (Q(-623091366), Q(0), Q(-6372)),
    (Q(870426), Q(0), Q(9)),
)


@dataclass(frozen=True)
class NormalizedSlice:
    raw_coefficients: tuple[Fraction, ...]
    normalized_coefficients: tuple[int, ...]
    removed_square_coefficients: tuple[Fraction, ...]
    ordinate_constant_scale: Fraction
    raw_degree: int
    normalized_degree: int
    genus: int
    factor_degrees_and_exponents: tuple[tuple[int, int], ...]

    def normalized_value(self, parameter: Fraction) -> Fraction:
        return poly_evaluate(tuple(Q(value) for value in self.normalized_coefficients), Q(parameter))

    def removed_square_value(self, parameter: Fraction) -> Fraction:
        return poly_evaluate(self.removed_square_coefficients, Q(parameter))

    def original_ordinate(self, parameter: Fraction, normalized_ordinate: Fraction) -> Fraction:
        return (
            self.removed_square_value(parameter)
            * self.ordinate_constant_scale
            * Q(normalized_ordinate)
        )


@dataclass(frozen=True)
class Slice:
    accidental_index: int
    source_point: tuple[Fraction, Fraction]
    slope: int
    intercept: Fraction
    normalized: NormalizedSlice

    @property
    def identifier(self) -> str:
        sign = "p" if self.slope >= 0 else "m"
        return f"a{self.accidental_index + 1:02d}_s{sign}{abs(self.slope):02d}"

    def x_value(self, parameter: Fraction) -> Fraction:
        return Q(self.slope) * Q(parameter) + self.intercept


@dataclass(frozen=True)
class SearchPlan:
    identifier: str
    kind: str
    slice_index: int
    polynomial: tuple[Fraction, ...]
    height_specification: str
    matrix: tuple[int, int, int, int] | None
    homogenized_degree: int


def trim_polynomial(coefficients: Sequence[Fraction]) -> tuple[Fraction, ...]:
    answer = [Q(value) for value in coefficients]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def polynomial_digest(points: Iterable[tuple[Fraction, Fraction]]) -> str:
    text = "\n".join(
        f"{rational_to_string(x_value)},{rational_to_string(y_value)}"
        for x_value, y_value in sorted(points)
    )
    return hashlib.sha256(text.encode()).hexdigest()


def slice_polynomial(slope: Fraction, intercept: Fraction) -> tuple[Fraction, ...]:
    """Return ``Q_T(slope*T+intercept)`` in ascending T-power."""

    linear = (Q(intercept), Q(slope))
    power = (Q(1),)
    answer = (Q(0),)
    for coefficient_polynomial in SECTION7_BIVARIATE_COEFFICIENTS:
        answer = poly_add(answer, poly_multiply(coefficient_polynomial, power))
        power = poly_multiply(power, linear)
    return trim_polynomial(answer)


def _sympy_polynomial(coefficients: Sequence[Fraction], symbol: Any) -> sp.Poly:
    expression = sum(
        sp.Rational(Q(value).numerator, Q(value).denominator) * symbol**index
        for index, value in enumerate(coefficients)
    )
    return sp.Poly(expression, symbol, domain=sp.QQ)


def _fraction_coefficients(polynomial: sp.Poly) -> tuple[Fraction, ...]:
    degree = int(polynomial.degree())
    answer = []
    for index in range(degree + 1):
        value = polynomial.nth(index)
        answer.append(Q(int(value.p), int(value.q)))
    return trim_polynomial(answer)


def _square_part_and_kernel(value: int) -> tuple[int, int]:
    if value <= 0:
        raise ValueError("integer content must be positive")
    square_part = 1
    kernel = 1
    for prime, exponent in sp.factorint(value).items():
        square_part *= int(prime) ** (int(exponent) // 2)
        if int(exponent) % 2:
            kernel *= int(prime)
    if square_part**2 * kernel != value:
        raise AssertionError("integer square-kernel decomposition failed")
    return square_part, kernel


def normalize_slice(coefficients: Sequence[Fraction]) -> NormalizedSlice:
    """Remove polynomial squares and normalize the residual by a square.

    The returned data satisfy the exact identity

    ``raw(T) = (square(T)*constant_scale)^2 * normalized(T)``.
    """

    raw = trim_polynomial(coefficients)
    symbol = sp.symbols("T")
    polynomial = _sympy_polynomial(raw, symbol)
    constant, factors = sp.factor_list(polynomial)
    square_expression = sp.Integer(1)
    residual_expression = constant
    factor_records = []
    for factor, exponent in factors:
        exponent = int(exponent)
        factor_records.append((int(factor.degree()), exponent))
        square_expression *= factor.as_expr() ** (exponent // 2)
        if exponent % 2:
            residual_expression *= factor.as_expr()
    square_polynomial = sp.Poly(square_expression, symbol, domain=sp.QQ)
    residual = sp.Poly(residual_expression, symbol, domain=sp.QQ)
    residual_coefficients = _fraction_coefficients(residual)
    denominator = 1
    for value in residual_coefficients:
        denominator = lcm(denominator, value.denominator)
    integer_coefficients = tuple(
        int(value * denominator**2) for value in residual_coefficients
    )
    content = 0
    for value in integer_coefficients:
        content = gcd(content, abs(value))
    if content == 0:
        raise ValueError("the slice polynomial vanished")
    primitive = tuple(value // content for value in integer_coefficients)
    integer_square, integer_kernel = _square_part_and_kernel(content)
    normalized = tuple(integer_kernel * value for value in primitive)
    normalized = tuple(int(value) for value in trim_polynomial(tuple(Q(v) for v in normalized)))
    square_coefficients = _fraction_coefficients(square_polynomial)
    constant_scale = Q(integer_square, denominator)
    normalized_degree = len(normalized) - 1
    genus = max(0, (normalized_degree - 1) // 2)
    answer = NormalizedSlice(
        raw_coefficients=raw,
        normalized_coefficients=normalized,
        removed_square_coefficients=square_coefficients,
        ordinate_constant_scale=constant_scale,
        raw_degree=len(raw) - 1,
        normalized_degree=normalized_degree,
        genus=genus,
        factor_degrees_and_exponents=tuple(sorted(factor_records)),
    )
    # Independent exact polynomial replay.
    right = poly_multiply(
        poly_multiply(square_coefficients, square_coefficients),
        tuple(Q(value) for value in normalized),
    )
    right = tuple(constant_scale**2 * value for value in right)
    if trim_polynomial(right) != raw:
        raise AssertionError("the exact square normalization identity failed")
    return answer


def homogenized_transform(
    coefficients: Sequence[Fraction],
    matrix: Sequence[int],
    *,
    total_degree: int,
) -> tuple[Fraction, ...]:
    """Return ``(cU+d)^D f((aU+b)/(cU+d))`` exactly."""

    if total_degree <= 0 or total_degree % 2:
        raise ValueError("the homogeneous search degree must be positive and even")
    if len(coefficients) - 1 > total_degree or len(matrix) != 4:
        raise ValueError("invalid polynomial degree or Mobius matrix")
    a_value, b_value, c_value, d_value = (int(value) for value in matrix)
    if a_value * d_value - b_value * c_value == 0:
        raise ValueError("the Mobius matrix is singular")
    answer = [Q(0) for _ in range(total_degree + 1)]
    for power, coefficient in enumerate(coefficients):
        for left_power in range(power + 1):
            left = (
                comb(power, left_power)
                * a_value**left_power
                * b_value ** (power - left_power)
            )
            for right_power in range(total_degree - power + 1):
                right = (
                    comb(total_degree - power, right_power)
                    * c_value**right_power
                    * d_value ** (total_degree - power - right_power)
                )
                answer[left_power + right_power] += Q(coefficient) * left * right
    return trim_polynomial(answer)


def map_transformed_point(
    point: tuple[Fraction, Fraction],
    matrix: Sequence[int],
    *,
    total_degree: int,
) -> tuple[Fraction, Fraction] | None:
    parameter, ordinate = (Q(value) for value in point)
    a_value, b_value, c_value, d_value = (int(value) for value in matrix)
    denominator = c_value * parameter + d_value
    if denominator == 0:
        return None
    return (
        Q(a_value * parameter + b_value, denominator),
        ordinate / denominator ** (total_degree // 2),
    )


def generic_quartic_points(parameter: Fraction) -> tuple[tuple[str, tuple[Fraction, Fraction]], ...]:
    parameter = Q(parameter)
    records = [
        (f"visible-{index:02d}", point)
        for index, point in enumerate(primitive_visible_points(CONSTRUCTION, parameter))
    ]
    records.extend(
        (f"linear-{section.label}", section.point(parameter))
        for section in SECTION7_LINEAR_COMPANION_SECTIONS
    )
    records.extend(
        (f"quadratic-{section.label}", section.point(parameter))
        for section in SECTION7_QUADRATIC_COMPANION_SECTIONS
    )
    if len(records) != 21 or len({point[0] for _, point in records}) != 21:
        raise AssertionError("the 21 generic section abscissas collided")
    return tuple(records)


def generic_labels_for_x(parameter: Fraction, x_value: Fraction) -> tuple[str, ...]:
    return tuple(
        label for label, point in generic_quartic_points(parameter) if point[0] == Q(x_value)
    )


def generic_labels_for_jacobian_image(
    parameter: Fraction, quartic_point: tuple[Fraction, Fraction]
) -> tuple[str, ...]:
    """Return generic labels with the same Jacobian sign pair as a point."""

    parameter = Q(parameter)
    image = quartic_point_to_short_jacobian(CONSTRUCTION, parameter, quartic_point)
    coefficients = short_jacobian_coefficients(CONSTRUCTION, parameter)
    if not point_on_short_curve(coefficients, image):
        raise AssertionError("a slice point missed the specialized Jacobian")
    labels = []
    for label, generic_point in generic_quartic_points(parameter):
        generic_image = quartic_point_to_short_jacobian(
            CONSTRUCTION, parameter, generic_point
        )
        if not point_on_short_curve(coefficients, generic_image):
            raise AssertionError("a generic section missed the specialized Jacobian")
        # Equal short-model X means the same point up to sign.  It therefore
        # cannot supply a new Mordell--Weil direction.
        if generic_image[0] == image[0]:
            labels.append(label)
    return tuple(labels)


def reconstruct_rank20_quartic_points(
    *,
    uniform_timeout: float,
    chart_timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    quartic, _, _ = exact_curve_data()
    raw, uniform_process = search_original_quartic(
        quartic,
        str(UNIFORM_HEIGHT),
        timeout=uniform_timeout,
        stack_bytes=stack_bytes,
    )
    if uniform_process["status"] != "completed":
        raise RuntimeError(f"the pinned uniform reconstruction did not complete: {uniform_process}")
    uniform = signless_quartic_points(raw)
    if len(raw) != EXPECTED_UNIFORM_SIGNED_POINTS or len(uniform) != EXPECTED_UNIFORM_ABSCISSAS:
        raise AssertionError("the pinned uniform point set changed")
    by_x = {point[0]: point for point in uniform}
    charts = optimized_cross_ratio_charts(tuple(point[0] for point in uniform), count=CHART_COUNT)
    chart_records = []
    for index, chart in enumerate(charts, 1):
        transformed = transform_binary_quartic(quartic, chart.matrix)
        raw_chart, process = search_original_quartic(
            transformed,
            str(CHART_HEIGHT),
            timeout=chart_timeout,
            stack_bytes=stack_bytes,
        )
        before = len(by_x)
        if process["status"] == "completed":
            for point in map_run_points(quartic, raw_chart, chart.matrix):
                by_x.setdefault(point[0], point)
        chart_records.append(
            {
                "id": chart.identifier,
                "matrix": list(chart.matrix),
                **process,
                "new_global_abscissas": len(by_x) - before,
            }
        )
        if index % 12 == 0:
            print(f"reconstruction charts {index}/{len(charts)} union={len(by_x)}", flush=True)
    points = tuple(sorted(by_x.values()))
    if len(points) != EXPECTED_TOTAL_ABSCISSAS:
        raise AssertionError("the exact reconstructed union did not reach 34 abscissas")
    if polynomial_digest(points) != EXPECTED_RECONSTRUCTION_SHA256:
        raise AssertionError("the exact reconstructed quartic-point digest changed")
    if any(point[1] ** 2 != quartic_value(quartic, point[0]) for point in points):
        raise AssertionError("a reconstructed point missed the rank-20 quartic")
    return points, {
        "uniform": uniform_process,
        "uniform_signed_point_count": len(raw),
        "uniform_signless_point_count": len(uniform),
        "chart_records": chart_records,
        "completed_chart_count": sum(record["status"] == "completed" for record in chart_records),
        "timed_out_chart_count": sum(record["status"] == "timeout" for record in chart_records),
        "reconstructed_point_count": len(points),
        "reconstructed_point_sha256": polynomial_digest(points),
    }


def decontaminate_rank20_points(
    points: Sequence[tuple[Fraction, Fraction]],
) -> tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    generic = generic_quartic_points(T0)
    generic_by_x = {point[0]: label for label, point in generic}
    accidental = tuple(point for point in points if point[0] not in generic_by_x)
    found_labels = tuple(label for label, point in generic if point[0] in {p[0] for p in points})
    if len(accidental) != 16 or len(found_labels) != 18:
        raise AssertionError("the generic/accidental rank-20 split changed")
    return accidental, {
        "declared_generic_section_count": len(generic),
        "generic_abscissas_present_in_reconstruction": len(found_labels),
        "generic_labels_present": list(found_labels),
        "generic_labels_absent": [label for label, _ in generic if label not in found_labels],
        "accidental_point_count": len(accidental),
    }


def build_slices(
    accidental: Sequence[tuple[Fraction, Fraction]],
) -> tuple[Slice, ...]:
    answer = []
    for accidental_index, point in enumerate(accidental):
        x_value, y_value = point
        for slope in SLOPES:
            intercept = x_value - slope * T0
            normalized = normalize_slice(slice_polynomial(Q(slope), intercept))
            if poly_evaluate(normalized.raw_coefficients, T0) != y_value**2:
                raise AssertionError("a slice failed to replay its source point")
            removed = normalized.removed_square_value(T0) * normalized.ordinate_constant_scale
            if removed == 0:
                raise AssertionError("a source point lies at a removed square zero")
            normalized_y = y_value / removed
            if normalized_y**2 != normalized.normalized_value(T0):
                raise AssertionError("a normalized slice lost its source point")
            answer.append(Slice(accidental_index, point, slope, intercept, normalized))
    if len(answer) != len(accidental) * len(SLOPES):
        raise AssertionError("the complete slope grid changed size")
    return tuple(answer)


def select_minimum_intercept_priority_slices(
    slices: Sequence[Slice],
) -> tuple[Slice, ...]:
    """Choose the lower-height one of the two genus-one slopes per accident."""

    selected = []
    accidental_indices = sorted({item.accidental_index for item in slices})
    for accidental_index in accidental_indices:
        options = tuple(
            item
            for item in slices
            if item.accidental_index == accidental_index
            and item.slope in PRIORITY_SLOPES
        )
        if len(options) != 2 or any(item.normalized.genus != 1 for item in options):
            raise AssertionError("an accidental point lost a +/-1 genus-one option")
        selected.append(
            min(
                options,
                key=lambda item: (
                    max(abs(item.intercept.numerator), item.intercept.denominator),
                    abs(item.intercept.numerator).bit_length()
                    + item.intercept.denominator.bit_length(),
                    -item.slope,
                ),
            )
        )
    return tuple(selected)


def validate_known_linear_sections() -> tuple[dict[str, Any], ...]:
    records = []
    visible = primitive_visible_points(CONSTRUCTION, T0)
    for index, point in enumerate(visible):
        matches = tuple(
            (Q(root), slope)
            for root in SECTION7_ROOTS
            for slope in (-1, 1)
            if point[0] == Q(root) + slope * T0
        )
        if len(matches) != 1:
            raise AssertionError("a visible section has an unexpected abscissa")
        root, slope = matches[0]
        normalized = normalize_slice(slice_polynomial(Q(slope), root))
        if normalized.genus != 0 or normalized.normalized_degree != 0:
            raise AssertionError("a visible section slice did not normalize to a square")
        records.append(
            {
                "label": f"visible-{index:02d}",
                "slope": str(slope),
                "intercept": str(root),
                "normalized_degree": normalized.normalized_degree,
                "genus": normalized.genus,
            }
        )
    for section in SECTION7_LINEAR_COMPANION_SECTIONS:
        normalized = normalize_slice(slice_polynomial(section.slope, section.intercept))
        if normalized.genus != 0 or normalized.normalized_degree != 0:
            raise AssertionError("a linear companion slice did not normalize to a square")
        records.append(
            {
                "label": f"linear-{section.label}",
                "slope": str(section.slope),
                "intercept": str(section.intercept),
                "normalized_degree": normalized.normalized_degree,
                "genus": normalized.genus,
            }
        )
    if len(records) != 18:
        raise AssertionError("the known linear-section validation count changed")
    return tuple(records)


def build_priority_search_plans(
    priority: Sequence[Slice],
    *,
    identity_height: int,
    chart_height: int,
    chart_shifts: Sequence[int],
    include_identity: bool = True,
) -> tuple[SearchPlan, ...]:
    plans = []
    for slice_index, item in enumerate(priority):
        degree = 2 * item.normalized.genus + 2
        if degree < max(2, item.normalized.normalized_degree):
            raise AssertionError("the homogenized search degree is too small")
        polynomial = tuple(Q(value) for value in item.normalized.normalized_coefficients)
        if include_identity:
            plans.append(
                SearchPlan(
                    f"{item.identifier}_identity",
                    "identity",
                    slice_index,
                    polynomial,
                    str(identity_height),
                    None,
                    degree,
                )
            )
        for shift in chart_shifts:
            matrix = centered_unimodular_matrix(T0, shift)
            plans.append(
                SearchPlan(
                    f"{item.identifier}_centered_{shift}",
                    "centered_mobius",
                    slice_index,
                    homogenized_transform(polynomial, matrix, total_degree=degree),
                    str(chart_height),
                    matrix,
                    degree,
                )
            )
    return tuple(plans)


def map_search_points(
    plan: SearchPlan,
    item: Slice,
    raw_points: Iterable[tuple[Fraction, Fraction]],
) -> tuple[tuple[tuple[Fraction, Fraction], ...], int]:
    mapped = []
    poles = 0
    for point in raw_points:
        if poly_evaluate(plan.polynomial, point[0]) != point[1] ** 2:
            raise AssertionError("PARI returned a point off a transformed slice")
        if plan.matrix is None:
            original = Q(point[0]), Q(point[1])
        else:
            original = map_transformed_point(
                point, plan.matrix, total_degree=plan.homogenized_degree
            )
            if original is None:
                poles += 1
                continue
        if original[1] ** 2 != item.normalized.normalized_value(original[0]):
            raise AssertionError("a chart point missed the normalized auxiliary curve")
        mapped.append(original)
    return tuple(mapped), poles


def rational_square_root(value: Fraction) -> Fraction | None:
    value = Q(value)
    if value < 0:
        return None
    numerator_root = isqrt(value.numerator)
    denominator_root = isqrt(value.denominator)
    if numerator_root**2 != value.numerator or denominator_root**2 != value.denominator:
        return None
    return Q(numerator_root, denominator_root)


def replay_pinned_identity_parameters(
    priority: Sequence[Slice],
) -> tuple[tuple[set[tuple[Fraction, Fraction]], ...], dict[str, Any]]:
    """Associate every pinned H=200000 parameter by exact square testing."""

    by_slice: list[set[tuple[Fraction, Fraction]]] = [set() for _ in priority]
    association_records = []
    for parameter in PINNED_IDENTITY_HEIGHT_200000_PARAMETERS:
        matches = []
        for slice_index, item in enumerate(priority):
            root = rational_square_root(item.normalized.normalized_value(parameter))
            if root is None:
                continue
            matches.append(slice_index)
            by_slice[slice_index].add((parameter, root))
            if root:
                by_slice[slice_index].add((parameter, -root))
        if not matches:
            raise AssertionError(f"pinned pilot parameter {parameter} missed every slice")
        association_records.append(
            {
                "T": str(parameter),
                "matching_priority_slice_indices_zero_based": matches,
                "exact_square_replay": True,
            }
        )
    # T0 was the known base point on each auxiliary curve and is replayed
    # separately rather than duplicated in the signless pilot list.
    for slice_index, item in enumerate(priority):
        root = rational_square_root(item.normalized.normalized_value(T0))
        if root is None:
            raise AssertionError("T0 ceased to be a point on a priority slice")
        by_slice[slice_index].add((T0, root))
        if root:
            by_slice[slice_index].add((T0, -root))
    return tuple(by_slice), {
        "height": 200_000,
        "signless_parameter_count_excluding_T0": len(
            PINNED_IDENTITY_HEIGHT_200000_PARAMETERS
        ),
        "parameters_excluding_T0": [
            str(value) for value in PINNED_IDENTITY_HEIGHT_200000_PARAMETERS
        ],
        "association_records": association_records,
        "T0_replayed_separately_on_every_priority_slice": True,
        "all_values_verified_by_exact_rational_square_tests": True,
        "enumeration_stdout_was_temporary_and_is_not_a_portable_input": True,
    }


def slice_manifest_record(item: Slice) -> dict[str, Any]:
    return {
        "id": item.identifier,
        "accidental_index_one_based": item.accidental_index + 1,
        "source_x": str(item.source_point[0]),
        "source_y": str(item.source_point[1]),
        "slope": item.slope,
        "intercept": str(item.intercept),
        "raw_degree": item.normalized.raw_degree,
        "normalized_degree": item.normalized.normalized_degree,
        "genus": item.normalized.genus,
        "raw_coefficients_ascending": [str(value) for value in item.normalized.raw_coefficients],
        "normalized_integer_coefficients_ascending": list(
            item.normalized.normalized_coefficients
        ),
        "removed_square_coefficients_ascending": [
            str(value) for value in item.normalized.removed_square_coefficients
        ],
        "ordinate_constant_scale": str(item.normalized.ordinate_constant_scale),
        "factor_degrees_and_exponents": [
            list(value) for value in item.normalized.factor_degrees_and_exponents
        ],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--uniform-timeout", type=float, default=40.0)
    parser.add_argument("--reconstruction-chart-timeout", type=float, default=3.0)
    parser.add_argument("--slice-timeout", type=float, default=30.0)
    parser.add_argument("--identity-height", type=int, default=200_000)
    parser.add_argument("--chart-height", type=int, default=50_000)
    parser.add_argument("--chart-shifts", default="0,-1")
    parser.add_argument("--rerun-identity", action="store_true")
    parser.add_argument("--proxy-trial-prime-bound", type=int, default=251)
    parser.add_argument("--proxy-frontier-count", type=int, default=32)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    return parser


def parse_shifts(text: str) -> tuple[int, ...]:
    values = tuple(int(value) for value in text.split(",") if value)
    if not values or len(set(values)) != len(values):
        raise ValueError("chart shifts must be a nonempty distinct integer list")
    return values


def main() -> None:
    args = build_parser().parse_args()
    if not 0 < args.uniform_timeout <= 60:
        raise SystemExit("--uniform-timeout must lie in (0,60]")
    if not 0 < args.reconstruction_chart_timeout <= 60 or not 0 < args.slice_timeout <= 60:
        raise SystemExit("per-process timeouts must lie in (0,60]")
    if min(
        args.identity_height,
        args.chart_height,
        args.proxy_trial_prime_bound,
        args.proxy_frontier_count,
    ) <= 0:
        raise SystemExit("heights, proxy bounds, and frontier size must be positive")
    if args.stack_bytes < 64_000_000:
        raise SystemExit("--stack-bytes must be at least 64MB")
    try:
        chart_shifts = parse_shifts(args.chart_shifts)
    except ValueError as error:
        raise SystemExit(str(error)) from error

    started = time.monotonic()
    points, reconstruction = reconstruct_rank20_quartic_points(
        uniform_timeout=args.uniform_timeout,
        chart_timeout=args.reconstruction_chart_timeout,
        stack_bytes=args.stack_bytes,
    )
    accidental, decontamination = decontaminate_rank20_points(points)
    known_validation = validate_known_linear_sections()
    slices = build_slices(accidental)
    all_genus_one = tuple(item for item in slices if item.slope in PRIORITY_SLOPES)
    if len(all_genus_one) != 32 or any(item.normalized.genus != 1 for item in all_genus_one):
        raise AssertionError("the 32 m=+/-1 slices are not all genus one")
    priority = select_minimum_intercept_priority_slices(slices)
    if len(priority) != 16:
        raise AssertionError("the minimum-intercept priority set must have 16 slices")
    classification_counts = Counter(
        (item.normalized.raw_degree, item.normalized.normalized_degree, item.normalized.genus)
        for item in slices
    )
    manifest_payload = [slice_manifest_record(item) for item in slices]
    manifest_sha256 = hashlib.sha256(
        json.dumps(manifest_payload, separators=(",", ":")).encode()
    ).hexdigest()
    print(
        f"slice classification: total={len(slices)} genus1={len(all_genus_one)} "
        f"priority_min_intercept={len(priority)} "
        f"classes={dict(classification_counts)}",
        flush=True,
    )

    plans = build_priority_search_plans(
        priority,
        identity_height=args.identity_height,
        chart_height=args.chart_height,
        chart_shifts=chart_shifts,
        include_identity=args.rerun_identity,
    )
    replayed, pinned_pilot = replay_pinned_identity_parameters(priority)
    mapped_by_slice: list[set[tuple[Fraction, Fraction]]] = [
        set(points) for points in replayed
    ]
    run_records = []
    for index, plan in enumerate(plans, 1):
        raw, process = search_original_quartic(
            plan.polynomial,
            plan.height_specification,
            timeout=args.slice_timeout,
            stack_bytes=args.stack_bytes,
        )
        mapped, poles = map_search_points(plan, priority[plan.slice_index], raw)
        mapped_by_slice[plan.slice_index].update(mapped)
        run_records.append(
            {
                "id": plan.identifier,
                "kind": plan.kind,
                "slice_index_zero_based": plan.slice_index,
                "height_specification": plan.height_specification,
                "matrix": None if plan.matrix is None else list(plan.matrix),
                **process,
                "mapped_exact_point_count_including_signs": len(mapped),
                "chart_pole_count": poles,
                "retried": False,
            }
        )
        if index % 16 == 0 or process["status"] != "completed":
            print(
                f"auxiliary searches {index}/{len(plans)} "
                f"status={process['status']}",
                flush=True,
            )

    parameter_records: dict[Fraction, dict[str, Any]] = {}
    slice_point_records = []
    for slice_index, (item, auxiliary_points) in enumerate(zip(priority, mapped_by_slice)):
        exact_affine_count = len(auxiliary_points)
        positive_rank_by_cardinality = exact_affine_count > 16
        signless_parameters = {point[0] for point in auxiliary_points}
        slice_point_records.append(
            {
                "slice_index_zero_based": slice_index,
                "id": item.identifier,
                "slope": item.slope,
                "intercept": str(item.intercept),
                "normalized_integer_coefficients_ascending": list(
                    item.normalized.normalized_coefficients
                ),
                "exact_distinct_affine_auxiliary_point_count": exact_affine_count,
                "exact_distinct_parameter_count": len(signless_parameters),
                "positive_rank_certified_by_point_count": positive_rank_by_cardinality,
                "positive_rank_argument": (
                    "a rational base point identifies the genus-one slice with its "
                    "Jacobian; more than 16 rational affine points exceed Mazur's "
                    "maximum rational torsion cardinality"
                    if positive_rank_by_cardinality
                    else None
                ),
            }
        )
        for parameter, normalized_ordinate in sorted(auxiliary_points):
            if abs(parameter) == T0:
                continue
            original_y = item.normalized.original_ordinate(parameter, normalized_ordinate)
            x_value = item.x_value(parameter)
            quartic = section7_primitive_quartic_coefficients(parameter)
            if original_y**2 != quartic_value(quartic, x_value):
                raise AssertionError("a recovered slice point missed the specialized quartic")
            quartic_labels = generic_labels_for_x(parameter, x_value)
            jacobian_labels = generic_labels_for_jacobian_image(
                parameter, (x_value, original_y)
            )
            labels = tuple(
                sorted(
                    {
                        *(f"quartic-x:{label}" for label in quartic_labels),
                        *(f"jacobian-sign-pair:{label}" for label in jacobian_labels),
                    }
                )
            )
            record = parameter_records.setdefault(
                parameter,
                {
                    "T": str(parameter),
                    "parameter_height": max(abs(parameter.numerator), parameter.denominator),
                    "generic_intersection_sources": set(),
                    "forced_non_generic_sources": set(),
                    "forced_points": {},
                },
            )
            source = item.identifier
            if labels:
                record["generic_intersection_sources"].add(
                    f"{source}:{','.join(labels)}"
                )
            else:
                record["forced_non_generic_sources"].add(source)
                record["forced_points"].setdefault(
                    str(x_value),
                    {
                        "x": str(x_value),
                        "y": str(original_y),
                        "exact_quartic_membership_checked": True,
                    },
                )

    parameter_output = []
    forced_parameters = []
    for parameter, record in sorted(parameter_records.items()):
        forced = bool(record["forced_non_generic_sources"])
        output = {
            "T": record["T"],
            "parameter_height": record["parameter_height"],
            "classification": (
                "forced_non_generic_quartic_point"
                if forced
                else "only_intersections_with_21_generic_sections"
            ),
            "generic_intersection_sources": sorted(record["generic_intersection_sources"]),
            "forced_non_generic_sources": sorted(record["forced_non_generic_sources"]),
            "forced_points": list(record["forced_points"].values()),
            "exact_all_21_generic_abscissas_checked": True,
            "exact_all_21_generic_jacobian_sign_pairs_checked": True,
        }
        if forced:
            try:
                proxy = conductor_radical_proxy(
                    parameter, trial_prime_bound=args.proxy_trial_prime_bound
                )
            except ValueError as error:
                proxy = {"status": "singular", "error": str(error)}
            else:
                proxy = {"status": "completed", **proxy}
            output["conductor_radical_proxy"] = proxy
            forced_parameters.append((parameter, output))
        parameter_output.append(output)

    def proxy_key(item: tuple[Fraction, dict[str, Any]]) -> tuple[Any, ...]:
        parameter, record = item
        proxy = record["conductor_radical_proxy"]
        return (
            proxy.get("log_radical_upper_proxy", float("inf")),
            max(abs(parameter.numerator), parameter.denominator),
            parameter,
        )

    forced_parameters.sort(key=proxy_key)
    proxy_frontier = [record for _, record in forced_parameters[: args.proxy_frontier_count]]
    completed_runs = sum(record["status"] == "completed" for record in run_records)
    artifact = {
        "schema_version": 1,
        "artifact_kind": "bounded_section7_accidental_linear_slice_construction",
        "status": "bounded_constructive_slice_search_complete",
        "claim_scope": {
            "exact": (
                "rank-20 point reconstruction, generic-section decontamination, "
                "slice equations and normalizations, all returned point maps, "
                "and forced-vs-generic classification"
            ),
            "bounded": (
                "only m=+/-1 auxiliary curves receive the declared PARI boxes; "
                "conductor values are proxies, not exact conductors"
            ),
            "rank_certificate": False,
        },
        "input_fiber": {
            "constructor_T0": str(T0),
            "roots": list(SECTION7_ROOTS),
            "certified_rank_lower_bound_from_separate_input": 20,
        },
        "reconstruction": reconstruction,
        "decontamination_at_T0": {
            **decontamination,
            "accidental_points": [
                {"x": str(point[0]), "y": str(point[1])} for point in accidental
            ],
        },
        "known_section_slice_validation": {
            "count": len(known_validation),
            "records": list(known_validation),
            "all_normalized_to_genus_zero_constant": True,
        },
        "complete_slice_classification": {
            "slope_range": [SLOPES[0], SLOPES[-1]],
            "accidental_point_count": len(accidental),
            "slice_count": len(slices),
            "classification_counts": [
                {
                    "raw_degree": key[0],
                    "normalized_degree": key[1],
                    "genus": key[2],
                    "count": count,
                }
                for key, count in sorted(classification_counts.items())
            ],
            "manifest_sha256": manifest_sha256,
            "manifest_stored_inline": False,
            "m_plus_minus_one_genus_one_slice_count": len(all_genus_one),
            "priority_minimum_intercept_slice_count": len(priority),
            "priority_selection": (
                "for each accidental point choose m in {-1,+1} minimizing "
                "projective height of n=x0-m*T0"
            ),
            "all_m_plus_minus_one_slices_genus_one": True,
        },
        "auxiliary_search_budget": {
            "identity_height": args.identity_height,
            "identity_enumeration_rerun": args.rerun_identity,
            "pinned_identity_pilot_replayed_exactly": True,
            "centered_chart_height": args.chart_height,
            "centered_chart_shifts": list(chart_shifts),
            "per_call_timeout_seconds": args.slice_timeout,
            "stack_bytes": args.stack_bytes,
            "declared_call_count": len(plans),
            "one_pass_no_retry": True,
        },
        "auxiliary_search": {
            "pinned_identity_pilot": pinned_pilot,
            "completed_call_count": completed_runs,
            "timed_out_call_count": sum(record["status"] == "timeout" for record in run_records),
            "other_noncompleted_call_count": len(plans)
            - completed_runs
            - sum(record["status"] == "timeout" for record in run_records),
            "run_records": run_records,
            "slice_point_records": slice_point_records,
            "positive_rank_slice_count_by_point_cardinality": sum(
                record["positive_rank_certified_by_point_count"]
                for record in slice_point_records
            ),
        },
        "new_parameter_decontamination": {
            "distinct_new_parameter_count_before_generic_filter": len(parameter_output),
            "generic_only_parameter_count": sum(
                record["classification"] == "only_intersections_with_21_generic_sections"
                for record in parameter_output
            ),
            "forced_non_generic_parameter_count": len(forced_parameters),
            "records": parameter_output,
        },
        "conductor_proxy_frontier": {
            "trial_prime_bound": args.proxy_trial_prime_bound,
            "retained_count": len(proxy_frontier),
            "selection": "smallest log radical upper proxy, then parameter height and T",
            "records": proxy_frontier,
            "exact_conductor_not_claimed": True,
        },
        "target": {
            "rank": 21,
            "strict_log_conductor_bound": "182.72",
            "hit": False,
            "reason": (
                "this artifact constructs and decontaminates candidates; exact "
                "conductor and rank certification are separate gated passes"
            ),
        },
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "actual_command": " ".join(
                shlex.quote(part) for part in [sys.executable, *sys.argv]
            ),
            "python": platform.python_version(),
            "sympy": sp.__version__,
            "pari_gp": pari_version(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "wall_seconds": time.monotonic() - started,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")
    print(
        f"new_T={len(parameter_output)} generic_only="
        f"{artifact['new_parameter_decontamination']['generic_only_parameter_count']} "
        f"forced={len(forced_parameters)} positive_rank_slices="
        f"{artifact['auxiliary_search']['positive_rank_slice_count_by_point_cardinality']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
