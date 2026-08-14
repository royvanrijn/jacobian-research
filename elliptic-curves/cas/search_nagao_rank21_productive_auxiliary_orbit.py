#!/usr/bin/env python3
"""Search a Mordell--Weil orbit forced by Nagao's published rank-21 fibre.

Nagao's record specialization uses the six-root construction
``(399,380,352,47,4,0)`` at constructor parameter ``T=14721/188``.  This
script first reconstructs exact quartic preimages of all 21 printed points.
Eleven preimages are not among the twelve visible Mestre sections.  The most
productive bounded slice through those points is

``x = -T + 57361/139``.

Its auxiliary equation is a genus-one quartic.  Rather than repeatedly
raising a naive ``hyperellratpoints`` box, we turn it into a pointed
generalized Weierstrass curve, find an eight-dimensional numerical basis from
the complete ``H=200000`` replay, saturate that basis at small primes, and
enumerate the exact coefficient box ``|c_i|<=2, sum |c_i|<=8``.  Every group
combination is inverted exactly to a specialization parameter.  The complete
orbit receives a small-prime discriminant-radical proxy; all plausible
survivors receive exact conductor calculations.

Numerical height rank is used only to choose generators.  Exact group-law
identities and conductor computations are reproducible computations, not a
proof that the auxiliary curve has rank exactly eight or that an unsearched
coefficient vector cannot produce a better specialization.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import isqrt, log
from pathlib import Path
import platform
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

import sympy as sp

from alternate_quartic_covers import short_add
from ek_k3 import primes_up_to, rational_square_root, rational_to_string
from extend_nagao_u42_frontier import saturate_exact_basis
from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK21_CONSTRUCTION,
    RANK21_CONSTRUCTOR_PARAMETER,
    RANK21_PUBLISHED_MODEL,
    RANK21_PUBLISHED_POINTS,
    primitive_quartic_coefficients,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    quartic_value,
    rank21_short_jacobian_coefficients,
    short_jacobian_coefficients,
)
from pari_bridge import minimal_curve_data, pari_version
from search_elkies_klagsbrun_rank30 import poly_add, poly_multiply
from search_extra_points import signless_quartic_points
from search_nagao_rank21_accidental_slices import NormalizedSlice, normalize_slice
from search_nagao_rank21_neighborhood import DISCRIMINANT_POLYNOMIAL
from search_nagao_rank21_t956_skew import search_original_quartic
from search_nagao_section7_auxiliary_jacobians import (
    translate_polynomial,
    weierstrass_add,
    weierstrass_multiply,
)
from triage_nagao_rank13_finalists import (
    height_matrix_replay,
    point_digest,
    point_on_short_curve,
    stable_height_rank,
)


Q = Fraction
T0 = RANK21_CONSTRUCTOR_PARAMETER
TARGET_LOG_CONDUCTOR = Decimal("182.72")
MODEL_CHANGE_SHORT_TO_MINIMAL = (
    Q(3, 8836),
    Q(15, 312299584),
    Q(3, 17672),
    Q(27, 1379739562112),
)
PRODUCTIVE_SOURCE_X = Q(8737649, 26132)
PRODUCTIVE_SLOPE = -1
PRODUCTIVE_INTERCEPT = Q(57361, 139)
ORBIT_MAXIMUM_ABSOLUTE_COEFFICIENT = 2
ORBIT_MAXIMUM_L1_NORM = 8
PROXY_PRIME_BOUND = 1000
PROXY_EXACT_CONDUCTOR_GATE = 190.0
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/"
    "elliptic_nagao_rank21_productive_auxiliary_orbit.json"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank21_productive_auxiliary_orbit.py"
)


# Ascending T-coefficients for the ascending x-coefficients (e,d,c,b,a).
# The script independently checks these values against the generic constructor
# at more points than were needed to interpolate them.
RANK21_QUARTIC_COEFFICIENT_POLYNOMIALS: tuple[tuple[Fraction, ...], ...] = (
    (Q(197804341504), Q(0), Q(-5788761232), Q(0), Q(4306112), Q(0), Q(-23)),
    (Q(-37842168672), Q(0), Q(98440030), Q(0), Q(-18124)),
    (Q(6945772145), Q(0), Q(-3814698), Q(0), Q(46)),
    (Q(-35986906), Q(0), Q(18124)),
    (Q(47342), Q(0), Q(-23)),
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def rational_record(value: Fraction) -> str:
    return rational_to_string(Q(value))


def rational_pair_record(point: tuple[Fraction, Fraction]) -> dict[str, str]:
    return {"x": rational_record(point[0]), "y": rational_record(point[1])}


def polynomial_value(coefficients: Sequence[Fraction], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * Q(value) + Q(coefficient)
    return answer


def polynomial_digest(values: Iterable[Fraction]) -> str:
    text = "\n".join(rational_record(value) for value in values)
    return hashlib.sha256(text.encode()).hexdigest()


def orbit_record_digest(
    records: Iterable[tuple[Fraction, tuple[int, ...]]],
) -> str:
    text = "\n".join(
        f"{rational_record(parameter)}|{','.join(map(str, vector))}"
        for parameter, vector in records
    )
    return hashlib.sha256(text.encode()).hexdigest()


def validate_bivariate_quartic() -> None:
    for parameter in (Q(-9), Q(-1), Q(1, 2), Q(1), Q(3, 2), Q(8), T0):
        replay = tuple(
            polynomial_value(coefficient_polynomial, parameter)
            for coefficient_polynomial in RANK21_QUARTIC_COEFFICIENT_POLYNOMIALS
        )
        expected = primitive_quartic_coefficients(RANK21_CONSTRUCTION, parameter)
        if replay != expected:
            raise AssertionError("the rank-21 bivariate quartic replay failed")


def slice_polynomial(slope: int, intercept: Fraction) -> tuple[Fraction, ...]:
    if slope not in (-1, 1):
        raise ValueError("only the degree-cancelling slopes +/-1 are supported")
    linear = (Q(intercept), Q(slope))
    power = (Q(1),)
    answer = (Q(0),)
    for coefficient_polynomial in RANK21_QUARTIC_COEFFICIENT_POLYNOMIALS:
        answer = poly_add(answer, poly_multiply(coefficient_polynomial, power))
        power = poly_multiply(power, linear)
    while len(answer) > 1 and answer[-1] == 0:
        answer = answer[:-1]
    return tuple(Q(value) for value in answer)


def minimum_to_short(point: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    u_value, r_value, s_value, t_value = MODEL_CHANGE_SHORT_TO_MINIMAL
    x_value, y_value = (Q(value) for value in point)
    answer = (
        u_value**2 * x_value + r_value,
        u_value**3 * y_value + s_value * u_value**2 * x_value + t_value,
    )
    if not point_on_short_curve(rank21_short_jacobian_coefficients(), answer):
        raise AssertionError("the published-to-short model change failed")
    return answer


def short_to_minimum(point: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    u_value, r_value, s_value, t_value = MODEL_CHANGE_SHORT_TO_MINIMAL
    x_value, y_value = (Q(value) for value in point)
    answer = (
        (x_value - r_value) / u_value**2,
        (y_value - s_value * (x_value - r_value) - t_value) / u_value**3,
    )
    if not point_on_short_curve(tuple(Q(value) for value in RANK21_PUBLISHED_MODEL), answer):
        raise AssertionError("the short-to-published model change failed")
    return answer


def short_negate(point: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return point[0], -point[1]


def rational_linear_roots(expression: sp.Expr, symbol: sp.Symbol) -> tuple[Fraction, ...]:
    roots: set[Fraction] = set()
    for factor, _ in sp.factor_list(expression)[1]:
        if sp.degree(factor, symbol) != 1:
            continue
        root = sp.solve(factor, symbol)[0]
        if root.is_Rational:
            roots.add(Q(int(root.p), int(root.q)))
    return tuple(sorted(roots))


def covariant_x_polynomial(target_x: Fraction) -> tuple[Fraction, ...]:
    e_value, d_value, c_value, b_value, a_value = primitive_quartic_coefficients(
        RANK21_CONSTRUCTION, T0
    )
    g_values = (
        b_value**2 / 16 - a_value * c_value / 6,
        b_value * c_value / 12 - a_value * d_value / 2,
        c_value**2 / 12 - b_value * d_value / 8 - a_value * e_value,
        c_value * d_value / 12 - b_value * e_value / 2,
        d_value**2 / 16 - c_value * e_value / 6,
    )
    return tuple(
        Q(target_x) * coefficient - 36 * g_value
        for coefficient, g_value in zip(
            (a_value, b_value, c_value, d_value, e_value),
            g_values,
            strict=True,
        )
    )


def invert_covariant_target(
    target: tuple[Fraction, Fraction],
) -> tuple[Fraction, Fraction]:
    symbol = sp.symbols("x")
    coefficients = covariant_x_polynomial(target[0])
    expression = sum(
        sp.Rational(value.numerator, value.denominator) * symbol ** (4 - index)
        for index, value in enumerate(coefficients)
    )
    quartic = primitive_quartic_coefficients(RANK21_CONSTRUCTION, T0)
    answers: list[tuple[Fraction, Fraction]] = []
    for x_value in rational_linear_roots(expression, symbol):
        y_value = rational_square_root(quartic_value(quartic, x_value))
        if y_value is None:
            continue
        for signed_y in {y_value, -y_value}:
            point = (x_value, signed_y)
            if quartic_point_to_short_jacobian(
                RANK21_CONSTRUCTION, T0, point
            ) == target:
                answers.append(point)
    if len(answers) != 1:
        raise ValueError(f"expected one rational covariant preimage, found {len(answers)}")
    return answers[0]


def reconstruct_with_convention(
    base_index: int, sign: int
) -> tuple[tuple[Fraction, Fraction], ...]:
    visible = primitive_visible_points(RANK21_CONSTRUCTION, T0)
    coefficients = rank21_short_jacobian_coefficients()
    offset = quartic_point_to_short_jacobian(
        RANK21_CONSTRUCTION, T0, visible[base_index]
    )
    answer = []
    for published in RANK21_PUBLISHED_POINTS:
        short_point = minimum_to_short(published)
        doubled = short_add(coefficients, short_point, short_point)
        if doubled is None:
            raise AssertionError("a published point doubled to zero")
        if sign == -1:
            doubled = short_negate(doubled)
        target = short_add(coefficients, offset, doubled)
        if target is None:
            raise AssertionError("a covariant target became infinity")
        answer.append(invert_covariant_target(target))
    return tuple(answer)


def select_reconstruction_convention() -> tuple[
    tuple[tuple[Fraction, Fraction], ...], tuple[dict[str, Any], ...]
]:
    visible = primitive_visible_points(RANK21_CONSTRUCTION, T0)
    visible_set = set(visible)
    visible_x = {point[0] for point in visible}
    trials: list[dict[str, Any]] = []
    successful: dict[tuple[int, int], tuple[tuple[Fraction, Fraction], ...]] = {}
    for base_index in range(len(visible)):
        for sign in (1, -1):
            try:
                preimages = reconstruct_with_convention(base_index, sign)
            except (ValueError, AssertionError) as error:
                trials.append(
                    {
                        "base_index_one_based": base_index + 1,
                        "sign": sign,
                        "status": "no-complete-rational-reconstruction",
                        "error": str(error),
                        "exact_oriented_generic_matches": 0,
                        "generic_abscissa_matches": 0,
                    }
                )
                continue
            successful[(base_index, sign)] = preimages
            trials.append(
                {
                    "base_index_one_based": base_index + 1,
                    "sign": sign,
                    "status": "completed",
                    "exact_oriented_generic_matches": sum(
                        point in visible_set for point in preimages
                    ),
                    "generic_abscissa_matches": sum(
                        point[0] in visible_x for point in preimages
                    ),
                }
            )
    maximum = max(int(record["generic_abscissa_matches"]) for record in trials)
    winners = [
        record for record in trials if int(record["generic_abscissa_matches"]) == maximum
    ]
    if len(winners) != 1:
        raise AssertionError("the published reconstruction convention is not unique")
    winner = winners[0]
    key = int(winner["base_index_one_based"]) - 1, int(winner["sign"])
    if key != (0, 1) or maximum != 10:
        raise AssertionError("the pinned rank-21 reconstruction convention changed")
    if int(winner["exact_oriented_generic_matches"]) != 6:
        raise AssertionError("the oriented generic match count changed")
    return successful[key], tuple(trials)


@dataclass(frozen=True)
class Slice:
    accidental_index: int
    source_point: tuple[Fraction, Fraction]
    slope: int
    intercept: Fraction
    normalized: NormalizedSlice

    @property
    def identifier(self) -> str:
        return f"a{self.accidental_index + 1:02d}_s{'p' if self.slope > 0 else 'm'}01"

    def raw_x(self, parameter: Fraction) -> Fraction:
        return Q(self.slope) * Q(parameter) + self.intercept


def build_slices(
    accidental: Sequence[tuple[Fraction, Fraction]],
) -> tuple[Slice, ...]:
    answer = []
    for accidental_index, source_point in enumerate(accidental):
        for slope in (-1, 1):
            intercept = source_point[0] - slope * T0
            normalized = normalize_slice(slice_polynomial(slope, intercept))
            if normalized.normalized_degree != 4 or normalized.genus != 1:
                raise AssertionError("an accidental +/-1 slice ceased to be genus one")
            raw_value = polynomial_value(normalized.raw_coefficients, T0)
            if raw_value != source_point[1] ** 2:
                raise AssertionError("an accidental slice lost its source point")
            answer.append(
                Slice(accidental_index, source_point, slope, intercept, normalized)
            )
    if len(answer) != 22:
        raise AssertionError("the accidental slice count changed")
    return tuple(answer)


def normalized_to_raw_point(
    item: Slice, point: tuple[Fraction, Fraction]
) -> tuple[Fraction, Fraction]:
    parameter, normalized_ordinate = point
    raw_ordinate = item.normalized.original_ordinate(parameter, normalized_ordinate)
    if raw_ordinate**2 != polynomial_value(item.normalized.raw_coefficients, parameter):
        raise AssertionError("a normalized slice point failed raw replay")
    return parameter, raw_ordinate


def visible_abscissas(parameter: Fraction) -> set[Fraction]:
    parameter = Q(parameter)
    return {
        root + sign * parameter
        for root in RANK21_CONSTRUCTION.roots
        for sign in (-1, 1)
    }


@dataclass(frozen=True)
class PointedQuartic:
    base_parameter: Fraction
    base_ordinate: Fraction
    quartic_coefficients: tuple[Fraction, ...]
    shifted_coefficients: tuple[Fraction, ...]
    weierstrass_coefficients: tuple[Fraction, ...]

    @classmethod
    def from_slice(cls, item: Slice) -> "PointedQuartic":
        removed = (
            item.normalized.removed_square_value(T0)
            * item.normalized.ordinate_constant_scale
        )
        if removed == 0:
            raise AssertionError("the productive source lies at a removed-square zero")
        base_ordinate = item.source_point[1] / removed
        quartic = tuple(Q(value) for value in item.normalized.normalized_coefficients)
        shifted = translate_polynomial(quartic, T0)
        if len(shifted) != 5 or shifted[0] != base_ordinate**2:
            raise AssertionError("the pointed auxiliary quartic lost its base point")
        _, d_value, c_value, b_value, a_value = shifted
        q_value = base_ordinate
        coefficients = (
            d_value / q_value,
            c_value - d_value**2 / (4 * q_value**2),
            2 * q_value * b_value,
            -4 * q_value**2 * a_value,
            a_value * (d_value**2 - 4 * q_value**2 * c_value),
        )
        return cls(T0, base_ordinate, quartic, shifted, coefficients)

    def quartic_value(self, parameter: Fraction) -> Fraction:
        return polynomial_value(self.quartic_coefficients, parameter)

    def forward(
        self, point: tuple[Fraction, Fraction]
    ) -> tuple[Fraction, Fraction]:
        parameter, ordinate = (Q(value) for value in point)
        if ordinate**2 != self.quartic_value(parameter):
            raise ValueError("the point missed the auxiliary quartic")
        u_value = parameter - self.base_parameter
        if u_value == 0:
            raise ValueError("the chosen base fibre maps exceptionally")
        _, d_value, c_value, _, _ = self.shifted_coefficients
        q_value = self.base_ordinate
        x_value = (
            2 * q_value * (ordinate + q_value) + d_value * u_value
        ) / u_value**2
        y_value = (
            4 * q_value**2 * (ordinate + q_value)
            + 2 * q_value * (d_value * u_value + c_value * u_value**2)
            - d_value**2 * u_value**2 / (2 * q_value)
        ) / u_value**3
        answer = x_value, y_value
        if not point_on_short_curve(self.weierstrass_coefficients, answer):
            raise AssertionError("the pointed-quartic forward map failed")
        return answer

    def inverse(
        self, point: tuple[Fraction, Fraction] | None
    ) -> tuple[Fraction, Fraction] | None:
        if point is None:
            return self.base_parameter, self.base_ordinate
        if not point_on_short_curve(self.weierstrass_coefficients, point):
            raise ValueError("the point missed the auxiliary Weierstrass curve")
        x_value, y_value = point
        if y_value == 0:
            return None
        _, d_value, c_value, _, _ = self.shifted_coefficients
        q_value = self.base_ordinate
        u_value = (
            4 * q_value**2 * (x_value + c_value) - d_value**2
        ) / (2 * q_value * y_value)
        if u_value == 0:
            return None
        ordinate = (x_value * u_value**2 - d_value * u_value) / (2 * q_value) - q_value
        answer = self.base_parameter + u_value, ordinate
        if ordinate**2 != self.quartic_value(answer[0]):
            raise AssertionError("the pointed-quartic inverse map failed")
        return answer


def exact_radical_proxy(
    parameter: Fraction, primes: Sequence[int]
) -> dict[str, Any]:
    parameter = abs(Q(parameter))
    numerator, denominator = parameter.numerator, parameter.denominator
    degree = len(DISCRIMINANT_POLYNOMIAL) - 1
    numerator_powers = [1]
    denominator_powers = [1]
    for _ in range(degree):
        numerator_powers.append(numerator_powers[-1] * numerator)
        denominator_powers.append(denominator_powers[-1] * denominator)
    discriminant = abs(
        sum(
            coefficient
            * numerator_powers[power]
            * denominator_powers[degree - power]
            for power, coefficient in enumerate(DISCRIMINANT_POLYNOMIAL)
        )
    )
    if discriminant == 0:
        raise ValueError("the orbit reached a singular specialization")
    remaining = discriminant
    savings = 0.0
    valuations = []
    for prime in primes:
        valuation = 0
        while remaining % prime == 0:
            remaining //= prime
            valuation += 1
        if valuation:
            valuations.append([prime, valuation])
            savings += (valuation - 1) * log(prime)
    raw_log = log(discriminant)
    return {
        "trial_prime_bound": max(primes),
        "raw_log_absolute_homogenized_discriminant": raw_log,
        "known_repeated_prime_log_savings": savings,
        "log_radical_upper_proxy": raw_log - savings,
        "small_prime_valuations": valuations,
        "unfactored_cofactor_decimal_digits_upper_bound": (
            (remaining.bit_length() * 30103) // 100000 + 1
        ),
    }


def enumerate_orbit(
    auxiliary: PointedQuartic,
    basis: Sequence[tuple[Fraction, Fraction]],
) -> tuple[dict[Fraction, dict[str, Any]], dict[str, Any]]:
    coefficient_choices = tuple(
        range(
            -ORBIT_MAXIMUM_ABSOLUTE_COEFFICIENT,
            ORBIT_MAXIMUM_ABSOLUTE_COEFFICIENT + 1,
        )
    )
    multiples = tuple(
        tuple(
            weierstrass_multiply(auxiliary.weierstrass_coefficients, point, scalar)
            for scalar in coefficient_choices
        )
        for point in basis
    )
    by_parameter: dict[Fraction, dict[str, Any]] = {}
    vector_count = 0

    def recurse(
        index: int,
        point: tuple[Fraction, Fraction] | None,
        l1_norm: int,
        vector: tuple[int, ...],
    ) -> None:
        nonlocal vector_count
        if index == len(basis):
            if l1_norm == 0:
                return
            vector_count += 1
            inverse = auxiliary.inverse(point)
            if inverse is None or inverse[0] == 0:
                return
            signed_parameter, normalized_ordinate = inverse
            canonical_parameter = abs(signed_parameter)
            raw_x = PRODUCTIVE_SLOPE * signed_parameter + PRODUCTIVE_INTERCEPT
            previous = by_parameter.get(canonical_parameter)
            candidate = {
                "signed_parameter": signed_parameter,
                "normalized_ordinate": normalized_ordinate,
                "forced_quartic_x_at_canonical_curve": raw_x,
                "coefficient_vector": vector,
            }
            if previous is None or vector < previous["coefficient_vector"]:
                by_parameter[canonical_parameter] = candidate
            return
        for scalar in coefficient_choices:
            new_l1 = l1_norm + abs(scalar)
            if new_l1 > ORBIT_MAXIMUM_L1_NORM:
                continue
            next_point = weierstrass_add(
                auxiliary.weierstrass_coefficients,
                point,
                multiples[index][scalar + ORBIT_MAXIMUM_ABSOLUTE_COEFFICIENT],
            )
            recurse(index + 1, next_point, new_l1, vector + (scalar,))

    recurse(0, None, 0, ())
    ordered = tuple(
        sorted(
            (
                parameter,
                tuple(record["coefficient_vector"]),
            )
            for parameter, record in by_parameter.items()
        )
    )
    return by_parameter, {
        "enumerated_nonzero_coefficient_vector_count": vector_count,
        "distinct_nonzero_canonical_parameter_count": len(by_parameter),
        "parameter_sha256": polynomial_digest(parameter for parameter, _ in ordered),
        "parameter_and_vector_sha256": orbit_record_digest(ordered),
    }


def conductor_attempt(
    parameter: Fraction, *, timeout: float, stack_bytes: int
) -> dict[str, Any]:
    try:
        data = minimal_curve_data(
            short_jacobian_coefficients(RANK21_CONSTRUCTION, parameter),
            timeout=timeout,
            stack_bytes=stack_bytes,
            local_primes=(2, 3, 5, 7, 11, 13, 17, 19, 23),
        )
    except (subprocess.TimeoutExpired, RuntimeError, ValueError) as error:
        return {
            "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
            "error": str(error)[:500],
            "retried": False,
        }
    return {
        "status": "completed",
        **data,
        "below_strict_log_conductor_target": (
            Decimal(data["log_conductor"]) < TARGET_LOG_CONDUCTOR
        ),
        "retried": False,
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slice-height", type=int, default=200_000)
    parser.add_argument("--slice-timeout", type=float, default=20.0)
    parser.add_argument("--height-timeout", type=float, default=30.0)
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--conductor-timeout", type=float, default=20.0)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--output", type=Path, default=root / DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.slice_height != 200_000:
        raise SystemExit("the canonical replay pins slice height 200000")
    if min(
        args.slice_timeout,
        args.height_timeout,
        args.saturation_timeout,
        args.conductor_timeout,
    ) <= 0 or max(
        args.slice_timeout,
        args.height_timeout,
        args.saturation_timeout,
        args.conductor_timeout,
    ) > 60:
        raise SystemExit("all one-shot timeouts must lie in (0,60]")
    if args.stack_bytes < 8_000_000:
        raise SystemExit("the PARI stack cap is too small")

    started = time.monotonic()
    validate_bivariate_quartic()
    short_points = tuple(minimum_to_short(point) for point in RANK21_PUBLISHED_POINTS)
    if tuple(short_to_minimum(point) for point in short_points) != RANK21_PUBLISHED_POINTS:
        raise AssertionError("the exact published model round trip failed")
    reconstruction, convention_trials = select_reconstruction_convention()
    visible = primitive_visible_points(RANK21_CONSTRUCTION, T0)
    visible_x = {point[0] for point in visible}
    accidental = tuple(point for point in reconstruction if point[0] not in visible_x)
    if len(accidental) != 11:
        raise AssertionError("the published accidental-preimage count changed")
    slices = build_slices(accidental)
    productive = next(
        item
        for item in slices
        if item.source_point[0] == PRODUCTIVE_SOURCE_X
        and item.slope == PRODUCTIVE_SLOPE
    )
    if productive.intercept != PRODUCTIVE_INTERCEPT:
        raise AssertionError("the productive slice intercept changed")

    run_records = []
    normalized_points_by_slice: dict[str, tuple[tuple[Fraction, Fraction], ...]] = {}
    forced_by_parameter: dict[Fraction, list[dict[str, Any]]] = {}
    for index, item in enumerate(slices, 1):
        polynomial = tuple(Q(value) for value in item.normalized.normalized_coefficients)
        raw, process_record = search_original_quartic(
            polynomial,
            str(args.slice_height),
            timeout=args.slice_timeout,
            stack_bytes=args.stack_bytes,
        )
        if process_record["status"] != "completed":
            raise RuntimeError(f"slice {item.identifier} did not complete")
        exact_points = tuple(sorted(set(raw)))
        normalized_points_by_slice[item.identifier] = exact_points
        new_count = 0
        for point in signless_quartic_points(exact_points):
            signed_parameter, raw_ordinate = normalized_to_raw_point(item, point)
            if abs(signed_parameter) == T0:
                continue
            canonical_parameter = abs(signed_parameter)
            raw_x = item.raw_x(signed_parameter)
            quartic = primitive_quartic_coefficients(
                RANK21_CONSTRUCTION, canonical_parameter
            )
            if raw_ordinate**2 != quartic_value(quartic, raw_x):
                raise AssertionError("a canonicalized slice point missed the family")
            if raw_x in visible_abscissas(canonical_parameter):
                continue
            forced_by_parameter.setdefault(canonical_parameter, []).append(
                {
                    "slice_id": item.identifier,
                    "signed_parameter_before_T_symmetry": rational_record(signed_parameter),
                    "forced_quartic_x": rational_record(raw_x),
                }
            )
            new_count += 1
        run_records.append(
            {
                "slice_id": item.identifier,
                "accidental_index_zero_based": item.accidental_index,
                "source_x": rational_record(item.source_point[0]),
                "slope": item.slope,
                "intercept": rational_record(item.intercept),
                "normalized_coefficients_ascending": list(
                    item.normalized.normalized_coefficients
                ),
                **process_record,
                "exact_forced_nonvisible_signless_incidence_count": new_count,
            }
        )
        print(
            f"slice {index}/{len(slices)} {item.identifier} "
            f"points={len(exact_points)} forced={new_count}",
            flush=True,
        )

    auxiliary = PointedQuartic.from_slice(productive)
    productive_points = normalized_points_by_slice[productive.identifier]
    auxiliary_points = tuple(
        sorted(
            {
                auxiliary.forward(point)
                for point in productive_points
                if point[0] != T0
            }
        )
    )
    height_runs = height_matrix_replay(
        auxiliary.weierstrass_coefficients,
        auxiliary_points,
        precisions=(72, 120),
        timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    auxiliary_rank = stable_height_rank(height_runs)
    if auxiliary_rank != 8:
        raise AssertionError("the productive auxiliary numerical rank changed")
    indices = tuple(height_runs[-1]["subset_indices_one_based"])
    selected = tuple(auxiliary_points[index - 1] for index in indices)
    saturated, saturation = saturate_exact_basis(
        auxiliary.weierstrass_coefficients,
        selected,
        prime_bound=20,
        timeout=args.saturation_timeout,
        stack_bytes=args.stack_bytes,
    )
    if len(saturated) != 8:
        raise AssertionError("small-prime saturation changed the auxiliary dimension")

    orbit, orbit_summary = enumerate_orbit(auxiliary, saturated)
    print(
        f"orbit vectors={orbit_summary['enumerated_nonzero_coefficient_vector_count']} "
        f"parameters={orbit_summary['distinct_nonzero_canonical_parameter_count']}",
        flush=True,
    )
    primes = primes_up_to(PROXY_PRIME_BOUND)
    proxy_records = []
    plausible_orbit_parameters: list[Fraction] = []
    for index, (parameter, record) in enumerate(sorted(orbit.items()), 1):
        try:
            proxy = exact_radical_proxy(parameter, primes)
        except ValueError:
            continue
        raw_x = Q(record["forced_quartic_x_at_canonical_curve"])
        generic = raw_x in visible_abscissas(parameter)
        record["proxy"] = proxy
        record["generic_visible_intersection"] = generic
        if (
            not generic
            and proxy["log_radical_upper_proxy"] < PROXY_EXACT_CONDUCTOR_GATE
        ):
            plausible_orbit_parameters.append(parameter)
        proxy_records.append(
            {
                "parameter": parameter,
                "coefficient_vector": tuple(record["coefficient_vector"]),
                "forced_quartic_x": raw_x,
                "generic_visible_intersection": generic,
                **proxy,
            }
        )
        if index % 10_000 == 0:
            print(
                f"orbit proxy {index}/{len(orbit)} "
                f"plausible_nonvisible={len(plausible_orbit_parameters)}",
                flush=True,
            )
    proxy_records.sort(
        key=lambda record: (
            record["log_radical_upper_proxy"],
            max(abs(record["parameter"].numerator), record["parameter"].denominator),
            record["parameter"],
        )
    )

    exact_conductor_parameters = sorted(
        set(forced_by_parameter) | set(plausible_orbit_parameters)
    )
    conductors = {}
    for parameter in exact_conductor_parameters:
        conductors[rational_record(parameter)] = conductor_attempt(
            parameter,
            timeout=args.conductor_timeout,
            stack_bytes=args.stack_bytes,
        )
        print(
            f"conductor T={rational_record(parameter)} "
            f"status={conductors[rational_record(parameter)]['status']}",
            flush=True,
        )
    hits = [
        parameter
        for parameter in exact_conductor_parameters
        if conductors[rational_record(parameter)].get(
            "below_strict_log_conductor_target", False
        )
    ]

    top_proxy_records = proxy_records[:64]
    artifact = {
        "schema_version": 1,
        "status": "bounded_auxiliary_group_orbit_complete",
        "goal": (
            "find E/Q with rank>=21 and log conductor<182.72, or rank>=30"
        ),
        "target_hit": False,
        "rank21_subthreshold_parameters_requiring_rank_triage": [
            rational_record(parameter) for parameter in hits
        ],
        "source": {
            "primary_source": PRIMARY_SOURCE,
            "root_tuple": [int(root) for root in RANK21_CONSTRUCTION.roots],
            "constructor_parameter": rational_record(T0),
            "published_parameter": "14721/376",
            "published_point_count": len(RANK21_PUBLISHED_POINTS),
            "model_change_short_to_published_u_r_s_t": [
                rational_record(value) for value in MODEL_CHANGE_SHORT_TO_MINIMAL
            ],
            "all_published_points_exactly_round_tripped": True,
        },
        "published_point_reconstruction": {
            "selected_base_index_one_based": 1,
            "selected_sign": 1,
            "exact_oriented_generic_matches": 6,
            "generic_abscissa_matches": 10,
            "accidental_preimage_count": len(accidental),
            "all_convention_trials": list(convention_trials),
            "accidental_preimages": [
                {
                    "index_zero_based": index,
                    "x": rational_record(point[0]),
                    "y": rational_record(point[1]),
                }
                for index, point in enumerate(accidental)
            ],
            "accidental_preimage_sha256": hashlib.sha256(
                "\n".join(
                    f"{rational_record(point[0])},{rational_record(point[1])}"
                    for point in accidental
                ).encode()
            ).hexdigest(),
        },
        "slice_search": {
            "height_bound": args.slice_height,
            "timeout_seconds_per_slice": args.slice_timeout,
            "slice_count": len(slices),
            "completed_slice_count": sum(
                record["status"] == "completed" for record in run_records
            ),
            "run_records": run_records,
            "distinct_forced_nonvisible_parameter_count": len(forced_by_parameter),
            "forced_nonvisible_parameters": [
                {
                    "constructor_parameter": rational_record(parameter),
                    "incidences": incidences,
                }
                for parameter, incidences in sorted(forced_by_parameter.items())
            ],
        },
        "productive_auxiliary_curve": {
            "slice_id": productive.identifier,
            "source_x": rational_record(productive.source_point[0]),
            "source_y": rational_record(productive.source_point[1]),
            "slope": productive.slope,
            "intercept": rational_record(productive.intercept),
            "normalized_quartic_coefficients_ascending": list(
                productive.normalized.normalized_coefficients
            ),
            "pointed_base_ordinate": rational_record(auxiliary.base_ordinate),
            "generalized_weierstrass_coefficients": [
                rational_record(value) for value in auxiliary.weierstrass_coefficients
            ],
            "H200000_signed_point_count": len(productive_points),
            "mapped_nonexceptional_auxiliary_point_count": len(auxiliary_points),
            "height_selection": {
                "stable_numerical_rank": auxiliary_rank,
                "precision_runs": list(height_runs),
                "selected_indices_one_based": list(indices),
                "selected_point_sha256": point_digest(selected),
                "scope_warning": "numerical generator selection only",
            },
            "small_prime_saturation": saturation,
            "saturated_basis": [rational_pair_record(point) for point in saturated],
        },
        "orbit": {
            "coefficient_range": [
                -ORBIT_MAXIMUM_ABSOLUTE_COEFFICIENT,
                ORBIT_MAXIMUM_ABSOLUTE_COEFFICIENT,
            ],
            "maximum_l1_norm": ORBIT_MAXIMUM_L1_NORM,
            **orbit_summary,
            "proxy_prime_bound": PROXY_PRIME_BOUND,
            "exact_conductor_gate": PROXY_EXACT_CONDUCTOR_GATE,
            "nonvisible_proxy_below_gate_count": len(plausible_orbit_parameters),
            "nonvisible_proxy_below_gate_parameters": [
                rational_record(parameter) for parameter in plausible_orbit_parameters
            ],
            "top_64_proxy_records": [
                {
                    **{
                        key: value
                        for key, value in record.items()
                        if key not in {"parameter", "forced_quartic_x"}
                    },
                    "parameter": rational_record(record["parameter"]),
                    "forced_quartic_x": rational_record(record["forced_quartic_x"]),
                    "coefficient_vector": list(record["coefficient_vector"]),
                }
                for record in top_proxy_records
            ],
        },
        "exact_conductors": conductors,
        "conclusion": {
            "target_hit": False,
            "strict_subtarget_parameter_count": len(hits),
            "statement": (
                "The complete declared slice and coefficient-box orbit produced "
                "no certified rank-21/subthreshold hit."
            ),
            "scope_warning": (
                "This is exhaustive only for the 22 H=200000 slice boxes and "
                "the declared rank-8 coefficient box; it is not an upper-rank proof."
            ),
        },
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "pari": pari_version(),
            "script_sha256": file_sha256(Path(__file__).resolve()),
            "wall_seconds": time.monotonic() - started,
            "all_subprocesses_foreground": True,
            "one_attempt_no_retries": True,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
