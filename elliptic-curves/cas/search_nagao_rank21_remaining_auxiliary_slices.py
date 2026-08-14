#!/usr/bin/env python3
"""Exhaust the 21 nonproductive auxiliary slices through Nagao's rank-21 fibre.

The published fibre for roots (399,380,352,47,4,0) occurs at constructor
parameter T0=14721/188.  This program reconstructs all 21 printed-point
preimages under the binary-quartic covariant, identifies the eleven
accidental preimages, and forms both degree-cancelling slopes through each.
The separately owned productive slice x=-T+57361/139 is explicitly excluded.

Every other slice receives exactly one hyperellratpoints run at H=200000.
Using its source point as origin, the normalized quartic is converted exactly
to a pointed generalized Weierstrass model.  The complete bounded point set is
mapped to that model, minimized inside PARI, and a basis is selected only when
72- and 120-digit height computations agree.  All nonzero ternary coefficient
vectors are then exhausted.  A quinary expansion is declared in advance but
is enabled only for a rank-at-least-four slice whose ternary search produces
at least 64 accepted incidences and a radical proxy below 250.

Generated specializations are checked exactly against the original quartic,
the complete bounded panel, singular fibres, visible sections, and visible
covariant sign-pairs.  Every proxy below 190 receives an exact conductor; the
two best remaining proxy leaders are also computed as a bounded calibration.
This is a finite experiment.  Numerical height ranks select generators but do
not prove the exact Mordell--Weil ranks of the auxiliary curves.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import itertools
import json
from math import log
from pathlib import Path
import platform
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

import sympy as sp

from alternate_quartic_covers import short_add
from ek_k3 import primes_up_to, rational_square_root, rational_to_string
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
from search_extra_points import gp_rational, gp_vector, run_gp
from search_nagao_rank21_accidental_slices import NormalizedSlice, normalize_slice
from search_nagao_rank21_t956_skew import search_original_quartic
from search_nagao_section7_auxiliary_jacobians import (
    translate_polynomial,
    weierstrass_add,
    weierstrass_multiply,
)
from triage_nagao_rank13_finalists import (
    parse_vecsmall,
    point_digest,
    point_on_short_curve,
    stable_height_rank,
)


Q = Fraction
if hasattr(sys, "set_int_max_str_digits"):
    # Exact auxiliary group sums routinely exceed Python's conservative
    # decimal-conversion guard; the search remains bounded by its vector box.
    sys.set_int_max_str_digits(0)
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
SLICE_HEIGHT = 200_000
HEIGHT_PRECISIONS = (72, 120)
TERNARY_ALPHABET = (-1, 0, 1)
QUINARY_ALPHABET = (-2, -1, 0, 1, 2)
EXPANSION_MINIMUM_RANK = 4
EXPANSION_MINIMUM_ACCEPTED_INCIDENCES = 64
EXPANSION_PROXY_GATE = 250.0
PROXY_PRIME_BOUND = 2_000
PROXY_EXACT_CONDUCTOR_GATE = 190.0
LEADER_EXACT_CONDUCTOR_COUNT = 2
EXPECTED_PREIMAGE_SHA256 = (
    "53ea5b29844b37bcb27ad47d9b2afe68c83bd3a44d1c8295e6d5d1aa8465d862"
)
EXPECTED_ACCIDENTAL_SHA256 = (
    "f5a019f38f3f8bfc99be4ff80356025bee4c2df58ea7a784746c0475cd89135b"
)
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/"
    "elliptic_nagao_rank21_remaining_auxiliary_slices.json"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank21_remaining_auxiliary_slices.py"
)


# Ascending T coefficients for the ascending x coefficients (e,d,c,b,a).
# They are evaluated against the generic constructor at seven independent T's.
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


def point_record(point: tuple[Fraction, Fraction]) -> dict[str, str]:
    return {"x": rational_record(point[0]), "y": rational_record(point[1])}


def polynomial_value(
    coefficients: Sequence[Fraction], value: Fraction
) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * Q(value) + Q(coefficient)
    return answer


def rational_digest(values: Iterable[Fraction]) -> str:
    payload = "\n".join(rational_record(value) for value in values)
    return hashlib.sha256(payload.encode()).hexdigest()


def point_sequence_digest(points: Iterable[tuple[Fraction, Fraction]]) -> str:
    payload = "\n".join(
        f"{rational_record(point[0])},{rational_record(point[1])}"
        for point in points
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def validate_bivariate_quartic() -> None:
    for parameter in (Q(-9), Q(-1), Q(1, 2), Q(1), Q(3, 2), Q(8), T0):
        replay = tuple(
            polynomial_value(coefficients, parameter)
            for coefficients in RANK21_QUARTIC_COEFFICIENT_POLYNOMIALS
        )
        if replay != primitive_quartic_coefficients(RANK21_CONSTRUCTION, parameter):
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


def minimum_to_short(
    point: tuple[Fraction, Fraction]
) -> tuple[Fraction, Fraction]:
    u_value, r_value, s_value, t_value = MODEL_CHANGE_SHORT_TO_MINIMAL
    x_value, y_value = (Q(value) for value in point)
    answer = (
        u_value**2 * x_value + r_value,
        u_value**3 * y_value + s_value * u_value**2 * x_value + t_value,
    )
    if not point_on_short_curve(rank21_short_jacobian_coefficients(), answer):
        raise AssertionError("the published-to-short model change failed")
    return answer


def short_to_minimum(
    point: tuple[Fraction, Fraction]
) -> tuple[Fraction, Fraction]:
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


def rational_linear_roots(
    expression: sp.Expr, symbol: sp.Symbol
) -> tuple[Fraction, ...]:
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
    target: tuple[Fraction, Fraction]
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
        raise ValueError(
            f"expected one rational covariant preimage, found {len(answers)}"
        )
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
        record
        for record in trials
        if int(record["generic_abscissa_matches"]) == maximum
    ]
    if len(winners) != 1:
        raise AssertionError("the published reconstruction convention is not unique")
    winner = winners[0]
    key = int(winner["base_index_one_based"]) - 1, int(winner["sign"])
    preimages = successful[key]
    if key != (0, 1) or maximum != 10:
        raise AssertionError("the pinned reconstruction convention changed")
    if int(winner["exact_oriented_generic_matches"]) != 6:
        raise AssertionError("the oriented generic match count changed")
    if point_sequence_digest(preimages) != EXPECTED_PREIMAGE_SHA256:
        raise AssertionError("the printed-point preimage digest changed")
    return preimages, tuple(trials)


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
    accidental: Sequence[tuple[Fraction, Fraction]]
) -> tuple[Slice, ...]:
    answer = []
    for accidental_index, source_point in enumerate(accidental):
        for slope in (-1, 1):
            intercept = source_point[0] - slope * T0
            normalized = normalize_slice(slice_polynomial(slope, intercept))
            if normalized.normalized_degree != 4 or normalized.genus != 1:
                raise AssertionError("an accidental +/-1 slice ceased to be genus one")
            if polynomial_value(normalized.raw_coefficients, T0) != source_point[1] ** 2:
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
    raw_ordinate = item.normalized.original_ordinate(
        parameter, normalized_ordinate
    )
    if raw_ordinate**2 != polynomial_value(
        item.normalized.raw_coefficients, parameter
    ):
        raise AssertionError("a normalized slice point failed raw replay")
    return parameter, raw_ordinate


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
            raise AssertionError("the source lies at a removed-square zero")
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

    def conjugate_base_image(self) -> tuple[Fraction, Fraction]:
        _, d_value, c_value, b_value, _ = self.shifted_coefficients
        q_value = self.base_ordinate
        discriminant_term = d_value**2 - 4 * q_value**2 * c_value
        answer = (
            discriminant_term / (4 * q_value**2),
            -d_value * discriminant_term / (4 * q_value**3) - 2 * q_value * b_value,
        )
        if not point_on_short_curve(self.weierstrass_coefficients, answer):
            raise AssertionError("the conjugate-base limit missed the auxiliary curve")
        return answer

    def forward(
        self, point: tuple[Fraction, Fraction]
    ) -> tuple[Fraction, Fraction] | None:
        parameter, ordinate = (Q(value) for value in point)
        if ordinate**2 != self.quartic_value(parameter):
            raise ValueError("the point missed the auxiliary quartic")
        if parameter == self.base_parameter:
            if ordinate == self.base_ordinate:
                return None
            if ordinate == -self.base_ordinate:
                return self.conjugate_base_image()
            raise AssertionError("the base fibre had an unexpected ordinate")
        u_value = parameter - self.base_parameter
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
        ordinate = (
            (x_value * u_value**2 - d_value * u_value) / (2 * q_value)
            - q_value
        )
        answer = self.base_parameter + u_value, ordinate
        if ordinate**2 != self.quartic_value(answer[0]):
            raise AssertionError("the pointed-quartic inverse map failed")
        return answer


def minimized_height_matrix_replay(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    precisions: tuple[int, ...],
    timeout: float,
    stack_bytes: int,
) -> tuple[dict[str, Any], ...]:
    """Replay heights after an exact PARI minimal-model point change.

    Each precision is a separate foreground process.  This is important for
    the second accidental preimage, whose exact coordinates have hundreds of
    digits but whose individual 72/120-digit runs remain below the cap.
    """

    if not points:
        raise ValueError("height replay needs at least one point")
    if any(not point_on_short_curve(coefficients, point) for point in points):
        raise AssertionError("a height-replay point missed its exact curve")
    curve = ",".join(gp_rational(Q(value)) for value in coefficients)
    point_vector = ",".join(gp_vector(point) for point in points)
    records = []
    for precision in precisions:
        program = "\n".join(
            (
                f"default(realprecision,{precision});",
                f"E0=ellinit([{curve}]);",
                "E=ellminimalmodel(E0,&V);",
                f"P0=[{point_vector}];",
                "P=vector(#P0,i,ellchangepoint(P0[i],V));",
                'print("MODEL_CHANGE ",V);',
                'print("MINIMAL_MODEL ",[E.a1,E.a2,E.a3,E.a4,E.a6]);',
                'print("ON_CURVE ",vecsum(vector(#P,i,ellisoncurve(E,P[i]))));',
                "H=ellheightmatrix(E,P);",
                "IX=matindexrank(H);",
                "K=vecextract(P,IX[2]);",
                "HK=ellheightmatrix(E,K);",
                'print("HEIGHT_BEGIN");',
                "print(matrank(H));",
                "print(IX[2]);",
                "print(matdet(HK));",
                "EV=mateigen(HK,1)[1];print(vecmin(EV));print(vecmax(EV));",
                'print("HEIGHT_END");',
                "quit",
            )
        ) + "\n"
        output, wall_seconds = run_gp(
            program, timeout=timeout, stack_bytes=stack_bytes
        )
        lines = [line.strip() for line in output.splitlines() if line.strip()]
        if f"ON_CURVE {len(points)}" not in lines:
            raise AssertionError("the minimal-model point change failed")
        start = lines.index("HEIGHT_BEGIN") + 1
        end = lines.index("HEIGHT_END")
        values = lines[start:end]
        records.append(
            {
                "decimal_precision": precision,
                "numerical_rank": int(values[0]),
                "subset_indices_one_based": parse_vecsmall(values[1]),
                "subset_height_determinant": values[2],
                "subset_smallest_eigenvalue": values[3],
                "subset_largest_eigenvalue": values[4],
                "minimal_model_change": next(
                    line.removeprefix("MODEL_CHANGE ")
                    for line in lines
                    if line.startswith("MODEL_CHANGE ")
                ),
                "minimal_model": next(
                    line.removeprefix("MINIMAL_MODEL ")
                    for line in lines
                    if line.startswith("MINIMAL_MODEL ")
                ),
                "gp_process_wall_seconds": wall_seconds,
            }
        )
    return tuple(records)


def coefficient_vectors(
    dimension: int, alphabet: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
    if dimension < 0 or dimension > 6:
        raise ValueError("the declared auxiliary dimension cap is six")
    choices = tuple(int(value) for value in alphabet)
    if 0 not in choices:
        raise ValueError("the coefficient alphabet must contain zero")
    return tuple(
        vector
        for vector in itertools.product(choices, repeat=dimension)
        if any(vector)
    )


def recursive_group_elements(
    coefficients: Sequence[Fraction],
    basis: Sequence[tuple[Fraction, Fraction]],
    alphabet: Sequence[int],
) -> tuple[tuple[tuple[int, ...], tuple[Fraction, Fraction] | None], ...]:
    states: list[tuple[tuple[int, ...], tuple[Fraction, Fraction] | None]] = [
        ((), None)
    ]
    choices = tuple(int(value) for value in alphabet)
    for basis_point in basis:
        multiples = {
            scalar: weierstrass_multiply(coefficients, basis_point, scalar)
            for scalar in choices
        }
        next_states = []
        for vector, point in states:
            for scalar in choices:
                next_states.append(
                    (
                        vector + (scalar,),
                        weierstrass_add(coefficients, point, multiples[scalar]),
                    )
                )
        states = next_states
    answer = tuple((vector, point) for vector, point in states if any(vector))
    if len(answer) != len(coefficient_vectors(len(basis), choices)):
        raise AssertionError("the recursive coefficient population changed")
    return answer


def homogenized_discriminant(parameter: Fraction) -> int:
    parameter = abs(Q(parameter))
    coefficients = RANK21_CONSTRUCTION.primitive_discriminant_polynomial
    if any(value.denominator != 1 for value in coefficients) or len(coefficients) != 21:
        raise AssertionError("the rank-21 discriminant normalization changed")
    numerator, denominator = parameter.numerator, parameter.denominator
    degree = len(coefficients) - 1
    value = sum(
        int(coefficient)
        * numerator**power
        * denominator ** (degree - power)
        for power, coefficient in enumerate(coefficients)
    )
    if value == 0:
        raise ValueError("singular rank-21 specialization")
    return value


def exact_radical_proxy(
    parameter: Fraction, primes: Sequence[int]
) -> dict[str, Any]:
    discriminant = abs(homogenized_discriminant(parameter))
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
            remaining.bit_length() * 30103 // 100000 + 1
        ),
    }


def canonical_visible_data(
    parameter: Fraction,
) -> tuple[set[Fraction], set[Fraction]]:
    points = primitive_visible_points(RANK21_CONSTRUCTION, Q(parameter))
    return (
        {point[0] for point in points},
        {
            quartic_point_to_short_jacobian(
                RANK21_CONSTRUCTION, Q(parameter), point
            )[0]
            for point in points
        },
    )


def expansion_is_warranted(
    *, stable_rank: int, accepted_incidences: int, best_proxy: float | None
) -> bool:
    return (
        stable_rank >= EXPANSION_MINIMUM_RANK
        and accepted_incidences >= EXPANSION_MINIMUM_ACCEPTED_INCIDENCES
        and best_proxy is not None
        and best_proxy < EXPANSION_PROXY_GATE
    )


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
            "status": (
                "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error"
            ),
            "error": str(error)[:500],
        }
    return {
        "status": "completed",
        **data,
        "below_strict_log_conductor_target": (
            Decimal(data["log_conductor"]) < TARGET_LOG_CONDUCTOR
        ),
    }


def write_artifact(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slice-height", type=int, default=SLICE_HEIGHT)
    parser.add_argument("--slice-timeout", type=float, default=20.0)
    parser.add_argument("--height-timeout", type=float, default=60.0)
    parser.add_argument("--conductor-timeout", type=float, default=30.0)
    parser.add_argument("--stack-bytes", type=int, default=768_000_000)
    parser.add_argument("--output", type=Path, default=root / DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.slice_height != SLICE_HEIGHT:
        raise SystemExit("the canonical replay pins slice height 200000")
    if min(args.slice_timeout, args.height_timeout, args.conductor_timeout) <= 0:
        raise SystemExit("timeouts must be positive")
    if max(args.slice_timeout, args.height_timeout, args.conductor_timeout) > 60:
        raise SystemExit("every one-shot timeout must be at most 60 seconds")
    if args.stack_bytes < 8_000_000:
        raise SystemExit("the PARI stack cap is too small")

    started = time.monotonic()
    validate_bivariate_quartic()
    short_points = tuple(minimum_to_short(point) for point in RANK21_PUBLISHED_POINTS)
    if tuple(short_to_minimum(point) for point in short_points) != RANK21_PUBLISHED_POINTS:
        raise AssertionError("the published-model exact round trip failed")
    reconstruction, convention_trials = select_reconstruction_convention()
    visible = primitive_visible_points(RANK21_CONSTRUCTION, T0)
    visible_x = {point[0] for point in visible}
    accidental = tuple(point for point in reconstruction if point[0] not in visible_x)
    if len(accidental) != 11 or point_sequence_digest(accidental) != EXPECTED_ACCIDENTAL_SHA256:
        raise AssertionError("the eleven accidental preimages changed")
    all_slices = build_slices(accidental)
    productive = tuple(
        item
        for item in all_slices
        if item.source_point[0] == PRODUCTIVE_SOURCE_X
        and item.slope == PRODUCTIVE_SLOPE
    )
    if len(productive) != 1 or productive[0].intercept != PRODUCTIVE_INTERCEPT:
        raise AssertionError("the separately owned productive slice changed")
    slices = tuple(item for item in all_slices if item != productive[0])
    if len(slices) != 21:
        raise AssertionError("the nonproductive slice panel changed")

    slice_records: list[dict[str, Any]] = []
    slice_runtime: dict[str, dict[str, Any]] = {}
    prior_parameters = {abs(T0)}
    for index, item in enumerate(slices, 1):
        polynomial = tuple(Q(value) for value in item.normalized.normalized_coefficients)
        raw_points, process = search_original_quartic(
            polynomial,
            str(args.slice_height),
            timeout=args.slice_timeout,
            stack_bytes=args.stack_bytes,
        )
        if process["status"] != "completed":
            raise RuntimeError(f"slice {item.identifier} did not complete")
        normalized_points = tuple(sorted(set(raw_points)))
        if (T0, item.source_point[1] / (
            item.normalized.removed_square_value(T0)
            * item.normalized.ordinate_constant_scale
        )) not in normalized_points:
            raise AssertionError("the bounded run omitted its defining source point")
        prior_parameters.update(abs(point[0]) for point in normalized_points)

        auxiliary = PointedQuartic.from_slice(item)
        mapped = []
        mapped_sources: dict[tuple[Fraction, Fraction], list[dict[str, str]]] = {}
        for quartic_point in normalized_points:
            auxiliary_point = auxiliary.forward(quartic_point)
            if auxiliary_point is None:
                continue
            mapped_sources.setdefault(auxiliary_point, []).append(
                {
                    "T": rational_record(quartic_point[0]),
                    "normalized_ordinate": rational_record(quartic_point[1]),
                }
            )
            if auxiliary_point not in mapped:
                mapped.append(auxiliary_point)
        mapped_points = tuple(mapped)
        height_runs = minimized_height_matrix_replay(
            auxiliary.weierstrass_coefficients,
            mapped_points,
            precisions=HEIGHT_PRECISIONS,
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
        )
        stable_rank = stable_height_rank(height_runs)
        indices = tuple(height_runs[-1]["subset_indices_one_based"])
        if stable_rank != len(indices) or stable_rank > 6:
            raise AssertionError("the stable auxiliary basis dimension changed")
        basis = tuple(mapped_points[position - 1] for position in indices)
        slice_runtime[item.identifier] = {
            "item": item,
            "auxiliary": auxiliary,
            "normalized_points": normalized_points,
            "basis": basis,
            "stable_rank": stable_rank,
        }
        record = {
            "id": item.identifier,
            "accidental_index_one_based": item.accidental_index + 1,
            "source_point": point_record(item.source_point),
            "slope": item.slope,
            "intercept": rational_record(item.intercept),
            "normalized_quartic_coefficients_ascending": list(
                item.normalized.normalized_coefficients
            ),
            "H200000_process": process,
            "H200000_signed_point_count": len(normalized_points),
            "H200000_distinct_parameter_count": len(
                {point[0] for point in normalized_points}
            ),
            "H200000_parameter_sha256": rational_digest(
                sorted({point[0] for point in normalized_points})
            ),
            "generalized_weierstrass_coefficients": [
                rational_record(value)
                for value in auxiliary.weierstrass_coefficients
            ],
            "pointed_base_ordinate": rational_record(auxiliary.base_ordinate),
            "conjugate_base_image": point_record(auxiliary.conjugate_base_image()),
            "mapped_nonidentity_point_count": len(mapped_points),
            "height_selection": {
                "stable_numerical_rank": stable_rank,
                "precision_runs": list(height_runs),
                "selected_indices_one_based": list(indices),
                "selected_point_sha256": point_digest(basis),
                "basis": [point_record(point) for point in basis],
                "scope_warning": "numerical generator selection only",
            },
        }
        slice_records.append(record)
        print(
            f"slice {index}/{len(slices)} {item.identifier} "
            f"points={len(normalized_points)} rank={stable_rank}",
            flush=True,
        )

    checkpoint = {
        "schema_version": 1,
        "status": "exact_H200000_auxiliary_checkpoint",
        "target_hit": False,
        "scope": {
            "root_tuple": [int(root) for root in RANK21_CONSTRUCTION.roots],
            "constructor_parameter": rational_record(T0),
            "excluded_productive_slice": {
                "source_x": rational_record(PRODUCTIVE_SOURCE_X),
                "slope": PRODUCTIVE_SLOPE,
                "intercept": rational_record(PRODUCTIVE_INTERCEPT),
            },
            "searched_slice_ids": [item.identifier for item in slices],
            "slice_height": SLICE_HEIGHT,
        },
        "published_reconstruction": {
            "preimage_sha256": EXPECTED_PREIMAGE_SHA256,
            "accidental_sha256": EXPECTED_ACCIDENTAL_SHA256,
            "all_convention_trials": list(convention_trials),
        },
        "slices": slice_records,
        "checkpoint": {
            "complete_H200000_slice_count": len(slice_records),
            "written_before_orbit_enumeration": True,
        },
    }
    write_artifact(args.output, checkpoint)
    print(f"checkpoint written {args.output}", flush=True)

    candidates: dict[Fraction, list[dict[str, Any]]] = {}
    slice_generation: dict[str, dict[str, Any]] = {}

    def enumerate_slice(item: Slice, alphabet: Sequence[int], phase: str) -> dict[str, Any]:
        runtime = slice_runtime[item.identifier]
        auxiliary: PointedQuartic = runtime["auxiliary"]
        basis = runtime["basis"]
        counts = {
            "coefficient_vectors": 0,
            "identity_relations": 0,
            "exceptional_inverse_images": 0,
            "zero_parameters": 0,
            "prior_H200000_or_T0_parameters": 0,
            "singular_parameters": 0,
            "visible_quartic_images": 0,
            "visible_covariant_sign_pairs": 0,
            "accepted_incidences": 0,
        }
        for vector, auxiliary_point in recursive_group_elements(
            auxiliary.weierstrass_coefficients, basis, alphabet
        ):
            if phase == "quinary-expansion" and all(
                scalar in TERNARY_ALPHABET for scalar in vector
            ):
                continue
            counts["coefficient_vectors"] += 1
            if auxiliary_point is None:
                counts["identity_relations"] += 1
                continue
            inverse = auxiliary.inverse(auxiliary_point)
            if inverse is None:
                counts["exceptional_inverse_images"] += 1
                continue
            signed_parameter, normalized_ordinate = inverse
            parameter = abs(signed_parameter)
            if parameter == 0:
                counts["zero_parameters"] += 1
                continue
            if parameter in prior_parameters:
                counts["prior_H200000_or_T0_parameters"] += 1
                continue
            try:
                homogenized_discriminant(parameter)
            except ValueError:
                counts["singular_parameters"] += 1
                continue
            _, raw_ordinate = normalized_to_raw_point(
                item, (signed_parameter, normalized_ordinate)
            )
            raw_x = item.raw_x(signed_parameter)
            canonical_quartic = primitive_quartic_coefficients(
                RANK21_CONSTRUCTION, parameter
            )
            if raw_ordinate**2 != quartic_value(canonical_quartic, raw_x):
                raise AssertionError("a canonicalized orbit point missed the family")
            if raw_ordinate == 0:
                counts["singular_parameters"] += 1
                continue
            visible_abscissas, visible_image_x = canonical_visible_data(parameter)
            if raw_x in visible_abscissas:
                counts["visible_quartic_images"] += 1
                continue
            direct_image = quartic_point_to_short_jacobian(
                RANK21_CONSTRUCTION,
                parameter,
                (raw_x, raw_ordinate),
            )
            if direct_image[0] in visible_image_x:
                counts["visible_covariant_sign_pairs"] += 1
                continue
            candidates.setdefault(parameter, []).append(
                {
                    "slice_id": item.identifier,
                    "phase": phase,
                    "coefficient_vector": list(vector),
                    "signed_parameter_before_even_symmetry": rational_record(
                        signed_parameter
                    ),
                    "forced_quartic_point": point_record((raw_x, raw_ordinate)),
                    "direct_jacobian_point": point_record(direct_image),
                }
            )
            counts["accepted_incidences"] += 1
        return counts

    for index, item in enumerate(slices, 1):
        counts = enumerate_slice(item, TERNARY_ALPHABET, "ternary")
        slice_generation[item.identifier] = {
            "ternary": counts,
            "stable_rank": slice_runtime[item.identifier]["stable_rank"],
        }
        print(
            f"ternary {index}/{len(slices)} {item.identifier} "
            f"vectors={counts['coefficient_vectors']} "
            f"accepted={counts['accepted_incidences']}",
            flush=True,
        )

    primes = primes_up_to(PROXY_PRIME_BOUND)
    proxies: dict[Fraction, dict[str, Any]] = {}
    for index, parameter in enumerate(sorted(candidates), 1):
        proxies[parameter] = exact_radical_proxy(parameter, primes)
        if index % 100 == 0:
            print(f"ternary proxy {index}/{len(candidates)}", flush=True)

    for item in slices:
        own_parameters = [
            parameter
            for parameter, incidences in candidates.items()
            if any(record["slice_id"] == item.identifier for record in incidences)
        ]
        best_proxy = (
            min(proxies[parameter]["log_radical_upper_proxy"] for parameter in own_parameters)
            if own_parameters
            else None
        )
        generation = slice_generation[item.identifier]
        warranted = expansion_is_warranted(
            stable_rank=generation["stable_rank"],
            accepted_incidences=generation["ternary"]["accepted_incidences"],
            best_proxy=best_proxy,
        )
        generation["ternary_best_proxy"] = best_proxy
        generation["quinary_expansion_warranted"] = warranted
        generation["quinary_expansion_rule"] = {
            "minimum_stable_rank": EXPANSION_MINIMUM_RANK,
            "minimum_accepted_ternary_incidences": EXPANSION_MINIMUM_ACCEPTED_INCIDENCES,
            "strict_ternary_best_proxy_gate": EXPANSION_PROXY_GATE,
        }
        if warranted:
            expansion_counts = enumerate_slice(
                item, QUINARY_ALPHABET, "quinary-expansion"
            )
            generation["quinary_expansion"] = expansion_counts
            print(
                f"expanded {item.identifier} "
                f"new_vectors={expansion_counts['coefficient_vectors']} "
                f"accepted={expansion_counts['accepted_incidences']}",
                flush=True,
            )

    for parameter in candidates:
        if parameter not in proxies:
            proxies[parameter] = exact_radical_proxy(parameter, primes)

    ordered_parameters = sorted(
        candidates,
        key=lambda parameter: (
            proxies[parameter]["log_radical_upper_proxy"],
            max(abs(parameter.numerator), parameter.denominator),
            parameter,
        ),
    )
    candidate_records = []
    maximum_pool_size = 0
    maximum_slice_diversity = 0
    for parameter in ordered_parameters:
        incidences = candidates[parameter]
        forced_x = {Q(record["direct_jacobian_point"]["x"]) for record in incidences}
        slice_diversity = len({record["slice_id"] for record in incidences})
        visible_direct_x = {
            quartic_point_to_short_jacobian(
                RANK21_CONSTRUCTION, parameter, point
            )[0]
            for point in primitive_visible_points(RANK21_CONSTRUCTION, parameter)
        }
        pool_size = len(visible_direct_x | forced_x)
        maximum_pool_size = max(maximum_pool_size, pool_size)
        maximum_slice_diversity = max(maximum_slice_diversity, slice_diversity)
        candidate_records.append(
            {
                "parameter": rational_record(parameter),
                "projective_height": max(
                    abs(parameter.numerator), parameter.denominator
                ),
                "incidence_count": len(incidences),
                "distinct_slice_count": slice_diversity,
                "visible_plus_forced_direct_image_x_count": pool_size,
                "source_slice_ids": sorted(
                    {record["slice_id"] for record in incidences}
                ),
                "log_radical_upper_proxy": proxies[parameter][
                    "log_radical_upper_proxy"
                ],
            }
        )

    plausible = [
        parameter
        for parameter in ordered_parameters
        if proxies[parameter]["log_radical_upper_proxy"]
        < PROXY_EXACT_CONDUCTOR_GATE
    ]
    conductor_parameters = list(plausible)
    for parameter in ordered_parameters[:LEADER_EXACT_CONDUCTOR_COUNT]:
        if parameter not in conductor_parameters:
            conductor_parameters.append(parameter)
    conductors = {}
    for parameter in conductor_parameters:
        record = conductor_attempt(
            parameter,
            timeout=args.conductor_timeout,
            stack_bytes=args.stack_bytes,
        )
        conductors[rational_record(parameter)] = record
        print(
            f"conductor T={rational_record(parameter)} status={record['status']}",
            flush=True,
        )
        if record.get("below_strict_log_conductor_target", False):
            print(
                f"ALERT strict-subtarget T={rational_record(parameter)} "
                f"logN={record['log_conductor']}",
                flush=True,
            )

    subtarget = [
        parameter
        for parameter in conductor_parameters
        if conductors[rational_record(parameter)].get(
            "below_strict_log_conductor_target", False
        )
    ]
    expanded_slice_ids = [
        identifier
        for identifier, record in slice_generation.items()
        if record["quinary_expansion_warranted"]
    ]
    incidence_payload = "\n".join(
        f"{rational_record(parameter)}|{record['slice_id']}|"
        f"{','.join(map(str, record['coefficient_vector']))}|{record['phase']}"
        for parameter in sorted(candidates)
        for record in sorted(
            candidates[parameter],
            key=lambda item: (
                item["slice_id"], item["phase"], item["coefficient_vector"]
            ),
        )
    )

    artifact = {
        "schema_version": 1,
        "status": "bounded_nonproductive_auxiliary_orbits_complete",
        "goal": "rank>=21 with log conductor<182.72, or rank>=30",
        "target_hit": False,
        "source": {
            "primary_source": PRIMARY_SOURCE,
            "root_tuple": [int(root) for root in RANK21_CONSTRUCTION.roots],
            "constructor_parameter": rational_record(T0),
            "published_parameter": "14721/376",
            "published_point_count": len(RANK21_PUBLISHED_POINTS),
            "model_change_short_to_published_u_r_s_t": [
                rational_record(value)
                for value in MODEL_CHANGE_SHORT_TO_MINIMAL
            ],
            "all_published_points_exactly_round_tripped": True,
        },
        "published_point_reconstruction": {
            "selected_base_index_one_based": 1,
            "selected_sign": 1,
            "exact_oriented_generic_matches": 6,
            "generic_abscissa_matches": 10,
            "preimage_count": len(reconstruction),
            "preimage_sha256": EXPECTED_PREIMAGE_SHA256,
            "preimages": [point_record(point) for point in reconstruction],
            "accidental_preimage_count": len(accidental),
            "accidental_preimage_sha256": EXPECTED_ACCIDENTAL_SHA256,
            "accidental_preimages": [point_record(point) for point in accidental],
            "all_convention_trials": list(convention_trials),
        },
        "scope": {
            "constructed_slope_count": len(all_slices),
            "searched_nonproductive_slice_count": len(slices),
            "searched_slice_ids": [item.identifier for item in slices],
            "excluded_productive_slice": {
                "slice_id": productive[0].identifier,
                "source_x": rational_record(PRODUCTIVE_SOURCE_X),
                "slope": PRODUCTIVE_SLOPE,
                "intercept": rational_record(PRODUCTIVE_INTERCEPT),
            },
            "slice_height": SLICE_HEIGHT,
            "base_coefficient_alphabet": list(TERNARY_ALPHABET),
            "conditional_expansion_alphabet": list(QUINARY_ALPHABET),
            "no_other_parameters_or_slices_searched": True,
        },
        "slices": [
            {
                **record,
                "generation": slice_generation[record["id"]],
            }
            for record in slice_records
        ],
        "decontamination": {
            "global_H200000_or_T0_parameter_count": len(prior_parameters),
            "global_H200000_or_T0_parameter_sha256": rational_digest(
                sorted(prior_parameters)
            ),
            "exact_filters": [
                "canonicalize by the exact even-family symmetry T -> -T",
                "exclude the union of every searched H=200000 parameter and T0",
                "exclude zero and exact primitive-discriminant zeros",
                "replay normalized ordinate on the original rank-21 quartic",
                "exclude every visible quartic abscissa",
                "exclude every visible covariant image up to sign by exact X",
                "deduplicate exact canonical parameters across slices",
            ],
            "all_filters_applied_before_proxy_and_conductor_work": True,
        },
        "generation": {
            "expanded_slice_ids": expanded_slice_ids,
            "unique_accepted_parameter_count": len(candidates),
            "accepted_incidence_count": sum(len(value) for value in candidates.values()),
            "parameter_sha256": rational_digest(sorted(candidates)),
            "parameter_incidence_sha256": hashlib.sha256(
                incidence_payload.encode()
            ).hexdigest(),
            "maximum_distinct_slice_count_at_one_parameter": maximum_slice_diversity,
            "maximum_visible_plus_forced_direct_image_x_count": maximum_pool_size,
            "rank_signal_at_least_18_in_constructed_exact_point_pools": (
                maximum_pool_size >= 18
            ),
            "all_candidate_records_ordered_by_proxy": candidate_records,
        },
        "proxy_filter": {
            "trial_prime_bound": PROXY_PRIME_BOUND,
            "strict_exact_conductor_gate": PROXY_EXACT_CONDUCTOR_GATE,
            "minimum_log_radical_upper_proxy": (
                proxies[ordered_parameters[0]]["log_radical_upper_proxy"]
                if ordered_parameters
                else None
            ),
            "below_gate_count": len(plausible),
            "below_gate_parameters": [
                rational_record(parameter) for parameter in plausible
            ],
            "top_64_full_proxy_records": [
                {
                    "parameter": rational_record(parameter),
                    **proxies[parameter],
                }
                for parameter in ordered_parameters[:64]
            ],
            "scope_warning": (
                "the discriminant radical is a conductor heuristic; it is not "
                "an equality or a rank certificate"
            ),
        },
        "exact_conductors": {
            "all_proxy_below_190_attempted": True,
            "calibration_leader_count": LEADER_EXACT_CONDUCTOR_COUNT,
            "attempted_parameter_count": len(conductor_parameters),
            "records": conductors,
            "strict_subtarget_parameters": [
                rational_record(parameter) for parameter in subtarget
            ],
        },
        "conclusion": {
            "target_hit": False,
            "strict_subtarget_parameter_count": len(subtarget),
            "rank_signal_at_least_18": maximum_pool_size >= 18,
            "statement": (
                "The declared nonproductive H=200000 auxiliary-slice panel "
                "and every warranted coefficient orbit produced no breakthrough."
            ),
            "scope_warning": (
                "This is exhaustive only for the 21 declared bounded slice runs "
                "and their declared coefficient alphabets."
            ),
        },
        "checkpoint": {
            "exact_H200000_checkpoint_written_before_orbit": True,
            "final_artifact_replaced_checkpoint_after_completion": True,
        },
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "pari": pari_version(),
            "script_sha256": file_sha256(Path(__file__).resolve()),
            "wall_seconds": time.monotonic() - started,
            "all_subprocesses_foreground": True,
            "maximum_one_shot_timeout_seconds": max(
                args.slice_timeout, args.height_timeout, args.conductor_timeout
            ),
            "one_attempt_per_declared_subprocess": True,
        },
    }
    if subtarget or maximum_pool_size >= 18:
        artifact["target_hit"] = False
    write_artifact(args.output, artifact)
    print(
        f"wrote {args.output} candidates={len(candidates)} "
        f"proxy<190={len(plausible)} subtarget={len(subtarget)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
