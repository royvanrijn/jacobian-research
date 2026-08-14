#!/usr/bin/env python3
"""Constructive genus-one slices through Fermigier's rank-22 near miss.

At the normalized record parameter ``T0=39508/39`` this script reconstructs
an exact quartic preimage of every one of Fermigier's 22 published points.
The reconstruction uses the degree-two covariant map in a transparent way.
If ``phi:C->J`` is that map and ``O`` is the first visible quartic point, the
chosen isomorphism ``psi:C->J`` is characterized by

``phi(Q) = phi(O) + 2*psi(Q)``.

Thus the preimage of a published point ``P`` is recovered by factoring the
quartic equation for ``phi(O)+2P``.  The convention is selected independently
from all 13 possible generic base points and both signs: it is the unique one
maximizing overlap with the published visible-section coordinates.

An independent, pinned ``H=10^6`` replay on the record quartic returns 27
abscissas: the 13 selected generic sections and 14 further abscissas.  One of
the latter is the ``T -> -T`` conjugate of Mestre's extra generic section, so
the artifact distinguishes 14 search-relative accidentals from 13 genuinely
non-generic specializations.  All 14 nevertheless define the requested 28
exact slices

``x = +/- T + n``.

Their leading degree-six terms cancel, leaving nonsingular quartics in ``T``.
Each receives one bounded ``H=200000`` ``hyperellratpoints`` search, with a foreground
process-group timeout and no retry.  Every returned point is mapped exactly;
all distinct new parameters receive an exact conductor attempt.  Only a
below-target candidate whose exact point pool has stable numerical rank at
least 21 can trigger small-prime saturation and an exact finite-reduction
independence attempt.

The record replay is also turned into a portable exact certificate: a stable
21-point numerical subset is used only to discover a lower-height basis via
``ellsaturation(E,P,2)``; exact finite-reduction images certify those 21
returned points together with published point P14 as independent.  This
proves rank at least 22 for Fermigier's benchmark, but its exact conductor has
``log(N)>182.72`` and is therefore explicitly not a target hit.

The slice searches remain bounded construction experiments.  Numerical
height ranks, root numbers, timeouts, and the absence of points outside the
declared boxes are not rank certificates.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
from typing import Any, Sequence

import sympy as sp

from alternate_quartic_covers import point_on_short_curve, short_add
from certify_nagao_rank17_frontier import exact_log_conductor_certificate
from ek_k3 import rational_square_root, rational_to_string
from extend_nagao_u42_frontier import saturate_exact_basis
from fermigier_mestre import FermigierMestreFamily, NORMALIZED_RECORD_PARAMETER
from mod2_reduction_independence import (
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from pari_bridge import minimal_curve_data, pari_version
from search_extra_points import (
    gp_rational,
    parse_point_vector,
    signless_quartic_points,
)
from triage_nagao_rank13_finalists import (
    height_matrix_replay,
    point_digest,
    stable_height_rank,
)
from verify_fermigier_benchmark import PUBLISHED_MODEL
from verify_fermigier_rank22_points import PUBLISHED_POINTS


Q = Fraction
T0 = NORMALIZED_RECORD_PARAMETER
TARGET_LOG_CONDUCTOR = Decimal("182.72")
MODEL_CHANGE_SHORT_TO_MINIMAL = (
    Q(14, 507),
    Q(49, 771147),
    Q(7, 507),
    Q(1372, 130323843),
)
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/elliptic_fermigier_rank22_accidental_slices.json"
)
DEFAULT_ESCALATION_CHECKPOINT = Path(
    "artifacts/generated-results/"
    "elliptic_fermigier_3115_3_h1000000_checkpoint.json"
)
ESCALATION_CHECKPOINT_SHA256 = (
    "fe4acebc95c02b61229107ef7220974e83696ded87f97629a8980cef3ba46471"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_fermigier_rank22_accidental_slices.py"
)

# These are exactly the abscissas returned by the pinned H=10^6 record-fiber
# replay after removing the 13 selected generic-section abscissas.  The order
# is fixed independently of PARI's output order.
RECORD_SEARCH_RELATIVE_ACCIDENTAL_X: tuple[Fraction, ...] = (
    Q(32609, 39),
    Q(256714, 39),
    Q(-51427, 195),
    Q(-2086, 195),
    Q(144932, 195),
    Q(618754, 195),
    Q(-790314, 871),
    Q(211132, 1131),
    Q(580286, 1209),
    Q(-433096, 2379),
    Q(464062, 2379),
    Q(719604, 9425),
    Q(-70102, 63375),
    Q(-74800, 98943),
)
CONJUGATE_EXTRA_SECTION_X = Q(144932, 195)


# Coefficients of R_T(x), first by ascending x-power and then ascending
# T-power.  They replay FermigierMestreFamily.quartic_coefficients exactly.
FERMIGIER_BIVARIATE_COEFFICIENTS: tuple[tuple[Fraction, ...], ...] = (
    (
        Q(18103855887324900),
        Q(0),
        Q(102302344648),
        Q(0),
        Q(-879500),
        Q(0),
        Q(1),
    ),
    (Q(-257843832010380), Q(0), Q(-650709150), Q(0), Q(1860)),
    (Q(1195214262641), Q(0), Q(1718550), Q(0), Q(-2)),
    (Q(-2051321790), Q(0), Q(-1860)),
    (Q(1149050), Q(0), Q(1)),
)


@dataclass(frozen=True)
class Slice:
    accidental_label: str
    source_point: tuple[Fraction, Fraction]
    slope: int
    intercept: Fraction
    coefficients: tuple[Fraction, ...]

    @property
    def identifier(self) -> str:
        slope = "p1" if self.slope == 1 else "m1"
        return f"{self.accidental_label.lower()}_{slope}"

    def x_value(self, parameter: Fraction) -> Fraction:
        return Q(self.slope) * Q(parameter) + self.intercept


def trim_polynomial(coefficients: Sequence[Fraction]) -> tuple[Fraction, ...]:
    answer = [Q(value) for value in coefficients]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def poly_add(*polynomials: Sequence[Fraction]) -> tuple[Fraction, ...]:
    length = max((len(polynomial) for polynomial in polynomials), default=1)
    return trim_polynomial(
        tuple(
            sum(
                (
                    Q(polynomial[index])
                    if index < len(polynomial)
                    else Q(0)
                )
                for polynomial in polynomials
            )
            for index in range(length)
        )
    )


def poly_multiply(*polynomials: Sequence[Fraction]) -> tuple[Fraction, ...]:
    answer = (Q(1),)
    for polynomial in polynomials:
        product = [Q(0)] * (len(answer) + len(polynomial) - 1)
        for left_degree, left in enumerate(answer):
            for right_degree, right in enumerate(polynomial):
                product[left_degree + right_degree] += Q(left) * Q(right)
        answer = trim_polynomial(product)
    return answer


def poly_evaluate(
    coefficients: Sequence[Fraction], value: Fraction
) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * Q(value) + Q(coefficient)
    return answer


def minimum_to_short(
    point: tuple[Fraction, Fraction]
) -> tuple[Fraction, Fraction]:
    """Invert PARI's pinned minimal-model change exactly."""

    u_value, r_value, s_value, t_value = MODEL_CHANGE_SHORT_TO_MINIMAL
    x_value, y_value = (Q(value) for value in point)
    answer = (
        u_value**2 * x_value + r_value,
        u_value**3 * y_value + s_value * u_value**2 * x_value + t_value,
    )
    coefficients = FermigierMestreFamily.coefficients(T0)
    if answer[1] ** 2 != (
        answer[0] ** 3 + coefficients[3] * answer[0] + coefficients[4]
    ):
        raise AssertionError("the inverse minimal-model change left the curve")
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
    a1, a2, a3, a4, a6 = (Q(value) for value in PUBLISHED_MODEL)
    if answer[1] ** 2 + a1 * answer[0] * answer[1] + a3 * answer[1] != (
        answer[0] ** 3 + a2 * answer[0] ** 2 + a4 * answer[0] + a6
    ):
        raise AssertionError("the minimal-model change left the curve")
    return answer


def short_negate(point: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return point[0], -point[1]


def short_sum(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction] | None],
) -> tuple[Fraction, Fraction] | None:
    answer = None
    for point in points:
        answer = short_add(coefficients, answer, point)
    return answer


def _rational_linear_roots(expression: sp.Expr, symbol: sp.Symbol) -> tuple[Fraction, ...]:
    roots: set[Fraction] = set()
    for factor, _ in sp.factor_list(expression)[1]:
        if sp.degree(factor, symbol) != 1:
            continue
        root = sp.solve(factor, symbol)[0]
        if root.is_Rational:
            roots.add(Q(int(root.p), int(root.q)))
    return tuple(sorted(roots))


def short_halves(
    coefficients: Sequence[Fraction],
    target: tuple[Fraction, Fraction] | None,
) -> tuple[tuple[Fraction, Fraction] | None, ...]:
    """Return every rational half of ``target`` and check it exactly.

    If ``target=(X,Y)`` and ``R=(r,s)``, the duplication formula gives

    ``r^4-4X*r^3-2A*r^2-(4AX+8B)r+A^2-4BX=0``.

    The quartic need not be squarefree when the target is 2-torsion, so this
    routine deliberately deduplicates linear roots rather than rejecting
    repeated factors.
    """

    coefficient_a, coefficient_b = (Q(value) for value in coefficients[3:])
    symbol = sp.symbols("r")
    if target is None:
        torsion_polynomial = (
            symbol**3
            + sp.Rational(coefficient_a.numerator, coefficient_a.denominator)
            * symbol
            + sp.Rational(coefficient_b.numerator, coefficient_b.denominator)
        )
        answers: list[tuple[Fraction, Fraction] | None] = [None]
        answers.extend(
            (root, Q(0))
            for root in _rational_linear_roots(torsion_polynomial, symbol)
        )
        return tuple(answers)

    target_x, target_y = (Q(value) for value in target)
    polynomial = (
        symbol**4
        - 4 * sp.Rational(target_x.numerator, target_x.denominator) * symbol**3
        - 2
        * sp.Rational(coefficient_a.numerator, coefficient_a.denominator)
        * symbol**2
        - sp.Rational(
            (4 * coefficient_a * target_x + 8 * coefficient_b).numerator,
            (4 * coefficient_a * target_x + 8 * coefficient_b).denominator,
        )
        * symbol
        + sp.Rational(
            (coefficient_a**2 - 4 * coefficient_b * target_x).numerator,
            (coefficient_a**2 - 4 * coefficient_b * target_x).denominator,
        )
    )
    answers_set: set[tuple[Fraction, Fraction]] = set()
    for root in _rational_linear_roots(polynomial, symbol):
        if target_y:
            ordinate = (
                root**3
                - 3 * target_x * root**2
                - coefficient_a * root
                - coefficient_a * target_x
                - 2 * coefficient_b
            ) / (2 * target_y)
            possible_ordinates = (ordinate,)
        else:
            square_root = rational_square_root(
                root**3 + coefficient_a * root + coefficient_b
            )
            possible_ordinates = (
                () if square_root is None else tuple({square_root, -square_root})
            )
        for ordinate in possible_ordinates:
            candidate = (root, ordinate)
            if point_on_short_curve(coefficients, candidate) and short_add(
                coefficients, candidate, candidate
            ) == target:
                answers_set.add(candidate)
    return tuple(sorted(answers_set))


def quartic_group_pullback(
    parameter: Fraction,
    point: tuple[Fraction, Fraction],
) -> tuple[Fraction, Fraction] | None:
    """Map a quartic point through the base-point isomorphism ``psi``.

    The first selected generic section is the origin.  On every fiber used
    below there is a unique rational half of ``phi(Q)-phi(O)``; ambiguity is
    rejected rather than silently choosing a torsion translate.
    """

    parameter = Q(parameter)
    coefficients = FermigierMestreFamily.coefficients(parameter)
    base = FermigierMestreFamily.known_quartic_points(parameter)[0]
    if point == base:
        return None
    image = FermigierMestreFamily.quartic_point_to_jacobian(parameter, point)
    base_image = FermigierMestreFamily.quartic_point_to_jacobian(parameter, base)
    difference = short_add(coefficients, image, short_negate(base_image))
    halves = short_halves(coefficients, difference)
    affine_halves = tuple(candidate for candidate in halves if candidate is not None)
    if len(affine_halves) != 1:
        raise AssertionError(
            f"expected a unique affine group pullback, found {len(affine_halves)}"
        )
    answer = affine_halves[0]
    if short_add(coefficients, answer, answer) != difference:
        raise AssertionError("the exact group pullback failed its doubling replay")
    return answer


def generic_group_seed_points(
    parameter: Fraction,
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Return 13 affine seeds in the based ``psi`` group coordinate.

    The twelve non-base selected sections give twelve pullbacks.  The other
    seed is ``D=psi(iota(O))=-phi(O)``, which makes choosing just one sheet of
    every quartic abscissa rank-safe because ``psi(iota(Q))=D-psi(Q)``.
    """

    parameter = Q(parameter)
    quartic_points = FermigierMestreFamily.known_quartic_points(parameter)
    base_image = FermigierMestreFamily.quartic_point_to_jacobian(
        parameter, quartic_points[0]
    )
    pullbacks = tuple(
        quartic_group_pullback(parameter, point) for point in quartic_points[1:]
    )
    if any(point is None for point in pullbacks):
        raise AssertionError("a non-base generic section pulled back to zero")
    answer = (short_negate(base_image),) + tuple(
        point for point in pullbacks if point is not None
    )
    if len(answer) != 13 or len({point[0] for point in answer}) != 13:
        raise AssertionError("generic group seeds collided modulo inverse")
    return answer


def _covariant_x_polynomial(
    target_x: Fraction,
) -> tuple[Fraction, ...]:
    """Return descending coefficients of X*R_T(x)-36*g_T(x)."""

    a_value, b_value, c_value, d_value, e_value = (
        FermigierMestreFamily.quartic_coefficients(T0)
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
    """Recover the unique rational affine quartic point above ``target``."""

    symbol = sp.symbols("x")
    coefficients = _covariant_x_polynomial(target[0])
    expression = sum(
        sp.Rational(value.numerator, value.denominator) * symbol ** (4 - index)
        for index, value in enumerate(coefficients)
    )
    rational_roots = _rational_linear_roots(expression, symbol)
    answers: list[tuple[Fraction, Fraction]] = []
    for x_value in rational_roots:
        z_value = rational_square_root(
            FermigierMestreFamily.quartic_value(T0, x_value)
        )
        if z_value is None:
            continue
        for signed_z in (z_value, -z_value):
            point = (x_value, signed_z)
            if FermigierMestreFamily.quartic_point_to_jacobian(T0, point) == target:
                answers.append(point)
    if len(answers) != 1:
        raise AssertionError(
            f"expected one rational covariant preimage, found {len(answers)}"
        )
    return answers[0]


def reconstruct_with_convention(
    base_index: int, sign: int
) -> tuple[tuple[Fraction, Fraction], ...]:
    generic = FermigierMestreFamily.known_quartic_points(T0)
    if base_index < 0 or base_index >= len(generic) or sign not in (-1, 1):
        raise ValueError("invalid reconstruction convention")
    coefficients = FermigierMestreFamily.coefficients(T0)
    offset = FermigierMestreFamily.quartic_point_to_jacobian(
        T0, generic[base_index]
    )
    answers = []
    for published in PUBLISHED_POINTS:
        short_point = minimum_to_short(published)
        doubled = short_add(coefficients, short_point, short_point)
        if doubled is None:
            raise AssertionError("a published point unexpectedly doubled to zero")
        if sign == -1:
            doubled = short_negate(doubled)
        target = short_add(coefficients, offset, doubled)
        if target is None:
            raise AssertionError("the covariant target unexpectedly became infinity")
        answers.append(invert_covariant_target(target))
    return tuple(answers)


def select_reconstruction_convention() -> tuple[
    int,
    int,
    tuple[tuple[Fraction, Fraction], ...],
    tuple[dict[str, Any], ...],
]:
    """Select the unique convention maximizing generic-section x-overlap."""

    generic = FermigierMestreFamily.known_quartic_points(T0)
    generic_x = {point[0] for point in generic}
    trials = []
    complete: dict[tuple[int, int], tuple[tuple[Fraction, Fraction], ...]] = {}
    for base_index in range(len(generic)):
        for sign in (1, -1):
            preimages = reconstruct_with_convention(base_index, sign)
            complete[(base_index, sign)] = preimages
            exact_matches = sum(point in set(generic) for point in preimages)
            abscissa_matches = sum(point[0] in generic_x for point in preimages)
            trials.append(
                {
                    "base_index_one_based": base_index + 1,
                    "sign": sign,
                    "exact_oriented_generic_matches": exact_matches,
                    "generic_abscissa_matches": abscissa_matches,
                }
            )
    maximum = max(int(trial["generic_abscissa_matches"]) for trial in trials)
    winners = [
        trial for trial in trials if trial["generic_abscissa_matches"] == maximum
    ]
    if len(winners) != 1:
        raise AssertionError("the reconstruction convention is not uniquely selected")
    winner = winners[0]
    key = (int(winner["base_index_one_based"]) - 1, int(winner["sign"]))
    if key != (0, 1) or maximum != 11:
        raise AssertionError("the pinned Fermigier reconstruction convention changed")
    return key[0], key[1], complete[key], tuple(trials)


def slice_polynomial(slope: int, intercept: Fraction) -> tuple[Fraction, ...]:
    if slope not in (-1, 1):
        raise ValueError("this construction only uses slopes +/-1")
    linear = (Q(intercept), Q(slope))
    x_power = (Q(1),)
    answer = (Q(0),)
    for coefficient_polynomial in FERMIGIER_BIVARIATE_COEFFICIENTS:
        answer = poly_add(
            answer,
            poly_multiply(coefficient_polynomial, x_power),
        )
        x_power = poly_multiply(x_power, linear)
    answer = trim_polynomial(answer)
    if len(answer) - 1 != 4:
        raise AssertionError("the +/-1 Fermigier slice did not become quartic")
    return answer


def published_accidental_points(
    reconstruction: Sequence[tuple[Fraction, Fraction]],
) -> tuple[tuple[str, tuple[Fraction, Fraction]], ...]:
    generic_x = {
        point[0] for point in FermigierMestreFamily.known_quartic_points(T0)
    }
    accidentals = tuple(
        (f"P{index}", point)
        for index, point in enumerate(reconstruction, start=1)
        if point[0] not in generic_x
    )
    if [label for label, _ in accidentals] != [
        "P6",
        *[f"P{index}" for index in range(13, 23)],
    ]:
        raise AssertionError("the accidental rank-22 point labels changed")
    return accidentals


def canonical_signless_points(
    points: Sequence[tuple[Fraction, Fraction]],
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Choose the positive ordinate for each abscissa, sorted exactly."""

    by_x: dict[Fraction, Fraction] = {}
    for x_value, z_value in points:
        x_value, z_value = Q(x_value), Q(z_value)
        if x_value in by_x and by_x[x_value] ** 2 != z_value**2:
            raise AssertionError("one quartic abscissa had inconsistent ordinates")
        by_x[x_value] = abs(z_value)
    return tuple(sorted(by_x.items()))


def record_search_relative_accidentals(
    record_points: Sequence[tuple[Fraction, Fraction]],
) -> tuple[tuple[str, tuple[Fraction, Fraction]], ...]:
    signless = canonical_signless_points(record_points)
    generic_x = {
        point[0] for point in FermigierMestreFamily.known_quartic_points(T0)
    }
    by_x = {point[0]: point for point in signless if point[0] not in generic_x}
    if set(by_x) != set(RECORD_SEARCH_RELATIVE_ACCIDENTAL_X):
        raise AssertionError("the pinned H=10^6 accidental abscissas changed")
    return tuple(
        (f"A{index:02d}", by_x[x_value])
        for index, x_value in enumerate(
            RECORD_SEARCH_RELATIVE_ACCIDENTAL_X, start=1
        )
    )


def build_slices(
    accidentals: Sequence[tuple[str, tuple[Fraction, Fraction]]],
) -> tuple[Slice, ...]:
    slices = []
    symbol = sp.symbols("T")
    for label, point in accidentals:
        for slope in (-1, 1):
            intercept = point[0] - slope * T0
            coefficients = slice_polynomial(slope, intercept)
            if poly_evaluate(coefficients, T0) != point[1] ** 2:
                raise AssertionError("a slice missed its source point")
            expression = sum(
                sp.Rational(value.numerator, value.denominator) * symbol**degree
                for degree, value in enumerate(coefficients)
            )
            if sp.discriminant(expression, symbol) == 0:
                raise AssertionError("a priority Fermigier slice is singular")
            slices.append(Slice(label, point, slope, intercept, coefficients))
    return tuple(slices)


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


def run_gp_once(
    program: str, *, timeout: float, stack_bytes: int
) -> tuple[str | None, dict[str, Any]]:
    if timeout <= 0 or timeout > 60:
        raise ValueError("a GP subprocess timeout must lie in (0,60]")
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


def gp_polynomial(coefficients: Sequence[Fraction]) -> str:
    return "+".join(
        f"{gp_rational(Q(coefficient))}*x^{degree}"
        for degree, coefficient in enumerate(coefficients)
        if coefficient
    ) or "0"


def search_polynomial(
    coefficients: Sequence[Fraction],
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    program = "\n".join(
        (
            f"Q={gp_polynomial(coefficients)};gettime();",
            f"R=hyperellratpoints(Q,{height_bound});",
            'print("PARI_MILLISECONDS ",gettime());',
            'print("POINTS_BEGIN");print(R);print("POINTS_END");',
            "quit",
        )
    ) + "\n"
    output, process_record = run_gp_once(
        program, timeout=timeout, stack_bytes=stack_bytes
    )
    record = {
        **process_record,
        "height_bound": height_bound,
        "timeout_seconds": timeout,
        "pari_stack_bytes": stack_bytes,
        "retried": False,
    }
    if output is None:
        return (), record
    marker = re.search(r"POINTS_BEGIN\n(.*?)\nPOINTS_END", output, re.DOTALL)
    milliseconds = re.search(r"^PARI_MILLISECONDS (\d+)$", output, re.MULTILINE)
    if marker is None or milliseconds is None:
        raise AssertionError("PARI omitted slice-search markers")
    points = parse_point_vector(marker.group(1))
    record.update(
        {
            "pari_milliseconds": int(milliseconds.group(1)),
            "signed_point_count": len(points),
            "distinct_parameter_count": len(signless_quartic_points(points)),
        }
    )
    return points, record


def search_slice(
    slice_data: Slice,
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    return search_polynomial(
        slice_data.coefficients,
        height_bound=height_bound,
        timeout=timeout,
        stack_bytes=stack_bytes,
    )


def search_record_quartic(
    *, height_bound: int, timeout: float, stack_bytes: int
) -> tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    # FermigierMestreFamily returns descending x-coefficients, whereas
    # gp_polynomial consumes ascending coefficients.
    return search_polynomial(
        tuple(reversed(FermigierMestreFamily.quartic_coefficients(T0))),
        height_bound=height_bound,
        timeout=timeout,
        stack_bytes=stack_bytes,
    )


def search_specialized_quartic(
    parameter: Fraction,
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]:
    return search_polynomial(
        tuple(
            reversed(FermigierMestreFamily.quartic_coefficients(Q(parameter)))
        ),
        height_bound=height_bound,
        timeout=timeout,
        stack_bytes=stack_bytes,
    )


def projective_height(value: Fraction) -> int:
    value = Q(value)
    return max(abs(value.numerator), value.denominator)


def point_record(point: tuple[Fraction, Fraction]) -> dict[str, str | bool]:
    return {
        "jacobian_x": rational_to_string(point[0]),
        "jacobian_y": rational_to_string(point[1]),
        "exact_membership_checked": True,
    }


def specialized_quartic_screen(
    parameter: Fraction,
    *,
    height_bound: int,
    search_timeout: float,
    height_timeout: float,
    precisions: tuple[int, ...],
    stack_bytes: int,
) -> tuple[dict[str, Any], tuple[tuple[Fraction, Fraction], ...]]:
    """Search one specialized quartic and replay the direct-image height rank."""

    parameter = Q(parameter)
    raw_points, search = search_specialized_quartic(
        parameter,
        height_bound=height_bound,
        timeout=search_timeout,
        stack_bytes=stack_bytes,
    )
    record: dict[str, Any] = {
        "parameter_t": rational_to_string(parameter),
        "search": search,
        "height_rank": None,
    }
    if search["status"] != "completed":
        return record, ()
    signless = canonical_signless_points(raw_points)
    seeds = FermigierMestreFamily.known_jacobian_points(parameter)
    pool = list(seeds)
    seen_image_x = {point[0] for point in seeds}
    new_records = []
    for quartic_point in signless:
        image = FermigierMestreFamily.quartic_point_to_jacobian(
            parameter, quartic_point
        )
        if image[0] in seen_image_x:
            continue
        seen_image_x.add(image[0])
        pool.append(image)
        new_records.append(
            {
                "quartic_x": rational_to_string(quartic_point[0]),
                "quartic_z": rational_to_string(quartic_point[1]),
                "direct_covariant_image": point_record(image),
            }
        )
    exact_pool = tuple(pool)
    runs = height_matrix_replay(
        FermigierMestreFamily.coefficients(parameter),
        exact_pool,
        precisions=precisions,
        timeout=height_timeout,
        stack_bytes=stack_bytes,
    )
    numerical_rank = stable_height_rank(runs)
    indices = tuple(runs[-1]["subset_indices_one_based"])
    selected = tuple(exact_pool[index - 1] for index in indices)
    record.update(
        {
            "signed_point_count": len(raw_points),
            "signless_abscissa_count": len(signless),
            "generic_direct_seed_count": len(seeds),
            "new_distinct_direct_images_modulo_sign": len(new_records),
            "new_points": new_records,
            "exact_pool_point_count": len(exact_pool),
            "height_rank": {
                "status": "completed",
                "stable_numerical_rank": numerical_rank,
                "precision_runs": list(runs),
                "selected_subset_indices_one_based": list(indices),
                "selected_point_sha256": point_digest(selected),
                "selected_points": [point_record(point) for point in selected],
                "scope_warning": "numerical triage only",
            },
        }
    )
    return record, selected


def load_escalation_checkpoint(
    path: Path,
    *,
    precisions: tuple[int, ...],
    height_timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    """Ingest and exactly replay the durable H=10^6 T=3115/3 checkpoint."""

    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != ESCALATION_CHECKPOINT_SHA256:
        raise AssertionError("the pinned 3115/3 escalation checkpoint hash changed")
    checkpoint = json.loads(raw)
    results = checkpoint.get("results", [])
    if len(results) != 1:
        raise AssertionError("the escalation checkpoint must contain one result")
    source = results[0]
    expected = {
        "t": "3115/3",
        "quartic_height_bound": 1_000_000,
        "signed_quartic_points_found": 134,
        "distinct_quartic_x_values": 67,
        "visible_section_x_values_found": 13,
        "new_x_values_beyond_visible_sections": 54,
        "stable_numerical_rank": 15,
    }
    if any(source.get(key) != value for key, value in expected.items()):
        raise AssertionError("the pinned 3115/3 escalation counts changed")
    selected = tuple(
        (Q(record["jacobian_x"]), Q(record["jacobian_y"]))
        for record in source["explicit_numerically_independent_subset"]
    )
    coefficients = FermigierMestreFamily.coefficients(Q(3115, 3))
    if len(selected) != 15 or any(
        not point_on_short_curve(coefficients, point) for point in selected
    ):
        raise AssertionError("the escalation subset failed exact membership")
    replay_runs = height_matrix_replay(
        coefficients,
        selected,
        precisions=precisions,
        timeout=height_timeout,
        stack_bytes=stack_bytes,
    )
    if stable_height_rank(replay_runs) != 15:
        raise AssertionError("the durable escalation subset lost numerical rank 15")
    return {
        "status": "completed checkpoint ingested and exact subset replayed",
        **expected,
        "checkpoint_path": str(path),
        "checkpoint_sha256": digest,
        "source_engine_script_sha256": checkpoint["script_sha256"],
        "source_reproducing_command": checkpoint["reproducing_command"],
        "source_software": checkpoint["software"],
        "exact_selected_point_count": len(selected),
        "exact_selected_point_sha256": point_digest(selected),
        "exact_selected_points": [point_record(point) for point in selected],
        "selected_subset_height_replay": list(replay_runs),
        "scope_warning": (
            "the bounded enumeration counts are pinned by the hashed durable "
            "checkpoint; this ingestion rechecks its 15 selected points exactly "
            "and replays their numerical height rank, but does not rerun H=10^6"
        ),
    }


def candidate_point_pools(
    slices: Sequence[Slice],
    search_results: Sequence[tuple[tuple[tuple[Fraction, Fraction], ...], dict[str, Any]]],
) -> tuple[
    dict[Fraction, tuple[tuple[Fraction, Fraction], ...]],
    dict[Fraction, list[dict[str, Any]]],
]:
    """Exact-map non-generic slice points and aggregate by the even fiber |T|."""

    quartic_by_parameter: dict[Fraction, dict[Fraction, dict[str, Any]]] = (
        defaultdict(dict)
    )
    for slice_data, (raw_points, _) in zip(slices, search_results, strict=True):
        for original_parameter, ordinate in canonical_signless_points(raw_points):
            normalized_parameter = abs(Q(original_parameter))
            if normalized_parameter in (Q(0), abs(T0)):
                continue
            if FermigierMestreFamily.discriminant_factor(normalized_parameter) == 0:
                continue
            x_value = slice_data.x_value(Q(original_parameter))
            if ordinate**2 != FermigierMestreFamily.quartic_value(
                normalized_parameter, x_value
            ):
                raise AssertionError("an exact slice point missed the normalized fiber")
            # R_T and E_T are even in T.  Include both orientations of the
            # extra Mestre section so that a T-sign conjugate is not counted
            # as a specialization-only point.
            generic_x = {
                point[0]
                for signed_parameter in (
                    normalized_parameter,
                    -normalized_parameter,
                )
                for point in FermigierMestreFamily.known_quartic_points(
                    signed_parameter
                )
            }
            if x_value in generic_x:
                continue
            existing = quartic_by_parameter[normalized_parameter].get(x_value)
            if existing is None:
                quartic_by_parameter[normalized_parameter][x_value] = {
                    "ordinate": abs(ordinate),
                    "slice_ids": {slice_data.identifier},
                    "signed_parameters": {Q(original_parameter)},
                }
            else:
                if Q(existing["ordinate"]) ** 2 != ordinate**2:
                    raise AssertionError("a shared slice point had inconsistent ordinates")
                existing["slice_ids"].add(slice_data.identifier)
                existing["signed_parameters"].add(Q(original_parameter))

    pools: dict[Fraction, tuple[tuple[Fraction, Fraction], ...]] = {}
    provenance: dict[Fraction, list[dict[str, Any]]] = {}
    for parameter, points_by_x in quartic_by_parameter.items():
        seeds = generic_group_seed_points(parameter)
        seen_direct_x = {
            point[0] for point in FermigierMestreFamily.known_jacobian_points(parameter)
        }
        seen_pullback_x = {point[0] for point in seeds}
        extras = []
        records = []
        for x_value, source in sorted(points_by_x.items()):
            ordinate = Q(source["ordinate"])
            image = FermigierMestreFamily.quartic_point_to_jacobian(
                parameter, (x_value, ordinate)
            )
            pullback = quartic_group_pullback(parameter, (x_value, ordinate))
            if pullback is None:
                raise AssertionError("a non-generic slice point pulled back to zero")
            direct_image_is_new = image[0] not in seen_direct_x
            seen_direct_x.add(image[0])
            pullback_is_new = pullback[0] not in seen_pullback_x
            if pullback_is_new:
                seen_pullback_x.add(pullback[0])
                extras.append(pullback)
            records.append(
                {
                    "quartic_x": rational_to_string(x_value),
                    "quartic_z": rational_to_string(ordinate),
                    "slice_ids": sorted(source["slice_ids"]),
                    "signed_slice_parameters": [
                        rational_to_string(value)
                        for value in sorted(source["signed_parameters"])
                    ],
                    "direct_covariant_image": point_record(image),
                    "direct_image_new_modulo_sign": direct_image_is_new,
                    "basepoint_group_pullback": point_record(pullback),
                    "group_pullback_new_modulo_inverse": pullback_is_new,
                    "exact_group_pullback_doubling_checked": True,
                }
            )
        pools[parameter] = tuple(seeds) + tuple(extras)
        provenance[parameter] = records
    return pools, provenance


def slice_parameter_incidences(
    slice_data: Slice,
    raw_points: Sequence[tuple[Fraction, Fraction]],
) -> list[dict[str, Any]]:
    """Retain every signless parameter incidence with its exact classification."""

    records = []
    for original_parameter, ordinate in canonical_signless_points(raw_points):
        normalized = abs(original_parameter)
        x_value = slice_data.x_value(original_parameter)
        parameter_sign = -1 if original_parameter < 0 else 1
        canonical_slope = slice_data.slope * parameter_sign
        if original_parameter and (
            Q(canonical_slope) * normalized + slice_data.intercept != x_value
        ):
            raise AssertionError("negative-parameter slice canonicalization failed")
        classification = "non-generic candidate"
        if normalized == 0:
            classification = "zero parameter excluded"
        elif normalized == abs(T0):
            classification = "record source fiber excluded"
        elif FermigierMestreFamily.discriminant_factor(normalized) == 0:
            classification = "singular fiber excluded"
        else:
            generic_x = {
                point[0]
                for signed_parameter in (normalized, -normalized)
                for point in FermigierMestreFamily.known_quartic_points(
                    signed_parameter
                )
            }
            if x_value in generic_x:
                classification = "generic-or-T-sign-conjugate collision excluded"
        records.append(
            {
                "signed_parameter_t": rational_to_string(original_parameter),
                "normalized_even_fiber_t": rational_to_string(normalized),
                "raw_slice_slope": slice_data.slope,
                "canonical_even_fiber_slope": canonical_slope,
                "slice_intercept": rational_to_string(slice_data.intercept),
                "quartic_x": rational_to_string(x_value),
                "quartic_z": rational_to_string(ordinate),
                "classification": classification,
            }
        )
    return records


def conductor_probe(
    parameter: Fraction, *, timeout: float, stack_bytes: int
) -> dict[str, Any]:
    try:
        data = minimal_curve_data(
            FermigierMestreFamily.coefficients(parameter),
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        return {
            "status": "completed",
            **data,
            "below_strict_log_conductor_target": (
                Decimal(data["log_conductor"]) < TARGET_LOG_CONDUCTOR
            ),
        }
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "timeout_seconds": timeout}
    except (RuntimeError, AssertionError, ValueError) as error:
        return {"status": "error", "error": str(error)[:1000]}


def finite_reduction_attempt(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    saturation_timeout: float,
    stack_bytes: int,
    certificate_prime_bound: int,
) -> dict[str, Any]:
    saturated, saturation = saturate_exact_basis(
        coefficients,
        points,
        prime_bound=20,
        timeout=saturation_timeout,
        stack_bytes=stack_bytes,
    )
    signatures = find_mod2_reduction_certificate(
        coefficients, saturated, prime_bound=certificate_prime_bound
    )
    exact_rank = combined_mod2_rank(signatures, len(saturated))
    certified = exact_rank == len(saturated)
    return {
        "status": "certified" if certified else "bounded-search-rank-deficient",
        "input_point_count": len(points),
        "saturated_point_count": len(saturated),
        "saturation": saturation,
        "saturated_point_sha256": point_digest(saturated),
        "certificate_prime_bound": certificate_prime_bound,
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
        "certified_algebraic_rank_lower_bound": len(saturated) if certified else None,
        "saturated_basis": [point_record(point) for point in saturated],
    }


def exact_mod2_certificate_record(
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
        "proof_note": (
            "full rank of the exact reduction images modulo 2 proves the "
            "listed rational points Z-independent by infinite descent on any "
            "putative integer relation"
            if certified
            else "the declared prime box did not prove independence"
        ),
    }


def certify_record_replay(
    record_points: Sequence[tuple[Fraction, Fraction]],
    *,
    precisions: tuple[int, ...],
    height_timeout: float,
    saturation_timeout: float,
    certificate_prime_bound: int,
    conductor_timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    """Produce the standalone exact rank-22 certificate at Fermigier's T0."""

    signless = canonical_signless_points(record_points)
    if len(record_points) != 54 or len(signless) != 27:
        raise AssertionError("the pinned record-fiber point counts changed")
    coefficients = FermigierMestreFamily.coefficients(T0)
    direct_images = tuple(
        FermigierMestreFamily.quartic_point_to_jacobian(T0, point)
        for point in signless
    )
    if len({point[0] for point in direct_images}) != len(direct_images):
        raise AssertionError("the record replay had duplicate direct images modulo sign")
    height_runs = height_matrix_replay(
        coefficients,
        direct_images,
        precisions=precisions,
        timeout=height_timeout,
        stack_bytes=stack_bytes,
    )
    numerical_rank = stable_height_rank(height_runs)
    if numerical_rank != 21:
        raise AssertionError("the pinned record-fiber numerical rank changed")
    indices = tuple(height_runs[-1]["subset_indices_one_based"])
    selected = tuple(direct_images[index - 1] for index in indices)
    saturated, saturation = saturate_exact_basis(
        coefficients,
        selected,
        prime_bound=2,
        timeout=saturation_timeout,
        stack_bytes=stack_bytes,
    )
    if len(saturated) != 21:
        raise AssertionError("2-saturation changed the record basis length")

    published_p14_short = minimum_to_short(PUBLISHED_POINTS[13])
    augmented = tuple(saturated) + (published_p14_short,)
    if len({point[0] for point in augmented}) != 22:
        raise AssertionError("P14 duplicated a saturated record point modulo sign")
    constructive_certificate = exact_mod2_certificate_record(
        coefficients,
        augmented,
        prime_bound=certificate_prime_bound,
    )
    if (
        constructive_certificate["status"] != "certified"
        or constructive_certificate["combined_exact_rank_over_F2"] != 22
    ):
        raise AssertionError("the pinned rank-22 finite-reduction certificate failed")
    published_short_points = tuple(minimum_to_short(point) for point in PUBLISHED_POINTS)
    published_certificate = exact_mod2_certificate_record(
        coefficients,
        published_short_points,
        prime_bound=certificate_prime_bound,
    )
    if (
        published_certificate["status"] != "certified"
        or published_certificate["combined_exact_rank_over_F2"] != 22
    ):
        raise AssertionError("the canonical published-point rank-22 certificate failed")
    conductor = conductor_probe(
        T0, timeout=conductor_timeout, stack_bytes=stack_bytes
    )
    if conductor.get("status") != "completed":
        raise AssertionError("the record conductor computation did not complete")
    if conductor["below_strict_log_conductor_target"]:
        raise AssertionError("the record curve unexpectedly met the conductor target")
    return {
        "record_signless_points": [
            {
                "quartic_x": rational_to_string(point[0]),
                "quartic_z": rational_to_string(point[1]),
                "direct_covariant_image": point_record(image),
            }
            for point, image in zip(signless, direct_images, strict=True)
        ],
        "direct_image_pool_count": len(direct_images),
        "height_rank": {
            "status": "completed",
            "stable_numerical_rank": numerical_rank,
            "precision_runs": list(height_runs),
            "selected_subset_indices_one_based": list(indices),
            "selected_point_sha256": point_digest(selected),
            "scope_warning": "numerical selection only; not the independence proof",
        },
        "two_saturation_discovery": saturation,
        "rank22_basis": [point_record(point) for point in augmented],
        "rank22_basis_construction": (
            "the exact 21-point output of ellsaturation(E,P,2), followed by "
            "published minimal-model point P14 transformed exactly to the short model"
        ),
        "published_augmenting_point_label": "P14",
        "published_augmenting_point_short": point_record(published_p14_short),
        "canonical_published_point_rank22_basis": [
            point_record(point) for point in published_short_points
        ],
        "canonical_published_point_finite_reduction_certificate": (
            published_certificate
        ),
        "constructive_saturated_plus_P14_finite_reduction_certificate": (
            constructive_certificate
        ),
        "conductor": conductor,
        "target_classification": (
            "unconditional rank >= 22 near miss; log conductor is above the "
            "strict 182.72 target"
        ),
    }


def parse_precisions(value: str) -> tuple[int, ...]:
    try:
        precisions = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("precisions must be integers") from error
    if len(precisions) < 2 or tuple(sorted(set(precisions))) != precisions:
        raise argparse.ArgumentTypeError("provide increasing distinct precisions")
    return precisions


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record-height", type=int, default=1_000_000)
    parser.add_argument("--record-timeout", type=float, default=60.0)
    parser.add_argument("--slice-height", type=int, default=200_000)
    parser.add_argument("--slice-timeout", type=float, default=10.0)
    parser.add_argument("--specialization-height", type=int, default=50_000)
    parser.add_argument("--specialization-timeout", type=float, default=15.0)
    parser.add_argument("--conductor-timeout", type=float, default=8.0)
    parser.add_argument("--height-timeout", type=float, default=20.0)
    parser.add_argument("--height-precisions", type=parse_precisions, default=(72, 120))
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=500)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--escalation-checkpoint",
        type=Path,
        default=root / DEFAULT_ESCALATION_CHECKPOINT,
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root / DEFAULT_OUTPUT,
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if (
        args.record_height != 1_000_000
        or args.slice_height != 200_000
        or args.specialization_height != 50_000
    ):
        raise SystemExit(
            "the canonical replay pins record H=1000000, slice H=200000, "
            "and specialization H=50000"
        )
    if (
        not 0 < args.record_timeout <= 60
        or not 0 < args.slice_timeout <= 60
        or not 0 < args.specialization_timeout <= 60
    ):
        raise SystemExit("all bounded search timeouts must lie in (0,60]")
    if min(
        args.conductor_timeout,
        args.height_timeout,
        args.saturation_timeout,
    ) <= 0 or max(
        args.conductor_timeout,
        args.height_timeout,
        args.saturation_timeout,
    ) > 60:
        raise SystemExit("all remaining timeouts must lie in (0,60]")
    if args.stack_bytes < 8_000_000 or args.certificate_prime_bound < 3:
        raise SystemExit("invalid stack or finite-reduction prime bound")

    # Pin the exact family/model transformation and every source point before
    # any bounded search.
    if tuple(short_to_minimum(minimum_to_short(point)) for point in PUBLISHED_POINTS) != PUBLISHED_POINTS:
        raise AssertionError("the pinned model transformation failed its round trip")
    base_index, sign, reconstruction, convention_trials = (
        select_reconstruction_convention()
    )
    generic = FermigierMestreFamily.known_quartic_points(T0)
    generic_x = {point[0] for point in generic}
    preimage_records = []
    for index, (published, preimage) in enumerate(
        zip(PUBLISHED_POINTS, reconstruction, strict=True), start=1
    ):
        short_point = minimum_to_short(published)
        offset = FermigierMestreFamily.quartic_point_to_jacobian(
            T0, generic[base_index]
        )
        doubled = short_add(
            FermigierMestreFamily.coefficients(T0), short_point, short_point
        )
        if doubled is None:
            raise AssertionError("a published point doubled to zero")
        if sign == -1:
            doubled = short_negate(doubled)
        target = short_add(FermigierMestreFamily.coefficients(T0), offset, doubled)
        if target is None or FermigierMestreFamily.quartic_point_to_jacobian(
            T0, preimage
        ) != target:
            raise AssertionError("a reconstructed preimage failed exact replay")
        group_pullback = quartic_group_pullback(T0, preimage)
        if group_pullback != short_point:
            raise AssertionError(
                "a reconstructed preimage failed exact basepoint pullback"
            )
        generic_index = next(
            (
                position
                for position, point in enumerate(generic, start=1)
                if point[0] == preimage[0]
            ),
            None,
        )
        preimage_records.append(
            {
                "label": f"P{index}",
                "published_minimal_point": {
                    "x": rational_to_string(published[0]),
                    "y": rational_to_string(published[1]),
                },
                "short_point": point_record(short_point),
                "covariant_target": point_record(target),
                "basepoint_group_pullback": point_record(group_pullback),
                "quartic_preimage": {
                    "x": rational_to_string(preimage[0]),
                    "z": rational_to_string(preimage[1]),
                    "exact_quartic_membership_checked": True,
                    "exact_covariant_target_checked": True,
                },
                "classification": (
                    "generic-section abscissa"
                    if preimage[0] in generic_x
                    else "accidental rational preimage"
                ),
                "generic_section_index_one_based": generic_index,
                "generic_ordinate_orientation": (
                    "same"
                    if preimage in generic
                    else "opposite"
                    if (preimage[0], -preimage[1]) in generic
                    else None
                ),
            }
        )

    record_points, record_search = search_record_quartic(
        height_bound=args.record_height,
        timeout=args.record_timeout,
        stack_bytes=args.stack_bytes,
    )
    if record_search["status"] != "completed":
        raise RuntimeError(
            f"the pinned H=10^6 record replay did not complete: {record_search}"
        )
    record_accidentals = record_search_relative_accidentals(record_points)
    if len(record_accidentals) != 14:
        raise AssertionError("the record replay did not produce 14 priority sources")
    conjugate_point = next(
        point
        for _, point in record_accidentals
        if point[0] == CONJUGATE_EXTRA_SECTION_X
    )
    standard_pullbacks = tuple(
        quartic_group_pullback(T0, point) for point in generic[1:]
    )
    if any(point is None for point in standard_pullbacks):
        raise AssertionError("a non-base generic section pulled back to zero")
    conjugate_pullback = quartic_group_pullback(T0, conjugate_point)
    if conjugate_pullback is None:
        raise AssertionError("the conjugate extra section pulled back to zero")
    relation_terms = (
        standard_pullbacks[5],
        short_negate(standard_pullbacks[6]),
        short_negate(standard_pullbacks[7]),
        standard_pullbacks[8],
        short_negate(standard_pullbacks[11]),
        conjugate_pullback,
    )
    if short_sum(FermigierMestreFamily.coefficients(T0), relation_terms) is not None:
        raise AssertionError("the exact conjugate-section dependency changed")
    conjugate_dependency = {
        "identity": "G6-G7-G8+G9-G12+G13=O",
        "definitions": (
            "G1..G12 are psi of selected generic sections 2..13; G13 is "
            "psi of the positive-ordinate T-sign-conjugate extra source"
        ),
        "checked_by_exact_short_weierstrass_addition": True,
        "scope": "record fiber T=39508/39 only",
    }
    record_certificate = certify_record_replay(
        record_points,
        precisions=args.height_precisions,
        height_timeout=args.height_timeout,
        saturation_timeout=args.saturation_timeout,
        certificate_prime_bound=args.certificate_prime_bound,
        conductor_timeout=args.conductor_timeout,
        stack_bytes=args.stack_bytes,
    )
    print(
        "record H=1000000 signed=54 x=27 stable_rank=21 "
        "exact_certified_rank=22",
        flush=True,
    )

    published_accidentals = published_accidental_points(reconstruction)
    published_labels_by_x = {
        point[0]: label for label, point in published_accidentals
    }
    source_records = []
    for label, point in record_accidentals:
        source_records.append(
            {
                "label": label,
                "quartic_x": rational_to_string(point[0]),
                "quartic_z": rational_to_string(point[1]),
                "source_classification": (
                    "T-sign-conjugate generic section"
                    if point[0] == CONJUGATE_EXTRA_SECTION_X
                    else "specialization-only accidental at the record fiber"
                ),
                "published_preimage_label_at_same_x": published_labels_by_x.get(
                    point[0]
                ),
                "exact_quartic_membership_checked": (
                    point[1] ** 2
                    == FermigierMestreFamily.quartic_value(T0, point[0])
                ),
            }
        )
    overlap_labels = sorted(
        (record["published_preimage_label_at_same_x"] for record in source_records),
        key=lambda label: int(label[1:]) if label is not None else 10**9,
    )
    overlap_labels = [label for label in overlap_labels if label is not None]
    if overlap_labels != ["P6", "P13", "P16", "P17", "P18", "P19"]:
        raise AssertionError("the uniform/published accidental overlap changed")

    slices = build_slices(record_accidentals)
    search_results = []
    slice_records = []
    for slice_data in slices:
        points, search = search_slice(
            slice_data,
            height_bound=args.slice_height,
            timeout=args.slice_timeout,
            stack_bytes=args.stack_bytes,
        )
        search_results.append((points, search))
        slice_records.append(
            {
                "slice_id": slice_data.identifier,
                "source_label": slice_data.accidental_label,
                "source_quartic_point": {
                    "x": rational_to_string(slice_data.source_point[0]),
                    "z": rational_to_string(slice_data.source_point[1]),
                },
                "equation": (
                    f"z^2=R_T({slice_data.slope}*T+"
                    f"({rational_to_string(slice_data.intercept)}))"
                ),
                "slope": slice_data.slope,
                "intercept": rational_to_string(slice_data.intercept),
                "auxiliary_coefficients_ascending": [
                    rational_to_string(value) for value in slice_data.coefficients
                ],
                "auxiliary_degree": len(slice_data.coefficients) - 1,
                "auxiliary_genus": 1,
                "source_parameter_replay": True,
                "source_classification": (
                    "T-sign-conjugate generic section"
                    if slice_data.source_point[0] == CONJUGATE_EXTRA_SECTION_X
                    else "specialization-only accidental at the record fiber"
                ),
                "parameter_incidences": slice_parameter_incidences(
                    slice_data, points
                ),
                "search": search,
            }
        )
        print(
            f"slice={slice_data.identifier} status={search['status']} "
            f"parameters={search.get('distinct_parameter_count', 0)}",
            flush=True,
        )

    pools, extra_provenance = candidate_point_pools(slices, search_results)
    candidate_records = []
    target_hits = []
    alternative_rank_hits = []
    for parameter in sorted(
        pools,
        key=lambda value: (
            -len(
                {
                    slice_id
                    for extra in extra_provenance[value]
                    for slice_id in extra["slice_ids"]
                }
            ),
            -sum(
                bool(extra["group_pullback_new_modulo_inverse"])
                for extra in extra_provenance[value]
            ),
            projective_height(value),
            value,
        ),
    ):
        pool = pools[parameter]
        extras = extra_provenance[parameter]
        slice_ids = sorted(
            {
                slice_id
                for extra in extras
                for slice_id in extra["slice_ids"]
            }
        )
        pullback_x = {
            extra["basepoint_group_pullback"]["jacobian_x"]
            for extra in extras
            if extra["basepoint_group_pullback"] is not None
        }
        new_direct_count = sum(
            bool(extra["direct_image_new_modulo_sign"]) for extra in extras
        )
        new_group_count = sum(
            bool(extra["group_pullback_new_modulo_inverse"]) for extra in extras
        )
        conductor = conductor_probe(
            parameter,
            timeout=args.conductor_timeout,
            stack_bytes=args.stack_bytes,
        )
        record: dict[str, Any] = {
            "parameter_t": rational_to_string(parameter),
            "projective_height": projective_height(parameter),
            "generic_seed_count": 13,
            "generic_seed_expected_numerical_rank": 12,
            "slice_ids": slice_ids,
            "distinct_slice_incidence_count": len(slice_ids),
            "distinct_forced_quartic_abscissa_count": len(extras),
            "distinct_group_pullback_classes_modulo_inverse": len(pullback_x),
            "new_group_pullback_classes_modulo_inverse": new_group_count,
            "new_direct_covariant_images_modulo_sign": new_direct_count,
            "exact_pool_point_count": len(pool),
            "forced_accidental_points": extras,
            "conductor_probe": conductor,
            "height_rank": {
                "status": "not-run",
                "reason": (
                    "single-slice yield with fewer than 21 exact pool points"
                ),
            },
            "finite_reduction_certificate": None,
        }
        if len(slice_ids) >= 2 or len(pool) >= 21:
            try:
                runs = height_matrix_replay(
                    FermigierMestreFamily.coefficients(parameter),
                    pool,
                    precisions=args.height_precisions,
                    timeout=args.height_timeout,
                    stack_bytes=args.stack_bytes,
                )
                numerical_rank = stable_height_rank(runs)
                indices = runs[-1]["subset_indices_one_based"]
                selected = tuple(pool[index - 1] for index in indices)
                record["height_rank"] = {
                    "status": "completed",
                    "stable_numerical_rank": numerical_rank,
                    "precision_runs": list(runs),
                    "selected_subset_indices_one_based": indices,
                    "selected_point_sha256": point_digest(selected),
                    "status_note": "numerical triage only",
                }
                below_target_rank_lane = (
                    numerical_rank >= 21
                    and conductor.get("status") == "completed"
                    and conductor["below_strict_log_conductor_target"]
                )
                alternative_rank_lane = numerical_rank >= 30
                if below_target_rank_lane or alternative_rank_lane:
                    certificate = finite_reduction_attempt(
                        FermigierMestreFamily.coefficients(parameter),
                        selected,
                        saturation_timeout=args.saturation_timeout,
                        stack_bytes=args.stack_bytes,
                        certificate_prime_bound=args.certificate_prime_bound,
                    )
                    record["finite_reduction_certificate"] = certificate
                    if certificate["status"] == "certified":
                        certified_rank = int(
                            certificate["certified_algebraic_rank_lower_bound"]
                        )
                        if below_target_rank_lane and certified_rank >= 21:
                            record["exact_log_conductor_bound"] = (
                                exact_log_conductor_certificate(
                                    int(conductor["conductor"])
                                )
                            )
                            target_hits.append(rational_to_string(parameter))
                        if alternative_rank_lane and certified_rank >= 30:
                            alternative_rank_hits.append(
                                rational_to_string(parameter)
                            )
            except (
                subprocess.TimeoutExpired,
                RuntimeError,
                AssertionError,
                ValueError,
            ) as error:
                record["height_rank"] = {
                    "status": (
                        "timeout"
                        if isinstance(error, subprocess.TimeoutExpired)
                        else "error"
                    ),
                    "error": str(error)[:1000],
                }
        candidate_records.append(record)
        print(
            f"T={parameter} slices={len(slice_ids)} x={len(extras)} "
            f"pool={len(pool)} "
            f"conductor={conductor['status']} "
            f"logN={conductor.get('log_conductor')}",
            flush=True,
        )

    completed_conductors = [
        record
        for record in candidate_records
        if record["conductor_probe"]["status"] == "completed"
    ]
    below_target = [
        record
        for record in completed_conductors
        if record["conductor_probe"]["below_strict_log_conductor_target"]
    ]
    specialization_screens = []
    for candidate in below_target:
        parameter = Q(candidate["parameter_t"])
        try:
            screen, selected = specialized_quartic_screen(
                parameter,
                height_bound=args.specialization_height,
                search_timeout=args.specialization_timeout,
                height_timeout=args.height_timeout,
                precisions=args.height_precisions,
                stack_bytes=args.stack_bytes,
            )
            numerical_rank = (
                int(screen["height_rank"]["stable_numerical_rank"])
                if screen["height_rank"] is not None
                else 0
            )
            screen["finite_reduction_certificate"] = None
            if numerical_rank >= 21:
                certificate = finite_reduction_attempt(
                    FermigierMestreFamily.coefficients(parameter),
                    selected,
                    saturation_timeout=args.saturation_timeout,
                    stack_bytes=args.stack_bytes,
                    certificate_prime_bound=args.certificate_prime_bound,
                )
                screen["finite_reduction_certificate"] = certificate
                if (
                    certificate["status"] == "certified"
                    and int(certificate["certified_algebraic_rank_lower_bound"])
                    >= 21
                ):
                    if candidate["parameter_t"] not in target_hits:
                        target_hits.append(candidate["parameter_t"])
                    screen["exact_log_conductor_bound"] = (
                        exact_log_conductor_certificate(
                            int(candidate["conductor_probe"]["conductor"])
                        )
                    )
        except (
            subprocess.TimeoutExpired,
            RuntimeError,
            AssertionError,
            ValueError,
        ) as error:
            screen = {
                "parameter_t": rational_to_string(parameter),
                "search": {
                    "status": (
                        "timeout"
                        if isinstance(error, subprocess.TimeoutExpired)
                        else "error"
                    ),
                    "error": str(error)[:1000],
                },
                "height_rank": None,
                "finite_reduction_certificate": None,
            }
        candidate["specialized_quartic_height_50000_screen"] = screen
        specialization_screens.append(screen)
        print(
            f"specialization T={parameter} H={args.specialization_height} "
            f"status={screen['search']['status']} "
            f"new={screen.get('new_distinct_direct_images_modulo_sign')} "
            f"rank={(screen.get('height_rank') or {}).get('stable_numerical_rank')}",
            flush=True,
        )

    h50_frontier = next(
        (
            screen
            for screen in specialization_screens
            if screen["parameter_t"] == "3115/3"
            and (screen.get("height_rank") or {}).get("stable_numerical_rank")
            == 15
        ),
        None,
    )
    if h50_frontier is None:
        raise AssertionError("the T=3115/3 H=50000 rank-15 frontier was not recovered")
    escalation_checkpoint = load_escalation_checkpoint(
        args.escalation_checkpoint,
        precisions=args.height_precisions,
        height_timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    h50_frontier["height_1000000_escalation"] = escalation_checkpoint
    multi_slice_intersections = [
        {
            "parameter_t": record["parameter_t"],
            "slice_ids": record["slice_ids"],
            "distinct_slice_incidence_count": record[
                "distinct_slice_incidence_count"
            ],
            "distinct_forced_quartic_abscissa_count": record[
                "distinct_forced_quartic_abscissa_count"
            ],
            "distinct_group_pullback_classes_modulo_inverse": record[
                "distinct_group_pullback_classes_modulo_inverse"
            ],
            "height_rank": record["height_rank"],
            "conductor_probe": record["conductor_probe"],
        }
        for record in candidate_records
        if record["distinct_slice_incidence_count"] >= 2
    ]
    incidence_histogram: dict[str, int] = defaultdict(int)
    for record in candidate_records:
        incidence_histogram[str(record["distinct_slice_incidence_count"])] += 1
    candidate_stable_ranks = [
        int(record["height_rank"]["stable_numerical_rank"])
        for record in candidate_records
        if record["height_rank"]["status"] == "completed"
    ]
    script_path = Path(__file__).resolve()
    command = " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv])
    artifact = {
        "schema_version": 2,
        "status": (
            "unconditional exact rank >= 22 certificate for Fermigier's above-target "
            "benchmark, exact published-point pullbacks, and a bounded 28-slice "
            "construction/conductor experiment"
        ),
        "record_parameter_normalized_T": rational_to_string(T0),
        "published_model": list(PUBLISHED_MODEL),
        "short_to_minimal_change": [
            rational_to_string(value) for value in MODEL_CHANGE_SHORT_TO_MINIMAL
        ],
        "preimage_map_convention": {
            "identity": "phi(Q)=phi(O)+2*psi(Q)",
            "base_generic_section_index_one_based": base_index + 1,
            "sign": sign,
            "unique_maximum_generic_abscissa_overlap": 11,
            "all_convention_trials": list(convention_trials),
            "mathematical_scope": (
                "these are preimages under the base-point isomorphism psi, not "
                "direct preimages of P under the degree-two covariant phi"
            ),
        },
        "published_point_preimages": preimage_records,
        "generic_sections": {
            "count": len(generic),
            "points": [
                {
                    "index_one_based": index,
                    "x": rational_to_string(point[0]),
                    "z": rational_to_string(point[1]),
                }
                for index, point in enumerate(generic, start=1)
            ],
            "published_preimages_on_generic_abscissas": sum(
                record["classification"] == "generic-section abscissa"
                for record in preimage_records
            ),
        },
        "published_accidental_preimages": {
            "count": sum(
                record["classification"] == "accidental rational preimage"
                for record in preimage_records
            ),
            "labels": [
                record["label"]
                for record in preimage_records
                if record["classification"] == "accidental rational preimage"
            ],
        },
        "record_fiber_height_1000000_replay": {
            "search": record_search,
            "signed_point_count_pinned": 54,
            "signless_abscissa_count_pinned": 27,
            "positive_parameter_generic_catalog_count": 13,
            "extra_vs_positive_parameter_catalog_count": 14,
            "T_sign_conjugate_generic_source_count": 1,
            "specialization_only_accidental_source_count": 13,
            "priority_source_count": len(record_accidentals),
            "priority_sources": source_records,
            "T_sign_conjugate_exact_dependency": conjugate_dependency,
            "overlap_with_published_preimage_labels": overlap_labels,
            "exact_rank_certificate": record_certificate,
        },
        "slices": slice_records,
        "candidate_conductor_screen": candidate_records,
        "multi_slice_intersections": multi_slice_intersections,
        "specialized_quartic_target_triage": {
            "height_50000_screen_parameters": [
                screen["parameter_t"] for screen in specialization_screens
            ],
            "screen_location": (
                "candidate_conductor_screen[*]."
                "specialized_quartic_height_50000_screen"
            ),
            "escalated_parameter_t": "3115/3",
            "escalation_reason": (
                "unique H=50000 frontier with stable numerical rank 15, a "
                "three-rank gain over the visible rank-12 baseline"
            ),
            "escalation_record_location": (
                "candidate_conductor_screen entry T=3115/3, within "
                "specialized_quartic_height_50000_screen."
                "height_1000000_escalation"
            ),
            "escalation_checkpoint_sha256": escalation_checkpoint[
                "checkpoint_sha256"
            ],
            "conclusion": (
                "the H=1000000 escalation remains stable numerical rank 15; "
                "no exact target-certificate trigger was reached"
            ),
        },
        "summary": {
            "published_points_reconstructed_exactly": len(preimage_records),
            "generic_abscissa_preimages": 11,
            "published_accidental_preimages": 11,
            "record_replay_signed_points": 54,
            "record_replay_signless_abscissas": 27,
            "record_replay_extra_vs_positive_catalog": 14,
            "record_replay_T_sign_conjugate_generic_sources": 1,
            "record_replay_specialization_accidental_sources": 13,
            "record_replay_exact_certified_rank_lower_bound": 22,
            "record_replay_meets_strict_conductor_target": False,
            "genus_one_slices": len(slices),
            "slices_from_specialization_accidental_sources": 26,
            "slices_from_T_sign_conjugate_generic_source": 2,
            "slice_searches_completed": sum(
                record["search"]["status"] == "completed"
                for record in slice_records
            ),
            "slice_search_timeouts": sum(
                record["search"]["status"] == "timeout"
                for record in slice_records
            ),
            "distinct_new_parameters": len(candidate_records),
            "canonical_candidate_parameters": [
                record["parameter_t"] for record in candidate_records
            ],
            "slice_incidence_multiplicity_histogram": dict(
                sorted(incidence_histogram.items(), key=lambda item: int(item[0]))
            ),
            "maximum_slice_incidence_count": max(
                (
                    record["distinct_slice_incidence_count"]
                    for record in candidate_records
                ),
                default=0,
            ),
            "multi_slice_parameter_count": len(multi_slice_intersections),
            "conductor_calls_completed": len(completed_conductors),
            "completed_conductors_below_strict_target": len(below_target),
            "specialized_quartic_height_50000_screens": len(
                specialization_screens
            ),
            "specialized_quartic_height_50000_completed": sum(
                screen["search"]["status"] == "completed"
                for screen in specialization_screens
            ),
            "specialized_quartic_height_50000_ranks": {
                screen["parameter_t"]: (
                    None
                    if screen["height_rank"] is None
                    else screen["height_rank"]["stable_numerical_rank"]
                )
                for screen in specialization_screens
            },
            "specialized_quartic_maximum_stable_numerical_rank": max(
                (
                    int(screen["height_rank"]["stable_numerical_rank"])
                    for screen in specialization_screens
                    if screen["height_rank"] is not None
                ),
                default=None,
            ),
            "specialized_quartic_frontier_parameter_t": "3115/3",
            "specialized_quartic_frontier_H1000000_rank": 15,
            "maximum_forced_accidental_abscissas": max(
                (
                    record["distinct_forced_quartic_abscissa_count"]
                    for record in candidate_records
                ),
                default=0,
            ),
            "maximum_forced_pool_numerical_rank_before_specialized_search": max(
                candidate_stable_ranks, default=None
            ),
            "maximum_exact_pool_point_count": max(
                (record["exact_pool_point_count"] for record in candidate_records),
                default=13,
            ),
            "finite_reduction_certificates_triggered": sum(
                record["finite_reduction_certificate"] is not None
                for record in candidate_records
            ),
            "target_hits": target_hits,
            "alternative_rank_hits": alternative_rank_hits,
        },
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "low_conductor_certified_hits": target_hits,
            "alternative_rank_certified_hits": alternative_rank_hits,
            "record_benchmark_classification": (
                "certified rank >= 22 but above strict log-conductor target"
            ),
        },
        "parameters": {
            "record_height_bound": args.record_height,
            "record_timeout_seconds": args.record_timeout,
            "slice_height_bound": args.slice_height,
            "slice_timeout_seconds_per_process": args.slice_timeout,
            "specialization_height_bound": args.specialization_height,
            "specialization_timeout_seconds_per_process": (
                args.specialization_timeout
            ),
            "escalation_checkpoint": str(args.escalation_checkpoint),
            "escalation_checkpoint_sha256": ESCALATION_CHECKPOINT_SHA256,
            "conductor_timeout_seconds_per_candidate": args.conductor_timeout,
            "height_precisions": list(args.height_precisions),
            "height_timeout_seconds_per_candidate": args.height_timeout,
            "saturation_timeout_seconds_per_candidate": args.saturation_timeout,
            "certificate_prime_bound": args.certificate_prime_bound,
            "pari_stack_bytes": args.stack_bytes,
            "retries": 0,
            "temporary_directory_dependency": None,
            "output": str(args.output),
        },
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
            "sympy": sp.__version__,
        },
        "reproducing_command": command,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    print(
        f"preimages=22 record_sources=14 slices={len(slices)} "
        f"candidates={len(candidate_records)} below_target={len(below_target)} "
        f"multi_slice={len(multi_slice_intersections)} hits={len(target_hits)} "
        f"alternative_hits={len(alternative_rank_hits)}"
    )


if __name__ == "__main__":
    main()
