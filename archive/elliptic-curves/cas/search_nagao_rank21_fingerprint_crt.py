#!/usr/bin/env python3
"""Outside-box local-fingerprint CRT/Gauss search in Nagao's rank-21 family.

This experiment deliberately does not clone the residue vector of the known
``T=6793/64`` rank-19 specialization.  At six good primes it forms the full
union of finite symbols satisfying declared trace thresholds, combines them
with a deterministic bounded beam, and penalizes the exact shortest height in
the corresponding CRT lattice.  Independent strata add discriminant-power
conditions at 13, 37, and 83.  Gauss-reduced representatives are required to
lie outside the earlier positive ``a<=10000, b<=100`` rectangle and a pinned
set of known leads is excluded exactly.

Selection scores through 200 and 2000 omit every prime used to construct a
trace or discriminant-power condition.  Thus the score tail is leakage-free
with respect to the local fingerprint.  Exact homogenized-discriminant
valuations and a trial-division radical proxy gate the expensive searches.
Quartic points are mapped and decontaminated exactly before two-precision
height triage.  Stable numerical rank at least 18 triggers an immediate exact
finite-reduction checkpoint and conductor replay.  Numerical rank and bounded
search failures are never promoted to theorems.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from math import gcd, log
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
from typing import Any, Iterable, Sequence

from crt_lattice import crt_pair, gauss_reduce
from ek_k3 import primes_up_to, rational_to_string
from multiple_root_lifting import (
    affine_variable_coefficients,
    fixed_divisor_valuation,
)
from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK21_CONSTRUCTION,
    short_jacobian_coefficients,
)
from pari_bridge import minimal_curve_data, pari_version
from search_nagao_rank21_mutations import (
    conductor_radical_proxy,
    homogenized_discriminant,
)
from search_nagao_rank21_neighborhood import DISCRIMINANT_POLYNOMIAL
from search_nagao_rank21_unbiased import (
    PointPool,
    build_residue_tables,
    finite_reduction_certificate,
    parallel_point_search,
    parallel_rank_replay,
    point_pool_record,
    pool_priority,
    projective_index,
    residue_score,
    residue_table,
)


Q = Fraction
TARGET_LOG_CONDUCTOR = Decimal("182.72")
TRACE_THRESHOLDS = (
    (11, -4),
    (19, -5),
    (31, -8),
    (41, -10),
    (47, -9),
    (59, -9),
)
POWER_PRIMES = (13, 37, 83)
CHEAP_POWER_PRIMES = (5, 7, 23)
PROXY_PRIMES = (2, 5, 7, 13, 23, 37, 83)
# This is the exact nine-prime omission used by the supplied decontaminated
# calibration.  The cheap 5/7/23 strata are bad-reduction conditions and add
# no Nagao summand, so retaining those primes preserves comparability.
SCORE_OMITTED_PRIMES = (
    tuple(prime for prime, _ in TRACE_THRESHOLDS) + POWER_PRIMES
)
DESIGN_PRIMES = SCORE_OMITTED_PRIMES
TRACE_BEAM_WIDTH = 6_000
SINGLE_POWER_WIDTH = 1_200
DOUBLE_POWER_WIDTH = 500
TRIPLE_POWER_WIDTH = 200
RADICAL_PROXY_LIMIT = Decimal("190")
B200_RETAIN_COUNT = 800
EXACT_POINT_RETAIN_COUNT = 96
EXACT_SCORE_CUTOFF = 2_000
B200_TAIL_FLOOR = 17.0
CALIBRATION_PARAMETER = Q(6793, 64)
CALIBRATION_RESIDUES_AND_TRACES = (
    (11, 8, -5),
    (19, 15, -6),
    (31, 2, -9),
    (41, 3, -11),
    (47, 7, -12),
    (59, 37, -14),
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank21_fingerprint_crt.py"
)


# The first seven are previous deep/high-rank Nagao-family leads.  The next
# two are earlier outside-box mutation leaders.  The remainder are the smooth-
# denominator population supplied by the independently completed exhaustive
# stratum.  Equality is exact after Fraction normalization.
KNOWN_LEAD_PARAMETERS = frozenset(
    abs(parameter)
    for parameter in {
        Q(531, 2),
        Q(956, 9),
        Q(1637, 12),
        Q(5777, 32),
        Q(3137, 72),
        Q(5783, 16),
        Q(6793, 64),
        Q(-7591, 204),
        Q(-12847, 38),
        Q(18839, 800),
        Q(34927, 675),
        Q(11225, 1176),
        Q(16276, 225),
        Q(81265, 1024),
        Q(92941, 128),
        Q(52595, 336),
        Q(81815, 512),
        Q(50641, 108),
        Q(5483, 625),
        Q(5, 576),
        Q(86969, 1458),
        Q(94117, 480),
        Q(93827, 1458),
        Q(85541, 700),
        Q(39307, 300),
        Q(55059, 1000),
        Q(98717, 240),
        Q(22417, 1680),
        Q(3589, 512),
    }
)


SMOOTH_STRATUM_PARAMETERS = (
    Q(18839, 800),
    Q(34927, 675),
    Q(11225, 1176),
    Q(16276, 225),
    Q(81265, 1024),
    Q(92941, 128),
    Q(52595, 336),
    Q(81815, 512),
    Q(50641, 108),
    Q(5483, 625),
    Q(5, 576),
    Q(86969, 1458),
    Q(94117, 480),
    Q(93827, 1458),
    Q(85541, 700),
    Q(39307, 300),
    Q(55059, 1000),
    Q(98717, 240),
    Q(22417, 1680),
    Q(3589, 512),
)


@dataclass(frozen=True)
class TraceSymbol:
    prime: int
    residue: int
    ellap: int
    contribution: float


@dataclass(frozen=True)
class PowerSymbol:
    prime: int
    modulus: int
    residue: int
    forced_valuation: int
    condition_kind: str = "affine_parameter"

    @property
    def label(self) -> str:
        return f"p{self.prime}-m{self.modulus}-r{self.residue}-v{self.forced_valuation}"


@dataclass(frozen=True)
class BeamState:
    residue: int
    modulus: int
    trace_symbols: tuple[TraceSymbol, ...]
    power_symbols: tuple[PowerSymbol, ...]
    representative: Fraction
    representative_height: int
    reduced_basis: tuple[tuple[int, int], tuple[int, int]]
    objective: float


@dataclass(frozen=True)
class GeneratedCandidate:
    parameter: Fraction
    stratum: str
    trace_symbols: tuple[TraceSymbol, ...]
    power_symbols: tuple[PowerSymbol, ...]
    crt_residue: int
    crt_modulus: int
    gauss_reduced_basis: tuple[tuple[int, int], tuple[int, int]]
    gauss_coefficients: tuple[int, int]
    radical_proxy: dict[str, Any]

    @property
    def identifier(self) -> str:
        sign = "m" if self.parameter.numerator < 0 else "p"
        return (
            f"fingerprint-{sign}{abs(self.parameter.numerator)}-"
            f"{self.parameter.denominator}"
        )

    @property
    def height(self) -> int:
        return max(abs(self.parameter.numerator), self.parameter.denominator)


@dataclass(frozen=True)
class FingerprintPrefilter:
    generated: GeneratedCandidate
    residue_score_b200: float
    residue_good_primes: int
    residue_bad_primes: int

    @property
    def numerator(self) -> int:
        return self.generated.parameter.numerator

    @property
    def denominator(self) -> int:
        return self.generated.parameter.denominator

    @property
    def parameter(self) -> Fraction:
        return self.generated.parameter

    @property
    def identifier(self) -> str:
        return self.generated.identifier


@dataclass(frozen=True)
class FingerprintExact:
    prefilter: FingerprintPrefilter
    exact_score_b2000: str
    exact_good_primes: int
    exact_bad_primes: int
    exact_last_prime: int
    exact_design_primes_omitted: int

    @property
    def parameter(self) -> Fraction:
        return self.prefilter.parameter

    @property
    def identifier(self) -> str:
        return self.prefilter.identifier


def trace_symbol_union(prime: int, threshold: int) -> tuple[TraceSymbol, ...]:
    """Return every finite good symbol with ``a_p <= threshold``."""

    symbols = tuple(
        TraceSymbol(
            prime,
            symbol.projective_index,
            int(symbol.ellap),
            symbol.contribution,
        )
        for symbol in residue_table(prime)
        if symbol.projective_index < prime
        and symbol.good_reduction
        and symbol.ellap is not None
        and symbol.ellap <= threshold
    )
    if not symbols:
        raise AssertionError("a declared trace union is empty")
    return symbols


def calibration_verification() -> tuple[dict[str, Any], ...]:
    """Check, but never select on, the rank-19 calibration fingerprint."""

    records = []
    for prime, expected_residue, expected_trace in CALIBRATION_RESIDUES_AND_TRACES:
        residue = projective_index(
            CALIBRATION_PARAMETER.numerator,
            CALIBRATION_PARAMETER.denominator,
            prime,
        )
        symbol = residue_table(prime)[residue]
        if residue != expected_residue or symbol.ellap != expected_trace:
            raise AssertionError("the rank-19 calibration fingerprint changed")
        records.append(
            {
                "prime": prime,
                "residue": residue,
                "ellap": symbol.ellap,
                "used_for_generation": False,
            }
        )
    return tuple(records)


def gauss_representatives(
    residue: int,
    modulus: int,
    *,
    coefficient_radius: int = 2,
    limit: int = 2,
) -> tuple[
    tuple[
        Fraction,
        tuple[int, int],
        tuple[tuple[int, int], tuple[int, int]],
    ],
    ...,
]:
    """Return the shortest verified primitive representatives of one class."""

    if modulus < 1 or not 0 <= residue < modulus:
        raise ValueError("invalid CRT class")
    if coefficient_radius < 1 or limit < 1:
        raise ValueError("invalid Gauss representative bounds")
    basis = gauss_reduce((modulus, 0), (residue, 1))
    found: dict[Fraction, tuple[int, int]] = {}
    for left in range(-coefficient_radius, coefficient_radius + 1):
        for right in range(-coefficient_radius, coefficient_radius + 1):
            if left == 0 and right == 0:
                continue
            numerator = left * basis[0][0] + right * basis[1][0]
            denominator = left * basis[0][1] + right * basis[1][1]
            if denominator == 0:
                continue
            parameter = Q(numerator, denominator)
            if parameter == 0 or gcd(parameter.denominator, modulus) != 1:
                continue
            if (parameter.numerator - residue * parameter.denominator) % modulus:
                continue
            prior = found.get(parameter)
            coefficients = (left, right)
            if prior is None or coefficients < prior:
                found[parameter] = coefficients
    ordered = sorted(
        found.items(),
        key=lambda item: (
            max(abs(item[0].numerator), item[0].denominator),
            abs(item[0].numerator),
            item[0].denominator,
            item[0].numerator,
            item[1],
        ),
    )
    if not ordered:
        raise AssertionError("a CRT lattice had no local-unit representative")
    return tuple((parameter, pair, basis) for parameter, pair in ordered[:limit])


def state_representatives(
    state: BeamState,
    *,
    coefficient_radius: int = 4,
    limit: int = 3,
) -> tuple[
    tuple[
        Fraction,
        tuple[int, int],
        tuple[tuple[int, int], tuple[int, int]],
    ],
    ...,
]:
    """Search a stored reduced basis, including the projective 2-adic one."""

    found: dict[Fraction, tuple[int, int]] = {}
    basis = state.reduced_basis
    for left in range(-coefficient_radius, coefficient_radius + 1):
        for right in range(-coefficient_radius, coefficient_radius + 1):
            if left == 0 and right == 0:
                continue
            numerator = left * basis[0][0] + right * basis[1][0]
            denominator = left * basis[0][1] + right * basis[1][1]
            if denominator == 0:
                continue
            parameter = Q(numerator, denominator)
            if parameter == 0 or not _satisfies_state(parameter, state):
                continue
            pair = (left, right)
            if parameter not in found or pair < found[parameter]:
                found[parameter] = pair
    ordered = sorted(
        found.items(),
        key=lambda item: (
            max(abs(item[0].numerator), item[0].denominator),
            abs(item[0].numerator),
            item[0].denominator,
            item[0].numerator,
            item[1],
        ),
    )
    if not ordered:
        raise AssertionError("a reduced beam lattice had no valid representative")
    return tuple((parameter, pair, basis) for parameter, pair in ordered[:limit])


def _state(
    residue: int,
    modulus: int,
    trace_symbols: tuple[TraceSymbol, ...],
    power_symbols: tuple[PowerSymbol, ...],
) -> BeamState:
    parameter, _, basis = gauss_representatives(residue, modulus, limit=1)[0]
    height = max(abs(parameter.numerator), parameter.denominator)
    trace_reward = sum(symbol.contribution for symbol in trace_symbols)
    power_savings = sum(
        (symbol.forced_valuation - 1) * log(symbol.prime)
        for symbol in power_symbols
    )
    # The discriminant is degree 20 in the primitive pair.  This objective is
    # a deliberately simple proxy, used only inside each bounded beam.
    objective = 20 * log(height) - power_savings - 2 * trace_reward
    return BeamState(
        residue,
        modulus,
        trace_symbols,
        power_symbols,
        parameter,
        height,
        basis,
        objective,
    )


def _retain_beam(states: Iterable[BeamState], width: int) -> tuple[BeamState, ...]:
    if width < 1:
        raise ValueError("beam width must be positive")
    deduplicated: dict[tuple[int, int], BeamState] = {}
    for state in states:
        key = state.residue, state.modulus
        previous = deduplicated.get(key)
        ordering = (
            state.objective,
            state.representative_height,
            state.residue,
            tuple((s.prime, s.residue) for s in state.trace_symbols),
            tuple((s.prime, s.modulus, s.residue) for s in state.power_symbols),
        )
        if previous is None:
            deduplicated[key] = state
        else:
            previous_ordering = (
                previous.objective,
                previous.representative_height,
                previous.residue,
                tuple((s.prime, s.residue) for s in previous.trace_symbols),
                tuple(
                    (s.prime, s.modulus, s.residue)
                    for s in previous.power_symbols
                ),
            )
            if ordering < previous_ordering:
                deduplicated[key] = state
    return tuple(
        sorted(
            deduplicated.values(),
            key=lambda state: (
                state.objective,
                state.representative_height,
                state.modulus,
                state.residue,
            ),
        )[:width]
    )


def build_trace_beam(
    *, width: int = TRACE_BEAM_WIDTH
) -> tuple[tuple[BeamState, ...], tuple[dict[str, Any], ...]]:
    states = (_state(0, 1, (), ()),)
    audit = []
    for prime, threshold in TRACE_THRESHOLDS:
        symbols = trace_symbol_union(prime, threshold)
        expanded = []
        for state in states:
            for symbol in symbols:
                residue, modulus = crt_pair(
                    state.residue,
                    state.modulus,
                    symbol.residue,
                    prime,
                )
                expanded.append(
                    _state(
                        residue,
                        modulus,
                        state.trace_symbols + (symbol,),
                        state.power_symbols,
                    )
                )
        states = _retain_beam(expanded, width)
        audit.append(
            {
                "prime": prime,
                "trace_threshold": threshold,
                "symbol_union_size": len(symbols),
                "expanded_state_count": len(expanded),
                "retained_state_count": len(states),
                "minimum_representative_height": min(
                    state.representative_height for state in states
                ),
                "maximum_retained_representative_height": max(
                    state.representative_height for state in states
                ),
            }
        )
    return states, tuple(audit)


def power_symbol_union(prime: int) -> tuple[PowerSymbol, ...]:
    """Return base multiple-root balls and exceptional stronger p^2 lifts."""

    if prime not in POWER_PRIMES:
        raise ValueError("undeclared power prime")
    options = []
    for residue in range(prime):
        base_valuation = fixed_divisor_valuation(
            affine_variable_coefficients(
                DISCRIMINANT_POLYNOMIAL, residue, prime
            ),
            prime,
        )
        if base_valuation < 2:
            continue
        options.append(PowerSymbol(prime, prime, residue, base_valuation))
        square = prime * prime
        for lift_index in range(prime):
            lifted_residue = residue + lift_index * prime
            lifted_valuation = fixed_divisor_valuation(
                affine_variable_coefficients(
                    DISCRIMINANT_POLYNOMIAL, lifted_residue, square
                ),
                prime,
            )
            if lifted_valuation > base_valuation:
                options.append(
                    PowerSymbol(
                        prime,
                        square,
                        lifted_residue,
                        lifted_valuation,
                    )
                )
    options.sort(
        key=lambda option: (
            option.modulus,
            option.residue,
            -option.forced_valuation,
        )
    )
    return tuple(options)


def cheap_power_symbol_union(prime: int) -> tuple[PowerSymbol, ...]:
    """Return automatic alternative unions at common/fixed bad primes."""

    pinned = {
        5: (5, (1, 2, 3, 4)),
        7: (
            343,
            (4, 59, 60, 66, 91, 161, 182, 252, 277, 283, 284, 339),
        ),
        23: (23, (2, 3, 4, 6, 8, 11, 12, 15, 17, 19, 20, 21)),
    }
    if prime not in pinned:
        raise ValueError("undeclared cheap power prime")
    modulus, residues = pinned[prime]
    answer = []
    for residue in residues:
        valuation = fixed_divisor_valuation(
            affine_variable_coefficients(
                DISCRIMINANT_POLYNOMIAL, residue, modulus
            ),
            prime,
        )
        if valuation < {5: 3, 7: 6, 23: 3}[prime]:
            raise AssertionError("an automatic cheap power ball changed")
        answer.append(PowerSymbol(prime, modulus, residue, valuation))
    if prime == 7:
        for residue in (1836, 4149, 12658, 14971):
            valuation = fixed_divisor_valuation(
                affine_variable_coefficients(
                    DISCRIMINANT_POLYNOMIAL, residue, 16807
                ),
                7,
            )
            if valuation < 6:
                raise AssertionError("an exceptional 7-adic ball changed")
            answer.append(PowerSymbol(7, 16807, residue, valuation))
    return tuple(answer)


def extend_power_beam(
    states: Sequence[BeamState],
    options: Sequence[PowerSymbol],
    *,
    width: int,
) -> tuple[BeamState, ...]:
    expanded = []
    for state in states:
        for option in options:
            if gcd(state.modulus, option.modulus) != 1:
                raise AssertionError("a power-prime stratum repeated a modulus")
            residue, modulus = crt_pair(
                state.residue,
                state.modulus,
                option.residue,
                option.modulus,
            )
            expanded.append(
                _state(
                    residue,
                    modulus,
                    state.trace_symbols,
                    state.power_symbols + (option,),
                )
            )
    return _retain_beam(expanded, width)


def build_even_denominator_beam(
    trace_states: Sequence[BeamState], *, width: int
) -> tuple[BeamState, ...]:
    """Intersect each odd trace lattice with odd ``a`` and ``16 | b``."""

    expanded = []
    symbol = PowerSymbol(2, 16, 0, 22, "denominator_divisibility")
    for state in trace_states:
        # a=M*k+16*r*c, b=16*c.  The odd-a primitive vectors are precisely
        # the desired projective 2-adic ball because the trace modulus is odd.
        basis = gauss_reduce(
            (state.modulus, 0), (16 * state.residue, 16)
        )
        provisional = BeamState(
            state.residue,
            state.modulus,
            state.trace_symbols,
            state.power_symbols + (symbol,),
            Q(1),
            1,
            basis,
            0.0,
        )
        try:
            parameter, _, _ = state_representatives(
                provisional, coefficient_radius=4, limit=1
            )[0]
        except AssertionError:
            continue
        height = max(abs(parameter.numerator), parameter.denominator)
        trace_reward = sum(item.contribution for item in state.trace_symbols)
        savings = 21 * log(2)
        expanded.append(
            BeamState(
                state.residue,
                state.modulus,
                state.trace_symbols,
                state.power_symbols + (symbol,),
                parameter,
                height,
                basis,
                20 * log(height) - savings - 2 * trace_reward,
            )
        )
    return _retain_beam(expanded, width)


def build_power_strata(
    trace_states: Sequence[BeamState],
    *,
    single_width: int = SINGLE_POWER_WIDTH,
    double_width: int = DOUBLE_POWER_WIDTH,
    triple_width: int = TRIPLE_POWER_WIDTH,
) -> dict[str, tuple[BeamState, ...]]:
    # ``double_width`` and ``triple_width`` remain explicit API parameters so
    # old small-bound test calls fail deterministically if the interface moves;
    # their exact-shortest branches were independently exhausted and are not
    # regenerated in this complementary lane.
    if double_width < 1 or triple_width < 1:
        raise ValueError("closed-branch widths must remain positive")
    options = {prime: power_symbol_union(prime) for prime in POWER_PRIMES}
    strata: dict[str, tuple[BeamState, ...]] = {"trace-only": tuple(trace_states)}
    for prime in POWER_PRIMES:
        strata[f"power-{prime}"] = extend_power_beam(
            trace_states, options[prime], width=single_width
        )
    for prime in CHEAP_POWER_PRIMES:
        strata[f"cheap-power-{prime}"] = extend_power_beam(
            trace_states,
            cheap_power_symbol_union(prime),
            width=single_width,
        )
    strata["even-denominator-v2-22"] = build_even_denominator_beam(
        trace_states, width=single_width
    )
    return strata


def in_old_positive_rectangle(parameter: Fraction) -> bool:
    parameter = Q(parameter)
    return (
        abs(parameter.numerator) <= 10_000
        and parameter.denominator <= 100
    )


def canonicalize_even_parameter(
    parameter: Fraction,
    state: BeamState,
    basis: tuple[tuple[int, int], tuple[int, int]],
) -> tuple[
    Fraction,
    tuple[TraceSymbol, ...],
    tuple[PowerSymbol, ...],
    int,
    tuple[tuple[int, int], tuple[int, int]],
]:
    """Canonicalize ``T`` to ``|T|`` and transport every local condition."""

    parameter = Q(parameter)
    if parameter > 0:
        return (
            parameter,
            state.trace_symbols,
            state.power_symbols,
            state.residue,
            basis,
        )
    traces = tuple(
        TraceSymbol(
            symbol.prime,
            (-symbol.residue) % symbol.prime,
            symbol.ellap,
            symbol.contribution,
        )
        for symbol in state.trace_symbols
    )
    powers = tuple(
        PowerSymbol(
            symbol.prime,
            symbol.modulus,
            (
                symbol.residue
                if symbol.condition_kind == "denominator_divisibility"
                else (-symbol.residue) % symbol.modulus
            ),
            symbol.forced_valuation,
            symbol.condition_kind,
        )
        for symbol in state.power_symbols
    )
    transported_basis = tuple((-x_value, y_value) for x_value, y_value in basis)
    return -parameter, traces, powers, (-state.residue) % state.modulus, transported_basis


def _satisfies_state(parameter: Fraction, state: BeamState) -> bool:
    for symbol in state.trace_symbols:
        if (
            projective_index(
                parameter.numerator, parameter.denominator, symbol.prime
            )
            != symbol.residue
        ):
            return False
    for symbol in state.power_symbols:
        if symbol.condition_kind == "denominator_divisibility":
            if (
                parameter.numerator % 2 != 1
                or parameter.denominator % symbol.modulus != 0
            ):
                return False
        elif not (
            (parameter.numerator - symbol.residue * parameter.denominator)
            % symbol.modulus
            == 0
            and gcd(parameter.denominator, symbol.modulus) == 1
        ):
            return False
    return True


def generate_candidates(
    strata: dict[str, tuple[BeamState, ...]],
    *,
    proxy_limit: Decimal = RADICAL_PROXY_LIMIT,
) -> tuple[tuple[GeneratedCandidate, ...], dict[str, Any]]:
    retained: dict[Fraction, GeneratedCandidate] = {}
    raw_representatives = 0
    excluded_old_box = 0
    excluded_known = 0
    excluded_proxy = 0
    singular = 0
    for stratum, states in strata.items():
        for state in states:
            shell = state_representatives(
                state,
                coefficient_radius=4,
                limit=12,
            )
            requested_indices = (
                (0, 2, 5)
                if stratum in ("trace-only", "even-denominator-v2-22")
                else (0, 3)
            )
            representatives = tuple(
                shell[index] for index in requested_indices if index < len(shell)
            )
            for parameter, pair, basis in representatives:
                raw_representatives += 1
                if not _satisfies_state(parameter, state):
                    raise AssertionError("primitive Gauss normalization lost a symbol")
                (
                    parameter,
                    trace_symbols,
                    power_symbols,
                    crt_residue,
                    basis,
                ) = canonicalize_even_parameter(parameter, state, basis)
                if in_old_positive_rectangle(parameter):
                    excluded_old_box += 1
                    continue
                if parameter in KNOWN_LEAD_PARAMETERS:
                    excluded_known += 1
                    continue
                try:
                    proxy = conductor_radical_proxy(parameter)
                except ValueError:
                    singular += 1
                    continue
                if Decimal(str(proxy["log_radical_upper_proxy"])) >= proxy_limit:
                    excluded_proxy += 1
                    continue
                candidate = GeneratedCandidate(
                    parameter,
                    stratum,
                    trace_symbols,
                    power_symbols,
                    crt_residue,
                    state.modulus,
                    basis,
                    pair,
                    proxy,
                )
                previous = retained.get(parameter)
                candidate_key = (
                    -sum(
                        symbol.forced_valuation - 1
                        for symbol in candidate.power_symbols
                    ),
                    candidate.radical_proxy["log_radical_upper_proxy"],
                    candidate.stratum,
                )
                if previous is None:
                    retained[parameter] = candidate
                else:
                    previous_key = (
                        -sum(
                            symbol.forced_valuation - 1
                            for symbol in previous.power_symbols
                        ),
                        previous.radical_proxy["log_radical_upper_proxy"],
                        previous.stratum,
                    )
                    if candidate_key < previous_key:
                        retained[parameter] = candidate
    candidates = tuple(
        sorted(
            retained.values(),
            key=lambda candidate: (
                candidate.radical_proxy["log_radical_upper_proxy"],
                candidate.height,
                candidate.identifier,
            ),
        )
    )
    digest = hashlib.sha256()
    for candidate in candidates:
        digest.update(
            (
                f"{candidate.parameter}|{candidate.stratum}|"
                f"{candidate.radical_proxy['log_radical_upper_proxy']!r}\n"
            ).encode()
        )
    stratum_survivors = {
        stratum: sum(candidate.stratum == stratum for candidate in candidates)
        for stratum in strata
    }
    return candidates, {
        "raw_gauss_representatives": raw_representatives,
        "excluded_by_old_positive_rectangle": excluded_old_box,
        "excluded_as_known_leads": excluded_known,
        "excluded_by_radical_proxy": excluded_proxy,
        "singular_specializations": singular,
        "exactly_deduplicated_proxy_survivors": len(candidates),
        "stratum_selected_provenance_counts": stratum_survivors,
        "survivor_stream_sha256": digest.hexdigest(),
    }


def leakage_free_tables(cutoff: int = 200) -> dict[int, Any]:
    tables = build_residue_tables(cutoff)
    return {prime: table for prime, table in tables.items() if prime not in DESIGN_PRIMES}


def b200_prefilter(
    candidates: Sequence[GeneratedCandidate],
    *,
    keep_count: int,
    tables: dict[int, Any],
) -> tuple[FingerprintPrefilter, ...]:
    scored = []
    for candidate in candidates:
        score, good, bad = residue_score(
            candidate.parameter.numerator,
            candidate.parameter.denominator,
            tables,
        )
        scored.append(FingerprintPrefilter(candidate, score, good, bad))
    scored.sort(
        key=lambda candidate: (
            -candidate.residue_score_b200,
            candidate.generated.radical_proxy["log_radical_upper_proxy"],
            candidate.generated.height,
            candidate.identifier,
        )
    )
    return tuple(scored[: min(keep_count, len(scored))])


def _gp_rational(value: Fraction) -> str:
    return f"({value.numerator}/{value.denominator})"


def exact_decontaminated_scores(
    candidates: Sequence[FingerprintPrefilter],
    *,
    cutoff: int,
    batch_size: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[FingerprintExact, ...]:
    """Use PARI for an exact B score with every design prime omitted."""

    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    if batch_size < 1 or cutoff < 5 or timeout <= 0 or stack_bytes < 8_000_000:
        raise ValueError("invalid exact-score bounds")
    omitted_test = "||".join(f"p=={prime}" for prime in DESIGN_PRIMES)
    last_prime = primes_up_to(cutoff)[-1]
    records: dict[str, dict[str, Any]] = {}
    for start in range(0, len(candidates), batch_size):
        batch = candidates[start : start + batch_size]
        commands = ["default(realprecision,80);"]
        for index, candidate in enumerate(batch):
            coefficients = short_jacobian_coefficients(
                RANK21_CONSTRUCTION, candidate.parameter
            )
            vector = ",".join(_gp_rational(value) for value in coefficients)
            commands.extend(
                (
                    f"E=ellminimalmodel(ellinit([{vector}]));",
                    "S=0;USED=0;BAD=0;OMITTED=0;",
                    (
                        f"forprime(p=5,{cutoff},if({omitted_test},OMITTED++,"
                        "if(valuation(E.disc,p)>0,BAD++,"
                        "A=ellap(E,p);S+=(2-A)/(p+1-A)*log(p);USED++)));"
                    ),
                    (
                        f'print("ROW|{index}|",S,"|",USED,"|",BAD,"|",'
                        f'OMITTED,"|{last_prime}");'
                    ),
                )
            )
        commands.append("quit")
        result = subprocess.run(
            [executable, "-q", "-s", str(stack_bytes)],
            input="\n".join(commands) + "\n",
            text=True,
            capture_output=True,
            timeout=timeout,
        )
        if result.returncode != 0 or "***" in result.stderr:
            raise RuntimeError(f"PARI/GP score failed: {result.stderr.strip()}")
        observed = 0
        for line in result.stdout.splitlines():
            if not line.startswith("ROW|"):
                continue
            _, index_text, score, used, bad, omitted, observed_last = line.split("|")
            candidate = batch[int(index_text)]
            records[candidate.identifier] = {
                "score": score,
                "used": int(used),
                "bad": int(bad),
                "omitted": int(omitted),
                "last": int(observed_last),
            }
            observed += 1
        if observed != len(batch):
            raise RuntimeError("PARI omitted one or more exact-score records")
    answer = tuple(
        FingerprintExact(
            candidate,
            records[candidate.identifier]["score"],
            records[candidate.identifier]["used"],
            records[candidate.identifier]["bad"],
            records[candidate.identifier]["last"],
            records[candidate.identifier]["omitted"],
        )
        for candidate in candidates
    )
    return tuple(
        sorted(
            answer,
            key=lambda candidate: (
                -Decimal(candidate.exact_score_b2000),
                candidate.prefilter.generated.radical_proxy[
                    "log_radical_upper_proxy"
                ],
                candidate.prefilter.generated.height,
                candidate.identifier,
            ),
        )
    )


def select_point_survivors(
    candidates: Sequence[FingerprintExact],
    *,
    keep_count: int,
    b200_floor: float = B200_TAIL_FLOOR,
) -> tuple[FingerprintExact, ...]:
    """Preserve the declared B200 tail, then fill by exact B2000 rank."""

    mandatory = {
        candidate.identifier
        for candidate in candidates
        if candidate.prefilter.residue_score_b200 >= b200_floor
    }
    if len(mandatory) > keep_count:
        raise ValueError("the declared B200 tail exceeds the point budget")
    selected = [candidate for candidate in candidates if candidate.identifier in mandatory]
    selected_ids = set(mandatory)
    for candidate in candidates:
        if len(selected) == keep_count:
            break
        if candidate.identifier not in selected_ids:
            selected.append(candidate)
            selected_ids.add(candidate.identifier)
    return tuple(selected)


def candidate_record(candidate: FingerprintExact) -> dict[str, Any]:
    generated = candidate.prefilter.generated
    return {
        "candidate_id": candidate.identifier,
        "constructor_parameter": rational_to_string(candidate.parameter),
        "height": generated.height,
        "stratum": generated.stratum,
        "trace_symbols": [symbol.__dict__ for symbol in generated.trace_symbols],
        "power_symbols": [
            {
                "prime": symbol.prime,
                "modulus": symbol.modulus,
                "residue": symbol.residue,
                "forced_homogenized_discriminant_valuation": (
                    symbol.forced_valuation
                ),
            }
            for symbol in generated.power_symbols
        ],
        "crt_residue": generated.crt_residue,
        "crt_modulus": generated.crt_modulus,
        "gauss_reduced_basis": [list(vector) for vector in generated.gauss_reduced_basis],
        "gauss_coefficients": list(generated.gauss_coefficients),
        "radical_proxy": generated.radical_proxy,
        "exact_homogenized_discriminant_valuations_2_5_7_13_23_37_83": {
            str(prime): integer_valuation(
                homogenized_discriminant(generated.parameter), prime
            )
            for prime in PROXY_PRIMES
        },
        "leakage_free_b200": {
            "score": candidate.prefilter.residue_score_b200,
            "good_primes_used": candidate.prefilter.residue_good_primes,
            "bad_primes_skipped": candidate.prefilter.residue_bad_primes,
            "design_primes_omitted": list(DESIGN_PRIMES),
        },
        "exact_pari_leakage_free_b2000": {
            "score": candidate.exact_score_b2000,
            "good_primes_used": candidate.exact_good_primes,
            "bad_primes_skipped": candidate.exact_bad_primes,
            "design_primes_omitted_count": candidate.exact_design_primes_omitted,
            "last_numerical_prime": candidate.exact_last_prime,
        },
    }


def integer_valuation(value: int, prime: int) -> int:
    if value == 0:
        raise ValueError("zero has no finite valuation")
    value = abs(value)
    valuation = 0
    while value % prime == 0:
        value //= prime
        valuation += 1
    return valuation


def smooth_stratum_closure() -> tuple[dict[str, Any], ...]:
    records = []
    for parameter in SMOOTH_STRATUM_PARAMETERS:
        proxy = conductor_radical_proxy(parameter)
        records.append(
            {
                "constructor_parameter": rational_to_string(parameter),
                "radical_proxy": proxy,
                "point_search_in_this_lane": False,
                "reason": (
                    "independently exhausted smooth-denominator stratum; its "
                    "best exact proxy/conductor cases are above target"
                ),
            }
        )
    minimum = min(record["radical_proxy"]["log_radical_upper_proxy"] for record in records)
    if abs(minimum - 191.86811185320965) > 1e-10:
        raise AssertionError("the supplied smooth-stratum proxy frontier changed")
    return tuple(records)


def parse_positive_ints(value: str) -> tuple[int, ...]:
    try:
        result = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from error
    if not result or any(item < 1 for item in result):
        raise argparse.ArgumentTypeError("all integers must be positive")
    return result


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace-beam-width", type=int, default=TRACE_BEAM_WIDTH)
    parser.add_argument("--b200-keep", type=int, default=B200_RETAIN_COUNT)
    parser.add_argument("--exact-point-keep", type=int, default=EXACT_POINT_RETAIN_COUNT)
    parser.add_argument("--exact-score-batch", type=int, default=200)
    parser.add_argument("--exact-score-timeout", type=float, default=60.0)
    parser.add_argument("--stage-heights", type=parse_positive_ints, default=(50_000, 250_000, 1_000_000))
    parser.add_argument("--stage-keeps", type=parse_positive_ints, default=(20, 5))
    parser.add_argument("--stage-timeouts", type=parse_positive_ints, default=(5, 20, 90))
    parser.add_argument("--stage-workers", type=parse_positive_ints, default=(4, 4, 2))
    parser.add_argument("--height-precisions", type=parse_positive_ints, default=(72, 120))
    parser.add_argument("--height-timeout", type=float, default=20.0)
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--conductor-timeout", type=float, default=30.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=1_000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts/generated-results/elliptic_nagao_rank21_fingerprint_crt.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.trace_beam_width < 100:
        raise SystemExit("--trace-beam-width must be at least 100")
    if len(args.stage_heights) != 3 or len(args.stage_keeps) != 2:
        raise SystemExit("this run requires three point heights and two keep counts")
    if len(args.stage_timeouts) != 3 or len(args.stage_workers) != 3:
        raise SystemExit("provide one timeout and worker count per point stage")
    if tuple(sorted(set(args.stage_heights))) != args.stage_heights:
        raise SystemExit("point heights must be strictly increasing")
    if args.stage_keeps[1] > args.stage_keeps[0]:
        raise SystemExit("point keep counts must be nonincreasing")
    if args.b200_keep < args.exact_point_keep or args.stage_keeps[0] > args.exact_point_keep:
        raise SystemExit("selection keep counts must be nonincreasing")
    if min(
        args.exact_score_timeout,
        *args.stage_timeouts,
        args.height_timeout,
        args.saturation_timeout,
        args.conductor_timeout,
    ) <= 0:
        raise SystemExit("all subprocess timeouts must be positive")
    if args.stack_bytes < 64_000_000 or args.certificate_prime_bound < 3:
        raise SystemExit("invalid stack or finite-reduction prime bound")

    calibration = calibration_verification()
    symbol_unions = {
        prime: trace_symbol_union(prime, threshold)
        for prime, threshold in TRACE_THRESHOLDS
    }
    trace_states, trace_audit = build_trace_beam(width=args.trace_beam_width)
    power_strata = build_power_strata(trace_states)
    candidates, generation_audit = generate_candidates(power_strata)
    if not candidates:
        raise SystemExit("the exact radical gate removed the entire population")
    print(
        f"generated proxy survivors={len(candidates)} "
        f"from raw={generation_audit['raw_gauss_representatives']}",
        flush=True,
    )

    tables = leakage_free_tables(200)
    prefiltered = b200_prefilter(
        candidates, keep_count=args.b200_keep, tables=tables
    )
    exactly_scored = exact_decontaminated_scores(
        prefiltered,
        cutoff=EXACT_SCORE_CUTOFF,
        batch_size=args.exact_score_batch,
        timeout=args.exact_score_timeout,
        stack_bytes=args.stack_bytes,
    )
    exact_survivors = select_point_survivors(
        exactly_scored, keep_count=args.exact_point_keep
    )

    calibration_generated = GeneratedCandidate(
        CALIBRATION_PARAMETER,
        "calibration-only",
        tuple(
            TraceSymbol(prime, residue, trace, 0.0)
            for prime, residue, trace in CALIBRATION_RESIDUES_AND_TRACES
        ),
        (),
        0,
        1,
        ((1, 0), (0, 1)),
        (0, 1),
        conductor_radical_proxy(CALIBRATION_PARAMETER),
    )
    calibration_prefilter = FingerprintPrefilter(
        calibration_generated,
        *residue_score(
            CALIBRATION_PARAMETER.numerator,
            CALIBRATION_PARAMETER.denominator,
            tables,
        ),
    )
    calibration_exact = exact_decontaminated_scores(
        (calibration_prefilter,),
        cutoff=EXACT_SCORE_CUTOFF,
        batch_size=1,
        timeout=args.exact_score_timeout,
        stack_bytes=args.stack_bytes,
    )[0]

    point_stages = []
    current = tuple(exact_survivors)
    checkpoint_targets: dict[str, tuple[PointPool, dict[str, Any]]] = {}
    exact_checkpoints: dict[str, dict[str, Any]] = {}
    conductor_replays: dict[str, dict[str, Any]] = {}
    final_pools: tuple[PointPool, ...] = ()
    final_ranks: dict[str, dict[str, Any]] = {}

    def checkpoint(
        pool: PointPool, rank_record: dict[str, Any], stage_height: int
    ) -> None:
        identifier = pool.candidate.identifier
        stable_rank = int(rank_record["stable_numerical_rank"])
        if stable_rank < 18:
            return
        prior = checkpoint_targets.get(identifier)
        if prior is not None and prior[0].height_bound >= pool.height_bound:
            return
        checkpoint_targets[identifier] = (pool, rank_record)
        try:
            certificate = finite_reduction_certificate(
                pool,
                rank_record,
                saturation_timeout=args.saturation_timeout,
                certificate_prime_bound=args.certificate_prime_bound,
                stack_bytes=args.stack_bytes,
            )
        except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
            certificate = {
                "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
                "error": str(error)[:500],
            }
        certificate["trigger_stage_height"] = stage_height
        certificate["trigger_stable_numerical_rank"] = stable_rank
        exact_checkpoints[identifier] = certificate
        try:
            conductor = {
                "status": "completed",
                **minimal_curve_data(
                    pool.coefficients,
                    timeout=args.conductor_timeout,
                    local_primes=POWER_PRIMES,
                    stack_bytes=args.stack_bytes,
                ),
            }
            conductor["below_strict_log_conductor_target"] = (
                Decimal(conductor["log_conductor"]) < TARGET_LOG_CONDUCTOR
            )
        except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
            conductor = {
                "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
                "error": str(error)[:500],
            }
        conductor["trigger_stage_height"] = stage_height
        conductor_replays[identifier] = conductor

    for stage_index, (height_bound, timeout, workers) in enumerate(
        zip(args.stage_heights, args.stage_timeouts, args.stage_workers), start=1
    ):
        pools = parallel_point_search(
            current,
            height_bound=height_bound,
            timeout=float(timeout),
            stack_bytes=args.stack_bytes,
            workers=workers,
        )
        ranks = parallel_rank_replay(
            pools,
            precisions=args.height_precisions,
            timeout=args.height_timeout,
            stack_bytes=args.stack_bytes,
            workers=workers,
        )
        ordered = tuple(
            sorted(
                pools,
                key=lambda pool: pool_priority(
                    pool, ranks[pool.candidate.identifier]
                ),
            )
        )
        for pool in ordered:
            record = ranks[pool.candidate.identifier]
            if record.get("status") == "completed":
                checkpoint(pool, record, height_bound)
        keep_count = (
            args.stage_keeps[stage_index - 1]
            if stage_index <= len(args.stage_keeps)
            else len(ordered)
        )
        retained = ordered[: min(keep_count, len(ordered))]
        point_stages.append(
            {
                "stage": stage_index,
                "quartic_naive_height_bound": height_bound,
                "population_searched": len(pools),
                "completed_point_searches": sum(pool.status == "completed" for pool in pools),
                "point_search_timeouts": sum(pool.status == "timeout" for pool in pools),
                "point_search_errors": sum(pool.status == "error" for pool in pools),
                "ranked_population": [
                    {
                        **point_pool_record(
                            pool,
                            ranks[pool.candidate.identifier],
                            include_points=(stage_index == 3),
                        ),
                        "fingerprint_generation": candidate_record(pool.candidate),
                    }
                    for pool in ordered
                ],
                "retained_candidate_ids": [pool.candidate.identifier for pool in retained],
            }
        )
        current = tuple(pool.candidate for pool in retained)
        final_pools, final_ranks = ordered, ranks
        print(
            f"point stage H={height_bound} population={len(pools)} "
            f"best_rank={max((int(r.get('stable_numerical_rank',0)) for r in ranks.values()),default=0)}",
            flush=True,
        )

    # Final leaders receive conductor replays even below the exact checkpoint
    # threshold, so a negative result still records the conductor frontier.
    for pool in final_pools[: min(5, len(final_pools))]:
        identifier = pool.candidate.identifier
        if identifier in conductor_replays:
            continue
        try:
            conductor = {
                "status": "completed",
                **minimal_curve_data(
                    pool.coefficients,
                    timeout=args.conductor_timeout,
                    local_primes=POWER_PRIMES,
                    stack_bytes=args.stack_bytes,
                ),
            }
            conductor["below_strict_log_conductor_target"] = (
                Decimal(conductor["log_conductor"]) < TARGET_LOG_CONDUCTOR
            )
        except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
            conductor = {
                "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
                "error": str(error)[:500],
            }
        conductor_replays[identifier] = conductor

    certified_hits = []
    for identifier, certificate in exact_checkpoints.items():
        conductor = conductor_replays.get(identifier, {})
        rank = certificate.get("certified_algebraic_rank_lower_bound")
        if certificate.get("status") != "certified" or rank is None:
            continue
        if int(rank) >= 30 or (
            int(rank) >= 21
            and conductor.get("below_strict_log_conductor_target") is True
        ):
            pool, _ = checkpoint_targets[identifier]
            certified_hits.append(
                {
                    "candidate_id": identifier,
                    "constructor_parameter": rational_to_string(pool.candidate.parameter),
                    "certified_rank_lower_bound": rank,
                    "conductor": conductor.get("conductor"),
                    "log_conductor": conductor.get("log_conductor"),
                }
            )

    smooth_closure = smooth_stratum_closure()
    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded outside-box fingerprint CRT/Gauss search; numerical "
            "height rank is triage only"
        ),
        "primary_source": PRIMARY_SOURCE,
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "certified_hits": certified_hits,
        },
        "local_fingerprint": {
            "selection_rule": "union of all finite good symbols with a_p at or below the threshold",
            "thresholds": {str(prime): threshold for prime, threshold in TRACE_THRESHOLDS},
            "symbol_unions": {
                str(prime): [symbol.__dict__ for symbol in symbols]
                for prime, symbols in symbol_unions.items()
            },
            "rank19_calibration_only": list(calibration),
            "calibration_parameter": rational_to_string(CALIBRATION_PARAMETER),
            "calibration_never_inserted_into_generation_or_selection": True,
        },
        "beam_generation": {
            "trace_beam_width": args.trace_beam_width,
            "objective": "20*log(shortest_Gauss_height)-power_savings-2*trace_reward",
            "trace_stages": list(trace_audit),
            "power_symbol_unions": {
                str(prime): [symbol.__dict__ for symbol in power_symbol_union(prime)]
                for prime in POWER_PRIMES
            },
            "cheap_power_symbol_unions": {
                str(prime): [
                    symbol.__dict__ for symbol in cheap_power_symbol_union(prime)
                ]
                for prime in CHEAP_POWER_PRIMES
            },
            "stratum_state_counts": {
                stratum: len(states) for stratum, states in power_strata.items()
            },
            "single_power_width": SINGLE_POWER_WIDTH,
            "closed_double_and_triple_power_widths_not_used": {
                "double": DOUBLE_POWER_WIDTH,
                "triple": TRIPLE_POWER_WIDTH,
            },
            "old_positive_rectangle_excluded": {
                "numerator_at_most": 10_000,
                "denominator_at_most": 100,
            },
            "known_leads_excluded": [
                rational_to_string(parameter)
                for parameter in sorted(KNOWN_LEAD_PARAMETERS)
            ],
            "radical_proxy_strict_upper_limit": str(RADICAL_PROXY_LIMIT),
            **generation_audit,
        },
        "closed_smooth_denominator_stratum": {
            "selection_role": "recorded exclusion; never point-searched here",
            "records": list(smooth_closure),
            "independent_exact_conductor_frontier": {
                "3589/512": {"log_conductor": "190.7695", "root_number": -1},
                "93827/1458": {"log_conductor": "204.6691", "root_number": -1},
                "86969/1458": {"status": "strict-20-second-timeout"},
            },
        },
        "closed_exact_shortest_multi_power_branches": {
            "selection_role": (
                "independently exhausted; this lane uses single/cheap/2-adic "
                "strata and non-short coefficient shells"
            ),
            "full_13_37_83": {
                "outside_box_gauss_representatives": 1638,
                "proxy_below_190_and_all_six_trace_thresholds": 3,
                "h50000_results": [
                    {"constructor_parameter": "2279/138", "stable_numerical_rank": 12},
                    {"constructor_parameter": "2621/218", "stable_numerical_rank": 11},
                    {"constructor_parameter": "444/203", "stable_numerical_rank": 11},
                ],
            },
            "six_double_power_pairs": {
                "outside_box_proxy_below_190": 1818,
                "all_six_trace_thresholds": 17,
                "decontaminated_b200_at_least_17": 168,
                "h50000_multiobjective_leaders": 25,
                "maximum_stable_numerical_rank": 12,
            },
        },
        "leakage_control": {
            "design_primes_omitted_from_b200_and_exact_b2000": list(DESIGN_PRIMES),
            "b200_primes_used": list(tables),
            "b200_population_scored": len(candidates),
            "b200_retained": len(prefiltered),
            "exact_b2000_population_rescored": len(exactly_scored),
            "exact_point_population_retained": len(exact_survivors),
            "b200_tail_floor_preserved_for_points": B200_TAIL_FLOOR,
            "b200_tail_population": sum(
                candidate.residue_score_b200 >= B200_TAIL_FLOOR
                for candidate in prefiltered
            ),
        },
        "b200_ranked_tail": [
            {
                "candidate_id": candidate.identifier,
                "constructor_parameter": rational_to_string(candidate.parameter),
                "score": candidate.residue_score_b200,
                "good_primes_used": candidate.residue_good_primes,
                "bad_primes_skipped": candidate.residue_bad_primes,
                "radical_proxy": candidate.generated.radical_proxy,
                "stratum": candidate.generated.stratum,
            }
            for candidate in prefiltered
        ],
        "exact_b2000_ranked_population": [
            candidate_record(candidate) for candidate in exactly_scored
        ],
        "calibration_scores": candidate_record(calibration_exact),
        "point_stages": point_stages,
        "exact_rank_checkpoints": exact_checkpoints,
        "conductor_replays": conductor_replays,
        "summary": {
            "generated_proxy_survivors": len(candidates),
            "b200_survivors": len(prefiltered),
            "exact_b2000_point_survivors": len(exact_survivors),
            "deep_population": len(final_pools),
            "maximum_stable_numerical_rank": max(
                (
                    int(record["stable_numerical_rank"])
                    for record in final_ranks.values()
                    if record.get("status") == "completed"
                ),
                default=0,
            ),
            "stable_rank_at_least_18_candidate_ids": sorted(checkpoint_targets),
            "certified_target_hit": bool(certified_hits),
        },
        "bounds": {
            "b200_keep_count": args.b200_keep,
            "exact_point_keep_count": args.exact_point_keep,
            "exact_score_batch_size": args.exact_score_batch,
            "exact_score_timeout_seconds_per_batch": args.exact_score_timeout,
            "quartic_naive_height_stages": list(args.stage_heights),
            "point_stage_keep_counts": list(args.stage_keeps),
            "search_timeout_seconds_per_candidate": list(args.stage_timeouts),
            "point_worker_counts": list(args.stage_workers),
            "height_precisions": list(args.height_precisions),
            "height_timeout_seconds": args.height_timeout,
            "saturation_timeout_seconds": args.saturation_timeout,
            "conductor_timeout_seconds": args.conductor_timeout,
            "certificate_prime_bound": args.certificate_prime_bound,
            "pari_stack_bytes_per_process": args.stack_bytes,
            "process_policy": (
                "bounded thread pools contain synchronous subprocess.run "
                "calls; every pool joins all workers and children"
            ),
        },
        "interpretation": (
            "Beam and score stages define a bounded rare-event search, not a "
            "rank estimate. Radical proxies are selection features, not "
            "conductor bounds. Only stored finite reductions certify rank."
        ),
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "invocation": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "reproducing_command": REPRODUCING_COMMAND,
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True, default=str) + "\n")
    print(f"wrote {args.output}", flush=True)
    print(json.dumps(artifact["summary"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
