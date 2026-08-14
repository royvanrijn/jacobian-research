#!/usr/bin/env python3
"""Bounded CRT/parity neighborhood search around Nagao's rank-20 section-7 fiber.

The certified specialization ``T=5081/47`` is calibration only.  This script
learns the family invariants, its exceptional local trace symbols, and the
discriminant root balls responsible for the calibration's small radical.  It
then searches a declared union of CRT/Gauss, near-calibration, and smooth even
denominator strata.  Every selected parameter is positive, primitive, outside
the previous ``a<=10000,b<=100`` box, and unequal to the calibration.

All primes used to design the population are omitted from the B=200 and exact
PARI B=2000 scores.  Exact conductor and root-number replays precede bounded
point searches.  The twelve visible sections and nine additional generic
sections are included as predeclared seeds; bounded-search yield is therefore
decontaminated against all twenty-one.  Numerical height rank is triage only.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import gcd, log
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

from certify_nagao_rank20_t5081 import (
    CONSTRUCTION,
    EXPECTED_CONDUCTOR,
    EXPECTED_MINIMAL_DISCRIMINANT,
    PARAMETER_T as CALIBRATION_PARAMETER,
    ROOTS,
)
from crt_lattice import crt_pair, gauss_reduce
from ek_k3 import legendre_symbol, primes_up_to, rational_to_string
from multiple_root_lifting import (
    affine_variable_coefficients,
    fixed_divisor_valuation,
)
from nagao_1994 import (
    PRIMARY_SOURCE,
    primitive_quartic_coefficients,
    primitive_visible_points,
    quartic_point_to_short_jacobian,
    quartic_value,
    short_jacobian_coefficients,
)
from pari_bridge import minimal_curve_data, pari_version
from search_extra_points import parse_point_vector, run_gp, signless_quartic_points
from search_nagao_rank21_unbiased import finite_reduction_certificate
from triage_nagao_rank13_finalists import (
    height_matrix_replay,
    point_digest,
    point_on_short_curve,
    quartic_gp_polynomial,
    stable_height_rank,
)


Q = Fraction
TARGET_LOG_CONDUCTOR = Decimal("182.72")
OLD_BOX_A_MAX = 10_000
OLD_BOX_B_MAX = 100
TRACE_DESIGN_COUNT = 5
TRACE_BEAM_WIDTH = 4_000
ROOT_BEAM_WIDTH = 500
PROXY_TRIAL_BOUND = 2_000
PROXY_LIMIT = Decimal("190")
MAX_PROXY_SURVIVORS = 4_000
B200_KEEP = 400
CONDUCTOR_KEEP = 48
POINT_KEEP = 32
EXACT_SCORE_CUTOFF = 2_000
SAVING_PRIMES = (2, 3, 5, 7, 13, 17, 23)
EXPECTED_TRACE_DESIGN = (
    (53, 10, -12, 12),
    (109, 71, -20, 12),
    (151, 31, -23, 14),
    (163, 11, -24, 12),
    (197, 62, -27, 6),
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank20_t5081_neighborhood.py"
)


# Each tuple is (x slope, x intercept, ascending y coefficients).  These are
# exact primitive-quartic sections, not search discoveries.
COMPANION_SECTION_DATA = (
    (Q(7, 27), Q(6920, 27), (Q(84770), Q(-18554923, 243), Q(-29974, 243), Q(680, 243))),
    (Q(-7, 27), Q(6920, 27), (Q(-84770), Q(-18554923, 243), Q(29974, 243), Q(680, 243))),
    (Q(17, 27), Q(5462, 27), (Q(5138284), Q(-6202747, 243), Q(-23222, 243), Q(440, 243))),
    (Q(-17, 27), Q(5462, 27), (Q(-5138284), Q(-6202747, 243), Q(23222, 243), Q(440, 243))),
    (Q(43, 27), Q(-4015, 27), (Q(94091525), Q(-236806588, 243), Q(756284, 243), Q(-1120, 243))),
    (Q(-43, 27), Q(-4015, 27), (Q(-94091525), Q(-236806588, 243), Q(-756284, 243), Q(-1120, 243))),
)

# Three further exact generic sections have x=m*T^2+k.  Their specializations
# lie in the rank-12 generic span, but they remain mandatory contaminants for
# any bounded point-yield count.
QUADRATIC_COMPANION_SECTION_DATA = (
    (
        Q(56, 5373),
        Q(1389190, 5373),
        (Q(0), Q(684218797630, 9623043), Q(0), Q(171853351, 9623043), Q(0), Q(3136, 9623043)),
    ),
    (
        Q(-22, 5373),
        Q(1389190, 5373),
        (Q(0), Q(638541742570, 9623043), Q(0), Q(-24743777, 9623043), Q(0), Q(484, 9623043)),
    ),
    (
        Q(-34, 5373),
        Q(1389190, 5373),
        (Q(0), Q(-631221199130, 9623043), Q(0), Q(-2763929, 9623043), Q(0), Q(1156, 9623043)),
    ),
)

# Coordinates in the declared generic rank-12 basis
# (visible sections 0..10, then the +7/27 linear companion).  These exact
# relations were obtained by symbolic elliptic-curve group-law replay over
# Q(T); hence the quadratic sections do not increase the generic baseline.
QUADRATIC_COMPANION_RELATIONS = (
    (0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 1),
    (0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1),
    (-1, -1, -1, 0, 0, -1, -1, -1, 0, 0, 0, 1),
)


def polynomial_add(*values: Sequence[Fraction]) -> tuple[Fraction, ...]:
    length = max(len(value) for value in values)
    answer = tuple(
        sum(
            (value[index] if index < len(value) else Q(0))
            for value in values
        )
        for index in range(length)
    )
    while len(answer) > 1 and answer[-1] == 0:
        answer = answer[:-1]
    return answer


def polynomial_multiply(*values: Sequence[Fraction]) -> tuple[Fraction, ...]:
    answer = (Q(1),)
    for value in values:
        product = [Q(0)] * (len(answer) + len(value) - 1)
        for left_index, left in enumerate(answer):
            for right_index, right in enumerate(value):
                product[left_index + right_index] += left * right
        answer = tuple(product)
    return answer


def polynomial_scale(value: Sequence[Fraction], scale: Fraction) -> tuple[Fraction, ...]:
    return tuple(Q(scale) * coefficient for coefficient in value)


def polynomial_value(coefficients: Sequence[Fraction], value: Fraction) -> Fraction:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * Q(value) + coefficient
    return answer


def interpolate_polynomial(
    xs: Sequence[Fraction], ys: Sequence[Fraction]
) -> tuple[Fraction, ...]:
    """Return the exact power-basis interpolant through distinct points."""

    if len(xs) != len(ys) or not xs or len(set(xs)) != len(xs):
        raise ValueError("interpolation requires equally many distinct abscissas")
    answer = (Q(0),)
    for index, (x_value, y_value) in enumerate(zip(xs, ys)):
        basis = (Q(1),)
        denominator = Q(1)
        for other_index, other_x in enumerate(xs):
            if index == other_index:
                continue
            basis = polynomial_multiply(basis, (-Q(other_x), Q(1)))
            denominator *= Q(x_value) - Q(other_x)
        answer = polynomial_add(
            answer, polynomial_scale(basis, Q(y_value) / denominator)
        )
    return answer


def learn_invariant_polynomials() -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Interpolate and independently check the degree-8/12 primitive invariants."""

    i_xs = tuple(Q(value) for value in range(1, 10))
    j_xs = tuple(Q(value) for value in range(1, 14))
    invariant_i = interpolate_polynomial(
        i_xs,
        tuple(CONSTRUCTION.primitive_binary_invariants(value)[0] for value in i_xs),
    )
    invariant_j = interpolate_polynomial(
        j_xs,
        tuple(CONSTRUCTION.primitive_binary_invariants(value)[1] for value in j_xs),
    )
    if len(invariant_i) != 9 or len(invariant_j) != 13:
        raise AssertionError("the learned section-7 invariant degrees changed")
    if any(value.denominator != 1 for value in invariant_i + invariant_j):
        raise AssertionError("the learned section-7 invariants are not integral")
    for parameter in (Q(14), Q(17, 3), CALIBRATION_PARAMETER):
        expected = CONSTRUCTION.primitive_binary_invariants(parameter)
        observed = (
            polynomial_value(invariant_i, parameter),
            polynomial_value(invariant_j, parameter),
        )
        if observed != expected:
            raise AssertionError("section-7 invariant interpolation failed")
    return (
        tuple(int(value) for value in invariant_i),
        tuple(int(value) for value in invariant_j),
    )


INVARIANT_I, INVARIANT_J = learn_invariant_polynomials()
DISCRIMINANT_POLYNOMIAL = tuple(
    int(value) for value in CONSTRUCTION.primitive_discriminant_polynomial
)
if len(DISCRIMINANT_POLYNOMIAL) != 21:
    raise AssertionError("the section-7 discriminant degree changed")


def homogeneous_value_mod(
    coefficients: Sequence[int], numerator: int, denominator: int, prime: int
) -> int:
    degree = len(coefficients) - 1
    return sum(
        coefficient
        * pow(numerator, power, prime)
        * pow(denominator, degree - power, prime)
        for power, coefficient in enumerate(coefficients)
    ) % prime


def homogenized_discriminant(parameter: Fraction) -> int:
    parameter = abs(Q(parameter))
    numerator, denominator = parameter.numerator, parameter.denominator
    degree = len(DISCRIMINANT_POLYNOMIAL) - 1
    value = sum(
        coefficient * numerator**power * denominator ** (degree - power)
        for power, coefficient in enumerate(DISCRIMINANT_POLYNOMIAL)
    )
    if value == 0:
        raise ValueError("singular section-7 specialization")
    return value


def integer_valuation(value: int, prime: int) -> int:
    value = abs(value)
    if value == 0:
        raise ValueError("zero has no finite valuation")
    answer = 0
    while value % prime == 0:
        value //= prime
        answer += 1
    return answer


def conductor_radical_proxy(
    parameter: Fraction, *, trial_prime_bound: int = PROXY_TRIAL_BOUND
) -> dict[str, Any]:
    discriminant = abs(homogenized_discriminant(parameter))
    remaining = discriminant
    savings = 0.0
    valuations = []
    for prime in primes_up_to(trial_prime_bound):
        valuation = 0
        while remaining % prime == 0:
            remaining //= prime
            valuation += 1
        if valuation:
            valuations.append([prime, valuation])
            savings += (valuation - 1) * log(prime)
    raw_log = log(discriminant)
    return {
        "trial_prime_bound": trial_prime_bound,
        "raw_log_absolute_homogenized_discriminant": raw_log,
        "known_repeated_prime_log_savings": savings,
        "log_radical_upper_proxy": raw_log - savings,
        "small_prime_valuations": valuations,
        "unfactored_cofactor_decimal_digits": len(str(remaining)),
    }


@dataclass(frozen=True)
class ResidueSymbol:
    prime: int
    residue: int
    ellap: int | None
    good_reduction: bool
    contribution: float


def projective_index(numerator: int, denominator: int, prime: int) -> int:
    if denominator % prime == 0:
        if numerator % prime == 0:
            raise ValueError("a nonprimitive pair has no projective reduction")
        return prime
    return numerator * pow(denominator, -1, prime) % prime


def residue_table(prime: int) -> tuple[ResidueSymbol, ...]:
    if prime < 5 or prime not in primes_up_to(prime):
        raise ValueError("residue tables require a prime at least five")
    symbols = []
    for index in range(prime + 1):
        numerator, denominator = (index, 1) if index < prime else (1, 0)
        coefficient_a = -27 * homogeneous_value_mod(
            INVARIANT_I, numerator, denominator, prime
        ) % prime
        coefficient_b = -27 * homogeneous_value_mod(
            INVARIANT_J, numerator, denominator, prime
        ) % prime
        if (4 * coefficient_a**3 + 27 * coefficient_b**2) % prime == 0:
            symbols.append(ResidueSymbol(prime, index, None, False, 0.0))
            continue
        trace = -sum(
            legendre_symbol(
                (x_value**3 + coefficient_a * x_value + coefficient_b) % prime,
                prime,
            )
            for x_value in range(prime)
        )
        contribution = (2 - trace) / (prime + 1 - trace) * log(prime)
        symbols.append(ResidueSymbol(prime, index, trace, True, contribution))
    return tuple(symbols)


def build_residue_tables(cutoff: int) -> dict[int, tuple[ResidueSymbol, ...]]:
    return {
        prime: residue_table(prime)
        for prime in primes_up_to(cutoff)
        if prime >= 5
    }


def learn_local_trace_fingerprint(
    tables: dict[int, tuple[ResidueSymbol, ...]], *, count: int = TRACE_DESIGN_COUNT
) -> tuple[tuple[int, int, int, int], ...]:
    """Select rare negative-trace symbols by information-weighted contribution."""

    ranked = []
    for prime, table in tables.items():
        if prime < 11:
            continue
        index = projective_index(
            CALIBRATION_PARAMETER.numerator,
            CALIBRATION_PARAMETER.denominator,
            prime,
        )
        symbol = table[index]
        if not symbol.good_reduction or symbol.ellap is None or symbol.ellap >= 0:
            continue
        union_size = sum(
            candidate.residue < prime
            and candidate.good_reduction
            and candidate.ellap is not None
            and candidate.ellap <= symbol.ellap
            for candidate in table
        )
        if union_size == 0:
            continue
        information_weight = symbol.contribution * log(prime / union_size)
        ranked.append(
            (
                -information_weight,
                prime,
                index,
                int(symbol.ellap),
                union_size,
            )
        )
    selected = sorted(ranked)[:count]
    answer = tuple(
        sorted((prime, index, trace, union) for _, prime, index, trace, union in selected)
    )
    if count == TRACE_DESIGN_COUNT and answer != EXPECTED_TRACE_DESIGN:
        raise AssertionError("the automatically learned trace fingerprint changed")
    return answer


def residue_score(
    parameter: Fraction, tables: dict[int, tuple[ResidueSymbol, ...]]
) -> tuple[float, int, int]:
    score = 0.0
    good = 0
    bad = 0
    for prime, table in tables.items():
        symbol = table[
            projective_index(parameter.numerator, parameter.denominator, prime)
        ]
        if symbol.good_reduction:
            score += symbol.contribution
            good += 1
        else:
            bad += 1
    return score, good, bad


@dataclass(frozen=True)
class RootBall:
    prime: int
    residue: int
    modulus: int
    forced_valuation: int


def root_ball_union(prime: int, target_valuation: int) -> tuple[RootBall, ...]:
    """Return disjoint affine balls on which v_p(D) is at least the target."""

    if prime < 2 or target_valuation < 1:
        raise ValueError("invalid root-ball request")
    active = [
        (residue, prime)
        for residue in range(prime)
        if polynomial_value(DISCRIMINANT_POLYNOMIAL, Q(residue)).numerator % prime == 0
    ]
    completed: list[RootBall] = []
    while active:
        residue, modulus = active.pop()
        forced = fixed_divisor_valuation(
            affine_variable_coefficients(
                DISCRIMINANT_POLYNOMIAL, residue, modulus
            ),
            prime,
        )
        if forced >= target_valuation:
            completed.append(RootBall(prime, residue, modulus, forced))
            continue
        if modulus > prime ** (target_valuation + 1):
            raise AssertionError("a discriminant root ball failed to stabilize")
        lifted_modulus = modulus * prime
        for lift in range(prime):
            lifted = residue + lift * modulus
            if (
                polynomial_value(DISCRIMINANT_POLYNOMIAL, Q(lifted)).numerator
                % lifted_modulus
                == 0
            ):
                active.append((lifted, lifted_modulus))
    return tuple(
        sorted(completed, key=lambda ball: (ball.modulus, ball.residue))
    )


def learn_discriminant_root_balls() -> dict[int, tuple[RootBall, ...]]:
    lead_discriminant = homogenized_discriminant(CALIBRATION_PARAMETER)
    targets = {
        prime: integer_valuation(lead_discriminant, prime)
        for prime in SAVING_PRIMES
    }
    expected_targets = {2: 6, 3: 14, 5: 7, 7: 6, 13: 2, 17: 2, 23: 2}
    if targets != expected_targets:
        raise AssertionError("the calibration discriminant valuations changed")
    answer = {
        prime: root_ball_union(prime, target)
        for prime, target in targets.items()
    }
    for prime, balls in answer.items():
        if not any(
            (CALIBRATION_PARAMETER.numerator - ball.residue * CALIBRATION_PARAMETER.denominator)
            % ball.modulus
            == 0
            for ball in balls
        ):
            raise AssertionError(f"the calibration left its p={prime} root-ball union")
    return answer


@dataclass(frozen=True)
class TraceChoice:
    prime: int
    residue: int
    ellap: int
    contribution: float


@dataclass(frozen=True)
class BeamState:
    residue: int
    modulus: int
    trace_choices: tuple[TraceChoice, ...]
    root_balls: tuple[RootBall, ...]
    representative: Fraction
    basis: tuple[tuple[int, int], tuple[int, int]]
    objective: float


def gauss_shell(
    residue: int,
    modulus: int,
    *,
    radius: int = 4,
    limit: int = 12,
) -> tuple[tuple[Fraction, tuple[int, int], tuple[tuple[int, int], tuple[int, int]]], ...]:
    if modulus < 1 or not 0 <= residue < modulus:
        raise ValueError("invalid CRT class")
    basis = gauss_reduce((modulus, 0), (residue, 1))
    found: dict[Fraction, tuple[int, int]] = {}
    for left in range(-radius, radius + 1):
        for right in range(-radius, radius + 1):
            if left == 0 and right == 0:
                continue
            numerator = left * basis[0][0] + right * basis[1][0]
            denominator = left * basis[0][1] + right * basis[1][1]
            if denominator == 0:
                continue
            parameter = abs(Q(numerator, denominator))
            if parameter == 0 or gcd(parameter.denominator, modulus) != 1:
                continue
            if (parameter.numerator - residue * parameter.denominator) % modulus and (
                parameter.numerator + residue * parameter.denominator
            ) % modulus:
                continue
            pair = (left, right)
            if parameter not in found or pair < found[parameter]:
                found[parameter] = pair
    ordered = sorted(
        found.items(),
        key=lambda item: (
            max(item[0].numerator, item[0].denominator),
            item[0].numerator,
            item[0].denominator,
            item[1],
        ),
    )
    return tuple((parameter, pair, basis) for parameter, pair in ordered[:limit])


def make_state(
    residue: int,
    modulus: int,
    traces: tuple[TraceChoice, ...],
    balls: tuple[RootBall, ...],
) -> BeamState:
    shell = gauss_shell(residue, modulus, radius=2, limit=1)
    if not shell:
        shell = gauss_shell(residue, modulus, radius=6, limit=1)
    if shell:
        representative, _, basis = shell[0]
    else:
        # The original lattice vector is always a valid affine unit class.
        # It is intentionally expensive in the height objective, but keeps a
        # valid state representable when a very short reduced-basis shell has
        # only non-unit denominators.
        representative = Q(residue, 1)
        basis = gauss_reduce((modulus, 0), (residue, 1))
    height = max(representative.numerator, representative.denominator)
    trace_reward = sum(choice.contribution for choice in traces)
    savings = sum((ball.forced_valuation - 1) * log(ball.prime) for ball in balls)
    return BeamState(
        residue,
        modulus,
        traces,
        balls,
        representative,
        basis,
        20 * log(height) - 2 * trace_reward - savings,
    )


def retain_beam(states: Iterable[BeamState], width: int) -> tuple[BeamState, ...]:
    retained: dict[tuple[int, int], BeamState] = {}
    for state in states:
        key = state.residue, state.modulus
        previous = retained.get(key)
        if previous is None or (
            state.objective,
            max(state.representative.numerator, state.representative.denominator),
            state.residue,
        ) < (
            previous.objective,
            max(previous.representative.numerator, previous.representative.denominator),
            previous.residue,
        ):
            retained[key] = state
    return tuple(
        sorted(
            retained.values(),
            key=lambda state: (
                state.objective,
                max(state.representative.numerator, state.representative.denominator),
                state.residue,
            ),
        )[:width]
    )


def trace_choice_union(
    prime: int, threshold: int, tables: dict[int, tuple[ResidueSymbol, ...]]
) -> tuple[TraceChoice, ...]:
    # Infinity is excluded here so every selected symbol has the affine
    # congruence a-rb=0 and can enter the two-dimensional CRT lattice.
    return tuple(
        TraceChoice(prime, symbol.residue, int(symbol.ellap), symbol.contribution)
        for symbol in tables[prime]
        if symbol.residue < prime
        and symbol.good_reduction
        and symbol.ellap is not None
        and symbol.ellap <= threshold
    )


def build_trace_beams(
    fingerprint: Sequence[tuple[int, int, int, int]],
    tables: dict[int, tuple[ResidueSymbol, ...]],
    *,
    width: int,
) -> tuple[tuple[BeamState, ...], tuple[BeamState, ...], tuple[dict[str, Any], ...]]:
    states = (make_state(0, 1, (), ()),)
    prefix = states
    audit = []
    for stage, (prime, _, threshold, _) in enumerate(fingerprint, start=1):
        choices = trace_choice_union(prime, threshold, tables)
        expanded = []
        for state in states:
            for choice in choices:
                residue, modulus = crt_pair(
                    state.residue, state.modulus, choice.residue, prime
                )
                expanded.append(
                    make_state(
                        residue,
                        modulus,
                        state.trace_choices + (choice,),
                        state.root_balls,
                    )
                )
        states = retain_beam(expanded, width)
        if stage == 3:
            prefix = states
        audit.append(
            {
                "prime": prime,
                "threshold": threshold,
                "affine_symbol_union_size": len(choices),
                "expanded": len(expanded),
                "retained": len(states),
            }
        )
    return states, prefix, tuple(audit)


def extend_root_beam(
    states: Sequence[BeamState], balls: Sequence[RootBall], *, width: int
) -> tuple[BeamState, ...]:
    expanded = []
    for state in states:
        for ball in balls:
            if gcd(state.modulus, ball.modulus) != 1:
                raise AssertionError("a local design prime was repeated")
            residue, modulus = crt_pair(
                state.residue, state.modulus, ball.residue, ball.modulus
            )
            expanded.append(
                make_state(
                    residue,
                    modulus,
                    state.trace_choices,
                    state.root_balls + (ball,),
                )
            )
    return retain_beam(expanded, width)


def build_beam_strata(
    full_trace: Sequence[BeamState],
    prefix_trace: Sequence[BeamState],
    root_balls: dict[int, tuple[RootBall, ...]],
    *,
    root_width: int,
) -> dict[str, tuple[BeamState, ...]]:
    strata: dict[str, tuple[BeamState, ...]] = {
        "trace-full": tuple(full_trace),
        "trace-prefix-3": tuple(prefix_trace),
    }
    for prime in SAVING_PRIMES:
        strata[f"trace-full-power-{prime}"] = extend_root_beam(
            full_trace, root_balls[prime], width=root_width
        )
    for label, primes in (
        ("power-3-5", (3, 5)),
        ("power-3-7", (3, 7)),
        ("power-5-7", (5, 7)),
        ("power-13-17-23", (13, 17, 23)),
        ("power-2-3-5-7", (2, 3, 5, 7)),
    ):
        states = tuple(prefix_trace)
        for prime in primes:
            states = extend_root_beam(states, root_balls[prime], width=root_width)
        strata[f"trace-prefix-{label}"] = states
    return strata


@dataclass(frozen=True)
class GeneratedCandidate:
    parameter: Fraction
    stratum: str
    conditions: tuple[str, ...]
    proxy: dict[str, Any]

    @property
    def identifier(self) -> str:
        return f"section7-{self.parameter.numerator}-{self.parameter.denominator}"

    @property
    def height(self) -> int:
        return max(self.parameter.numerator, self.parameter.denominator)


def outside_prior_box(parameter: Fraction) -> bool:
    parameter = abs(Q(parameter))
    return not (
        parameter.numerator <= OLD_BOX_A_MAX
        and parameter.denominator <= OLD_BOX_B_MAX
    )


def smooth_even_denominators() -> tuple[int, ...]:
    values = set()
    for power_two in range(4, 11):
        for power_three in range(0, 7):
            base = 2**power_two * 3**power_three
            for extra in (1, 5, 7):
                denominator = base * extra
                if 101 <= denominator <= 2_000:
                    values.add(denominator)
    return tuple(sorted(values))


def nearest_integer(value: Fraction) -> int:
    quotient, remainder = divmod(Q(value).numerator, Q(value).denominator)
    return quotient + (2 * remainder >= Q(value).denominator)


def raw_parameter_strata(
    beam_strata: dict[str, tuple[BeamState, ...]],
) -> tuple[tuple[Fraction, str, tuple[str, ...]], ...]:
    raw: list[tuple[Fraction, str, tuple[str, ...]]] = []
    for stratum, states in beam_strata.items():
        for state in states:
            shell = gauss_shell(state.residue, state.modulus, radius=4, limit=12)
            indices = (0, 2, 5) if "power" not in stratum else (0, 3)
            conditions = tuple(
                [
                    f"a_{choice.prime}<={choice.ellap}@{choice.residue}"
                    for choice in state.trace_choices
                ]
                + [
                    f"v{ball.prime}>={ball.forced_valuation}@{ball.residue}mod{ball.modulus}"
                    for ball in state.root_balls
                ]
            )
            for index in indices:
                if index < len(shell):
                    raw.append((shell[index][0], stratum, conditions))

    # A rational neighborhood outside the old denominator box.  It is not
    # conditioned on the learned trace primes and thus provides an independent
    # check against overfitting the CRT design.
    for denominator in range(101, 321):
        center = nearest_integer(CALIBRATION_PARAMETER * denominator)
        for offset in range(-24, 25):
            raw.append((Q(center + offset, denominator), "near-lead-shell", ()))

    # Explicitly retain odd-numerator/even-smooth-denominator fibers.  The
    # offsets are deliberately wider than the generic near-lead shell.
    for denominator in smooth_even_denominators():
        center = nearest_integer(CALIBRATION_PARAMETER * denominator)
        for offset in range(-64, 65, 4):
            numerator = center + offset
            if numerator % 2 == 0:
                numerator += 1
            raw.append(
                (
                    Q(numerator, denominator),
                    "smooth-even-denominator",
                    ("odd numerator; denominator 2^u*3^v times optional 5 or 7",),
                )
            )
    return tuple(raw)


def generate_candidates(
    beam_strata: dict[str, tuple[BeamState, ...]],
    *,
    proxy_limit: Decimal,
    max_survivors: int,
) -> tuple[tuple[GeneratedCandidate, ...], dict[str, Any]]:
    raw = raw_parameter_strata(beam_strata)
    retained: dict[Fraction, GeneratedCandidate] = {}
    excluded_box = 0
    excluded_lead = 0
    singular = 0
    proxy_rejected = 0
    for parameter, stratum, conditions in raw:
        parameter = abs(Q(parameter))
        if parameter == 0:
            continue
        if not outside_prior_box(parameter):
            excluded_box += 1
            continue
        if parameter == CALIBRATION_PARAMETER:
            excluded_lead += 1
            continue
        try:
            proxy = conductor_radical_proxy(parameter)
        except ValueError:
            singular += 1
            continue
        if Decimal(str(proxy["log_radical_upper_proxy"])) >= proxy_limit:
            proxy_rejected += 1
            continue
        candidate = GeneratedCandidate(parameter, stratum, conditions, proxy)
        previous = retained.get(parameter)
        if previous is None or (
            proxy["log_radical_upper_proxy"], stratum, conditions
        ) < (
            previous.proxy["log_radical_upper_proxy"],
            previous.stratum,
            previous.conditions,
        ):
            retained[parameter] = candidate

    ordered_all = sorted(
        retained.values(),
        key=lambda candidate: (
            candidate.proxy["log_radical_upper_proxy"],
            candidate.height,
            candidate.identifier,
        ),
    )
    # Preserve a small exact quota from every construction stratum before the
    # global proxy fill, so a single dense neighborhood cannot erase CRT arms.
    selected: dict[Fraction, GeneratedCandidate] = {}
    strata = sorted({candidate.stratum for candidate in ordered_all})
    quota = max(1, min(40, max_survivors // max(1, len(strata))))
    for stratum in strata:
        for candidate in (item for item in ordered_all if item.stratum == stratum):
            selected[candidate.parameter] = candidate
            if sum(item.stratum == stratum for item in selected.values()) >= quota:
                break
    for candidate in ordered_all:
        if len(selected) >= max_survivors:
            break
        selected.setdefault(candidate.parameter, candidate)
    answer = tuple(
        sorted(
            selected.values(),
            key=lambda candidate: (
                candidate.proxy["log_radical_upper_proxy"],
                candidate.height,
                candidate.identifier,
            ),
        )
    )
    if any(candidate.parameter == CALIBRATION_PARAMETER for candidate in answer):
        raise AssertionError("the calibration leaked into the generated population")
    if any(not outside_prior_box(candidate.parameter) for candidate in answer):
        raise AssertionError("an old-box parameter leaked into the population")
    digest = hashlib.sha256()
    for candidate in answer:
        digest.update(
            f"{candidate.parameter}|{candidate.stratum}|{candidate.proxy['log_radical_upper_proxy']!r}\n".encode()
        )
    return answer, {
        "raw_parameter_records": len(raw),
        "excluded_old_positive_box": excluded_box,
        "excluded_calibration": excluded_lead,
        "singular": singular,
        "proxy_rejected": proxy_rejected,
        "deduplicated_below_proxy_before_cap": len(ordered_all),
        "retained_after_stratified_cap": len(answer),
        "stratum_counts": {
            stratum: sum(candidate.stratum == stratum for candidate in answer)
            for stratum in strata
        },
        "survivor_stream_sha256": digest.hexdigest(),
    }


@dataclass(frozen=True)
class PrefilterCandidate:
    generated: GeneratedCandidate
    score_b200: float
    good_primes: int
    bad_primes: int

    @property
    def parameter(self) -> Fraction:
        return self.generated.parameter

    @property
    def identifier(self) -> str:
        return self.generated.identifier


@dataclass(frozen=True)
class ExactCandidate:
    prefilter: PrefilterCandidate
    score_b2000: str
    good_primes: int
    bad_primes: int
    omitted_primes: int
    last_prime: int

    @property
    def parameter(self) -> Fraction:
        return self.prefilter.parameter

    @property
    def identifier(self) -> str:
        return self.prefilter.identifier


def prefilter_candidates(
    candidates: Sequence[GeneratedCandidate],
    *,
    tables: dict[int, tuple[ResidueSymbol, ...]],
    keep: int,
) -> tuple[PrefilterCandidate, ...]:
    scored = tuple(
        PrefilterCandidate(candidate, *residue_score(candidate.parameter, tables))
        for candidate in candidates
    )
    return tuple(
        sorted(
            scored,
            key=lambda candidate: (
                -candidate.score_b200,
                candidate.generated.proxy["log_radical_upper_proxy"],
                candidate.generated.height,
                candidate.identifier,
            ),
        )[: min(keep, len(scored))]
    )


def gp_rational(value: Fraction) -> str:
    return f"({value.numerator}/{value.denominator})"


def exact_decontaminated_scores(
    candidates: Sequence[PrefilterCandidate],
    *,
    omitted_primes: Sequence[int],
    cutoff: int,
    batch_size: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[ExactCandidate, ...]:
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    if batch_size < 1 or timeout <= 0:
        raise ValueError("invalid exact-score bounds")
    omitted = tuple(sorted(set(omitted_primes)))
    omitted_test = "||".join(f"p=={prime}" for prime in omitted if prime >= 5)
    if not omitted_test:
        omitted_test = "0"
    last_prime = primes_up_to(cutoff)[-1]
    records: dict[str, dict[str, Any]] = {}
    for start in range(0, len(candidates), batch_size):
        batch = candidates[start : start + batch_size]
        commands = ["default(realprecision,80);"]
        for index, candidate in enumerate(batch):
            coefficients = short_jacobian_coefficients(
                CONSTRUCTION, candidate.parameter
            )
            vector = ",".join(gp_rational(value) for value in coefficients)
            commands.extend(
                (
                    f"E=ellminimalmodel(ellinit([{vector}]));",
                    "S=0;USED=0;BAD=0;OMITTED=0;",
                    (
                        f"forprime(p=5,{cutoff},if({omitted_test},OMITTED++,"
                        "if(valuation(E.disc,p)>0,BAD++,"
                        "A=ellap(E,p);S+=(2-A)/(p+1-A)*log(p);USED++)));"
                    ),
                    f'print("ROW|{index}|",S,"|",USED,"|",BAD,"|",OMITTED);',
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
            raise RuntimeError(f"PARI exact score failed: {result.stderr.strip()}")
        observed = 0
        for line in result.stdout.splitlines():
            if not line.startswith("ROW|"):
                continue
            _, index_text, score, used, bad, omitted_count = line.split("|")
            candidate = batch[int(index_text)]
            records[candidate.identifier] = {
                "score": score,
                "used": int(used),
                "bad": int(bad),
                "omitted": int(omitted_count),
            }
            observed += 1
        if observed != len(batch):
            raise RuntimeError("PARI omitted exact-score output")
    answer = tuple(
        ExactCandidate(
            candidate,
            records[candidate.identifier]["score"],
            records[candidate.identifier]["used"],
            records[candidate.identifier]["bad"],
            records[candidate.identifier]["omitted"],
            last_prime,
        )
        for candidate in candidates
    )
    return tuple(
        sorted(
            answer,
            key=lambda candidate: (
                -Decimal(candidate.score_b2000),
                candidate.prefilter.generated.proxy["log_radical_upper_proxy"],
                candidate.prefilter.generated.height,
                candidate.identifier,
            ),
        )
    )


def select_conductor_population(
    candidates: Sequence[ExactCandidate], *, keep: int
) -> tuple[ExactCandidate, ...]:
    selected: dict[str, ExactCandidate] = {}
    for candidate in candidates[: min(24, len(candidates))]:
        selected[candidate.identifier] = candidate
    by_proxy = sorted(
        candidates,
        key=lambda candidate: (
            candidate.prefilter.generated.proxy["log_radical_upper_proxy"],
            -Decimal(candidate.score_b2000),
            candidate.identifier,
        ),
    )
    for candidate in by_proxy[: min(16, len(by_proxy))]:
        selected[candidate.identifier] = candidate
    for stratum in sorted(
        {candidate.prefilter.generated.stratum for candidate in candidates}
    ):
        leader = next(
            candidate
            for candidate in candidates
            if candidate.prefilter.generated.stratum == stratum
        )
        selected[leader.identifier] = leader
    for candidate in candidates:
        if len(selected) >= keep:
            break
        selected.setdefault(candidate.identifier, candidate)
    return tuple(
        sorted(
            selected.values(),
            key=lambda candidate: (
                -Decimal(candidate.score_b2000), candidate.identifier
            ),
        )[:keep]
    )


@dataclass(frozen=True)
class ConductorReplay:
    candidate: ExactCandidate
    status: str
    data: dict[str, Any]
    error: str | None = None


def replay_one_conductor(
    candidate: ExactCandidate,
    *,
    timeout: float,
    stack_bytes: int,
) -> ConductorReplay:
    try:
        data = minimal_curve_data(
            short_jacobian_coefficients(CONSTRUCTION, candidate.parameter),
            timeout=timeout,
            local_primes=SAVING_PRIMES,
            stack_bytes=stack_bytes,
        )
        data["below_strict_log_conductor_target"] = (
            Decimal(data["log_conductor"]) < TARGET_LOG_CONDUCTOR
        )
        return ConductorReplay(candidate, "completed", data)
    except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
        return ConductorReplay(
            candidate,
            "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
            {},
            str(error)[:500],
        )


def parallel_conductor_replay(
    candidates: Sequence[ExactCandidate],
    *,
    timeout: float,
    stack_bytes: int,
    workers: int,
) -> tuple[ConductorReplay, ...]:
    by_id: dict[str, ConductorReplay] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                replay_one_conductor,
                candidate,
                timeout=timeout,
                stack_bytes=stack_bytes,
            ): candidate.identifier
            for candidate in candidates
        }
        for future in as_completed(futures):
            by_id[futures[future]] = future.result()
    return tuple(by_id[candidate.identifier] for candidate in candidates)


def select_point_population(
    replays: Sequence[ConductorReplay], *, keep: int
) -> tuple[ExactCandidate, ...]:
    selected: dict[str, ExactCandidate] = {}
    completed = [replay for replay in replays if replay.status == "completed"]
    root_minus_feasible = sorted(
        (
            replay
            for replay in completed
            if replay.data.get("below_strict_log_conductor_target") is True
            and int(replay.data["root_number"]) == -1
        ),
        key=lambda replay: (
            -Decimal(replay.candidate.score_b2000), replay.candidate.identifier
        ),
    )
    for replay in root_minus_feasible:
        selected[replay.candidate.identifier] = replay.candidate
        if len(selected) >= keep:
            break
    # Exact conductor-feasible score leaders of either parity remain eligible.
    feasible = sorted(
        (
            replay
            for replay in completed
            if replay.data.get("below_strict_log_conductor_target") is True
        ),
        key=lambda replay: (
            -Decimal(replay.candidate.score_b2000), replay.candidate.identifier
        ),
    )
    for replay in feasible[:12]:
        selected.setdefault(replay.candidate.identifier, replay.candidate)
    # Preserve four pure score leaders regardless of conductor/parity outcome.
    for replay in replays[:4]:
        selected.setdefault(replay.candidate.identifier, replay.candidate)
    return tuple(
        sorted(
            selected.values(),
            key=lambda candidate: (
                candidate.identifier not in {
                    replay.candidate.identifier for replay in root_minus_feasible
                },
                -Decimal(candidate.score_b2000),
                candidate.identifier,
            ),
        )[:keep]
    )


def companion_quartic_points(parameter: Fraction) -> tuple[tuple[Fraction, Fraction], ...]:
    parameter = Q(parameter)
    quartic = primitive_quartic_coefficients(CONSTRUCTION, parameter)
    points = []
    for slope, intercept, y_coefficients in COMPANION_SECTION_DATA:
        x_value = slope * parameter + intercept
        y_value = polynomial_value(y_coefficients, parameter)
        if y_value**2 != quartic_value(quartic, x_value):
            raise AssertionError("a generic companion section missed the quartic")
        points.append((x_value, y_value))
    return tuple(points)


def quadratic_companion_quartic_points(
    parameter: Fraction,
) -> tuple[tuple[Fraction, Fraction], ...]:
    parameter = Q(parameter)
    quartic = primitive_quartic_coefficients(CONSTRUCTION, parameter)
    points = []
    for quadratic_slope, intercept, y_coefficients in QUADRATIC_COMPANION_SECTION_DATA:
        x_value = quadratic_slope * parameter**2 + intercept
        y_value = polynomial_value(y_coefficients, parameter)
        if y_value**2 != quartic_value(quartic, x_value):
            raise AssertionError("a quadratic generic companion missed the quartic")
        points.append((x_value, y_value))
    return tuple(points)


def exact_predeclared_seeds(
    parameter: Fraction,
) -> tuple[
    tuple[tuple[Fraction, Fraction], ...],
    tuple[tuple[Fraction, Fraction], ...],
    tuple[Fraction, ...],
]:
    visible = primitive_visible_points(CONSTRUCTION, parameter)
    companions = companion_quartic_points(parameter) + quadratic_companion_quartic_points(parameter)
    quartic_points = visible + companions
    if len(visible) != 12 or len(companions) != 9:
        raise AssertionError("the section-7 seed counts changed")
    coefficients = short_jacobian_coefficients(CONSTRUCTION, parameter)
    images = tuple(
        quartic_point_to_short_jacobian(CONSTRUCTION, parameter, point)
        for point in quartic_points
    )
    if any(not point_on_short_curve(coefficients, point) for point in images):
        raise AssertionError("a predeclared section missed the Jacobian")
    # Special fibers may collide.  Retain one exact sign representative per x.
    by_x: dict[Fraction, tuple[Fraction, Fraction]] = {}
    quartic_by_x: dict[Fraction, tuple[Fraction, Fraction]] = {}
    for quartic_point, image in zip(quartic_points, images):
        quartic_by_x.setdefault(quartic_point[0], quartic_point)
        by_x.setdefault(image[0], image)
    return tuple(quartic_by_x.values()), tuple(by_x.values()), coefficients


def bounded_quartic_points(
    parameter: Fraction,
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> tuple[tuple[tuple[Fraction, Fraction], ...], float, int]:
    quartic = primitive_quartic_coefficients(CONSTRUCTION, parameter)
    program = "\n".join(
        (
            f"Q={quartic_gp_polynomial(quartic)};",
            "gettime();",
            f"R=hyperellratpoints(Q,{height_bound});",
            'print("PARI_MILLISECONDS ",gettime());',
            'print("POINTS ",R);',
            "quit",
        )
    ) + "\n"
    output, wall_seconds = run_gp(
        program, timeout=timeout, stack_bytes=stack_bytes
    )
    match = re.search(r"PARI_MILLISECONDS (\d+)", output)
    if match is None or "POINTS " not in output:
        raise AssertionError("PARI omitted bounded point output")
    return (
        parse_point_vector(output.split("POINTS ", 1)[1]),
        wall_seconds,
        int(match.group(1)),
    )


def map_new_points(
    parameter: Fraction,
    raw_points: Sequence[tuple[Fraction, Fraction]],
    seed_quartic: Sequence[tuple[Fraction, Fraction]],
    seed_images: Sequence[tuple[Fraction, Fraction]],
) -> tuple[tuple[tuple[Fraction, Fraction], ...], int, int]:
    coefficients = short_jacobian_coefficients(CONSTRUCTION, parameter)
    seed_quartic_x = {point[0] for point in seed_quartic}
    seen_image_x = {point[0] for point in seed_images}
    new_images = []
    predeclared_returned = 0
    zero_ordinates = 0
    for quartic_point in signless_quartic_points(tuple(raw_points)):
        if quartic_point[1] == 0:
            zero_ordinates += 1
            continue
        if quartic_point[0] in seed_quartic_x:
            predeclared_returned += 1
        image = quartic_point_to_short_jacobian(
            CONSTRUCTION, parameter, quartic_point
        )
        if not point_on_short_curve(coefficients, image):
            raise AssertionError("a searched quartic point missed the Jacobian")
        if image[0] in seen_image_x:
            continue
        seen_image_x.add(image[0])
        new_images.append(image)
    return tuple(new_images), predeclared_returned, zero_ordinates


@dataclass(frozen=True)
class PointPool:
    candidate: ExactCandidate
    height_bound: int
    status: str
    signed_points: int
    signless_points: int
    predeclared_abscissas_returned: int
    new_images: tuple[tuple[Fraction, Fraction], ...]
    seed_images: tuple[tuple[Fraction, Fraction], ...]
    coefficients: tuple[Fraction, ...]
    wall_seconds: float
    pari_milliseconds: int
    error: str | None = None

    @property
    def pool(self) -> tuple[tuple[Fraction, Fraction], ...]:
        return self.seed_images + self.new_images


def search_one_pool(
    candidate: ExactCandidate,
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
) -> PointPool:
    started = time.monotonic()
    try:
        seed_quartic, seed_images, coefficients = exact_predeclared_seeds(
            candidate.parameter
        )
        raw, wall, milliseconds = bounded_quartic_points(
            candidate.parameter,
            height_bound=height_bound,
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        new_images, returned, _ = map_new_points(
            candidate.parameter, raw, seed_quartic, seed_images
        )
        return PointPool(
            candidate,
            height_bound,
            "completed",
            len(raw),
            len({point[0] for point in raw}),
            returned,
            new_images,
            seed_images,
            coefficients,
            wall,
            milliseconds,
        )
    except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
        return PointPool(
            candidate,
            height_bound,
            "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
            0,
            0,
            0,
            (),
            (),
            (),
            time.monotonic() - started,
            0,
            str(error)[:500],
        )


def parallel_point_search(
    candidates: Sequence[ExactCandidate],
    *,
    height_bound: int,
    timeout: float,
    stack_bytes: int,
    workers: int,
) -> tuple[PointPool, ...]:
    by_id: dict[str, PointPool] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                search_one_pool,
                candidate,
                height_bound=height_bound,
                timeout=timeout,
                stack_bytes=stack_bytes,
            ): candidate.identifier
            for candidate in candidates
        }
        for future in as_completed(futures):
            by_id[futures[future]] = future.result()
    return tuple(by_id[candidate.identifier] for candidate in candidates)


def rank_one_pool(
    pool: PointPool,
    *,
    precisions: tuple[int, ...],
    timeout: float,
    stack_bytes: int,
) -> dict[str, Any]:
    if pool.status != "completed":
        return {"status": "not-run", "reason": "point search did not complete"}
    try:
        runs = height_matrix_replay(
            pool.coefficients,
            pool.pool,
            precisions=precisions,
            timeout=timeout,
            stack_bytes=stack_bytes,
        )
        rank = stable_height_rank(runs)
        indices = tuple(runs[-1]["subset_indices_one_based"])
        selected = tuple(pool.pool[index - 1] for index in indices)
        return {
            "status": "completed",
            "stable_numerical_rank": rank,
            "precision_runs": list(runs),
            "selected_subset_indices_one_based": list(indices),
            "selected_points": selected,
            "selected_point_sha256": point_digest(selected),
        }
    except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
        return {
            "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
            "error": str(error)[:500],
        }


def parallel_rank_replay(
    pools: Sequence[PointPool],
    *,
    precisions: tuple[int, ...],
    timeout: float,
    stack_bytes: int,
    workers: int,
) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                rank_one_pool,
                pool,
                precisions=precisions,
                timeout=timeout,
                stack_bytes=stack_bytes,
            ): pool.candidate.identifier
            for pool in pools
        }
        for future in as_completed(futures):
            records[futures[future]] = future.result()
    return records


def pool_priority(pool: PointPool, rank: dict[str, Any], root_minus_ids: set[str]) -> tuple[Any, ...]:
    numerical_rank = int(rank.get("stable_numerical_rank", -1))
    return (
        rank.get("status") != "completed",
        -numerical_rank,
        pool.candidate.identifier not in root_minus_ids,
        -len(pool.new_images),
        -Decimal(pool.candidate.score_b2000),
        pool.candidate.identifier,
    )


def parse_positive_ints(value: str) -> tuple[int, ...]:
    try:
        answer = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from error
    if not answer or any(item < 1 for item in answer):
        raise argparse.ArgumentTypeError("integers must be positive")
    return answer


def generated_record(candidate: GeneratedCandidate) -> dict[str, Any]:
    return {
        "candidate_id": candidate.identifier,
        "constructor_parameter_T": rational_to_string(candidate.parameter),
        "height": candidate.height,
        "stratum": candidate.stratum,
        "local_conditions": list(candidate.conditions),
        "radical_proxy": candidate.proxy,
        "exact_homogenized_discriminant_valuations": {
            str(prime): integer_valuation(
                homogenized_discriminant(candidate.parameter), prime
            )
            for prime in SAVING_PRIMES
        },
    }


def exact_record(candidate: ExactCandidate) -> dict[str, Any]:
    return {
        **generated_record(candidate.prefilter.generated),
        "leakage_free_b200": {
            "score": candidate.prefilter.score_b200,
            "good_primes_used": candidate.prefilter.good_primes,
            "bad_primes_skipped": candidate.prefilter.bad_primes,
        },
        "exact_pari_leakage_free_b2000": {
            "score": candidate.score_b2000,
            "good_primes_used": candidate.good_primes,
            "bad_primes_skipped": candidate.bad_primes,
            "design_primes_omitted_count": candidate.omitted_primes,
            "last_prime": candidate.last_prime,
        },
    }


def point_pool_record(
    pool: PointPool, rank: dict[str, Any], *, include_selected_points: bool
) -> dict[str, Any]:
    rank_record = {key: value for key, value in rank.items() if key != "selected_points"}
    record = {
        **exact_record(pool.candidate),
        "point_search": {
            "height_bound": pool.height_bound,
            "status": pool.status,
            "signed_points": pool.signed_points,
            "distinct_quartic_abscissas": pool.signless_points,
            "predeclared_21_section_abscissas_returned": (
                pool.predeclared_abscissas_returned
            ),
            "new_distinct_jacobian_sign_pairs_beyond_21_sections": len(pool.new_images),
            "new_image_sha256": point_digest(pool.new_images),
            "exact_mapping_and_deduplication": True,
            "wall_seconds": pool.wall_seconds,
            "pari_milliseconds": pool.pari_milliseconds,
        },
        "height_rank": rank_record,
    }
    if pool.error is not None:
        record["point_search"]["error"] = pool.error
    if include_selected_points:
        record["exact_selected_points"] = [
            {
                "jacobian_x": rational_to_string(point[0]),
                "jacobian_y": rational_to_string(point[1]),
                "exact_membership_checked": True,
            }
            for point in rank.get("selected_points", ())
        ]
    return record


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace-beam-width", type=int, default=TRACE_BEAM_WIDTH)
    parser.add_argument("--root-beam-width", type=int, default=ROOT_BEAM_WIDTH)
    parser.add_argument("--max-proxy-survivors", type=int, default=MAX_PROXY_SURVIVORS)
    parser.add_argument("--b200-keep", type=int, default=B200_KEEP)
    parser.add_argument("--conductor-keep", type=int, default=CONDUCTOR_KEEP)
    parser.add_argument("--point-keep", type=int, default=POINT_KEEP)
    parser.add_argument("--exact-score-batch", type=int, default=100)
    parser.add_argument("--exact-score-timeout", type=float, default=60.0)
    parser.add_argument("--conductor-timeout", type=float, default=20.0)
    parser.add_argument("--conductor-workers", type=int, default=4)
    parser.add_argument(
        "--stage-heights", type=parse_positive_ints, default=(50_000, 250_000, 1_000_000)
    )
    parser.add_argument("--stage-keeps", type=parse_positive_ints, default=(8, 2))
    parser.add_argument("--stage-timeouts", type=parse_positive_ints, default=(5, 20, 90))
    parser.add_argument("--stage-workers", type=parse_positive_ints, default=(4, 4, 2))
    parser.add_argument("--height-precisions", type=parse_positive_ints, default=(72, 120))
    parser.add_argument("--height-timeout", type=float, default=25.0)
    parser.add_argument("--saturation-timeout", type=float, default=30.0)
    parser.add_argument("--certificate-prime-bound", type=int, default=1_000)
    parser.add_argument("--stack-bytes", type=int, default=512_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts/generated-results/elliptic_nagao_rank20_t5081_neighborhood.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if min(
        args.trace_beam_width,
        args.root_beam_width,
        args.max_proxy_survivors,
        args.b200_keep,
        args.conductor_keep,
        args.point_keep,
        args.exact_score_batch,
        args.conductor_workers,
    ) < 1:
        raise SystemExit("all population and worker bounds must be positive")
    if not (
        args.max_proxy_survivors >= args.b200_keep >= args.conductor_keep >= args.point_keep
    ):
        raise SystemExit("population caps must be nonincreasing")
    if len(args.stage_heights) != 3 or len(args.stage_keeps) != 2:
        raise SystemExit("this lane requires three point heights and two keep counts")
    if len(args.stage_timeouts) != 3 or len(args.stage_workers) != 3:
        raise SystemExit("provide one timeout and worker count per point stage")
    if tuple(sorted(set(args.stage_heights))) != args.stage_heights:
        raise SystemExit("point heights must be strictly increasing")
    if tuple(sorted(set(args.height_precisions))) != args.height_precisions:
        raise SystemExit("height precisions must be strictly increasing")
    if args.stage_keeps[1] > args.stage_keeps[0] or args.stage_keeps[0] > args.point_keep:
        raise SystemExit("point-stage populations must be nonincreasing")
    if min(
        args.exact_score_timeout,
        args.conductor_timeout,
        *args.stage_timeouts,
        args.height_timeout,
        args.saturation_timeout,
    ) <= 0:
        raise SystemExit("subprocess timeouts must be positive")
    if args.stack_bytes < 64_000_000 or args.certificate_prime_bound < 5:
        raise SystemExit("invalid PARI stack or certificate-prime bound")

    all_tables = build_residue_tables(200)
    fingerprint = learn_local_trace_fingerprint(all_tables)
    root_balls = learn_discriminant_root_balls()
    trace_primes = tuple(prime for prime, _, _, _ in fingerprint)
    design_primes = tuple(sorted(set(trace_primes) | set(SAVING_PRIMES)))
    score_tables = {
        prime: table
        for prime, table in all_tables.items()
        if prime not in design_primes
    }

    full_trace, prefix_trace, trace_audit = build_trace_beams(
        fingerprint, all_tables, width=args.trace_beam_width
    )
    beam_strata = build_beam_strata(
        full_trace,
        prefix_trace,
        root_balls,
        root_width=args.root_beam_width,
    )
    candidates, generation_audit = generate_candidates(
        beam_strata,
        proxy_limit=PROXY_LIMIT,
        max_survivors=args.max_proxy_survivors,
    )
    if not candidates:
        raise SystemExit("the exact proxy gate removed the population")
    print(
        f"section7 generated={len(candidates)} raw={generation_audit['raw_parameter_records']} "
        f"proxy_pre_cap={generation_audit['deduplicated_below_proxy_before_cap']}",
        flush=True,
    )

    prefiltered = prefilter_candidates(
        candidates, tables=score_tables, keep=args.b200_keep
    )
    exact = exact_decontaminated_scores(
        prefiltered,
        omitted_primes=design_primes,
        cutoff=EXACT_SCORE_CUTOFF,
        batch_size=args.exact_score_batch,
        timeout=args.exact_score_timeout,
        stack_bytes=args.stack_bytes,
    )
    conductor_population = select_conductor_population(
        exact, keep=args.conductor_keep
    )
    conductor_replays = parallel_conductor_replay(
        conductor_population,
        timeout=args.conductor_timeout,
        stack_bytes=args.stack_bytes,
        workers=args.conductor_workers,
    )
    conductor_by_id = {
        replay.candidate.identifier: replay for replay in conductor_replays
    }
    completed_conductors = [
        replay for replay in conductor_replays if replay.status == "completed"
    ]
    print(
        f"section7 exact_B2000={len(exact)} conductors={len(completed_conductors)}/"
        f"{len(conductor_replays)} subtarget={sum(r.data.get('below_strict_log_conductor_target') is True for r in completed_conductors)} "
        f"root_minus={sum(int(r.data.get('root_number', 0)) == -1 for r in completed_conductors)}",
        flush=True,
    )

    point_population = select_point_population(
        conductor_replays, keep=args.point_keep
    )
    if not point_population:
        raise SystemExit("conductor triage produced no point population")

    # Calibration is processed separately and never consumes a selection slot.
    calibration_prefilter = PrefilterCandidate(
        GeneratedCandidate(
            CALIBRATION_PARAMETER,
            "calibration-only",
            (),
            conductor_radical_proxy(CALIBRATION_PARAMETER),
        ),
        *residue_score(CALIBRATION_PARAMETER, score_tables),
    )
    calibration_exact = exact_decontaminated_scores(
        (calibration_prefilter,),
        omitted_primes=design_primes,
        cutoff=EXACT_SCORE_CUTOFF,
        batch_size=1,
        timeout=args.exact_score_timeout,
        stack_bytes=args.stack_bytes,
    )[0]
    calibration_seed_quartic, calibration_seed_images, calibration_coefficients = (
        exact_predeclared_seeds(CALIBRATION_PARAMETER)
    )
    calibration_height_runs = height_matrix_replay(
        calibration_coefficients,
        calibration_seed_images,
        precisions=args.height_precisions,
        timeout=args.height_timeout,
        stack_bytes=args.stack_bytes,
    )
    calibration_seed_rank = stable_height_rank(calibration_height_runs)
    if calibration_seed_rank != 12:
        raise AssertionError("the section-7 generic seed baseline changed")

    root_minus_ids = {
        replay.candidate.identifier
        for replay in completed_conductors
        if int(replay.data["root_number"]) == -1
        and replay.data.get("below_strict_log_conductor_target") is True
    }
    stages = []
    current = tuple(point_population)
    checkpoints: dict[str, dict[str, Any]] = {}
    checkpoint_pools: dict[str, PointPool] = {}
    final_pools: tuple[PointPool, ...] = ()
    final_ranks: dict[str, dict[str, Any]] = {}

    for stage_index, (height, search_timeout, workers) in enumerate(
        zip(args.stage_heights, args.stage_timeouts, args.stage_workers), start=1
    ):
        pools = parallel_point_search(
            current,
            height_bound=height,
            timeout=float(search_timeout),
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
                    pool, ranks[pool.candidate.identifier], root_minus_ids
                ),
            )
        )
        for pool in ordered:
            rank = ranks[pool.candidate.identifier]
            if rank.get("status") != "completed" or int(rank["stable_numerical_rank"]) < 19:
                continue
            try:
                certificate = finite_reduction_certificate(
                    pool,
                    rank,
                    saturation_timeout=args.saturation_timeout,
                    certificate_prime_bound=args.certificate_prime_bound,
                    stack_bytes=args.stack_bytes,
                )
            except (subprocess.TimeoutExpired, RuntimeError, AssertionError, ValueError) as error:
                certificate = {
                    "status": "timeout" if isinstance(error, subprocess.TimeoutExpired) else "error",
                    "error": str(error)[:500],
                }
            certificate["trigger_height"] = height
            certificate["trigger_stable_numerical_rank"] = rank["stable_numerical_rank"]
            certificate["exact_selected_points"] = [
                [rational_to_string(point[0]), rational_to_string(point[1])]
                for point in rank.get("selected_points", ())
            ]
            checkpoints[pool.candidate.identifier] = certificate
            checkpoint_pools[pool.candidate.identifier] = pool
            replay = conductor_by_id.get(pool.candidate.identifier)
            if replay is None or replay.status != "completed":
                replay = replay_one_conductor(
                    pool.candidate,
                    timeout=args.conductor_timeout,
                    stack_bytes=args.stack_bytes,
                )
                conductor_by_id[pool.candidate.identifier] = replay

        keep_count = (
            args.stage_keeps[stage_index - 1]
            if stage_index <= len(args.stage_keeps)
            else len(ordered)
        )
        retained = ordered[: min(keep_count, len(ordered))]
        stages.append(
            {
                "stage": stage_index,
                "quartic_naive_height_bound": height,
                "population_searched": len(pools),
                "completed_point_searches": sum(pool.status == "completed" for pool in pools),
                "point_search_timeouts": sum(pool.status == "timeout" for pool in pools),
                "point_search_errors": sum(pool.status == "error" for pool in pools),
                "ranked_population": [
                    point_pool_record(
                        pool,
                        ranks[pool.candidate.identifier],
                        include_selected_points=(
                            stage_index == 3
                            or int(ranks[pool.candidate.identifier].get("stable_numerical_rank", 0)) >= 19
                        ),
                    )
                    for pool in ordered
                ],
                "retained_candidate_ids": [pool.candidate.identifier for pool in retained],
            }
        )
        current = tuple(pool.candidate for pool in retained)
        final_pools, final_ranks = ordered, ranks
        print(
            f"section7 H={height} population={len(pools)} "
            f"best_rank={max((int(rank.get('stable_numerical_rank', 0)) for rank in ranks.values()), default=0)}",
            flush=True,
        )

    certified_hits = []
    for identifier, certificate in checkpoints.items():
        replay = conductor_by_id.get(identifier)
        rank = certificate.get("certified_algebraic_rank_lower_bound")
        if certificate.get("status") != "certified" or rank is None:
            continue
        if int(rank) >= 30 or (
            int(rank) >= 21
            and replay is not None
            and replay.status == "completed"
            and replay.data.get("below_strict_log_conductor_target") is True
        ):
            pool = checkpoint_pools[identifier]
            certified_hits.append(
                {
                    "candidate_id": identifier,
                    "constructor_parameter_T": rational_to_string(pool.candidate.parameter),
                    "certified_rank_lower_bound": rank,
                    "conductor": replay.data.get("conductor") if replay else None,
                    "log_conductor": replay.data.get("log_conductor") if replay else None,
                    "root_number": replay.data.get("root_number") if replay else None,
                }
            )

    lead_discriminant = homogenized_discriminant(CALIBRATION_PARAMETER)
    lead_minimal_valuations = {
        str(prime): integer_valuation(EXPECTED_MINIMAL_DISCRIMINANT, prime)
        for prime in (*SAVING_PRIMES, 47)
    }
    bad_calibration = []
    for prime, table in all_tables.items():
        index = projective_index(
            CALIBRATION_PARAMETER.numerator,
            CALIBRATION_PARAMETER.denominator,
            prime,
        )
        if not table[index].good_reduction:
            bad_calibration.append(prime)
    if tuple(bad_calibration) != (5, 7, 13, 17, 23, 47):
        raise AssertionError("the calibration bad-prime fingerprint changed")

    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded section-7 parameter-neighborhood/parity search; numerical "
            "height ranks are triage only"
        ),
        "primary_source": PRIMARY_SOURCE,
        "target": {
            "rank_at_least": 21,
            "strict_log_conductor_upper_bound": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
            "certified_hits": certified_hits,
        },
        "family": {
            "roots_in_source_order": list(ROOTS),
            "primitive_quartic_square_scale": rational_to_string(CONSTRUCTION.quartic_square_scale),
            "learned_invariant_i_ascending": list(INVARIANT_I),
            "learned_invariant_j_ascending": list(INVARIANT_J),
            "primitive_discriminant_polynomial_ascending": list(DISCRIMINANT_POLYNOMIAL),
            "parameter_sign_canonicalization": "T -> abs(T); invariant/discriminant family is even",
        },
        "calibration_only": {
            "constructor_parameter_T": rational_to_string(CALIBRATION_PARAMETER),
            "excluded_before_every_population_selection": True,
            "exact_conductor": str(EXPECTED_CONDUCTOR),
            "exact_log_conductor": "174.249816228548",
            "exact_root_number": 1,
            "exact_minimal_discriminant": str(EXPECTED_MINIMAL_DISCRIMINANT),
            "homogenized_discriminant_valuations": {
                str(prime): integer_valuation(lead_discriminant, prime)
                for prime in SAVING_PRIMES
            },
            "minimal_discriminant_valuations_at_saving_and_denominator_primes": lead_minimal_valuations,
            "bad_primes_through_200": bad_calibration,
            "radical_proxy": calibration_prefilter.generated.proxy,
            "leakage_free_b200": {
                "score": calibration_prefilter.score_b200,
                "good_primes_used": calibration_prefilter.good_primes,
                "bad_primes_skipped": calibration_prefilter.bad_primes,
            },
            "exact_pari_leakage_free_b2000": exact_record(calibration_exact)["exact_pari_leakage_free_b2000"],
            "predeclared_seed_sections": {
                "visible": 12,
                "linear_generic_companions": 6,
                "quadratic_generic_companions": 3,
                "quadratic_companions_symbolically_dependent_over_QT": True,
                "quadratic_companion_relations_in_visible0_through_10_plus_linear_7_over_27_basis": [
                    list(relation) for relation in QUADRATIC_COMPANION_RELATIONS
                ],
                "distinct_jacobian_sign_pairs_at_calibration": len(calibration_seed_images),
                "stable_height_rank_72_120": calibration_seed_rank,
                "precision_runs": list(calibration_height_runs),
            },
        },
        "local_design": {
            "trace_selection_rule": (
                "top five calibration symbols by contribution*log(p/affine-union-size), "
                "with negative trace and 11<=p<=200"
            ),
            "learned_trace_fingerprint": [
                {
                    "prime": prime,
                    "calibration_residue": residue,
                    "calibration_ellap": trace,
                    "affine_union_size": union,
                    "union_threshold": trace,
                }
                for prime, residue, trace, union in fingerprint
            ],
            "discriminant_saving_primes": list(SAVING_PRIMES),
            "automatic_disjoint_root_ball_unions": {
                str(prime): [ball.__dict__ for ball in balls]
                for prime, balls in root_balls.items()
            },
            "all_design_primes_omitted_from_scores": list(design_primes),
            "trace_beam_audit": list(trace_audit),
            "beam_stratum_state_counts": {
                stratum: len(states) for stratum, states in beam_strata.items()
            },
        },
        "population": {
            "positive_primitive_parameters_only": True,
            "excluded_prior_box": {
                "numerator_at_most": OLD_BOX_A_MAX,
                "denominator_at_most": OLD_BOX_B_MAX,
            },
            "proxy_strict_upper_limit": str(PROXY_LIMIT),
            "max_proxy_survivors": args.max_proxy_survivors,
            "smooth_even_denominators": list(smooth_even_denominators()),
            **generation_audit,
        },
        "leakage_free_scoring": {
            "b200_population_scored": len(candidates),
            "b200_retained": len(prefiltered),
            "exact_b2000_rescored": len(exact),
            "score_primes_through_200": list(score_tables),
            "design_primes_omitted": list(design_primes),
            "b200_ranked_tail": [
                {
                    "candidate_id": candidate.identifier,
                    "constructor_parameter_T": rational_to_string(candidate.parameter),
                    "score": candidate.score_b200,
                    "stratum": candidate.generated.stratum,
                    "radical_proxy": candidate.generated.proxy["log_radical_upper_proxy"],
                }
                for candidate in prefiltered[:100]
            ],
            "exact_b2000_ranked_population": [exact_record(candidate) for candidate in exact],
        },
        "conductor_parity_triage": {
            "population_selected": len(conductor_population),
            "selection_union": "24 score leaders + 16 proxy leaders + every stratum leader, capped",
            "completed": len(completed_conductors),
            "subtarget": sum(
                replay.data.get("below_strict_log_conductor_target") is True
                for replay in completed_conductors
            ),
            "subtarget_root_minus": sum(
                replay.data.get("below_strict_log_conductor_target") is True
                and int(replay.data["root_number"]) == -1
                for replay in completed_conductors
            ),
            "records": [
                {
                    "candidate_id": replay.candidate.identifier,
                    "constructor_parameter_T": rational_to_string(replay.candidate.parameter),
                    "status": replay.status,
                    **replay.data,
                    **({"error": replay.error} if replay.error else {}),
                }
                for replay in conductor_by_id.values()
            ],
        },
        "point_triage": {
            "selected_population": [candidate.identifier for candidate in point_population],
            "selection_priority": "subtarget root -1, then subtarget score leaders, plus four unconditional score leaders",
            "predeclared_sections_decontaminated": 21,
            "height_precisions": list(args.height_precisions),
            "stages": stages,
        },
        "exact_checkpoints_stable_numerical_rank_at_least_19": checkpoints,
        "final_frontier": [
            {
                "candidate_id": pool.candidate.identifier,
                "constructor_parameter_T": rational_to_string(pool.candidate.parameter),
                "stable_numerical_rank": final_ranks[pool.candidate.identifier].get("stable_numerical_rank"),
                "new_points_beyond_21_predeclared": len(pool.new_images),
                "exact_score_b2000": pool.candidate.score_b2000,
                "conductor": (
                    conductor_by_id[pool.candidate.identifier].data
                    if pool.candidate.identifier in conductor_by_id
                    else None
                ),
            }
            for pool in final_pools
        ],
        "bounds_and_caveats": {
            "bounded_population_only": True,
            "bounded_quartic_search_only": True,
            "numerical_height_rank_is_not_a_rank_certificate": True,
            "conductor_proxy_is_not_a_conductor_bound": True,
            "target_hit_requires_exact_finite_reduction_and_conductor_replay": True,
            "all_subprocesses_synchronous_with_finite_timeouts": True,
        },
        "reproducibility": {
            "command": REPRODUCING_COMMAND,
            "argv": [shlex.join(sys.argv)],
            "python": platform.python_version(),
            "pari": pari_version(),
            "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        f"wrote {args.output} certified_hits={len(certified_hits)} "
        f"checkpoints={len(checkpoints)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
