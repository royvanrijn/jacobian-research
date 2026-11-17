#!/usr/bin/env python3
"""Search the marked multi-pole triangular peak-reduction conjecture.

This is an exact bounded experiment, not a proof.  A state is a pair of
rational functions on P^1 together with an ordered set of marked boundary
points and a passive conductor pairing.  The peak is the sum of the pole
orders of both functions at every marked boundary point.

The search asks for a state with no lowering complete polynomial shear in
either orientation but with a globally lowering alternating Jung word.
Such a state and word would disprove the proposed peak-reduction statement.
All arithmetic is over Q and every reported local initial coefficient is
computed exactly.
"""
from __future__ import annotations

import argparse
import functools
import itertools
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import sympy as sp


t = sp.Symbol("t")
U = sp.Symbol("U")
Boundary = sp.Rational | None
Polynomial = tuple[Fraction, ...]


def _trim(poly: Polynomial) -> Polynomial:
    end = len(poly)
    while end > 1 and poly[end - 1] == 0:
        end -= 1
    return poly[:end]


def _poly_add(first: Polynomial, second: Polynomial) -> Polynomial:
    length = max(len(first), len(second))
    return _trim(
        tuple(
            (first[index] if index < len(first) else Fraction(0))
            + (second[index] if index < len(second) else Fraction(0))
            for index in range(length)
        )
    )


def _poly_scale(poly: Polynomial, scalar: Fraction) -> Polynomial:
    return _trim(tuple(scalar * coefficient for coefficient in poly))


def _poly_multiply(first: Polynomial, second: Polynomial) -> Polynomial:
    result = [Fraction(0)] * (len(first) + len(second) - 1)
    for first_degree, first_coefficient in enumerate(first):
        for second_degree, second_coefficient in enumerate(second):
            result[first_degree + second_degree] += (
                first_coefficient * second_coefficient
            )
    return _trim(tuple(result))


def _poly_divmod(
    numerator: Polynomial, denominator: Polynomial
) -> tuple[Polynomial, Polynomial]:
    if denominator == (Fraction(0),):
        raise ZeroDivisionError
    remainder = list(numerator)
    quotient = [Fraction(0)] * max(
        1, len(numerator) - len(denominator) + 1
    )
    while len(remainder) >= len(denominator) and any(remainder):
        degree = len(remainder) - len(denominator)
        coefficient = remainder[-1] / denominator[-1]
        quotient[degree] += coefficient
        for index, value in enumerate(denominator):
            remainder[degree + index] -= coefficient * value
        while len(remainder) > 1 and remainder[-1] == 0:
            remainder.pop()
    return _trim(tuple(quotient)), _trim(tuple(remainder))


def _poly_gcd(first: Polynomial, second: Polynomial) -> Polynomial:
    current, following = first, second
    while following != (Fraction(0),):
        _, remainder = _poly_divmod(current, following)
        current, following = following, remainder
    return _poly_scale(current, Fraction(1, 1) / current[-1])


def _poly_evaluate(poly: Polynomial, point: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(poly):
        result = result * point + coefficient
    return result


def _poly_linear_multiplicity(
    poly: Polynomial, point: Fraction
) -> tuple[int, Polynomial]:
    multiplicity = 0
    current = poly
    divisor = (-point, Fraction(1))
    while _poly_evaluate(current, point) == 0:
        current, remainder = _poly_divmod(current, divisor)
        assert remainder == (Fraction(0),)
        multiplicity += 1
    return multiplicity, current


@dataclass(frozen=True)
class RationalFunction:
    """Reduced element of Q(t), stored as ascending coefficient tuples."""

    numerator: Polynomial
    denominator: Polynomial

    @staticmethod
    def make(
        numerator: Polynomial, denominator: Polynomial = (Fraction(1),)
    ) -> "RationalFunction":
        """Construct a quotient without cancelling its boundary factors.

        Search states keep a common-denominator presentation.  Local orders
        subtract numerator and denominator multiplicities exactly, so gcd
        cancellation is unnecessary in the hot path and would dominate the
        degree-27 runs.
        """
        numerator = _trim(numerator)
        denominator = _trim(denominator)
        if denominator == (Fraction(0),):
            raise ZeroDivisionError
        if numerator == (Fraction(0),):
            return RationalFunction((Fraction(0),), (Fraction(1),))
        return RationalFunction(numerator, denominator)

    @staticmethod
    def integer(value: int | Fraction) -> "RationalFunction":
        return RationalFunction.make((Fraction(value),))

    def __add__(
        self, other: "RationalFunction" | int | Fraction
    ) -> "RationalFunction":
        if not isinstance(other, RationalFunction):
            other = RationalFunction.integer(other)
        return RationalFunction.make(
            _poly_add(
                _poly_multiply(self.numerator, other.denominator),
                _poly_multiply(other.numerator, self.denominator),
            ),
            _poly_multiply(self.denominator, other.denominator),
        )

    __radd__ = __add__

    def __neg__(self) -> "RationalFunction":
        return RationalFunction(
            _poly_scale(self.numerator, Fraction(-1)),
            self.denominator,
        )

    def __sub__(
        self, other: "RationalFunction" | int | Fraction
    ) -> "RationalFunction":
        return self + (-other if isinstance(other, RationalFunction) else -other)

    def __mul__(
        self, other: "RationalFunction" | int | Fraction
    ) -> "RationalFunction":
        if not isinstance(other, RationalFunction):
            other = RationalFunction.integer(other)
        return RationalFunction.make(
            _poly_multiply(self.numerator, other.numerator),
            _poly_multiply(self.denominator, other.denominator),
        )

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> "RationalFunction":
        if exponent < 0:
            return RationalFunction.make(
                _poly_power(self.denominator, -exponent),
                _poly_power(self.numerator, -exponent),
            )
        return RationalFunction.make(
            _poly_power(self.numerator, exponent),
            _poly_power(self.denominator, exponent),
        )

    def to_sympy(self) -> sp.Expr:
        def convert(poly: Polynomial) -> sp.Expr:
            return sum(
                sp.Rational(value.numerator, value.denominator) * t**degree
                for degree, value in enumerate(poly)
            )

        return sp.cancel(convert(self.numerator) / convert(self.denominator))

    def __str__(self) -> str:
        return str(self.to_sympy())


def equivalent(
    first: RationalFunction, second: RationalFunction
) -> bool:
    """Exact equality in Q(t), independent of quotient presentation."""
    return _poly_multiply(
        first.numerator, second.denominator
    ) == _poly_multiply(second.numerator, first.denominator)


def _poly_power(poly: Polynomial, exponent: int) -> Polynomial:
    result = (Fraction(1),)
    current = poly
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = _poly_multiply(result, current)
        remaining >>= 1
        if remaining:
            current = _poly_multiply(current, current)
    return result


@dataclass(frozen=True)
class LocalInitial:
    """Order, pole order, and initial coefficient in the chosen parameter."""

    order: int
    pole: int
    residue: Fraction


@dataclass(frozen=True)
class Factor:
    """The triangular factor target <- target + polynomial(base)."""

    target: int
    polynomial: sp.Expr

    @property
    def degree(self) -> int:
        return int(sp.Poly(self.polynomial, U).degree())

    def inverse(self) -> "Factor":
        return Factor(self.target, sp.expand(-self.polynomial))


@functools.lru_cache(maxsize=None)
def local_initial(
    expression: RationalFunction, boundary: Boundary
) -> LocalInitial:
    """Return the exact valuation and initial coefficient at one boundary."""
    numerator = expression.numerator
    denominator = expression.denominator
    if numerator == (Fraction(0),):
        raise ValueError("the zero function has no local initial form")

    if boundary is None:
        order = len(denominator) - len(numerator)
        residue = numerator[-1] / denominator[-1]
    else:
        rational_boundary = Fraction(int(boundary.p), int(boundary.q))
        num_order, num_unit = _poly_linear_multiplicity(
            numerator, rational_boundary
        )
        den_order, den_unit = _poly_linear_multiplicity(
            denominator, rational_boundary
        )
        order = num_order - den_order
        residue = (
            _poly_evaluate(num_unit, rational_boundary)
            / _poly_evaluate(den_unit, rational_boundary)
        )

    return LocalInitial(order, max(0, -order), residue)


@functools.lru_cache(maxsize=None)
def boundary_ledger(
    first: RationalFunction,
    second: RationalFunction,
    boundaries: tuple[Boundary, ...],
) -> tuple[tuple[LocalInitial, LocalInitial], ...]:
    return tuple(
        (local_initial(first, point), local_initial(second, point))
        for point in boundaries
    )


@functools.lru_cache(maxsize=None)
def peak(
    first: RationalFunction,
    second: RationalFunction,
    boundaries: tuple[Boundary, ...],
) -> int:
    return sum(
        first_local.pole + second_local.pole
        for first_local, second_local in boundary_ledger(
            first, second, boundaries
        )
    )


def residue_character(
    first: RationalFunction,
    second: RationalFunction,
    boundaries: tuple[Boundary, ...],
) -> tuple[tuple[str, str], ...]:
    """The two exact initial-form coefficients at each marked boundary."""
    return tuple(
        (str(first_local.residue), str(second_local.residue))
        for first_local, second_local in boundary_ledger(
            first, second, boundaries
        )
    )


def _evaluate_polynomial(
    polynomial: sp.Expr, base: RationalFunction
) -> RationalFunction:
    """Evaluate a small polynomial without creating a nested substitution."""
    terms = sp.Poly(polynomial, U, domain=sp.QQ).terms()
    return sum(
        (
            Fraction(int(coefficient.p), int(coefficient.q))
            * base**degree[0]
            for degree, coefficient in terms
        ),
        start=RationalFunction.integer(0),
    )


@functools.lru_cache(maxsize=None)
def apply_factor(
    first: RationalFunction,
    second: RationalFunction,
    factor: Factor,
) -> tuple[RationalFunction, RationalFunction]:
    base = second if factor.target == 0 else first
    changed = (
        (first if factor.target == 0 else second)
        + _evaluate_polynomial(factor.polynomial, base)
    )
    if changed.numerator == (Fraction(0),):
        # A constant translation is invisible to every marked pole and
        # restores a nonzero restriction.
        changed = RationalFunction.integer(1)
    return (changed, second) if factor.target == 0 else (first, changed)


def _forced_cancellations(
    base: RationalFunction,
    target: RationalFunction,
    boundaries: tuple[Boundary, ...],
) -> tuple[tuple[int, sp.Expr], ...]:
    """Candidate terms forced by a tied leading pole at some boundary."""
    candidates: set[tuple[int, sp.Expr]] = set()
    for base_local, target_local in zip(
        (
            local_initial(base, point)
            for point in boundaries
        ),
        (
            local_initial(target, point)
            for point in boundaries
        ),
        strict=True,
    ):
        if base_local.pole == 0 or target_local.pole == 0:
            continue
        if target_local.pole % base_local.pole:
            continue
        degree = target_local.pole // base_local.pole
        if degree < 1:
            continue
        coefficient = -target_local.residue / base_local.residue**degree
        candidates.add((degree, coefficient))
    return tuple(
        sorted(
            candidates,
            key=lambda item: (
                -item[0],
                item[1],
            ),
        )
    )


def lowering_complete_shears(
    first: RationalFunction,
    second: RationalFunction,
    boundaries: tuple[Boundary, ...],
) -> tuple[Factor, ...]:
    """Exhaust endpoint-relevant complete lowering polynomial shears."""
    initial_peak = peak(first, second, boundaries)
    found: dict[tuple[int, sp.Expr], Factor] = {}

    for target_index, (base, target) in enumerate(
        ((second, first), (first, second))
    ):
        visited: set[tuple[RationalFunction, int]] = set()

        def explore(
            remainder: RationalFunction,
            maximum_degree: int | None,
            polynomial: sp.Expr,
        ) -> None:
            state = (remainder, maximum_degree or 0)
            if state in visited:
                return
            visited.add(state)
            for degree, coefficient in _forced_cancellations(
                base, remainder, boundaries
            ):
                if maximum_degree is not None and degree >= maximum_degree:
                    continue
                sympy_coefficient = sp.Rational(
                    coefficient.numerator, coefficient.denominator
                )
                next_polynomial = sp.expand(
                    polynomial + sympy_coefficient * U**degree
                )
                next_remainder = remainder + coefficient * base**degree
                if next_remainder.numerator == (Fraction(0),):
                    next_remainder = RationalFunction.integer(1)
                if target_index == 0:
                    candidate = (next_remainder, second)
                else:
                    candidate = (first, next_remainder)
                if peak(*candidate, boundaries) < initial_peak:
                    factor = Factor(target_index, next_polynomial)
                    found[(target_index, next_polynomial)] = factor
                explore(
                    next_remainder,
                    degree,
                    next_polynomial,
                )

        explore(target, None, sp.Integer(0))

    return tuple(
        sorted(
            found.values(),
            key=lambda factor: (
                factor.degree,
                factor.target,
                sp.default_sort_key(factor.polynomial),
            ),
        )
    )


def _boundary_name(point: Boundary) -> str:
    return "infinity" if point is None else str(point)


def trace_word(
    first: RationalFunction,
    second: RationalFunction,
    boundaries: tuple[Boundary, ...],
    pairing: tuple[tuple[int, int], ...],
    word: tuple[Factor, ...],
) -> dict[str, object]:
    """Return the complete exact per-factor marked ledger."""
    states: list[dict[str, object]] = []
    current = (first, second)
    for index in range(len(word) + 1):
        ledger = boundary_ledger(*current, boundaries)
        states.append(
            {
                "index": index,
                "first": str(current[0]),
                "second": str(current[1]),
                "peak": peak(*current, boundaries),
                "boundary_valuations": [
                    {
                        "boundary": _boundary_name(point),
                        "orders": [row[0].order, row[1].order],
                        "poles": [row[0].pole, row[1].pole],
                    }
                    for point, row in zip(boundaries, ledger, strict=True)
                ],
                "residue_character": [
                    {
                        "boundary": _boundary_name(point),
                        "initial_coefficients": [
                            str(row[0].residue),
                            str(row[1].residue),
                        ],
                    }
                    for point, row in zip(boundaries, ledger, strict=True)
                ],
                "conductor_pairing": [list(pair) for pair in pairing],
            }
        )
        if index < len(word):
            current = apply_factor(*current, word[index])

    return {
        "polydegree": [factor.degree for factor in word],
        "factors": [
            {
                "target": factor.target,
                "polynomial": str(factor.polynomial),
                "peak_before": states[index]["peak"],
                "peak_after": states[index + 1]["peak"],
            }
            for index, factor in enumerate(word)
        ],
        "states": states,
    }


def coefficient_polynomials(
    maximum_degree: int,
    complete: bool,
    include_linear: bool,
) -> tuple[sp.Expr, ...]:
    """Normalized bounded factor polynomials in increasing degree."""
    polynomials: list[sp.Expr] = []
    for degree in range(1 if include_linear else 2, maximum_degree + 1):
        if complete:
            lower_rows: Iterable[tuple[int, ...]] = itertools.product(
                (-1, 0, 1), repeat=degree - 1
            )
        else:
            lower_rows = ((0,) * (degree - 1),)
        for leading in (-1, 1):
            for lower in lower_rows:
                polynomial = leading * U**degree
                polynomial += sum(
                    lower[index - 1] * U**index
                    for index in range(1, degree)
                )
                polynomials.append(sp.expand(polynomial))
    return tuple(polynomials)


def verify_local_increment_inequality(
    maximum_pole: int = 12,
    maximum_degree: int = 6,
) -> int:
    """Bounded replay of the uniform two-factor pole-change lemma."""

    def possible_outputs(
        target_pole: int, base_pole: int, degree: int
    ) -> tuple[int, ...]:
        polynomial_pole = degree * base_pole
        if target_pole != polynomial_pole:
            return (max(target_pole, polynomial_pole),)
        # Equal leading poles may cancel to any lower nonnegative order.
        return tuple(range(target_pole + 1))

    cases = 0
    for a, b, degree, next_degree in itertools.product(
        range(maximum_pole + 1),
        range(maximum_pole + 1),
        range(2, maximum_degree + 1),
        range(2, maximum_degree + 1),
    ):
        for a_prime in possible_outputs(a, b, degree):
            for b_prime in possible_outputs(
                b, a_prime, next_degree
            ):
                assert b_prime - b >= a_prime - a
                cases += 1
    return cases


def seed_functions(
    atoms: tuple[RationalFunction, ...],
    signed: bool,
    maximum_terms: int | None,
) -> tuple[RationalFunction, ...]:
    coefficients = (-1, 0, 1) if signed else (0, 1)
    rows = []
    for row in itertools.product(coefficients, repeat=len(atoms)):
        if not any(row):
            continue
        if maximum_terms is not None and sum(value != 0 for value in row) > (
            maximum_terms
        ):
            continue
        rows.append(
            sum(
                (c * atom for c, atom in zip(row, atoms)),
                start=RationalFunction.integer(0),
            )
        )
    return tuple(dict.fromkeys(rows))


def alternating_words(
    length: int,
    polynomials: tuple[sp.Expr, ...],
) -> Iterable[tuple[Factor, ...]]:
    for degrees_and_coefficients in itertools.product(
        polynomials, repeat=length
    ):
        for initial_target in (0, 1):
            yield tuple(
                Factor((initial_target + index) % 2, polynomial)
                for index, polynomial in enumerate(
                    degrees_and_coefficients
                )
            )


def orbit_endpoint_search(
    functions: tuple[RationalFunction, ...],
    boundaries: tuple[Boundary, ...],
    pairing: tuple[tuple[int, int], ...],
    polynomials: tuple[sp.Expr, ...],
    maximum_length: int,
) -> dict[str, object]:
    """Seek a terminal endpoint whose inverse word lowers globally.

    Generating the candidate as a forward orbit reaches states far outside
    the bounded seed-function box.  If the final forward factor increased
    height, its inverse is already a lowering complete shear, so only
    nonincreasing last factors need the expensive terminal test.
    """
    summary: dict[str, object] = {
        "base_pairs": 0,
        "forward_words_tested": 0,
        "endpoints_above_base": 0,
        "nonincreasing_last_factors": 0,
        "complete_shear_terminal_endpoints": 0,
        "counterexample": None,
    }
    base_pairs = [
        (first, second)
        for first, second in itertools.product(functions, repeat=2)
        if not equivalent(first, second)
        and not equivalent(first, -second)
    ]
    summary["base_pairs"] = len(base_pairs)

    for length in range(2, maximum_length + 1):
        words = sorted(
            alternating_words(length, polynomials),
            key=lambda word: (
                int(sp.prod(factor.degree for factor in word)),
                tuple(factor.degree for factor in word),
                tuple(
                    sp.default_sort_key(factor.polynomial)
                    for factor in word
                ),
                word[0].target,
            ),
        )
        for word in words:
            inverse_word = tuple(
                factor.inverse() for factor in reversed(word)
            )
            for base_first, base_second in base_pairs:
                summary["forward_words_tested"] = int(
                    summary["forward_words_tested"]
                ) + 1
                states = [(base_first, base_second)]
                for factor in word:
                    states.append(apply_factor(*states[-1], factor))
                endpoint = states[-1]
                base_peak = peak(base_first, base_second, boundaries)
                endpoint_peak = peak(*endpoint, boundaries)
                if endpoint_peak <= base_peak:
                    continue
                summary["endpoints_above_base"] = int(
                    summary["endpoints_above_base"]
                ) + 1
                if endpoint_peak > peak(*states[-2], boundaries):
                    continue
                summary["nonincreasing_last_factors"] = int(
                    summary["nonincreasing_last_factors"]
                ) + 1

                lowering = lowering_complete_shears(
                    *endpoint, boundaries
                )
                if lowering:
                    continue
                summary["complete_shear_terminal_endpoints"] = int(
                    summary["complete_shear_terminal_endpoints"]
                ) + 1

                recovered = endpoint
                for factor in inverse_word:
                    recovered = apply_factor(*recovered, factor)
                if not (
                    equivalent(recovered[0], base_first)
                    and equivalent(recovered[1], base_second)
                ):
                    continue
                summary["counterexample"] = {
                    "forward_trace": trace_word(
                        base_first,
                        base_second,
                        boundaries,
                        pairing,
                        word,
                    ),
                    "globally_lowering_inverse_trace": trace_word(
                        *endpoint,
                        boundaries,
                        pairing,
                        inverse_word,
                    ),
                    "endpoint_complete_lowering_shears": [],
                }
                return summary
    return summary


def search_chart(
    name: str,
    functions: tuple[RationalFunction, ...],
    boundaries: tuple[Boundary, ...],
    pairing: tuple[tuple[int, int], ...],
    polynomials: tuple[sp.Expr, ...],
    maximum_length: int,
    scan_all: bool,
    search_orbit_endpoints: bool,
    orbit_only: bool,
) -> dict[str, object]:
    """Search in length/polydegree order and return the first counterexample."""
    summary: dict[str, object] = {
        "chart": name,
        "boundaries": [_boundary_name(point) for point in boundaries],
        "conductor_pairing": [list(pair) for pair in pairing],
        "seed_function_count": len(functions),
        "ordered_seed_pair_count": len(functions) ** 2,
        "locally_terminal_seed_pairs": 0,
        "words_tested_from_terminal_pairs": 0,
        "globally_lowering_words": 0,
        "all_words_tested": 0,
        "all_globally_lowering_words": 0,
        "delayed_globally_lowering_words": 0,
        "delayed_words_with_initial_complete_reduction": 0,
        "delayed_reduced_jung_words": 0,
        "delayed_polydegree_counts": {},
        "peak_delta_sequences": {},
        "largest_delayed_peak_excess": None,
        "largest_delayed_peak_trace": None,
        "largest_polydegree_delayed_trace": None,
        "counterexample": None,
        "orbit_endpoint_search": None,
    }
    if orbit_only:
        summary["orbit_endpoint_search"] = orbit_endpoint_search(
            functions,
            boundaries,
            pairing,
            polynomials,
            maximum_length,
        )
        summary["counterexample"] = summary["orbit_endpoint_search"][
            "counterexample"
        ]
        return summary

    terminal_pairs: list[tuple[RationalFunction, RationalFunction]] = []
    all_pairs: list[
        tuple[
            RationalFunction,
            RationalFunction,
            tuple[Factor, ...],
        ]
    ] = []
    for first, second in itertools.product(functions, repeat=2):
        if first == second or first == -second:
            continue
        lowering_shears = lowering_complete_shears(
            first, second, boundaries
        )
        all_pairs.append((first, second, lowering_shears))
        if not lowering_shears:
            terminal_pairs.append((first, second))
    summary["locally_terminal_seed_pairs"] = len(terminal_pairs)

    for length in range(2, maximum_length + 1):
        words = sorted(
            alternating_words(length, polynomials),
            key=lambda word: (
                sp.prod(factor.degree for factor in word),
                tuple(factor.degree for factor in word),
                tuple(
                    sp.default_sort_key(factor.polynomial)
                    for factor in word
                ),
                word[0].target,
            ),
        )
        for first, second in terminal_pairs:
            initial_peak = peak(first, second, boundaries)
            for word in words:
                summary["words_tested_from_terminal_pairs"] = int(
                    summary["words_tested_from_terminal_pairs"]
                ) + 1
                current = (first, second)
                for factor in word:
                    current = apply_factor(*current, factor)
                if peak(*current, boundaries) >= initial_peak:
                    continue
                summary["globally_lowering_words"] = int(
                    summary["globally_lowering_words"]
                ) + 1
                summary["counterexample"] = {
                    "initial_complete_lowering_shears": [],
                    "trace": trace_word(
                        first,
                        second,
                        boundaries,
                        pairing,
                        word,
                    ),
                }
                return summary

    if scan_all:
        delta_counts: dict[str, int] = {}
        delayed_polydegree_counts: dict[str, int] = {}
        largest_excess = -1
        largest_polydegree = -1
        for length in range(2, maximum_length + 1):
            words = sorted(
                alternating_words(length, polynomials),
                key=lambda word: (
                    sp.prod(factor.degree for factor in word),
                    tuple(factor.degree for factor in word),
                    tuple(
                        sp.default_sort_key(factor.polynomial)
                        for factor in word
                    ),
                    word[0].target,
                ),
            )
            for first, second, lowering_shears in all_pairs:
                initial_peak = peak(first, second, boundaries)
                for word in words:
                    summary["all_words_tested"] = int(
                        summary["all_words_tested"]
                    ) + 1
                    current = (first, second)
                    peaks = [initial_peak]
                    for factor in word:
                        current = apply_factor(*current, factor)
                        peaks.append(peak(*current, boundaries))
                    if peaks[-1] >= initial_peak:
                        continue
                    summary["all_globally_lowering_words"] = int(
                        summary["all_globally_lowering_words"]
                    ) + 1
                    deltas = tuple(value - initial_peak for value in peaks)
                    delta_key = ",".join(str(value) for value in deltas)
                    delta_counts[delta_key] = delta_counts.get(delta_key, 0) + 1
                    if peaks[1] < initial_peak:
                        continue
                    summary["delayed_globally_lowering_words"] = int(
                        summary["delayed_globally_lowering_words"]
                    ) + 1
                    polydegree = tuple(factor.degree for factor in word)
                    polydegree_key = ",".join(
                        str(value) for value in polydegree
                    )
                    delayed_polydegree_counts[polydegree_key] = (
                        delayed_polydegree_counts.get(polydegree_key, 0) + 1
                    )
                    if all(value >= 2 for value in polydegree):
                        summary["delayed_reduced_jung_words"] = int(
                            summary["delayed_reduced_jung_words"]
                        ) + 1
                    if lowering_shears:
                        summary[
                            "delayed_words_with_initial_complete_reduction"
                        ] = int(
                            summary[
                                "delayed_words_with_initial_complete_reduction"
                            ]
                        ) + 1
                    excess = max(peaks) - initial_peak
                    if excess > largest_excess:
                        largest_excess = excess
                        local_reductions = []
                        for factor in lowering_shears:
                            reduced = apply_factor(first, second, factor)
                            local_reductions.append(
                                {
                                    "target": factor.target,
                                    "polynomial": str(factor.polynomial),
                                    "peak_before": initial_peak,
                                    "peak_after": peak(
                                        *reduced, boundaries
                                    ),
                                }
                            )
                        summary["largest_delayed_peak_excess"] = excess
                        summary["largest_delayed_peak_trace"] = {
                            "initial_complete_reductions": local_reductions,
                            "word": trace_word(
                                first,
                                second,
                                boundaries,
                                pairing,
                                word,
                            ),
                        }
                    polydegree_product = int(sp.prod(polydegree))
                    if polydegree_product > largest_polydegree:
                        largest_polydegree = polydegree_product
                        summary["largest_polydegree_delayed_trace"] = (
                            trace_word(
                                first,
                                second,
                                boundaries,
                                pairing,
                                word,
                            )
                        )
        summary["peak_delta_sequences"] = dict(
            sorted(delta_counts.items())
        )
        summary["delayed_polydegree_counts"] = dict(
            sorted(delayed_polydegree_counts.items())
        )
    if search_orbit_endpoints:
        summary["orbit_endpoint_search"] = orbit_endpoint_search(
            functions,
            boundaries,
            pairing,
            polynomials,
            maximum_length,
        )
        orbit_counterexample = summary["orbit_endpoint_search"][
            "counterexample"
        ]
        if orbit_counterexample is not None:
            summary["counterexample"] = orbit_counterexample
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-degree", type=int, default=3)
    parser.add_argument("--max-length", type=int, default=3)
    parser.add_argument(
        "--complete-factors",
        action="store_true",
        help="include all lower coefficients in {-1,0,1}",
    )
    parser.add_argument(
        "--include-linear",
        action="store_true",
        help="include affine triangular factors of degree one",
    )
    parser.add_argument(
        "--signed-seeds",
        action="store_true",
        help="use seed coefficients in {-1,0,1}, not just {0,1}",
    )
    parser.add_argument(
        "--extended-seeds",
        action="store_true",
        help="include pole orders two in the seed basis",
    )
    parser.add_argument(
        "--four-pole",
        action="store_true",
        help="also test the marked boundaries -1, 0, 1, infinity",
    )
    parser.add_argument(
        "--max-seed-terms",
        type=int,
        default=None,
        help="limit the number of nonzero basis terms in each seed function",
    )
    parser.add_argument(
        "--scan-all",
        action="store_true",
        help="also aggregate every globally lowering word, not only terminals",
    )
    parser.add_argument(
        "--orbit-endpoints",
        action="store_true",
        help="generate high-complexity endpoints and test their inverse words",
    )
    parser.add_argument(
        "--orbit-only",
        action="store_true",
        help="run only the generated-orbit endpoint search",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "artifacts/generated-results/marked_multipole_peak_search.json"
        ),
    )
    args = parser.parse_args()

    local_increment_cases = verify_local_increment_inequality()
    factor_polynomials = coefficient_polynomials(
        args.max_degree, args.complete_factors, args.include_linear
    )
    inverse_t = RationalFunction.make(
        (Fraction(1),), (Fraction(0), Fraction(1))
    )
    inverse_t_squared = inverse_t**2
    inverse_t_minus_one = RationalFunction.make(
        (Fraction(1),), (Fraction(-1), Fraction(1))
    )
    inverse_t_minus_one_squared = inverse_t_minus_one**2
    inverse_t_plus_one = RationalFunction.make(
        (Fraction(1),), (Fraction(1), Fraction(1))
    )
    inverse_t_plus_one_squared = inverse_t_plus_one**2
    constant_one = RationalFunction.integer(1)
    parameter_t = RationalFunction.make((Fraction(0), Fraction(1)))
    parameter_t_squared = parameter_t**2
    two_pole_atoms = (
        (
            inverse_t_squared,
            inverse_t,
            constant_one,
            parameter_t,
            parameter_t_squared,
        )
        if args.extended_seeds
        else (inverse_t, constant_one, parameter_t)
    )
    three_pole_atoms = (
        (
            inverse_t_squared,
            inverse_t,
            inverse_t_minus_one_squared,
            inverse_t_minus_one,
            constant_one,
            parameter_t,
            parameter_t_squared,
        )
        if args.extended_seeds
        else (inverse_t, inverse_t_minus_one, parameter_t)
    )
    four_pole_atoms = (
        (
            inverse_t_plus_one_squared,
            inverse_t_plus_one,
            inverse_t_squared,
            inverse_t,
            inverse_t_minus_one_squared,
            inverse_t_minus_one,
            parameter_t,
            parameter_t_squared,
        )
        if args.extended_seeds
        else (
            inverse_t_plus_one,
            inverse_t,
            inverse_t_minus_one,
            parameter_t,
        )
    )
    charts = [
        (
            "two-pole Laurent",
            two_pole_atoms,
            (sp.Rational(0), None),
            ((0, 1),),
        ),
        (
            "three-pole rational",
            three_pole_atoms,
            (sp.Rational(0), sp.Rational(1), None),
            ((0, 1),),
        ),
    ]
    if args.four_pole:
        charts.append(
            (
                "four-pole rational",
                four_pole_atoms,
                (
                    sp.Rational(-1),
                    sp.Rational(0),
                    sp.Rational(1),
                    None,
                ),
                ((0, 2), (1, 3)),
            )
        )

    results = []
    for name, atoms, boundaries, pairing in charts:
        functions = seed_functions(
            atoms, args.signed_seeds, args.max_seed_terms
        )
        result = search_chart(
            name,
            functions,
            boundaries,
            pairing,
            factor_polynomials,
            args.max_length,
            args.scan_all,
            args.orbit_endpoints or args.orbit_only,
            args.orbit_only,
        )
        results.append(result)
        if result["counterexample"] is not None:
            break

    payload = {
        "status": (
            "counterexample"
            if any(result["counterexample"] for result in results)
            else "no counterexample in bounded search"
        ),
        "scope": {
            "max_factor_degree": args.max_degree,
            "max_word_length": args.max_length,
            "factor_lower_coefficients": (
                "all in {-1,0,1}"
                if args.complete_factors
                else "zero (monomial factors)"
            ),
            "linear_triangular_factors": args.include_linear,
            "seed_coefficients": (
                "{-1,0,1}" if args.signed_seeds else "{0,1}"
            ),
            "extended_seed_basis": args.extended_seeds,
            "maximum_seed_terms": args.max_seed_terms,
            "four_pole_chart": args.four_pole,
            "word_order": (
                "length, total polydegree, polydegree tuple, coefficients, "
                "initial orientation"
            ),
            "verified_local_increment_cases": local_increment_cases,
        },
        "invariant_tested": (
            "a globally lowering alternating word cannot start at a state "
            "with no complete lowering triangular shear"
        ),
        "charts": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
