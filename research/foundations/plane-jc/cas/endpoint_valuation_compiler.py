"""Exact residue-aware endpoint reduction for two-ended Laurent charts.

The module exhausts both strictly decreasing monomial shears and complete
polynomial shears assembled from descending-degree cancellation chains.
Prefixes of a lowering polynomial shear may preserve or increase height.
The resulting strict-descent normal form is not claimed to be minimal under
the full polynomial automorphism group; that requires a separate marked
multi-pole peak-reduction theorem for alternating triangular directions.
"""
from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class EndpointInitial:
    """Pole order and leading residue at t=0 or t=infinity."""

    pole: int
    residue: sp.Expr | None


@dataclass(frozen=True)
class ShearStep:
    """One target <- target - coefficient*base^degree reduction."""

    target: int
    degree: int
    coefficient: sp.Expr
    height_before: int
    height_after: int


@dataclass(frozen=True)
class PolynomialShear:
    """One target <- target - polynomial(base) strict reduction."""

    target: int
    polynomial: sp.Expr
    terms: tuple[tuple[int, sp.Expr], ...]
    height_before: int
    height_after: int


@dataclass(frozen=True)
class ReductionResult:
    """A monotone elementary normal form and its reduction certificate."""

    first: sp.Expr
    second: sp.Expr
    height: int
    steps: tuple[ShearStep, ...]


@dataclass(frozen=True)
class PolynomialReductionResult:
    """A strict-descent normal form using complete polynomial shears."""

    first: sp.Expr
    second: sp.Expr
    height: int
    steps: tuple[PolynomialShear, ...]


def _laurent_coefficients(
    expression: sp.Expr, parameter: sp.Symbol
) -> dict[int, sp.Expr]:
    coefficients: dict[int, sp.Expr] = {}
    for term in sp.Add.make_args(sp.expand(expression)):
        coefficient, exponent = term.as_coeff_exponent(parameter)
        if parameter in coefficient.free_symbols:
            raise ValueError("expression is not a Laurent polynomial")
        integer_exponent = int(exponent)
        if exponent != integer_exponent:
            raise ValueError("Laurent exponents must be integral")
        coefficients[integer_exponent] = sp.simplify(
            coefficients.get(integer_exponent, 0) + coefficient
        )
    return {
        exponent: coefficient
        for exponent, coefficient in coefficients.items()
        if coefficient != 0
    }


def endpoint_initials(
    expression: sp.Expr, parameter: sp.Symbol
) -> tuple[EndpointInitial, EndpointInitial]:
    """Return initial pole data at (t=0, t=infinity)."""
    coefficients = _laurent_coefficients(expression, parameter)
    if not coefficients:
        raise ValueError("the zero Laurent polynomial has no endpoint initial")

    minimum = min(coefficients)
    maximum = max(coefficients)
    zero_pole = max(0, -minimum)
    infinity_pole = max(0, maximum)
    return (
        EndpointInitial(
            zero_pole,
            coefficients[minimum] if zero_pole else None,
        ),
        EndpointInitial(
            infinity_pole,
            coefficients[maximum] if infinity_pole else None,
        ),
    )


def pole_matrix(
    first: sp.Expr, second: sp.Expr, parameter: sp.Symbol
) -> tuple[tuple[int, int], tuple[int, int]]:
    """Rows are the pole pairs at t=0 and t=infinity."""
    first_initials = endpoint_initials(first, parameter)
    second_initials = endpoint_initials(second, parameter)
    return (
        (first_initials[0].pole, second_initials[0].pole),
        (first_initials[1].pole, second_initials[1].pole),
    )


def pole_height(
    first: sp.Expr, second: sp.Expr, parameter: sp.Symbol
) -> int:
    """Sum of the four nonnegative entries of the two-endpoint pole matrix."""
    return sum(sum(row) for row in pole_matrix(first, second, parameter))


def _candidate_shears(
    base: sp.Expr,
    target: sp.Expr,
    parameter: sp.Symbol,
) -> tuple[tuple[int, sp.Expr], ...]:
    """Finite degrees/residues that can cancel a current target pole."""
    base_initials = endpoint_initials(base, parameter)
    target_initials = endpoint_initials(target, parameter)
    candidates: set[tuple[int, sp.Expr]] = set()

    for base_initial, target_initial in zip(
        base_initials, target_initials, strict=True
    ):
        if base_initial.pole == 0 or target_initial.pole == 0:
            continue
        if target_initial.pole % base_initial.pole:
            continue
        degree = target_initial.pole // base_initial.pole
        if degree < 1:
            continue
        assert base_initial.residue is not None
        assert target_initial.residue is not None
        coefficient = sp.simplify(
            target_initial.residue / base_initial.residue**degree
        )
        candidates.add((degree, coefficient))

    return tuple(
        sorted(candidates, key=lambda item: (item[0], sp.default_sort_key(item[1])))
    )


def reducing_shears(
    first: sp.Expr,
    second: sp.Expr,
    parameter: sp.Symbol,
) -> tuple[tuple[sp.Expr, sp.Expr, ShearStep], ...]:
    """Enumerate every elementary shear that strictly lowers pole height."""
    before = pole_height(first, second, parameter)
    moves: list[tuple[sp.Expr, sp.Expr, ShearStep]] = []

    for target_index, (base, target) in enumerate(
        ((second, first), (first, second))
    ):
        for degree, coefficient in _candidate_shears(
            base, target, parameter
        ):
            reduced_target = sp.expand(target - coefficient * base**degree)
            if reduced_target == 0:
                # A constant translation, invisible to endpoint poles,
                # replaces the zero restriction by the nonzero constant 1.
                reduced_target = sp.Integer(1)
            if target_index == 0:
                new_first, new_second = reduced_target, second
            else:
                new_first, new_second = first, reduced_target
            after = pole_height(new_first, new_second, parameter)
            if after >= before:
                continue
            moves.append(
                (
                    new_first,
                    new_second,
                    ShearStep(
                        target=target_index,
                        degree=degree,
                        coefficient=coefficient,
                        height_before=before,
                        height_after=after,
                    ),
                )
            )

    return tuple(
        sorted(
            moves,
            key=lambda move: (
                move[2].height_after,
                pole_matrix(move[0], move[1], parameter),
                sp.default_sort_key(move[0]),
                sp.default_sort_key(move[1]),
            ),
        )
    )


def _polynomial_shears_one_orientation(
    base: sp.Expr,
    target: sp.Expr,
    parameter: sp.Symbol,
    target_index: int,
) -> tuple[PolynomialShear, ...]:
    """Exhaust lowering P(base) shears via descending cancellation degrees."""
    initial_height = pole_height(base, target, parameter)
    found: dict[tuple[sp.Expr, sp.Expr], PolynomialShear] = {}
    visited: set[tuple[sp.Expr, int]] = set()

    def explore(
        remainder: sp.Expr,
        maximum_degree: int | None,
        polynomial: sp.Expr,
        terms: tuple[tuple[int, sp.Expr], ...],
    ) -> None:
        state = (sp.expand(remainder), maximum_degree or 0)
        if state in visited:
            return
        visited.add(state)

        for degree, coefficient in _candidate_shears(
            base, remainder, parameter
        ):
            if maximum_degree is not None and degree >= maximum_degree:
                continue
            next_remainder = sp.expand(
                remainder - coefficient * base**degree
            )
            if next_remainder == 0:
                next_remainder = sp.Integer(1)
            next_polynomial = sp.expand(
                polynomial + coefficient * sp.Symbol("_U") ** degree
            )
            next_terms = terms + ((degree, coefficient),)
            next_height = pole_height(base, next_remainder, parameter)

            if next_height < initial_height:
                key = (next_polynomial, next_remainder)
                found[key] = PolynomialShear(
                    target=target_index,
                    polynomial=next_polynomial,
                    terms=next_terms,
                    height_before=initial_height,
                    height_after=next_height,
                )

            explore(
                next_remainder,
                degree,
                next_polynomial,
                next_terms,
            )

    explore(target, None, sp.Integer(0), ())
    return tuple(
        sorted(
            found.values(),
            key=lambda move: (
                move.height_after,
                len(move.terms),
                move.terms,
                sp.default_sort_key(move.polynomial),
            ),
        )
    )


def reducing_polynomial_shears(
    first: sp.Expr,
    second: sp.Expr,
    parameter: sp.Symbol,
) -> tuple[tuple[sp.Expr, sp.Expr, PolynomialShear], ...]:
    """Exhaust strict reductions by one triangular polynomial in either direction."""
    moves: list[tuple[sp.Expr, sp.Expr, PolynomialShear]] = []

    for target_index, (base, target) in enumerate(
        ((second, first), (first, second))
    ):
        for move in _polynomial_shears_one_orientation(
            base, target, parameter, target_index
        ):
            evaluated_polynomial = sp.Integer(0)
            for degree, coefficient in move.terms:
                evaluated_polynomial += coefficient * base**degree
            reduced_target = sp.expand(target - evaluated_polynomial)
            if reduced_target == 0:
                reduced_target = sp.Integer(1)
            if target_index == 0:
                new_first, new_second = reduced_target, second
            else:
                new_first, new_second = first, reduced_target
            assert pole_height(new_first, new_second, parameter) == (
                move.height_after
            )
            moves.append((new_first, new_second, move))

    return tuple(
        sorted(
            moves,
            key=lambda item: (
                item[2].height_after,
                len(item[2].terms),
                pole_matrix(item[0], item[1], parameter),
                sp.default_sort_key(item[2].polynomial),
            ),
        )
    )


def monotone_polynomial_reduce(
    first: sp.Expr,
    second: sp.Expr,
    parameter: sp.Symbol,
) -> PolynomialReductionResult:
    """Iterate best complete polynomial shears while total height decreases."""
    current_first = sp.expand(first)
    current_second = sp.expand(second)
    steps: list[PolynomialShear] = []

    while True:
        moves = reducing_polynomial_shears(
            current_first, current_second, parameter
        )
        if not moves:
            break
        current_first, current_second, step = moves[0]
        steps.append(step)

    height = pole_height(current_first, current_second, parameter)
    assert all(
        earlier.height_after == later.height_before
        for earlier, later in zip(steps, steps[1:])
    )
    assert all(step.height_after < step.height_before for step in steps)
    return PolynomialReductionResult(
        first=current_first,
        second=current_second,
        height=height,
        steps=tuple(steps),
    )


def monotone_reduce(
    first: sp.Expr,
    second: sp.Expr,
    parameter: sp.Symbol,
) -> ReductionResult:
    """Greedily choose the least-height elementary reduction until terminal."""
    current_first = sp.expand(first)
    current_second = sp.expand(second)
    steps: list[ShearStep] = []

    while True:
        moves = reducing_shears(
            current_first, current_second, parameter
        )
        if not moves:
            break
        current_first, current_second, step = moves[0]
        steps.append(step)

    height = pole_height(current_first, current_second, parameter)
    assert all(
        earlier.height_after == later.height_before
        for earlier, later in zip(steps, steps[1:])
    )
    assert all(step.height_after < step.height_before for step in steps)
    return ReductionResult(
        first=current_first,
        second=current_second,
        height=height,
        steps=tuple(steps),
    )
