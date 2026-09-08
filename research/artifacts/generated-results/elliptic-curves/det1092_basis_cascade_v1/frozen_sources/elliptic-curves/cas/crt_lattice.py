#!/usr/bin/env python3
"""Exact CRT and two-dimensional lattice reduction utilities.

For a residue r modulo M, rational reconstruction searches the lattice

    L(r,M) = {(a,b) in Z^2 : a - r*b == 0 (mod M)}

with basis (M,0), (r,1).  In dimension two, exact Gauss reduction is stronger
and simpler than invoking a general-purpose floating-point LLL routine.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd, log
from typing import Iterable, Sequence


Vector = tuple[int, int]


def dot(left: Vector, right: Vector) -> int:
    return left[0] * right[0] + left[1] * right[1]


def norm_squared(vector: Vector) -> int:
    return dot(vector, vector)


def nearest_integer(numerator: int, denominator: int) -> int:
    """Round an exact rational to the nearest integer (ties toward zero).

    Choosing zero on half-integer ties gives Gauss reduction a deterministic
    boundary convention and avoids a two-cycle between equally short vectors.
    """

    if denominator <= 0:
        raise ValueError("denominator must be positive")
    sign = -1 if numerator < 0 else 1
    quotient, remainder = divmod(abs(numerator), denominator)
    if 2 * remainder > denominator:
        quotient += 1
    return sign * quotient


def gauss_reduce(first: Vector, second: Vector) -> tuple[Vector, Vector]:
    """Return a Gauss-reduced basis of a rank-two integral lattice."""

    if first[0] * second[1] - first[1] * second[0] == 0:
        raise ValueError("the lattice basis is singular")
    b1, b2 = first, second
    while True:
        if norm_squared(b2) < norm_squared(b1):
            b1, b2 = b2, b1
        coefficient = nearest_integer(dot(b1, b2), norm_squared(b1))
        if coefficient == 0:
            return b1, b2
        b2 = (b2[0] - coefficient * b1[0], b2[1] - coefficient * b1[1])


def crt_pair(
    residue_left: int, modulus_left: int, residue_right: int, modulus_right: int
) -> tuple[int, int]:
    """Combine two congruences with coprime positive moduli."""

    if modulus_left <= 0 or modulus_right <= 0:
        raise ValueError("CRT moduli must be positive")
    if gcd(modulus_left, modulus_right) != 1:
        raise ValueError("this search uses CRT only for coprime moduli")
    residue_left %= modulus_left
    residue_right %= modulus_right
    step = (
        (residue_right - residue_left)
        * pow(modulus_left, -1, modulus_right)
        % modulus_right
    )
    modulus = modulus_left * modulus_right
    return (residue_left + modulus_left * step) % modulus, modulus


@dataclass(frozen=True)
class RationalRepresentative:
    numerator: int
    denominator: int
    height: int


def short_rational_representatives(
    residue: int,
    modulus: int,
    *,
    coefficient_radius: int = 8,
    limit: int = 8,
) -> tuple[RationalRepresentative, ...]:
    """Enumerate short primitive rational representatives of r modulo M.

    Only denominators invertible modulo M are retained.  Dividing a lattice
    vector by its gcd then preserves the congruence because that gcd is also a
    unit modulo M.
    """

    if modulus <= 1:
        raise ValueError("the rational-reconstruction modulus must exceed one")
    residue %= modulus
    reduced = gauss_reduce((modulus, 0), (residue, 1))
    answers: dict[tuple[int, int], RationalRepresentative] = {}
    for left in range(-coefficient_radius, coefficient_radius + 1):
        for right in range(-coefficient_radius, coefficient_radius + 1):
            if left == 0 and right == 0:
                continue
            numerator = left * reduced[0][0] + right * reduced[1][0]
            denominator = left * reduced[0][1] + right * reduced[1][1]
            if denominator == 0 or gcd(denominator, modulus) != 1:
                continue
            common = gcd(abs(numerator), abs(denominator))
            numerator //= common
            denominator //= common
            if denominator < 0:
                numerator = -numerator
                denominator = -denominator
            if gcd(denominator, modulus) != 1:
                continue
            if (numerator - residue * denominator) % modulus != 0:
                continue
            height = max(abs(numerator), denominator)
            answers[(numerator, denominator)] = RationalRepresentative(
                numerator, denominator, height
            )
    return tuple(
        sorted(
            answers.values(),
            key=lambda item: (item.height, abs(item.numerator), item.denominator),
        )[:limit]
    )


@dataclass(frozen=True)
class ConstraintChoice:
    prime: int
    modulus: int
    residue: int
    kind: str
    label: str
    local_score: float = 0.0


@dataclass(frozen=True)
class CRTState:
    residue: int
    modulus: int
    choices: tuple[ConstraintChoice, ...]
    local_score: float
    representative: RationalRepresentative | None
    objective: float


def beam_combine(
    groups: Sequence[Sequence[ConstraintChoice]],
    *,
    beam_width: int,
    height_weight: float,
    coefficient_radius: int = 8,
) -> tuple[CRTState, ...]:
    """Heuristically combine one choice from each prime group.

    This routine is deliberately non-exhaustive whenever ``beam_width`` is
    smaller than the expanded state population. The true minimum height is
    nondecreasing along a fixed branch, because its feasible sets shrink.
    Different branches can change their relative ordering, so a state
    discarded for its current height can still lead to the best final
    representative. The bounded representative enumeration supplies no
    certified minimum. Callers may use the result to rank candidates,
    but must not use beam survival as a mathematical exclusion or completeness
    certificate.  ``test_beam_width_one_counterexample`` records an exact
    width-one false negative.
    """

    if beam_width < 1:
        raise ValueError("beam width must be positive")
    states = (CRTState(0, 1, (), 0.0, None, 0.0),)
    used_primes: set[int] = set()
    for group in groups:
        if not group:
            raise ValueError("constraint groups must be nonempty")
        primes = {choice.prime for choice in group}
        if len(primes) != 1:
            raise ValueError("each constraint group must belong to one prime")
        prime = next(iter(primes))
        if prime in used_primes:
            raise ValueError("a prime may occur in only one independent CRT group")
        used_primes.add(prime)

        expanded: dict[tuple[int, int], CRTState] = {}
        for state in states:
            for choice in group:
                residue, modulus = crt_pair(
                    state.residue, state.modulus, choice.residue, choice.modulus
                )
                representatives = short_rational_representatives(
                    residue,
                    modulus,
                    coefficient_radius=coefficient_radius,
                    limit=1,
                )
                if not representatives:
                    continue
                representative = representatives[0]
                local_score = state.local_score + choice.local_score
                objective = local_score - height_weight * log(max(2, representative.height))
                candidate = CRTState(
                    residue,
                    modulus,
                    state.choices + (choice,),
                    local_score,
                    representative,
                    objective,
                )
                key = (residue, modulus)
                previous = expanded.get(key)
                if previous is None or candidate.objective > previous.objective:
                    expanded[key] = candidate
        states = tuple(
            sorted(
                expanded.values(),
                key=lambda state: (
                    -state.objective,
                    state.representative.height if state.representative else 0,
                    state.residue,
                ),
            )[:beam_width]
        )
        if not states:
            raise RuntimeError("beam pruning removed every CRT state")
    return states


def lattice_determinant(vectors: Iterable[Vector]) -> int:
    first, second = tuple(vectors)
    return first[0] * second[1] - first[1] * second[0]
