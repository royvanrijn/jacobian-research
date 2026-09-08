#!/usr/bin/env python3
"""Dependency-free policy helpers for the R17 Kummer quotient collector.

This module contains no number-field arithmetic.  It makes the search policy
and the two quotient projections independently testable under ordinary
CPython; the Sage collector is responsible for every ideal identity.
"""

from __future__ import annotations

from dataclasses import dataclass
import random
from typing import Iterable, Mapping, Sequence


GENERIC_POINT_COUNT = 17
STRATEGIES = ("single", "pair", "sparse")


class BinaryRows:
    """Incremental row space over GF(2), using packed integer rows."""

    def __init__(self) -> None:
        self.pivots: dict[int, int] = {}
        self.rows: list[int] = []

    def reduce(self, row: int) -> int:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in self.pivots:
                break
            row ^= self.pivots[pivot]
        return row

    def add(self, row: int) -> bool:
        original = row
        row = self.reduce(row)
        if not row:
            return False
        self.pivots[row.bit_length() - 1] = row
        self.rows.append(original)
        return True

    def free_column(self, width: int, start: int = 0) -> int | None:
        if width < 0 or not (0 <= start <= max(width, 0)):
            raise ValueError("invalid free-column search interval")
        for offset in range(width):
            index = (start + offset) % width
            if self.reduce(1 << index):
                return index
        return None

    def free_column_combination(
        self, width: int, start: int = 0, count: int = 1
    ) -> tuple[int, ...]:
        """Choose a nonzero quotient class supported on several free columns.

        Columns are scanned cyclically.  A column is eligible only when its
        unit vector is unresolved, and it is retained only when the cumulative
        mask remains nonzero modulo the current row space.
        """

        if count <= 0:
            raise ValueError("the requested target-column count must be positive")
        if width < 0 or not (0 <= start <= max(width, 0)):
            raise ValueError("invalid free-column search interval")
        selected = []
        mask = 0
        for offset in range(width):
            index = (start + offset) % width
            unit = 1 << index
            if not self.reduce(unit):
                continue
            candidate = mask ^ unit
            if not self.reduce(candidate):
                continue
            selected.append(index)
            mask = candidate
            if len(selected) == count:
                break
        return tuple(selected)

    @property
    def rank(self) -> int:
        return len(self.pivots)


@dataclass(frozen=True)
class CompanionTerm:
    """One known half-ideal and its short signed exponent."""

    index: int
    label: str
    exponent: int
    role: str

    @property
    def parity(self) -> int:
        return abs(self.exponent) & 1


def parse_strategies(value: str) -> tuple[str, ...]:
    strategies = tuple(item.strip() for item in value.split(",") if item.strip())
    if not strategies:
        raise ValueError("at least one companion strategy is required")
    unknown = sorted(set(strategies) - set(STRATEGIES))
    if unknown:
        raise ValueError(f"unknown companion strategies: {', '.join(unknown)}")
    return strategies


def _weighted_sample_without_replacement(
    rng: random.Random,
    population: Sequence[int],
    count: int,
    weights: Mapping[int, int],
) -> list[int]:
    """Sample distinct indices with reproducible positive integral weights."""

    remaining = list(population)
    chosen: list[int] = []
    for _ in range(min(count, len(remaining))):
        total = sum(weights[index] for index in remaining)
        ticket = rng.randrange(total)
        subtotal = 0
        for offset, index in enumerate(remaining):
            subtotal += weights[index]
            if ticket < subtotal:
                chosen.append(index)
                remaining.pop(offset)
                break
    return chosen


def select_companion_terms(
    *,
    rng: random.Random,
    attempt: int,
    labels: Sequence[str],
    generic_point_count: int = GENERIC_POINT_COUNT,
    strategies: Sequence[str] = STRATEGIES,
    sparse_min: int = 3,
    sparse_max: int = 6,
    exponent_radius: int = 1,
    exceptional_weight: int = 4,
    signed_exponents: bool = False,
) -> tuple[str, tuple[CompanionTerm, ...]]:
    """Choose an ``I_i``, ``I_i I_j``, or short signed ideal product.

    Single and pair trials use exponent +1 exactly.  Sparse trials use
    positive exponents up to the declared radius, or nonzero signed exponents
    when requested.  Exceptional point ideals
    receive a larger sampling weight, but every known ideal remains reachable.
    """

    if attempt <= 0:
        raise ValueError("attempt numbers start at one")
    if not labels:
        return strategies[(attempt - 1) % len(strategies)], ()
    if not (0 <= generic_point_count <= len(labels)):
        raise ValueError("invalid generic point count")
    if sparse_min < 0 or sparse_min > sparse_max:
        raise ValueError("invalid sparse companion interval")
    if exponent_radius <= 0 or exceptional_weight <= 0:
        raise ValueError("exponent radius and exceptional weight must be positive")
    if not strategies or any(strategy not in STRATEGIES for strategy in strategies):
        raise ValueError("invalid companion strategy sequence")

    strategy = strategies[(attempt - 1) % len(strategies)]
    if strategy == "single":
        count = 1
    elif strategy == "pair":
        count = 2
    else:
        count = rng.randint(sparse_min, sparse_max)

    population = list(range(len(labels)))
    weights = {
        index: (1 if index < generic_point_count else exceptional_weight)
        for index in population
    }
    indices = _weighted_sample_without_replacement(
        rng, population, count, weights
    )
    terms = []
    for index in indices:
        exponent = 1
        if strategy == "sparse":
            choices = tuple(range(1, exponent_radius + 1))
            if signed_exponents:
                choices = tuple(range(-exponent_radius, 0)) + choices
            exponent = rng.choice(choices)
        terms.append(
            CompanionTerm(
                index=index,
                label=labels[index],
                exponent=exponent,
                role=(
                    "generic_MW17"
                    if index < generic_point_count
                    else "known_exceptional"
                ),
            )
        )
    return strategy, tuple(terms)


def projection_masks(
    *,
    base_mask: int,
    exceptional_parity_mask: int,
    factor_base_width: int,
) -> tuple[int, int]:
    """Return rows modulo generic MW17 and modulo the full known subgroup."""

    if base_mask < 0 or exceptional_parity_mask < 0 or factor_base_width < 0:
        raise ValueError("packed masks and widths must be nonnegative")
    if base_mask >> factor_base_width:
        raise ValueError("base mask exceeds the factor-base width")
    generic_mask = base_mask | (exceptional_parity_mask << factor_base_width)
    full_known_mask = base_mask
    return generic_mask, full_known_mask


def independent_rank(rows: Iterable[int]) -> int:
    space = BinaryRows()
    for row in rows:
        space.add(row)
    return space.rank
