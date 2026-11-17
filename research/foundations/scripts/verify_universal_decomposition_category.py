#!/usr/bin/env python3
"""Degree-only regressions for Programme 4.

These finite checks locate the first degrees at which arithmetic permits
competing atomic words.  They do not assert that one Keller map realizes two
such words; the geometric and polynomial-sandwich gates are separate.
"""

from __future__ import annotations

from functools import lru_cache
from math import prod


@lru_cache(maxsize=None)
def ordered_atomic_words(degree: int) -> tuple[tuple[int, ...], ...]:
    """Return ordered nontrivial words with factors at least three."""

    words: set[tuple[int, ...]] = set()
    for first in range(3, degree + 1):
        if degree % first:
            continue
        quotient = degree // first
        if quotient >= 3:
            words.add((first, quotient))
            for suffix in ordered_atomic_words(quotient):
                words.add((first,) + suffix)
    return tuple(sorted(words))


def unordered_patterns(degree: int) -> set[tuple[int, ...]]:
    return {tuple(sorted(word)) for word in ordered_atomic_words(degree)}


def word_lengths(degree: int) -> set[int]:
    return {len(word) for word in ordered_atomic_words(degree)}


def first_degree(predicate, limit: int = 200) -> int:
    for degree in range(3, limit + 1):
        if predicate(degree):
            return degree
    raise AssertionError(f"no witness through {limit}")


for degree in range(3, 201):
    for word in ordered_atomic_words(degree):
        assert len(word) >= 2
        assert all(factor >= 3 for factor in word)
        assert prod(word) == degree

# The first degree supporting any decomposable Keller example is 9.
assert first_degree(lambda degree: bool(ordered_atomic_words(degree))) == 9

# Reversing 3 and 4 gives the first distinct ordered two-factor words.
assert first_degree(
    lambda degree: len(
        {word for word in ordered_atomic_words(degree) if len(word) == 2}
    )
    > 1
) == 12
assert (3, 4) in ordered_atomic_words(12)
assert (4, 3) in ordered_atomic_words(12)

# Degree 24 is the first place where the unordered two-factor multisets can
# differ: {3,8} and {4,6}.
assert first_degree(lambda degree: len(unordered_patterns(degree)) > 1) == 24
assert {(3, 8), (4, 6)} <= unordered_patterns(24)

# Degree 27 is the first place where nontrivial atomic words can have
# different lengths: (3,9) versus (3,3,3).
assert first_degree(lambda degree: len(word_lengths(degree)) > 1) == 27
assert {2, 3} <= word_lengths(27)

print("PASS: degree 9 is the first arithmetic decomposable degree")
print("PASS: degree 12 is the first ordered adjacent-swap laboratory")
print("PASS: degree 24 is the first competing unordered factor-degree laboratory")
print("PASS: degree 27 is the first competing factor-length laboratory")
print("NOTE: these are degree-word possibilities, not common-map factorizations")
