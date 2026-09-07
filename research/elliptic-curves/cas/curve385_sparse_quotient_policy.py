#!/usr/bin/env python3
"""Deterministic sparse quotient-word policy for the curve-385 rank-32 search."""

from __future__ import annotations

from hashlib import sha256
from itertools import combinations
import json
from typing import Any, Iterable, Sequence


POLICY_DOMAIN = "curve385-sparse-rank32-v1"
STAGE_SPECS = (
    ("natural-weight-1", "natural", (1,)),
    ("natural-weight-2", "natural", (2,)),
    ("alternate-a-weight-at-most-2", "alternate-a", (1, 2)),
    ("alternate-b-weight-at-most-2", "alternate-b", (1, 2)),
    ("natural-weight-3", "natural", (3,)),
    ("natural-weight-4", "natural", (4,)),
)


def canonical_hash(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256(encoded).hexdigest()


def gf2_rank(words: Iterable[int]) -> int:
    pivots: list[int] = []
    for word in words:
        reduced = int(word)
        for pivot in pivots:
            reduced = min(reduced, reduced ^ pivot)
        if not reduced:
            continue
        pivots.append(reduced)
        pivots.sort(reverse=True)
    return len(pivots)


def quotient_basis_words(bit_count: int, basis_id: str) -> tuple[int, ...]:
    if bit_count < 1:
        raise ValueError("quotient dimension must be positive")
    if basis_id == "natural":
        return tuple(1 << index for index in range(bit_count))
    if basis_id not in {"alternate-a", "alternate-b"}:
        raise ValueError(f"unknown quotient basis {basis_id!r}")
    mask = (1 << bit_count) - 1
    answer: list[int] = []
    counter = 0
    while len(answer) < bit_count:
        digest = sha256(
            f"{POLICY_DOMAIN}:{bit_count}:{basis_id}:{counter}".encode()
        ).digest()
        counter += 1
        candidate = int.from_bytes(digest, "big") & mask
        if candidate and gf2_rank(answer + [candidate]) > len(answer):
            answer.append(candidate)
    if gf2_rank(answer) != bit_count:
        raise ArithmeticError("deterministic quotient basis is singular")
    return tuple(answer)


def physical_words(
    basis_words: Sequence[int], weights: Iterable[int]
) -> tuple[int, ...]:
    weights = tuple(sorted(set(map(int, weights))))
    if not weights or weights[0] < 1 or weights[-1] > len(basis_words):
        raise ValueError("invalid quotient-weight shell")
    answer = set()
    for weight in weights:
        for indices in combinations(range(len(basis_words)), weight):
            word = 0
            for index in indices:
                word ^= int(basis_words[index])
            if not word:
                raise ArithmeticError("a sparse word vanished in a quotient basis")
            answer.add(word)
    return tuple(sorted(answer))


def stage_plan(bit_count: int, old_class_count: int = 43) -> list[dict[str, Any]]:
    seen: set[int] = set()
    answer = []
    for index, (stage_id, basis_id, weights) in enumerate(STAGE_SPECS, start=1):
        basis = quotient_basis_words(bit_count, basis_id)
        selected = physical_words(basis, weights)
        new_words = tuple(word for word in selected if word not in seen)
        seen.update(selected)
        answer.append(
            {
                "index": index,
                "id": stage_id,
                "basis_id": basis_id,
                "basis_words_in_natural_coordinates": list(basis),
                "basis_words_binary": [f"{word:0{bit_count}b}" for word in basis],
                "basis_words_sha256": canonical_hash(list(basis)),
                "basis_weight_shell": list(weights),
                "selected_physical_word_count_before_deduplication": len(selected),
                "new_physical_words": list(new_words),
                "new_physical_words_sha256": canonical_hash(list(new_words)),
                "new_physical_word_count": len(new_words),
                "new_chart_count": old_class_count * len(new_words),
                "cumulative_physical_word_count": len(seen),
                "cumulative_chart_count": old_class_count * len(seen),
            }
        )
    return answer


def validate_stage_plan(
    plan: Sequence[dict[str, Any]], bit_count: int, old_class_count: int = 43
) -> None:
    expected = stage_plan(bit_count, old_class_count)
    if list(plan) != expected:
        raise ArithmeticError("sparse quotient stage plan differs from policy")
    seen: set[int] = set()
    for expected_index, row in enumerate(plan, start=1):
        if row["index"] != expected_index or not row["new_physical_words"]:
            raise ArithmeticError("sparse quotient stage order is invalid")
        new_words = set(map(int, row["new_physical_words"]))
        if seen.intersection(new_words):
            raise ArithmeticError("a physical quotient word occurs in two stages")
        seen.update(new_words)
        if row["new_chart_count"] != old_class_count * len(new_words):
            raise ArithmeticError("a sparse stage has the wrong chart count")
