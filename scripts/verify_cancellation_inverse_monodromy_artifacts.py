#!/usr/bin/env python3
"""Exact permutation replay of stored cancellation monodromy searches."""

from __future__ import annotations

import json
import math
from pathlib import Path

from sympy.combinatorics import Permutation, PermutationGroup

from search_cancellation_inverse_monodromy import exact_model


ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = (
    ROOT / "artifacts/cancellation_inverse_monodromy_small.json",
    ROOT / "artifacts/cancellation_inverse_monodromy_6_4.json",
    ROOT / "artifacts/cancellation_inverse_monodromy_6_5.json",
)


def main() -> None:
    for m in range(1, 7):
        for r in range(1, 6):
            model = exact_model(m, r)
            assert model.degree == r * (m + 1) + 1
            assert model.critical_value_polynomial.degree() == m + 1

    seen: set[tuple[int, int]] = set()
    for path in ARTIFACTS:
        records = json.loads(path.read_text())
        for record in records:
            m = int(record["m"])
            r = int(record["r"])
            degree = int(record["degree"])
            assert (m, r) not in seen
            seen.add((m, r))
            assert degree == r * (m + 1) + 1
            assert len(record["generators_zero_based"]) == m + 1
            assert all(
                lengths == [r + 1] for lengths in record["cycle_lengths"]
            )

            generators = [
                Permutation(permutation)
                for permutation in record["generators_zero_based"]
            ]
            group = PermutationGroup(*generators)
            assert group.is_transitive()
            assert group.is_primitive()
            expected_order = math.factorial(degree)
            expected_name = f"S_{degree}"
            if r % 2 == 0:
                expected_order //= 2
                expected_name = f"A_{degree}"
            assert group.order() == expected_order
            assert record["classification"] == expected_name

    assert len(seen) == 11
    assert (6, 4) in seen and (6, 5) in seen
    print(
        "PASS: exact critical-value separation on the 6x5 grid and exact "
        "replay of 11 stored primitive A_N/S_N branch-cycle tuples"
    )


if __name__ == "__main__":
    main()
