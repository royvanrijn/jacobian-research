#!/usr/bin/env python3
"""Verify the literature-derived geometric-degree and generic-fiber sieve."""

from __future__ import annotations

from fractions import Fraction
from math import ceil


def chau_semigroup_candidates(limit: int) -> tuple[int, ...]:
    values = {
        5 * r + 3 * s
        for r in range(limit + 1)
        for s in range(limit + 1)
        if r + s >= 1 and 5 * r + 3 * s <= limit
    }
    return tuple(sorted(value for value in values if value >= 6))


def audit_degree_sieve() -> None:
    # Makar-Limanov: N < (n-m) n b0 / (n(a0+b0)-a0).
    m, n = 15, 60
    a0, b0 = 3, 5
    upper = Fraction((n - m) * n * b0, n * (a0 + b0) - a0)
    assert upper == Fraction(1500, 53)
    assert upper > 28
    assert upper < 29

    candidates = chau_semigroup_candidates(28)
    assert candidates == (6, *range(8, 29))
    assert 7 not in candidates
    assert tuple(value for value in candidates if value >= 12) == tuple(range(12, 29))


def audit_generic_fiber_identities() -> None:
    # A_j is moved-sheet degree and F_j is residue degree over a target
    # component with normalization degrees (3 k_j, 5 k_j).
    for degree in range(6, 29):
        for weighted_a in range(1, 18):
            for weighted_f in range(weighted_a + 1):
                contact_different = weighted_a - weighted_f
                for p in range(1, 8):
                    for q in range(1, 8):
                        chi_p = degree - 3 * weighted_a
                        chi_q = degree - 5 * weighted_a
                        punctures_p = p + 3 * weighted_f
                        punctures_q = q + 5 * weighted_f
                        twice_genus_p = 2 - punctures_p - chi_p
                        twice_genus_q = 2 - punctures_q - chi_q

                        assert twice_genus_p == (
                            2 - degree - p + 3 * contact_different
                        )
                        assert twice_genus_q == (
                            2 - degree - q + 5 * contact_different
                        )

                        # Chau's degree-weighted Euler identity for
                        # (deg P,deg Q)=(75,125).
                        assert 75 * (chi_q - degree) == 125 * (chi_p - degree)

                        if twice_genus_p >= 0:
                            assert contact_different >= ceil((degree + p - 2) / 3)
                        if twice_genus_q >= 0:
                            assert contact_different >= ceil((degree + q - 2) / 5)


def main() -> None:
    audit_degree_sieve()
    audit_generic_fiber_identities()
    print(
        "PASS: F2 geometric degree is in {6,8,...,28}; "
        "generic-fiber Euler/genus identities retain the unknown pole counts"
    )


if __name__ == "__main__":
    main()
