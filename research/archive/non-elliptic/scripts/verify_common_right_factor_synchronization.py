#!/usr/bin/env python3
"""Audit the common-right-factor synchronization top-jet theorem.

The proof itself is uniform.  This checker replays its triangular
reconstruction in the degree-30 and degree-42 frontier patterns and checks
the associated cut census.
"""

from __future__ import annotations

import itertools
import math

import sympy as sp


X = sp.symbols("x")


def proper_cuts(degree: int) -> tuple[int, ...]:
    return tuple(
        cut
        for cut in range(2, degree)
        if degree % cut == 0
    )


def decorated_incomparable_pairs(
    degree: int,
) -> tuple[tuple[int, int, int], ...]:
    answer = []
    for left, right in itertools.combinations(proper_cuts(degree), 2):
        if right % left == 0 or left % right == 0:
            continue
        common_right_degree = math.gcd(degree // left, degree // right)
        if common_right_degree > 1:
            answer.append((left, right, common_right_degree))
    return tuple(answer)


def reconstruct_right_factor(composite: sp.Expr, m: int, r: int) -> sp.Expr:
    """Recover a monic original degree-r right factor from the top jet."""

    degree = m * r
    source = sp.Poly(composite, X)
    unknowns = sp.symbols(f"q1:{r}")
    candidate = X**r + sum(
        unknowns[power - 1] * X**power
        for power in range(1, r)
    )
    leading_power = sp.Poly(sp.expand(candidate**m), X)
    solved: dict[sp.Symbol, sp.Expr] = {}
    for drop in range(1, r):
        variable = unknowns[r - drop - 1]
        equation = sp.expand(
            leading_power.nth(degree - drop).subs(
                solved, simultaneous=True
            )
            - source.nth(degree - drop)
        )
        linear = sp.Poly(equation, variable)
        assert linear.degree() == 1
        assert linear.nth(1) == m
        solved[variable] = sp.cancel(-linear.nth(0) / m)
    return sp.expand(candidate.subs(solved, simultaneous=True))


def audit_top_jet(m: int, r: int) -> None:
    coefficients = sp.symbols(f"rho_{m}_{r}_1:{r}")
    right_factor = X**r + sum(
        coefficients[power - 1] * X**power
        for power in range(1, r)
    )
    # Lower outer terms have degree at most (m-1)r and cannot affect the
    # top r-1 coefficients.  A symbolic leading power is therefore enough.
    composite_top = sp.expand(right_factor**m)
    reconstructed = reconstruct_right_factor(composite_top, m, r)
    assert sp.expand(reconstructed - right_factor) == 0

    linear = sp.symbols(f"lambda_{m}_{r}")
    reconstructed_after_hessian_shift = reconstruct_right_factor(
        composite_top + linear * X,
        m,
        r,
    )
    assert (
        sp.expand(reconstructed_after_hessian_shift - right_factor) == 0
    )


def audit_bad_characteristic_boundary() -> None:
    epsilon, z = sp.symbols("epsilon z")
    outer = z**2 + z
    first = X**2 + epsilon * X
    second = X**2
    difference = sp.expand(
        outer.subs(z, first) - outer.subs(z, second) - epsilon * X
    )
    quotient_basis = sp.groebner(
        [epsilon**2],
        X,
        epsilon,
        modulus=2,
    )
    assert quotient_basis.reduce(difference)[1] == 0
    assert quotient_basis.reduce(first - second)[1] != 0


EXPECTED = {
    30: ((2, 3, 5), (2, 5, 3), (3, 5, 2)),
    42: ((2, 3, 7), (2, 7, 3), (3, 7, 2)),
}


def main() -> None:
    audited_shapes = set()
    for degree, expected in EXPECTED.items():
        actual = decorated_incomparable_pairs(degree)
        assert actual == expected
        for left, right, common_right_degree in actual:
            # The top-jet theorem only sees the total outer degree.
            shape = (degree // common_right_degree, common_right_degree)
            if shape not in audited_shapes:
                audit_top_jet(*shape)
                audited_shapes.add(shape)
            print(
                f"degree {degree}: cuts {{{left},{right}}}, "
                f"common right degree {common_right_degree}"
            )
    print(
        "PASS: the degree-30 and degree-42 decorated-pair census is exact"
    )
    print(
        "PASS: every occurring terminal factor is unchanged by a linear "
        "Hessian lift"
    )
    audit_bad_characteristic_boundary()
    print(
        "PASS: the characteristic-two dual-number boundary counterexample "
        "is exact"
    )


if __name__ == "__main__":
    main()
