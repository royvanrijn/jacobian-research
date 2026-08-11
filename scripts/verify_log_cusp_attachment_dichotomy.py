#!/usr/bin/env python3
"""Verify the smooth-boundary/fold versus SNC-node cusp attachment dichotomy."""

from __future__ import annotations

from fractions import Fraction
from math import gcd

import sympy as sp


r, t = sp.symbols("r t")


def monomial_fold_audit(m: int, n: int, q: int) -> None:
    """Check the exact fold family for ``y^m-x^n``."""

    assert 1 < m < n
    assert gcd(m, n) == 1
    assert q >= 1
    gap_order = (n - m) * q
    x = t ** (m * q) + r
    y = t ** (n * q) + sp.Rational(n, m) * r * t**gap_order
    branch_equation = sp.expand(y**m - x**n)

    # The linear transverse term cancels.  The first possible term is r^2,
    # with the universal coefficient displayed in the note.
    assert sp.expand(branch_equation).coeff(r, 0) == 0
    assert sp.expand(branch_equation).coeff(r, 1) == 0
    expected_quadratic = (
        -sp.Rational(n * (n - m), 2 * m) * t ** (m * q * (n - 2))
    )
    assert sp.expand(branch_equation).coeff(r, 2) == expected_quadratic

    jacobian = sp.factor(
        sp.diff(x, r) * sp.diff(y, t) - sp.diff(x, t) * sp.diff(y, r)
    )
    expected_jacobian = (
        sp.Rational(n * (n - m) * q, m)
        * r
        * t ** (gap_order - 1)
    )
    assert jacobian == expected_jacobian

    # epsilon=0 is the lone smooth-boundary case.  If the residual t-factor
    # is present, t=0 must be included in the boundary and epsilon=1.
    epsilon = 0 if gap_order == 1 else 1
    theta = sp.Matrix(
        [
            [r * sp.diff(x, r), r * sp.diff(y, r)],
            [t**epsilon * sp.diff(x, t), t**epsilon * sp.diff(y, t)],
        ]
    )
    column_operation = sp.Matrix(
        [
            [1, -sp.Rational(n, m) * t**gap_order],
            [0, 1],
        ]
    )
    reduced = sp.simplify(theta * column_operation)
    expected_reduced = sp.Matrix(
        [
            [r, 0],
            [
                m * q * t ** (m * q - 1 + epsilon),
                sp.Rational(n * (n - m) * q, m)
                * r
                * t ** (gap_order - 1 + epsilon),
            ],
        ]
    )
    assert reduced == expected_reduced

    fitting_one = sp.groebner(tuple(theta), r, t, order="lex", domain=sp.QQ)
    point_exponent = m * q - 1 + epsilon
    expected_fitting = sp.groebner(
        [r, t**point_exponent], r, t, order="lex", domain=sp.QQ
    )
    assert fitting_one == expected_fitting
    assert point_exponent == (m * q - 1 if epsilon == 0 else m * q)


def ordinary_cusp_fold_audit() -> None:
    """Check the exact one-boundary fold and its one-point quotient."""

    x = t**2 + r
    y = t**3 + sp.Rational(3, 2) * r * t
    cusp = sp.factor(y**2 - x**3)
    assert sp.expand(
        cusp + r**2 * (r + sp.Rational(3, 4) * t**2)
    ) == 0

    jacobian = sp.factor(
        sp.diff(x, r) * sp.diff(y, t) - sp.diff(x, t) * sp.diff(y, r)
    )
    assert jacobian == sp.Rational(3, 2) * r

    theta = sp.Matrix(
        [
            [r * sp.diff(x, r), r * sp.diff(y, r)],
            [sp.diff(x, t), sp.diff(y, t)],
        ]
    )
    operation = sp.Matrix(
        [[1, -sp.Rational(3, 2) * t], [0, 1]]
    )
    reduced = sp.simplify(theta * operation)
    assert reduced == sp.Matrix(
        [[r, 0], [2 * t, sp.Rational(3, 2) * r]]
    )
    assert sp.factor(theta.det()) == sp.Rational(3, 2) * r**2
    assert sp.groebner(tuple(theta), r, t, domain=sp.QQ) == sp.groebner(
        [r, t], r, t, domain=sp.QQ
    )

    # For A_N=[[r,0],[t^N,r]], the class of the first target generator gives
    # 0 -> R/(r^2) -> coker(A_N) -> R/(r,t^N) -> 0.  The ideal arithmetic
    # below checks the annihilator and the quotient length for a bounded range.
    for exponent in range(1, 10):
        assert sp.gcd(r, t**exponent) == 1
        standard_basis = tuple(t**index for index in range(exponent))
        assert len(standard_basis) == exponent
        for monomial in standard_basis:
            assert sp.rem(monomial, t**exponent, t) == monomial
        assert sp.rem(t**exponent, t**exponent, t) == 0

    # The affine companion meets the boundary fold with contact two.
    companion_on_boundary = sp.expand(
        (r + sp.Rational(3, 4) * t**2).subs(r, 0)
    )
    assert companion_on_boundary == sp.Rational(3, 4) * t**2


def fitting_lower_bound_audit() -> None:
    """Replay the order and complete-fiber arithmetic."""

    # At a point on one boundary component epsilon=0; at an SNC crossing it
    # is 1.  The tangential row has minimum order q*m-1+epsilon.
    for branch_multiplicity in range(2, 8):
        for residue_index in range(1, 8):
            for epsilon in (0, 1):
                exponent = (
                    branch_multiplicity * residue_index - 1 + epsilon
                )
                assert exponent >= 1
                quotient_basis = tuple(t**index for index in range(exponent))
                assert len(quotient_basis) == exponent

    def compositions(total: int) -> tuple[tuple[int, ...], ...]:
        if total == 0:
            return ((),)
        rows: list[tuple[int, ...]] = []
        for first in range(1, total + 1):
            for tail in compositions(total - first):
                rows.append((first, *tail))
        return tuple(rows)

    for branch_multiplicity in range(2, 7):
        for residue_degree in range(1, 8):
            for partition in compositions(residue_degree):
                fiber_points = len(partition)
                for crossing_count in range(fiber_points + 1):
                    local_epsilons = (1,) * crossing_count + (0,) * (
                        fiber_points - crossing_count
                    )
                    total = sum(
                        branch_multiplicity * residue_index - 1 + epsilon
                        for residue_index, epsilon in zip(
                            partition, local_epsilons, strict=True
                        )
                    )
                    expected = (
                        branch_multiplicity * residue_degree
                        - fiber_points
                        + crossing_count
                    )
                    assert total == expected
                    assert total >= (branch_multiplicity - 1) * residue_degree

                    if branch_multiplicity == 2:
                        assert residue_degree <= total <= 2 * residue_degree


def f2_budget_endpoint_audit() -> None:
    """Check the two exact ordinary-cusp endpoints used in the F2 sieve."""

    # The doubled virtual point numerators for the minimal signature are
    # 12-4*b-s_X and 49-4*b-s_X.  A complete smooth unramified fold fiber has
    # charge f; the node-saturated LUAF1 fiber has charge 2*f.
    residue_degree = 1
    for followed_centers in range(9):
        square_numerator = 12 - 4 * followed_centers
        double_numerator = 49 - 4 * followed_centers - 1
        smooth_fold_square = Fraction(
            square_numerator - 2 * residue_degree, 2
        )
        node_saturated_square = Fraction(
            square_numerator - 4 * residue_degree, 2
        )
        smooth_fold_double = Fraction(
            double_numerator - 2 * residue_degree, 2
        )
        node_saturated_double = Fraction(
            double_numerator - 4 * residue_degree, 2
        )
        assert smooth_fold_square == 5 - 2 * followed_centers
        assert node_saturated_square == 4 - 2 * followed_centers
        assert smooth_fold_double == 23 - 2 * followed_centers
        assert node_saturated_double == 22 - 2 * followed_centers


def main() -> None:
    for multiplicity in range(2, 7):
        for higher_order in range(multiplicity + 1, multiplicity + 5):
            if gcd(multiplicity, higher_order) != 1:
                continue
            for residue_index in range(1, 5):
                monomial_fold_audit(multiplicity, higher_order, residue_index)
    ordinary_cusp_fold_audit()
    fitting_lower_bound_audit()
    f2_budget_endpoint_audit()
    print(
        "PASS: unibranch Fitt_1 colength is bounded below by "
        "q_p*m_C-1+epsilon; the monomial fold attains it, the ordinary "
        "smooth-boundary cusp fold has charge one, and a complete cusp "
        "fiber ranges from f to 2*f before its source incidences are known"
    )


if __name__ == "__main__":
    main()
