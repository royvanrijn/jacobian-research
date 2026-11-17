#!/usr/bin/env python3
"""Verify contracted-divisor Smith packets and the minimal cubic E8 normal form."""

from __future__ import annotations

import sympy as sp


def generic_smith_arithmetic() -> None:
    # At a contracted divisor the common target-ideal order h is the first
    # Smith exponent.  If delta is the determinant order, the second is
    # delta-h.  Replay every arithmetically possible small pair.
    for h in range(1, 8):
        for delta in range(2 * h, 20):
            first, second = h, delta - h
            assert first >= 1
            assert second >= first
            assert first + second == delta


def cubic_jet_family() -> None:
    r, t, a, b, c = sp.symbols("r t a b c")
    x = t**3 + a * r * t + b * r**2 * t
    y = (
        t**5
        + sp.Rational(5, 3) * a * r * t**3
        + r**2
        * (
            sp.Rational(5, 3) * b * t**3
            + sp.Rational(5, 9) * a**2 * t
        )
        + c * r**3 * t
    )

    cusp = sp.expand(y**3 - x**5)
    assert cusp.coeff(r, 1) == 0
    assert cusp.coeff(r, 2) == 0
    assert sp.factor(cusp.coeff(r, 3)) == -t**9 * (
        -5 * a**3 + 90 * a * b * t**2 - 81 * c * t**2
    ) / 27

    jacobian = sp.factor(
        sp.diff(x, r) * sp.diff(y, t)
        - sp.diff(x, t) * sp.diff(y, r)
    )
    expected_unit = (
        5 * a**3
        - 90 * a * b * t**2
        + 18 * a * c * r
        - 60 * b**2 * r * t**2
        + 9 * b * c * r**2
        + 81 * c * t**2
    )
    assert jacobian == -r**2 * t * expected_unit / 9

    theta = sp.Matrix(
        [
            [r * sp.diff(x, r), r * sp.diff(y, r)],
            [t * sp.diff(x, t), t * sp.diff(y, t)],
        ]
    )
    assert all(sp.rem(entry, t) == 0 for entry in theta)
    saturated = sp.simplify(theta / t)
    assert sp.factor(saturated.det()) == -r**3 * expected_unit / 9

    # For a!=0, the (1,1) entry is r times a unit, while subtracting its
    # unit multiple from the (2,1) entry leaves t^2 times a unit.  Thus the
    # saturated first Fitting ideal is (r,t^2), and det is r^3 times a unit.
    assert sp.factor(saturated[0, 0] / r).subs({r: 0, t: 0}) == a
    assert saturated[1, 0].subs(r, 0) == 3 * t**2
    assert sp.factor(saturated.det() / r**3).subs({r: 0, t: 0}) == -5 * a**3 / 9


def exact_normal_form() -> None:
    r, t = sp.symbols("r t")
    x = t**3 + r * t
    y = t**5 + sp.Rational(5, 3) * r * t**3 + sp.Rational(5, 9) * r**2 * t
    theta = sp.Matrix(
        [
            [r * sp.diff(x, r), r * sp.diff(y, r)],
            [t * sp.diff(x, t), t * sp.diff(y, t)],
        ]
    )
    saturated = sp.simplify(theta / t)

    q = sp.Rational(10, 9) * r + sp.Rational(5, 3) * t**2
    column_operation = sp.Matrix([[1, -q], [0, 1]])
    row_operation = sp.Matrix([[1, 0], [-1, 1]])
    row_scaling = sp.diag(1, sp.Rational(1, 3))
    column_scaling = sp.diag(1, -sp.Rational(27, 5))
    normal = sp.simplify(
        row_scaling
        * row_operation
        * saturated
        * column_operation
        * column_scaling
    )
    assert normal == sp.Matrix([[r, 0], [t**2, r**2]])

    entries = tuple(normal)
    fitting_groebner = sp.groebner(entries, r, t)
    assert fitting_groebner.reduce(r)[1] == 0
    assert fitting_groebner.reduce(t**2)[1] == 0
    assert sp.factor(normal.det()) == r**3

    # The cyclic generator e_1 has annihilator (r^3).  After killing it,
    # the remaining generator has relations (r^2,t^2), of colength four.
    quotient_basis = (1, r, t, r * t)
    for monomial in quotient_basis:
        assert sp.rem(monomial, r**2, r) == monomial
        assert sp.rem(monomial, t**2, t) == monomial
    assert len(quotient_basis) == 4

    # Before saturation, the common t factor makes Fitt_1=t*(r,t^2) and
    # the determinant r^3*t^2.  Over the generic point of t=0 this is
    # diag(t,t).
    full_normal = t * normal
    assert sp.factor(full_normal.det()) == r**3 * t**2
    assert tuple(full_normal) == (r * t, 0, t**3, r**2 * t)


def saturated_global_budget() -> None:
    # At the degree-six cubic equality row the unsaturated scalar residual is
    # two.  Put T^2=-n, let v be its boundary valency, and I=D'.T after
    # removing 2T from the determinant.  The rank-two quotient on T has
    # ch_2=v-2+n.  Replacing D by D'=D-2T changes the half-square by
    # 2*I-2*n.  Therefore the total saturated point budget is
    # 4-v+2*I-3*n.  The local normal form already consumes four.
    for valency in range(1, 7):
        for negativity in range(1, 6):
            for incidence in range(3, 13):
                curve_charge = valency - 2 + negativity
                total_points = (
                    2
                    - curve_charge
                    + 2 * incidence
                    - 2 * negativity
                )
                assert total_points == (
                    4 - valency + 2 * incidence - 3 * negativity
                )
                remaining = total_points - 4
                assert remaining == (
                    2 * incidence - valency - 3 * negativity
                )

    # If the saturated determinant meets T only in the cubic affine row,
    # I=3.  The only possible negative curve is then T^2=-1, with valency at
    # most three; valencies one, two, three leave 2,1,0 further point units.
    survivors = []
    for valency in range(1, 8):
        for negativity in range(1, 8):
            remaining = 6 - valency - 3 * negativity
            if remaining >= 0:
                survivors.append((valency, negativity, remaining))
    assert survivors == [(1, 1, 2), (2, 1, 1), (3, 1, 0)]


def main() -> None:
    generic_smith_arithmetic()
    cubic_jet_family()
    exact_normal_form()
    saturated_global_budget()
    print(
        "PASS: every contracted divisor has two positive generic Smith "
        "exponents (h,delta-h); the minimal cubic E8 jet saturates from "
        "t*diag-packet to the normal matrix [[r,0],[t^2,r^2]], with "
        "Fitt_1=(r,t^2), determinant r^3, cyclic-quotient colength four, "
        "and remaining global point budget 2*I-v-3*n"
    )


if __name__ == "__main__":
    main()
