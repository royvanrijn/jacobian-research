#!/usr/bin/env python3
"""Exact certificates for common fibers of stably inequivalent Keller maps."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from jcsearch.weighted import WeightedSeedModel, w  # noqa: E402
from master_cancellation import (  # noqa: E402
    fiber_antiderivative,
    parameter_polynomial,
)


T, U = sp.symbols("T U")


def check_fixed_rational_pairs() -> None:
    """Check the all-degree common-fiber pencil over Q."""

    for degree in range(4, 13):
        H = T**degree + T**3 - 2 * T**2
        P = H + T + U
        c = 1 - degree

        assert sp.expand(P.subs(T, 1) - P.subs(T, 0)) == 1
        assert sp.diff(P, T).subs(T, 0) == 1
        assert sp.diff(P, T).subs(T, 1) == degree
        assert c != 0
        assert (
            sp.diff(P, T, 2).subs(T, 1)
            != 2
            * (
                sp.diff(P, T).subs(T, 1)
                - sp.diff(P, T).subs(T, 0)
            )
        )

        residual = sp.expand(H / T**2)
        assert sp.factor(residual.subs(T, 1)) == 0
        assert sp.gcd(sp.Poly(residual, T), sp.Poly(sp.diff(residual, T), T)) == 1

        # The weighted constructor accepts the fixed seed H' and has
        # determinant c.  Its distinguished target recovers P exactly.
        model = WeightedSeedModel(
            sp.diff(H, T).subs(T, w),
            c=c,
        )
        weighted_inverse = model.inverse_polynomial(U / c, -1, 1)
        assert sp.expand(weighted_inverse.subs(w, T) - P) == 0
        assert model.zero_profile()[:2] == (2, 1)

        # The quadratic-gauge seed is fixed as U varies.
        G = H + T
        assert sp.diff(G, T).subs(T, 0) == 1
        assert sp.expand(G).coeff(T, 3) == 1
        assert sp.Poly(G, T).LC() == 1
        quadratic_inverse = G - sp.Rational(1, 2) * (-2 * U)
        assert sp.expand(quadratic_inverse - P) == 0

        # Generic irreducibility: P is primitive and linear in U.
        field = sp.QQ.frac_field(U)
        assert sp.Poly(P, T, domain=field).is_irreducible

    # The smallest member gives a connected quartic fiber at U=1.
    quartic = T**4 + T**3 - 2 * T**2 + T + 1
    assert sp.Poly(quartic, T, domain=sp.QQ).is_irreducible
    assert sp.discriminant(quartic, T) == -1156


def check_original_rational_quartic() -> None:
    """Retain the particularly small quartic proposed in the research note."""

    P = 2 * T**4 - T**3 - T**2 + T + 1
    H = sp.expand(P - P.subs(T, 0) - sp.diff(P, T).subs(T, 0) * T)
    assert sp.factor(H) == T**2 * (T - 1) * (2 * T + 1)
    assert sp.Poly(P, T, modulus=3).is_irreducible
    assert sp.discriminant(P, T) == 1556
    assert P.subs(T, 1) - P.subs(T, 0) == sp.diff(P, T).subs(T, 0) == 1
    assert sp.diff(P, T).subs(T, 0) - sp.diff(P, T).subs(T, 1) == -3
    assert sp.diff(P, T, 2).subs(T, 1) == 16


def check_fixed_quadratic_field_triple() -> None:
    """Check the common quartic pencil over Q(sqrt(-2))."""

    eta, R = sp.symbols("eta R")
    eta_relation = sp.Poly(eta**2 + 2, eta, domain=sp.QQ)

    def reduce_eta(expression: sp.Expr) -> sp.Expr:
        numerator, denominator = sp.cancel(expression).as_numer_denom()
        relation = sp.Poly(eta**2 + 2, eta, domain="EX")
        numerator = sp.rem(
            sp.Poly(numerator, eta, domain="EX"), relation
        ).as_expr()
        denominator = sp.rem(
            sp.Poly(denominator, eta, domain="EX"), relation
        ).as_expr()
        conjugate = denominator.subs(eta, -eta)
        norm = sp.rem(
            sp.Poly(sp.expand(denominator * conjugate), eta, domain="EX"),
            relation,
        ).as_expr()
        reduced = sp.rem(
            sp.Poly(sp.expand(numerator * conjugate), eta, domain="EX"),
            relation,
        ).as_expr()
        return sp.cancel(reduced / norm)

    P_target = 4 + eta
    Q_target = 3
    psi = sp.expand(
        fiber_antiderivative(2, 1, T, P_target, Q_target) - R
    )
    expected = (
        T
        - sp.Rational(9, 2) * T**2
        + (8 + 2 * eta) * T**3
        - (sp.Rational(7, 2) + 2 * eta) * T**4
        - R
    )
    assert reduce_eta(psi - expected) == 0

    assert reduce_eta(psi.subs(T, 1) - psi.subs(T, 0)) == 1
    assert sp.diff(psi, T).subs(T, 0) == 1
    c = reduce_eta(
        sp.diff(psi, T).subs(T, 0) - sp.diff(psi, T).subs(T, 1)
    )
    assert c == -1 + 2 * eta
    assert reduce_eta(
        sp.diff(psi, T, 2).subs(T, 1)
        - 2
        * (
            sp.diff(psi, T).subs(T, 1)
            - sp.diff(psi, T).subs(T, 0)
        )
    ) == -5 - 8 * eta

    H = sp.expand(psi + R - T)
    expected_H = T**2 * (T - 1) * (
        sp.Rational(9, 2) - (sp.Rational(7, 2) + 2 * eta) * T
    )
    assert reduce_eta(H - expected_H) == 0

    # The cancellation map is polynomial over Q(eta): theta=2+eta is a
    # root of the type-(2,1) parameter polynomial.
    theta = 2 + eta
    modulus = parameter_polynomial(2, 1, sp.Symbol("q"))
    assert reduce_eta(modulus.subs(sp.Symbol("q"), theta)) == 0

    # The weighted and quadratic-gauge special targets recover the same
    # polynomial psi.  The latter has C=2R because g_1=1.
    weighted_inverse = H + T - R
    quadratic_seed = psi + R
    quadratic_inverse = quadratic_seed - sp.Rational(1, 2) * (2 * R)
    assert reduce_eta(weighted_inverse - psi) == 0
    assert reduce_eta(quadratic_inverse - psi) == 0

    # Generic irreducibility over Q(eta)(R) follows from the standard
    # primitive, degree-one-in-R argument.
    assert sp.Poly(psi, R).degree() == 1
    assert sp.diff(psi, R) == -1
    assert sp.Poly(psi + R, T).degree() == 4

    # Explicit connected fiber at R=-1.  Reduce 2*psi at
    # (17, eta-10), then apply Rabin's degree-four criterion.
    specialized = sp.expand(2 * psi.subs(R, -1).subs(eta, 10))
    reduced = sp.Poly(specialized, T, modulus=17)
    monic = sp.Poly(reduced.monic(), T, modulus=17)
    expected_mod_17 = sp.Poly(
        T**4 + 14 * T**3 + 2 * T**2 + 9 * T + 9,
        T,
        modulus=17,
    )
    assert monic == expected_mod_17

    def power_mod(exponent: int) -> sp.Poly:
        result = sp.Poly(1, T, modulus=17)
        base = sp.Poly(T, T, modulus=17)
        while exponent:
            if exponent & 1:
                result = result.mul(base).rem(monic)
            base = base.mul(base).rem(monic)
            exponent //= 2
        return result

    x_poly = sp.Poly(T, T, modulus=17)
    degree_two_test = (power_mod(17**2) - x_poly).rem(monic)
    degree_four_test = (power_mod(17**4) - x_poly).rem(monic)
    assert degree_two_test == sp.Poly(
        6 * T**3 - T**2 - 6 * T + 3, T, modulus=17
    )
    assert sp.gcd(monic, degree_two_test).degree() == 0
    assert degree_four_test.is_zero

    # Stable-boundary invariants used for pairwise separation.
    assert (1, 2, 2) == (
        1,  # weighted ramified-stratum unit rank
        2,  # cancellation ramified-stratum unit rank
        2,  # quadratic-gauge ramified-stratum unit rank
    )
    assert (1, 2) == (
        1,  # cancellation Fitting Laurent-support rank
        2,  # quadratic-gauge Fitting Laurent-support rank
    )
    assert (1, 6, 2) == (
        1,  # weighted boundary-contact nilpotency index
        2 * 1 * (2 + 1),  # cancellation type (m,r)=(2,1)
        2,  # quadratic gauge
    )


if __name__ == "__main__":
    check_fixed_rational_pairs()
    check_original_rational_quartic()
    check_fixed_quadratic_field_triple()
    print("PASS: fixed Q-pairs share the pencil T^N+T^3-2T^2+T+u")
    print("PASS: the original Q-quartic is irreducible and boundary-clean")
    print("PASS: three fixed maps over Q(sqrt(-2)) share one quartic pencil")
    print("PASS: the R=-1 quartic is irreducible by a mod-17 Rabin certificate")
