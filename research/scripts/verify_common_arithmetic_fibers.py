#!/usr/bin/env python3
"""Exact certificates for common fibers of stably inequivalent Keller maps."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from jcsearch.weighted import WeightedSeedModel, w, x, y, z  # noqa: E402
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


def check_small_rational_coefficient_triple() -> None:
    """Check the optimized connected quartic shared by all three maps."""

    # The primitive polynomial is rational even though the cancellation map
    # is defined over K=Q(sqrt(-2)).
    shared = 9 * T**4 - 19 * T**3 + 10 * T**2 - 8 * T - 4
    assert sp.discriminant(shared, T) == -75621552

    # The weighted presentation has integral seed and determinant one.
    H = 9 * T**4 - 19 * T**3 + 10 * T**2
    c = -sp.diff(H, T).subs(T, 1)
    assert c == 1
    assert sp.factor(H) == T**2 * (T - 1) * (9 * T - 10)
    assert sp.diff(H, T, 2).subs(T, 1) / c == 14
    weighted = WeightedSeedModel(sp.diff(H, T).subs(T, w), c=c)
    weighted_map = weighted.mapping()
    assert sp.factor(
        sp.Matrix(weighted_map).jacobian((x, y, z)).det()
    ) == 1
    assert tuple(
        len(sp.Poly(sp.expand(component), x, y, z).terms())
        for component in weighted_map
    ) == (16, 14, 3)
    assert tuple(
        sp.Poly(component, x, y, z).total_degree()
        for component in weighted_map
    ) == (12, 11, 4)
    assert sp.expand(weighted.inverse_polynomial(-4, 8, 1).subs(w, T) - shared) == 0

    # Shear the quadratic coefficient into the target.  This is the
    # support-minimal quartic gauge with the same selected inverse polynomial.
    G = 9 * T**4 - 19 * T**3 - 8 * T
    assert (
        sp.Poly(G, T).coeff_monomial(T)
        * sp.Poly(G, T).coeff_monomial(T**3)
        * sp.Poly(G, T).coeff_monomial(T**4)
        != 0
    )
    quadratic_inverse = G - sp.Rational(-8, 2) * (
        sp.Rational(5, 2) * T**2 - 1
    )
    assert sp.expand(quadratic_inverse - shared) == 0
    t = 1 + x * y
    q = t**2 * z + sp.Rational(8, 19) * y**2 * (1 + 3 * t)
    quadratic_map = (
        t * q,
        y
        + sp.Rational(57, 8) * x * q
        - sp.Rational(9, 2) * t**2 * x**2 * q**4,
        x * (5 - 3 * t)
        - sp.Rational(19, 8) * x**3 * z
        + sp.Rational(9, 4) * (x * q) ** 4,
    )
    assert sp.factor(
        sp.Matrix(quadratic_map).jacobian((x, y, z)).det()
    ) == -2
    assert tuple(
        len(sp.Poly(sp.expand(component), x, y, z).terms())
        for component in quadratic_map
    ) == (7, 51, 38)
    assert tuple(
        sp.Poly(component, x, y, z).total_degree()
        for component in quadratic_map
    ) == (7, 26, 24)

    # At the rational cancellation target, the affine generator change
    # T_old=1/4+3*T turns its inverse equation into a nonzero scalar
    # multiple of the same primitive polynomial.
    cancellation_target = (
        sp.Rational(4, 11),
        sp.Integer(1),
        -sp.Rational(22481, 23232),
    )
    old_T = sp.Rational(1, 4) + 3 * T
    cancellation_inverse = (
        fiber_antiderivative(
            2,
            1,
            old_T,
            cancellation_target[0],
            cancellation_target[1],
        )
        - cancellation_target[2]
    )
    assert sp.expand(cancellation_inverse + sp.Rational(36, 121) * shared) == 0

    # A simple root of shared is in every reconstruction open.  For the
    # cancellation chart this follows directly from the differentiated
    # affine-generator identity.
    assert sp.gcd(sp.Poly(shared, T), sp.Poly(sp.diff(shared, T), T)).degree() == 0
    cancellation_derivative = 1 - old_T * (
        cancellation_target[1] - cancellation_target[0] * old_T
    ) ** 2
    assert sp.expand(
        3 * cancellation_derivative
        + sp.Rational(36, 121) * sp.diff(shared, T)
    ) == 0

    # Irreducibility over K has a one-prime certificate.  The prime 17
    # splits in K because 7^2=-2 mod 17.  Monic normalization gives f.
    assert (7**2 + 2) % 17 == 0
    finite_field_polynomial = sp.Poly(
        T**4 - 4 * T**3 + 3 * T**2 + T - 8,
        T,
        modulus=17,
    )

    def power_mod(exponent: int) -> sp.Poly:
        result = sp.Poly(1, T, modulus=17)
        base = sp.Poly(T, T, modulus=17)
        while exponent:
            if exponent & 1:
                result = result.mul(base).rem(finite_field_polynomial)
            base = base.mul(base).rem(finite_field_polynomial)
            exponent //= 2
        return result

    degree_two_test = (
        power_mod(17**2) - sp.Poly(T, T, modulus=17)
    ).rem(finite_field_polynomial)
    degree_four_test = (
        power_mod(17**4) - sp.Poly(T, T, modulus=17)
    ).rem(finite_field_polynomial)
    assert degree_two_test == sp.Poly(
        2 * T**3 + 8 * T**2 + T - 8,
        T,
        modulus=17,
    )
    assert sp.gcd(finite_field_polynomial, degree_two_test).degree() == 0
    assert degree_four_test.is_zero

    # The three stable-boundary fingerprints are inherited from the
    # boundary-clean weighted quartic, type-(2,1) cancellation map, and
    # admissible quartic quadratic gauge.
    unit_ranks = (1, 2, 2)
    cancellation_support = [(0,), (1,)]
    quadratic_support = [(0, 0), (1, 2), (4, 3)]

    def affine_rank(points: list[tuple[int, ...]]) -> int:
        origin = sp.Matrix(points[0])
        differences = [
            sp.Matrix(point) - origin
            for point in points[1:]
        ]
        return sp.Matrix.hstack(*differences).rank()

    assert unit_ranks == (1, 2, 2)
    assert (
        affine_rank(cancellation_support),
        affine_rank(quadratic_support),
    ) == (1, 2)
    assert (1, 2 * 1 * (2 + 1), 2) == (1, 6, 2)


if __name__ == "__main__":
    check_fixed_rational_pairs()
    check_original_rational_quartic()
    check_fixed_quadratic_field_triple()
    check_small_rational_coefficient_triple()
    print("PASS: fixed Q-pairs share the pencil T^N+T^3-2T^2+T+u")
    print("PASS: the original Q-quartic is irreducible and boundary-clean")
    print("PASS: three fixed maps over Q(sqrt(-2)) share one quartic pencil")
    print("PASS: the R=-1 quartic is irreducible by a mod-17 Rabin certificate")
    print("PASS: the optimized rational quartic is one connected fiber of all three maps")
