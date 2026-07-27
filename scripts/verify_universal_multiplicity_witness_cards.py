#!/usr/bin/env python3
"""Exact degree 4, 5, and 6 universal-multiplicity witness cards."""

from __future__ import annotations

import sympy as sp


T, W, S = sp.symbols("T W S")


def irreducible_mod_prime(polynomial: sp.Expr, prime: int, degree: int) -> None:
    """Check a one-factor irreducibility certificate modulo ``prime``."""

    _, factors = sp.factor_list(polynomial, T, modulus=prime)
    assert len(factors) == 1
    factor, multiplicity = factors[0]
    assert multiplicity == 1
    assert sp.Poly(factor, T, modulus=prime).degree() == degree


# ---------------------------------------------------------------------------
# Degree four: one connected quartic field, three weighted presentations.
# ---------------------------------------------------------------------------

quartic = T**4 - 3 * T**2 - 1
irreducible_mod_prime(quartic, 7, 4)
assert sp.discriminant(quartic, T) == -2704

# The trace-zero generator has Tr(eta^2)=6.  These are three rational points
# of e^2+2u^2=3, equivalently Tr(eta^2)=2e^2+4u^2.
quartic_rows = (
    {
        "e": sp.Rational(5, 3),
        "u": sp.Rational(-1, 3),
        "alpha": sp.Rational(-7, 10),
        "target": (sp.Rational(107, 3430), sp.Rational(5, 49), 1),
    },
    {
        "e": sp.Rational(5, 3),
        "u": sp.Rational(1, 3),
        "alpha": sp.Rational(-3, 10),
        "target": (sp.Rational(107, 270), sp.Rational(-5, 9), 1),
    },
    {
        "e": sp.Integer(-1),
        "u": sp.Integer(-1),
        "alpha": sp.Rational(1, 2),
        "target": (sp.Rational(-3, 2), sp.Integer(-1), 1),
    },
)

quartic_alphas: list[sp.Expr] = []
for row in quartic_rows:
    e = row["e"]
    u = row["u"]
    step = sp.factor(e - 2 * u)
    alpha = row["alpha"]
    assert e**2 + 2 * u**2 == 3
    assert e != 0 and step != 0
    assert sp.factor(u / e - sp.Rational(1, 2)) == alpha

    inverse = sp.factor(
        -quartic.subs(T, u + step * W) / (2 * step**3 * e)
    )
    constant = sp.factor(inverse.subs(W, 0))
    linear = sp.factor(sp.diff(inverse, W).subs(W, 0))
    seed = sp.factor(inverse - constant - linear * W)
    expected_seed = sp.factor(
        W**2 * (W - 1) * (alpha * W - alpha - 1)
    )
    assert sp.factor(seed - expected_seed) == 0
    assert (constant, -linear, 1) == row["target"]

    # Exact-double, distinct primitive roots, Hessian-clean, marked-point
    # clean.  The additional alpha=1 exclusion belongs to the declared
    # boundary-clean weighted open.
    assert alpha not in (0, 1, -1, 2)
    assert 4 * alpha**2 + 4 * alpha + 3 != 0
    assert sp.diff(seed, W).subs(W, 1) == -1
    quartic_alphas.append(alpha)

assert len(set(quartic_alphas)) == 3


# ---------------------------------------------------------------------------
# Degrees five and six: three translations of one connected field each.
# ---------------------------------------------------------------------------

quintic = T**5 + T**3 + 1
sextic = T**6 + T**4 + 1
irreducible_mod_prime(quintic, 2, 5)
irreducible_mod_prime(sextic, 7, 6)
assert sp.discriminant(quintic, T) == 3233
assert sp.discriminant(sextic, T) == -61504

expected_quintic = {
    1: {
        "coefficients": (8, 13, 11, 5, 1),
        "target_c": sp.Rational(-3, 4),
        "invariant": sp.Rational(64, 171875),
    },
    2: {
        "coefficients": (92, 86, 41, 10, 1),
        "target_c": sp.Rational(-41, 46),
        "invariant": sp.Rational(529, 2562500),
    },
    3: {
        "coefficients": (432, 279, 91, 15, 1),
        "target_c": sp.Rational(-271, 216),
        "invariant": sp.Rational(256, 1421875),
    },
}

expected_sextic = {
    1: {
        "coefficients": (10, 21, 24, 16, 6, 1),
        "target_c": sp.Rational(-3, 5),
        "invariant": sp.Rational(4, 9),
    },
    2: {
        "coefficients": (224, 264, 168, 61, 12, 1),
        "target_c": sp.Rational(-81, 112),
        "invariant": sp.Rational(61, 144),
    },
    3: {
        "coefficients": (1566, 1269, 552, 136, 18, 1),
        "target_c": sp.Rational(-811, 783),
        "invariant": sp.Rational(34, 81),
    },
}


def verify_quadratic_card(
    polynomial: sp.Expr,
    degree: int,
    expected_rows: dict[int, dict[str, object]],
) -> None:
    invariant_values: list[sp.Expr] = []
    for shift, expected_row in expected_rows.items():
        translated = sp.Poly(
            sp.expand(polynomial.subs(T, shift + S) - polynomial.subs(T, shift)),
            S,
        )
        coefficients = tuple(
            translated.coeff_monomial(S**index)
            for index in range(1, degree + 1)
        )
        assert coefficients == expected_row["coefficients"]
        assert all(coefficient != 0 for coefficient in coefficients[2:])
        g1 = coefficients[0]
        target_c = sp.factor(-2 * polynomial.subs(T, shift) / g1)
        assert target_c == expected_row["target_c"]

        inverse = sp.expand(translated.as_expr() - g1 * target_c / 2)
        assert sp.expand(inverse - polynomial.subs(T, shift + S)) == 0

        normalized = tuple(
            sp.factor(coefficient / g1) for coefficient in coefficients
        )
        if degree == 5:
            invariant = sp.factor(
                normalized[4] ** 5
                / (normalized[2] * normalized[3] ** 6)
            )
        else:
            invariant = sp.factor(
                normalized[3] * normalized[5] / normalized[4] ** 2
            )
        assert invariant == expected_row["invariant"]
        invariant_values.append(invariant)

    assert len(set(invariant_values)) == len(expected_rows)


verify_quadratic_card(quintic, 5, expected_quintic)
verify_quadratic_card(sextic, 6, expected_sextic)

print("PASS: connected quartic witness has three clean weighted presentations")
print("PASS: connected quintic witness has three distinct quadratic invariants")
print("PASS: connected sextic witness has three distinct top-jet invariants")
print("PASS: every displayed inverse polynomial defines the fixed field")
