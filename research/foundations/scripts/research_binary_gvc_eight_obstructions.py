#!/usr/bin/env python3
"""Exact certificates for the reduced binary GVC return obstructions.

Corollary 5.11 of ``BINARY_GVC_UNIFORM_FACE_TERMINATION.md`` leaves six
quadratic--cubic and two double-quadratic fourth-moment expressions.
The default command evaluates their signs exactly on a bounded endpoint
window.  ``--prove-h00`` verifies Lemma 5.12, while
``--prove-negative-corners`` verifies Lemma 5.13.  Only the two proof
modes just named and ``--prove-three-more`` are unbounded arithmetic
results for those obstructions.  The ``--verify-opposite-packet`` mode
exactly checks the coefficient identity and a bounded window of the
central-binomial inequalities used in the separate unbounded packet
theorem.  A bounded sign table and ``--explore-remaining-cones`` are
evidence.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import comb, factorial

import sympy as sp


QC_TYPES = ((0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1))
DQ_TYPES = ((0, 1), (0, 2))


def factorial_ratio(n: int, multiplier: int) -> int:
    return factorial(multiplier * n) // factorial(n) ** multiplier


def normalized_weights(
    r: int, s: int
) -> tuple[dict[tuple[int, int], Fraction], dict[int, Fraction]]:
    """Divide W_(m,k) by the harmless common factor (r!)^m."""

    w = {
        (m, k): Fraction(
            comb(m, k)
            * (-1) ** k
            * factorial_ratio(r, m - k)
            * factorial_ratio(s, k)
        )
        for m in range(2, 5)
        for k in range(m + 1)
    }
    endpoint = {
        m: sum(Fraction(comb(m, k)) * w[m, k] for k in range(m + 1))
        for m in range(2, 5)
    }
    return w, endpoint


def quadratic_cubic(
    w: dict[tuple[int, int], Fraction | sp.Expr],
    endpoint: dict[int, Fraction | sp.Expr],
    q: int,
    ell: int,
) -> Fraction | sp.Expr:
    u = -endpoint[2] / (2 * w[2, q])
    s_q = w[3, q] + w[3, q + 1]
    v = -(endpoint[3] + 6 * s_q * u) / (3 * w[3, ell])
    a_q = w[4, q] + 2 * w[4, q + 1] + w[4, q + 2]
    return (
        endpoint[4]
        + 12 * a_q * u
        + 6 * w[4, 2 * q] * u**2
        + 12 * (w[4, ell] + w[4, ell + 1]) * v
    )


def double_quadratic(
    w: dict[tuple[int, int], Fraction | sp.Expr],
    endpoint: dict[int, Fraction | sp.Expr],
    q: int,
    q_prime: int,
) -> tuple[Fraction | sp.Expr, Fraction | sp.Expr]:
    a = 2 * w[2, q]
    b = 2 * w[2, q_prime]
    c = 6 * (w[3, q] + w[3, q + 1])
    d = 6 * (w[3, q_prime] + w[3, q_prime + 1])
    determinant = a * d - b * c
    if determinant == 0:
        return determinant, Fraction(0)
    u = (-endpoint[2] * d + b * endpoint[3]) / determinant
    v = (-a * endpoint[3] + c * endpoint[2]) / determinant
    a_q = w[4, q] + 2 * w[4, q + 1] + w[4, q + 2]
    a_q_prime = (
        w[4, q_prime]
        + 2 * w[4, q_prime + 1]
        + w[4, q_prime + 2]
    )
    obstruction = (
        endpoint[4]
        + 12 * a_q * u
        + 6 * w[4, 2 * q] * u**2
        + 12 * a_q_prime * v
        + 6 * w[4, 2 * q_prime] * v**2
        + 12 * w[4, q + q_prime] * u * v
    )
    return determinant, obstruction


def sign(value: Fraction) -> int:
    return (value > 0) - (value < 0)


def symbolic_data() -> tuple[
    dict[str, sp.Symbol],
    dict[str, sp.Expr],
]:
    c_r, c_s, t_r, t_s, q_r, q_s = sp.symbols(
        "C_r C_s T_r T_s Q_r Q_s"
    )
    w = {
        (2, 0): c_r,
        (2, 1): -2,
        (2, 2): c_s,
        (3, 0): t_r,
        (3, 1): -3 * c_r,
        (3, 2): 3 * c_s,
        (3, 3): -t_s,
        (4, 0): q_r,
        (4, 1): -4 * t_r,
        (4, 2): 6 * c_r * c_s,
        (4, 3): -4 * t_s,
        (4, 4): q_s,
    }
    endpoint = {
        m: sum(sp.binomial(m, k) * w[m, k] for k in range(m + 1))
        for m in range(2, 5)
    }
    expressions = {
        f"H{kind}": sp.cancel(quadratic_cubic(w, endpoint, *kind))
        for kind in QC_TYPES
    }
    for kind in DQ_TYPES:
        determinant, obstruction = double_quadratic(w, endpoint, *kind)
        expressions[f"Delta{kind}"] = sp.factor(determinant)
        expressions[f"D{kind}"] = sp.cancel(obstruction)
    symbols = {
        "C_r": c_r,
        "C_s": c_s,
        "T_r": t_r,
        "T_s": t_s,
        "Q_r": q_r,
        "Q_s": q_s,
    }
    return symbols, expressions


def prove_h00() -> None:
    """Verify Lemma 5.12's all-order coefficient-positive certificate."""

    symbols, expressions = symbolic_data()
    c_r, c_s = symbols["C_r"], symbols["C_s"]
    t_r, t_s = symbols["T_r"], symbols["T_s"]
    q_r, q_s = symbols["Q_r"], symbols["Q_s"]
    l_r, l_s, m_r, m_s = sp.symbols("L_r L_s M_r M_s")
    ratio_substitution = {
        t_r: c_r * l_r,
        t_s: c_s * l_s,
        q_r: c_r * l_r * m_r,
        q_s: c_s * l_s * m_s,
    }
    assert sp.expand(
        expressions["Delta(0, 1)"].subs(ratio_substitution)
        - 12 * c_r * (2 * l_r - 3 * c_r + 3 * c_s - 6)
    ) == 0
    assert sp.expand(
        expressions["Delta(0, 2)"].subs(ratio_substitution)
        + 12 * c_r * c_s * (l_r + l_s - 6)
    ) == 0

    n = sp.symbols("n", integer=True, positive=True)
    next_s = {
        c_s: 2 * (2 * n + 1) * c_s / (n + 1),
        t_s: (
            (3 * n + 1)
            * (3 * n + 2)
            * (3 * n + 3)
            * t_s
            / (n + 1) ** 3
        ),
        q_s: (
            sp.prod(4 * n + j for j in range(1, 5))
            * q_s
            / (n + 1) ** 4
        ),
    }
    difference = sp.cancel(
        expressions["H(0, 0)"].subs(next_s, simultaneous=True)
        - expressions["H(0, 0)"]
    )
    difference_numerator = sp.together(difference).as_numer_denom()[0]

    x_r, y_r, z_r, x_s, y_s, z_s, n_0 = sp.symbols(
        "x_r y_r z_r x_s y_s z_s n_0"
    )
    cone = {
        c_r: 2 + x_r,
        l_r: sp.Rational(3, 2) * (2 + x_r) + y_r,
        m_r: sp.Rational(4, 3)
        * (sp.Rational(3, 2) * (2 + x_r) + y_r)
        + z_r,
        c_s: 2 + x_s,
        l_s: sp.Rational(3, 2) * (2 + x_s) + y_s,
        m_s: sp.Rational(4, 3)
        * (sp.Rational(3, 2) * (2 + x_s) + y_s)
        + z_s,
        n: 1 + n_0,
    }
    difference_polynomial = sp.Poly(
        sp.expand(
            difference_numerator.subs(ratio_substitution).subs(cone)
        ),
        x_r,
        y_r,
        z_r,
        x_s,
        y_s,
        z_s,
        n_0,
    )
    assert all(
        coefficient >= 0
        for coefficient in difference_polynomial.coeffs()
    )
    assert difference_polynomial.eval(
        {
            x_r: 0,
            y_r: 0,
            z_r: 0,
            x_s: 0,
            y_s: 0,
            z_s: 0,
            n_0: 0,
        }
    ) > 0

    base_numerator = sp.together(
        expressions["H(0, 0)"].subs({c_s: 2, t_s: 6, q_s: 24})
    ).as_numer_denom()[0]
    base_cone = {
        c_r: 6 + x_r,
        t_r: (6 + x_r)
        * (sp.Rational(5, 2) * (6 + x_r) + y_r),
        q_r: (6 + x_r)
        * (sp.Rational(5, 2) * (6 + x_r) + y_r)
        * (
            sp.Rational(28, 15)
            * (sp.Rational(5, 2) * (6 + x_r) + y_r)
            + z_r
        ),
    }
    base_polynomial = sp.Poly(
        sp.expand(base_numerator.subs(base_cone)),
        x_r,
        y_r,
        z_r,
    )
    assert all(coefficient >= 0 for coefficient in base_polynomial.coeffs())
    assert base_polynomial.eval({x_r: 0, y_r: 0, z_r: 0}) > 0
    assert expressions["H(0, 0)"].subs(
        {c_r: 2, t_r: 6, q_r: 24, c_s: 2, t_s: 6, q_s: 24}
    ) == 0

    print(
        "PASS determinants and H_(0,0): 19-term positive base and "
        f"{len(difference_polynomial.terms())}-term positive "
        "forward-difference certificate"
    )


def prove_negative_corners() -> None:
    """Verify Lemma 5.13's coupled-ratio cone certificates."""

    symbols, expressions = symbolic_data()
    c_r, c_s = symbols["C_r"], symbols["C_s"]
    t_r, t_s = symbols["T_r"], symbols["T_s"]
    q_r, q_s = symbols["Q_r"], symbols["Q_s"]
    l_r, l_s, m_r, m_s = sp.symbols("L_r L_s M_r M_s")
    ratio_substitution = {
        t_r: c_r * l_r,
        t_s: c_s * l_s,
        q_r: c_r * l_r * m_r,
        q_s: c_s * l_s * m_s,
    }

    # Exact one-step inequalities.  Together with C_1=2, L_1=3,
    # M_1=4, these prove C>=2, L/C>=3/2, and M/L>=4/3.
    n = sp.symbols("n", integer=True, positive=True)
    c_step = 2 * (2 * n + 1) / (n + 1)
    l_step = (
        3 * (3 * n + 1) * (3 * n + 2)
        / (2 * (n + 1) * (2 * n + 1))
    )
    m_step = (
        4 * (4 * n + 1) * (4 * n + 2) * (4 * n + 3)
        / (3 * (n + 1) * (3 * n + 1) * (3 * n + 2))
    )
    assert sp.factor(c_step - 3) == (n - 1) / (n + 1)
    assert sp.factor(l_step / c_step - sp.Rational(5, 3)) == (
        (n - 1) * (n + 2) / (12 * (2 * n + 1) ** 2)
    )
    assert sp.factor(m_step / l_step - sp.Rational(7, 5)) == (
        (n - 1)
        * (n + 2)
        * (17 * n**2 + 17 * n + 6)
        / (45 * (3 * n + 1) ** 2 * (3 * n + 2) ** 2)
    )

    # D_(0,2) is invariant under swapping the two endpoints.
    endpoint_swap = {
        c_r: c_s,
        c_s: c_r,
        t_r: t_s,
        t_s: t_r,
        q_r: q_s,
        q_s: q_r,
    }
    assert sp.cancel(
        expressions["D(0, 2)"]
        - expressions["D(0, 2)"].xreplace(endpoint_swap)
    ) == 0

    x, y, z, a, b, d = sp.symbols("x y z a b d")
    lower_c = 2 + x
    lower_l = sp.Rational(3, 2) * lower_c + y
    lower_m = sp.Rational(4, 3) * lower_l + z
    gap_c = 3 + a
    gap_l = sp.Rational(5, 3) * gap_c + b
    gap_m = sp.Rational(7, 5) * gap_l + d
    ordered_cones = {
        "s>r": {
            c_r: lower_c,
            l_r: lower_l,
            m_r: lower_m,
            c_s: lower_c * gap_c,
            l_s: lower_l * gap_l,
            m_s: lower_m * gap_m,
        },
        "r>s": {
            c_s: lower_c,
            l_s: lower_l,
            m_s: lower_m,
            c_r: lower_c * gap_c,
            l_r: lower_l * gap_l,
            m_r: lower_m * gap_m,
        },
    }
    expected_terms = {
        ("H(0, 3)", "s>r"): 266,
        ("H(0, 3)", "r>s"): 361,
        ("D(0, 2)", "s>r"): 2236,
        ("D(0, 2)", "r>s"): 2236,
    }
    term_counts: dict[tuple[str, str], int] = {}
    for name in ("H(0, 3)", "D(0, 2)"):
        numerator, denominator = sp.together(
            -expressions[name].subs(ratio_substitution)
        ).as_numer_denom()
        if name == "H(0, 3)":
            assert denominator == 2 * c_r
        else:
            assert sp.factor(
                denominator
                - 6 * c_r * c_s * (l_r + l_s - 6) ** 2
            ) == 0
        for orientation, cone in ordered_cones.items():
            polynomial = sp.Poly(
                sp.expand(numerator.subs(cone)),
                x,
                y,
                z,
                a,
                b,
                d,
            )
            assert all(
                coefficient >= 0 for coefficient in polynomial.coeffs()
            )
            assert polynomial.eval(
                {x: 0, y: 0, z: 0, a: 0, b: 0, d: 0}
            ) > 0
            term_counts[name, orientation] = len(polynomial.terms())

    assert term_counts == expected_terms
    print(
        "PASS H_(0,3)<0 and D_(0,2)<0 for all unequal positive "
        "endpoints: coupled-ratio coefficient counts "
        f"{term_counts}"
    )


def prove_three_more() -> None:
    """Verify the ordered-tail and factorial-product cone closures."""

    symbols, expressions = symbolic_data()
    c_r, c_s = symbols["C_r"], symbols["C_s"]
    t_r, t_s = symbols["T_r"], symbols["T_s"]
    q_r, q_s = symbols["Q_r"], symbols["Q_s"]
    l_r, l_s, m_r, m_s = sp.symbols("L_r L_s M_r M_s")
    ratio_substitution = {
        t_r: c_r * l_r,
        t_s: c_s * l_s,
        q_r: c_r * l_r * m_r,
        q_s: c_s * l_s * m_s,
    }

    n = sp.symbols("n", integer=True, positive=True)
    c_step = 2 * (2 * n + 1) / (n + 1)
    l_over_c_step = (
        3 * (3 * n + 1) * (3 * n + 2) / (4 * (2 * n + 1) ** 2)
    )
    m_over_l_step = (
        8
        * (4 * n + 1)
        * (4 * n + 2)
        * (4 * n + 3)
        * (2 * n + 1)
        / (9 * (3 * n + 1) ** 2 * (3 * n + 2) ** 2)
    )
    assert sp.factor(
        c_step.subs(n, n + 1)
        - c_step
        - 2 / ((n + 1) * (n + 2))
    ) == 0
    assert sp.factor(
        l_over_c_step.subs(n, n + 1)
        - l_over_c_step
        - 3
        * (n + 1)
        / (2 * (2 * n + 1) ** 2 * (2 * n + 3) ** 2)
    ) == 0
    assert sp.factor(
        m_over_l_step.subs(n, n + 1)
        - m_over_l_step
        - 64
        * (n + 1)
        * (18 * n**4 + 72 * n**3 + 103 * n**2 + 62 * n + 15)
        / (
            9
            * (3 * n + 1) ** 2
            * (3 * n + 2) ** 2
            * (3 * n + 4) ** 2
            * (3 * n + 5) ** 2
        )
    ) == 0
    product_ratio = sp.factor(
        (
            c_step * l_over_c_step
            * c_step * l_over_c_step * m_over_l_step
        )
        / c_step**3
    )
    assert sp.factor(
        product_ratio
        - 1
        - (8 * n**2 + 7 * n + 1) / (2 * (2 * n + 1) ** 3)
    ) == 0

    x, y, z, a, b, d = sp.symbols("x y z a b d")

    def numerator(name: str, multiplier: int) -> sp.Expr:
        signed, denominator = sp.together(
            multiplier * expressions[name].subs(ratio_substitution)
        ).as_numer_denom()
        expected_denominator = {
            "H(0, 1)": 6 * c_r,
            "H(0, 2)": 6 * c_r,
            "H(1, 0)": 4,
            "H(1, 1)": 12,
            "D(0, 1)": (
                6
                * c_r
                * (3 * c_r - 3 * c_s - 2 * l_r + 6) ** 2
            ),
        }[name]
        assert sp.factor(denominator - expected_denominator) == 0
        return signed

    def certify(polynomial: sp.Poly, expected_terms: int) -> None:
        assert len(polynomial.terms()) == expected_terms
        assert all(coefficient >= 0 for coefficient in polynomial.coeffs())
        assert polynomial.eval(
            {variable: 0 for variable in polynomial.gens}
        ) > 0

    def tail_polynomial(
        name: str,
        multiplier: int,
        orientation: str,
        base: int,
    ) -> sp.Poly:
        base_c = sp.Integer(comb(2 * base, base))
        base_l = sp.Integer(comb(3 * base, base))
        base_m = sp.Integer(comb(4 * base, base))
        base_n = sp.Integer(base)
        lower_c = base_c + x
        lower_l = base_l * lower_c / base_c + y
        lower_m = base_m * lower_l / base_l + z
        gap_c = c_step.subs(n, base_n) + a
        gap_l_over_c = l_over_c_step.subs(n, base_n) + b
        gap_m_over_l = m_over_l_step.subs(n, base_n) + d
        gap_l = gap_c * gap_l_over_c
        gap_m = gap_l * gap_m_over_l
        if orientation == "s>r":
            cone = {
                c_r: lower_c,
                l_r: lower_l,
                m_r: lower_m,
                c_s: lower_c * gap_c,
                l_s: lower_l * gap_l,
                m_s: lower_m * gap_m,
            }
        else:
            cone = {
                c_s: lower_c,
                l_s: lower_l,
                m_s: lower_m,
                c_r: lower_c * gap_c,
                l_r: lower_l * gap_l,
                m_r: lower_m * gap_m,
            }
        return sp.Poly(
            sp.expand(numerator(name, multiplier).subs(cone)),
            x,
            y,
            z,
            a,
            b,
            d,
        )

    def fixed_ray_polynomial(
        name: str,
        multiplier: int,
        fixed_side: str,
        fixed: int,
        moving_base: int,
    ) -> sp.Poly:
        fixed_c = sp.Integer(comb(2 * fixed, fixed))
        fixed_l = sp.Integer(comb(3 * fixed, fixed))
        fixed_m = sp.Integer(comb(4 * fixed, fixed))
        moving_c0 = sp.Integer(comb(2 * moving_base, moving_base))
        moving_l0 = sp.Integer(comb(3 * moving_base, moving_base))
        moving_m0 = sp.Integer(comb(4 * moving_base, moving_base))
        moving_c = moving_c0 + x
        moving_l = moving_l0 * moving_c / moving_c0 + y
        moving_m = moving_m0 * moving_l / moving_l0 + z
        if fixed_side == "r":
            cone = {
                c_r: fixed_c,
                l_r: fixed_l,
                m_r: fixed_m,
                c_s: moving_c,
                l_s: moving_l,
                m_s: moving_m,
            }
        else:
            cone = {
                c_s: fixed_c,
                l_s: fixed_l,
                m_s: fixed_m,
                c_r: moving_c,
                l_r: moving_l,
                m_r: moving_m,
            }
        return sp.Poly(
            sp.expand(numerator(name, multiplier).subs(cone)),
            x,
            y,
            z,
        )

    tail_certificates = (
        ("H(0, 1)", 1, "s>r", 2, 456),
        ("H(0, 1)", -1, "r>s", 4, 570),
        ("H(0, 2)", 1, "s>r", 6, 456),
        ("H(1, 1)", 1, "s>r", 2, 340),
        ("H(1, 1)", 1, "r>s", 6, 340),
        # These two oriented certificates narrow the final two formulas.
        ("H(1, 0)", 1, "s>r", 2, 340),
        ("D(0, 1)", 1, "s>r", 2, 1500),
    )
    for name, multiplier, orientation, base, terms in tail_certificates:
        certify(
            tail_polynomial(name, multiplier, orientation, base),
            terms,
        )

    lower_c = 2 + x
    lower_l = sp.Rational(3, 2) * lower_c + y
    lower_m = sp.Rational(4, 3) * lower_l + z
    gap_c = 3 + a
    gap_l = sp.Rational(5, 3) * gap_c + b
    gap_m = sp.Rational(7, 5) * gap_l + d
    h02_reverse = sp.Poly(
        sp.expand(
            numerator("H(0, 2)", -1).subs(
                {
                    c_s: lower_c,
                    l_s: lower_l,
                    m_s: lower_m,
                    c_r: lower_c * gap_c,
                    l_r: lower_l * gap_l,
                    m_r: lower_m * gap_m,
                }
            )
        ),
        x,
        y,
        z,
        a,
        b,
        d,
    )
    certify(h02_reverse, 361)

    fixed_certificates = (
        ("H(0, 1)", 1, "r", 1, 3, 14),
        ("H(0, 1)", -1, "s", 1, 4, 19),
        ("H(0, 1)", -1, "s", 2, 4, 19),
        ("H(0, 1)", -1, "s", 3, 4, 19),
        ("H(0, 2)", 1, "r", 1, 6, 14),
        ("H(0, 2)", 1, "r", 2, 6, 14),
        ("H(0, 2)", 1, "r", 3, 6, 14),
        ("H(0, 2)", 1, "r", 4, 6, 14),
        ("H(0, 2)", 1, "r", 5, 7, 14),
        ("H(1, 1)", 1, "r", 1, 3, 14),
        ("H(1, 1)", 1, "s", 1, 6, 14),
        ("H(1, 1)", 1, "s", 2, 6, 14),
        ("H(1, 1)", 1, "s", 3, 6, 14),
        ("H(1, 1)", 1, "s", 4, 6, 14),
        ("H(1, 1)", 1, "s", 5, 6, 14),
        ("H(1, 0)", 1, "r", 1, 3, 14),
        ("H(1, 0)", -1, "s", 1, 4, 14),
        ("H(1, 0)", -1, "s", 2, 4, 14),
        ("H(1, 0)", -1, "s", 3, 4, 14),
        ("D(0, 1)", 1, "r", 1, 4, 24),
        ("D(0, 1)", -1, "s", 1, 10, 43),
        ("D(0, 1)", -1, "s", 2, 10, 43),
        ("D(0, 1)", -1, "s", 3, 10, 43),
    )
    for (
        name,
        multiplier,
        fixed_side,
        fixed,
        moving_base,
        terms,
    ) in fixed_certificates:
        certify(
            fixed_ray_polynomial(
                name,
                multiplier,
                fixed_side,
                fixed,
                moving_base,
            ),
            terms,
        )

    # On the final wedge r>s>=4, the broad linear cone misses the
    # factorial correlation L_s*M_s >= kappa_4*C_s^3.  Its exact ratio
    # is increasing by the product-ratio identity above.
    product_slack = sp.symbols("product_slack")
    base = 4
    base_c = sp.Integer(comb(2 * base, base))
    base_l = sp.Integer(comb(3 * base, base))
    base_m = sp.Integer(comb(4 * base, base))
    lower_c = base_c + x
    lower_l = base_l * lower_c / base_c + y
    product_constant = base_l * base_m / base_c**3
    lower_m = (
        product_constant * lower_c**3 + product_slack
    ) / lower_l
    gap_c = c_step.subs(n, base) + a
    gap_l = gap_c * (l_over_c_step.subs(n, base) + b)
    gap_m = gap_l * (m_over_l_step.subs(n, base) + d)
    product_cone = {
        c_s: lower_c,
        l_s: lower_l,
        m_s: lower_m,
        c_r: lower_c * gap_c,
        l_r: lower_l * gap_l,
        m_r: lower_m * gap_m,
    }
    for name, terms in (("H(1, 0)", 408), ("D(0, 1)", 1692)):
        numerator(name, -1)
        product_numerator, product_denominator = sp.together(
            -expressions[name]
            .subs(ratio_substitution)
            .subs(product_cone)
        ).as_numer_denom()
        denominator_polynomial = sp.Poly(
            sp.expand(product_denominator),
            x,
            y,
            product_slack,
            a,
            b,
            d,
        )
        assert all(
            coefficient >= 0
            for coefficient in denominator_polynomial.coeffs()
        )
        assert denominator_polynomial.eval(
            {
                x: 0,
                y: 0,
                product_slack: 0,
                a: 0,
                b: 0,
                d: 0,
            }
        ) > 0
        certify(
            sp.Poly(
                sp.expand(product_numerator),
                x,
                y,
                product_slack,
                a,
                b,
                d,
            ),
            terms,
        )

    # Exact finite complements of the tail and fixed-ray cones.
    h01_finite = ((1, 2), (2, 1), (3, 1), (3, 2))
    h02_finite = tuple(
        (r, s)
        for r in range(1, 6)
        for s in range(r + 1, 7)
        if not (r == 5 and s == 6)
    ) + ((5, 6),)
    h11_finite = ((1, 2),) + tuple(
        (r, s) for r in range(2, 6) for s in range(1, r)
    )
    for name, pairs in (
        ("H(0, 1)", h01_finite),
        ("H(0, 2)", h02_finite),
        ("H(1, 1)", h11_finite),
    ):
        q, ell = {
            "H(0, 1)": (0, 1),
            "H(0, 2)": (0, 2),
            "H(1, 1)": (1, 1),
        }[name]
        for r, s in pairs:
            w, endpoint = normalized_weights(r, s)
            assert quadratic_cubic(w, endpoint, q, ell) != 0

    # The same certificates leave only r>s, s>=4 for H10 and D01.
    for r, s in ((1, 2), (2, 1), (3, 1), (3, 2)):
        w, endpoint = normalized_weights(r, s)
        assert quadratic_cubic(w, endpoint, 1, 0) != 0
    for r, s in ((1, 2), (1, 3)):
        w, endpoint = normalized_weights(r, s)
        assert double_quadratic(w, endpoint, 0, 1)[1] != 0
    for s in range(1, 4):
        for r in range(s + 1, 10):
            w, endpoint = normalized_weights(r, s)
            assert double_quadratic(w, endpoint, 0, 1)[1] != 0

    print(
        "PASS H_(0,1), H_(0,2), H_(1,0), H_(1,1), and D_(0,1) "
        "are nonzero for all unequal endpoints; the final wedge uses "
        "408- and 1692-term factorial-product cones"
    )


def bounded_sign_table(limit: int) -> None:
    qc_signs = {kind: set() for kind in QC_TYPES}
    dq_signs = {kind: set() for kind in DQ_TYPES}
    determinant_signs = {kind: set() for kind in DQ_TYPES}
    for r in range(1, limit + 1):
        for s in range(1, limit + 1):
            if r == s:
                continue
            w, endpoint = normalized_weights(r, s)
            for kind in QC_TYPES:
                qc_signs[kind].add(sign(quadratic_cubic(w, endpoint, *kind)))
            for kind in DQ_TYPES:
                determinant, obstruction = double_quadratic(
                    w, endpoint, *kind
                )
                determinant_signs[kind].add(sign(determinant))
                dq_signs[kind].add(sign(obstruction))

    print(f"exact bounded sign table for 1 <= r,s <= {limit}, r != s")
    print(f"quadratic--cubic: {qc_signs}")
    print(f"double--quadratic determinants: {determinant_signs}")
    print(f"double--quadratic obstructions: {dq_signs}")
    assert all(0 not in signs for signs in qc_signs.values())
    assert all(0 not in signs for signs in dq_signs.values())
    assert all(0 not in signs for signs in determinant_signs.values())
    print(
        "STATUS: exact bounded regression; the three proof modes "
        "independently close all eight obstructions"
    )


def positive_compositions_of_three(total: int):
    """Yield the ordered positive triples summing to ``total``."""

    for first in range(1, total - 1):
        for second in range(1, total - first):
            yield first, second, total - first - second


def verify_opposite_packet(limit: int) -> None:
    """Replay the exact identities behind Theorem 7.5.

    The loop is deliberately only a bounded regression.  The theorem's
    unbounded inequalities follow from strict central-binomial
    supermultiplicativity, as recorded in the canonical note.
    """

    a_1, b_1, r_a, r_b, u = sp.symbols("A_1 B_1 R_A R_B u")
    v = -a_1 * u / b_1
    second_coefficient = (
        a_1**2 * r_a * u**2
        + a_1 * b_1 * u * v
        + b_1**2 * r_b * v**2
    )
    expected = a_1**2 * u**2 * (r_a + r_b - 1)
    assert sp.expand(second_coefficient - expected) == 0

    central = lambda n: comb(2 * n, n)
    checked_profiles = 0
    for degree in range(3, limit + 1):
        for profile in positive_compositions_of_three(degree):
            profile_product = 1
            for part in profile:
                profile_product *= central(part)
            assert central(degree) > profile_product
            for endpoint_order in range(1, limit + 1):
                ratio = Fraction(
                    central(degree * endpoint_order),
                    central(degree) * profile_product,
                )
                if endpoint_order == 1:
                    assert ratio == Fraction(1, profile_product)
                    assert ratio <= Fraction(1, 8)
                else:
                    assert ratio > 1
                checked_profiles += 1

    print(
        "PASS opposite three-by-three packet coefficient identity; "
        f"checked {checked_profiles} exact profile/order ratios through "
        f"d,r <= {limit}"
    )
    print(
        "STATUS: bounded regression for Theorem 7.5; its unbounded proof "
        "uses strict central-binomial supermultiplicativity"
    )


def explore_remaining_cones() -> None:
    """Report coefficient signs in the two basic ordered endpoint cones."""

    symbols, expressions = symbolic_data()
    c_r, c_s = symbols["C_r"], symbols["C_s"]
    t_r, t_s = symbols["T_r"], symbols["T_s"]
    q_r, q_s = symbols["Q_r"], symbols["Q_s"]
    l_r, l_s, m_r, m_s = sp.symbols("L_r L_s M_r M_s")
    ratio_substitution = {
        t_r: c_r * l_r,
        t_s: c_s * l_s,
        q_r: c_r * l_r * m_r,
        q_s: c_s * l_s * m_s,
    }
    x, y, z, a, b, d = sp.symbols("x y z a b d")
    lower_c = 2 + x
    lower_l = sp.Rational(3, 2) * lower_c + y
    lower_m = sp.Rational(4, 3) * lower_l + z
    gap_c = 3 + a
    gap_l = sp.Rational(5, 3) * gap_c + b
    gap_m = sp.Rational(7, 5) * gap_l + d
    cones = {
        "s>r": {
            c_r: lower_c,
            l_r: lower_l,
            m_r: lower_m,
            c_s: lower_c * gap_c,
            l_s: lower_l * gap_l,
            m_s: lower_m * gap_m,
        },
        "r>s": {
            c_s: lower_c,
            l_s: lower_l,
            m_s: lower_m,
            c_r: lower_c * gap_c,
            l_r: lower_l * gap_l,
            m_r: lower_m * gap_m,
        },
    }
    for name in (
        "H(0, 1)",
        "H(0, 2)",
        "H(1, 0)",
        "H(1, 1)",
        "D(0, 1)",
    ):
        for multiplier in (1, -1):
            numerator, denominator = sp.together(
                multiplier * expressions[name].subs(ratio_substitution)
            ).as_numer_denom()
            for orientation, cone in cones.items():
                polynomial = sp.Poly(
                    sp.expand(numerator.subs(cone)),
                    x,
                    y,
                    z,
                    a,
                    b,
                    d,
                )
                negative = sum(bool(c < 0) for c in polynomial.coeffs())
                print(
                    name,
                    f"sign={multiplier:+}",
                    orientation,
                    f"negative={negative}/{len(polynomial.terms())}",
                    f"denominator={sp.factor(denominator)}",
                )

    tail_tests = (
        ("H(0, 1)", 1, "s>r", 2),
        ("H(0, 1)", -1, "r>s", 4),
        ("H(0, 2)", 1, "s>r", 6),
        ("H(1, 0)", 1, "s>r", 2),
        ("H(1, 0)", -1, "r>s", 4),
        ("H(1, 0)", -1, "r>s", 5),
        ("H(1, 0)", -1, "r>s", 6),
        ("H(1, 0)", -1, "r>s", 10),
        ("H(1, 1)", 1, "s>r", 2),
        ("H(1, 1)", 1, "r>s", 2),
        ("H(1, 1)", 1, "r>s", 3),
        ("H(1, 1)", 1, "r>s", 4),
        ("H(1, 1)", 1, "r>s", 6),
        ("H(1, 1)", 1, "r>s", 10),
        ("D(0, 1)", 1, "s>r", 2),
        ("D(0, 1)", -1, "r>s", 4),
        ("D(0, 1)", -1, "r>s", 6),
        ("D(0, 1)", -1, "r>s", 10),
    )
    for name, multiplier, orientation, base in tail_tests:
        base_c = sp.Integer(comb(2 * base, base))
        base_l = sp.Integer(comb(3 * base, base))
        base_m = sp.Integer(comb(4 * base, base))
        base_n = sp.Integer(base)
        tail_c = base_c + x
        tail_l = base_l * tail_c / base_c + y
        tail_m = base_m * tail_l / base_l + z
        step_c = 2 * (2 * base_n + 1) / (base_n + 1) + a
        step_l_over_c = (
            3
            * (3 * base_n + 1)
            * (3 * base_n + 2)
            / (4 * (2 * base_n + 1) ** 2)
            + b
        )
        step_m_over_l = (
            8
            * (4 * base_n + 1)
            * (4 * base_n + 2)
            * (4 * base_n + 3)
            * (2 * base_n + 1)
            / (
                9
                * (3 * base_n + 1) ** 2
                * (3 * base_n + 2) ** 2
            )
            + d
        )
        tail_gap_l = step_c * step_l_over_c
        tail_gap_m = tail_gap_l * step_m_over_l
        if orientation == "s>r":
            tail_cone = {
                c_r: tail_c,
                l_r: tail_l,
                m_r: tail_m,
                c_s: tail_c * step_c,
                l_s: tail_l * tail_gap_l,
                m_s: tail_m * tail_gap_m,
            }
        else:
            tail_cone = {
                c_s: tail_c,
                l_s: tail_l,
                m_s: tail_m,
                c_r: tail_c * step_c,
                l_r: tail_l * tail_gap_l,
                m_r: tail_m * tail_gap_m,
            }
        numerator = sp.together(
            multiplier * expressions[name].subs(ratio_substitution)
        ).as_numer_denom()[0]
        polynomial = sp.Poly(
            sp.expand(numerator.subs(tail_cone)),
            x,
            y,
            z,
            a,
            b,
            d,
        )
        negative = sum(bool(c < 0) for c in polynomial.coeffs())
        print(
            name,
            f"tail sign={multiplier:+}",
            orientation,
            f"base>={base}",
            f"negative={negative}/{len(polynomial.terms())}",
        )
        if 0 < negative <= 10:
            print(
                "  negative terms=",
                [
                    term
                    for term in polynomial.terms()
                    if term[1] < 0
                ],
            )

    fixed_ray_tests = (
        ("H(0, 1)", 1, "r", 1, 3),
        ("H(0, 1)", -1, "s", 1, 4),
        ("H(0, 1)", -1, "s", 2, 4),
        ("H(0, 1)", -1, "s", 3, 4),
        ("H(0, 2)", 1, "r", 1, 6),
        ("H(0, 2)", 1, "r", 2, 6),
        ("H(0, 2)", 1, "r", 3, 6),
        ("H(0, 2)", 1, "r", 4, 6),
        ("H(0, 2)", 1, "r", 5, 7),
        ("H(1, 0)", 1, "r", 1, 3),
        ("H(1, 0)", -1, "s", 1, 4),
        ("H(1, 0)", -1, "s", 2, 4),
        ("H(1, 0)", -1, "s", 3, 4),
        ("H(1, 1)", 1, "r", 1, 3),
        ("H(1, 1)", 1, "s", 1, 6),
        ("H(1, 1)", 1, "s", 2, 6),
        ("H(1, 1)", 1, "s", 3, 6),
        ("H(1, 1)", 1, "s", 4, 6),
        ("H(1, 1)", 1, "s", 5, 6),
        ("D(0, 1)", 1, "r", 1, 3),
        ("D(0, 1)", 1, "r", 1, 4),
        ("D(0, 1)", 1, "r", 1, 6),
        ("D(0, 1)", 1, "r", 1, 10),
        ("D(0, 1)", -1, "s", 1, 4),
        ("D(0, 1)", -1, "s", 2, 4),
        ("D(0, 1)", -1, "s", 3, 4),
        ("D(0, 1)", -1, "s", 1, 10),
        ("D(0, 1)", -1, "s", 2, 10),
        ("D(0, 1)", -1, "s", 3, 10),
    )
    for name, multiplier, fixed_side, fixed, moving_base in fixed_ray_tests:
        fixed_c = sp.Integer(comb(2 * fixed, fixed))
        fixed_l = sp.Integer(comb(3 * fixed, fixed))
        fixed_m = sp.Integer(comb(4 * fixed, fixed))
        moving_c0 = sp.Integer(comb(2 * moving_base, moving_base))
        moving_l0 = sp.Integer(comb(3 * moving_base, moving_base))
        moving_m0 = sp.Integer(comb(4 * moving_base, moving_base))
        moving_c = moving_c0 + x
        moving_l = moving_l0 * moving_c / moving_c0 + y
        moving_m = moving_m0 * moving_l / moving_l0 + z
        if fixed_side == "r":
            fixed_ray_cone = {
                c_r: fixed_c,
                l_r: fixed_l,
                m_r: fixed_m,
                c_s: moving_c,
                l_s: moving_l,
                m_s: moving_m,
            }
            moving_side = "s"
        else:
            fixed_ray_cone = {
                c_s: fixed_c,
                l_s: fixed_l,
                m_s: fixed_m,
                c_r: moving_c,
                l_r: moving_l,
                m_r: moving_m,
            }
            moving_side = "r"
        numerator = sp.together(
            multiplier * expressions[name].subs(ratio_substitution)
        ).as_numer_denom()[0]
        polynomial = sp.Poly(
            sp.expand(numerator.subs(fixed_ray_cone)),
            x,
            y,
            z,
        )
        negative = sum(bool(c < 0) for c in polynomial.coeffs())
        print(
            name,
            f"fixed {fixed_side}={fixed}",
            f"{moving_side}>={moving_base}",
            f"sign={multiplier:+}",
            f"negative={negative}/{len(polynomial.terms())}",
        )

    n_0 = sp.symbols("n_0")
    n = sp.symbols("n", integer=True, positive=True)
    next_c_ratio = 2 * (2 * n + 1) / (n + 1)
    next_t_ratio = (
        (3 * n + 1)
        * (3 * n + 2)
        * (3 * n + 3)
        / (n + 1) ** 3
    )
    next_q_ratio = (
        sp.prod(4 * n + j for j in range(1, 5)) / (n + 1) ** 4
    )
    next_r = {
        c_r: c_r * next_c_ratio,
        t_r: t_r * next_t_ratio,
        q_r: q_r * next_q_ratio,
    }
    h10_difference = sp.cancel(
        expressions["H(1, 0)"].subs(next_r, simultaneous=True)
        - expressions["H(1, 0)"]
    )
    difference_numerator = sp.together(
        -h10_difference.subs(ratio_substitution)
    ).as_numer_denom()[0]
    base = 3
    base_c = sp.Integer(comb(2 * base, base))
    base_l = sp.Integer(comb(3 * base, base))
    base_m = sp.Integer(comb(4 * base, base))
    tail_c = base_c + x
    tail_l = base_l * tail_c / base_c + y
    tail_m = base_m * tail_l / base_l + z
    base_n = sp.Integer(base)
    step_c = 2 * (2 * base_n + 1) / (base_n + 1) + a
    step_l_over_c = (
        3
        * (3 * base_n + 1)
        * (3 * base_n + 2)
        / (4 * (2 * base_n + 1) ** 2)
        + b
    )
    step_m_over_l = (
        8
        * (4 * base_n + 1)
        * (4 * base_n + 2)
        * (4 * base_n + 3)
        * (2 * base_n + 1)
        / (
            9
            * (3 * base_n + 1) ** 2
            * (3 * base_n + 2) ** 2
        )
        + d
    )
    difference_cone = {
        c_s: tail_c,
        l_s: tail_l,
        m_s: tail_m,
        c_r: tail_c * step_c,
        l_r: tail_l * step_c * step_l_over_c,
        m_r: tail_m
        * step_c
        * step_l_over_c
        * step_m_over_l,
        n: 4 + n_0,
    }
    difference_polynomial = sp.Poly(
        sp.expand(difference_numerator.subs(difference_cone)),
        x,
        y,
        z,
        a,
        b,
        d,
        n_0,
    )
    print(
        "H(1, 0) negative r-successor for s>=3:",
        f"negative={sum(bool(c < 0) for c in difference_polynomial.coeffs())}"
        f"/{len(difference_polynomial.terms())}",
    )

    boundary_substitution = {
        c_r: c_s * next_c_ratio,
        t_r: t_s * next_t_ratio,
        q_r: q_s * next_q_ratio,
    }
    boundary_numerator = sp.together(
        -expressions["H(1, 0)"]
        .subs(boundary_substitution, simultaneous=True)
        .subs(ratio_substitution)
    ).as_numer_denom()[0]
    boundary_cone = {
        c_s: tail_c,
        l_s: tail_l,
        m_s: tail_m,
        n: 3 + n_0,
    }
    boundary_polynomial = sp.Poly(
        sp.expand(boundary_numerator.subs(boundary_cone)),
        x,
        y,
        z,
        n_0,
    )
    print(
        "H(1, 0) negative adjacent boundary for s>=3:",
        f"negative={sum(bool(c < 0) for c in boundary_polynomial.coeffs())}"
        f"/{len(boundary_polynomial.terms())}",
    )

    strong_base = 4
    strong_c0 = sp.Integer(comb(2 * strong_base, strong_base))
    strong_l0 = sp.Integer(comb(3 * strong_base, strong_base))
    strong_m0 = sp.Integer(comb(4 * strong_base, strong_base))
    strong_c1 = sp.Integer(comb(2 * (strong_base + 1), strong_base + 1))
    strong_l1 = sp.Integer(comb(3 * (strong_base + 1), strong_base + 1))
    strong_l_slope = (
        strong_l1
        - strong_l0 * strong_c1 / strong_c0
    ) / (strong_c1 - strong_c0)
    strong_c = strong_c0 + x
    strong_l = (
        strong_l0 * strong_c / strong_c0 + strong_l_slope * x + y
    )
    strong_m = strong_m0 * strong_l / strong_l0 + z
    strong_n = sp.Integer(strong_base)
    strong_gap_c = (
        2 * (2 * strong_n + 1) / (strong_n + 1) + a
    )
    strong_gap_l = strong_gap_c * (
        3
        * (3 * strong_n + 1)
        * (3 * strong_n + 2)
        / (4 * (2 * strong_n + 1) ** 2)
        + b
    )
    strong_gap_m = strong_gap_l * (
        8
        * (4 * strong_n + 1)
        * (4 * strong_n + 2)
        * (4 * strong_n + 3)
        * (2 * strong_n + 1)
        / (
            9
            * (3 * strong_n + 1) ** 2
            * (3 * strong_n + 2) ** 2
        )
        + d
    )
    strong_cone = {
        c_s: strong_c,
        l_s: strong_l,
        m_s: strong_m,
        c_r: strong_c * strong_gap_c,
        l_r: strong_l * strong_gap_l,
        m_r: strong_m * strong_gap_m,
    }
    for name in ("H(1, 0)", "D(0, 1)"):
        strong_numerator = sp.together(
            -expressions[name].subs(ratio_substitution)
        ).as_numer_denom()[0]
        strong_polynomial = sp.Poly(
            sp.expand(strong_numerator.subs(strong_cone)),
            x,
            y,
            z,
            a,
            b,
            d,
        )
        print(
            name,
            "strong L(C) secant wedge:",
            "negative="
            f"{sum(bool(c < 0) for c in strong_polynomial.coeffs())}"
            f"/{len(strong_polynomial.terms())}",
        )
        if name == "H(1, 0)":
            print(
                [
                    term
                    for term in strong_polynomial.terms()
                    if term[1] < 0
                ]
            )

    diagonal_a, diagonal_b, diagonal_d = sp.symbols(
        "diagonal_a diagonal_b diagonal_d"
    )
    diagonal_c_ratio = 1 + diagonal_a
    diagonal_l_ratio = diagonal_c_ratio + diagonal_b
    diagonal_m_ratio = diagonal_l_ratio + diagonal_d
    diagonal_cone = {
        c_s: strong_c,
        l_s: strong_l,
        m_s: strong_m,
        c_r: strong_c * diagonal_c_ratio,
        l_r: strong_l * diagonal_l_ratio,
        m_r: strong_m * diagonal_m_ratio,
    }
    for name in ("H(1, 0)", "D(0, 1)"):
        diagonal_numerator = sp.together(
            -expressions[name].subs(ratio_substitution)
        ).as_numer_denom()[0]
        diagonal_polynomial = sp.Poly(
            sp.expand(diagonal_numerator.subs(diagonal_cone)),
            x,
            y,
            z,
            diagonal_a,
            diagonal_b,
            diagonal_d,
        )
        print(
            name,
            "diagonal-ratio wedge:",
            "negative="
            f"{sum(bool(c < 0) for c in diagonal_polynomial.coeffs())}"
            f"/{len(diagonal_polynomial.terms())}",
        )

    basic_l = strong_l0 * strong_c / strong_c0 + y
    basic_m = strong_m0 * basic_l / strong_l0 + z
    basic_cone = {
        c_s: strong_c,
        l_s: basic_l,
        m_s: basic_m,
        c_r: strong_c * strong_gap_c,
        l_r: basic_l * strong_gap_l,
        m_r: basic_m * strong_gap_m,
    }
    h10_numerator = sp.together(
        expressions["H(1, 0)"].subs(ratio_substitution)
    ).as_numer_denom()[0]
    m_r_coefficient = sp.diff(h10_numerator, m_r)
    m_s_coefficient = sp.diff(h10_numerator, m_s)
    dominant_block = m_r_coefficient * m_r + m_s_coefficient * m_s
    remainder_block = sp.expand(h10_numerator - dominant_block)
    for label, block in (
        ("-P", -dominant_block),
        ("-(P+9R)", -(dominant_block + 9 * remainder_block)),
    ):
        block_polynomial = sp.Poly(
            sp.expand(block.subs(basic_cone)),
            x,
            y,
            z,
            a,
            b,
            d,
        )
        print(
            "H(1, 0)",
            label,
            "block cone:",
            "negative="
            f"{sum(bool(c < 0) for c in block_polynomial.coeffs())}"
            f"/{len(block_polynomial.terms())}",
        )
    positive_remainder_majorant = (
        9 * c_r**3 * c_s
        + 18 * c_r**2 * c_s**2
        + 72 * c_r**2 * c_s
        + 9 * c_r * c_s**3
        + 72 * c_r * c_s**2
        + 192 * c_r * l_r
        + 576 * c_r
        + 288 * c_s**2
        + 64 * c_s * l_s
    )
    majorant_polynomial = sp.Poly(
        sp.expand(
            (-dominant_block - positive_remainder_majorant).subs(
                basic_cone
            )
        ),
        x,
        y,
        z,
        a,
        b,
        d,
    )
    print(
        "H(1, 0) dominant block minus positive remainder majorant:",
        "negative="
        f"{sum(bool(c < 0) for c in majorant_polynomial.coeffs())}"
        f"/{len(majorant_polynomial.terms())}",
    )
    print(
        [
            term
            for term in majorant_polynomial.terms()
            if term[1] < 0
        ]
    )

    product_slack = sp.symbols("product_slack")
    product_constant = strong_l0 * strong_m0 / strong_c0**3
    product_l = strong_l0 * strong_c / strong_c0 + y
    product_m = (
        product_constant * strong_c**3 + product_slack
    ) / product_l
    product_cone = {
        c_s: strong_c,
        l_s: product_l,
        m_s: product_m,
        c_r: strong_c * strong_gap_c,
        l_r: product_l * strong_gap_l,
        m_r: product_m * strong_gap_m,
    }
    for name in ("H(1, 0)", "D(0, 1)"):
        product_numerator = sp.together(
            -expressions[name]
            .subs(ratio_substitution)
            .subs(product_cone)
        ).as_numer_denom()[0]
        product_polynomial = sp.Poly(
            sp.expand(product_numerator),
            x,
            y,
            product_slack,
            a,
            b,
            d,
        )
        print(
            name,
            "factorial-product cone:",
            "negative="
            f"{sum(bool(c < 0) for c in product_polynomial.coeffs())}"
            f"/{len(product_polynomial.terms())}",
        )

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=40)
    parser.add_argument("--prove-h00", action="store_true")
    parser.add_argument("--prove-negative-corners", action="store_true")
    parser.add_argument("--prove-three-more", action="store_true")
    parser.add_argument("--verify-opposite-packet", action="store_true")
    parser.add_argument("--explore-remaining-cones", action="store_true")
    arguments = parser.parse_args()
    if arguments.prove_h00:
        prove_h00()
    elif arguments.prove_negative_corners:
        prove_negative_corners()
    elif arguments.prove_three_more:
        prove_three_more()
    elif arguments.verify_opposite_packet:
        verify_opposite_packet(arguments.limit)
    elif arguments.explore_remaining_cones:
        explore_remaining_cones()
    else:
        bounded_sign_table(arguments.limit)


if __name__ == "__main__":
    main()
