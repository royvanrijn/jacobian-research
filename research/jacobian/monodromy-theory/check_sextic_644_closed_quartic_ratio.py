"""Exact first-weight and binomial comparisons for the closed quartic ratio.

Run with .venv/bin/python from the repository root. The checked
identities are the finite coefficients used in
SEXTIC_KELLER_644_CLOSED_QUARTIC_RATIO.md.
"""

from fractions import Fraction

from sympy import Poly, Rational, factor, series, symbols


def check_candidate_weights():
    weights = {
        Fraction(transverse, 4 - x_degree)
        for x_degree in range(3)
        for transverse in range(1, 4 - x_degree)
    }
    assert weights == {
        Fraction(1, 2),
        Fraction(1, 3),
        Fraction(2, 3),
        Fraction(1, 4),
        Fraction(3, 4),
    }
    for weight in weights:
        r, s = weight.numerator, weight.denominator
        assert weight < Fraction(4, 5)
        for n in range(1, 4 // s + 1):
            assert (n * r) % 4 != 0
            assert 1 + n * (s - r) < 11
            assert 1 - n * s < 0
        if s in (3, 4):
            test_order = 3 if s == 3 else 2
            assert test_order * (s - r) < 11
            assert 6 - test_order * s < 0
            assert (test_order * r) % 4 != 0
            if s == 3 and r == 2:
                assert 6 - 2 * s == 0


def check_binomial_coefficients():
    t, linear, quadratic, correction = symbols("t L M D")
    short_edge = Poly(
        series((1 + linear*t)**Rational(3, 2), t, 0, 4).removeO(), t
    )
    assert short_edge.nth(2) == 3*linear**2/8
    assert short_edge.nth(3) == -linear**3/16

    edge = 1 + linear*t + quadratic*t**2
    expansion = series(
        edge**Rational(3, 2) + correction*t**4*edge**Rational(-1, 2),
        t,
        0,
        7,
    ).removeO()
    coefficients = Poly(expansion, t)
    delta = linear**2 - 4*quadratic
    expected_correction = -3*delta**2/128
    assert factor(coefficients.nth(4) - correction + expected_correction) == 0
    assert factor(
        coefficients.nth(6).subs(correction, expected_correction) + delta**3/512
    ) == 0


if __name__ == "__main__":
    check_candidate_weights()
    check_binomial_coefficients()
    print("Checked the first weights and quartic binomial coefficients.")
