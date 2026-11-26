"""Exact polynomial identities used in the written hyperelliptic octic proof."""

from sympy import Poly, Rational, diff, expand, integrate, symbols


def check_degree_six_generator():
    variable, quadratic, linear, constant = symbols("T p s r")
    cubic = variable**3 + quadratic * variable**2 + linear * variable + constant
    octic_root = (
        variable**4
        + Rational(3, 2) * quadratic * variable**3
        + Rational(3, 8) * (quadratic**2 + 4 * linear) * variable**2
        - (quadratic**3 - 12 * quadratic * linear - 24 * constant) * variable / 16
        + Rational(3, 128)
        * (quadratic**4 - 8 * quadratic**2 * linear + 32 * quadratic * constant
           + 16 * linear**2)
    )
    difference = Poly(variable * octic_root**2 - cubic**3, variable)
    assert difference.degree() == 4
    discriminant = quadratic**2 - 4 * linear
    assert expand(
        difference.nth(4) - Rational(3, 128)
        * discriminant * (quadratic * discriminant + 16 * constant)
    ) == 0

    residual = Poly(difference.as_expr() - difference.nth(3) * cubic, variable)
    factor = (
        25 * quadratic**5 - 128 * quadratic**3 * linear
        + 408 * quadratic**2 * constant + 112 * quadratic * linear**2
        - 480 * linear * constant
    )
    assert expand(residual.nth(2) + discriminant * factor / 1024) == 0
    assert expand(residual.nth(1).subs(linear, quadratic**2 / 4)) == 0
    assert expand(
        factor.subs(constant, -quadratic * discriminant / 16)
        + quadratic * discriminant**2 / 2
    ) == 0
    assert expand(
        octic_root.subs({quadratic: 0, constant: 0})
        - variable**4 - Rational(3, 2) * linear * variable**2
        - Rational(3, 8) * linear**2
    ) == 0


def check_constant_conductor_and_periods():
    coefficients = [Rational(1)]
    for index in range(4):
        coefficients.append(coefficients[-1] * (Rational(11, 2 * index + 3) - 1))
    assert coefficients == [1, Rational(8, 3), Rational(16, 5),
                            Rational(64, 35), Rational(128, 315)]

    variable, endpoint, linear = symbols("u r s")
    for index in range(5):
        odd_generator = variable**(2 * index + 1)
        canonical = integrate(-diff(variable**2, variable) * odd_generator,
                              (variable, -endpoint, endpoint)) / endpoint
        assert expand(canonical + 4 * endpoint**(2 * index + 2) / (2 * index + 3)) == 0
        cubic_in_square = variable**6 + linear * variable**2
        alternative = integrate(-diff(cubic_in_square, variable) * odd_generator,
                                (variable, -endpoint, endpoint)) / endpoint
        predicted = -4 * (3 * endpoint**(2 * index + 6) / (2 * index + 7)
                          + linear * endpoint**(2 * index + 2) / (2 * index + 3))
        assert expand(alternative - predicted) == 0


if __name__ == "__main__":
    check_degree_six_generator()
    check_constant_conductor_and_periods()
    print("Verified octic coefficient identities and both canonical action formulas.")
