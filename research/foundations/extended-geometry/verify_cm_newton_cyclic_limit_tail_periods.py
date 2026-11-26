#!/usr/bin/env python3
"""Verify the local rescaling and resonant coefficient of the tail examples."""

from sympy import Rational, series, simplify, sqrt, symbols


def check_correction_scale():
    angular_scale = Rational(9, 7)
    height_scale = -Rational(2, 7)
    assert 18*height_scale + 4*angular_scale == 0
    assert 1 + 17*height_scale + 4*angular_scale == Rational(9, 7)
    assert 2 + 16*height_scale + 2*angular_scale == 0
    assert 1 + 12*height_scale + angular_scale == -Rational(8, 7)

    z = symbols("z")
    quadratic = 1 + z*z
    log_q = 8*z/quadratic - 1/z
    assert simplify(z*quadratic.diff(z) - Rational(2, 3)*z*log_q*quadratic
                    - Rational(2, 3)*(1-4*z*z)) == 0


def check_resonance():
    tau = symbols("tau")
    discriminant_root = (2+tau)**2 * sqrt(1+2*tau*tau) / (1+tau)
    truncated = 4 + 5*tau*tau - tau**3
    assert series(discriminant_root, tau, 0, 7).removeO() == (
        truncated-2*tau**5+Rational(7, 2)*tau**6
    )
    square_remainder = (2+tau)**4/(1+2*tau*tau)*(1-truncated**2/discriminant_root**2)
    assert series(square_remainder, tau, 0, 7).removeO() == -16*tau**5-4*tau**6

    def log_derivative_at_zero(function):
        return function.diff(tau).subs(tau, 0)/function.subs(tau, 0)

    even_part = (1+2*tau*tau)/(1+tau)
    square_center = truncated/even_part
    a_log = log_derivative_at_zero(square_center**4*(2+tau)**2/even_part)
    beta_log = log_derivative_at_zero(square_center**18/even_part)
    assert (a_log, beta_log) == (6, 19)
    assert a_log-(Rational(5, 9)*beta_log+Rational(4, 9)*Rational(1, 4))/2 == Rational(2, 3)

    V, beta_0, r_0 = symbols("V beta_0 r_0")
    square = beta_0*V**32-r_0*V**5
    assert simplify(10 - V*square.diff(V)/(2*square) - (
        Rational(15, 2) - Rational(27, 2)*beta_0*V**27/(beta_0*V**27-r_0)
    )) == 0
    assert Rational(15, 54) == Rational(5, 18)


if __name__ == "__main__":
    check_correction_scale()
    check_resonance()
    print("Exact escaping residue scalings and resonant coefficient checked")
