"""Verify the coefficient and leading-shape tests for (14,4,11)."""

from itertools import product

from sympy import Poly, Rational, binomial, cancel, diff, expand, factorial, gcd, invert
from sympy import rem, resultant, symbols


def quartic_coefficients(a, b, c, power, maximum):
    coefficients = [0] * (maximum + 1)
    for orders in product(range(maximum // 2 + 1), range(maximum // 3 + 1),
                          range(maximum // 4 + 1)):
        weight = sum(w * i for w, i in zip((2, 3, 4), orders))
        if weight > maximum:
            continue
        count = sum(orders)
        term = binomial(power, count) * factorial(count)
        for parameter, order in zip((a, b, c), orders):
            term *= parameter**order / factorial(order)
        coefficients[weight] += term
    return coefficients


def check_quartic_tetradecic():
    z, a, b, delta, alpha, beta, gamma = symbols("z a b Delta alpha beta gamma")
    c = (a**2 - delta) / 4
    quartic = z**4 + a*z**2 + b*z + c
    root_powers = [(14, Rational(7, 2), 1), (6, Rational(3, 2), alpha / 2),
                   (2, Rational(1, 2), beta / 2), (-2, Rational(-1, 2), gamma / 2)]
    root_terms = {i: 0 for i in range(-3, 15)}
    for highest, power, prefactor in root_powers:
        coefficients = quartic_coefficients(a, b, c, power, highest + 3)
        for i, coefficient in enumerate(coefficients):
            root_terms[highest - i] += prefactor * coefficient
    larger = sum(root_terms[i] * z**i for i in range(15))
    first = 35*delta**3 - 280*a*b**2*delta + 192*alpha*delta + 56*b**4 - 512*beta
    second = (2240*a**2*b**4 + 3920*a*b**2*delta**2 - 3072*alpha*a*b**2
              - 3456*alpha*delta**2 - 665*delta**4 + 8192*beta*delta + 16384*gamma)
    gap_factor = 56*a*b**2 - 64*alpha - 35*delta**2
    assert expand(2048*root_terms[-1] + b*first) == 0
    assert expand(32768*root_terms[-2] - 20*delta*first - second) == 0
    assert expand(4096*root_terms[-3] - a*b*first - 2*b**3*gap_factor) == 0

    S, P = symbols("S P")
    differences = [0, 1]
    for i in range(2, 15):
        differences.append(expand(S*differences[-1] - P*differences[-2]))
    divided = sum(root_terms[i]*differences[i] for i in range(1, 15))
    divided = divided.subs(P, (S**3 + a*S + b) / (2*S))
    beta_relation = (first + 512*beta) / 512
    k = S**2 + 2*a
    sum_polynomial = (8*S**5*k**3 + 28*delta*S**3*k**2 - 8*b*S**4*k**2
                      - 28*b*delta*S**2*k - 16*b**2*S**3*k - gap_factor*S*k
                      + 8*b**3*S**2 - 28*b**2*delta*S + b*gap_factor)
    assert expand(-512*divided.subs(beta, beta_relation) / S**2 - sum_polynomial) == 0

    y, X = symbols("kappa X")
    shape = 64*y**2 + 112*y - 19
    assert gcd(shape, (8*y-1)*(8*y-5)) == 1
    delta_shape = rem((5*y - Rational(5, 8))*invert(y**2, shape, y), shape, y)
    b_square_shape = rem(y*delta_shape**2, shape, y)
    assert delta_shape == (1600*y + 2040) / 361
    assert b_square_shape == (70400*y + 32000) / 6859
    leading_quartic = quartic.subs(a, 1)
    leading_larger = larger.subs({a: 1, alpha: 0, beta: 0})
    quartic_prime = diff(leading_quartic, z)
    larger_remainder = rem(diff(leading_larger, z), quartic_prime, z)
    critical_resultant = resultant(quartic_prime, larger_remainder, z)
    reduced = 0
    for (b_degree, delta_degree), coefficient in Poly(critical_resultant, b, delta).terms():
        assert b_degree % 2 == 1
        reduced += rem(coefficient*b_square_shape**(b_degree // 2)
                       * delta_shape**delta_degree, shape, y)
    desired = -Rational(46648000000000000, 104127350297911241532841) * (
        20161770792104*y - 6158358133417)
    assert expand(reduced - desired) == 0
    assert gcd(reduced, shape) == 1

    beta_only = first.subs({alpha: 0, a: y*delta**2/b**2})
    beta_only = beta_only.subs(delta**3, -8192*beta / (35*shape))
    expected = 56*(X**2 - Rational(64, 7)*beta*(64*y**2-16*y-3)/shape)
    assert cancel(beta_only.subs(b**4, X**2) - expected) == 0
    assert gcd(shape, diff(shape, y)) == 1
    numerator = (8*y + 1)*(8*y - 3)
    assert expand(numerator) == 64*y**2 - 16*y - 3
    assert gcd(shape, numerator) == 1


if __name__ == "__main__":
    check_quartic_tetradecic()
    print("Verified the (14,4,11) tail coefficients, pair sums and leading forms.")
