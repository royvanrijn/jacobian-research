"""Exact examples of prescribed moments on covers with zero Frobenius return."""

from collections import defaultdict
from fractions import Fraction
from math import comb, gcd


def root_polynomial(roots, modulus=None):
    coefficients = [1]
    for root, multiplicity in roots:
        for _ in range(multiplicity):
            updated = [0] * (len(coefficients) + 1)
            for index, value in enumerate(coefficients):
                updated[index] -= root * value
                updated[index + 1] += value
            coefficients = updated if modulus is None else [v % modulus for v in updated]
    return coefficients


def parameter_choice(prime, moment, b, lower, upper):
    step = 2 * moment * prime
    ell = 0
    while True:
        numerator = moment * prime ** ell
        source = 1 + step * (int((numerator / upper - 1) // step) + 1)
        if lower < Fraction(numerator, source) < upper:
            break
        ell += 1

    iterate = 1
    while True:
        numerator = source * prime ** iterate - 1
        if numerator % moment == 0:
            order = numerator // moment
            if order > max(2 * source, b) and lower < Fraction(prime ** (iterate + ell), order) < upper:
                return source, ell, iterate, order
        iterate += 1


def coherent_dimension(character, order, branches):
    residues = sum(character * a % order for a in branches)
    assert residues % order == 0
    return residues // order - 1


def check_return_map(prime, moment, finite=None, interval=None):
    if finite is None:
        finite = tuple((root, 1 if root == 1 else 2) for root in range(1, prime))
    b = sum(a for _, a in finite) - 1
    if interval is None:
        lower = Fraction(prime + 1, 2 * prime)
        interval = (lower, lower + Fraction(1, b * prime))
    source, ell, iterate, order = parameter_choice(prime, moment, b, *interval)
    k = prime ** (iterate + ell)
    target = prime * k % order
    assert k < order and gcd(k, order) == 1
    assert sum(a for _, a in finite) == b + 1
    roots = ((0, order - b), *finite)
    branches = tuple(a for _, a in roots) + (order - 1,)
    assert sum(branches) == 2 * order and gcd(order, *branches) == 1
    assert coherent_dimension(k, order, branches) == coherent_dimension(target, order, branches) == 1
    fiber = tuple((root, prime * (a * k % order) // order) for root, a in roots)
    assert root_polynomial(fiber, prime)[prime - 1] == 0

    polynomial = root_polynomial(((root, a * moment) for root, a in finite), prime)
    selected = polynomial[b * moment]
    assert selected != 0 and source * prime ** iterate % order == 1
    print(f"p={prime}, m={moment}, source={source}, iterate={iterate}, "
          f"connected order={order}, return zero, selected entry={selected}.")
    return source, iterate, order, selected, finite


def general_branch_data(prime, moment):
    power = prime
    while power <= moment:
        power *= prime
    constant = 2 if prime == 3 and moment % 3 == 1 else 1
    inverse = constant * pow(moment, -1, power) % power
    x = next(index for index in range(1, prime - 1)
             if 1 <= -inverse * index % prime <= index)
    finite = tuple((root, power - inverse if root == 1 else power)
                   for root in range(1, prime))
    b = sum(a for _, a in finite) - 1
    lower = Fraction(x, prime) + Fraction(1, prime * power)
    # These are the possible breaks for either coherent dimension and
    # for all the finite exponents of the one-step map.
    divisors = (prime, prime * b, prime * power, prime * (power - inverse))
    upper = min(Fraction((lower * d) // 1 + 1, d) for d in divisors)
    return finite, (lower, upper)


def laurent_lift(order, source, finite):
    polynomial = root_polynomial(finite)
    b = len(polynomial) - 2
    result = defaultdict(Fraction)
    for degree, value in enumerate(polynomial):
        angle = degree - b
        if angle >= 0:
            result[angle, order - 2 * source * angle] += value
        else:
            pole = -angle
            for k in range(pole + 1):
                result[angle, order + 2 * source * k] += (
                    value * (-1) ** (pole + k) * comb(pole, k)
                )
    return result


def second_moment(phase):
    return sum(a * b / (j + ell + 1)
               for (i, j), a in phase.items()
               for (k, ell), b in phase.items()
               if i + k == 0 and (j + ell) % 2 == 0)


def check_off_congruence_precision(prime, source, iterate, order, selected, finite):
    assert (prime, source, iterate, order) == (5, 81, 1, 202)
    phase = laurent_lift(order, source, finite)
    lower_height = order % source
    phase[0, lower_height] += Fraction(7, 2)

    # The off-congruence term t^172 meets the leading constant angle at
    # interval denominator 375, which has valuation three.
    forbidden_height = 172
    assert forbidden_height < order and (forbidden_height - order) % source != 0
    assert order + forbidden_height + 1 == 3 * prime ** 3
    phase[0, forbidden_height] = Fraction(prime ** 3)

    difference = source * prime ** iterate * second_moment(phase) - selected
    assert difference.denominator % prime != 0 and difference.numerator % prime == 0
    phase[0, forbidden_height] = Fraction(prime ** 2)
    normalized = source * prime ** iterate * second_moment(phase)
    assert normalized.denominator % prime != 0
    assert normalized.numerator * pow(normalized.denominator, -1, prime) % prime == 2
    print("Off-congruence t^172: precision 5^3 preserves residue 3; 5^2 changes it to 2.")


def main():
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    prime = 5
    check_off_congruence_precision(prime, *check_return_map(prime, 2))
    check_return_map(5, 26)
    check_return_map(7, 3)
    check_return_map(11, 5)
    for prime, moment, expected in ((3, 4, 2), (3, 5, 1), (5, 7, 1), (5, 11, 1)):
        finite, interval = general_branch_data(prime, moment)
        assert check_return_map(prime, moment, finite, interval)[3] == expected


if __name__ == "__main__":
    main()
