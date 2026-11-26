"""Replay selected Frobenius certificates at moment orders divisible by p."""

from collections import defaultdict
from fractions import Fraction
from math import comb, gcd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from verify_cm_newton_prescribed_moment_eigenlines import (
    coherent_dimension,
    general_branch_data,
    laurent_lift,
    root_polynomial,
)


def valuation(integer, prime):
    exponent = 0
    while integer % prime == 0:
        integer //= prime
        exponent += 1
    return exponent


def parameters(prime, moment, minimum_gap=0):
    k = moment - 1
    finite, (lo, hi) = general_branch_data(prime, k)
    b = sum(exponent for _, exponent in finite) - 1
    step = 2 * k * prime
    shift = 1
    while True:
        numerator = k * prime ** shift
        lower = numerator / (k + hi)
        upper = numerator / (k + lo)
        source = 1 + step * (((lower - 1) // step) + 1)
        if lower < source < upper:
            break
        shift += 1

    iterate = 1
    while True:
        numerator = source * prime ** iterate - 1
        if numerator % k == 0:
            order = numerator // k
            exponent = iterate + shift
            radial = prime ** exponent - source * prime ** iterate
            if (order > max(2 * prime, b, source) and exponent > valuation(moment, prime)
                    and 0 < radial < order and lo < Fraction(radial + 1, order) < hi
                    and order - radial > minimum_gap):
                return finite, b, source, iterate, order, exponent, radial
        iterate += 1


def multiply(first, second):
    output = defaultdict(int)
    for (i, j), a in first.items():
        for (k, ell), b in second.items():
            output[i + k, j + ell] += a * b
    return {term: value for term, value in output.items() if value}


def moment_with_radial(phase, power, radial_height, radial_coefficient):
    powers = [{(0, 0): 1}]
    for _ in range(power):
        powers.append(multiply(powers[-1], phase))
    result = Fraction(0)
    for radial_count in range(power + 1):
        multiplier = comb(power, radial_count) * radial_coefficient ** radial_count
        shift = radial_count * radial_height
        result += sum((Fraction(multiplier * value, j + shift + 1)
                       for (i, j), value in powers[power - radial_count].items()
                       if i == 0 and (j + shift) % 2 == 0), Fraction(0))
    return result


def residue(value, prime):
    assert value.denominator % prime != 0
    return value.numerator * pow(value.denominator, -1, prime) % prime


def derivative_coefficient(phase, k, height):
    power = {(0, 0): 1}
    for _ in range(k):
        power = multiply(power, phase)
    return sum(j * value * power.get((-i, height - j), 0)
               for (i, j), value in phase.items())


def check_auxiliary_cover(prime, moment, data):
    finite, b, source, iterate, order, exponent, radial = data
    k = moment - 1
    roots = ((0, order - b), *finite)
    branches = tuple(a for _, a in roots) + (order - 1,)
    assert sum(branches) == 2 * order and gcd(order, *branches) == 1
    assert gcd(order, prime) == 1 and source * prime ** iterate % order == 1
    selected = root_polynomial(((root, a * k) for root, a in finite), prime)[b * k]
    assert selected != 0

    # A separate map in the character-one return orbit is zero.
    character = prime ** exponent % order
    assert character == radial + 1
    target = prime * character % order
    assert coherent_dimension(character, order, branches) == 1
    assert coherent_dimension(target, order, branches) == 1
    factor = tuple((root, prime * (character * a % order) // order) for root, a in roots)
    assert root_polynomial(factor, prime)[prime - 1] == 0
    return selected


def lift_phase(order, prime, finite):
    # The preceding checker uses weight j+2*source*i for its lift.
    # Our weight is j+2*p*i.
    return {term: int(value) for term, value in laurent_lift(order, prime, finite).items()
            if value}


def check_leading_variant(prime, moment):
    data = parameters(prime, moment)
    finite, b, _, _, order, exponent, radial = data
    selected = check_auxiliary_cover(prime, moment, data)
    phase = lift_phase(order, prime, finite)
    phase[0, 1] = 7  # same weight residue as order
    phase[0, 2] = prime  # a p-small off-congruence term
    d = valuation(moment, prime)
    normalized = prime ** (exponent - d) * moment_with_radial(phase, moment, radial, 1)
    expected = moment // prime ** d * selected % prime
    assert residue(normalized, prime) == expected != 0
    print(f"p={prime}, m={moment}, connected order={order}, radial height={radial}, "
          f"leading degree={order + (2 * prime - 1) * b}, residue={expected}.")


def check_exception_and_dominant_variant():
    prime = moment = 3
    finite, b, _, _, order, exponent, radial = data = parameters(prime, moment, minimum_gap=25)
    selected = check_auxiliary_cover(prime, moment, data)
    assert selected == 1 and root_polynomial(finite, prime)[b] == 2
    assert order - radial > moment * (2 * prime * b + 1)
    phase = lift_phase(order, prime, finite)
    phase[0, order - 1] = 1
    phase[0, radial + 1] = 2
    normalized = prime ** (exponent - 1)
    assert derivative_coefficient(phase, 2, prime ** exponent - 1) % prime == 1
    assert residue(normalized * moment_with_radial(phase, moment, radial, 1), prime) == 0
    assert residue(normalized * moment_with_radial(phase, moment, radial, 2), prime) == 1
    print(f"Order {order}, arbitrary lower terms: radial coefficients 1 and 2 give residues 0 and 1.")

    height_unit = order * (moment - 1) + 1
    higher_exponent = 2
    while prime ** (higher_exponent - 1) <= moment * height_unit:
        higher_exponent += 1
    dominant_height = prime ** higher_exponent - height_unit
    normalized = prime ** (higher_exponent - 1)
    assert residue(normalized * moment_with_radial(phase, moment, dominant_height, 1), prime) == 1
    print(f"Same arbitrary lower terms, dominant radial height {dominant_height}: residue 1.")


def main():
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    for prime, moment in ((3, 3), (5, 5), (3, 9)):
        check_leading_variant(prime, moment)
    check_exception_and_dominant_variant()


if __name__ == "__main__":
    main()
