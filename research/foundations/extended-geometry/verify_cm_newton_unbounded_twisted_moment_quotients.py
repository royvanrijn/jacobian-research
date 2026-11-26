#!/usr/bin/env python3
"""Replay twisted quotient dimensions, certifying classes, and prime-power coefficients."""

from functools import cache
from math import comb, gcd


def require(condition, failure):
    if not condition:
        raise AssertionError(failure)


def binomial_mod(top, bottom, prime):
    if bottom < 0 or bottom > top:
        return 0
    result = 1
    while top:
        top, top_digit = divmod(top, prime)
        bottom, bottom_digit = divmod(bottom, prime)
        result = result * comb(top_digit, bottom_digit) % prime
    return result


def direct_coefficient(prime, conductor, e, f, c, iterate):
    """Take the coefficient directly in U, without translating the roots."""
    order = (prime**iterate - 1) // conductor
    carries = {0: 1}
    for _ in range(iterate):
        order, digit = divmod(order, prime)
        a, b = e * digit, f * digit
        left = [binomial_mod(a, i, prime) * (-1)**((a - i) % 2)
                for i in range(a + 1)]
        right = [binomial_mod(b, i, prime) * pow(-c, b - i, prime) % prime
                 for i in range(b + 1)]

        @cache
        def local_coefficient(index):
            return sum(left[i] * right[index - i]
                       for i in range(max(0, index - b), min(a, index) + 1)) % prime

        next_carries = {}
        for carry, coefficient in carries.items():
            first_degree = (prime - 1 - carry) % prime
            for degree in range(first_degree, a + b + 1, prime):
                outgoing = (carry + degree) // prime
                total = next_carries.get(outgoing, 0) + coefficient * local_coefficient(degree)
                if total % prime:
                    next_carries[outgoing] = total % prime
                else:
                    next_carries.pop(outgoing, None)
        carries = next_carries
    return carries.get(0, 0)


def check_case(multiplicity, prime, c, expected_residue):
    conductor = 7 * multiplicity
    e, f = multiplicity + 1, 6 * multiplicity
    require(multiplicity % 14 == 5 and c % 4 == 2, "parameters")
    a, inverse = 1 + 3 * multiplicity, 1 + 2 * multiplicity
    residue_class = prime % conductor
    require(residue_class in (a, inverse), "prime progression")
    require(gcd(conductor, e) == 1 and c % prime not in (0, 1), "branch parameters")
    orbit = tuple(pow(prime, j, conductor) for j in range(3))
    expected_orbit = (1, a, inverse) if residue_class == a else (1, inverse, a)
    require(orbit == expected_orbit, "orbit")

    def dimensions(twist):
        exponents = (twist * multiplicity, e, f, -twist * multiplicity - e - f)
        return tuple(sum((index * power) % conductor for power in exponents) // conductor - 1
                     for index in orbit)

    require(dimensions(0) == (1, 1, 1), "selected orbit")
    require(dimensions(1) == (1, 1, 1), "toric quotient orbit")
    require(dimensions(3) == (1, 2, 2), "nonordinary orbit")
    order = (prime**3 - 1) // conductor
    direct = direct_coefficient(prime, conductor, e, f, c, len(orbit))
    translated = pow(1 - c, order, prime) * binomial_mod(f * order, order, prime) % prime
    require(direct == translated == expected_residue and direct != 0, "moment residue")
    print(f"n={conductor}, p={prime}, c={c}: order {order}, residue {direct}, "
          f"quotient dimensions {dimensions(3)}")


def cyclic_orbit(residue_class, conductor):
    orbit = [1]
    index = residue_class
    while index != 1:
        orbit.append(index)
        index = index * residue_class % conductor
    return orbit


def check_classes_at_35():
    conductor, multiplicity = 35, 5
    units = [i for i in range(1, conductor) if gcd(i, conductor) == 1]
    classes = [
        residue_class for residue_class in units
        if all(i + multiplicity * (i % 7) < conductor
               for i in cyclic_orbit(residue_class, conductor))
    ]
    require(classes == [1, 4, 8, 9, 11, 16, 22, 29], "certifying classes at n=35")

    def dimension(twist, index):
        exponents = (twist * multiplicity, 6, 30, -twist * multiplicity - 36)
        return sum((index * exponent) % conductor for exponent in exponents) // conductor - 1

    for residue_class in units:
        mismatch = any(
            dimension(twist, i) != dimension(twist, i * residue_class % conductor)
            for i in units for twist in range(7)
        )
        require(mismatch == (residue_class != 1), "quotient dimensions at n=35")
    print(f"n={conductor}: certifying classes {classes}, {len(classes) - 1} nonsplit")

    for prime, multiple, expected_residue in (
        (29, 1, 1), (43, 1, 14), (79, 1, 67), (3, 1, 0), (29, 2, 1)
    ):
        iterate = multiple * len(cyclic_orbit(prime % conductor, conductor))
        q = prime**iterate
        order = (q - 1) // conductor
        direct = direct_coefficient(prime, conductor, 6, 30, 2, iterate)
        translated = pow(-1, order, prime) * binomial_mod(30 * order, order, prime) % prime
        require(direct == translated == expected_residue, "prime-power residue at n=35")
        require(bool(direct) == (prime % conductor in classes), "class prediction at n=35")
        print(f"n={conductor}, p={prime}, iterate={iterate}: order {order}, residue {direct}")


if __name__ == "__main__":
    check_case(5, 191, 2, 182)
    check_case(19, 191, 6, 186)
    check_case(5, 151, 2, 23)
    check_case(19, 571, 6, 46)
    check_classes_at_35()
