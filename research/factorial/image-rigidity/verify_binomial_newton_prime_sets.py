#!/usr/bin/env python3
"""Verify the finite certificates for the fixed-prime Newton polygon obstruction."""

from fractions import Fraction
from math import comb, factorial, isqrt, prod

A = 6
D = 14
B = A + D
PRIME_CUTOFF = 3 * A + 2 * D

# A row gives p, the lower order, and (slope numerator, denominator, width).
SEEDS = (
    (2, 6, ((27, 2, 2),)),
    (3, 15, ((13, 2, 2),)),
    (5, 7, ((4, 1, 2),)),
    (7, 2, ((2, 1, 2),)),
    (11, 3, ((3, 2, 2),)),
    (13, 5, ((1, 1, 2),)),
    (17, 2, ((1, 1, 2),)),
    (19, 2, ((1, 1, 2),)),
    (23, 4, ((0, 1, 1), (1, 1, 1))),
    (29, 5, ((0, 1, 1), (1, 1, 1))),
    (31, 4, ((1, 2, 2),)),
    (37, 12, ((1, 2, 2),)),
    (41, 7, ((0, 1, 2),)),
    (43, 8, ((0, 1, 2),)),
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def valuation(integer, prime):
    exponent = 0
    while integer % prime == 0:
        exponent += 1
        integer //= prime
    return exponent


def factorial_valuation(argument, prime):
    exponent = 0
    while argument:
        argument //= prime
        exponent += argument
    return exponent


def moment_coefficients(order):
    denominator = factorial(A * order)
    return [
        comb(order, index) * (factorial(A * order + D * index) // denominator)
        for index in range(order + 1)
    ]


def check_reverse_tail(order, prime, slope, minimum):
    difference = slope - Fraction(D, prime - 1)
    support_valuation = valuation(B, prime)

    # Beyond this power, the lower bound increases from one band to the next.
    power = 1
    exponent = 0
    while (difference * power - 2 * exponent - support_valuation <= minimum
           or difference * (prime - 1) * power < 2):
        power *= prime
        exponent += 1

    endpoint = B * order
    first_canceled_index = endpoint // D + 1
    fixed = (
        factorial_valuation(order, prime)
        - factorial_valuation(endpoint, prime)
        - support_valuation
    )
    count = 0
    for index in range(max(order + 1, first_canceled_index), power):
        extended_valuation = (
            fixed
            + factorial_valuation(index - 1 - order, prime)
            - factorial_valuation(index, prime)
            - factorial_valuation(D * index - 1 - endpoint, prime)
        )
        require(
            extended_valuation + slope * index > minimum,
            f"reverse tail ties or undercuts slope {slope} at p={prime}, m={order}, j={index}",
        )
        count += 1
    return count


def certify(seeds=SEEDS):
    primes = tuple(
        candidate
        for candidate in range(2, PRIME_CUTOFF + 1)
        if all(candidate % divisor for divisor in range(2, isqrt(candidate) + 1))
    )
    require(tuple(row[0] for row in seeds) == primes, "incomplete prime coverage")

    edge_count = 0
    tail_count = 0
    for prime, lower_order, selections in seeds:
        require(lower_order >= 2, f"invalid seed order at p={prime}")
        require(sum(row[2] for row in selections) == 2, f"invalid widths at p={prime}")
        slopes = []
        for numerator, denominator, width in selections:
            slope = Fraction(numerator, denominator)
            require(width > 0 and width % slope.denominator == 0,
                    f"incompatible width for slope {slope} at p={prime}")
            require(slope not in slopes, f"duplicate slope {slope} at p={prime}")
            require(slope != Fraction(D, prime - 1), f"mean slope at p={prime}")
            slopes.append(slope)

            for order in (lower_order, lower_order + 1):
                values = [valuation(coefficient, prime)
                          for coefficient in moment_coefficients(order)]
                weighted = [value - slope * k for k, value in enumerate(values)]
                minimum = min(weighted)
                first = weighted.index(minimum)
                last = len(weighted) - 1 - weighted[::-1].index(minimum)
                require(last - first >= width,
                        f"insufficient width for slope {slope} at p={prime}, m={order}")
                edge_count += 1

                if slope > Fraction(D, prime - 1):
                    reverse_minimum = minimum - values[-1] + slope * order
                    tail_count += check_reverse_tail(order, prime, slope, reverse_minimum)

    return len(primes), edge_count, tail_count


def check_residual_example():
    order = 146603499986155091679122137738
    require(order.bit_length() == 97 and order + 1 < 2**99, "residual order bound")
    examples = (
        (79, 625251304, 4, 78, 58, 20, 1, 0),
        (83, 3566867498, 4, 76, 67, 1, 21, 17),
        (97, 579419484714, 5, 77, 31, 51, 51, 0),
    )
    compared = 0
    for prime, remainder, degree, lower, linear, upper, root, lift in examples:
        require(order % prime**(degree + 1) == remainder, "wrong base order")
        stop = prime**2 // 20
        require(19 * (order + 1) < prime**18, "wrong factorial upper bound")
        require(Fraction(prime - 21, prime - 1) * (stop + 1) - 36 > 1,
                "reverse tail can contribute to the example")

        residuals = []
        for m, expected in ((order, {0: 1, 2 * degree: lower}),
                            (order + 1, {0: 1, 2: linear, 2 * degree + 2: upper})):
            exponent = 0
            unit = 1
            edge = {0: 1}
            next_digit = {}
            for j in range(1, stop + 1):
                numerator = m - j + 1
                denominator = j * prod(range(19 * m - 10 * j + 1,
                                              19 * m - 10 * (j - 1) + 1))
                num_v = valuation(numerator, prime)
                den_v = valuation(denominator, prime)
                exponent += num_v - den_v
                unit = unit * (numerator // prime**num_v) * pow(
                    denominator // prime**den_v, -1, prime
                ) % prime
                weight = 2 * exponent + j
                require(weight >= 0, f"edge undercut at p={prime}, m={m}, j={j}")
                if weight == 0:
                    edge[j] = unit
                if weight == 1:
                    next_digit[j] = unit
                compared += 1
            require(edge == expected, f"unexpected edge residual at p={prime}")
            edge_value = sum(c * pow(root, j // 2, prime) for j, c in edge.items())
            require(edge_value % prime == 0, f"residual does not vanish at p={prime}")
            derivative = sum(j * c * pow(root, j // 2, prime) for j, c in edge.items())
            correction = sum(c * pow(root, (j - 1) // 2, prime)
                             for j, c in next_digit.items())
            residuals.append((derivative, correction))

        require((upper * pow(lower, -1, prime) - linear) * root % prime == 1,
                f"could not reduce the residual gcd to degree one at p={prime}")
        (derivative_lower, correction_lower), (derivative_upper, correction_upper) = residuals
        determinant = derivative_lower * correction_upper - derivative_upper * correction_lower
        require(determinant * pow(linear * root, -1, prime) % prime == lift,
                f"incorrect first lifting value at p={prime}")
    return len(examples), compared


def check_lifting_constants():
    numerators = []
    for delta in (-3, -1, 1, 3):
        for i in range(max(1, delta + 1), min(9, 9 + delta) + 1):
            m = Fraction(-i, 19)
            ell = Fraction(-delta, 10)
            r = 10 + delta - i
            tau = 38 * (m - ell) / ((ell + 1) * (ell + 2) * (delta - i))
            b0 = m / prod(range(-i - 9, -i + 1))
            b1 = (m + 1) / prod(range(10 - i, 20 - i))
            w0 = (m - ell) / ((ell + 1) * prod(range(r - 19, r - 9)))
            w1 = (m - 1) / (3 * prod(range(-i - 10, -i)))
            w2 = (m - ell - 1) / ((ell + 3) * prod(range(r - 20, r - 10)))
            value = -ell * ((tau - 1) * b1 + w1 - tau * w2)
            value -= (2 - (ell + 2) * tau) * (b0 - w0)
            numerators.append(abs(value.numerator))
    require(len(numerators) == 28 and all(numerators), "lifting constants can vanish")
    require(max(numerators) == 418112777, "unexpected lifting numerator bound")
    return len(numerators)


if __name__ == "__main__":
    primes, edges, tail_values = certify()
    print(f"Verified {primes} primes, {edges} seed edges, and {tail_values} reverse tail values.")
    primes, compared = check_residual_example()
    constants = check_lifting_constants()
    print(f"Verified the residual and first lifting digit at {primes} primes using "
          f"{compared} reverse coefficients, and {constants} exact lifting constants.")
