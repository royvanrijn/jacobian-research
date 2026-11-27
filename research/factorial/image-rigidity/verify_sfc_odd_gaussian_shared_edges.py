#!/usr/bin/env python3
"""Exact controls for the two-layer Hermite shared-edge obstruction.

The unrestricted theorem is proved in the accompanying note. The finite
checks here audit the digit valuation formula, the complete near-face
inventory, coefficient cancellations, perturbation margins and controls.
No rational extended gcd or enumeration of a longer start prefix is used.
"""

import json
import platform
from fractions import Fraction as F
from math import comb, factorial, isqrt

from verify_sfc_odd_gaussian_prime_divisors import (
    factorial_tables, gcd_degree, initial_polynomial,
    positive_root_valuations, pure_hermite, verify_first,
)


def is_prime(p):
    return p >= 2 and all(p % j for j in range(2, isqrt(p) + 1))


def digit_sum(value, p):
    result = 0
    while value:
        value, digit = divmod(value, p)
        result += digit
    return result


def layer_certificate(p, e, k, d, even, odd, corrupt=False):
    q = (2 * k + 1) * p ** e
    power = p ** (d - 1)
    left, right = (power + 1) // 2, (p * power + 1) // 2
    middle = left + power
    lam = F(power + 1, power * (p - 1))
    alpha = lam / 2
    gamma = even[0][left] + left * lam
    delta = F(2, p - 1)
    threshold = gamma + delta

    def low_layer(coefficients, shift):
        return {i + shift: F(v) + (i + shift) * lam
                for i, v in enumerate(coefficients[0])
                if F(v) + (i + shift) * lam <= threshold}

    expected_even = {left: gamma, right: gamma, middle: threshold}
    if threshold >= 0:
        expected_even[0] = F(0)
    assert low_layer(even, 0) == expected_even
    assert low_layer(odd, 1) == {
        left: gamma, right: gamma, middle: threshold,
        left + 1: gamma + lam, right + 1: gamma + lam,
    }
    assert low_layer(even, 1) == {left + 1: gamma + lam, right + 1: gamma + lam}

    assert even[0][left] == odd[0][left - 1]
    kappa = even[1][left] * pow(odd[1][left - 1], -1, p) % p
    assert kappa == pow(2, (2 - left) % (p - 1), p)
    # Exact c = E_(left+1)/E_left - O_left/O_(left-1).
    # The first ratio has positive valuation and the second is a unit.
    assert even[0][left + 1] > even[0][left]
    assert odd[0][left] == odd[0][left - 1]
    c = -odd[1][left] * pow(odd[1][left - 1], -1, p) % p
    assert c == p - 2
    if corrupt:
        kappa = (kappa + 1) % p
        assert kappa

    residual_layer = {}
    for i in range(len(even[0]) + 1):
        terms = []
        if i < len(even[0]):
            terms.append((even[0][i], even[1][i]))
        if i:
            terms.append((even[0][i - 1], -c * even[1][i - 1]))
            terms.append((odd[0][i - 1], -kappa * odd[1][i - 1]))
        minimum = min(v for v, _ in terms)
        weight = minimum + i * lam
        if weight > threshold:
            continue
        residue = sum(unit for v, unit in terms if v == minimum) % p
        if residue:
            residual_layer[i] = (weight, residue)
        else:
            # One higher integral valuation is already beyond the retained
            # layer; no unrecorded higher-precision coefficient is inferred.
            assert weight + 1 > threshold

    expected = {middle: (threshold, -even[1][middle] % p)}
    if threshold >= 0:
        expected[0] = (F(0), 1)
    assert residual_layer == expected
    assert threshold != 0
    minimum_weight = min(weight for weight, _ in residual_layer.values())
    lowest = [i for i, (weight, _) in residual_layer.items() if weight == minimum_weight]
    assert len(lowest) == 1
    beta = F(e, 2) - sum((F(1, p ** j) for j in range(1, e + 1)), F(0)) / 2
    margin = gamma + beta - alpha - minimum_weight
    assert beta > alpha and margin > 0
    # This was a genuine collision in the previous leading-initial test.
    assert gcd_degree(initial_polynomial(*even, alpha),
                      initial_polynomial(*odd, alpha), p) > 0
    return dict(edge=d, alpha=str(alpha), beta=str(beta), gamma=str(gamma),
                residual_weight=str(minimum_weight), lowest_index=lowest[0],
                perturbation_margin=str(margin))


def verify_case(p, e, k):
    assert is_prime(p) and p > 4 * k + 2 and e >= 2
    assert pow(2, (p - 1) // 2, p) == 1
    q = (2 * k + 1) * p ** e
    values, units = factorial_tables(2 * q + 4, p)
    even = pure_hermite(q + 1, p, values, units)
    odd = pure_hermite(q + 2, p, values, units)
    for i in range(len(even[0])):
        indicators = sum(i % p ** j == (p ** j + 1) // 2 for j in range(1, e + 1))
        assert even[0][i] == -values[i] - indicators
        assert values[i] == (i - digit_sum(i, p)) // (p - 1)
    rows = [layer_certificate(p, e, k, d, even, odd) for d in range(2, e + 1)]
    return dict(prime=p, exponent=e, cofactor=2 * k + 1, start=q, layers=rows)


def perturbation_coefficients():
    # Independently generate every coefficient of the quadratic phase at
    # representative unit orders and check the exact pure-coefficient ratio.
    count = 0
    for m in (2, 3, 6, 9, 14, 50, 51):
        def moment(power):
            out = F(1)
            for j in range(1, power, 2):
                out *= F(j, m)
            return out
        for i in range(m + 1):
            for j in range(m - i + 1):
                linear = m - i - j
                if linear % 2:
                    continue
                actual = F(factorial(m), factorial(i) * factorial(j) * factorial(linear))
                actual *= moment(linear + 2 * j)
                pure = comb(m, i + j) * moment(linear)
                ratio = F(comb(i + j, j))
                for r in range(1, j + 1):
                    ratio *= F(linear + 2 * r - 1, m)
                assert actual == pure * ratio
                count += 1
    return count


def exact_start_49():
    def coefficients(m, odd):
        return [F((2 * m) ** i * factorial(25),
                  factorial(25 - i) * factorial(2 * i + odd))
                for i in range(26)]
    even, odd = coefficients(50, 0), coefficients(51, 1)
    kappa = even[4] / odd[3]
    c = even[5] / even[4] - odd[4] / odd[3]
    assert kappa == F(34375000, 132651) and c == F(-47, 6)
    weights = []
    for i in range(27):
        value = even[i] if i < 26 else F(0)
        if i:
            value -= c * even[i - 1] + kappa * odd[i - 1]
        if value:
            numerator, denominator = abs(value.numerator), value.denominator
            valuation = 0
            while numerator % 7 == 0:
                numerator //= 7
                valuation += 1
            while denominator % 7 == 0:
                denominator //= 7
                valuation -= 1
            weights.append((i, valuation + F(4 * i, 21)))
    assert min(weight for _, weight in weights) == 0
    assert [i for i, weight in weights if weight == 0] == [0]


def main():
    # Different primes, shared-edge depths and nontrivial dominant cofactors.
    inputs = [(7, 2, 0), (7, 3, 0), (7, 4, 0), (7, 5, 0),
              (17, 2, 0), (17, 3, 0), (23, 2, 0), (31, 2, 0),
              (41, 2, 0), (7, 3, 1), (17, 2, 2), (23, 2, 4)]
    cases = [verify_case(*case) for case in inputs]
    # The older leading test does not close the p=7 prime at start 49.
    assert not verify_first(49, 7)
    values, units = factorial_tables(102, 7)
    even, odd = [pure_hermite(m, 7, values, units) for m in (50, 51)]
    try:
        layer_certificate(7, 2, 0, 2, even, odd, corrupt=True)
    except AssertionError:
        pass
    else:
        raise AssertionError("A corrupted syzygy was accepted")
    assert verify_first(25, 5)  # Existing nonresidue branch is a positive control.
    coefficient_checks = perturbation_coefficients()
    exact_start_49()
    # A separate final-slope exception remains: S_2 is proportional to
    # Y^2+20Y+60, and S_2(Y),S_2(2Y) share Y=2 modulo 13.
    residual_first, residual_second = [8, 7, 1], [8, 1, 4]
    assert gcd_degree(residual_first, residual_second, 13) == 1
    assert all(sum(c * 2 ** i for i, c in enumerate(poly)) % 13 == 0
               for poly in (residual_first, residual_second))
    print(json.dumps(dict(
        status="PASS", scope="finite exact controls of the unrestricted written shared-edge theorem",
        python=platform.python_version(),
        cases=cases, shared_edges=sum(len(c["layers"]) for c in cases),
        perturbation_coefficient_checks=coefficient_checks,
        controls="old leading collision retained; corrupted kappa rejected; nonresidue case passes",
        exact_rational_control="start 49 at prime 7, full 27-coefficient syzygy has unique constant least term",
        remaining_final_slope_example="q=5*13^e; S_2(Y) and S_2(2Y) share Y=2 modulo 13, not an actual Gaussian common-zero certificate",
        assurance="internal replay; not independent mathematical review or formal verification",
    ), indent=2))


if __name__ == "__main__":
    main()
