#!/usr/bin/env python3
"""Exact inequalities and finite ordinary moment pairs for support {4,14}."""
import hashlib
import json
import math
import platform
from fractions import Fraction as Q
from pathlib import Path

import flint


A, D, GAP = 4, 14, 10
FIRST, STOP = 13, 4000
PRIMES = (1000003, 1000033)
LAST = D * (STOP - 1)
E_LOW, E_HIGH = Q(163, 60), Q(87, 32)


def check(condition, label):
    if not condition:
        raise ValueError(label)


def exp_upper(x):
    """Rational envelope for exp(x) when x is an integer or half-integer."""
    x = Q(x)
    integer = x.numerator // x.denominator
    remainder = x - integer
    check(remainder in (0, Q(1, 2)), "Exponential bound input")
    return (E_HIGH if integer >= 0 else E_LOW) ** integer * (
        1 if remainder == 0 else Q(5, 3)
    )


def low_envelope(x):
    return (Q(x) / A) ** A * exp_upper(A - x)


def upper_envelope(x):
    return (Q(x) / D) ** D * exp_upper(D - x)


def upper_ratio(x):
    return Q(6, 5) * (Q(D, A) ** A) * E_LOW ** (-GAP) * (D / Q(x)) ** GAP


def low_ratio(x):
    return Q(5, 4) * (Q(A, D) ** A) * E_HIGH ** GAP * (Q(x) / D) ** GAP


def verify_bounds():
    # From the exponential series, 163/60 is its sum through degree five.
    partial = sum(Q(1, math.factorial(k)) for k in range(7))
    check(partial + Q(1, math.factorial(7)) * Q(8, 7) < E_HIGH, "Upper e")
    check(E_HIGH < Q(25, 9), "Upper exp(1/2)")
    check(upper_ratio(Q(55, 4)) < Q(1, 100), "Upper saddle disk")
    check(low_ratio(Q(65, 16)) < Q(1, 1200), "Lower saddle disk")
    check(Q(1, 100) * Q(41, 4) < Q(1, 4), "Upper Rouche margin")
    check(Q(1, 1200) * Q(161, 16) < Q(1, 16), "Lower Rouche margin")

    upper_diff = Q(1, 4) * (126 + Q(13, 4)) / (14 * Q(55, 4) ** 2)
    lower_diff = Q(1, 16) * (44 + Q(3, 16)) / (4 * Q(63, 16) ** 2)
    check(upper_diff < min(Q(1, 14) - Q(1, 20), Q(1, 10) - Q(1, 14)),
          "Upper saddle curvature")
    check(lower_diff < Q(1, 20), "Lower saddle curvature")

    for name, base, x, v, ratio, target in (
        ("upper", 14, 12, Q(1, 4), upper_ratio(12), Q(1, 30)),
        ("lower", 4, 3, Q(1, 16), low_ratio(Q(501, 100)), Q(1, 8)),
    ):
        check(1 - 3 * (v / x) ** 2 > Q(99, 100), name + " real inverse square")
        check(ratio < (Q(1, 25) if base == 14 else Q(1, 150)), name + " strip ratio")
        y = ratio / (1 - ratio)
        end = 16 if base == 14 else 5
        check(Q(99, 100) * base - 110 * y - 100 * y ** 2 > end ** 2 * target,
              name + " strip curvature")
    check(25 + Q(1, 16) ** 2 < Q(501, 100) ** 2, "Lower strip modulus")

    for name, base, x, ratio, bound, radius, half_curvature, decay in (
        ("upper", 14, Q(53, 4), upper_ratio(Q(53, 4)), Q(1, 250),
         Q(1, 2), Q(1, 40), Q(1, 45)),
        ("lower", 4, Q(59, 16), low_ratio(Q(69, 16)), Q(1, 25),
         Q(1, 4), Q(1, 10), Q(1, 12)),
    ):
        y = ratio / (1 - ratio)
        check(2 * base + 1320 * y + 3300 * y ** 2 + 2000 * y ** 3
              < 6 * bound * x ** 3, name + " cubic remainder")
        check(half_curvature - bound * radius > decay, name + " local decay")

    # Each mixed envelope is bounded by separate monotone endpoint values.
    for l, r in ((9, 10), (10, 11), (11, 12), (16, 16)):
        check(upper_envelope(r) + Q(6, 5) * low_envelope(l) < Q(189, 200),
              "Upper remote interval")
    for l, r in ((3, 3), (5, 6), (6, 7), (7, 8), (8, 9)):
        check(low_envelope(l) + Q(5, 4) * upper_envelope(r) < Q(93, 100),
              "Lower remote interval")
    check(Q(7, 1296) < Q(1, 180), "Upper imaginary correction")
    check(Q(7, 1024) < Q(1, 140), "Lower imaginary correction")
    check(Q(189, 200) < 1 - Q(1, 180) - Q(1, 80) - Q(1, 30),
          "Upper remote height loss")
    check(Q(93, 100) < 1 - Q(1, 140) - Q(1, 800) - Q(1, 20),
          "Lower remote height loss")

    for l, r in ((4, 5), (5, 6), (6, 7), (7, 8), (8, 9)):
        bound = (Q(1) if l == 4 else low_envelope(l)) * Q(4, 5) + upper_envelope(r)
        check(bound < Q(82, 100), "Small parameter opposite packet")
    for l, r in ((9, 12), (12, 16), (16, 16)):
        bound = low_envelope(l) + Q(5, 6) * (
            Q(1) if l == 12 else upper_envelope(r)
        )
        check(bound < Q(9, 10), "Large parameter opposite packet")
    check(Q(82, 100) < Q(9, 10) * (1 - Q(1, 80)), "Small parameter decay")
    check(Q(9, 10) < Q(92, 100) * (1 - Q(1, 800)), "Large parameter decay")

    check(E_HIGH ** 5 / (1 - Q(1, 4)) < 200, "Upper remote amplitude")
    check(E_HIGH ** 4 / (1 - Q(1, 16)) < 60, "Lower remote amplitude")
    check(E_HIGH ** 2 / (1 - Q(1, 4)) < 10, "Upper central amplitude")
    check(E_HIGH / (1 - Q(1, 16)) < 3, "Lower central amplitude")
    m = STOP - 1
    check(m > 63 ** 2 and m < 64 ** 2 and m // 256 >= 15, "Error budget input")
    upper_error = ((Q(5, 3) * (Q(81, 10) + 45) + Q(1290, 2 ** 15)) / (7 * 63)
                   + Q(800 * 64, 7 * 2 ** 100))
    lower_error = ((Q(4, 3) * (Q(144, 25) + 12) + Q(240, 2 ** 15)) / (4 * 63)
                   + Q(600 * 64, 4 * 2 ** 100))
    check(max(upper_error, lower_error) < Q(1, 4), "Weighted Gaussian errors")
    check(3 * 64 * Q(92, 100) ** STOP < 1, "Opposite packet dominance")

    saddle_u = (Q(55, 56) ** 14) * Q(99, 100)
    saddle_l = Q(6, 5) * (Q(65, 64) ** 4) * Q(1201, 1200) * E_LOW ** -10
    check(saddle_u > 4 * saddle_l, "Stationary value separation")


def factorial_tables(p):
    factorial = [1] * (LAST + 1)
    inv_factorial = [1] * (LAST + 1)
    for k in range(1, LAST + 1):
        factorial[k] = k * factorial[k - 1] % p
    inv_factorial[LAST] = pow(factorial[LAST], -1, p)
    for k in range(LAST, 0, -1):
        inv_factorial[k - 1] = k * inv_factorial[k] % p
    return factorial, inv_factorial


def moment(m, p, tables):
    factorial, inv_factorial = tables
    shape = A * m
    normalization = factorial[m] * inv_factorial[shape] % p
    polynomial = flint.nmod_poly([
        normalization * factorial[shape + GAP * k]
        * inv_factorial[k] * inv_factorial[m - k] % p
        for k in range(m + 1)
    ], p)
    check(polynomial.degree() == m and polynomial[0] == 1, "Actual moment degree")
    return polynomial


def verify_finite():
    for p in PRIMES:
        check(p > LAST and all(p % k for k in range(2, math.isqrt(p) + 1)),
              "Degree-preserving prime")
    tables = [factorial_tables(p) for p in PRIMES]
    for p, table in zip(PRIMES, tables):
        for m in (1, 2, 3, 12, 13):
            actual = moment(m, p, table)
            expected = [
                (math.comb(m, k) * (math.factorial(A * m + GAP * k) // math.factorial(A * m))) % p
                for k in range(m + 1)
            ]
            check(list(actual) == expected, "Literal coefficient control")
    previous = moment(FIRST - 1, PRIMES[0], tables[0])
    stream = hashlib.sha256()
    fallback = []
    count = 0
    for n in range(FIRST, STOP):
        used = PRIMES[0]
        upper = moment(n, used, tables[0])
        gcd, x, y = previous.xgcd(upper)
        if gcd.degree() != 0:
            fallback.append({"terminal": n, "first_gcd_degree": gcd.degree()})
            used = PRIMES[1]
            lower_at_second = moment(n - 1, used, tables[1])
            upper_at_second = moment(n, used, tables[1])
            gcd, x, y = lower_at_second.xgcd(upper_at_second)
            check(gcd == 1 and x * lower_at_second + y * upper_at_second == 1,
                  "Fallback Bezout identity at N=" + str(n))
        else:
            check(gcd == 1 and x * previous + y * upper == 1,
                  "Bezout identity at N=" + str(n))
        stream.update(f"{n},{used}\n".encode())
        previous = upper
        count += 1
    check(count == STOP - FIRST, "Complete finite range")
    return {
        "first_terminal": FIRST,
        "last_terminal": STOP - 1,
        "terminal_pairs": count,
        "maximum_factorial_argument": LAST,
        "primes": PRIMES,
        "fallbacks": fallback,
        "case_stream_sha256": stream.hexdigest(),
    }


def main():
    verify_bounds()
    result = {
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python": platform.python_version(),
        "python_flint": flint.__version__,
        "analytic_rational_inequalities": "verified",
        "finite_bezout_replay": verify_finite(),
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
