#!/usr/bin/env python3
"""Exact checks for the integer-Beta radial Mathieu counterexample family.

For positive integers s,t, let B_{s,t} be normalized integration against
U^(s-1)(1-U)^(t-1) on [0,1], and put

    f=(1-T^-1)((1-U)+UT).

The theorem checked here says, for m >= s+t-1,

    Phi(f^m) = 0,
    Phi(T^-t f^m)
      = (-1)^(m+t) binom(s+t-1,t) / binom(m+t,t).

The calculation is dependency-free and uses exact rational arithmetic.
"""

from fractions import Fraction
from math import comb, factorial


def beta_bernstein_moment(s: int, t: int, m: int, k: int) -> Fraction:
    """B_{s,t}(binom(m,k) U^k (1-U)^(m-k))."""
    return Fraction(
        comb(m, k)
        * factorial(s + t - 1)
        * factorial(k + s - 1)
        * factorial(m - k + t - 1),
        factorial(s - 1)
        * factorial(t - 1)
        * factorial(m + s + t - 1),
    )


def circuit_moment(s: int, t: int, m: int, angular_shift: int = 0) -> Fraction:
    """Phi(T^-angular_shift f^m), expanded in the angular variable."""
    total = Fraction(0)
    for j in range(m + 1):
        k = j + angular_shift
        if 0 <= k <= m:
            total += (
                (-1) ** j
                * comb(m, j)
                * beta_bernstein_moment(s, t, m, k)
            )
    return total


def mixed_closed_form(s: int, t: int, m: int) -> Fraction:
    return Fraction(
        (-1) ** (m + t) * comb(s + t - 1, t),
        comb(m + t, t),
    )


def polynomial_weight_moment(
    coefficients: list[int], t: int, m: int, angular_shift: int = 0
) -> Fraction:
    """Integral CT(T^-shift f^m) (1-U)^(t-1) h(U) dU."""
    total = Fraction(0)
    for r, coefficient in enumerate(coefficients):
        if not coefficient:
            continue
        normalization = Fraction(
            factorial(r) * factorial(t - 1),
            factorial(r + t),
        )
        total += (
            coefficient
            * normalization
            * circuit_moment(r + 1, t, m, angular_shift)
        )
    return total


def main() -> None:
    for s in range(1, 7):
        for t in range(1, 7):
            threshold = s + t - 1

            # The pure identity starts exactly after the finite-difference
            # degree s+t-2.
            assert circuit_moment(s, t, threshold - 1) != 0
            for m in range(threshold, threshold + 12):
                assert circuit_moment(s, t, m) == 0
                mixed = circuit_moment(s, t, m, angular_shift=t)
                assert mixed == mixed_closed_form(s, t, m)
                assert mixed != 0

            # Endpoint prime ratios of a_n=(s)_t/(n+s)_t are 1 mod p.
            for p in (17, 19, 23, 29, 31):
                if p <= s + t:
                    continue
                for n in range(3):
                    for j in range(n + 1, 4):
                        numerator = 1
                        denominator = 1
                        for h in range(t):
                            numerator *= n * p + s + h
                            denominator *= j * p + s + h
                        assert numerator % p == denominator % p

    # General polynomial density
    # w(U)=(1-U)^(t-1)h(U).  The mixed boundary term depends only on h(1).
    polynomial_cases = [
        ([2, -3, 5], 1),
        ([1, 0, 4, -2], 2),
        ([-3, 7, 1, 0, 2], 4),
    ]
    for coefficients, t in polynomial_cases:
        degree = len(coefficients) - 1
        threshold = degree + t
        endpoint_value = sum(coefficients)
        assert endpoint_value != 0
        for m in range(threshold, threshold + 12):
            assert polynomial_weight_moment(coefficients, t, m) == 0
            mixed = polynomial_weight_moment(
                coefficients, t, m, angular_shift=t
            )
            expected = Fraction(
                (-1) ** (m + t) * endpoint_value,
                t * comb(m + t, t),
            )
            assert mixed == expected
            assert mixed != 0

    print("PASS integer-Beta pure finite differences for 1 <= s,t <= 6")
    print("PASS integer-Beta mixed closed form for 12 moments past threshold")
    print("PASS integer-Beta endpoint ratios are units, not prime-separated")
    print("PASS arbitrary polynomial-density boundary formula")


if __name__ == "__main__":
    main()
