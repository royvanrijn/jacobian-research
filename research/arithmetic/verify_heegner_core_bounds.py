#!/usr/bin/env python3
"""Exact rational certificates for the Heegner split-prime core sieve."""

import hashlib
import json
from bisect import bisect_right
from fractions import Fraction
from itertools import product
from math import isqrt, prod
from pathlib import Path

from verify_primitive_norm_minus11 import cubert_floor, divisor_frontier, prime_flags

if not __debug__:
    raise RuntimeError("Run this certificate checker with Python assertions active.")

PRIME_LIMIT = 2**16
ROUNDING = 2**64
VALUATION_CLASSES = (0, 1, 2)  # 2 denotes valuation at least two.


def ceil_scaled(numerator, denominator):
    return (ROUNDING * numerator + denominator - 1) // denominator


def sieve_envelopes(p, inert_primes):
    classes = (p - 1) // 2
    progression_tail = classes * (
        Fraction(1, PRIME_LIMIT**2) + Fraction(1, 2 * p * PRIME_LIMIT)
    )
    mass = Fraction(
        sum(ceil_scaled(2, q * (q - 1)) for q in inert_primes), ROUNDING
    ) + Fraction(2 * PRIME_LIMIT, PRIME_LIMIT - 1) * progression_tail
    suffix = [0] * (len(inert_primes) + 1)
    for i in reversed(range(len(inert_primes))):
        suffix[i] = suffix[i + 1] + ceil_scaled(1, inert_primes[i] ** 2)
    return mass, [Fraction(x, ROUNDING) + progression_tail for x in suffix]


def delta_envelope_sixth(odd_primes, valuations):
    bound = Fraction(3, 2) ** 6
    for q in odd_primes:
        if q >= 64:
            break
        valuation = valuations.get(q)
        exponents = {
            None: range(11), 0: (0,), 1: (1,), 2: range(2, 11)
        }[valuation]
        bound *= max(
            Fraction((a + 1)**6 * q**(5*a), ((q**(a + 1) - 1) // (q - 1))**6)
            for a in exponents
        )
    return bound


def upward_sixth_root(sixth_power):
    hundredths = 1
    while Fraction(hundredths, 100)**6 < sixth_power:
        hundredths += 1
    return Fraction(hundredths, 100)


def inert_19_profiles(inert_primes):
    profiles = []
    for v3, v13 in product(VALUATION_CLASSES, repeat=2):
        valuations = {3: v3, 13: v13}
        bound = delta_envelope_sixth(inert_primes, valuations)
        profiles.append((valuations, upward_sixth_root(bound)))
    return profiles


def check_19_intervals(split_primes, inert_primes, mass, tails):
    profiles = inert_19_profiles(inert_primes)
    assert all(bound < 5 for _, bound in profiles)
    # A boundary target is 4q or 20q for one split prime q > 2^12.
    profiles.append(({}, Fraction(3, 2)))

    frontier, _ = divisor_frontier(split_primes, 2**72 - 1)
    frontier_indices = [size for size, _ in frontier]
    first_cutoff = bisect_right(inert_primes, 13)
    worst = (Fraction(0), None)
    count = 0
    for k in range(30, 72):
        lower = 2**k
        _, divisor_bound = frontier[bisect_right(frontier_indices, 2*lower - 1) - 1]
        for valuations, delta_bound in profiles:
            profile_mass = mass
            error_adjustment = Fraction(0)
            for q, valuation in valuations.items():
                if valuation:
                    profile_mass -= Fraction(2, q * (q - 1))
                    error_adjustment -= 1
                if valuation == 2:
                    profile_mass += Fraction(1, q*q)
                    error_adjustment += Fraction(1, q)
            error = delta_bound / cubert_floor(lower)
            for i in range(first_cutoff, len(inert_primes) + 1):
                upper = (
                    profile_mass + Fraction(24, 5) * divisor_bound * tails[i]
                    + 3 * error * (Fraction(8, 5) + i + error_adjustment)
                )
                if upper < Fraction(4, 5):
                    count += 1
                    if upper > worst[0]:
                        worst = upper, {
                            "k": k, "valuations": valuations,
                            "P": inert_primes[i - 1], "M": divisor_bound
                        }
                    break
            else:
                raise AssertionError(f"uncovered bin {k}, valuations {valuations}")
    assert count == 420

    # Beyond these bins take P=n^(1/4) and use the nine residue progressions.
    infinite_bound = (
        mass + Fraction(24, 5) * 7 * (Fraction(9, 38*64) + Fraction(9, 2**24))
        + 15 * (Fraction(9, 38*64) + Fraction(53, 5*2**24))
    )
    assert infinite_bound < Fraction(4, 5)
    return {
        "intervals": [30, 71], "profile_checks": count,
        "worst_interval": worst[1], "worst_upper_bound": str(worst[0]),
        "infinite_upper_bound_at_2^72": str(infinite_bound),
        "inert_delta_profiles": [
            {"valuations": valuations, "constant": str(bound)}
            for valuations, bound in profiles if valuations
        ],
    }



SHARP_PROFILES = {
    # p: (initial k, terminal k, odd lattice linear/square-root coefficients,
    #     tail bound for E*n**(1/3), bound for M/n**(1/6), largest candidate)
    43: (32, 72, Fraction(18, 25), Fraction(4), 7, 4, 74),
    67: (35, 84, Fraction(29, 50), Fraction(4), 9, 3, 138),
    163: (43, 84, Fraction(37, 100), Fraction(7, 2), 14, 2, 992),
}


def first_minimum(lo, hi, value):
    while lo < hi:
        mid = (lo + hi) // 2
        if value(mid) <= value(mid + 1):
            hi = mid
        else:
            lo = mid + 1
    return lo


def check_sharp_intervals(p, split_primes, inert_primes, mass, tails):
    start, stop, linear, root, delta_max, divisor_max, largest = SHARP_PROFILES[p]
    classes = (p - 1) // 2
    w = Fraction(6, p - 1)
    nu = 2*w / (p + 1)
    limit = 1 - Fraction(4, p + 1)
    assert (14*linear)**2 * p > 66**2  # linear > 3*pi/(2*sqrt(p))
    assert ((root - 3)/6)**2 * p > 1  # root > (3/2)*(2+4/sqrt(p))
    assert delta_envelope_sixth(inert_primes, {}) < delta_max**6
    assert prod(max(Fraction((a + 1)**6, q**a) for a in range(11))
                for q in split_primes if q < 64) < divisor_max**6

    boundary_start = 34 if p == 163 else start
    q_floor = 2**boundary_start // largest
    assert q_floor > largest
    odd_part_bound = Fraction(1)
    for m in range(1, largest//2 + 1, 2):
        divisors = [d for d in range(1, m + 1) if m % d == 0]
        odd_part_bound = max(
            odd_part_bound, Fraction(len(divisors)**6 * m**5, sum(divisors)**6)
        )
    boundary_sixth = Fraction(3, 2)**6 * odd_part_bound * Fraction(2**6, q_floor)
    boundary_constant = upward_sixth_root(boundary_sixth)
    assert boundary_constant <= delta_max

    harmonic = [0] * (len(inert_primes) + 1)
    for i in reversed(range(len(inert_primes))):
        harmonic[i] = harmonic[i + 1] + ceil_scaled(1, inert_primes[i])
    harmonic = [Fraction(factor, ROUNDING) for factor in harmonic]
    selected = inert_primes[:(6 if p == 163 else 4)]
    profiles = []
    for exponents in product(VALUATION_CLASSES, repeat=len(selected)):
        valuations = dict(zip(selected, exponents, strict=True))
        delta_sixth = delta_envelope_sixth(inert_primes, valuations)
        delta_constant = upward_sixth_root(delta_sixth)
        profile_mass, error_offset = mass, 2 - w - 2*nu
        shell_factor = Fraction(1)
        for q, exponent in valuations.items():
            shell_factor *= sum(Fraction(1, q**a) for a in range(exponent + 1))
            if exponent:
                profile_mass -= Fraction(2, q*(q-1))
                error_offset -= 1
            if exponent == 2:
                profile_mass += Fraction(1, q*q)
                error_offset += Fraction(1, q)
        profiles.append((valuations, delta_constant, profile_mass, error_offset, shell_factor))
    profiles.append(({}, boundary_constant, mass, 2-w-2*nu, Fraction(1)))

    frontier, _ = divisor_frontier(split_primes, 2**stop - 1)
    frontier_indices = [size for size, _ in frontier]
    worst = (Fraction(0), None)
    count = 0
    for k in range(min(start, boundary_start), stop):
        lower = 2**k
        _, M = frontier[bisect_right(frontier_indices, 2*lower - 1) - 1]
        # sqrt(n) < 2**ceil((k+1)/2). Beyond R use residue progressions and log(2)<7/10.
        harmonic_tail = classes * (
            Fraction(1, PRIME_LIMIT) + Fraction(7*max(0, (k+2)//2 - 16), 20*p)
        )
        square_root_lower = isqrt(lower)
        concentration_tail = [
            Fraction(2*M, 3*w) * (
                linear * t + (root*(h + harmonic_tail) + Fraction(1, 2))/square_root_lower
                + Fraction(2*classes, lower)
            )
            for t, h in zip(tails, harmonic, strict=True)
        ]
        for valuations, delta, profile_mass, offset, shell_factor in profiles:
            if valuations and k < start:
                continue
            E = delta / cubert_floor(lower)

            def margin(i):
                return profile_mass + concentration_tail[i]/shell_factor + (offset + i)*E/w

            i = first_minimum(len(selected), len(inert_primes), margin)
            upper = margin(i)
            count += 1
            if upper > worst[0]:
                worst = upper, dict(k=k, valuations=valuations, P=inert_primes[i-1], M=M)
            if upper >= limit:
                raise AssertionError((p, k, valuations, float(upper), float(limit)))

    # Above the last bin take P=n^(1/4)/2. With x=n^(1/12),
    # log(sqrt(n)/P)=log(2*x^3) <= 1+3*x.
    assert stop % 12 == 0
    x = 2**(stop // 12)
    rate = Fraction(classes, 2*p)
    infinite = (
        mass + Fraction(2*divisor_max, 3*w) * (
            linear*(2*rate/x + Fraction(4*classes, x**4))
            + root*(Fraction(2*classes, x**7) + rate*(1+3*x)/x**4)
            + 2*rate/x**4 + Fraction(2*classes, x**10)
        ) + Fraction(delta_max, w) * (rate/(2*x) + Fraction(classes+2, x**4))
    )
    assert infinite < limit
    return dict(inert_intervals=[start, stop-1], boundary_intervals=[boundary_start, stop-1],
                profile_checks=count, worst_interval=worst[1],
                worst_upper_bound=str(worst[0]), infinite_upper_bound=str(infinite))


def main():
    flags = prime_flags(PRIME_LIMIT)
    primes = [p for p in range(2, PRIME_LIMIT + 1) if flags[p]]
    assert 12**6 < 3 * 11**6
    assert 9**6 > 2 * 8**6 and 10**6 < 2 * 9**6
    assert delta_envelope_sixth(primes[1:], {}) < 16**6
    mass_upper_bounds = {}
    sharper = {}
    for p in (19, 43, 67, 163):
        split_primes = [q for q in primes if pow(q, (p - 1) // 2, p) == 1]
        inert_primes = [q for q in primes[1:] if pow(q, (p - 1) // 2, p) == p - 1]
        divisor_bound_sixth = prod(
            max(Fraction((a + 1)**6, q**a) for a in range(11))
            for q in split_primes if q < 64
        )
        assert divisor_bound_sixth < 7**6
        mass, tails = sieve_envelopes(p, inert_primes)
        assert mass < Fraction(27, 50)
        mass_upper_bounds[p] = str(mass)
        if p == 19:
            assert mass < Fraction(359, 1000)
            exact_19 = check_19_intervals(split_primes, inert_primes, mass, tails)

        else:
            sharper[p] = check_sharp_intervals(p, split_primes, inert_primes, mass, tails)

    common_upper_bound = (
        Fraction(27, 50) + Fraction(216, 5) * 7 * (Fraction(1, 4*2**10) + Fraction(81, 2**40))
        + 27 * 16 * (Fraction(1, 4*2**10) + Fraction(83, 2**40))
    )
    assert common_upper_bound < Fraction(4, 5)
    paths = ("verify_heegner_core_bounds.py", "verify_primitive_norm_minus11.py")
    print(json.dumps({
        "checker_inputs_sha256": {
            name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in paths
        },
        "common_onset_power_of_two": 120,
        "common_upper_bound": str(common_upper_bound),
        "inert_mass_upper_bounds": mass_upper_bounds,
        "discriminant_minus_19": exact_19,
        "sharp_reduced_family_certificates": sharper,
    }, indent=2))


if __name__ == "__main__":
    main()
