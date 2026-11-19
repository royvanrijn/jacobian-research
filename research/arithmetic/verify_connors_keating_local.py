#!/usr/bin/env python3
"""Exact rational checks of the Connors-Keating local pair densities.

Finite residue tests include exact corrections for unresolved valuation tails;
these are finite arithmetic checks, not a proof of the global asymptotic.
"""
from fractions import Fraction as Q
import json

PRIMES = (3, 7, 11, 19)
SHIFTS = range(1, 65)

def valuation(n, p):
    assert n > 0
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

def norm_flags(p, k):
    modulus = p ** k
    flags = bytearray(modulus)
    flags[0] = 1
    for residue in range(1, modulus):
        unit = residue
        v = 0
        while unit % p == 0:
            unit //= p
            v += 1
        if p == 2:
            flags[residue] = (k - v <= 1 or unit % 4 == 1)
        else:
            flags[residue] = (v % 2 == 0)
    return flags

def predicted_odd_prime(p, h):
    j = valuation(h, p)
    return Q(p, p + 1) * (1 - Q(1, p ** (j + 1)))

def square_valuation_partition(p, h):
    """Sum disjoint exact-even-valuation cells, with geometric tails."""
    j = valuation(h, p)
    m = j // 2
    if j % 2:
        return sum((Q(1, p ** (2 * a)) * (1 - Q(1, p))
                    for a in range(m + 1)), Q(0))
    equal_below = sum((Q(1, p ** (2 * a)) * (1 - Q(1, p))
                       for a in range(m)), Q(0))
    equal_at = Q(1, p ** (2 * m)) * (1 - Q(2, p))
    unequal_tail = (2 * (1 - Q(1, p)) * Q(1, p ** (2 * (m + 1)))
                    / (1 - Q(1, p * p)))
    return equal_below + equal_at + unequal_tail

def predicted_two(h):
    j = valuation(h, 2)
    w = Q(1) if j == 0 else 2 - Q(3, 2 ** j)
    return w / 4

def dyadic_partition(h):
    j = valuation(h, 2)
    if j == 0:
        return Q(1, 4)
    result = Q(1, 8)
    for _ in range(2, j + 1):
        result = Q(1, 4) + result / 2
    return result

def check():
    residue_checks = valuation_checks = singleton_checks = 0
    moduli = {}
    for p in PRIMES:
        max_j = max(valuation(h, p) for h in SHIFTS)
        depths = (max_j + 2, max_j + 3)
        moduli[str(p)] = [p ** k for k in depths]
        for k in depths:
            modulus = p ** k
            flags = norm_flags(p, k)
            # Conditional probability that a multiple of p^k has even valuation.
            unresolved_acceptance = Q(p if k % 2 == 0 else 1, p + 1)
            singleton = Q(sum(flags), modulus) - (1 - unresolved_acceptance) / modulus
            assert singleton == Q(p, p + 1), (p, k, singleton)
            singleton_checks += 1
            for h in SHIFTS:
                finite = Q(sum(flags[a] and flags[(a + h) % modulus]
                               for a in range(modulus)), modulus)
                j = valuation(h, p)
                assert k > j
                # Only residues 0 and -h remain unresolved; both contribute
                # only when the other, fixed valuation j is even.
                correction = (2 * (1 - unresolved_acceptance) / modulus
                              if j % 2 == 0 else Q(0))
                exact = finite - correction
                expected = predicted_odd_prime(p, h)
                assert exact == expected, (p, k, h, finite, correction, expected)
                residue_checks += 1
        for h in SHIFTS:
            assert square_valuation_partition(p, h) == predicted_odd_prime(p, h), (p, h)
            valuation_checks += 1
    max_j = max(valuation(h, 2) for h in SHIFTS)
    depths = (max_j + 3, max_j + 4)
    moduli['2'] = [2 ** k for k in depths]
    for k in depths:
        modulus = 2 ** k
        flags = norm_flags(2, k)
        assert Q(sum(flags), modulus) - Q(1, modulus) == Q(1, 2)
        singleton_checks += 1
        for h in SHIFTS:
            finite = Q(sum(flags[a] and flags[(a + h) % modulus]
                           for a in range(modulus)), modulus)
            assert k >= valuation(h, 2) + 3
            # The neighborhoods of 0 and -h have opposite fixed unit signs;
            # exactly one contributes an excess 1/modulus in the finite model.
            exact = finite - Q(1, modulus)
            assert exact == predicted_two(h), (2, k, h, finite, predicted_two(h))
            residue_checks += 1
    for h in SHIFTS:
        assert dyadic_partition(h) == predicted_two(h), h
        valuation_checks += 1
    result = dict(status='PASS', shifts='1..64', odd_primes=list(PRIMES),
                  residue_pair_checks=residue_checks,
                  exact_valuation_partition_checks=valuation_checks,
                  singleton_checks=singleton_checks, moduli=moduli,
                  arithmetic='fractions.Fraction throughout',
                  global_asymptotic_tested=False)
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    check()
