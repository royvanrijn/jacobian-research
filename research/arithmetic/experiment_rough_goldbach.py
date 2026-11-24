#!/usr/bin/env python3
"""Exact finite checks for the rough Goldbach and discriminant calculations."""

from fractions import Fraction
from math import gcd, isqrt


def least_prime_factors(limit):
    least = list(range(limit + 1))
    for prime in range(2, isqrt(limit) + 1):
        if least[prime] == prime:
            for multiple in range(prime * prime, limit + 1, prime):
                if least[multiple] == multiple:
                    least[multiple] = prime
    return least


def squarefree_data(least):
    omega = bytearray(len(least))
    squarefree = bytearray(len(least))
    squarefree[1] = 1
    for number in range(2, len(least)):
        prime = least[number]
        cofactor = number // prime
        omega[number] = omega[cofactor] + 1
        squarefree[number] = squarefree[cofactor] and cofactor % prime != 0
    return omega, squarefree


def integer_root(number, degree):
    low, high = 0, 1 << ((number.bit_length() + degree - 1) // degree)
    while high - low > 1:
        midpoint = (low + high) // 2
        if midpoint**degree <= number:
            low = midpoint
        else:
            high = midpoint
    return low


def coefficients(target, least, omega, squarefree, exponent):
    numerator, denominator = exponent
    # The note uses strict roughness. All reported counts are ordered.
    cutoff = integer_root(target**numerator, denominator)
    max_factors = (denominator - 1) // numerator
    counts = [[0] * (max_factors + 1) for _ in range(max_factors + 1)]
    for left in range(cutoff + 1 + cutoff % 2, target // 2, 2):
        right = target - left
        if (
            least[left] > cutoff
            and least[right] > cutoff
            and squarefree[left]
            and squarefree[right]
            and gcd(left, right) == 1
        ):
            i, j = omega[left], omega[right]
            counts[i][j] += 1
            counts[j][i] += 1
    return counts


def by_degree(counts):
    return [
        sum(counts[i][degree - i] for i in range(1, len(counts))
            if 1 <= degree - i < len(counts))
        for degree in range(2, 2 * len(counts) - 1)
    ]


def trial_prime_factors(number):
    factors = []
    for divisor in range(2, isqrt(number) + 1):
        while number % divisor == 0:
            factors.append(divisor)
            number //= divisor
    if number > 1:
        factors.append(number)
    return factors


def crosscheck_small_target(target, counts, exponent):
    numerator, denominator = exponent
    reference = [0] * len(by_degree(counts))
    for left in range(2, target - 1):
        factors = trial_prime_factors(left) + trial_prime_factors(target - left)
        if len(factors) == len(set(factors)) and min(factors)**denominator > target**numerator:
            reference[len(factors) - 2] += 1
    assert reference == by_degree(counts)


def parity_counts(counts):
    prime_prime, prime_semiprime = counts[1][1:]
    semiprime_prime, semiprime_semiprime = counts[2][1:]
    assert prime_semiprime == semiprime_prime
    total = prime_prime + 2 * prime_semiprime + semiprime_semiprime
    marginal = semiprime_semiprime - prime_prime
    correlation = prime_prime - 2 * prime_semiprime + semiprime_semiprime
    assert total - 2 * marginal + correlation == 4 * prime_prime
    return total, marginal, correlation


def prime_support(number, least):
    support = []
    while number > 1:
        prime = least[number]
        support.append(prime)
        number //= prime
    return support


def discriminant_check(target, least, squarefree):
    cutoff = integer_root(target, 5)
    supports = []
    residues = {}
    for left in range(cutoff + 1 + cutoff % 2, target // 2, 2):
        right = target - left
        if (
            min(least[left], least[right]) > cutoff
            and squarefree[left] and squarefree[right] and gcd(left, right) == 1
        ):
            support = prime_support(left, least) + prime_support(right, least)
            local = Fraction(1)
            for prime in support:
                if prime not in residues:
                    residues[prime] = {-i * i % prime for i in range(1, prime)}
                local *= Fraction(prime, 2 * (prime + 1))
            supports.append((support, local))
    assert supports
    for scale in (1024, 4096, 16384, 65536):
        discriminants = [d for d in range(scale + 3, 2 * scale, 4) if squarefree[d]]
        probabilities = [
            Fraction(sum(all(d % q in residues[q] for q in support)
                         for d in discriminants), len(discriminants))
            for support, _ in supports
        ]
        errors = []
        for fields in (1, 2, 4):
            corrected = sum((known / (local * 2**len(support)))**fields
                            for (support, local), known in zip(supports, probabilities))
            independent = sum(Fraction(1, 2**(fields * len(support)))
                              for support, _ in supports)
            errors.append(f"{float(corrected / independent - 1):+.3%}")
        print("discriminants", scale, len(discriminants), *errors)


def main():
    least = least_prime_factors(1_000_000)
    omega, squarefree = squarefree_data(least)
    cases = []
    for target in range(100, 5001, 2):
        counts = coefficients(target, least, omega, squarefree, (1, 3))
        total, marginal, correlation = parity_counts(counts)
        if total - 2 * sum(counts[2]) <= 0:
            cases.append((target, by_degree(counts), marginal, correlation))
    print("one-mark certificate nonpositive", len(cases), "of", 2451, "first", cases[:3])
    exponents = [(1, 5), (9, 40), (1, 4), (1, 3)]
    for target in (1000, 10000, 100000, 999998, 1000000):
        print("target", target)
        for exponent in exponents:
            counts = coefficients(target, least, omega, squarefree, exponent)
            degree_counts = by_degree(counts)
            if target == 1000:
                crosscheck_small_target(target, counts, exponent)
            print("  roughness", exponent, "by degree", degree_counts)
            if target == 1_000_000 and exponent == (1, 5):
                total = sum(degree_counts)
                first_moment = sum(degree * count
                                   for degree, count in enumerate(degree_counts, 2))
                print("  factor-count certificate", 3 * total - first_moment)
                for t in (Fraction(1, 4), Fraction(1, 2)):
                    extraction = sum(count * t**i for i, count in enumerate(degree_counts))
                    print("  extraction at", t, "lower bound", extraction - t * total)
            if exponent == (1, 3):
                total, marginal, correlation = parity_counts(counts)
                print("  parity cells", counts[1][1:], counts[2][1:],
                      "marginal and correlation over R",
                      f"{marginal / total:+.4f}", f"{correlation / total:+.4f}")
    discriminant_check(1000, least, squarefree)


if __name__ == "__main__":
    main()
