#!/usr/bin/env python3
"""Replay catalog and fixed-coefficient searches for Erdős problem 939."""

import argparse
import gzip
import hashlib
import json
from math import gcd, isqrt


def fourth_root(n):
    return isqrt(isqrt(n))


def least_prime_factors(limit):
    factors = list(range(limit + 1))
    for prime in range(2, isqrt(limit) + 1):
        if factors[prime] == prime:
            for n in range(prime * prime, limit + 1, prime):
                if factors[n] == n:
                    factors[n] = prime
    return factors


def factor(n, least_prime):
    while n > 1:
        prime = least_prime[n]
        exponent = 0
        while n % prime == 0:
            n //= prime
            exponent += 1
        yield prime, exponent


def coefficients(limit):
    """Yield the canonical multipliers d5^5 d6^6 d7^7 and their radicals."""
    # The smallest contribution of a prime to a nontrivial multiplier is p^5.
    low, high = 0, 1 << ((limit.bit_length() + 4) // 5)
    while high - low > 1:
        middle = (low + high) // 2
        if middle**5 <= limit:
            low = middle
        else:
            high = middle
    factors = least_prime_factors(low)
    primes = [p for p in range(2, low + 1) if factors[p] == p]

    def extend(value, radical, start):
        yield value, radical
        for index in range(start, len(primes)):
            prime = primes[index]
            if value * prime**5 > limit:
                break
            for exponent in (5, 6, 7):
                next_value = value * prime**exponent
                if next_value <= limit:
                    yield from extend(next_value, radical * prime, index + 1)

    yield from extend(1, 1, 0)


def full_values(limit):
    return {
        coefficient * base**4
        for coefficient, _ in coefficients(limit)
        for base in range(1, fourth_root(limit // coefficient) + 1)
    }


def enumeration_digest(bound):
    summands = hashlib.sha256()
    targets = hashlib.sha256()
    count = 0
    values = sorted(full_values(2 * bound))
    for value in values:
        word = value.to_bytes(16, "little")
        targets.update(word)
        if value <= bound:
            summands.update(word)
            count += 1
    return {
        "summands": count,
        "targets": len(values),
        "summands_sha256": summands.hexdigest(),
        "targets_sha256": targets.hexdigest(),
    }


def scan_catalog(path, limit, expected_rows):
    full = full_values(limit)
    counts = [0, 0, 0]
    hits = []
    with gzip.open(path, "rb") as rows:
        for row in rows:
            c, a = map(int, row.split())
            if not 0 < a < c <= limit:
                raise ValueError(f"out-of-range row {(c, a)}")
            counts[0] += 1
            if c in full:
                counts[1] += 1
                b = c - a
                if a in full and b in full and gcd(a, b) == 1:
                    hits.append((a, b, c))
                if a in full or b in full:
                    counts[2] += 1
    if counts[0] != expected_rows:
        raise ValueError(f"expected {expected_rows} rows, read {counts[0]}")
    return {
        "height_limit": limit,
        "four_full_values": len(full),
        "rows": counts[0],
        "four_full_targets": counts[1],
        "four_full_targets_with_a_four_full_summand": counts[2],
        "hits": hits,
    }


def locally_allowed(a, b, c, prime):
    modulus = 16 if prime == 2 else prime
    units = {pow(x, 4, modulus) for x in range(1, modulus) if x % prime}
    signed = (a, b, -c)
    for i in range(3):
        # The omitted coordinate may vanish modulo p. The others are units.
        j, k = (i + 1) % 3, (i + 2) % 3
        if signed[j] % prime and signed[k] % prime:
            ratio = -signed[j] * pow(signed[k], -1, modulus) % modulus
            if ratio in units:
                return True
    if any(coefficient % prime == 0 for coefficient in signed):
        return False
    # Divide through by the fourth power of the first coordinate.
    inverse_c = pow(c, -1, modulus)
    return any((a + b * value) * inverse_c % modulus in units for value in units)


def ratio_roots(a, b, c, z, least_prime, local_roots):
    exponents = dict(factor(c, least_prime))
    for prime, exponent in factor(z, least_prime):
        exponents[prime] = exponents.get(prime, 0) + 4 * exponent
    roots, modulus = [0], 1
    for prime, exponent in exponents.items():
        prime_power = prime**exponent
        local = local_roots(a, b, prime, exponent)
        if not local:
            return [], modulus
        inverse = pow(modulus, -1, prime_power)
        roots = [
            root + modulus * ((other - root) * inverse % prime_power)
            for root in roots
            for other in local
        ]
        modulus *= prime_power
    return roots, modulus


def reconstruct(root, modulus, x_bound, y_bound):
    """Recover a reduced positive x/y if 2*x_bound*y_bound < modulus."""
    previous, remainder = modulus, root
    previous_y, y = 0, 1
    while remainder and abs(y) <= y_bound:
        if remainder <= x_bound and y > 0 and gcd(remainder, y) == 1:
            return remainder, y
        quotient, next_remainder = divmod(previous, remainder)
        previous, remainder = remainder, next_remainder
        previous_y, y = y, previous_y - quotient * y
    return None


def search_coefficients(limit, base_limit):
    from functools import cache

    from sympy.ntheory.residue_ntheory import sqrt_mod

    @cache
    def local_roots(a, b, prime, exponent):
        if a % prime == 0 or b % prime == 0:
            return ()
        modulus = prime**exponent
        ratio = -b * pow(a, -1, modulus) % modulus
        return tuple(
            x
            for square in sqrt_mod(ratio, modulus, all_roots=True)
            for x in sqrt_mod(square, modulus, all_roots=True)
        )

    bank = sorted(coefficients(limit))
    least_prime = least_prime_factors(max(limit, base_limit))
    screened = locally_surviving = targets = root_targets = 0
    hits = []
    for index, (a, rad_a) in enumerate(bank):
        for b, rad_b in bank[index:]:
            if gcd(rad_a, rad_b) != 1:
                continue
            for c, rad_c in bank:
                if gcd(rad_a * rad_b, rad_c) != 1:
                    continue
                screened += 1
                primes = {2, 3, 5, 7, 11, 13, 17, 29, 41}
                for radical in (rad_a, rad_b, rad_c):
                    primes.update(prime for prime, _ in factor(radical, least_prime))
                if not all(locally_allowed(a, b, c, p) for p in sorted(primes)):
                    continue
                locally_surviving += 1
                for z in range(rad_c, base_limit + 1, rad_c):
                    modulus = c * z**4
                    x_bound = fourth_root((modulus - 1) // a)
                    y_bound = fourth_root((modulus - 1) // b)
                    if not x_bound or not y_bound:
                        continue
                    if gcd(z, rad_a * rad_b) != 1:
                        continue
                    targets += 1
                    roots, recovered_modulus = ratio_roots(a, b, c, z, least_prime, local_roots)
                    if roots:
                        if recovered_modulus != modulus or 2 * x_bound * y_bound >= modulus:
                            raise ArithmeticError("invalid reconstruction bounds")
                        root_targets += 1
                    for root in roots:
                        pair = reconstruct(root, modulus, x_bound, y_bound)
                        if pair is None:
                            continue
                        x, y = pair
                        first, second = a * x**4, b * y**4
                        if first + second == modulus and gcd(first, second) == 1:
                            hits.append((min(first, second), max(first, second), modulus))
    return {
        "coefficient_limit": limit,
        "coefficient_count": len(bank),
        "target_base_limit": base_limit,
        "coprime_coefficient_triples": screened,
        "locally_surviving_coefficient_triples": locally_surviving,
        "target_values_considered": targets,
        "target_values_with_modular_roots": root_targets,
        "hits": sorted(set(hits)),
    }


def positive_int(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("expected a positive integer")
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    catalog = commands.add_parser("catalog")
    catalog.add_argument("path")
    catalog.add_argument("limit", type=positive_int)
    catalog.add_argument("expected_rows", type=positive_int)
    enumeration = commands.add_parser("enumeration")
    enumeration.add_argument("max_summand", type=positive_int)
    fixed = commands.add_parser("fixed")
    fixed.add_argument("coefficient_limit", type=positive_int)
    fixed.add_argument("target_base_limit", type=positive_int)
    args = parser.parse_args()
    if args.command == "catalog":
        result = scan_catalog(args.path, args.limit, args.expected_rows)
    elif args.command == "enumeration":
        result = enumeration_digest(args.max_summand)
    else:
        result = search_coefficients(args.coefficient_limit, args.target_base_limit)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
