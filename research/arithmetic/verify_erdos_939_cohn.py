#!/usr/bin/env python3
"""Check full triples and possible Cohn seeds in pinned abc and S-unit tables."""

import argparse
import gzip
import hashlib
import json
from math import gcd, isqrt

from verify_erdos_939 import full_values, least_prime_factors


ABC_ARCHIVES = {
    "abc-small": (
        "8d40b46d6962aa14b2f4723eafcbaa68644a786a5b35db0d706a22ce51fd912c",
        14_482_065,
        1,
        10**18 - 1,
    ),
    "abc-additional": (
        "a7e6b9c755a156f44407fb944d425faf688e9a1398ca66015f8ca6afd0979c2b",
        9_345_651,
        10**18,
        2**63 - 1,
    ),
}
S_UNIT_TABLES = {
    "primes-53": (
        "80f3cd4e2535285f9a5c5811d9646b60523f7caef6fb0506202af761a743da98",
        1_663_449,
    ),
    "radical-1e7": (
        "c28a87a2015c987d6665647117094cab4e9072a7d23301c5d89f892d3d2b2fb3",
        432_408,
    ),
}


def scan_abc(path, minimum, maximum):
    full = full_values(maximum)
    count = degenerate = 0
    hits = []
    with gzip.open(path, "rb") as rows:
        for row in rows:
            c, a = map(int, row.split())
            b = c - a
            if not (minimum <= c <= maximum and a > 0 and b > 0 and gcd(a, b) == 1):
                raise ValueError(f"invalid triple {(a, b, c)}")
            count += 1
            if a == b:
                degenerate += 1
                continue
            if c % 2 == 0:
                even, odd1, odd2 = c, a, b
            elif a % 2 == 0:
                even, odd1, odd2 = a, b, c
            else:
                even, odd1, odd2 = b, a, c
            if odd1 in full and odd2 in full and even // (even & -even) in full:
                hits.append((a, b, c))
    return {"rows": count, "degenerate_seeds": degenerate, "possible_seeds": hits}


def scan_s_units(path, support):
    cap = 10**7
    least_prime = least_prime_factors(isqrt(cap))
    primes = [p for p in range(2, len(least_prime)) if p == least_prime[p]]
    first_sixteen = primes[:16]
    count = degenerate = height = 0
    full = []
    possible_seeds = []
    with open(path, "rt", encoding="utf-8") as rows:
        for row in rows:
            if row.startswith("#") or not row.strip():
                continue
            rad, a, plus, b, equal, c, _ = row.split()
            r, a, b, c = int(rad.removesuffix(":")), int(a), int(b), int(c)
            if not (
                plus == "+" and equal == "=" and 0 < a <= b < c == a + b
                and gcd(a, b) == 1 and r > 1
            ):
                raise ValueError(f"invalid triple {row}")
            if support == "primes-53":
                available = first_sixteen
            else:
                if r > cap:
                    raise ValueError(f"radical out of range: {row}")
                available = []
                remaining = r
                for prime in primes:
                    if prime * prime > remaining:
                        break
                    if remaining % prime == 0:
                        remaining //= prime
                        available.append(prime)
                        if remaining % prime == 0:
                            raise ValueError(f"nonsquarefree radical: {row}")
                if remaining > 1:
                    available.append(remaining)
            tracked_radical = 1
            defect_at_two = False
            odd_defect = False
            for value in (a, b, c):
                for prime in available:
                    if value == 1:
                        break
                    if value % prime:
                        continue
                    tracked_radical *= prime
                    exponent = 0
                    while value % prime == 0:
                        value //= prime
                        exponent += 1
                    if exponent < 4:
                        if prime == 2:
                            defect_at_two = True
                        else:
                            odd_defect = True
                if value != 1:
                    raise ValueError(f"incomplete support: {row}")
            if tracked_radical != r:
                raise ValueError(f"incorrect radical: {row}")
            count += 1
            height = max(height, c)
            degenerate += a == b
            if not odd_defect:
                if a != b:
                    possible_seeds.append((a, b, c))
                if not defect_at_two:
                    full.append((a, b, c))
    return {
        "rows": count,
        "maximum_height": height,
        "degenerate_seeds": degenerate,
        "four_full_triples": full,
        "possible_seeds": possible_seeds,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", choices=(*ABC_ARCHIVES, *S_UNIT_TABLES))
    parser.add_argument("path")
    args = parser.parse_args()
    if args.table in ABC_ARCHIVES:
        digest, expected_rows, minimum, maximum = ABC_ARCHIVES[args.table]
    else:
        digest, expected_rows = S_UNIT_TABLES[args.table]
    with open(args.path, "rb") as source:
        if hashlib.file_digest(source, "sha256").hexdigest() != digest:
            raise ValueError(f"unexpected SHA-256 for {args.path}")
    if args.table in ABC_ARCHIVES:
        result = scan_abc(args.path, minimum, maximum)
    else:
        result = scan_s_units(args.path, args.table)
    if result["rows"] != expected_rows:
        raise ValueError(f"expected {expected_rows} rows, read {result['rows']}")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
