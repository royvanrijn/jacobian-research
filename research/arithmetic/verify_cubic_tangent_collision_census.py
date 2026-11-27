"""With Sage's Python, count exact tangent-field collisions up to a center bound."""

import argparse
from math import isqrt

from sage.all import pari


def census(limit):
    prime = bytearray(b"\1") * (2 * limit + 1)
    prime[:2] = b"\0\0"
    for q in range(2, isqrt(2 * limit) + 1):
        if prime[q]:
            prime[q * q : 2 * limit + 1 : q] = b"\0" * ((2 * limit - q * q) // q + 1)

    first = {}
    total = 0
    repeated = 0
    for p in range(7, limit + 1, 2):
        if not prime[p]:
            continue
        for h in range(6, p - 3, 6):
            if not prime[p - h] or not prime[p + h]:
                continue
            radicand = h**4 + 6 * p**2 * h**2 - 3 * p**4
            d = int(pari(radicand).core())
            total += 1
            if d in first:
                repeated += 1
                print(f"Field {d} occurs at both {first[d]} and {(p, h)}")
            else:
                first[d] = (p, h)
    return total, repeated


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("limit", nargs="?", type=int, default=45000)
    args = parser.parse_args()
    if args.limit < 1:
        parser.error("the center bound must be positive")
    total, repeated = census(args.limit)
    if args.limit == 45000:
        assert (total, repeated) == (1365099, 0)
    print(f"Centers <= {args.limit}: {total} triangles, {repeated} repeat events")
