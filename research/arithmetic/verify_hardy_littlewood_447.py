#!/usr/bin/env python3
"""Verify the Engelsma 447-tuple and exclude a bounded range of translates."""

from math import isqrt
from pathlib import Path


def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def extend_residues(residues, modulus, allowed, prime):
    inverse = pow(modulus, -1, prime)
    return [
        r + modulus * (((a - r) * inverse) % prime)
        for r in residues
        for a in allowed
    ]


def verify(offsets):
    if len(offsets) != 447 or sorted(set(offsets)) != offsets:
        raise ValueError("expected 447 distinct increasing offsets")
    if (offsets[0], offsets[-1]) != (0, 3158):
        raise ValueError("unexpected endpoints")
    width = offsets[-1] + 1
    if sum(is_prime(n) for n in range(2, width + 1)) != 446:
        raise ValueError("unexpected initial prime count")

    allowed = {
        p: set(range(p)) - {(-h) % p for h in offsets}
        for p in range(2, len(offsets) + 1)
        if is_prime(p)
    }
    if not all(allowed.values()):
        raise ValueError("tuple is inadmissible")

    # Handle starts for which divisibility by a tested prime can be equality.
    if any(all(is_prime(n + h) for h in offsets) for n in range(1, max(allowed) + 1)):
        raise ValueError("small prime translate found")

    anchor, modulus = 0, 1
    for p, classes in allowed.items():
        if len(classes) == 1:
            (anchor,) = extend_residues([anchor], modulus, classes, p)
            modulus *= p

    # Every remaining start is anchor + modulus * q, with q >= 0.
    # Test one full period of the first eight nonsingleton constraints.
    quotient_residues, period = [0], 1
    remaining = [(p, classes) for p, classes in allowed.items() if len(classes) > 1]
    for p, classes in remaining[:8]:
        inverse = pow(modulus, -1, p)
        quotient_classes = {(inverse * (a - anchor)) % p for a in classes}
        quotient_residues = extend_residues(quotient_residues, period, quotient_classes, p)
        period *= p

    for q in quotient_residues:
        n = anchor + modulus * q
        if all(n % p in classes for p, classes in remaining[8:]):
            raise ValueError(f"locally surviving start {n}")

    print(f"447 admissible offsets in width {width}, with pi({width}) = 446")
    print(f"Forced start: {anchor} modulo {modulus}")
    print(f"Tested {len(quotient_residues)} quotient residues in [0, {period})")
    print(f"No prime translate with 1 <= start < {anchor + modulus * period}")


if __name__ == "__main__":
    path = Path(__file__).with_name("admissible_447_3158.txt")
    verify(list(map(int, path.read_text(encoding="ascii").split())))
