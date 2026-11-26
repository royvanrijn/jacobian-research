#!/usr/bin/env python3
"""Verify uniform local covering of factorial subsums at jumps through size four."""

from math import prod

BASE = 5


def valuation_of_factorial(n: int) -> int:
    valuation = 0
    while n:
        n //= BASE
        valuation += n
    return valuation


def subsums_at_boundary(n: int, width: int, modulus: int) -> int:
    """Return a bit set of residues after dividing by the initial factorial."""
    mask = (1 << modulus) - 1
    reachable = weight = 1
    start = n - width
    for index in range(start, n):
        if index != start:
            weight = weight * index % modulus
        if weight:
            rotated = ((reachable << weight) & mask) | (reachable >> (modulus - weight))
            reachable |= rotated
    return reachable


def covers_split_block(n: int, start: int, split: int) -> bool:
    """Check each low residue with a separate suffix bit set."""
    low_modulus = BASE ** (valuation_of_factorial(split) - valuation_of_factorial(start))
    high_modulus = BASE ** (valuation_of_factorial(n) - valuation_of_factorial(split))
    full_modulus = low_modulus * high_modulus
    ratio = prod(range(start + 1, split + 1))
    inverse_unit = pow(ratio // low_modulus, -1, high_modulus)

    prefix_sums = [0]
    weight = 1
    for index in range(start, split):
        if index != start:
            weight = weight * index % full_modulus
        prefix_sums += [(value + weight) % full_modulus for value in prefix_sums]
    shifts = [set() for _ in range(low_modulus)]
    for value in prefix_sums:
        high, low = divmod(value, low_modulus)
        shifts[low].add(high * inverse_unit % high_modulus)

    mask = (1 << high_modulus) - 1
    missing = mask ^ subsums_at_boundary(n, n - split, high_modulus)
    for admissible_shifts in shifts:
        residual = mask
        for shift in admissible_shifts:
            rotated = ((missing << shift) & mask) | (missing >> (high_modulus - shift))
            residual &= rotated
            if not residual:
                break
        if residual:
            return False
    return True


def verify() -> None:
    for jump, width in ((1, 5), (2, 10), (3, 20), (4, 25)):
        precision = jump + valuation_of_factorial(width - 1)
        modulus = BASE**precision
        step = BASE**jump
        # At precision one the entire residue class is represented by zero.
        classes = [0] if precision == jump else [
            residue for residue in range(step, modulus, step)
            if residue % (BASE * step)
        ]
        mask = (1 << modulus) - 1
        for residue in classes:
            if subsums_at_boundary(residue, width, modulus) != mask:
                raise AssertionError((jump, width, residue))
        print(f"jump={jump}, width={width}, precision={precision}, classes={len(classes)}")

    # These are all size-five boundaries before the first jump seven.
    width = 35
    modulus = BASE ** (5 + valuation_of_factorial(width - 1))
    mask = (1 << modulus) - 1
    for boundary in range(BASE**5, BASE**7, BASE**5):
        if boundary % BASE**6 == 0:
            continue
        if subsums_at_boundary(boundary, width, modulus) != mask:
            raise AssertionError((5, width, boundary))
    print(f"jump=5, width={width}, precision=12, boundaries=20")

    # Four jump-six blocks use a split to avoid a 6.1 billion bit array.
    for boundary in range(BASE**6, BASE**7, BASE**6):
        if not covers_split_block(boundary, boundary - 40, boundary - 25):
            raise AssertionError((6, 40, boundary))
    if covers_split_block(BASE**6, BASE**6 - 35, BASE**6 - 25):
        raise AssertionError("short split block negative control")
    print("jump=6, width=40, precision=14, boundaries=4")

    # The shorter 15-term block at N=125 misses four residues.
    if subsums_at_boundary(125, 15, BASE**5).bit_count() != BASE**5 - 4:
        raise AssertionError("short block negative control")


if __name__ == "__main__":
    verify()
