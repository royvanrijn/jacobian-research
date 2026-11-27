"""Check the block alphabets used for the dyadic cancellation dimension."""

from collections import Counter
from math import prod

from verify_erdos404_local_blocks import subsums_at_boundary


def alphabet_product() -> int:
    """Return the product of the alphabet sizes over one full period."""
    width = 16
    cap = 5
    counts: Counter[int] = Counter()
    for extra in range(cap + 1):
        modulus = 1 << (15 + extra)
        block_period = modulus // width
        # b = 16k. For extra < cap, v_2(k+1) is exactly extra.
        # At the cap all higher valuations use the same modulus.
        step = 1 << (extra if extra == cap else extra + 1)
        repeats = 1 << (cap - extra)
        for k in range((1 << extra) - 1, block_period, step):
            boundary = width * (k + 1)
            size = subsums_at_boundary(boundary, width, modulus).bit_count()
            counts[size] += repeats
    assert sum(counts.values()) == 1 << 16
    return prod(size**count for size, count in counts.items())


if __name__ == "__main__":
    log_floor = alphabet_product().bit_length() - 1
    assert 100_000 * log_floor > 91_598 * 16 * (1 << 16)
    print(f"floor(log2(product))={log_floor}; dimension > 0.91598")
