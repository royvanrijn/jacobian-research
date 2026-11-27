"""Exact 2-adic factorial-sum excess for the stated finite set of starts."""

from collections.abc import Iterable


def v2(n: int) -> int:
    return (n & -n).bit_length() - 1


def factorial_valuation(n: int) -> int:
    return n - n.bit_count()


def excess_at_failure(start: int, cutoff: int) -> tuple[int, int] | None:
    """Return the first failure index and exact excess, or None if cutoff lifts."""
    start_valuation = factorial_valuation(start)
    cutoff_valuation = factorial_valuation(cutoff)
    precision = cutoff_valuation - start_valuation
    mask = (1 << precision) - 1
    unit = 1
    # Divide all sums by the odd part of start! and by the current power of 2.
    states = {1}

    # The two factorials in a block have the same valuation.
    for n in range(start, cutoff, 2):
        jump = v2(n + 2)
        odd_weight = ((n + 1) * unit) & mask
        options = (0, odd_weight) if n == start else (0, unit, odd_weight, unit + odd_weight)
        next_states: set[int] = set()
        low_mask = (1 << jump) - 1
        best_low_valuation = 0
        for value in states:
            for offset in options:
                candidate = value + offset
                low = candidate & low_mask
                if low == 0:
                    next_states.add((candidate & mask) >> jump)
                else:
                    best_low_valuation = max(best_low_valuation, v2(low))
        if not next_states:
            excess = factorial_valuation(n) - start_valuation + best_low_valuation
            return n + 2, excess

        states = next_states
        mask >>= jump
        unit = (odd_weight * ((n + 2) >> jump)) & mask
    return None


def checked_profile(starts: Iterable[int]) -> dict[int, tuple[int, int]]:
    profile = {}
    for start in starts:
        cutoff = 18_448 if start == 34 else 4_096
        answer = excess_at_failure(start, cutoff)
        if answer is None:
            raise RuntimeError(f"Start {start} survives the specified cutoff")
        profile[start] = answer
    return profile


if __name__ == "__main__":
    profile = checked_profile(range(2, 81, 2))
    for start, (boundary, excess) in profile.items():
        print(f"{start} {boundary} {excess}")
    assert profile[34] == (18_448, 18_412)
    assert max(excess for start, (_, excess) in profile.items() if start != 34) < 18_412

