"""Replay the necessary characteristic-sequence enumeration at degree twenty."""

from math import gcd


def semigroup_contains(generators, target):
    represented = [True] + [False] * target
    for degree in range(1, target + 1):
        represented[degree] = any(
            generator <= degree and represented[degree - generator]
            for generator in generators
        )
    return represented[target]


def characteristic_sequences(conductor):
    def extend(prefix, previous_gcd, current_gcd, consumed):
        if current_gcd == 1:
            if consumed == conductor:
                yield prefix
            return
        for next_gcd in range(1, current_gcd):
            if current_gcd % next_gcd:
                continue
            quotient = current_gcd // next_gcd
            maximum = next_gcd + (conductor - consumed) // (quotient - 1)
            for entry in range(next_gcd + 1, maximum + 1):
                if gcd(current_gcd, entry) != next_gcd:
                    continue
                if prefix[-1] * previous_gcd <= entry * current_gcd:
                    continue
                if not semigroup_contains(prefix, quotient * entry):
                    continue
                cost = (quotient - 1) * (entry - next_gcd)
                yield from extend(prefix + (entry,), current_gcd, next_gcd, consumed + cost)

    for divisor in range(1, conductor // 2 + 1):
        for lower in range(2, conductor // divisor + 2):
            upper_limit = conductor // (divisor * (lower - 1)) + 1
            for upper in range(lower + 1, upper_limit + 1):
                if gcd(lower, upper) != 1:
                    continue
                initial_cost = divisor * (upper - 1) * (lower - 1)
                first_degree, second_degree = divisor * upper, divisor * lower
                yield from extend(
                    (first_degree, second_degree), first_degree, divisor, initial_cost
                )


if __name__ == "__main__":
    assert set(characteristic_sequences(18)) == {
        (19, 2), (10, 3), (7, 4), (10, 4, 11), (14, 4, 7), (8, 6, 7),
        (10, 6, 3), (9, 6, 7), (15, 6, 4), (15, 10, 3), (12, 8, 6, 7),
    }
    expected = {
        (6, 5), (11, 3), (21, 2),
        (8, 6, 9), (10, 4, 13), (10, 6, 5), (14, 4, 9),
        (9, 6, 8), (15, 6, 5), (21, 6, 2), (12, 8, 5), (21, 14, 2),
        (12, 8, 6, 9), (12, 8, 10, 5),
    }
    assert set(characteristic_sequences(20)) == expected
    assert set(characteristic_sequences(22)) == {
        (23, 2), (8, 6, 11), (10, 4, 15), (14, 4, 11), (15, 10, 4), (12, 8, 6, 11),
    }
    print("Verified the necessary characteristic sequences at conductors 18, 20 and 22.")
