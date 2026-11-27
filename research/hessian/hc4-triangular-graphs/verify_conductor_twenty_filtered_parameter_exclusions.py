"""Replay the Laurent degree and numerical-semigroup comparisons."""


def largest_nongap(generators, limit):
    reachable = {0}
    for degree in range(limit + 1):
        if degree in reachable:
            reachable.update(degree + value for value in generators if degree + value <= limit)
    return max(reachable)


def check_square_root_case(n, m, r, expected):
    e = 2*n // m
    assert e*m == 2*n
    assert len({m*i % n for i in range(e)}) == e
    assert m-1+r-n < 0
    degrees = [
        (index, m*index+m-1-n)
        for index in range(e-2, -1, -1)
        if m*index+m-1-n >= 0
    ]
    comparisons = [(index, degree, largest_nongap((n, m, r), degree+1))
                   for index, degree in degrees]
    assert comparisons == expected
    assert all(largest-1 < degree for _, degree, largest in comparisons)


def check_cubic_case(n, m, r, expected):
    e = 3*n // m
    assert e*m == 3*n
    assert len({m*i % n for i in range(e)}) == e
    assert m-1+r < 2*n
    candidates = (
        [("a", index, m*index+m-1-n) for index in range(2*n//m+1)]
        + [("b", index, m*index+m-1-2*n) for index in range(e-1)]
    )
    degrees = sorted(
        (entry for entry in candidates if entry[2] >= 0), key=lambda entry: -entry[2]
    )
    assert len({degree for _, _, degree in degrees}) == len(degrees)
    comparisons = [(name, index, degree, largest_nongap((n, m, r), degree+1))
                   for name, index, degree in degrees]
    assert comparisons == expected
    assert all(largest-1 < degree for _, _, degree, largest in comparisons)


if __name__ == "__main__":
    check_square_root_case(14, 4, 9, [(5, 9, 9), (4, 5, 4), (3, 1, 0)])
    check_square_root_case(15, 6, 5, [(3, 8, 6), (2, 2, 0)])
    check_square_root_case(21, 6, 2, [(5, 14, 14), (4, 8, 8), (3, 2, 2)])
    check_cubic_case(8, 6, 9, [("a", 2, 9, 9), ("a", 1, 3, 0), ("b", 2, 1, 0)])
    check_cubic_case(10, 6, 5, [("a", 3, 13, 12), ("a", 2, 7, 6), ("b", 3, 3, 0), ("a", 1, 1, 0)])
    print("Verified the five conductor-twenty Laurent and semigroup comparisons.")
