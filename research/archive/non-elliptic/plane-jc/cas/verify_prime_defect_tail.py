#!/usr/bin/env python3
"""Verify the uniform arithmetic behind the omitted GGHV tail census.

For a current corner C=(A,B) and a proper proportional continuation
C'=(a,b), put

    g=gcd(A-a,B-b),  (p,q)=((A-a)/g,(B-b)/g).

The auxiliary-F divisibility used in GGV/GGHV reduces to

    B*p - A*q  |  A-B.

The script enumerates all positive continuations for the two relevant corners
(21,8) and (24,7), and verifies the exact k=1 parallelism conclusion for the
(24,7) tail.
"""

from __future__ import annotations

from math import gcd


def proper_continuations(A: int, B: int):
    result = []
    for p in range(1, A + 1):
        for q in range(0, B + 1):
            if gcd(p, q) != 1:
                continue
            defect = B * p - A * q
            if defect <= 0 or (A - B) % defect:
                continue
            for scale in range(1, A + B + 1):
                a = A - scale * p
                b = B - scale * q
                if a <= 0 or b < 0:
                    break
                result.append((a, b, p, q, scale, defect))
    return result


def determinant_direct(a: int, b: int, k: int) -> int:
    # en(P)=(-k,0), en(Q)=(k+1,1)
    return -2 * a + (5 * k + 2) * b - k


def determinant_swapped_twice(a: int, b: int, k: int) -> int:
    # Twice the determinant for the swapped endpoint assignment.  This keeps
    # the calculation integral on the affine line 2a=7b-1.
    return 2 * (3 * a - (5 * k + 3) * b + k)


def main() -> None:
    expected_21_8 = [
        (13, 5, 8, 3, 1, 1),
        (5, 2, 8, 3, 2, 1),
        (1, 1, 20, 7, 1, 13),
    ]
    expected_24_7 = [
        (17, 5, 7, 2, 1, 1),
        (10, 3, 7, 2, 2, 1),
        (3, 1, 7, 2, 3, 1),
        (1, 1, 23, 6, 1, 17),
    ]
    assert proper_continuations(21, 8) == expected_21_8
    assert proper_continuations(24, 7) == expected_24_7

    proper = [(17, 5), (10, 3), (3, 1)]
    for a, b in proper:
        assert 2 * a == 7 * b - 1
        for k in (1, 2):
            direct = determinant_direct(a, b, k)
            swapped_twice = determinant_swapped_twice(a, b, k)
            assert direct == (k - 1) * (5 * b - 1)
            assert swapped_twice == -(2 * k - 3) * (5 * b - 1)
        assert determinant_direct(a, b, 1) == 0
        assert determinant_swapped_twice(a, b, 1) != 0
        assert determinant_direct(a, b, 2) != 0
        assert determinant_swapped_twice(a, b, 2) != 0

    print("TAIL_21_8=(13,5),(5,2); diagonal=(1,1)")
    print("TAIL_24_7=(17,5),(10,3),(3,1); diagonal=(1,1)")
    print("PRIMITIVE_24_7_STEP=(7,2)")
    print("SUCCESSOR_NORMAL=(-2,7)")
    print("PARALLEL_ENDPOINTS=en(P)=(-1,0),en(Q)=(2,1)")
    print("PRIME_DEFECT_TAIL_PASS")


if __name__ == "__main__":
    main()
