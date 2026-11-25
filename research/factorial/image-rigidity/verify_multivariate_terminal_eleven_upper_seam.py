#!/usr/bin/env python3
"""Exact terminal-eleven Cayley bridge from outer mean .94 to .95."""

from fractions import Fraction as Q
from sys import argv

import verify_multivariate_terminal_eleven_stable_upper as stable


RAYS = {
    50: (119137, 108477, 88273, 81314, 61691, 54840, 37590, 30097, 14473),
    55: (119515, 108787, 88459, 81473, 61780, 54913, 37629, 30126, 14485),
    60: (119813, 109103, 88649, 81635, 61871, 55629, 37669, 30156, 14498),
}
EXPECTED_CELLS = {
    (10, 50): (734, 65, 161, 62, 104, 56, 67, 37, 23),
    (10, 55): (1256, 65, 169, 62, 133, 60, 79, 40, 23),
    (11, 50): (51463, 81, 763, 223, 249, 74, 105, 35, 32),
    (11, 55): (112115, 85, 1394, 250, 435, 78, 179, 36, 42),
}


def main():
    order, ray = map(int, argv[1:]) if len(argv) > 1 else (0, -1)
    stable.require(order in (0, *stable.J) and -1 <= ray < 9, "replay selection")
    stable.require(RAYS[50] == stable.RAYS[50], "exact .95 seam")
    stable.LOWER_DEFICIT = Q(3, 50)
    stable.LAMBDA_TOP = Q(2449490, 10**7)
    stable.QMAX = Q(61, 1000)
    stable.RATIO = Q(100919, 100000), Q(101145, 100000)
    stable.C = 1 + stable.QMAX/(2*(1-stable.QMAX))
    stable.RAYS = RAYS
    stable.terms.cache_clear()
    stable.raw_error.cache_clear()
    stable.model.cache_clear()
    stable.verify_constants()
    total = 0
    for (current, left), counts in EXPECTED_CELLS.items():
        if order and current != order:
            continue
        for index, receipt in enumerate(counts):
            if ray >= 0 and ray != index:
                continue
            count = stable.check_ray(current, index, left, left+5)
            stable.require(count == receipt, ("ray count", current, left, index, count))
            total += count
            print(current, left, left+5, index, count, flush=True)
    print(f"Verified {total} exact projection cells.")


if __name__ == "__main__":
    main()
