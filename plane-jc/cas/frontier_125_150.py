#!/usr/bin/env python3
"""Exact bounded regression of the published <=150 candidate tables.

This is not a reimplementation of the complete-chain generator.  It encodes
the accepted output in Section 7 of arXiv:1708.07936, recomputes coordinate
degrees from A0 and (m,n), checks the printed maxima, quotients by reversal,
and emits precisely the distinct pairs in the requested [125,150] window.
"""

from dataclasses import dataclass
from math import gcd


@dataclass(frozen=True)
class Case:
    source: str
    a0: tuple[int, int]
    chain: tuple[str, ...]
    mn: tuple[int, int]
    printed_max: int

    @property
    def degrees(self) -> tuple[int, int]:
        scale = sum(self.a0)
        return self.mn[0] * scale, self.mn[1] * scale

    @property
    def normalized(self) -> tuple[int, int]:
        return tuple(sorted(self.degrees))


CASES = (
    Case("F2", (5, 20), ("(7/5,2)",), (3, 5), 125),
    Case("F7", (6, 15), ("(7/3,4)",), (2, 7), 147),
    Case("F8", (6, 15), ("(8/3,5)",), (3, 7), 147),
    Case("F9", (7, 21), ("(11/7,2)",), (3, 5), 140),
    Case("F11", (7, 21), ("(13/7,3)",), (2, 5), 140),
    Case("F24", (8, 24), ("(14/4,6)", "(5/4,0)", "(19/8,3)"), (3, 4), 128),
    Case("table-L1", (7, 35), ("(19/7,5)",), (2, 3), 126),
    Case("table-L1a", (7, 42), ("(13/7,6)",), (3, 2), 147),
    Case("table-L1b", (7, 42), ("(13/7,6)",), (2, 3), 147),
    Case("table-L1", (8, 28), ("(7/4,3)",), (3, 4), 144),
    Case("table-L1a", (9, 36), ("(17/9,4)",), (3, 2), 135),
    Case("table-L1b", (9, 36), ("(17/9,4)",), (2, 3), 135),
    Case("table-L1", (11, 33), ("(19/4,8)",), (2, 3), 132),
    Case("table-L1", (12, 33), ("(11/3,8)",), (2, 3), 135),
    Case("table-L2", (8, 40), ("(8,28)", "(11/4,7)"), (3, 2), 144),
    Case("table-L2", (9, 36), ("(9,24)", "(11/3,8)"), (2, 3), 135),
    Case("table-L2a", (10, 40), ("(16/5,6)", "(23/10,3)"), (3, 2), 150),
    Case("table-L2b", (10, 40), ("(18/5,8)", "(8/5,3)"), (3, 2), 150),
    Case("table-L2", (12, 30), ("(16/3,10)", "(11/6,3)"), (3, 2), 126),
    Case("table-L2a", (12, 36), ("(12,33)", "(11/3,8)"), (2, 3), 144),
    Case("table-L2b", (12, 36), ("(9,24)", "(11/3,8)"), (2, 3), 144),
    Case("table-L2c", (12, 36), ("(21/4,9)", "(19/4,8)"), (2, 3), 144),
    Case("table-L2d", (12, 36), ("(21/4,9)", "(12/4,5)"), (2, 3), 144),
    Case("table-L3", (12, 36), ("(12,30)", "(16/3,10)", "(11/6,3)"), (3, 2), 144),
)


def main() -> None:
    by_pair: dict[tuple[int, int], list[Case]] = {}
    for case in CASES:
        assert gcd(*case.mn) == 1
        assert max(case.degrees) == case.printed_max, case
        assert 125 <= case.printed_max <= 150
        by_pair.setdefault(case.normalized, []).append(case)

    expected = {
        (75, 125), (84, 126), (96, 128), (88, 132), (90, 135),
        (56, 140), (84, 140), (96, 144), (108, 144), (42, 147),
        (63, 147), (98, 147), (100, 150),
    }
    assert set(by_pair) == expected
    print("PAIR|GCD|RATIO|REALIZATIONS")
    for pair in sorted(by_pair, key=lambda p: (p[1], p[0])):
        g = gcd(*pair)
        ratio = f"{pair[0] // g}:{pair[1] // g}"
        print(f"{pair[0]},{pair[1]}|{g}|{ratio}|{len(by_pair[pair])}")
    print("FRONTIER_125_150_PASS")


if __name__ == "__main__":
    main()
