#!/usr/bin/env python3
"""Exact retained complements for shared terminal eleven on .44–.76."""

from ctypes import CDLL
from fractions import Fraction as Q
from functools import lru_cache
from math import prod
from pathlib import Path
from subprocess import run
from sys import argv
from tempfile import TemporaryDirectory

from verify_multivariate_terminal_eleven_retained_profiles import (
    CAP, normalized_native, replay_python, verify_interval, verify_parameters,
)
from verify_multivariate_terminal_six_private_angles import require
from verify_multivariate_terminal_ten_retained_mid_bands import bracket_tail


BANDS = (
    (Q(44, 100), Q(62, 100), Q(9, 100)),
    (Q(62, 100), Q(70, 100), Q(3, 40)),
    (Q(70, 100), Q(76, 100), Q(11, 200)),
)


@lru_cache(maxsize=None)
def row_square(base, gap):
    return prod(Q(10 * base + h, 10 * base + 9 * gap + h) for h in range(1, gap + 1))


def retained_cores(lower, cap):
    gap = 1
    threshold = 1 + 9 * cap / (1 + cap)
    candidates = []
    while threshold**(-gap) >= lower**2:
        base = 1
        while 10 * base + 1 < gap / cap:
            if Q(gap, 10 * base + 1) <= CAP:
                squared = row_square(base, gap)
                if squared >= lower**2:
                    candidates.extend(((base, gap, sign), squared) for sign in (-1, 1))
            base += 1
        gap += 1

    def collect(index, rows, squared):
        if index == len(candidates):
            if rows and rows <= tuple(sorted((a, d, -sign) for a, d, sign in rows)):
                yield rows, squared
            return
        row, row_squared = candidates[index]
        require(0 < row_squared < 1, ("row squared mean", row))
        while squared >= lower**2:
            yield from collect(index + 1, rows, squared)
            rows = tuple(sorted((*rows, row)))
            squared *= row_squared

    return tuple(sorted(collect(0, (), Q(1))))


def verify(replay):
    verify_parameters()
    leaves = cells = 0
    for (lower, upper, cap), expected in zip(BANDS, (4, 1, 1), strict=True):
        cores = retained_cores(lower, cap)
        require(len(cores) == expected, ("retained list count", lower, cap, cores))
        for rows, squared in cores:
            outer = bracket_tail(lower, upper, squared)
            count, radial = verify_interval(rows, False, outer, cap, replay)
            leaves += count
            cells += radial
            print(lower, upper, rows, outer, count, radial, flush=True)
    print(f"Verified {leaves} outer cells and {cells} radial cells.", flush=True)


def main():
    require(argv[1:] in ([], ["--native"]), "optional --native flag")
    if "--native" not in argv:
        verify(replay_python)
        return
    with TemporaryDirectory(prefix="sfc-eleven-retained-mid-") as directory:
        library = Path(directory) / "bernstein.so"
        source = Path(__file__).with_name("verify_multivariate_terminal_ten_shared_light.cpp")
        run(("c++", "-std=c++17", "-O3", "-fPIC", "-shared", str(source), "-o", str(library)),
            check=True)
        verify(normalized_native(CDLL(str(library))))


if __name__ == "__main__":
    main()
