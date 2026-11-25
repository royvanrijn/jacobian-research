#!/usr/bin/env python3
"""Exact retained complements for shared terminal thirteen."""

from argparse import ArgumentParser
from ctypes import CDLL
from fractions import Fraction as Q
from functools import lru_cache, partial
from math import prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

from verify_multivariate_terminal_thirteen_retained_profiles import (
    CAP, LOWER, coefficient_model, refined_model, replay_python, scaled_native, verify_interval,
    verify_parameters,
)
from verify_multivariate_terminal_six_private_angles import require
from verify_multivariate_terminal_ten_retained_mid_bands import bracket_tail


BANDS = {
    "lower": (Q(9,50), Q(27,100), Q(3,40), CAP),
    "lower-middle": (Q(27,100), Q(43,100), Q(3,50), CAP),
    "middle": (Q(43,100), Q(3,5), Q(1,20), CAP),
    "middle-low-tight": (Q(47,100), Q(1,2), Q(1,22), Q(1,20)),
    "middle-high-tight": (Q(13,25), Q(11,20), Q(1,22), Q(1,20)),
    "upper-middle": (Q(3,5), Q(33,50), Q(1,30), CAP),
    "bridge-left": (Q(33,50), Q(18,25), Q(1,30), CAP),
    "bridge-left-tight": (Q(33,50), Q(18,25), Q(1,36), Q(1,30)),
    "bridge-right": (Q(18,25), Q(19,25), Q(1,34), CAP),
    "cayley-middle": (Q(19,25), Q(171,200), Q(1,55), CAP),
    "cayley-middle-tight": (Q(171,200), Q(22,25), Q(1,60), CAP),
    "cayley-upper": (Q(22,25), Q(23,25), Q(1,70), CAP),
    "cayley-upper-tight": (Q(23,25), Q(19,20), Q(1,90), CAP),
}


@lru_cache(maxsize=None)
def row_square(base, gap):
    return prod(Q(LOWER * base + h, LOWER * base + (LOWER-1) * gap + h) for h in range(1, gap + 1))


def retained_cores(lower, cap, limit=CAP, *, canonical=True):
    gap = 1
    threshold = 1 + (LOWER-1) * cap / (1 + cap)
    candidates = []
    while threshold**(-gap) >= lower**2:
        base = 1
        while LOWER * base + 1 < gap / cap:
            if Q(gap, LOWER * base + 1) <= limit:
                squared = row_square(base, gap)
                if squared >= lower**2:
                    candidates.extend(((base, gap, sign), squared) for sign in (-1, 1))
            base += 1
        gap += 1
    candidates.sort()

    def collect(index, rows, squared):
        if index == len(candidates):
            if rows and (not canonical or rows <= tuple(sorted((a, d, -sign) for a, d, sign in rows))):
                yield rows, squared
            return
        row, row_squared = candidates[index]
        require(0 < row_squared < 1, ("row squared mean", row))
        while squared >= lower**2:
            yield from collect(index + 1, rows, squared)
            rows += (row,)
            squared *= row_squared

    return tuple(sorted(collect(0, (), Q(1))))


PRECISE_MODEL = partial(coefficient_model, precise=True, monotone=True)
REFINED_BANDS = {"bridge-left", "bridge-left-tight", "bridge-right",
                 "middle-low-tight", "middle-high-tight"}


def residual_cells(rows, residual, cap):
    """Retain an additional row when the first retained row leaves a broad tail."""
    if rows != ((1,1,-1),) or cap != Q(1,20):
        yield rows, residual, cap
        return
    smaller_cap = Q(1,30)
    yield rows, residual, smaller_cap
    for additional, squared in retained_cores(residual[0], smaller_cap, cap, canonical=False):
        yield tuple(sorted(rows+additional)), bracket_tail(*residual,squared), smaller_cap


def verify(replay, name, selection):
    verify_parameters()
    outer_sum = radial_sum = disk_sum = 0
    for band in [name] if name else BANDS:
        lower, upper, cap, limit = BANDS[band]
        model = refined_model if band in REFINED_BANDS else PRECISE_MODEL
        cores = retained_cores(lower, cap, limit)
        indices = list(range(len(cores))) if not selection else selection
        require(all(0 <= i < len(cores) for i in indices), ("valid list indices", band))
        outer_cells = radial_cells = disk_cells = 0
        for index in indices:
            rows, squared = cores[index]
            residual = bracket_tail(lower, upper, squared)
            leaves = cells = disks = 0
            for core, tail, tail_cap in residual_cells(rows, residual, cap):
                l, r, d = verify_interval(core, ((), ()), tail, tail_cap, replay,
                                          model_builder=model, grouped_depth=0,
                                          maximum_grouped_width=Q(1,200))
                leaves += l
                cells += r
                disks += d
            outer_cells += leaves
            radial_cells += cells
            disk_cells += disks
            print(band, index, rows, residual, leaves, cells, disks, flush=True)
        print(f"Verified {band}. {len(indices)} of {len(cores)} lists, "
              f"{outer_cells} outer cells, {radial_cells} radial cells, and {disk_cells} outer disk covers.", flush=True)
        outer_sum += outer_cells
        radial_sum += radial_cells
        disk_sum += disk_cells
    if not name:
        print(f"Verified {outer_sum} outer cells, {radial_sum} radial cells, and {disk_sum} outer disk covers.")


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--native", action="store_true")
    parser.add_argument("band", choices=BANDS, nargs="?")
    parser.add_argument("indices", type=int, nargs="*")
    args = parser.parse_args()
    require(args.band or not args.indices, "a band is required when listing indices")
    if not args.native:
        verify(replay_python, args.band, args.indices)
        return
    with TemporaryDirectory(prefix="sfc-thirteen-retained-bands-") as directory:
        source = Path(__file__).with_name("verify_multivariate_terminal_ten_shared_light.cpp")
        library = Path(directory) / "bernstein.so"
        run(("c++", "-std=c++17", "-O3", "-fPIC", "-shared", str(source), "-o", str(library)),
            check=True)
        verify(scaled_native(CDLL(str(library)), compress=True), args.band, args.indices)


if __name__ == "__main__":
    main()
