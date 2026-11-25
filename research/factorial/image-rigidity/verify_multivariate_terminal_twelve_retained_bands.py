#!/usr/bin/env python3
"""Exact retained complements for shared terminal twelve."""

from argparse import ArgumentParser
from ctypes import CDLL
from fractions import Fraction as Q
from functools import lru_cache, partial
from math import prod
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

from verify_multivariate_terminal_twelve_retained_profiles import (
    CAP, ORDERS, anchor_angles, bounds, cap_at, coefficient_model,
    enclosing_annulus, replay_python, scaled_native, verify_parameters, verify_rays,
)
from verify_multivariate_terminal_six_private_angles import require
from verify_multivariate_terminal_ten_retained_mid_bands import bracket_tail, core_factor



BANDS = {
    "low-tight": (Q(7,20), Q(2,5), Q(17,200), CAP),
    "low-narrow": (Q(9,25), Q(2,5), Q(3,40), Q(17,200)),
    "low-reduced": (Q(9,25), Q(2,5), Q(3,50), Q(3,40)),
    "middle-lower": (Q(2,5), Q(14,25), Q(3,40), CAP),
    "middle-lower-short": (Q(23,50), Q(12,25), Q(13,200), Q(3,40)),
    "middle-lower-tight": (Q(49,100), Q(14,25), Q(11,200), Q(3,40)),
    "middle-lower-narrow": (Q(1,2), Q(14,25), Q(1,20), Q(11,200)),
    "middle": (Q(14,25), Q(7,10), Q(1,20), CAP),
    "middle-narrow": (Q(33,50), Q(7,10), Q(1,25), Q(1,20)),
    "middle-upper": (Q(7,10), Q(19,25), Q(3,80), CAP),
    "middle-tight": (Q(7,10), Q(19,25), Q(3,100), Q(3,80)),
    "cayley-middle": (Q(19,25), Q(22,25), Q(1,55), CAP),
    "cayley-middle-tight": (Q(171,200), Q(22,25), Q(1,60), Q(1,55)),
    "cayley-upper": (Q(22,25), Q(19,20), Q(1,70), CAP),
    "cayley-upper-tight": (Q(23,25), Q(19,20), Q(1,90), CAP),
}

EXPECTED_LISTS = (2, 5, 1, 6, 1, 1, 1, 2, 4, 2, 1, 14, 1, 3, 4)


@lru_cache(maxsize=None)
def row_square(base, gap):
    return prod(Q(11 * base + h, 11 * base + 10 * gap + h) for h in range(1, gap + 1))


def retained_cores(lower, cap, limit=CAP):
    gap = 1
    threshold = 1 + 10 * cap / (1 + cap)
    candidates = []
    while threshold**(-gap) >= lower**2:
        base = 1
        while 11 * base + 1 < gap / cap:
            if Q(gap, 11 * base + 1) <= limit:
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




PRECISE_MODEL = partial(coefficient_model, precise=True, monotone=True)


def verify_core(rows, outer, declared, replay):
    core = {order: core_factor(order, rows) for order in ORDERS}
    global_annulus = enclosing_annulus(bounds(core, outer, cap_at(outer[0], declared)))
    pending = [(outer, 0)]
    leaves = cells = 0
    while pending:
        means, depth = pending.pop()
        midpoint = sum(means) / 2
        count = None
        if means[1] - means[0] <= Q(1, 200):
            cap = cap_at(means[0], declared)
            targets = anchor_angles(rows, False, midpoint)
            if targets is not None:
                box = bounds(core, means, cap)
                if global_annulus or enclosing_annulus(box):
                    count = verify_rays(core, box, means, cap, targets, True, replay,
                                        model_builder=PRECISE_MODEL)
        if count is not None:
            leaves += 1
            cells += count
            continue
        require(depth < 14, ("unresolved retained cell", rows, means))
        pending.extend((((means[0], midpoint), depth + 1), ((midpoint, means[1]), depth + 1)))
    return leaves, cells



def verify(replay, name, selection):
    verify_parameters()
    outer_sum = radial_sum = 0
    for band in [name] if name else BANDS:
        lower, upper, cap, limit = BANDS[band]
        cores = retained_cores(lower, cap, limit)
        expected = dict(zip(BANDS, EXPECTED_LISTS, strict=True))[band]
        require(len(cores) == expected, ("retained list count", band, len(cores)))
        indices = list(range(len(cores))) if not selection else selection
        require(all(0 <= i < len(cores) for i in indices), ("valid list indices", band))
        outer_cells = radial_cells = 0
        for index in indices:
            rows, squared = cores[index]
            residual = bracket_tail(lower, upper, squared)
            leaves, cells = verify_core(rows, residual, cap, replay)
            outer_cells += leaves
            radial_cells += cells
            print(band, index, rows, residual, leaves, cells, flush=True)
        print(f"Verified {band}. {len(indices)} of {len(cores)} lists, "
              f"{outer_cells} outer cells and {radial_cells} radial cells.", flush=True)
        outer_sum += outer_cells
        radial_sum += radial_cells
    if not name:
        print(f"Verified {outer_sum} outer cells and {radial_sum} radial cells in total.")


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
    with TemporaryDirectory(prefix="sfc-twelve-retained-bands-") as directory:
        source = Path(__file__).with_name("verify_multivariate_terminal_ten_shared_light.cpp")
        library = Path(directory) / "bernstein.so"
        run(("c++", "-std=c++17", "-O3", "-fPIC", "-shared", str(source), "-o", str(library)),
            check=True)
        verify(scaled_native(CDLL(str(library)), compress=True), args.band, args.indices)


if __name__ == "__main__":
    main()
