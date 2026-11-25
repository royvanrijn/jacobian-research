#!/usr/bin/env python3
"""Exact retained complements for tighter terminal-thirteen light caps."""

from argparse import ArgumentParser
from ctypes import CDLL
from fractions import Fraction as Q
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

from verify_multivariate_terminal_thirteen_retained_bands import bracket_tail, retained_cores
from verify_multivariate_terminal_thirteen_retained_profiles import (
    refined_model, scaled_native, verify_interval, verify_parameters,
)
from verify_multivariate_terminal_six_private_angles import require


# Outer endpoints, lower and upper caps, and expected retained lists.
BANDS = {
    "53-54": (Q(53,100), Q(27,50), Q(1,30), Q(1,22), 10),
    "54-55": (Q(27,50), Q(11,20), Q(1,30), Q(1,22), 10),
    "56-57": (Q(14,25), Q(57,100), Q(1,40), Q(1,20), 67),
    "57-60": (Q(57,100), Q(3,5), Q(1,40), Q(1,20), 49),
    "60-62": (Q(3,5), Q(31,50), Q(1,40), Q(1,30), 14),
    "62-63": (Q(31,50), Q(63,100), Q(1,40), Q(1,30), 13),
    "66-67": (Q(33,50), Q(67,100), Q(1,50), Q(1,30), 51),
    "67-68": (Q(67,100), Q(17,25), Q(1,50), Q(1,30), 38),
}


def verify(name, selection, replay):
    lower, upper, cap, limit, expected = BANDS[name]
    lists = retained_cores(lower, cap, limit)
    require(len(lists) == expected, ("retained list count", name, len(lists)))
    indices = selection if selection else range(len(lists))
    require(all(0 <= index < len(lists) for index in indices), "retained list indices")
    outer_cells = radial_cells = disk_cells = 0
    for index in indices:
        rows, squared = lists[index]
        residual = bracket_tail(lower, upper, squared)
        outer, radial, disks = verify_interval(
            rows, ((), ()), residual, cap, replay, model_builder=refined_model,
            grouped_depth=0, maximum_grouped_width=Q(1,200))
        outer_cells += outer
        radial_cells += radial
        disk_cells += disks
        print(name, index, rows, outer, radial, disks, flush=True)
    print(f"Verified {name}, {len(indices)} of {len(lists)} lists, {outer_cells} outer cells, "
          f"{radial_cells} radial cells, and {disk_cells} outer disk covers.", flush=True)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("band", choices=BANDS, nargs="?")
    parser.add_argument("indices", type=int, nargs="*")
    args = parser.parse_args()
    require(args.band or not args.indices, "a band is required for list selection")
    verify_parameters()
    with TemporaryDirectory(prefix="sfc-thirteen-retained-fallbacks-") as directory:
        source = Path(__file__).with_name("verify_multivariate_terminal_ten_shared_light.cpp")
        library = Path(directory) / "bernstein.so"
        run(("c++", "-std=c++17", "-O3", "-fPIC", "-shared", str(source), "-o", str(library)),
            check=True)
        replay = scaled_native(CDLL(str(library)), compress=True)
        for name in (args.band,) if args.band else BANDS:
            verify(name, args.indices, replay)


if __name__ == "__main__":
    main()
