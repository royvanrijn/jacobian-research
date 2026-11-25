#!/usr/bin/env python3
"""Exact replay for high and tightened retained bands at terminal eleven."""

from ctypes import CDLL
from fractions import Fraction as Q
from functools import partial
from pathlib import Path
from subprocess import run
from sys import argv
from tempfile import TemporaryDirectory

from verify_multivariate_terminal_eleven_retained_mid_bands import retained_cores
from verify_multivariate_terminal_eleven_retained_profiles import (
    ANNULUS, ORDERS, anchor_angles, bounds, cap_at, verify_parameters,
    verify_rays, replay_python,
)
from verify_multivariate_terminal_eleven_shared_light import coefficient_model, scaled_native
from verify_multivariate_terminal_six_private_angles import require
from verify_multivariate_terminal_ten_retained_mid_bands import bracket_tail, core_factor


HIGH_BAND = (Q(22,25), Q(19,20), Q(1,70), Q(1,8))
UPPER_MIDDLE_BAND = (Q(19,25), Q(22,25), Q(1,55), Q(1,8))
TIGHTENED_BANDS = (
    (Q(21,50), Q(11,25), Q(93,1000), Q(1,8)),
    (Q(57,100), Q(16,25), Q(3,50), Q(9,100)),
    (Q(73,100), Q(19,25), Q(9,200), Q(11,200)),
)
HIGH_LISTS = tuple(((base,1,-1),) for base in range(3,7))
TIGHTENED_LISTS = ((((2,2,-1),),), (((3,2,-1),),), (((2,1,-1),),))
PRECISE_MODEL = partial(coefficient_model, precise=True, monotone=True)


def verify_core(rows, outer, declared, replay):
    core = {order: core_factor(order, rows) for order in ORDERS}
    inner, upper = ANNULUS
    entire = bounds(core, outer, cap_at(outer[0], declared))
    for order in ORDERS:
        require(all(a > inner*d and b < upper*c for (a,b),(c,d) in
                    zip(entire[order], entire[order][1:])),
                ("Enestrom annulus", rows, order, entire[order]))

    pending = [(outer, 0)]
    outer_cells = radial_cells = 0
    while pending:
        means, depth = pending.pop()
        midpoint = sum(means) / 2
        count = None
        if means[1] - means[0] <= Q(1, 200):
            cap = cap_at(means[0], declared)
            targets = anchor_angles(rows, False, midpoint)
            if targets is not None:
                box = bounds(core, means, cap)
                count = verify_rays(core, box, means, cap, targets, True, replay,
                                    model_builder=PRECISE_MODEL)
        if count is not None:
            outer_cells += 1
            radial_cells += count
            continue
        require(depth < 14, ("unresolved high retained cell", rows, means))
        pending.extend((((means[0], midpoint), depth + 1), ((midpoint, means[1]), depth + 1)))
    return outer_cells, radial_cells


def verify(replay, selected, family):
    verify_parameters()
    bands, lists = {
        "high": ((HIGH_BAND,), (HIGH_LISTS,)),
        "upper-middle": ((UPPER_MIDDLE_BAND,), (37,)),
        "tightened": (TIGHTENED_BANDS, TIGHTENED_LISTS),
    }[family]
    outer_cells = radial_cells = 0
    for (lower, upper, declared, limit), expected in zip(bands, lists, strict=True):
        cores = tuple((rows,squared) for rows,squared in retained_cores(lower, declared)
                      if all(Q(d,10*a+1) <= limit for a,d,_ in rows))
        found = len(cores) if isinstance(expected, int) else tuple(rows for rows, _ in cores)
        require(found == expected, ("retained lists", cores))
        for index, (rows, squared) in enumerate(cores):
            key = index if family == "upper-middle" else rows[0][0]
            if selected and key not in selected:
                continue
            outer = bracket_tail(lower, upper, squared)
            count, radial = verify_core(rows, outer, declared, replay)
            outer_cells += count
            radial_cells += radial
            print(lower, upper, rows, outer, count, radial, flush=True)
    print(f"Verified {outer_cells} outer cells and {radial_cells} radial cells.", flush=True)


def main():
    switches = ("--tightened", "--upper-middle")
    selected = tuple(map(int, (x for x in argv[1:] if x not in ("--native", *switches))))
    active = [switch.removeprefix("--") for switch in switches if switch in argv]
    require(len(active) <= 1, "at most one alternative family")
    family = active[0] if active else "high"
    valid = {"high": {3,4,5,6}, "upper-middle": set(range(37)), "tightened": set()}[family]
    require(set(selected) <= valid, ("valid selections", family))
    if "--native" not in argv:
        verify(replay_python, selected, family)
        return
    with TemporaryDirectory(prefix="sfc-eleven-retained-high-") as directory:
        library = Path(directory) / "bernstein.so"
        source = Path(__file__).with_name("verify_multivariate_terminal_ten_shared_light.cpp")
        run(("c++", "-std=c++17", "-O3", "-fPIC", "-shared", str(source), "-o", str(library)),
            check=True)
        verify(scaled_native(CDLL(str(library)), compress=True), selected, family)


if __name__ == "__main__":
    main()
