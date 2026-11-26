#!/usr/bin/env python3
"""Exact terminal-fourteen separation on the remaining small outer band."""

from argparse import ArgumentParser
from fractions import Fraction as Q
from tempfile import TemporaryDirectory

import verify_multivariate_terminal_fourteen_shared_low_u as disks
import verify_multivariate_terminal_thirteen_shared_light as light
from verify_multivariate_terminal_fourteen_retained_profiles import (
    CAP, build_projector, configure_replay, retained,
)


LIGHT_CAP = Q(2, 25)
RAY_BANDS = (
    (20, 24, (1, 111, 0)),
    (24, 28, (1, 116, 0)),
    (28, 32, (1, 124, 0)),
)
DISK_BANDS = ((32, 36, 18), (36, 40, 9))


def replay_rays(rows, means, replay):
    return retained.verify_interval(
        rows, ((), ()), means, LIGHT_CAP, replay,
        model_builder=retained.refined_model, grouped_depth=0,
        maximum_grouped_width=Q(1, 200),
    )


def verify_light_rays(replay):
    for start, stop, expected in RAY_BANDS:
        counts = replay_rays((), (Q(start, 1000), Q(stop, 1000)), replay)
        retained.require(counts == expected, ("empty core ray counts", start, counts))
        print("Light rays", start, stop, counts, flush=True)


def verify_retained_rays(replay):
    # The enumerator binds its order on import.
    from verify_multivariate_terminal_thirteen_retained_bands import bracket_tail, retained_cores

    lower, upper = Q(1, 50), Q(1, 25)
    lists = retained_cores(lower, LIGHT_CAP, CAP)
    outer = radial = 0
    for index, (rows, squared) in enumerate(lists):
        means = bracket_tail(lower, upper, squared)
        cells, subdivisions, disk_cells = replay_rays(rows, means, replay)
        retained.require(disk_cells == 0, ("retained ray replay", rows, disk_cells))
        outer += cells
        radial += subdivisions
        print("Retained", index, rows, cells, subdivisions, flush=True)
    counts = len(lists), outer, radial
    retained.require(counts == (43, 66, 6903), ("retained ray counts", counts))
    print("Retained totals", counts, flush=True)


def verify_light_disks():
    disks.ctx.prec = 112
    light.exact_cap = disks.exact_cap
    light.ordinary_ratio_upper = disks.ordinary_ratio_upper
    for start, stop, expected in DISK_BANDS:
        counts = disks.check_band(light, Q(start, 1000), Q(stop, 1000), LIGHT_CAP)
        retained.require(counts == expected, ("light disk counts", start, counts))
        print("Light disks", start, stop, counts, flush=True)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--part", choices=("all", "light", "retained"), default="all")
    args = parser.parse_args()
    configure_replay()
    light.FLOOR = Q(1, 150)
    with TemporaryDirectory(prefix="sfc-fourteen-corner-") as directory:
        replay = build_projector(directory)
        if args.part in ("all", "light"):
            verify_light_rays(replay)
        if args.part in ("all", "retained"):
            verify_retained_rays(replay)
    if args.part in ("all", "light"):
        verify_light_disks()


if __name__ == "__main__":
    main()
