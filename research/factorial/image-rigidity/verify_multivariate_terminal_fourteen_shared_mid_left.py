#!/usr/bin/env python3
"""Exact terminal-fourteen shared separation on [.25,.30]."""

from argparse import ArgumentParser
from fractions import Fraction as Q
from tempfile import TemporaryDirectory

from flint import ctx

from verify_multivariate_terminal_fourteen_retained_profiles import (
    build_projector, configure_replay, retained,
)
from verify_multivariate_terminal_fourteen_shared_mid import verify_disks, verify_retained


def verify_tightened(replay):
    from verify_multivariate_terminal_thirteen_retained_bands import bracket_tail, retained_cores

    lo, hi, cap, limit = Q(7, 25), Q(3, 10), Q(1, 20), Q(1, 18)
    lists = retained_cores(lo, cap, limit)
    retained.require(len(lists) == 2, ("tightened retained list count", len(lists)))
    outer = radial = 0
    for rows, squared in lists:
        residual = bracket_tail(lo, hi, squared)
        l, r, d = retained.verify_interval(
            rows, ((), ()), residual, cap, replay,
            model_builder=retained.refined_model, grouped_depth=0,
            maximum_grouped_width=Q(1, 200))
        retained.require(d == 0, "tightened retained ray replay")
        outer += l
        radial += r
        print("Tightened retained", rows, l, r, flush=True)
    retained.require((outer, radial) == (24, 3380), ("tightened retained counts", outer, radial))
    print("Verified two tightened lists, 24 outer and 3380 radial cells.", flush=True)


def main():
    import verify_multivariate_terminal_thirteen_shared_light as light

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--part", choices=("all", "direct", "retained"), default="all")
    args = parser.parse_args()
    configure_replay()
    if args.part in ("all", "direct"):
        verify_disks(light, 25, Q(1, 18), (332, 468, 944))
        verify_disks(light, 28, Q(1, 20), (416, 473))
    if args.part in ("all", "retained"):
        ctx.prec = 192
        light.FLOOR = Q(1, 150)
        with TemporaryDirectory(prefix="sfc-fourteen-shared-mid-left-") as directory:
            replay = build_projector(directory)
            verify_retained(25, 30, Q(1, 18), (31, 790, 111384), replay)
            verify_tightened(replay)


if __name__ == "__main__":
    main()
