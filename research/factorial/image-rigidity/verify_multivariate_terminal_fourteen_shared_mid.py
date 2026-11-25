#!/usr/bin/env python3
"""Exact paired separation for terminal-fourteen shared profiles on [.30,.45]."""

from argparse import ArgumentParser
from fractions import Fraction as Q
from tempfile import TemporaryDirectory

from flint import arb, ctx

from verify_multivariate_terminal_fourteen_retained_profiles import (
    CAP, build_projector, configure_replay, retained,
)
from verify_multivariate_terminal_fourteen_shared_disks import check_band


# Lower hundredth, lower-order reciprocal cap, disk counts per hundredth,
# and retained list, outer cell and radial cell counts.
BANDS = {
    30: (Q(1, 22), (368, 360, 307, 307, 354), (25, 768, 115196)),
    35: (Q(1, 25), (153, 330, 399, 416, 402), (21, 616, 98164)),
    40: (Q(1, 25), (368, 336, 344, 376, 403), (10, 294, 51701)),
}


def verify_disks(model, start, cap, counts):
    ctx.prec = 112
    model.FLOOR = Q(1, 80)
    radii = tuple(arb(k)/2**p for p in range(20, 4, -1) for k in range(4, 8))
    for index, expected in enumerate(counts, start=start):
        count = check_band(model, Q(index, 100), Q(index+1, 100), cap, radii)
        retained.require(count == expected, ("shared disk count", index, count, expected))
        print("Direct", index, index+1, cap, count, flush=True)
        model.coefficient_model.cache_clear()
        model.ordinary_ratio_upper.cache_clear()


def verify_retained(start, stop, cap, expected, replay):
    # Import only after binding the enumerator to lower order thirteen.
    from verify_multivariate_terminal_thirteen_retained_bands import bracket_tail, retained_cores

    lower, upper = Q(start, 100), Q(stop, 100)
    lists = retained_cores(lower, cap, CAP)
    retained.require(len(lists) == expected[0], ("retained list count", start, len(lists)))
    outer = radial = 0
    for rows, squared in lists:
        residual = bracket_tail(lower, upper, squared)
        l, r, d = retained.verify_interval(
            rows, ((), ()), residual, cap, replay,
            model_builder=retained.refined_model, grouped_depth=0,
            maximum_grouped_width=Q(1, 200))
        retained.require(d == 0, "middle retained ray replay")
        outer += l
        radial += r
        print("Retained", start, rows, l, r, flush=True)
    retained.require((len(lists), outer, radial) == expected,
                     ("retained cell counts", start, outer, radial))
    print(f"Verified {start}–{stop}: {len(lists)} retained lists, "
          f"{outer} outer cells and {radial} radial cells.", flush=True)


def main():
    import verify_multivariate_terminal_thirteen_shared_light as light

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("start", nargs="?", type=int, choices=BANDS)
    parser.add_argument("--part", choices=("all", "direct", "retained"), default="all")
    args = parser.parse_args()
    configure_replay()
    starts = (args.start,) if args.start is not None else BANDS
    if args.part in ("all", "direct"):
        for start in starts:
            cap, counts, _ = BANDS[start]
            verify_disks(light, start, cap, counts)
    if args.part in ("all", "retained"):
        ctx.prec = 192
        light.FLOOR = Q(1, 150)
        with TemporaryDirectory(prefix="sfc-fourteen-shared-mid-") as directory:
            replay = build_projector(directory)
            for start in starts:
                cap, counts, expected = BANDS[start]
                verify_retained(start, start+len(counts), cap, expected, replay)


if __name__ == "__main__":
    main()
