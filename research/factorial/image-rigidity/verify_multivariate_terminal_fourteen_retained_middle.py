#!/usr/bin/env python3
"""Exact terminal-fourteen retained complements of the middle light band."""

from argparse import ArgumentParser
from fractions import Fraction as Q
from tempfile import TemporaryDirectory

from verify_multivariate_terminal_fourteen_retained_profiles import (
    CAP, build_projector, configure_replay, retained,
)


BANDS = {
    "64-72": (Q(16, 25), Q(18, 25), 126, 2073, 485724),
    "72-74": (Q(18, 25), Q(37, 50), 29, 219, 66515),
    "74-76": (Q(37, 50), Q(19, 25), 18, 132, 46323),
    "76-81": (Q(19, 25), Q(81, 100), 13, 218, 65795),
}
MIDDLE_CAP = Q(1, 60)


def verify_band(name, replay):
    # The band enumerator imports the lower order after configuration.
    from verify_multivariate_terminal_thirteen_retained_bands import bracket_tail, retained_cores

    lower, upper, expected_lists, expected_outer, expected_radial = BANDS[name]
    lists = retained_cores(lower, MIDDLE_CAP, CAP)
    retained.require(len(lists) == expected_lists, ("middle retained lists", name, len(lists)))
    outer = radial = 0
    for rows, squared in lists:
        residual = bracket_tail(lower, upper, squared)
        l, r, d = retained.verify_interval(
            rows, ((), ()), residual, MIDDLE_CAP, replay,
            model_builder=retained.refined_model, grouped_depth=0,
            maximum_grouped_width=Q(1, 200))
        retained.require(d == 0, "middle ray replay")
        outer += l
        radial += r
        print(name, rows, l, r, flush=True)
    retained.require((outer, radial) == (expected_outer, expected_radial),
                     ("middle cell counts", name, outer, radial))
    print(f"Verified {name}: {len(lists)} lists, {outer} outer and {radial} radial cells.",
          flush=True)


def main():
    import verify_multivariate_terminal_thirteen_shared_light as light

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("band", choices=BANDS, nargs="?")
    args = parser.parse_args()
    configure_replay()
    light.FLOOR = Q(1, 150)
    with TemporaryDirectory(prefix="sfc-fourteen-retained-middle-") as directory:
        replay = build_projector(directory)
        for name in (args.band,) if args.band else BANDS:
            verify_band(name, replay)


if __name__ == "__main__":
    main()
