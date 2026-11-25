#!/usr/bin/env python3
"""Exact terminal-fourteen retained complements of the lower-middle light bands."""

from argparse import ArgumentParser
from fractions import Fraction as Q
from tempfile import TemporaryDirectory

from verify_multivariate_terminal_fourteen_retained_bridge import verify_centered_disks
from verify_multivariate_terminal_fourteen_retained_profiles import (
    CAP, ORDERS, build_projector, configure_replay, retained,
)


BANDS = {
    "50-55": (Q(1, 2), Q(11, 20), Q(1, 40), 66, 1207, 213668, 0),
    "55-60": (Q(11, 20), Q(3, 5), Q(1, 40), 32, 553, 116072, 0),
    "60-63": (Q(3, 5), Q(63, 100), Q(1, 40), 17, 253, 51093, 64),
    "63-64": (Q(63, 100), Q(16, 25), Q(1, 44), 28, 116, 28151, 0),
}


def verify_singleton(rows, means, declared):
    endpoints = (Q(511, 625), Q(21, 25), Q(4293, 5000))
    retained.require(means == (endpoints[0], endpoints[-1]), ("lower singleton coverage", means))
    powers = {order: retained.fixed_powers(order, rows, ((), ())) for order in ORDERS}
    for shard in zip(endpoints[:-1], endpoints[1:], strict=True):
        pending = [(shard, 0)]
        accepted = 0
        while pending:
            cell, depth = pending.pop()
            cap = retained.cap_at(cell[0], declared)
            if cell[1]-cell[0] <= Q(1, 1000) and verify_centered_disks(powers, cell, cap):
                accepted += 1
                continue
            retained.require(depth < 17, ("unresolved lower singleton", cell, cap))
            midpoint = sum(cell)/2
            pending.extend((((cell[0], midpoint), depth+1), ((midpoint, cell[1]), depth+1)))
        retained.require(accepted == 32, ("lower singleton disk count", shard, accepted))
        print("Centered lower singleton shard", shard, accepted, flush=True)
    return 64, 0, 64


def verify_band(name, replay):
    # The band enumerator imports the lower order after configuration.
    from verify_multivariate_terminal_thirteen_retained_bands import bracket_tail, retained_cores

    lower, upper, cap, expected_lists, expected_outer, expected_radial, expected_disks = BANDS[name]
    lists = retained_cores(lower, cap, CAP)
    retained.require(len(lists) == expected_lists, ("lower-middle retained lists", name, len(lists)))
    outer = radial = disks = 0
    for rows, squared in lists:
        residual = bracket_tail(lower, upper, squared)
        if name == "60-63" and rows == ((1, 1, -1),):
            l, r, d = verify_singleton(rows, residual, cap)
        else:
            l, r, d = retained.verify_interval(
                rows, ((), ()), residual, cap, replay,
                model_builder=retained.refined_model, grouped_depth=0,
                maximum_grouped_width=Q(1, 200))
            retained.require(d == 0, "lower-middle ray replay")
        outer += l
        radial += r
        disks += d
        print(name, rows, l, r, d, flush=True)
    retained.require((outer, radial, disks) == (expected_outer, expected_radial, expected_disks),
                     ("lower-middle cell counts", name, outer, radial, disks))
    print(f"Verified {name}: {len(lists)} lists, {outer} outer cells, "
          f"{radial} radial cells, and {disks} disk covers.", flush=True)


def main():
    import verify_multivariate_terminal_thirteen_shared_light as light

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("band", choices=BANDS, nargs="?")
    args = parser.parse_args()
    configure_replay()
    light.FLOOR = Q(1, 150)
    with TemporaryDirectory(prefix="sfc-fourteen-retained-low-mid-") as directory:
        replay = build_projector(directory)
        for name in (args.band,) if args.band else BANDS:
            verify_band(name, replay)


if __name__ == "__main__":
    main()
