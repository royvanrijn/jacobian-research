#!/usr/bin/env python3
"""Terminal-fourteen retained complements on five low outer-mean bands."""

from argparse import ArgumentParser
from fractions import Fraction as Q
from tempfile import TemporaryDirectory

from verify_multivariate_terminal_fourteen_retained_profiles import build_projector, configure_replay
from verify_multivariate_terminal_fourteen_shared_mid import verify_retained


BANDS = {
    10: (11, Q(1, 9), (1, 4, 489)),
    13: (15, Q(1, 10), (1, 16, 2064)),
    15: (16, Q(1, 11), (3, 32, 4239)),
    16: (17, Q(3, 40), (4, 48, 6474)),
    18: (20, Q(3, 40), (4, 96, 12805)),
}


def main():
    import verify_multivariate_terminal_thirteen_shared_light as light

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("start", nargs="?", type=int, choices=BANDS)
    args = parser.parse_args()
    configure_replay()
    light.FLOOR = Q(1, 150)
    starts = (args.start,) if args.start is not None else BANDS
    with TemporaryDirectory(prefix="sfc-fourteen-retained-low-") as directory:
        replay = build_projector(directory)
        for start in starts:
            stop, cap, expected = BANDS[start]
            verify_retained(start, stop, cap, expected, replay)


if __name__ == "__main__":
    main()
