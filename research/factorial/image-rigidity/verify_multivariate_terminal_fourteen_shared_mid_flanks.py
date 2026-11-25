#!/usr/bin/env python3
"""Exact terminal-fourteen shared separation adjoining the middle bands."""

from argparse import ArgumentParser
from fractions import Fraction as Q
from tempfile import TemporaryDirectory

from flint import ctx

from verify_multivariate_terminal_fourteen_retained_profiles import (
    build_projector, configure_replay,
)
from verify_multivariate_terminal_fourteen_shared_mid import verify_disks, verify_retained


def main():
    import verify_multivariate_terminal_thirteen_shared_light as light

    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--part", choices=("all", "direct", "retained"), default="all")
    args = parser.parse_args()
    configure_replay()
    if args.part in ("all", "direct"):
        verify_disks(light, 45, Q(1, 30), (350, 384, 416, 773, 717))
    if args.part in ("all", "retained"):
        ctx.prec = 192
        light.FLOOR = Q(1, 150)
        with TemporaryDirectory(prefix="sfc-fourteen-shared-mid-flanks-") as directory:
            replay = build_projector(directory)
            verify_retained(45, 50, Q(1, 30), (39, 915, 157316), replay)


if __name__ == "__main__":
    main()
