#!/usr/bin/env python3
"""Exact terminal-fourteen replay for private and large-reciprocal profiles."""

from argparse import ArgumentParser
from ctypes import CDLL
from fractions import Fraction as Q
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

import verify_multivariate_terminal_thirteen_retained_profiles as retained
import verify_multivariate_terminal_thirteen_shared_light as light
from verify_multivariate_terminal_fourteen_high_curvature import CASES, ORDERS


CAP = Q(1, 8)
TAIL_EXPONENT = Q(21, 26)
TAIL_BASE = Q(2, 5)
ATTENUATION = 280
REFERENCE_RATIO = Q(169, 168)
PRIVATE_ONE = ((1,), ())
RETENTION_CASES = tuple((a, b, cutoff) for a, b, cutoff, *_ in CASES
                        if 0 < cutoff and sum(a) + sum(b) <= 3)
LIST_COUNTS = (278, 89, 1, 1, 1, 12, 10, 6, 1, 1)
# Outer endpoints and residual reciprocal cap for the empty one-private list.
EMPTY_BANDS = (
    (Q(0), Q(3, 10), CAP),
    (Q(3, 10), Q(2, 5), Q(1, 20)),
    (Q(2, 5), Q(1, 2), Q(1, 20)),
    (Q(1, 2), Q(3, 5), Q(1, 20)),
    (Q(3, 5), Q(3, 4), Q(1, 20)),
    (Q(3, 4), Q(1), CAP),
)


def verify_parameters():
    retained.ctx.prec = 192
    lo, hi = ORDERS
    _, log_outer_upper = retained.log_bounds(1 / TAIL_BASE)
    retained.require(378 * TAIL_EXPONENT * log_outer_upper < ATTENUATION,
                     "factorial tail attenuation")
    retained.require(retained.log_bounds(Q(281, 112))[0] / log_outer_upper
                     > 1 + Q(1, ATTENUATION), "Gaussian tail attenuation")
    retained.require(Q(1, lo * (lo - 1)) > Q(hi, lo * ATTENUATION),
                     "Gaussian endpoint exponent")


def skip_direct_disks(*_):
    return False


def configure_shared_model():
    """Bind the shared formulas before their first evaluation in a replay process."""
    cached = (light.ratio_bounds, light.independent_exponents, light.ordinary_ratio_upper,
              light.constants, light.leading_skew_direction, light.coefficient_model)
    retained.require(all(method.cache_info().currsize == 0 for method in cached),
                     "fresh order-dependent shared-model caches")
    light.LOWER, light.UPPER = ORDERS
    light.ORDERS = ORDERS
    light.REFERENCE_RATIO = REFERENCE_RATIO
    return light


def configure_replay():
    configure_shared_model()
    retained.require(retained.core_coefficients.cache_info().currsize == 0,
                     "fresh retained coefficient cache")
    retained.LOWER, retained.UPPER = ORDERS
    retained.ORDERS = ORDERS
    retained.CAP = CAP
    retained.SECTOR_ORDERS = (14, 13) * 6 + (14,)
    retained.CASES = RETENTION_CASES
    retained.TAIL_EXPONENT = TAIL_EXPONENT
    retained.TAIL_BASE = TAIL_BASE
    retained.ATTENUATION = ATTENUATION
    retained.verify_parameters = verify_parameters
    retained.verify_disks = skip_direct_disks
    verify_parameters()


def build_projector(directory):
    # The endpoint products are bounded by the original input guards. The
    # additional shift at 15 slots bounds the transformed entries by 2^82.
    direction_bound = (1 << 62) + (1 << 66)
    control_bound = 15 * direction_bound * (1 << (82 - 48)) + 1
    midpoint_bound = control_bound * (1 << 14)
    retained.require((2 * midpoint_bound + 1) * 256 < 1 << 127,
                     "15-slot signed accumulation bound")
    original = Path(__file__).with_name("verify_multivariate_terminal_ten_shared_light.cpp")
    source = Path(directory) / "bernstein.cpp"
    library = Path(directory) / "bernstein.so"
    body = original.read_text()
    replacements = (
        ("max_coefficients=14", "max_coefficients=15"),
        ("n==14 ? 83 : 84", "n==15 ? 82 : (n==14 ? 83 : 84)"),
        ("(n==13 ? 8 : 13)", "(n==13 ? 8 : (n==14 ? 13 : 17))"),
    )
    for before, after in replacements:
        retained.require(body.count(before) == 1, ("native source substitution", before))
        body = body.replace(before, after)
    source.write_text(body)
    run(("c++", "-std=c++17", "-O3", "-fPIC", "-shared", str(source), "-o", str(library)),
        check=True)
    return retained.scaled_native(CDLL(str(library)), compress=True)


def verify_empty_one(replay):
    from verify_multivariate_terminal_thirteen_retained_bands import bracket_tail, retained_cores

    total_outer = total_radial = list_count = 0
    strength_cutoff = RETENTION_CASES[1][2]
    for lower, upper, cap in EMPTY_BANDS:
        if not lower:
            lower = retained.powq(TAIL_BASE, TAIL_EXPONENT * strength_cutoff)[0]
        additional = retained_cores(lower, cap, CAP, canonical=False)
        outer = radial = 0
        for rows, squared in (((), Q(1)), *additional):
            residual = bracket_tail(lower, upper, squared) if rows else (lower, upper)
            l, r, d = retained.verify_interval(
                rows, PRIVATE_ONE, residual, cap, replay, model_builder=retained.refined_model,
                grouped_depth=0, maximum_grouped_width=Q(1, 200))
            retained.require(d == 0, "ray-only replay")
            outer += l
            radial += r
        lists = 1 + len(additional)
        list_count += lists
        total_outer += outer
        total_radial += radial
        print(f"Empty core band {lower} to {upper} at cap {cap}: {lists} lists, "
              f"{outer} outer and {radial} radial cells.", flush=True)
    print(f"Verified empty one-private core, {list_count} list occurrences, "
          f"{total_outer} outer and {total_radial} radial cells.", flush=True)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("all", "shared", "one-private", "other-private"),
                        nargs="?", default="all")
    args = parser.parse_args()
    configure_replay()
    for case, expected in zip(RETENTION_CASES, LIST_COUNTS, strict=True):
        found = len(retained.cores(case))
        retained.require(found == expected, ("retained list count", case, found, expected))
    with TemporaryDirectory(prefix="sfc-fourteen-retained-") as directory:
        replay = build_projector(directory)
        if args.case in ("all", "shared"):
            retained.verify((0,), None, replay)
        if args.case in ("all", "one-private"):
            retained.verify((1,), tuple(range(1, LIST_COUNTS[1])), replay)
            verify_empty_one(replay)
        if args.case in ("all", "other-private"):
            retained.verify(range(2, len(RETENTION_CASES)), None, replay)


if __name__ == "__main__":
    main()
