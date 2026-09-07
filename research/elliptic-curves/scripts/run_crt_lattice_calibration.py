#!/usr/bin/env python3
"""Generate the fixed CRT--lattice calibration manifest."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.calibration import build_calibration, dump_calibration  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        help="fresh JSON output path (existing files are never overwritten)",
    )
    parser.add_argument("--maximum-height", type=int, default=262144)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = build_calibration(maximum_height=args.maximum_height)
    best = result["best_candidate"]
    print(
        "CALIBRATION "
        f"t={best['numerator']}/{best['denominator']} "
        f"M={best['crt_modulus']} "
        f"rank_bounds={best['pari_gp']['rank_bounds']} "
        f"conductor={best['pari_gp']['conductor']}"
    )
    if args.output is not None:
        dump_calibration(result, args.output)
        print(f"WROTE {args.output}")


if __name__ == "__main__":
    main()
