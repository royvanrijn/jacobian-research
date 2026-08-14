#!/usr/bin/env python3
"""Generate the first CRT seed in the Fermigier adapter family."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.fermigier_seed import (  # noqa: E402
    build_fermigier_seed,
    dump_fermigier_seed,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--maximum-height", type=int, default=2**23)
    args = parser.parse_args()
    result = build_fermigier_seed(maximum_height=args.maximum_height)
    best = result["best_seed"]
    print(
        "FERMIGIER_SEED "
        f"u={best['numerator']}/{best['denominator']} "
        f"M={best['crt_modulus']} height={best['height']}"
    )
    if args.output is not None:
        dump_fermigier_seed(result, args.output)
        print(f"WROTE {args.output}")


if __name__ == "__main__":
    main()
