#!/usr/bin/env python3
"""Generate the exact rank-at-least-29 calibration certificate."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.record_rank29 import build_rank29_manifest  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--maximum-reduction-prime", type=int, default=600)
    args = parser.parse_args()
    manifest = build_rank29_manifest(
        maximum_reduction_prime=args.maximum_reduction_prime
    )
    certificate = manifest["independence_certificate"]
    print(
        "E29_INDEPENDENCE "
        f"relation_prime={certificate['relation_prime']} "
        f"rows={len(certificate['rows'])}"
    )
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
            handle.write("\n")
        print(f"WROTE {args.output}")


if __name__ == "__main__":
    main()
