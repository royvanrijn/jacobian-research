#!/usr/bin/env python3
"""Generate the eight-stage global low-degree census ledgers."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from jcsearch.low_degree_pipeline import (  # noqa: E402
    MANIFEST_FILENAME,
    build_low_degree_census,
)


DEFAULT_OUTPUT = ROOT / "artifacts/generated-results"


def serialized(payload: object) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--max-degree", type=int, default=7)
    parser.add_argument("--target-profile", default="7,6,4")
    parser.add_argument("--max-nonlinear-support", type=int, default=6)
    parser.add_argument("--primes", default="11,13,17")
    args = parser.parse_args()

    target = tuple(int(value) for value in args.target_profile.split(","))
    primes = tuple(int(value) for value in args.primes.split(","))
    artifacts = build_low_degree_census(
        max_degree=args.max_degree,
        target_profile=target,
        support_bound=args.max_nonlinear_support,
        primes=primes,
        progress=lambda message: print(message, flush=True),
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for filename, payload in artifacts.items():
        content = serialized(payload)
        path = args.output_dir / filename
        path.write_text(content)
        hashes[filename] = hashlib.sha256(content.encode()).hexdigest()
        print(f"wrote {path}")

    manifest = {
        "schema": "global-low-degree-census.manifest.v1",
        "generator": "scripts/compile_global_low_degree_census.py",
        "parameters": {
            "max_degree": args.max_degree,
            "target_profile": list(target),
            "max_nonlinear_support": args.max_nonlinear_support,
            "primes": list(primes),
        },
        "stage_sha256": hashes,
        "claim_boundary": (
            "Complete for invariant profile enumeration, the dense quadratic collision "
            f"ideal, and exact supports of size at most {args.max_nonlinear_support} "
            f"inside raw coordinate degree {args.max_degree}. The cardinality-unbounded "
            "support census remains open."
        ),
    }
    manifest_content = serialized(manifest)
    manifest_path = args.output_dir / MANIFEST_FILENAME
    manifest_path.write_text(manifest_content)
    print(f"wrote {manifest_path}")


if __name__ == "__main__":
    main()
