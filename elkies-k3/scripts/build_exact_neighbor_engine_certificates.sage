#!/usr/bin/env sage
"""Build or verify the pinned q80 and H3 exact-neighbor certificates."""

import argparse
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
load(str(HERE / "exact_neighbor_examples.sage"))


DEFAULT_DIRECTORY = ROOT / "artifacts/generated-results"


def write_or_verify(path, certificate):
    canonical = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if path.exists():
        if path.read_text() != canonical:
            raise SystemExit(
                f"refusing to replace pinned certificate with different content: {path}"
            )
        return "verified"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(canonical)
    return "written"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--q80-output",
        type=Path,
        default=DEFAULT_DIRECTORY / "elkies-k3-exact-neighbor-q80-first-q4.json",
    )
    parser.add_argument(
        "--h3-output",
        type=Path,
        default=DEFAULT_DIRECTORY / "elkies-k3-exact-neighbor-h3-d13-q24.json",
    )
    arguments = parser.parse_args()

    for factory, output in (
        (q80_first_q4_example, arguments.q80_output),
        (h3_d13_q24_example, arguments.h3_output),
    ):
        example = factory()
        _, certificate = run_example(example)
        action = write_or_verify(output, certificate)
        print(
            "EXACTNEIGHBOR|example={}|artifact={}|sha256={}|action={}|status=PASS".format(
                example["name"], output, certificate["certificate_sha256"], action
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
