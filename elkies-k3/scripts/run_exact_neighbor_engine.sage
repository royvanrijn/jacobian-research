#!/usr/bin/env sage
"""Run the exact neighbor engine from a versioned JSON input."""

import argparse
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
load(str(HERE / "exact_neighbor_engine.sage"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    data = json.loads(arguments.input.read_text())
    result, certificate = run_neighbor_input(data)
    digest = write_neighbor_certificate(arguments.output, certificate)
    print(
        "EXACTNEIGHBOR|input={}|output={}|degree={}|root_data={}|"
        "minimization={}|sha256={}|status=PASS".format(
            arguments.input,
            arguments.output,
            certificate["divisor"]["old_fiber_degree"],
            tuple(certificate["child"]["root_data"]),
            result["minimization_status"],
            digest,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
