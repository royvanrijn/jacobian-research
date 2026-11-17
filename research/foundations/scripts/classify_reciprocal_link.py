#!/usr/bin/env python3
"""Display exact reciprocal-link certificates for the built-in examples.

For new candidates, import ``ReciprocalLinkCandidate`` and
``classify_reciprocal_link`` from ``jcsearch.reciprocal``; the example below
is also the smallest complete input template.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from jcsearch.reciprocal import (
    classify_boundary_reconstruction,
    classify_reciprocal_link,
    masuda_hidden_cover_example,
    masuda_plinth_example,
    standard_cancellation_example,
)


def stringify(value):
    if isinstance(value, sp.Basic):
        return str(value)
    if isinstance(value, tuple):
        return [stringify(item) for item in value]
    if isinstance(value, list):
        return [stringify(item) for item in value]
    if isinstance(value, dict):
        return {str(key): stringify(item) for key, item in value.items()}
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "example",
        choices=("cancellation", "masuda", "masuda-hidden"),
        help="built-in exact candidate to classify",
    )
    args = parser.parse_args()

    if args.example == "cancellation":
        certificate = classify_reciprocal_link(standard_cancellation_example())
    elif args.example == "masuda":
        certificate = classify_boundary_reconstruction(masuda_plinth_example())
    else:
        certificate = classify_boundary_reconstruction(masuda_hidden_cover_example())
    print(json.dumps(stringify(asdict(certificate)), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
