#!/usr/bin/env python3
"""Build the recovered Fermigier score/search corpus before model fitting."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
ECSEARCH = ROOT / "elliptic-curves/ecsearch"
if str(ECSEARCH) not in sys.path:
    sys.path.insert(0, str(ECSEARCH))

from fermigier_labelled_corpus import build_corpus  # noqa: E402


DEFAULT_OUTPUT = (
    ROOT / "artifacts/local/elliptic-curves/fermigier-labelled-corpus-v1.jsonl.gz"
)
DEFAULT_SUMMARY = (
    ROOT / "artifacts/local/elliptic-curves/fermigier-labelled-corpus-v1-summary.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--split-salt", default="fermigier-labelled-corpus-v1")
    args = parser.parse_args()
    summary = build_corpus(
        ROOT,
        output=args.output,
        summary_output=args.summary,
        split_salt=args.split_salt,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

