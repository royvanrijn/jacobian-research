#!/usr/bin/env python3
"""Reject unresolved references, citations, and duplicate labels in LaTeX logs."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

_FATAL_PATTERNS: tuple[tuple[str, re.Pattern[str]], ...] = (
    ("undefined reference", re.compile(r"LaTeX Warning: Reference .* undefined")),
    ("undefined citation", re.compile(r"LaTeX Warning: Citation .* undefined")),
    ("undefined references", re.compile(r"LaTeX Warning: There were undefined references")),
    ("undefined citations", re.compile(r"LaTeX Warning: There were undefined citations")),
    ("multiply-defined label", re.compile(r"LaTeX Warning: Label .* multiply defined")),
    ("multiply-defined labels", re.compile(r"LaTeX Warning: There were multiply-defined labels")),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", nargs="+", type=Path, help="LaTeX .log files to inspect")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    failures: list[str] = []

    for path in args.logs:
        if not path.is_file():
            failures.append(f"{path}: log file is missing")
            continue

        for line_number, line in enumerate(path.read_text(errors="replace").splitlines(), start=1):
            for label, pattern in _FATAL_PATTERNS:
                if pattern.search(line):
                    failures.append(f"{path}:{line_number}: {label}: {line.strip()}")

    if failures:
        print("LaTeX log audit failed:", file=sys.stderr)
        for failure in failures:
            print(f"  {failure}", file=sys.stderr)
        return 1

    print("LaTeX log audit passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
