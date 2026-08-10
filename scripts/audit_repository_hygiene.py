#!/usr/bin/env python3
"""Audit active-note identities and repository-only generated debris."""

from __future__ import annotations

import re
import subprocess
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXCLUDED_PARTS = {
    ".cache",
    ".git",
    ".idea",
    ".lake",
    ".venv",
    "archive",
    "output",
    "tmp",
}
DEFINITION_PATTERNS = (
    re.compile(
        r"\*\*(?:Theorem|Corollary|Proposition|Lemma)\s+`?"
        r"([A-Z][A-Z0-9-]*\d[A-Z0-9-]*)"
    ),
    re.compile(r"\*\*`?([A-Z][A-Z0-9-]*\d[A-Z0-9-]*)`?\.\*\*"),
)


definitions: dict[str, set[Path]] = defaultdict(set)
for document in sorted(ROOT.rglob("*.md")):
    relative = document.relative_to(ROOT)
    if EXCLUDED_PARTS.intersection(relative.parts):
        continue
    for line in document.read_text().splitlines():
        for pattern in DEFINITION_PATTERNS:
            match = pattern.search(line)
            if match:
                definitions[match.group(1)].add(relative)
                break

duplicates = {
    identifier: paths
    for identifier, paths in definitions.items()
    if len(paths) > 1
}
if duplicates:
    details = []
    for identifier, paths in sorted(duplicates.items()):
        details.append(f"{identifier}: {', '.join(map(str, sorted(paths)))}")
    raise SystemExit(
        "headline theorem identifiers are defined in multiple active notes:\n"
        + "\n".join(details)
    )

tracked_ignored: list[str] = []
if (ROOT / ".git").exists():
    result = subprocess.run(
        ["git", "ls-files", "-ci", "--exclude-standard"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    tracked_ignored = result.stdout.splitlines()
if tracked_ignored:
    raise SystemExit(
        "ignored generated files are still tracked:\n" + "\n".join(tracked_ignored)
    )

print(
    f"PASS repository hygiene: {len(definitions)} active headline identifiers are "
    "file-unique and no ignored generated files are tracked"
)
