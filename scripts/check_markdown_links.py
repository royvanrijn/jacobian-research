#!/usr/bin/env python3
"""Fail if a local Markdown document reference does not exist."""

from pathlib import Path
import re


root = Path(__file__).resolve().parents[1]
excluded_parts = {".git", ".venv", ".cache", ".idea"}
documents = sorted(
    path for path in root.rglob("*.md")
    if not excluded_parts.intersection(path.relative_to(root).parts)
)

missing = []
checked = 0
for document in documents:
    text = document.read_text()
    targets = re.findall(r"\[[^]]+\]\(([^)]+)\)", text)
    for target in targets:
        if "://" in target or target.startswith(("#", "mailto:")):
            continue
        if "." not in target and "/" not in target:
            continue
        path_text = target.split("#", 1)[0]
        candidate = document.parent / path_text
        checked += 1
        if not candidate.exists():
            missing.append(f"{document.relative_to(root)} -> {target}")

if missing:
    raise SystemExit("missing Markdown references:\n" + "\n".join(missing))

print(f"PASS local Markdown references ({checked} checked in {len(documents)} files)")
