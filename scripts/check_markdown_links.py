#!/usr/bin/env python3
"""Fail if a local Markdown document reference does not exist."""

from pathlib import Path
import re


root = Path(__file__).resolve().parents[1]
documents = [root / "README.md", root / "results" / "README.md"]
documents.extend(sorted((root / "notes").glob("*.md")))

missing = []
checked = 0
for document in documents:
    text = document.read_text()
    targets = re.findall(r"\[[^]]+\]\(([^)]+)\)", text)
    for target in targets:
        if "://" in target or target.startswith(("#", "mailto:")):
            continue
        # Avoid interpreting TeX such as ``[X_i,X_j](F_k)`` as a file link.
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
