#!/usr/bin/env python3
"""Fail if a source-document reference is broken.

Links into the generated/local artifact trees are reproducibility declarations:
some outputs are intentionally absent from a lightweight checkout.  Validate
their paths when present and report absent outputs without treating them as
broken source documentation.
"""

from pathlib import Path
import re


root = Path(__file__).resolve().parents[1]
excluded_parts = {".git", ".venv", ".cache", ".idea", ".lake"}
documents = sorted(
    path for path in root.rglob("*.md")
    if not excluded_parts.intersection(path.relative_to(root).parts)
)

missing = []
absent_generated_outputs = []
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
            try:
                relative_candidate = candidate.resolve().relative_to(root.resolve())
            except ValueError:
                relative_candidate = None
            if relative_candidate is not None and (
                relative_candidate.parts[:2] == ("artifacts", "generated-results")
                or relative_candidate.parts[:2] == ("artifacts", "local")
            ):
                absent_generated_outputs.append(
                    f"{document.relative_to(root)} -> {target}"
                )
                continue
            missing.append(f"{document.relative_to(root)} -> {target}")

if missing:
    raise SystemExit("missing Markdown references:\n" + "\n".join(missing))

print(
    f"PASS local Markdown references ({checked} checked in {len(documents)} files; "
    f"{len(absent_generated_outputs)} documented generated outputs absent locally)"
)
