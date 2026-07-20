#!/usr/bin/env python3
"""Fail if a relative Markdown link in README.md does not exist."""
from pathlib import Path
import re

root = Path(__file__).resolve().parents[1]
readme = (root / "README.md").read_text()
links = re.findall(r"\[[^]]+\]\(([^)]+)\)", readme)
missing = [target for target in links
           if "://" not in target and not (root / target).exists()]
if missing:
    raise SystemExit("missing README links: " + ", ".join(missing))
print(f"PASS README relative links ({len(links)} checked)")
