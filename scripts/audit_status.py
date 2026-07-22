#!/usr/bin/env python3
"""Audit THEOREMS.yml and its generated STATUS.md summary."""

from render_status import ROOT, load_index, render, validate_index


index = load_index()
validate_index(index)

status_path = ROOT / "STATUS.md"
assert status_path.is_file(), "STATUS.md is missing"
assert status_path.read_text() == render(index), (
    "STATUS.md is stale; run python3 scripts/render_status.py"
)

for duplicate in ("CLAIMS.md", "RESEARCH_STATUS.md"):
    assert not (ROOT / duplicate).exists(), f"{duplicate} duplicates THEOREMS.yml"

readme = (ROOT / "README.md").read_text()
assert "THEOREMS.yml" in readme, "README.md must link the theorem source of truth"
assert "CLAIMS.md" not in readme and "RESEARCH_STATUS.md" not in readme

print("PASS theorem index: dependency graph is valid and STATUS.md is current")
