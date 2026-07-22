#!/usr/bin/env python3
"""Audit the single mathematical-status authority and its generated view."""

from render_status import ROOT, load_index, render, validate_index


index = load_index()
validate_index(index)
status_path = ROOT / "STATUS.md"
assert status_path.is_file(), "STATUS.md is missing"
assert status_path.read_text() == render(index), (
    "STATUS.md is stale; run python3 scripts/render_status.py"
)

for duplicate in ("THEOREMS.yml", "CLAIMS.md", "RESEARCH_STATUS.md"):
    assert not (ROOT / duplicate).exists(), f"{duplicate} duplicates MATH_STATUS.json"

readme = (ROOT / "README.md").read_text()
assert "MATH_STATUS.json" in readme, "README.md must link the status authority"
assert "THEOREMS.yml" not in readme, "README.md links the retired theorem index"

print("PASS mathematical status: one typed authority, valid graph, current rendering")
