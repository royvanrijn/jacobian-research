#!/usr/bin/env python3
"""Audit STATUS.md as the repository's sole mathematical-status ledger."""

from pathlib import Path


root = Path(__file__).resolve().parents[1]
status_path = root / "STATUS.md"
assert status_path.is_file(), "STATUS.md is missing"
assert not (root / "CLAIMS.md").exists(), "CLAIMS.md duplicates STATUS.md"
assert not (root / "RESEARCH_STATUS.md").exists(), (
    "RESEARCH_STATUS.md duplicates STATUS.md"
)

status = status_path.read_text()
header = "| Result | Mathematical status | External review | Location |"
assert status.count(header) == 1, "STATUS.md must contain the canonical four-column table"

rows = []
for line in status.splitlines():
    if not line.startswith("|") or line == header or line.startswith("|---"):
        continue
    cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
    assert len(cells) == 4, f"status row does not have four columns: {line}"
    rows.append(cells)

expected = [
    "Foundational Keller map and rational collision",
    "Exact symplectic cotangent lift and Weyl quantization",
    "Cubic marked-root realization",
    "Exact cubic image, fibers, and nonproperness",
    "Weighted marked-root family and `S_n` monodromy",
    "Generic discriminant geometry",
    "Weighted image and boundary theorems",
    "Full-contact omission and uniqueness",
    "Contact strata and dimensions",
    "Contact-atom principle",
    "Exceptional components and closure order",
    "Component normalizations",
    "Degree-twelve local singularity",
    "Effective finite-field Chebotarev law",
    "Explicit quartic weighted model",
    "External quartic-island classification",
    "Stable normal-form consequences",
    "Dicritical compactification",
    "Marked-point dimension barrier",
    "Cancellation construction",
    "Cancellation-parameter arithmetic",
    "Boundary distinction",
    "Rigidity within the current ansatz",
    "Degreewise stable-multiplicity theorem",
    "Degree-five stable-moduli theorem",
]
actual = [cells[0] for cells in rows]
assert actual == expected, "STATUS.md result rows are missing, duplicated, or out of order"

canonical_degreewise = "[Canonical paper](papers/marked-root-multiplicity/main.tex)"
assert status.count(canonical_degreewise) == 1, (
    "the degreewise theorem must have exactly one canonical source"
)

for result, mathematical, review, location in rows:
    assert mathematical, f"{result} has no mathematical status"
    assert review, f"{result} has no external-review entry"
    assert "](" in location, f"{result} has no linked location"

readme = (root / "README.md").read_text()
assert "The foundational map, cubic marked-root model, exact image theorem" in readme
assert "The canonical statement and proof are in the standalone paper" in readme
assert "| Claim | Description | Status |" not in readme
assert "CLAIMS.md" not in readme and "RESEARCH_STATUS.md" not in readme

print("PASS status ledger: one four-column source covers the named active results")
