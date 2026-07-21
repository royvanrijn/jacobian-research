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

expected = [f"C{i:02d}" for i in range(1, 17)] + [
    "Cancellation construction",
    "Cancellation-parameter arithmetic",
    "Boundary distinction",
    "Rigidity within the current ansatz",
]
actual = [cells[0].split(" — ", 1)[0] for cells in rows]
assert actual == expected, "STATUS.md result rows are missing, duplicated, or out of order"

for result, mathematical, review, location in rows:
    assert mathematical, f"{result} has no mathematical status"
    assert review, f"{result} has no external-review entry"
    assert "](" in location, f"{result} has no linked location"

readme = (root / "README.md").read_text()
assert "C01–C04 are the stable core. See [STATUS.md](STATUS.md)" in readme
assert "| Claim | Description | Status |" not in readme
assert "CLAIMS.md" not in readme and "RESEARCH_STATUS.md" not in readme

print("PASS status ledger: one four-column source covers C01--C16 and cancellation")
