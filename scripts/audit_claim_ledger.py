#!/usr/bin/env python3
"""Audit the concise four-status claim ledger."""

from pathlib import Path
import re


root = Path(__file__).resolve().parents[1]
ledger = (root / "CLAIMS.md").read_text()

allowed = {
    "Verified independently",
    "Proved internally",
    "Computational evidence",
    "Conjectural / incomplete",
}

rows = []
for line in ledger.splitlines():
    match = re.match(r"^\| C(\d{2}) \|", line)
    if not match:
        continue
    cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
    assert len(cells) == 5, f"C{match.group(1)} does not have five columns"
    rows.append((int(match.group(1)), cells[2]))

active_claims = list(range(1, 17)) + [24]
assert [number for number, _ in rows] == active_claims, (
    "CLAIMS.md must contain active claims C01--C16 and C24 in order; "
    "the archived transfer programme must not reappear as claim rows"
)
for number, status in rows:
    assert status in allowed, f"C{number:02d} has unsupported status {status!r}"

research_status = (root / "RESEARCH_STATUS.md").read_text()
assert re.search(r"^\| C24 \|.*\| \*\*Active\*\* \|", research_status, re.MULTILINE)
assert "Proved construction with incomplete global classification" in research_status
print("PASS claim ledger: active C01--C16 and C24 use approved statuses; transfer rows absent")
