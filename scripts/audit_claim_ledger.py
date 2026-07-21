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

assert [number for number, _ in rows] == list(range(1, 25)), (
    "CLAIMS.md must contain exactly C01--C24 in order"
)
for number, status in rows:
    assert status in allowed, f"C{number:02d} has unsupported status {status!r}"

assert "C24 | **Active**" in (root / "RESEARCH_STATUS.md").read_text()
print("PASS claim ledger: C01--C24 use only the four approved statuses")
