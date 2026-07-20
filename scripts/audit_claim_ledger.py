"""Audit the one-to-one vocabulary connecting README claims to CLAIMS.md."""

from pathlib import Path
import re

root = Path(__file__).resolve().parents[1]
ledger = (root / "CLAIMS.md").read_text()
readme = (root / "README.md").read_text()

ledger_ids = re.findall(r'<a id="(c\d{2})"></a>C\d{2}', ledger)
assert ledger_ids, "CLAIMS.md contains no claim rows"
assert len(ledger_ids) == len(set(ledger_ids)), "duplicate claim ID in CLAIMS.md"

readme_links = re.findall(r'\[Claim C(\d{2})\]\(CLAIMS\.md#(c\d{2})\)', readme)
assert readme_links, "README contains no claim-ledger links"
for label, anchor in readme_links:
    assert anchor == f"c{label}", f"claim label C{label} points to #{anchor}"
    assert anchor in ledger_ids, f"README points to missing ledger row {anchor}"

linked_ids = {anchor for _, anchor in readme_links}
unlinked = sorted(set(ledger_ids)-linked_ids)
assert not unlinked, f"ledger rows not linked from README: {unlinked}"
assert "[Claims " not in readme, "one README link must name exactly one claim row"

table_rows = [line for line in ledger.splitlines()
              if re.match(r"^\| <a id=", line)]
assert len(table_rows) == len(ledger_ids)
for row in table_rows:
    assert row.count("|") == 9, "claim row does not have exactly eight columns"

print(f"PASS claim ledger: {len(ledger_ids)} unique rows, all linked from README")
