#!/usr/bin/env python3
"""Independently replay the frozen extraction with bundled openpyxl, read-only."""

import ast
import hashlib
import json
from pathlib import Path
from openpyxl import load_workbook

HERE = Path(__file__).resolve().parent


def main():
    panel = json.loads((HERE / "lift_panel.json").read_text())
    path = HERE / panel["source"]
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest == panel["sha256"]
    workbook = load_workbook(path, read_only=True, data_only=False, keep_links=False)
    sheet = workbook.worksheets[0]
    assert sheet.title == panel["sheet"]
    found = []
    for number, row in enumerate(sheet.iter_rows(min_col=1, max_col=10, values_only=True), 1):
        if number == 1:
            assert list(row) == [panel["headers"][column] for column in "ABCDEFGHIJ"]
        elif str(row[2]) == "5" and str(row[3]) == "0":
            found.append({"name": row[0], "gauss_code": ast.literal_eval(row[1]),
                          "source_row": number, "source_range": f"A{number}:J{number}",
                          "reported_wirtinger_number": 5,
                          "reported_maximal_rank_finite_coxeter_hit": False})
            if len(found) == 16:
                break
    workbook.close()
    assert found == panel["cases"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    report = {"status": "PASS", "independent_reader": "openpyxl read_only=True; extraction producer uses standard-library ZIP/XML",
              "cases": len(found), "source_sha256": digest, "workbook_unchanged": True,
              "panel_sha256": hashlib.sha256((HERE / "lift_panel.json").read_bytes()).hexdigest(),
              "reader_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE / "panel_source_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
