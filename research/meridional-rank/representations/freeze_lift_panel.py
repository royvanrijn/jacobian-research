#!/usr/bin/env python3
"""Read the fixed author workbook; freeze cases without running a knot search.

Use the bundled Python runtime. Only standard-library ZIP/XML readers are
used; the workbook is neither modified nor exported.
"""

import ast
import hashlib
import itertools
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

HERE = Path(__file__).resolve().parent
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def main():
    path = HERE / "sources/wirt-hm/all_data_A.xlsx"
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = json.loads((path.parent / "manifest.json").read_text())
    assert next(x for x in manifest["files"] if x["path"] == path.name)["sha256"] == digest
    selected, headers, rows_seen = [], None, 0
    with zipfile.ZipFile(path) as archive:
        strings = []
        for event, element in ET.iterparse(archive.open("xl/sharedStrings.xml"), events=["end"]):
            if element.tag == NS + "si":
                strings.append("".join(e.text or "" for e in element.iter(NS + "t")))
                element.clear()
        for event, element in ET.iterparse(archive.open("xl/worksheets/sheet1.xml"), events=["end"]):
            if element.tag != NS + "row":
                continue
            row = {}
            for cell in element.findall(NS + "c"):
                assert cell.find(NS + "f") is None, "formula needs separate interpretation"
                value = cell.find(NS + "v")
                if value is None:
                    continue
                assert cell.get("t") != "e", "source cell has an Excel error"
                row["".join(x for x in cell.get("r") if x.isalpha())] = (
                    strings[int(value.text)] if cell.get("t") == "s" else value.text)
            number = int(element.get("r"))
            if number == 1:
                headers = row
                assert headers["C"] == "Wirtinger Number"
                assert headers["D"] == "Any Homomorphism? (1=Y,0=N)"
            else:
                rows_seen += 1
                if row.get("C") == "5" and row.get("D") == "0":
                    code = ast.literal_eval(row["B"])
                    assert len(code) % 2 == 0 and all(isinstance(v, int) and v for v in code)
                    selected.append({"name": row["A"], "gauss_code": code,
                                     "source_row": number, "source_range": f"A{number}:J{number}",
                                     "reported_wirtinger_number": 5,
                                     "reported_maximal_rank_finite_coxeter_hit": False})
                    if len(selected) == 16:
                        break
            element.clear()
    assert len(selected) == 16
    receipt = {"source": "sources/wirt-hm/all_data_A.xlsx", "sha256": digest,
               "commit": manifest["commit"], "sheet": "Sheet1", "headers": headers,
               "selection": "First 16 source-order rows with Wirtinger Number=5 and Any Homomorphism=0; frozen before testing",
               "rows_read_before_panel_complete": rows_seen,
               "search_limits": {"degrees": [3, 4], "seed_count": 5,
                                 "max_seed_trials_per_case": 5000, "max_seconds_per_case": 15,
                                 "no_adaptive_replacements": True}, "cases": selected}
    target = HERE / "lift_panel.json"
    if target.exists():
        assert json.loads(target.read_text()) == receipt, "frozen panel changed"
    else:
        target.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": "FROZEN", "cases": [x["name"] for x in selected],
                      "source_rows": [x["source_row"] for x in selected]}, indent=2))


if __name__ == "__main__":
    main()
