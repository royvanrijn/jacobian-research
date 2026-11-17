#!/usr/bin/env python3
"""Freeze first sixteen currently listed uncertain cases, before lift testing."""

import ast
import hashlib
import json
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def main():
    source = HERE / "sources/current13/data/still_open.json"
    names = json.loads(source.read_text())
    assert len(names) == len(set(names)) == 293
    wanted = names[:16]
    normalized = {name.replace("_", ""): name for name in wanted}
    workbook = HERE / "sources/wirt-hm/all_data_A.xlsx"
    rows = {}
    with zipfile.ZipFile(workbook) as z:
        strings = []
        for _, e in ET.iterparse(z.open("xl/sharedStrings.xml"), events=["end"]):
            if e.tag == NS + "si":
                strings.append("".join(t.text or "" for t in e.iter(NS + "t")))
                e.clear()
        for _, e in ET.iterparse(z.open("xl/worksheets/sheet1.xml"), events=["end"]):
            if e.tag != NS + "row":
                continue
            row = {}
            for c in e.findall(NS + "c"):
                v = c.find(NS + "v")
                if v is not None:
                    assert c.find(NS + "f") is None and c.get("t") != "e"
                    row["".join(x for x in c.get("r") if x.isalpha())] = (
                        strings[int(v.text)] if c.get("t") == "s" else v.text)
            if row.get("A") in normalized:
                name = normalized[row["A"]]
                assert name not in rows
                rows[name] = {"name": name, "workbook_name": row["A"],
                              "source_row": int(e.get("r")), "gauss_code": ast.literal_eval(row["B"]),
                              "reported_wirtinger_number": int(row["C"]),
                              "reported_maximal_finite_coxeter_hit": bool(int(row["D"]))}
                if len(rows) == 16:
                    break
            e.clear()
    assert set(rows) == set(wanted)
    manifest = json.loads((HERE / "sources/current13/manifest.json").read_text())
    record = {"selection": "First 16 names in pinned still_open.json, no adaptive replacement; join by removing only underscore against pinned BKM workbook",
              "current13_commit": manifest["commit"],
              "sources": {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in [source, workbook]},
              "search_limits": {"degrees": [3, 4], "seed_count": 4,
                                "max_seed_trials_per_case": 1000, "max_seconds_per_case": 20,
                                "coefficient_fields": "All characteristics, exactly by integral Smith form",
                                "no_adaptive_replacements": True},
              "claim_boundary": "Listed uncertain in this pinned release, not an exhaustive current literature audit. No new case until independent replay and prior-resolution audit.",
              "cases": [rows[name] for name in wanted]}
    out = HERE / "current13_panel.json"
    if out.exists():
        assert json.loads(out.read_text()) == record
    else:
        out.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({"status": "FROZEN", "cases": wanted,
                      "reported_upper_bounds": [x["reported_wirtinger_number"] for x in record["cases"]]}))


if __name__ == "__main__":
    main()
