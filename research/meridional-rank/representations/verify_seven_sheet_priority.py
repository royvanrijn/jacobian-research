#!/usr/bin/env python3
"""Check the exact two-source comparison for 24 seven-sheet equalities.

Default execution checks retained snapshots, not the external workbook.
--workbook PATH additionally hashes and re-extracts the pinned original XLSX
read-only.  --output records the actual scope of this invocation; without it
no existing verification receipt is overwritten.
"""

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def git_blob(data):
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()


def positive_names(inputs):
    names = [
        case["name"] for case in inputs["cases"]
        if any(action["expected_H1_F3_dimension"] == 1 for action in case["actions"])
    ]
    assert len(names) == len(set(names)) == 24
    return names


def check_snapshot(audit, inputs, still_open_bytes):
    names = positive_names(inputs)
    source = audit["still_open_source"]
    assert sha256(still_open_bytes) == source["sha256"], "still-open source hash differs"
    assert git_blob(still_open_bytes) == source["git_blob"], "still-open Git blob differs"
    still_open = json.loads(still_open_bytes)
    assert len(still_open) == len(set(still_open)) == 293
    assert set(names).issubset(still_open), "positive name absent from retained still-open list"
    rows = audit["rows"]
    assert [row["name"] for row in rows] == names, "retained rows do not match all 24 positives"
    assert len({row["row"] for row in rows}) == 24
    assert audit["upstream_current_blob"] == audit["downloaded_git_blob"]
    assert audit["upstream_current_blob"] == "3882c6da34c9057f8f13da31afb9d8973b7158cb"
    assert audit["bytes"] == 23554211
    assert audit["sha256"] == "514b12e118fa681520da5702a5ed4307ce7cf93b6b9f7b5e53f265ca7b01647b"
    assert set(audit["column_headers"]) == set("ABCDEFGHIJ")
    for row in rows:
        number = str(row["row"])
        cells = row["cached_cells"]
        assert set(cells) == {column + number for column in "ABCDEFGHIJ"}
        assert cells["A" + number] == row["name"].replace("_", ""), "workbook row name differs"
        assert cells["C" + number] == "4", "recorded Wirtinger number is not four"
        assert all(cells[column + number] == "0" for column in "DEGI"), "a maximal Coxeter hit is recorded"
        assert all(cells[column + number] == "N/A" for column in "FHJ")
        assert not row["formula_cells"], "a retained comparison cell is formula-driven"
    return names


def reextract_workbook(path, audit):
    """Read source ZIP members directly, preserving all cached cell text."""
    data = path.read_bytes()
    assert len(data) == audit["bytes"], "workbook byte length differs"
    assert sha256(data) == audit["sha256"], "workbook SHA256 differs"
    assert git_blob(data) == audit["downloaded_git_blob"], "workbook Git blob differs"
    wanted = {row["row"]: row for row in audit["rows"]}
    found = {}
    with ZipFile(path) as archive:
        shared = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        strings = ["".join(item.itertext()) for item in shared]
        for _, row in ET.iterparse(archive.open("xl/worksheets/sheet1.xml"), events=("end",)):
            if row.tag.rsplit("}", 1)[-1] != "row":
                continue
            number = int(row.attrib["r"])
            if number == 1 or number in wanted:
                cells, formulas = {}, []
                for cell in row:
                    value = cell.find("s:v", NS)
                    if value is not None:
                        cells[cell.attrib["r"]] = (
                            strings[int(value.text)] if cell.attrib.get("t") == "s" else value.text
                        )
                    if cell.find("s:f", NS) is not None:
                        formulas.append(cell.attrib["r"])
                if number == 1:
                    assert cells == {column + "1": text for column, text in audit["column_headers"].items()}
                else:
                    retained = wanted[number]
                    assert cells == retained["cached_cells"], "re-extracted workbook cells differ"
                    assert formulas == retained["formula_cells"], "re-extracted formulas differ"
                    found[number] = True
            row.clear()
            if len(found) == len(wanted):
                break
    assert set(found) == set(wanted), "workbook source rows are incomplete"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, help="optional original pinned all_data_A.xlsx")
    parser.add_argument("--output", type=Path, help="write an explicit receipt for this invocation")
    args = parser.parse_args()
    audit_path = HERE / "seven_sheet_priority_audit.json"
    audit = json.loads(audit_path.read_text())
    input_path = HERE / audit["positive_inputs"]
    inputs = json.loads(input_path.read_text())
    source_path = HERE / audit["still_open_source"]["path"]
    source_bytes = source_path.read_bytes()
    names = check_snapshot(audit, inputs, source_bytes)
    rejected = []
    bad = deepcopy(audit)
    first = bad["rows"][0]
    first["cached_cells"]["D" + str(first["row"])] = "1"
    try:
        check_snapshot(bad, inputs, source_bytes)
    except AssertionError as error:
        assert "Coxeter hit" in str(error)
        rejected.append("changed Coxeter hit rejected")
    else:
        raise AssertionError("corrupt source hit accepted")
    bad = deepcopy(audit)
    bad["rows"][0]["name"] = "13a_1"
    try:
        check_snapshot(bad, inputs, source_bytes)
    except AssertionError as error:
        assert "all 24 positives" in str(error)
        rejected.append("changed comparison name rejected")
    else:
        raise AssertionError("corrupt comparison name accepted")
    if args.workbook:
        reextract_workbook(args.workbook, audit)
    paths = [Path(__file__).resolve(), audit_path, input_path, source_path]
    report = {
        "status": "PASS",
        "mode": "pinned workbook re-extraction and retained snapshot checks" if args.workbook else "retained snapshot consistency only",
        "workbook_reextracted": bool(args.workbook),
        "positive_names_checked": names,
        "positive_count": len(names),
        "still_open_source_count": 293,
        "zero_recorded_maximal_coxeter_hits": 24,
        "negative_controls": rejected,
        "python_version": sys.version.split()[0],
        "scope": audit["scope"],
        "input_sha256": {path.relative_to(HERE).as_posix(): sha256(path.read_bytes()) for path in paths},
    }
    if args.workbook:
        report["workbook_sha256"] = audit["sha256"]
        report["workbook_git_blob"] = audit["downloaded_git_blob"]
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
