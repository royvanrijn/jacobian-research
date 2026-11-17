#!/usr/bin/env python3
"""Seal the new primary group-list input and precisely delimit its use."""

import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent


def main():
    url = "https://brauer.maths.qmul.ac.uk/Atlas/v3/alt/A8/"
    snapshot = HERE / "sources/atlas-a8.html"
    if not snapshot.exists():
        snapshot.write_bytes(urlopen(url, timeout=30).read())
    text = snapshot.read_text()
    assert "20160" in text and "L" in text and "GF(2)" in text
    record = {"id": "ATLASA8", "title": "ATLAS: Alternating group A8, Linear group L4(2)",
        "authors": "ATLAS of Finite Group Representations", "year": 2026, "url": url,
        "snapshot": "sources/atlas-a8.html", "location": "Group order, standard-generator heading and characteristic-two representation list",
        "audit": "Primary group order20160 and identification A8=L4(2) read, together with the two dimension4 characteristic2 models in the list. No full representation classification or model transport is imported. Conjugacy-class sizes210 and105 are separately proved and enumerated in this programme.",
        "scope": "Identification of the full alternating-group control and binary linear group. Distinct meridian classes prove a precise difference from the pinned binary rank-one search. No higher-field exclusion, new knot case or global priority conclusion."}
    catalogue = json.loads((HERE / "sources.json").read_text())
    known = {r["id"]: r for r in catalogue["sources"]}
    if record["id"] in known:
        assert known[record["id"]] == record
    else:
        catalogue["sources"].append(record)
        (HERE / "sources.json").write_text(json.dumps(catalogue, indent=2, ensure_ascii=False) + "\n")
    assert len(catalogue["sources"]) == len({r["id"] for r in catalogue["sources"]}) >= 40
    pinned = (HERE / "sources/current13/search/glsearch.py").read_text()
    assert "I + v f" in pinned and "seed vectors e_1..e_4" in pinned
    files = ["audit_double_transposition_sources.py", "sources/atlas-a8.html", "sources.json",
             "sources/current13/search/glsearch.py"]
    report = {"status": "PASS scoped double-transposition source audit", "new_record": record,
        "catalogue_distinct_ids": len(catalogue["sources"]), "frozen_main_manifest_unchanged": True,
        "imports": "Retained Cornwell torus presentation/meridian convention and Rost-Zieschang known torus MRC. BK20 and JP25 finite permutation/cycle quotient precedents credited. No new known-family or global priority audit completed.",
        "pinned_external_coverage": "Binary rank-one meridians are transvections of class size105 in GL4(2); a full A8 double-transposition quotient has meridian class size210. Conjugacy-class size is invariant under any group isomorphism. This coverage distinction does not challenge any valid deposited certificate or classify higher-field representations.",
        "input_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in files}}
    (HERE / "double_transposition_source_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "catalogue": len(catalogue["sources"])}))


if __name__ == "__main__":
    main()
