#!/usr/bin/env python3
"""Read pinned KnotInfo wheels; retain selected source columns without edits."""

import ast
import csv
import hashlib
import io
import json
from pathlib import Path
import urllib.request
import zipfile

HERE = Path(__file__).resolve().parent
RELEASES = [
    {"version": "2026.10.5", "url": "https://files.pythonhosted.org/packages/b3/46/d3acf7ecdd0d9dfd51edb6b12f15b8d997277567aa727a1250f2a68d6403/database_knotinfo-2026.10.5-py3-none-any.whl",
     "sha256": "7cb8cd135ddfdf439d1b799bb4597ed74bcee7b6e3a880b169966f428944a1ca"},
    {"version": "2026.11.2", "url": "https://files.pythonhosted.org/packages/ec/db/4e62ca9f13d5cd1d53ef5cd2aaeda2d654f0d62d7cc2f66d8ca96a097fb6/database_knotinfo-2026.11.2-py3-none-any.whl",
     "sha256": "5a91e8b77d24fb7d90a4a7048074acb97ca81cde205873983344a64458ad88a3"},
]


def main():
    certs = json.loads((HERE / "sources/current13/data/certificates.json").read_text())
    open_names = json.loads((HERE / "sources/current13/data/still_open.json").read_text())
    certified = {c["name"]: c for c in certs}
    wanted = set(certified) | set(open_names)
    assert len(wanted) == 536 and len(set(certified) & set(open_names)) == 0
    cache = HERE / "tmp/sources"
    cache.mkdir(parents=True, exist_ok=True)
    retained = HERE / "sources/knotinfo"
    retained.mkdir(exist_ok=True)
    csv.field_size_limit(10**9)
    results = []
    for release in RELEASES:
        archive_path = cache / release["url"].rsplit("/", 1)[1]
        if not archive_path.exists():
            with urllib.request.urlopen(release["url"], timeout=40) as response:
                archive_path.write_bytes(response.read())
        assert hashlib.sha256(archive_path.read_bytes()).hexdigest() == release["sha256"]
        with zipfile.ZipFile(archive_path) as archive:
            member = next(p for p in archive.namelist() if p.endswith("/knotinfo_data_complete.csv"))
            with archive.open(member) as stream:
                csv_hash = hashlib.file_digest(stream, "sha256").hexdigest()
            with archive.open(member) as stream:
                reader = csv.reader(io.TextIOWrapper(stream, encoding="utf-8-sig", newline=""), delimiter="|")
                header = next(reader)
                columns = {name: i for i, name in enumerate(header)}
                selected, uncertain = {}, []
                for row in reader:
                    name = row[columns["name"]]
                    bridge = row[columns["bridge_index"]]
                    if name.startswith("13") and bridge == "[3,4]":
                        uncertain.append(name)
                    if name in wanted:
                        selected[name] = {"name": name, "bridge_index": bridge,
                                          "pd_notation": ast.literal_eval(row[columns["pd_notation"]])}
        assert set(selected) == wanted
        exact_pd_matches = sum(selected[name]["pd_notation"] == cert["pd"] for name, cert in certified.items())
        subset = {"release": release, "csv_member": member, "csv_sha256": csv_hash,
                  "columns_retained": ["name", "bridge_index", "pd_notation"],
                  "selection": "All536 names in the pinned author243 certificates plus293 still-open list; no replacement",
                  "records": [selected[name] for name in sorted(wanted)]}
        path = retained / (release["version"] + "_selected.json")
        if path.exists():
            assert json.loads(path.read_text()) == subset
        else:
            path.write_text(json.dumps(subset, indent=2) + "\n")
        result = {"version": release["version"], "selected": len(selected),
                  "all_thirteen_crossing_3_or_4_entries": len(uncertain),
                  "author536_exactly_match_uncertain_set": set(uncertain) == wanted,
                  "author243_PD_exact_matches": exact_pd_matches,
                  "author243_still_table_uncertain": sum(selected[name]["bridge_index"] == "[3,4]" for name in certified),
                  "author293_still_table_uncertain": sum(selected[name]["bridge_index"] == "[3,4]" for name in open_names),
                  "subset": str(path.relative_to(HERE)), "subset_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        results.append(result)
    report = {"status": "PASS read-only snapshot audit", "results": results,
              "input_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest()
                               for p in ["audit_knotinfo_snapshots.py", "sources/current13/data/certificates.json",
                                         "sources/current13/data/still_open.json"]},
              "claim_boundary": "The two archived releases are checked, not live KnotInfo on every later date. A table uncertainty entry does not negate the separate certificate proof. The selected source columns are extracted read-only; original archives/CSV bytes were not edited."}
    (HERE / "knotinfo_snapshot_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
