#!/usr/bin/env python3
"""Check every October certificate's named diagram in a second census.

The exact dart-map verifier is independent of the census software.
No numerical isometry or invariant-only identity is accepted.
"""

import hashlib
import importlib.metadata as metadata
import json
from pathlib import Path
import warnings

warnings.filterwarnings("ignore", message="Plink failed to import tkinter.*")
from spherogram import Link
import snappy_manifolds
from verify_diagram_audit import map_witness, plane_map

HERE = Path(__file__).resolve().parent


def main():
    source = HERE / "sources/current13/data/certificates.json"
    certificates = json.loads(source.read_text())
    receipt = json.loads((HERE / "current13_certificate_verification.json").read_text())
    assert receipt["status"] == "PASS" and len(certificates) == receipt["certificates"] == 243
    for path, digest in receipt["input_sha256"].items():
        assert hashlib.sha256((HERE / path).read_bytes()).hexdigest() == digest
    versions = {name: metadata.version(name) for name in ["spherogram", "snappy_manifolds"]}
    assert versions == {"spherogram": "2.4.1", "snappy_manifolds": "1.4"}
    inputs = ["sources/current13/data/certificates.json", "current13_certificate_verification.json",
              "verify_diagram_audit.py", "audit_current13_names.py"]
    report = {"status": "INCOMPLETE", "versions": versions,
              "input_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in inputs},
              "records": [], "failures": [], "census_data": [],
              "claim_boundary": "Exact planar-map matches identify the certified knot types up to mirror in a separately packaged HT census. The already-checked pseudo-reflection and four-seed arguments therefore apply to these named types. Current KnotInfo uncertainty entries, publication priority and the full paper are separate questions."}
    package = Path(snappy_manifolds.__file__).parent
    for path in sorted(package.rglob("*.sqlite")):
        report["census_data"].append({"file": str(path.relative_to(package)), "bytes": path.stat().st_size,
                                     "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    output = HERE / "current13_name_audit.json"
    for c in certificates:
        name = c["name"]
        lookup = "K" + name.replace("_", "")
        link = Link(lookup)
        pd = [list(row) for row in link.PD_code()]
        plane_map(pd)
        try:
            witness = map_witness(c["pd"], pd)
        except AssertionError as error:
            report["failures"].append({"name": name, "reason": str(error)})
        else:
            report["records"].append({"name": name, "lookup": lookup, "pd": pd,
                                       "DT": link.DT_code(), "source_to_census": witness})
        output.write_text(json.dumps(report, indent=2) + "\n")
    report["status"] = "COMPLETE"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "exact_named_diagram_matches": len(report["records"]),
                      "failures": report["failures"]}, indent=2))


if __name__ == "__main__":
    main()
