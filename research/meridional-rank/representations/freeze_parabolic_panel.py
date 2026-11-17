#!/usr/bin/env python3
"""Freeze the same current13 targets with oriented primary PD data."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    panel = json.loads((HERE / "current13_panel.json").read_text())
    source = "sources/knotinfo/2026.11.2_selected.json"
    rows = {r["name"]: r for r in json.loads((HERE / source).read_text())["records"]}
    cases = []
    for c in panel["cases"]:
        row = rows[c["name"]]
        assert row["bridge_index"] == "[3,4]"
        cases.append({"name": c["name"], "pd": row["pd_notation"], "reported_bridge_interval": [3, 4]})
    plan = {"status": "FROZEN_BEFORE_PARABOLIC_TEST", "selection": "All same16 names from current13_panel.json, unchanged order and no replacements",
            "motivation": "An SL2(Fp) rank-one transvection base may have restricted affine-root lift dimension4 even when no S3/S4 transposition base exists. Evaluation of the lift space on arbitrary meridional generators gives mr>=dimension, without requiring the single base image to have rank4. This differs from the October nondegenerate Sp4 fixed-hyperplane certificates.",
            "sources": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in [source, "current13_panel.json"]},
            "limits": {"primes": [3, 5, 7], "seed_count": 4, "max_seconds_per_case": 25,
                       "max_seed_trials_per_case": 1000, "no_adaptive_replacements": True},
            "cases": cases,
            "claim_boundary": "Complete enumerations certify only these fields and this transvection conjugacy class, up to SL2 conjugation. No-hit does not exclude other primes or quotient types. A sharp hit requires independent matrix/seed/name replay and prior-resolution audit."}
    out = HERE / "parabolic_panel.json"
    if out.exists():
        assert json.loads(out.read_text()) == plan
    else:
        out.write_text(json.dumps(plan, indent=2) + "\n")
    print(json.dumps({"status": "FROZEN", "targets": len(cases), "primes": plan["limits"]["primes"]}))


if __name__ == "__main__":
    main()
