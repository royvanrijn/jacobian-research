#!/usr/bin/env sage-python
"""Correct the cohort-aware nontrivial-cover counts in the immutable v2 runs."""

from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts" / "generated-results" / "elliptic-curves"
CONTROL = ART / "curve302_strict_cover_control_benchmark_v2.json"
MISSING = ART / "curve302_strict_cover_missing_panel_v2.json"
OUTPUT = ART / "curve302_strict_cover_benchmark_v2_metric_audit_v1.json"


def read(path: Path):
    return json.loads(path.read_text())


def digest(path: Path):
    return sha256(path.read_bytes()).hexdigest()


def exact_count(payload, cohort):
    rows = payload["cases"]
    assert len(rows) == 8 and all(row["cohort"] == cohort for row in rows)
    # The zero member only means the trivial cover in H.  R+0 is R, so every
    # member of the affine missing coset is nontrivial.
    nontrivial = [row for row in rows if cohort != "recovered_control" or row["member"] != 0]
    return sum(bool(row["exact_rational_cover_points"]) for row in nontrivial), len(nontrivial)


def main():
    controls, missing = read(CONTROL), read(MISSING)
    assert controls["status"] == missing["status"] == "PASS"
    control_count, control_total = exact_count(controls, "recovered_control")
    missing_count, missing_total = exact_count(missing, "missing_strict_coset")
    assert (control_count, control_total) == (7, 7)
    assert (missing_count, missing_total) == (8, 8)
    record = {
        "schema": "elliptic-curves.curve302-strict-cover-benchmark-v2-metric-audit.v1",
        "status": "PASS",
        "bindings": {str(path.relative_to(ROOT)): digest(path) for path in (CONTROL, MISSING, Path(__file__))},
        "cohort_aware_counts": {
            "recovered_nontrivial_H_controls": {"solved": control_count, "total": control_total},
            "nontrivial_missing_R_plus_H_coset": {"solved": missing_count, "total": missing_total},
        },
        "correction": "The v2 comparison field applies member!=0 to both cohorts. That convention is correct for H but excludes R+0 in the affine R+H panel. This audit supplies the cohort-aware count without rewriting either immutable v2 artifact.",
    }
    with OUTPUT.open("x") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print("PASS", record["cohort_aware_counts"], flush=True)


if __name__ == "__main__":
    main()
