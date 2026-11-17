#!/usr/bin/env python3
"""Freeze a small diagram-identity/family audit, without testing identities."""

import hashlib
from itertools import product
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    candidates = []
    thirds = [[1, 3], [2, 3]]
    for fractions in product(thirds, repeat=5):
        candidates.append({"integer_twists": 0, "fractions": list(fractions)})
    for position in range(1, 5):
        for fifth in [[2, 5], [3, 5]]:
            for rest in product(thirds, repeat=3):
                fractions = [[1, 2]]
                remaining = iter(rest)
                for i in range(1, 5):
                    fractions.append(fifth if i == position else next(remaining))
                candidates.append({"integer_twists": 0, "fractions": fractions})
    for rest in product(thirds, repeat=4):
        candidates.append({"integer_twists": 1, "fractions": [[1, 2], *rest]})
    assert len(candidates) == 112
    sharp = HERE / "sharp_five_certificates.json"
    panel = json.loads(sharp.read_text())
    plan = {"date": "2026-11-16", "status": "FROZEN_BEFORE_FAMILY_TEST",
            "target_source": sharp.name, "target_sha256": hashlib.sha256(sharp.read_bytes()).hexdigest(),
            "targets": [c["name"] for c in panel["certificates"]],
            "motivation": "At five rational tangles and fifteen crossings, positive alternating candidates have all five length3 tangles; or one length2 and one length4 with three length3; or one extra integer crossing with one length2 and four length3. Nonintegral reduced slopes with at most one even denominator motivate this finite family. This is not an exhaustive Montesinos recognizer or a crossing-number theorem.",
            "candidate_count": len(candidates), "candidates": candidates,
            "method": "Numerator closure of sum of the specified rational tangles and integer tangle. Exact full signed-Gauss/planar-map matching up to traversal choices and global reflections. No numerical isometry or invariant-only acceptance.",
            "limits": {"max_candidates": 112, "no_adaptive_replacements": True},
            "claim_boundary": "A match supplies an explicit five-noninteger-tangle Montesinos presentation up to mirror. A miss says nothing about family membership. Prior MRC theorem scope must be sourced separately."}
    target = HERE / "montesinos_audit_plan.json"
    if target.exists():
        assert json.loads(target.read_text()) == plan
    else:
        target.write_text(json.dumps(plan, indent=2) + "\n")
    print(json.dumps({"status": "FROZEN", "targets": len(plan["targets"]), "candidates": len(candidates)}))


if __name__ == "__main__":
    main()
