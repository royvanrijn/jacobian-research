#!/usr/bin/env python3
"""Preserve all frozen candidates; check their 32 marked half-turn variants."""

import argparse
import hashlib
import json
from pathlib import Path

from audit_diagram_identity import Link, RationalTangle, normalize, pd_walk

HERE = Path(__file__).resolve().parent


def digest(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", action="store_true")
    args = parser.parse_args()
    plan_name = "montesinos_rotation_plan.json"
    original = json.loads((HERE / "montesinos_audit_plan.json").read_text())
    plan = {"status": "FROZEN_BEFORE_ROTATION_TEST", "original_plan_sha256": digest("montesinos_audit_plan.json"),
            "targets": original["targets"], "candidate_count": 112, "rotation_masks": list(range(32)),
            "max_diagrams": 3584, "no_adaptive_replacements": True,
            "motivation": "The first exact map test found seven matches. Keep every original fraction tuple and target; for each tangle independently allow a 180-degree rotation of its marked square. Such a rotation preserves rational-tangle slope (two applications of x->-1/x) and is a different diagram convention for the same Montesinos presentation. This is a fixed second-stage identity audit, not a replacement search.",
            "acceptance": "Exact full signed-Gauss planar-map identity up to global spatial reflections, with PD and transforms retained. Misses do not exclude Montesinos membership."}
    path = HERE / plan_name
    if args.freeze:
        if path.exists():
            assert json.loads(path.read_text()) == plan
        else:
            path.write_text(json.dumps(plan, indent=2) + "\n")
        print(json.dumps({"status": "FROZEN", "diagrams": 3584}))
        return
    assert json.loads(path.read_text()) == plan
    source = json.loads((HERE / "diagram_identity_audit.json").read_text())
    keys = {normalize(pd_walk(r["census_pd"]), True)[0]: r["name"] for r in source["records"]}
    output = HERE / "montesinos_rotation_results.json"
    hashes = {p: digest(p) for p in [plan_name, "montesinos_audit_plan.json", "diagram_identity_audit.json",
                                   "search_montesinos_rotations.py", "audit_diagram_identity.py"]}
    if output.exists():
        report = json.loads(output.read_text())
        assert report["input_sha256"] == hashes
    else:
        report = {"status": "INCOMPLETE", "input_sha256": hashes, "completed_candidates": 0,
                  "diagrams_checked": 0, "matches": {}, "acceptance": plan["acceptance"]}
    for index, c in enumerate(original["candidates"]):
        if index < report["completed_candidates"]:
            continue
        for mask in plan["rotation_masks"]:
            tangle = RationalTangle(c["integer_twists"])
            for i, (a, b) in enumerate(c["fractions"]):
                term = RationalTangle(a, b)
                if mask >> i & 1:
                    term = term.rotate(2)
                tangle = tangle + term
            link = tangle.numerator_closure()
            if len(link.link_components) == 1:
                pd = [list(crossing) for crossing in link.PD_code()]
                key, change = normalize(pd_walk(pd), True)
                name = keys.get(key)
                if name is not None and name not in report["matches"]:
                    report["matches"][name] = {"candidate_index": index, "half_turn_mask": mask,
                                               **c, "pd": pd, "full_transform": change}
            report["diagrams_checked"] += 1
        report["completed_candidates"] = index + 1
        output.write_text(json.dumps(report, indent=2) + "\n")
    report["status"] = "COMPLETE"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": "COMPLETE", "diagrams": report["diagrams_checked"],
                      "matches": [{"name": name, "fractions": c["fractions"], "mask": c["half_turn_mask"]}
                                  for name, c in sorted(report["matches"].items())]}, indent=2))


if __name__ == "__main__":
    main()
