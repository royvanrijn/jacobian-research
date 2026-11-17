#!/usr/bin/env python3
"""New coefficient-field audit of the retained five-seed panel, not discovery."""

import hashlib
import json
from pathlib import Path

from search_integer_lifts import test

HERE = Path(__file__).resolve().parent


def main():
    panel = json.loads((HERE / "lift_panel.json").read_text())
    plan = {"source_panel": "lift_panel.json", "source_sha256": hashlib.sha256((HERE / "lift_panel.json").read_bytes()).hexdigest(),
            "scope": "All-characteristic extension on the original sixteen cases, preserving all cases and reporting no discovery claim",
            "limits": {**panel["search_limits"], "max_seconds_per_case": 20,
                       "coefficient_fields": "All characteristics via exact integer Smith form"}}
    plan_path = HERE / "old_panel_integer_plan.json"
    if plan_path.exists():
        assert json.loads(plan_path.read_text()) == plan
    else:
        plan_path.write_text(json.dumps(plan, indent=2) + "\n")
    inputs = ["old_panel_integer_plan.json", "lift_panel.json", "integer_lifts.py",
              "search_integer_lifts.py", "search_old_panel_integer.py", "search_lift_panel.py", "signed_lifts.py"]
    hashes = {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in inputs}
    output = HERE / "old_panel_integer_results.json"
    if output.exists():
        report = json.loads(output.read_text())
        assert report["input_sha256"] == hashes
    else:
        report = {"input_sha256": hashes, "scope": plan["scope"], "results": []}
    done = {r["name"] for r in report["results"]}
    for case in panel["cases"]:
        if case["name"] in done:
            continue
        result = test(case, plan["limits"])
        report["results"].append(result)
        tmp = output.with_suffix(".tmp")
        tmp.write_text(json.dumps(report, indent=2) + "\n")
        tmp.replace(output)
        print(json.dumps({"name": case["name"], "status": result["status"], "sharp": result.get("sharp"),
                          "dimensions": [{"degree": d["degree"], "counts": d["maximum_over_all_characteristics_counts"]}
                                         for d in result["by_degree"]]}), flush=True)


if __name__ == "__main__":
    main()
