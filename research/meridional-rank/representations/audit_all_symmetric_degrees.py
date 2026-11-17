#!/usr/bin/env python3
"""Exact all-degree/all-field scope audit on the unchanged current13 panel.

The four-seed certificate bounds any transitive transposition quotient
degree by five. Exhaust S2..S5, with integral Smith forms for every base.
"""

from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path

from integer_lifts import relator_matrix, smith_certificate
from parabolic_lifts import oriented_diagram
from search_lift_panel import conjugation_table, extend_colouring
from signed_lifts import sign_potentials

HERE = Path(__file__).resolve().parent


def main():
    panel = json.loads((HERE / "parabolic_panel.json").read_text())
    old = json.loads((HERE / "parabolic_panel_results.json").read_text())
    plan = {"selection": "All same16 current13 cases and oriented source PDs, unchanged order and no replacement",
            "source_panel_sha256": hashlib.sha256((HERE / "parabolic_panel.json").read_bytes()).hexdigest(),
            "degrees": [2, 3, 4, 5], "fields": "Every characteristic, via exact integer Smith form",
            "motivation": "Four meridional seeds give w<=4. A connected graph of w transpositions has at most w+1 vertices, so every transitive symmetric quotient has degree<=5. This exhausts the missing S2 and S5 cases as well as replaying S3/S4 on the independent PDs.",
            "max_normalized_assignments_per_case": 1244, "no_adaptive_replacements": True}
    plan_path = HERE / "all_symmetric_degrees_plan.json"
    if not plan_path.exists():
        plan_path.write_text(json.dumps(plan, indent=2) + "\n")
        print("FROZEN all_symmetric_degrees_plan.json; rerun to execute")
        return
    assert json.loads(plan_path.read_text()) == plan
    inputs = ["all_symmetric_degrees_plan.json", "parabolic_panel.json", "parabolic_panel_results.json",
              "audit_all_symmetric_degrees.py", "integer_lifts.py", "parabolic_lifts.py",
              "search_lift_panel.py", "signed_lifts.py"]
    hashes = {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in inputs}
    output = HERE / "all_symmetric_degrees_results.json"
    if output.exists():
        report = json.loads(output.read_text())
        assert report["input_sha256"] == hashes
    else:
        report = {"input_sha256": hashes, "results": [],
                  "claim_boundary": "Uniform graph-generator argument covers r>=6; exact enumeration covers r2..5 and every field. This concerns the specified one-scalar symmetric reflection lift bound, not arbitrary quotient representations or nonabelian coefficients."}
    done = {r["name"] for r in report["results"]}
    for case, previous in zip(panel["cases"], old["results"]):
        assert case["name"] == previous["name"] and previous["status"] == "COMPLETE"
        if case["name"] in done:
            continue
        n = len(case["pd"])
        crossings = [c[:3] for c in oriented_diagram(case["pd"])]
        words = [[o+1,a+1,o+1,b+1] for o,a,b in crossings]
        result = {"name": case["name"], "seeds": previous["seeds"], "seed_trace": previous["seed_trace"],
                  "bridge_upper_bound": 4, "degrees": [], "status": "COMPLETE"}
        for degree in plan["degrees"]:
            edges, table = conjugation_table(degree)
            certificates, counts = [], Counter()
            for tail in product(range(len(edges)), repeat=3):
                assignment = (0,) + tail
                if sign_potentials(degree, [edges[i] for i in assignment], 0)[0] != 1:
                    continue
                labels = extend_colouring(n,crossings,previous["seeds"],assignment,table)
                if labels is None:
                    continue
                roots = [edges[i] for i in labels]
                matrix = relator_matrix(degree,roots,words)
                smith = smith_certificate(matrix,n)
                counts[smith["maximum_dimension_all_characteristics"]] += 1
                certificates.append({"assignment": assignment, "arc_labels": labels, "smith": smith})
            result["degrees"].append({"degree": degree, "normalized_assignments": len(edges)**3,
                                       "all_field_bound_counts": dict(sorted(counts.items())),
                                       "certificates": certificates})
        result["maximum_scalar_lift_bound_all_degrees_fields"] = max(int(b) for d in result["degrees"] for b in d["all_field_bound_counts"])
        report["results"].append(result)
        output.write_text(json.dumps(report,indent=2)+"\n")
        print(json.dumps({"name":case["name"], "maximum":result["maximum_scalar_lift_bound_all_degrees_fields"],
                          "bases":{d["degree"]:len(d["certificates"]) for d in result["degrees"]}}),flush=True)


if __name__ == "__main__":
    main()
