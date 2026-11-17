#!/usr/bin/env python3
"""Exact all-characteristic lift dimensions for the frozen 13-crossing panel."""

from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path
import time

from integer_lifts import relator_matrix, smith_certificate
from search_lift_panel import diagram, seed_certificate, conjugation_table, extend_colouring
from signed_lifts import sign_potentials

HERE = Path(__file__).resolve().parent


def test(case, limits):
    start = time.monotonic()
    deadline = start + limits["max_seconds_per_case"]
    n, crossings = diagram(case["gauss_code"])
    seeds, trace, trials = seed_certificate(n, crossings, limits["seed_count"],
                                           limits["max_seed_trials_per_case"], deadline)
    result = {"name": case["name"], "status": "INCOMPLETE", "seed_trials": trials,
              "by_degree": []}
    if seeds is None:
        result["reason"] = "no seed certificate within cap"
        return result
    result.update({"bridge_upper_bound": len(seeds), "seeds": seeds, "seed_trace": trace})
    words = [[o + 1, a + 1, o + 1, b + 1] for o, a, b in crossings]
    for r in limits["degrees"]:
        edges, table = conjugation_table(r)
        counts, certificates = Counter(), []
        for tail in product(range(len(edges)), repeat=len(seeds) - 1):
            if time.monotonic() > deadline:
                result["reason"] = "time cap"
                return result
            assignment = (0,) + tail
            if sign_potentials(r, [edges[x] for x in assignment], 0)[0] != 1:
                continue
            labels = extend_colouring(n, crossings, seeds, assignment, table)
            if labels is None:
                continue
            transpositions = [edges[x] for x in labels]
            matrix = relator_matrix(r, transpositions, words)
            smith = smith_certificate(matrix, n)
            bound = smith["maximum_dimension_all_characteristics"]
            assert bound <= len(seeds)
            counts[bound] += 1
            certificates.append({"seed_assignment": assignment, "arc_labels": labels,
                                 "smith": smith})
        result["by_degree"].append({"degree": r, "normalized_seed_assignments": len(edges) ** (len(seeds) - 1),
                                    "maximum_over_all_characteristics_counts": dict(sorted(counts.items())),
                                    "certificates": certificates})
    result["status"] = "COMPLETE"
    result["sharp"] = any(c["smith"]["maximum_dimension_all_characteristics"] == len(seeds)
                          for d in result["by_degree"] for c in d["certificates"])
    result["seconds"] = round(time.monotonic() - start, 6)
    return result


def main():
    panel_path = HERE / "current13_panel.json"
    panel = json.loads(panel_path.read_text())
    hashes = {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest()
              for p in ["current13_panel.json", "integer_lifts.py", "search_integer_lifts.py",
                        "search_lift_panel.py", "signed_lifts.py"]}
    output = HERE / "current13_integer_results.json"
    if output.exists():
        report = json.loads(output.read_text())
        assert report["input_sha256"] == hashes
    else:
        report = {"input_sha256": hashes, "results": [], "claim_boundary": panel["claim_boundary"]}
    done = {r["name"] for r in report["results"]}
    for case in panel["cases"]:
        if case["name"] in done:
            continue
        result = test(case, panel["search_limits"])
        report["results"].append(result)
        tmp = output.with_suffix(".tmp")
        tmp.write_text(json.dumps(report, indent=2) + "\n")
        tmp.replace(output)
        print(json.dumps({"name": case["name"], "status": result["status"], "sharp": result.get("sharp"),
                          "dimensions": [{"degree": d["degree"], "counts": d["maximum_over_all_characteristics_counts"]}
                                         for d in result["by_degree"]]}), flush=True)


if __name__ == "__main__":
    main()
