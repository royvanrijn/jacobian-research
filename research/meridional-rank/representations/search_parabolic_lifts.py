#!/usr/bin/env python3
"""Checkpoint every target in the fixed small parabolic-root lift panel."""

from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path
import time

from parabolic_lifts import roots_and_tables, oriented_diagram, extend, lift_kernel
from search_lift_panel import seed_certificate

HERE = Path(__file__).resolve().parent


def test(case, limits):
    started = time.monotonic()
    deadline = started + limits["max_seconds_per_case"]
    cr = oriented_diagram(case["pd"])
    n = len(case["pd"])
    seeds, trace, trials = seed_certificate(n, [c[:3] for c in cr], limits["seed_count"],
                                           limits["max_seed_trials_per_case"], deadline)
    result = {"name": case["name"], "status": "INCOMPLETE", "by_prime": [], "seed_trials": trials}
    if seeds is None:
        result["reason"] = "seed cap"
        return result
    result.update({"seeds": seeds, "seed_trace": trace, "bridge_upper_bound": 4})
    for p in limits["primes"]:
        roots, matrices, forward, backward = roots_and_tables(p)
        counts, best, accepted, tested = Counter(), None, 0, 0
        for tail in product(range(len(roots)), repeat=3):
            if time.monotonic() > deadline:
                result["reason"] = "time cap"
                return result
            assignment = (0,) + tail
            tested += 1
            labels = extend(n, cr, seeds, assignment, forward, backward)
            if labels is None:
                continue
            # Skip entirely parallel root families (commuting base images).
            if all((roots[0][0]*roots[i][1]-roots[0][1]*roots[i][0]) % p == 0 for i in labels):
                continue
            accepted += 1
            basis = lift_kernel(cr, labels, roots, matrices, p)
            bound = len(basis)
            assert bound <= 4
            counts[bound] += 1
            if best is None or bound > best["bound"]:
                best = {"bound": bound, "arc_roots": [roots[i] for i in labels],
                        "arc_matrices": [matrices[i] for i in labels], "basis": basis}
        result["by_prime"].append({"prime": p, "normalized_assignments": tested,
                                   "nonparallel_base_colourings": accepted,
                                   "restricted_lift_dimension_counts": dict(sorted(counts.items())),
                                   "best": best})
    result["status"] = "COMPLETE"
    result["sharp"] = any(r["best"] is not None and r["best"]["bound"] == 4 for r in result["by_prime"])
    result["seconds"] = round(time.monotonic()-started,6)
    return result


def main():
    panel = json.loads((HERE / "parabolic_panel.json").read_text())
    hashes = {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in
              ["parabolic_panel.json", "parabolic_lifts.py", "search_parabolic_lifts.py", "verify_integer_lifts.py", "search_lift_panel.py"]}
    out = HERE / "parabolic_panel_results.json"
    if out.exists():
        report = json.loads(out.read_text())
        assert report["input_sha256"] == hashes
    else:
        report = {"input_sha256": hashes, "results": [], "claim_boundary": panel["claim_boundary"]}
    done = {r["name"] for r in report["results"]}
    for case in panel["cases"]:
        if case["name"] in done:
            continue
        result = test(case,panel["limits"])
        report["results"].append(result)
        tmp = out.with_suffix(".tmp")
        tmp.write_text(json.dumps(report,indent=2)+"\n")
        tmp.replace(out)
        print(json.dumps({"name": case["name"], "status": result["status"], "sharp": result.get("sharp"),
                          "bounds": [{"prime": x["prime"], "counts": x["restricted_lift_dimension_counts"]}
                                     for x in result["by_prime"]]}),flush=True)


if __name__ == "__main__":
    main()
