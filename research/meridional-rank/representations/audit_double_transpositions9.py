#!/usr/bin/env python3
"""Degree-nine completion of the frozen double-transposition graph gate."""

from collections import Counter
import hashlib
import json
from pathlib import Path
import time

from double_transpositions import model, support_connected, tuple_representatives
from parabolic_lifts import oriented_diagram
from search_lift_panel import extend_colouring

HERE = Path(__file__).resolve().parent


def main():
    started = time.monotonic()
    plan = json.loads((HERE / "double_transposition_degree9_plan.json").read_text())
    inputs = ["double_transposition_degree9_plan.json", "double_transposition_plan.json",
              "audit_double_transpositions9.py", "double_transpositions.py", "parabolic_panel.json",
              "parabolic_panel_results.json", "parabolic_lifts.py", "search_lift_panel.py",
              "double_transposition_results.json"]
    hashes = {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in inputs}
    path = HERE / "double_transposition9_results.json"
    report = {"status": "RUNNING", "input_sha256": hashes, "results": []}
    if path.exists():
        report = json.loads(path.read_text())
        assert report["input_sha256"] == hashes
    elements, table, centralizer, actions = model(9)
    representatives, counts = tuple_representatives(actions)
    assert len(representatives) <= plan["limits"]["tuple_orbit_cap"]
    connected = [(t, s, w) for t, s, w in representatives
                 if support_connected([elements[i] for i in t])]
    matchings = [[(0, 1), (2, 3)], [(1, 2), (4, 5)], [(3, 4), (6, 7)], [(5, 6), (7, 8)]]
    control = []
    for pairs in matchings:
        p = list(range(9))
        for a, b in pairs:
            p[a], p[b] = b, a
        control.append(tuple(p))
    assert support_connected(control)
    report["graph_control"] = control
    report["tuple_catalogue"] = {"class_size": len(elements), "centralizer_size": len(centralizer),
        "prefix_counts": counts, "four_tuple_orbits": len(representatives),
        "raw_normalized_tuples": sum(w for _, _, w in representatives),
        "transitive_orbits": len(connected),
        "transitive_raw_normalized_tuples": sum(w for _, _, w in connected),
        "stabilizer_histogram": dict(sorted(Counter(s for _, s, _ in representatives).items()))}
    catalogue = {"elements": elements, "centralizer": centralizer,
        "representatives": [{"tuple": t, "stabilizer": s, "raw_weight": w,
            "transitive": support_connected([elements[i] for i in t])} for t, s, w in representatives],
        "counts": report["tuple_catalogue"]}
    catalogue_path = HERE / "double_transposition9_catalogue.json"
    catalogue_path.write_text(json.dumps(catalogue, separators=(",", ":")) + "\n")
    report["catalogue_sha256"] = hashlib.sha256(catalogue_path.read_bytes()).hexdigest()
    panel = json.loads((HERE / "parabolic_panel.json").read_text())["cases"]
    old = json.loads((HERE / "parabolic_panel_results.json").read_text())["results"]
    assert [c["name"] for c in panel] == json.loads((HERE / "double_transposition_plan.json").read_text())["targets"]
    done = {r["name"] for r in report["results"]}
    for case, previous in zip(panel, old):
        if case["name"] in done:
            continue
        if time.monotonic() - started > plan["limits"]["seconds_total"]:
            break
        tick = time.monotonic()
        crossings = [c[:3] for c in oriented_diagram(case["pd"])]
        successes, tested = [], 0
        for assignment, stabilizer, weight in connected:
            if time.monotonic() - tick > plan["limits"]["seconds_per_target"]:
                break
            tested += 1
            labels = extend_colouring(len(crossings), crossings, previous["seeds"], assignment, table)
            if labels is None:
                continue
            assert all(table[labels[o]][labels[a]] == labels[b] for o, a, b in crossings)
            successes.append({"seed_assignment": assignment, "arc_labels": labels,
                "arc_permutations": [elements[i] for i in labels], "stabilizer": stabilizer,
                "raw_weight": weight, "meridional_rank_lower_bound": 4})
        row = {"name": case["name"], "status": "COMPLETE" if tested == len(connected) else "CAPPED_INCOMPLETE",
            "seeds": previous["seeds"], "seed_trace": previous["seed_trace"], "bridge_upper_bound": 4,
            "transitive_seed_orbits_tested": tested, "certificates": successes,
            "seconds": round(time.monotonic() - tick, 6)}
        report["results"].append(row)
        path.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps({"name": row["name"], "status": row["status"], "tested": tested,
                          "certificates": len(successes)}), flush=True)
    report["status"] = "COMPLETE" if len(report["results"]) == 16 and all(r["status"] == "COMPLETE" for r in report["results"]) else "CAPPED_INCOMPLETE"
    report["sharp_names"] = [r["name"] for r in report["results"] if r["certificates"]]
    report["seconds_total"] = round(time.monotonic() - started, 6)
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ["status", "tuple_catalogue", "sharp_names", "seconds_total"]}, indent=2))


if __name__ == "__main__":
    main()
