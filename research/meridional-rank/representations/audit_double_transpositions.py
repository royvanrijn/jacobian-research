#!/usr/bin/env python3
"""Frozen sixteen-target degree-eight double-transposition search."""

from collections import Counter
import hashlib
import json
from pathlib import Path
import time

from double_transpositions import (IDENTITY, compose, generated_group, model,
                                  power, support_connected, tuple_representatives)
from parabolic_lifts import oriented_diagram
from search_lift_panel import extend_colouring

HERE = Path(__file__).resolve().parent


def main():
    started = time.monotonic()
    plan = json.loads((HERE / "double_transposition_plan.json").read_text())
    assert plan["frozen_before_execution"]
    inputs = ["double_transposition_plan.json", "double_transpositions.py",
              "audit_double_transpositions.py", "parabolic_panel.json",
              "parabolic_panel_results.json", "parabolic_lifts.py", "search_lift_panel.py",
              "sources/current13/search/glsearch.py"]
    hashes = {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in inputs}
    path = HERE / "double_transposition_results.json"
    report = {"status": "RUNNING", "input_sha256": hashes, "results": []}
    if path.exists():
        report = json.loads(path.read_text())
        assert report["input_sha256"] == hashes
    elements, table, centralizer, actions = model()
    representatives, level_counts = tuple_representatives(actions)
    assert len(representatives) < plan["limits"]["tuple_orbit_cap"]
    connected = [(t, s, w) for t, s, w in representatives
                 if support_connected([elements[i] for i in t])]
    report["tuple_catalogue"] = {"class_size": len(elements), "centralizer_size": len(centralizer),
        "prefix_counts": level_counts, "four_tuple_orbits": len(representatives),
        "raw_normalized_tuples": sum(w for _, _, w in representatives),
        "transitive_orbits": len(connected),
        "transitive_raw_normalized_tuples": sum(w for _, _, w in connected),
        "stabilizer_histogram": dict(sorted(Counter(s for _, s, _ in representatives).items()))}
    report["controls"] = []
    for control in plan["controls"]:
        y = tuple(control["y"])
        q, exponent = (7, 2) if control["name"] == "T(4,7)" else (15, 4)
        assert power(y, q) == IDENTITY
        found = []
        for i, m in enumerate(elements):
            x = compose(power(y, exponent), m)
            if power(x, 4) != IDENTITY or not support_connected([x, y]):
                continue
            group = generated_group([x, y])
            found.append({"meridian_index": i, "x": x, "y": y, "m": m,
                          "image_order": len(group)})
        report["controls"].append({"name": control["name"], "pairs_checked": 210,
            "transitive_pairs": found, "image_orders": dict(Counter(f["image_order"] for f in found))})
    assert any(c["transitive_pairs"] for c in report["controls"]), "no positive control"
    # The canonical tuple catalogue is retained independently of results.
    catalogue = {"elements": elements, "centralizer": centralizer,
                 "representatives": [{"tuple": t, "stabilizer": s, "raw_weight": w,
                    "transitive": support_connected([elements[i] for i in t])}
                     for t, s, w in representatives], "counts": report["tuple_catalogue"]}
    catalogue_path = HERE / "double_transposition_catalogue.json"
    catalogue_path.write_text(json.dumps(catalogue, separators=(",", ":")) + "\n")
    report["catalogue_sha256"] = hashlib.sha256(catalogue_path.read_bytes()).hexdigest()
    panel = json.loads((HERE / "parabolic_panel.json").read_text())["cases"]
    old = json.loads((HERE / "parabolic_panel_results.json").read_text())["results"]
    assert [c["name"] for c in panel] == plan["targets"]
    done = {r["name"] for r in report["results"]}
    for case, previous in zip(panel, old):
        if case["name"] in done:
            continue
        if time.monotonic() - started > plan["limits"]["seconds_total"]:
            break
        assert case["name"] == previous["name"] and previous["status"] == "COMPLETE"
        tick = time.monotonic()
        crossings = [c[:3] for c in oriented_diagram(case["pd"])]
        tested, successes = 0, []
        for assignment, stabilizer, weight in connected:
            if time.monotonic() - tick > plan["limits"]["seconds_per_target"]:
                break
            tested += 1
            labels = extend_colouring(len(crossings), crossings, previous["seeds"], assignment, table)
            if labels is None:
                continue
            assert all(table[labels[o]][labels[a]] == labels[b] for o, a, b in crossings)
            group = generated_group([elements[i] for i in assignment])
            successes.append({"seed_assignment": assignment, "arc_labels": labels,
                              "arc_permutations": [elements[i] for i in labels],
                              "stabilizer": stabilizer, "raw_weight": weight,
                              "image_order": len(group), "meridional_rank_lower_bound": 4})
        result = {"name": case["name"], "status": "COMPLETE" if tested == len(connected) else "CAPPED_INCOMPLETE",
            "seeds": previous["seeds"], "seed_trace": previous["seed_trace"],
            "bridge_upper_bound": 4, "transitive_seed_orbits_tested": tested,
            "certificates": successes, "seconds": round(time.monotonic() - tick, 6)}
        report["results"].append(result)
        path.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps({k: result[k] for k in ["name", "status", "transitive_seed_orbits_tested", "seconds"]}
                         | {"certificates": len(successes)}), flush=True)
    report["status"] = "COMPLETE" if len(report["results"]) == 16 and all(r["status"] == "COMPLETE" for r in report["results"]) else "CAPPED_INCOMPLETE"
    report["sharp_names"] = [r["name"] for r in report["results"] if r["certificates"]]
    report["seconds_total"] = round(time.monotonic() - started, 6)
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ["status", "tuple_catalogue", "sharp_names", "seconds_total"]}, indent=2))


if __name__ == "__main__":
    main()
