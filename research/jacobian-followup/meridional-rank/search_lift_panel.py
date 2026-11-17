#!/usr/bin/env python3
"""Bounded test of the frozen panel; retain every case and checkpoint result.

The unsigned crossing orientation is sufficient for these involution-valued
quotients. Knot names and Gauss-code provenance are the authors' table.
"""

from bisect import bisect_right
from collections import Counter
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
import time

from signed_lifts import nullspace, relator_rows, sign_potentials, solve

HERE = Path(__file__).resolve().parent


def diagram(code):
    n = len(code) // 2
    assert sorted(code) == list(range(-n, 0)) + list(range(1, n + 1))
    negative = [i for i, value in enumerate(code) if value < 0]
    under = {-code[pos]: i for i, pos in enumerate(negative)}
    over = {value: (bisect_right(negative, pos) - 1) % n
            for pos, value in enumerate(code) if value > 0}
    crossings = [(over[c], (under[c] - 1) % n, under[c]) for c in range(1, n + 1)]
    return n, crossings


def seed_certificate(n, crossings, count, trial_cap, deadline):
    for trial, seeds in enumerate(combinations(range(n), count), 1):
        if trial > trial_cap or time.monotonic() > deadline:
            return None, None, trial
        seen, trace = set(seeds), []
        changed = True
        while changed:
            changed = False
            for o, a, b in crossings:
                if o in seen and ((a in seen) ^ (b in seen)):
                    old, new = (a, b) if a in seen else (b, a)
                    seen.add(new)
                    trace.append([o, old, new])
                    changed = True
        if len(seen) == n:
            return list(seeds), trace, trial
    return None, None, trial


def conjugation_table(r):
    edges = list(combinations(range(r), 2))
    indices = {edge: i for i, edge in enumerate(edges)}
    table = []
    for a, b in edges:
        def image(v):
            return b if v == a else a if v == b else v
        table.append([indices[tuple(sorted((image(p), image(q))))] for p, q in edges])
    return edges, table


def extend_colouring(n, crossings, seeds, assignment, table):
    labels = [None] * n
    for seed, label in zip(seeds, assignment):
        labels[seed] = label
    changed = True
    while changed:
        changed = False
        for o, a, b in crossings:
            if labels[o] is None:
                continue
            if labels[a] is not None:
                value = table[labels[o]][labels[a]]
                if labels[b] is None:
                    labels[b] = value
                    changed = True
                elif labels[b] != value:
                    return None
            elif labels[b] is not None:
                labels[a] = table[labels[o]][labels[b]]
                changed = True
    assert all(x is not None for x in labels), "seed propagation did not cover diagram"
    return labels


def test_case(case, limits):
    started = time.monotonic()
    deadline = started + limits["max_seconds_per_case"]
    n, crossings = diagram(case["gauss_code"])
    seeds, trace, trials = seed_certificate(n, crossings, limits["seed_count"],
                                           limits["max_seed_trials_per_case"], deadline)
    result = {"name": case["name"], "crossings": n, "source_row": case["source_row"],
              "seed_trials": trials, "status": "INCOMPLETE", "by_degree": []}
    if seeds is None:
        result["reason"] = "upper certificate not found within cap"
        return result
    result.update({"bridge_upper_bound": len(seeds), "seed_strands": seeds, "seed_propagation": trace})
    relators = [[o + 1, a + 1, o + 1, b + 1] for o, a, b in crossings]
    for r in limits["degrees"]:
        edges, table = conjugation_table(r)
        distribution, best, attempts = Counter(), None, 0
        for tail in product(range(len(edges)), repeat=len(seeds) - 1):
            if time.monotonic() > deadline:
                result["reason"] = f"time cap reached in symmetric degree {r}"
                return result
            assignment = (0,) + tail  # normalize first seed to (0,1)
            attempts += 1
            if sign_potentials(r, [edges[x] for x in assignment], 0)[0] != 1:
                continue
            labels = extend_colouring(n, crossings, seeds, assignment, table)
            if labels is None:
                continue
            transpositions = [edges[x] for x in labels]
            data = {"label": case["name"], "degree": r, "transpositions": transpositions,
                    "relators": relators}
            bound = solve(data)["meridional_rank_lower_bound"]
            assert bound <= len(seeds), "lower bound contradicts independent upper certificate"
            distribution[bound] += 1
            if best is None or bound > best["bound"]:
                basis = nullspace(relator_rows(r, transpositions, relators), n)
                best = {"bound": bound, "presentation": data, "lift_basis_bits": basis}
        result["by_degree"].append({"degree": r, "assignments_tested": attempts,
                                    "surjective_colouring_bound_counts": dict(sorted(distribution.items())),
                                    "best_certificate": best})
    result["status"] = "COMPLETE"
    result["seconds"] = round(time.monotonic() - started, 6)
    result["sharp_counting_certificate"] = any(x["best_certificate"] is not None
        and x["best_certificate"]["bound"] == len(seeds) for x in result["by_degree"])
    return result


def main():
    panel_path = HERE / "lift_panel.json"
    panel = json.loads(panel_path.read_text())
    output_path = HERE / "lift_panel_results.json"
    plan_hash = hashlib.sha256(panel_path.read_bytes()).hexdigest()
    source_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if output_path.exists():
        report = json.loads(output_path.read_text())
        assert report["panel_sha256"] == plan_hash and report["search_sha256"] == source_hash
    else:
        report = {"panel_sha256": plan_hash, "search_sha256": source_hash,
                  "selection_scope": panel["selection"], "results": [],
                  "claim_boundary": "A hit needs independent certificate replay and prior-resolution audit before a new-knot claim."}
    done = {x["name"] for x in report["results"]}
    for case in panel["cases"]:
        if case["name"] in done:
            continue
        result = test_case(case, panel["search_limits"])
        report["results"].append(result)
        temporary = output_path.with_suffix(".tmp")
        temporary.write_text(json.dumps(report, indent=2) + "\n")
        temporary.replace(output_path)
        print(json.dumps({"name": result["name"], "status": result["status"],
                          "sharp": result.get("sharp_counting_certificate"),
                          "degrees": [{"degree": x["degree"], "bounds": x["surjective_colouring_bound_counts"]}
                                      for x in result["by_degree"]]}), flush=True)


if __name__ == "__main__":
    main()
