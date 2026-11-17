#!/usr/bin/env python3
"""Check whether each retained first A5 core map extends by its deck involution.

This is a structural check on existing witnesses, not a new colouring search.
An implementing double transposition gives an H3 reflection base; an odd
transposition would contradict the separately certified absence of S5 bases.
"""

import hashlib
from itertools import permutations
import json
from pathlib import Path

from parabolic_lifts import oriented_diagram
from search_nonabelian_core import a5_tables

HERE = Path(__file__).resolve().parent


def multiply(a,b):
    return tuple(a[b[i]] for i in range(5))


def inverse(a):
    return tuple(a.index(i) for i in range(5))


def main():
    panel = json.loads((HERE/"parabolic_panel.json").read_text())
    prior = json.loads((HERE/"nonabelian_core_results.json").read_text())
    group,_ = a5_tables()
    identity = tuple(range(5))
    results = []
    for case,old in zip(panel["cases"],prior["results"]):
        assert case["name"] == old["name"]
        witness = old["first_nonconstant_colouring"]
        if witness is None:
            results.append({"name":case["name"],"status":"NO_NONCONSTANT_FIRST_MAP"})
            continue
        images = [group[i] for i in witness["arc_labels"]]
        implementers = [c for c in permutations(range(5))
                        if all(multiply(multiply(c,a),inverse(c)) == inverse(a) for a in images)]
        assert len(implementers) <= 1
        lifts = []
        for c in implementers:
            assert multiply(c,c) == identity
            meridians = [multiply(a,c) for a in images]
            assert all(multiply(m,m) == identity for m in meridians)
            for over,a,b,_ in oriented_diagram(case["pd"]):
                assert meridians[b] == multiply(multiply(meridians[over],meridians[a]),meridians[over])
            parity = sum(c[i]>c[j] for i in range(5) for j in range(i+1,5))%2
            moved = sum(c[i] != i for i in range(5))
            assert parity == 0 and moved == 4
            assert all(sum(m[i] != i for i in range(5)) == 4 for m in meridians)
            lifts.append({"deck_implementer":c,"parity":parity,"moved_points":moved,
                          "meridian_A5_components":meridians,"H3_base_available":True})
        results.append({"name":case["name"],"status":"COMPLETE","implementers":lifts})
    inputs = ["audit_core_involutions.py","parabolic_panel.json","nonabelian_core_results.json",
              "search_nonabelian_core.py","parabolic_lifts.py"]
    report = {"results":results,"input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs},
              "scope":"Exact inversion-equivariance audit on existing first maps only. An H3 base has rank3, not automatically a sharp-four certificate. No alternate witnesses are tested."}
    (HERE/"core_involution_audit.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"targets":len(results),"first_maps":sum(r["status"]=="COMPLETE" for r in results),
                      "H3_bases":sum(len(r.get("implementers",[])) for r in results)},indent=2))


if __name__ == "__main__":
    main()
