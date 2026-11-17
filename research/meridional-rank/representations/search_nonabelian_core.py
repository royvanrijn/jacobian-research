#!/usr/bin/env python3
"""Fixed A5 core-quandle audit on the unchanged current13 panel.

First invocation freezes the complete plan; rerun executes it.
All accepted normalized colourings count Hom(pi1(Sigma2(K)),A5)
by the classical core-group theorem, whose credit is retained.
"""

import hashlib
from itertools import permutations, product
import json
from pathlib import Path
import time

from parabolic_lifts import oriented_diagram

HERE = Path(__file__).resolve().parent


def a5_tables():
    group = [p for p in permutations(range(5)) if sum(p[i]>p[j] for i in range(5) for j in range(i+1,5)) % 2 == 0]
    identity = tuple(range(5))
    assert group[0] == identity and len(group) == 60
    index = {p:i for i,p in enumerate(group)}
    def mul(a,b):
        return tuple(a[b[j]] for j in range(5))
    inverses = [index[tuple(p.index(j) for j in range(5))] for p in group]
    core = [[index[mul(mul(b,group[inverses[i]]),b)] for i in range(60)] for b in group]
    assert all(core[b][core[b][a]] == a for a in range(60) for b in range(60))
    return group,core


def main():
    panel = json.loads((HERE/"parabolic_panel.json").read_text())
    prior = json.loads((HERE/"parabolic_panel_results.json").read_text())
    plan = {"targets":[c["name"] for c in panel["cases"]],
            "source_panel_sha256":hashlib.sha256((HERE/"parabolic_panel.json").read_bytes()).hexdigest(),
            "coefficient_group":"A5, even permutations of five letters", "group_order":60,
            "quandle":"Core(A5), incoming a under overstrand b maps to b*a^(-1)*b",
            "normalized_first_seed":"identity", "assignments_per_case":60**3,
            "max_seconds_per_case":20,"no_adaptive_replacements":True,
            "motivation":"The all-field abelian lift bound stops at2 on these16 cases. The nonabelian core-group theorem identifies normalized colourings with Hom(pi1(double branched cover),A5). More than60^2 such homomorphisms forces cover group rank>=3, hence meridional rank>=4, matching the retained four-seed upper bound. This is a different coefficient group, not a larger prime sweep.",
            "sharp_threshold_strictly_greater_than":60**2,
            "claim_boundary":"A positive needs independent group/count/seed/name replay and prior-resolution audit. A no-hit excludes only this coefficient group and does not bound all nonabelian quotients."}
    plan_path=HERE/"nonabelian_core_plan.json"
    if not plan_path.exists():
        plan_path.write_text(json.dumps(plan,indent=2)+"\n")
        print("FROZEN nonabelian_core_plan.json; rerun to execute")
        return
    assert json.loads(plan_path.read_text())==plan
    inputs=["nonabelian_core_plan.json","parabolic_panel.json","parabolic_panel_results.json",
            "search_nonabelian_core.py","parabolic_lifts.py"]
    hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}
    out=HERE/"nonabelian_core_results.json"
    if out.exists():
        report=json.loads(out.read_text())
        assert report["input_sha256"]==hashes
    else:
        report={"input_sha256":hashes,"results":[],"claim_boundary":plan["claim_boundary"]}
    group,table=a5_tables()
    done={r["name"] for r in report["results"]}
    for case,old in zip(panel["cases"],prior["results"]):
        if case["name"] in done:
            continue
        started=time.monotonic()
        crossings=[c[:3] for c in oriented_diagram(case["pd"])]
        count,tested,best=0,0,None
        complete=True
        for tail in product(range(60),repeat=3):
            if time.monotonic()-started>plan["max_seconds_per_case"]:
                complete=False
                break
            labels=[None]*len(crossings)
            for s,value in zip(old["seeds"],(0,)+tail):
                labels[s]=value
            for o,a,b in old["seed_trace"]:
                labels[b]=table[labels[o]][labels[a]]
            tested+=1
            if all(labels[b]==table[labels[o]][labels[a]] for o,a,b in crossings):
                count+=1
                if best is None and any(labels):
                    best={"seed_assignment":(0,)+tail,"arc_labels":labels}
        result={"name":case["name"],"status":"COMPLETE" if complete else "INCOMPLETE",
                "normalized_assignments_tested":tested,"normalized_colourings":count,
                "sharp_four_count_certificate":count>60**2,
                "first_nonconstant_colouring":best,"seeds":old["seeds"],"seed_trace":old["seed_trace"],
                "seconds":round(time.monotonic()-started,6)}
        report["results"].append(result)
        out.write_text(json.dumps(report,indent=2)+"\n")
        print(json.dumps({"name":case["name"],"status":result["status"],"normalized_colourings":count,
                          "sharp_four":result["sharp_four_count_certificate"]}),flush=True)


if __name__=="__main__":
    main()
