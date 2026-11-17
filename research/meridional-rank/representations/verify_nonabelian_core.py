#!/usr/bin/env python3
"""Replay A5-core counts using PSL2(F5), without permutation arithmetic."""

from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path

from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent


def psl_tables():
    def canonical(a):
        return min(tuple(x%5 for x in a),tuple(-x%5 for x in a))
    identity=canonical((1,0,0,1))
    elements=sorted({canonical(a) for a in product(range(5),repeat=4) if (a[0]*a[3]-a[1]*a[2])%5==1})
    elements.remove(identity)
    elements.insert(0,identity)
    assert len(elements)==60
    index={a:i for i,a in enumerate(elements)}
    def multiply(a,b):
        return canonical((a[0]*b[0]+a[1]*b[2],a[0]*b[1]+a[1]*b[3],
                          a[2]*b[0]+a[3]*b[2],a[2]*b[1]+a[3]*b[3]))
    table=[[index[multiply(a,b)] for b in elements] for a in elements]
    inverse=[index[canonical((a[3],-a[1],-a[2],a[0]))] for a in elements]
    assert all(table[i][inverse[i]]==table[inverse[i]][i]==0 for i in range(60))
    orders=[]
    for i in range(60):
        value=i
        order=1
        while value:
            value=table[value][i]
            order+=1
        orders.append(order)
    assert Counter(orders)=={1:1,2:15,3:20,5:24}
    core=[[table[table[b][inverse[a]]][b] for a in range(60)] for b in range(60)]
    # Classical trefoil normalized core colourings: arbitrary order-3 image.
    assert sum(core[0][core[i][0]]==i for i in range(60))==21
    return core


def main():
    panel=json.loads((HERE/"parabolic_panel.json").read_text())
    report=json.loads((HERE/"nonabelian_core_results.json").read_text())
    for path,digest in report["input_sha256"].items():
        assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest
    core=psl_tables()
    results=[]
    for case,result in zip(panel["cases"],report["results"]):
        assert result["name"]==case["name"] and result["status"]=="COMPLETE"
        crossings=[c[:3] for c in oriented_diagram(case["pd"])]
        seen=set(result["seeds"])
        assert len(seen)==4
        for o,a,b in result["seed_trace"]:
            assert o in seen and a in seen and b not in seen
            assert (o,a,b) in crossings or (o,b,a) in crossings
            seen.add(b)
        assert len(seen)==len(crossings)
        count=0
        for values in product(range(60),repeat=3):
            labels={s:a for s,a in zip(result["seeds"],(0,)+values)}
            valid=True
            for o,a,b in result["seed_trace"]:
                labels[b]=core[labels[o]][labels[a]]
                if any(labels[y]!=core[labels[x]][labels[z]] for x,z,y in crossings if x in labels and z in labels and y in labels):
                    valid=False
                    break
            if valid:
                count+=1
        assert count==result["normalized_colourings"] and count<=3600
        assert result["normalized_assignments_tested"]==60**3
        results.append({"name":case["name"],"normalized_A5_core_colourings":count,"sharp_four":False})
    inputs=["verify_nonabelian_core.py","nonabelian_core_results.json","nonabelian_core_plan.json","parabolic_panel.json","parabolic_lifts.py"]
    receipt={"status":"PASS","targets":len(results),"results":results,
             "independent_coefficient_group":"PSL2(F5) matrices modulo scalar +/-I; order and element-order distribution checked",
             "known_control":"Trefoil:21 normalized core colourings = Hom(C3,A5)",
             "scope":"Separate matrix coefficient group and partial-consistency propagation reproduce counts. PD parser is shared. No new sharp case and no exclusion of other coefficient groups.",
             "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}}
    (HERE/"nonabelian_core_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({"status":"PASS","targets":len(results),"counts":dict(Counter(r["normalized_A5_core_colourings"] for r in results)),"new_sharp_cases":0},indent=2))


if __name__=="__main__":
    main()
