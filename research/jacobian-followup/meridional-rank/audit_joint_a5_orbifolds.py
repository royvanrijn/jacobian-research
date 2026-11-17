#!/usr/bin/env python3
"""Frozen four-cone orbifold controls for same-source A5 kernel bounds.

Every epimorphism of a Seifert double-cover group to centerless A5 kills
the central fibre. Thus these are per-group counts, never pooled knots.
"""

import hashlib
from itertools import combinations_with_replacement
import json
from pathlib import Path

from joint_a5 import group_model,certify_distinct_same_source

HERE=Path(__file__).resolve().parent


def main():
    patterns=[[even,*odd] for even in [2,6,10,30] for odd in combinations_with_replacement([3,5,15],3)]
    plan={"status":"FROZEN_BEFORE_ORBIFOLD_COUNTS","patterns":patterns,"coefficient_group":"A5",
          "presentation":"Gamma(a,b,c,d)=<x1,x2,x3,x4 | x1^a=x2^b=x3^c=x4^d=x1*x2*x3*x4=1>",
          "selection":"All40 denominator-divisibility types with one even entry and three odd entries divisible by3 or5; no pattern replacement",
          "motivation":"Find an actual known Montesinos/pretzel family demonstrating the user's same-group20-kernel criterion. A5 element orders are1,2,3,5, so the listed types cover every nontrivial odd-cone option and every even-cone option relevant to a four-cone epimorphism.",
          "limits":{"patterns":40,"generating_triple_orbits_per_pattern":1668,"maximum_checks":66720,"no_adaptive_replacements":True},
          "input_sha256":{"joint_a5_catalogue.json":hashlib.sha256((HERE/"joint_a5_catalogue.json").read_bytes()).hexdigest()},
          "boundary":"These are different source orbifold groups with separate counts. A successful pattern is pulled back separately along one knot cover's quotient to that orbifold. The known Montesinos family theorem already settles the corresponding MRC cases; no novel knot resolution is asserted."}
    path=HERE/"joint_a5_orbifold_plan.json"
    if not path.exists():
        path.write_text(json.dumps(plan,indent=2)+"\n")
        print("FROZEN joint_a5_orbifold_plan.json; rerun to execute")
        return
    assert json.loads(path.read_text())==plan
    catalogue=json.loads((HERE/"joint_a5_catalogue.json").read_text())
    elements,table,inverses,orders,automorphisms,core=group_model()
    results=[]
    for pattern in patterns:
        source="orbifold-Gamma("+",".join(map(str,pattern))+")"
        accepted=[]
        for index,values in enumerate(catalogue["generating_triple_orbits"]):
            fourth=inverses[table[table[values[0]][values[1]]][values[2]]]
            if all(exponent%orders[value]==0 for exponent,value in zip(pattern,[*values,fourth])):
                accepted.append({"source_group":source,"orbit_index":index,
                                 "generator_images":values,"fourth_image":fourth})
        if len(accepted)>=20:
            assert certify_distinct_same_source(accepted[:20],source,table,automorphisms)==20
        results.append({"cone_orders":pattern,"source_group":source,"A5_quotient_kernels":len(accepted),
                        "rank_at_least3":len(accepted)>=20,"representatives":accepted})
    inputs=["audit_joint_a5_orbifolds.py","joint_a5_orbifold_plan.json","joint_a5_catalogue.json","joint_a5.py"]
    output={"status":"COMPLETE","results":results,"input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs},
            "scope":plan["boundary"]}
    (HERE/"joint_a5_orbifold_results.json").write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps({"patterns":40,"at_least20_kernel_patterns":sum(r["rank_at_least3"] for r in results),
                      "counts":[{"cone_orders":r["cone_orders"],"kernels":r["A5_quotient_kernels"]} for r in results]},indent=2))


if __name__=="__main__":
    main()
