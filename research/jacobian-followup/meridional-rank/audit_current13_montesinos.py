#!/usr/bin/env python3
"""Execute the frozen current13 three/four-tangle diagram identity audit."""

from collections import Counter
import hashlib
import importlib.metadata as metadata
import json
from pathlib import Path
import time

from audit_diagram_identity import RationalTangle,normalize,pd_walk

HERE=Path(__file__).resolve().parent


def main():
    plan=json.loads((HERE/"current13_montesinos_plan.json").read_text())
    for path,digest in plan["input_sha256"].items():
        assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest
    panel=json.loads((HERE/plan["source_panel"]).read_text())
    assert [c["name"] for c in panel["cases"]]==plan["targets"]
    versions={name:metadata.version(name) for name in ["snappy","spherogram","snappy_15_knots","snappy_manifolds"]}
    assert versions=={"snappy":"3.3.2","spherogram":"2.4.1","snappy_15_knots":"1.2.1","snappy_manifolds":"1.4"}
    keys={}
    for case in panel["cases"]:
        key=normalize(pd_walk(case["pd"]),True)[0]
        keys.setdefault(key,[]).append(case["name"])
    inputs=["audit_current13_montesinos.py","current13_montesinos_plan.json",plan["source_panel"],
            plan["source_determinants"],"audit_diagram_identity.py","diagram_audit_requirements.txt"]
    hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}
    path=HERE/"current13_montesinos_results.json"
    if path.exists():
        report=json.loads(path.read_text())
        assert report["input_sha256"]==hashes
    else:
        report={"status":"IN_PROGRESS","input_sha256":hashes,"versions":versions,
                "completed_candidates":0,"diagrams_checked":0,"knot_diagrams":0,"matches":{},
                "claim_boundary":plan["claim_boundary"]}
    started=time.monotonic()
    path.write_text(json.dumps(report,indent=2)+"\n")
    for index,candidate in enumerate(plan["candidates"]):
        if index<report["completed_candidates"]:
            continue
        if time.monotonic()-started>plan["limits"]["maximum_run_seconds"]:
            report["status"]="CAPPED_INCOMPLETE"
            break
        for mask in range(2**candidate["branches"]):
            tangle=RationalTangle(candidate["integer_twists"])
            for i,(a,b) in enumerate(candidate["fractions"]):
                term=RationalTangle(a,b)
                if mask>>i&1:
                    term=term.rotate(2)
                tangle=tangle+term
            link=tangle.numerator_closure()
            assert len(link.crossings)==13
            report["diagrams_checked"]+=1
            if len(link.link_components)!=1:
                continue
            report["knot_diagrams"]+=1
            pd=[list(c) for c in link.PD_code()]
            key,change=normalize(pd_walk(pd),True)
            for name in keys.get(key,[]):
                if name not in report["matches"]:
                    report["matches"][name]={"candidate_index":index,"half_turn_mask":mask,
                                             **candidate,"pd":pd,"full_transform":change}
        report["completed_candidates"]=index+1
        if (index+1)%plan["limits"]["checkpoint_every_candidates"]==0:
            path.write_text(json.dumps(report,indent=2)+"\n")
            print(json.dumps({"completed":index+1,"diagrams":report["diagrams_checked"],"matches":list(report["matches"])}),flush=True)
    else:
        report["status"]="COMPLETE"
    report["elapsed_this_run_seconds"]=round(time.monotonic()-started,6)
    path.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":report["status"],"candidates":report["completed_candidates"],"diagrams":report["diagrams_checked"],
                      "matches":[{"name":name,"branches":c["branches"],"fractions":c["fractions"],"integer":c["integer_twists"]}
                                 for name,c in sorted(report["matches"].items())]},indent=2))


if __name__=="__main__":
    main()
