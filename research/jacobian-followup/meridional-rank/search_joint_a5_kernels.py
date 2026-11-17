#!/usr/bin/env python3
"""Exact per-knot A5 quotient-kernel audit, checkpointed after every case."""

from collections import Counter
import hashlib
import json
from pathlib import Path
import time

from joint_a5 import group_model,core_extension,certify_distinct_same_source
from parabolic_lifts import oriented_diagram
from search_lift_panel import seed_certificate

HERE=Path(__file__).resolve().parent


def main():
    panel=json.loads((HERE/"joint_a5_panel.json").read_text())
    for path,digest in panel["input_sha256"].items():
        assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest
    catalogue=json.loads((HERE/"joint_a5_catalogue.json").read_text())
    elements,table,inverse,orders,automorphisms,core=group_model()
    assert [list(p) for p in elements]==catalogue["elements"]
    triples=catalogue["generating_triple_orbits"]
    old={r["name"]:r for r in json.loads((HERE/"nonabelian_core_results.json").read_text())["results"]}
    inputs=["search_joint_a5_kernels.py","joint_a5.py","joint_a5_panel.json","joint_a5_catalogue.json",
            "parabolic_lifts.py","search_lift_panel.py","nonabelian_core_results.json"]
    hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}
    path=HERE/"joint_a5_results.json"
    if path.exists():
        report=json.loads(path.read_text())
        assert report["input_sha256"]==hashes
    else:
        report={"status":"IN_PROGRESS","input_sha256":hashes,"results":[],"claim_boundary":panel["claim_boundary"]}
    done={r["name"] for r in report["results"]}
    started=time.monotonic()
    for case in panel["cases"]:
        if case["name"] in done:
            continue
        if time.monotonic()-started>panel["limits"]["seconds_per_run"]:
            report["status"]="CAPPED_INCOMPLETE"
            break
        case_started=time.monotonic()
        deadline=case_started+panel["limits"]["seconds_per_case"]
        crossings=[c[:3] for c in oriented_diagram(case["pd"])]
        n=len(crossings)
        seeds,trace,seed_trials=seed_certificate(n,crossings,4,panel["limits"]["seed_trials_per_case"],deadline)
        if seeds is None:
            result={"name":case["name"],"status":"NO_FOUR_SEED_CERTIFICATE_WITHIN_CAP","seed_trials":seed_trials}
        else:
            source_hash=hashlib.sha256(json.dumps({"name":case["name"],"pd":case["pd"],"normalized_arc":seeds[0]},sort_keys=True).encode()).hexdigest()
            source_group=case["name"]+":double-cover:"+source_hash
            maps=[]
            tested=0
            for index,values in enumerate(triples):
                if time.monotonic()>deadline:
                    break
                labels=core_extension(n,crossings,seeds,trace,values,core)
                tested+=1
                if labels is not None:
                    maps.append({"source_group":source_group,"orbit_index":index,
                                 "generator_images":values,"arc_images":labels})
            count=len(maps)
            result={"name":case["name"],"status":"COMPLETE" if tested==1668 else "CAPPED_INCOMPLETE",
                    "source_group":source_group,"normalized_arc":seeds[0],"seeds":seeds,"seed_trace":trace,
                    "seed_trials":seed_trials,"orbit_trials":tested,"A5_quotient_kernels":count,
                    "epimorphism_representatives":maps,"sharp_four":count>=20,
                    "seconds":round(time.monotonic()-case_started,6)}
            if case["name"] in old and tested==1668:
                assert old[case["name"]]["normalized_colourings"]==1+120*count
                result["old_first16_Hom_count_consistent"]=True
            if count>=20:
                assert certify_distinct_same_source(maps[:20],source_group,table,automorphisms)==20
        report["results"].append(result)
        path.write_text(json.dumps(report,indent=2)+"\n")
        if result.get("sharp_four") or len(report["results"])%25==0:
            print(json.dumps({"completed":len(report["results"]),"name":case["name"],"kernels":result.get("A5_quotient_kernels"),"sharp_four":result.get("sharp_four",False)}),flush=True)
    else:
        report["status"]="COMPLETE" if all(r["status"]=="COMPLETE" for r in report["results"]) else "INCOMPLETE_CASES"
    report["counts_per_knot_histogram"]=dict(sorted(Counter(r["A5_quotient_kernels"] for r in report["results"] if r["status"]=="COMPLETE").items()))
    report["sharp_four_names"]=[r["name"] for r in report["results"] if r.get("sharp_four")]
    path.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":report["status"],"cases":len(report["results"]),
                      "kernel_counts_per_knot_histogram":report["counts_per_knot_histogram"],
                      "sharp_four_names":report["sharp_four_names"]},indent=2))


if __name__=="__main__":
    main()
