#!/usr/bin/env python3
"""Freeze then test all A5 kernel classes on the retained sixteen-knot panel.

The new rank gate uses the homology gain, rather than Schreier's absolute
homology threshold. Old maps and homology evidence are preserved unchanged.
"""

import hashlib
import json
from pathlib import Path
import time

from audit_core_subgroup_homology import covering_matrix
from integer_lifts import smith_certificate
from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent


def prime_divisors(value):
    result=[]
    divisor=2
    while divisor*divisor<=value:
        if value%divisor==0:
            result.append(divisor)
            while value%divisor==0:
                value//=divisor
        divisor+=1
    if value>1:
        result.append(value)
    return result


def gains(base_free,base_torsion,cover_free,cover_torsion):
    primes=sorted({p for d in base_torsion+cover_torsion for p in prime_divisors(d)}-{5})
    generic=next(p for p in [2,3,7,11,13,17,19,23,29,31,37,41,43] if p not in primes)
    rows=[]
    for prime in primes:
        base=base_free+sum(d%prime==0 for d in base_torsion)
        cover=cover_free+sum(d%prime==0 for d in cover_torsion)
        assert cover>=base
        rows.append({"prime":prime,"b1_G":base,"b1_H":cover,"gain":cover-base})
    generic_gain=cover_free-base_free
    assert generic_gain>=0
    maximum=max([generic_gain]+[r["gain"] for r in rows])
    return {"exceptional_primes":rows,"all_other_primes_except5":{"sample_prime":generic,"b1_G":base_free,
            "b1_H":cover_free,"gain":generic_gain},"maximum_gain_all_primes_except5":maximum,
            "sharp_four_via_augmentation":maximum>=5}


def homology(images,words,base_free,base_torsion):
    matrix,tree=covering_matrix(images,words)
    columns=len(images)*5-4
    smith=smith_certificate(matrix,columns)
    free=columns-len(smith["diagonal"])
    torsion=[d for d in smith["diagonal"] if d>1]
    return {"generator_images":images,"presentation_words":words,"tree_edges":tree,
            "rows":len(matrix),"columns":columns,"smith":smith,"H1_H_free_rank":free,
            "H1_H_torsion_factors":torsion,"H1_G_free_rank":base_free,"H1_G_torsion_factors":base_torsion,
            **gains(base_free,base_torsion,free,torsion)}


def main():
    panel=json.loads((HERE/"joint_a5_panel.json").read_text())
    kernels=json.loads((HERE/"joint_a5_results.json").read_text())
    old_panel=json.loads((HERE/"parabolic_panel.json").read_text())
    catalogue=json.loads((HERE/"joint_a5_catalogue.json").read_text())
    base=json.loads((HERE/"branched_cover_verification.json").read_text())
    names=[c["name"] for c in old_panel["cases"]]
    assert names==[c["name"] for c in panel["cases"][:16]]
    source_inputs=["joint_a5_panel.json","joint_a5_results.json","joint_a5_verification.json",
                   "joint_a5_catalogue.json","parabolic_panel.json","branched_cover_verification.json"]
    selected=[{"name":r["name"],"source_group":r["source_group"],
               "orbit_indices":[m["orbit_index"] for m in r["epimorphism_representatives"]]}
              for r in kernels["results"][:16]]
    assert len(selected)==16 and sum(len(r["orbit_indices"]) for r in selected)==16
    plan={"status":"FROZEN_BEFORE_AUGMENTATION_TEST","targets":names,"selected_kernels":selected,
          "selection":"All A5 quotient-kernel classes on exactly the retained first16 knots, with no adaptive map or target replacement",
          "motivation":"For p!=5, b1(H;Fp)-b1(G;Fp)>=5 forces rank G>=3 via the four-dimensional augmentation module; this is stronger than the previous absolute Schreier threshold7",
          "cover_degree":5,"maximum_covers":16,"maximum_matrix_rows":70,"maximum_matrix_columns":61,
          "limits":{"total_seconds":60,"seconds_per_cover":10,"checkpoint_after_each_cover":True},
          "characteristics":"All primes except5, classified by integral Smith factors; no blind prime sweep",
          "controls":["F2","C3*C3","C3*C3*C3"],
          "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in source_inputs},
          "scope":"One representative per Aut(A5) orbit suffices: all order12 A4 subgroups are conjugate in A5 and Aut(A5) permutes them. This does not test other finite targets, actions or cover degrees."}
    plan_path=HERE/"augmentation_cover_plan.json"
    if not plan_path.exists():
        plan_path.write_text(json.dumps(plan,indent=2)+"\n")
        print("FROZEN augmentation_cover_plan.json; rerun to execute")
        return
    assert json.loads(plan_path.read_text())==plan
    inputs=source_inputs+["audit_augmentation_covers.py","augmentation_cover_plan.json",
                         "audit_core_subgroup_homology.py","integer_lifts.py","parabolic_lifts.py"]
    hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}
    output=HERE/"augmentation_cover_results.json"
    if output.exists():
        report=json.loads(output.read_text())
        assert report["input_sha256"]==hashes
    else:
        report={"status":"RUNNING","results":[],"input_sha256":hashes}
    done={(r["name"],r["orbit_index"]) for r in report["results"]}
    base_by_name={r["name"]:r["degrees"][0] for r in base["all_degree_field_panel"]}
    started=time.monotonic()
    for case,result in zip(panel["cases"][:16],kernels["results"][:16]):
        assert case["name"]==result["name"]
        crossings=oriented_diagram(case["pd"])
        root=result["normalized_arc"]
        words=[[o+1,-(a+1),o+1,-(b+1)] for o,a,b,_ in crossings]+[[root+1]]
        for mapping in result["epimorphism_representatives"]:
            key=(case["name"],mapping["orbit_index"])
            if key in done:
                continue
            assert time.monotonic()-started<60
            before=time.monotonic()
            images=[catalogue["elements"][i] for i in mapping["arc_images"]]
            row=base_by_name[case["name"]]
            data=homology(images,words,row["H1_free_rank"],row["H1_torsion_factors"])
            assert time.monotonic()-before<10 and data["rows"]==70 and data["columns"]==61
            data.update({"name":case["name"],"orbit_index":mapping["orbit_index"],"source_group":mapping["source_group"]})
            report["results"].append(data)
            output.write_text(json.dumps(report,indent=2)+"\n")
            print(json.dumps({k:data[k] for k in ["name","orbit_index","H1_H_torsion_factors","maximum_gain_all_primes_except5","sharp_four_via_augmentation"]}),flush=True)
    images=[(1,2,0,3,4),(0,1,3,4,2)]
    controls=[]
    for label,control_images,words,bf,bt in [
        ("F2",images,[],2,[]),
        ("C3*C3",images,[[1]*3,[2]*3],0,[3,3]),
        ("C3*C3*C3",images+[images[0]],[[1]*3,[2]*3,[3]*3],0,[3,3,3])]:
        controls.append({"group":label,**homology(control_images,words,bf,bt)})
    assert len(report["results"])==16
    report.update({"status":"COMPLETE","controls":controls,
                   "skipped_no_A5_quotient":[r["name"] for r in selected if not r["orbit_indices"]]})
    output.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":"COMPLETE","covers":16,"new_sharp_cases":sum(r["sharp_four_via_augmentation"] for r in report["results"])}))


if __name__=="__main__":
    main()
