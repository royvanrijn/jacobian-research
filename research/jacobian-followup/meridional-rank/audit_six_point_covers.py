#!/usr/bin/env python3
"""Frozen A5 six-point cover audit on every retained first16 kernel class."""

import hashlib
import json
from pathlib import Path
import time

from audit_core_subgroup_homology import covering_matrix
from integer_lifts import smith_certificate
from joint_a5 import group_model
from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent


def prime_divisors(value):
    result=set()
    d=2
    while d*d<=value:
        if value%d==0:
            result.add(d)
            while value%d==0:
                value//=d
        d+=1
    if value>1:
        result.add(value)
    return result


def six_point_action(catalogue):
    elements,table,inverses,orders,_,_=group_model()
    subgroups=sorted(tuple(r["elements"]) for r in catalogue["subgroups"] if len(r["elements"])==5)
    assert len(subgroups)==6
    action=[]
    for g in range(60):
        action.append(tuple(subgroups.index(tuple(sorted(table[table[g][h]][inverses[g]] for h in H))) for H in subgroups))
    assert len(set(action))==60
    assert {p[0] for p in action}==set(range(6))
    assert all(action[table[a][b]]==tuple(action[a][action[b][i]] for i in range(6)) for a in range(60) for b in range(60))
    assert sum(p[0]==0 for p in action)==10
    character=[sum(p[i]==i for i in range(6))-1 for p in action]
    old_character=[sum(p[i]==i for i in range(5))-1 for p in elements]
    assert sum(x*x for x in character)==60 and sum(x*y for x,y in zip(character,old_character))==0
    assert all(character[i]=={1:5,2:1,3:-1,5:0}[orders[i]] for i in range(60))
    return elements,action,subgroups


def homology(images,words,base_free,base_torsion):
    matrix,tree=covering_matrix(images,words)
    columns=6*len(images)-5
    smith=smith_certificate(matrix,columns)
    free=columns-len(smith["diagonal"])
    torsion=[d for d in smith["diagonal"] if d>1]
    primes=sorted(set().union(*(prime_divisors(d) for d in torsion+base_torsion))-{2,3})
    generic=next(p for p in [5,7,11,13,17,19,23,29,31,37,41,43] if p not in primes)
    checks=[]
    for prime in primes+[generic]:
        bg=base_free+sum(d%prime==0 for d in base_torsion)
        bh=free+sum(d%prime==0 for d in torsion)
        assert bh>=bg
        checks.append({"prime":prime,"b1_G":bg,"b1_H":bh,"gain":bh-bg,"generic":prime==generic})
    maximum=max(r["gain"] for r in checks)
    return {"generator_images":images,"presentation_words":words,"tree_edges":tree,
            "rows":len(matrix),"columns":columns,"smith":smith,"H1_H_free_rank":free,
            "H1_H_torsion_factors":torsion,"H1_G_free_rank":base_free,"H1_G_torsion_factors":base_torsion,
            "field_dimensions":checks,"maximum_gain_all_primes_except2_3":maximum,
            "sharp_four_via_augmentation":maximum>=6}


def main():
    source_paths=["joint_a5_catalogue.json","joint_a5_panel.json","joint_a5_results.json","joint_a5_verification.json",
                  "parabolic_panel.json","branched_cover_verification.json","augmentation_cover_verification.json"]
    sources={p:json.loads((HERE/p).read_text()) for p in source_paths}
    panel=sources["joint_a5_panel.json"]["cases"][:16]
    kernel=sources["joint_a5_results.json"]["results"][:16]
    assert [c["name"] for c in panel]==[c["name"] for c in sources["parabolic_panel.json"]["cases"]]
    selected=[{"name":r["name"],"orbit_indices":[m["orbit_index"] for m in r["epimorphism_representatives"]]} for r in kernel]
    assert sum(len(r["orbit_indices"]) for r in selected)==16
    plan={"status":"FROZEN_BEFORE_SIX_POINT_AUDIT","selection":selected,
          "action":"Conjugation of A5 on its six Sylow5 subgroups; point stabilizer is D5 of order10",
          "motivation":"The dimension5 augmentation character is irreducible over C and orthogonal to the old dimension4 character. Index6 allows characteristic5, unavailable for the old split index5 module. This tests a structurally different coefficient module, not a larger unmotivated degree sweep.",
          "sharp_four_gain_threshold":6,"characteristics":"All primes other than2,3, by integral Smith data",
          "limits":{"covers":16,"seconds_total":60,"seconds_per_cover":10,"rows":84,"columns":73,"checkpoint_after_each_cover":True},
          "controls":["F2 index6 F7","C3*C3 index6 F3","C3*C3*C3 index6 F7"],
          "scope":"All A5-induced index6 covers on these same targets. No claim about other degree6 monodromy such as A6. No adaptive replacements.",
          "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in source_paths}}
    path=HERE/"six_point_cover_plan.json"
    if not path.exists():
        path.write_text(json.dumps(plan,indent=2)+"\n")
        print("FROZEN six_point_cover_plan.json; rerun to execute")
        return
    assert json.loads(path.read_text())==plan
    inputs=source_paths+["audit_six_point_covers.py","six_point_cover_plan.json","audit_core_subgroup_homology.py",
                         "integer_lifts.py","joint_a5.py","parabolic_lifts.py"]
    hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in inputs}
    output=HERE/"six_point_cover_results.json"
    report=json.loads(output.read_text()) if output.exists() else {"status":"RUNNING","results":[],"input_sha256":hashes}
    assert report["input_sha256"]==hashes
    elements,action,subgroups=six_point_action(sources["joint_a5_catalogue.json"])
    base={r["name"]:r["degrees"][0] for r in sources["branched_cover_verification.json"]["all_degree_field_panel"]}
    done={(r["name"],r["orbit_index"]) for r in report["results"]}
    start=time.monotonic()
    for case,old in zip(panel,kernel):
        assert case["name"]==old["name"]
        words=[[o+1,-(a+1),o+1,-(b+1)] for o,a,b,_ in oriented_diagram(case["pd"])]+[[old["normalized_arc"]+1]]
        for mapping in old["epimorphism_representatives"]:
            if (case["name"],mapping["orbit_index"]) in done:
                continue
            assert time.monotonic()-start<60
            before=time.monotonic()
            images=[action[i] for i in mapping["arc_images"]]
            row=homology(images,words,base[case["name"]]["H1_free_rank"],base[case["name"]]["H1_torsion_factors"])
            assert row["rows"]==84 and row["columns"]==73 and time.monotonic()-before<10
            row.update({"name":case["name"],"source_group":mapping["source_group"],"orbit_index":mapping["orbit_index"]})
            report["results"].append(row)
            output.write_text(json.dumps(report,indent=2)+"\n")
            print(json.dumps({k:row[k] for k in ["name","orbit_index","H1_H_free_rank","H1_H_torsion_factors","maximum_gain_all_primes_except2_3","sharp_four_via_augmentation"]}),flush=True)
    index={tuple(p):i for i,p in enumerate(elements)}
    images=[action[index[p]] for p in [(1,2,0,3,4),(0,1,3,4,2)]]
    controls=[]
    for label,im,words,bf,bt in [("F2",images,[],2,[]),("C3*C3",images,[[1]*3,[2]*3],0,[3,3]),
                               ("C3*C3*C3",images+[images[0]],[[1]*3,[2]*3,[3]*3],0,[3,3,3])]:
        controls.append({"group":label,**homology(im,words,bf,bt)})
    assert len(report["results"])==16
    report.update({"status":"COMPLETE","six_sylow5_subgroups":subgroups,"A5_six_point_action":action,
                   "controls":controls,"character_inner_products":{"V5_with_itself":1,"V5_with_V4":0},
                   "skipped_no_A5_map":[r["name"] for r in selected if not r["orbit_indices"]]})
    output.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":"COMPLETE","covers":16,"new_sharp_cases":sum(r["sharp_four_via_augmentation"] for r in report["results"])}))


if __name__=="__main__":
    main()
