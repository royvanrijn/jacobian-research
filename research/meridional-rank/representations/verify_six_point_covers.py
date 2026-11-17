#!/usr/bin/env python3
"""Independent projective-line PSL2(F5) replay of the new degree6 action.

Also checks which cover cohomology classes descend to the knot orbifold.
The PD parser and Smith replay are shared with earlier work and disclosed.
"""

from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path

from parabolic_lifts import oriented_diagram
from verify_core_subgroup_homology import fox_cover_matrix,inverse,multiply
from verify_integer_lifts import nullspace,replay

HERE=Path(__file__).resolve().parent


def prime_factors(value):
    found=set();divisor=2
    while divisor*divisor<=value:
        if value%divisor==0:
            found.add(divisor)
            while value%divisor==0:
                value//=divisor
        divisor+=1
    if value>1:
        found.add(value)
    return found


def projective_model():
    def mm(a,b):
        return ((a[0]*b[0]+a[1]*b[2])%5,(a[0]*b[1]+a[1]*b[3])%5,
                (a[2]*b[0]+a[3]*b[2])%5,(a[2]*b[1]+a[3]*b[3])%5)
    def normal(a):
        return min(a,tuple(-x%5 for x in a))
    identity=normal((1,0,0,1))
    matrices=sorted({normal(a) for a in product(range(5),repeat=4) if (a[0]*a[3]-a[1]*a[2])%5==1})
    matrices.remove(identity);matrices.insert(0,identity)
    index={a:i for i,a in enumerate(matrices)}
    table=[[index[normal(mm(a,b))] for b in matrices] for a in matrices]
    inverses=[next(j for j in range(60) if table[i][j]==0) for i in range(60)]
    action=[]
    for a,b,c,d in matrices:
        images=[]
        for x in range(6):
            numerator,denominator=(a,c) if x==5 else ((a*x+b)%5,(c*x+d)%5)
            images.append(numerator*pow(denominator,-1,5)%5 if denominator else 5)
        action.append(tuple(images))
    assert len(matrices)==len(set(action))==60
    assert all(action[table[a][b]]==multiply(action[a],action[b]) for a in range(60) for b in range(60))
    return table,inverses,action


def coefficient_identification(catalogue,producer_action,table,action):
    elements=[tuple(p) for p in catalogue["elements"]]
    index={p:i for i,p in enumerate(elements)}
    a,b=[index[p] for p in [(1,2,0,3,4),(0,1,3,4,2)]]
    transitions=[[index[multiply(p,elements[g])] for g in [a,b]] for p in elements]
    found=None
    for u,v in product(range(1,60),repeat=2):
        mapping={0:0};queue=[0];failed=False
        while queue and not failed:
            old=queue.pop()
            for g,new_target in zip(transitions[old],[table[mapping[old]][u],table[mapping[old]][v]]):
                if g in mapping:
                    if mapping[g]!=new_target:
                        failed=True;break
                else:
                    mapping[g]=new_target;queue.append(g)
        if not failed and len(mapping)==len(set(mapping.values()))==60:
            found=[mapping[i] for i in range(60)];break
    assert found is not None
    assert all(found[index[multiply(p,q)]]==table[found[i]][found[j]] for i,p in enumerate(elements) for j,q in enumerate(elements))
    sheet_map=None
    for first in range(6):
        mapping={0:first};queue=[0];failed=False
        while queue and not failed:
            point=queue.pop()
            for g in [a,b]:
                old=producer_action[g][point]
                new=action[found[g]][mapping[point]]
                if old in mapping:
                    if mapping[old]!=new:
                        failed=True;break
                else:
                    mapping[old]=new;queue.append(old)
        if not failed and len(mapping)==len(set(mapping.values()))==6:
            sheet_map=[mapping[i] for i in range(6)];break
    assert sheet_map is not None
    assert all(sheet_map[producer_action[g][i]]==action[found[g]][sheet_map[i]] for g in range(60) for i in range(6))
    return found,sheet_map


def action_matrix(p):
    n=len(p)-1
    return [[int(j==p[i])-int(j==p[-1]) for i in range(n)] for j in range(n)]


def cocycle_dimensions(images,words,prime):
    dimension=len(images[0])-1
    rows=[]
    for word in words:
        prefix=tuple(range(dimension+1));block=[[0]*(dimension*len(images)) for _ in range(dimension)]
        for letter in word:
            i=abs(letter)-1
            if letter<0:
                prefix=multiply(prefix,inverse(images[i]))
            a=action_matrix(prefix)
            for j in range(dimension):
                for k in range(dimension):
                    block[j][dimension*i+k]+=(1 if letter>0 else -1)*a[j][k]
            if letter>0:
                prefix=multiply(prefix,images[i])
        assert prefix==tuple(range(dimension+1))
        rows.extend(block)
    principal=[]
    for p in images:
        a=action_matrix(p)
        principal.extend([[a[j][k]-int(j==k) for k in range(dimension)] for j in range(dimension)])
    assert all(sum(row[j]*principal[j][k] for j in range(len(principal)))%prime==0 for row in rows for k in range(dimension))
    invariants=len(nullspace(principal,dimension,prime))
    z=len(nullspace(rows,dimension*len(images),prime))
    return {"prime":prime,"invariants":invariants,"principal_dimension":dimension-invariants,
            "Z1_dimension":z,"H1_dimension":z-dimension+invariants}


def cohomology(images,words,prime):
    dimensions=cocycle_dimensions(images,words,prime)
    assert dimensions["invariants"]==0
    return dimensions["H1_dimension"],dimensions["Z1_dimension"]


def main():
    report=json.loads((HERE/"six_point_cover_results.json").read_text())
    plan=json.loads((HERE/"six_point_cover_plan.json").read_text())
    catalogue=json.loads((HERE/"joint_a5_catalogue.json").read_text())
    panel=json.loads((HERE/"joint_a5_panel.json").read_text())["cases"][:16]
    kernel=json.loads((HERE/"joint_a5_results.json").read_text())["results"][:16]
    for receipt in [report,plan]:
        for path,digest in receipt["input_sha256"].items():
            assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest
    assert report["status"]=="COMPLETE" and len(report["results"])==16
    table,inverses,action=projective_model()
    identification,sheet_map=coefficient_identification(catalogue,report["A5_six_point_action"],table,action)
    expected={(r["name"],m["orbit_index"]):(case,r,m) for case,r in zip(panel,kernel) for m in r["epimorphism_representatives"]}
    assert set(expected)=={(r["name"],r["orbit_index"]) for r in report["results"]}
    results=[]
    for row in report["results"]:
        case,source,mapping=expected[(row["name"],row["orbit_index"])]
        cr=oriented_diagram(case["pd"])
        assert row["presentation_words"]==[[o+1,-(a+1),o+1,-(b+1)] for o,a,b,_ in cr]+[[source["normalized_arc"]+1]]
        assert row["source_group"]==source["source_group"]==mapping["source_group"]
        matrix_indices=[identification[i] for i in mapping["arc_images"]]
        images=[action[i] for i in matrix_indices]
        assert all(images[g][sheet_map[i]]==sheet_map[row["generator_images"][g][i]] for g in range(13) for i in range(6))
        moved_tree=[6*(e//6)+sheet_map[e%6] for e in row["tree_edges"]]
        raw,columns=fox_cover_matrix(images,row["presentation_words"],moved_tree)
        retained=[i for i in range(78) if i not in row["tree_edges"]]
        raw_columns={e:i for i,e in enumerate(e for e in range(78) if e not in moved_tree)}
        ordered=[[raw[6*word+sheet_map[s]][raw_columns[6*(e//6)+sheet_map[e%6]]] for e in retained]
                 for word in range(14) for s in range(6)]
        diagonal=replay(ordered,row["smith"])
        free=columns-len(diagonal);torsion=[d for d in diagonal if d>1]
        assert columns==row["columns"]==73 and len(ordered)==row["rows"]==84
        assert free==row["H1_H_free_rank"]==0 and torsion==row["H1_H_torsion_factors"] and len(torsion)==1
        d=row["H1_G_torsion_factors"][0]
        assert row["H1_G_free_rank"]==0 and torsion[0] in [4*d,28*d] and d%7!=0
        exceptional=(prime_factors(d)|prime_factors(torsion[0]))-{2,3}
        assert {c["prime"] for c in row["field_dimensions"] if not c["generic"]}==exceptional
        generic=[c for c in row["field_dimensions"] if c["generic"]]
        assert len(generic)==1 and generic[0]["prime"] not in exceptional and generic[0]["gain"]==0
        deck=[t for t in range(1,60) if table[t][t]==0 and
              all(table[table[t][g]][t]==inverses[g] for g in matrix_indices)]
        assert len(deck)==1
        knot_images=[action[table[g][deck[0]]] for g in matrix_indices]
        orbifold_words=[[o+1,a+1,o+1,-(b+1)] for o,a,b,_ in cr]+[[i+1,i+1] for i in range(13)]
        checks=[]
        for item in row["field_dimensions"]:
            p=item["prime"]
            assert p>=5 and all(p%k for k in range(2,p))
            bg=sum(v%p==0 for v in row["H1_G_torsion_factors"])
            bh=sum(v%p==0 for v in torsion)
            h1,z=cohomology(images,row["presentation_words"],p)
            assert h1==bh-bg==item["gain"] and bg==item["b1_G"] and bh==item["b1_H"]
            descended,orbifold_z=cohomology(knot_images,orbifold_words,p)
            assert 0<=descended<=h1
            checks.append({**item,"cover_H1_V5":h1,"cover_Z1_V5":z,"orbifold_H1_V5":descended,
                           "orbifold_Z1_V5":orbifold_z})
        assert max(r["gain"] for r in checks)==row["maximum_gain_all_primes_except2_3"]<=1
        assert all(c["gain"]==int(c["prime"]==7 and torsion[0]==28*d) for c in checks)
        assert not row["sharp_four_via_augmentation"]
        modular_v5=[cocycle_dimensions(images,row["presentation_words"],p) for p in [2,3]]
        assert [(c["invariants"],c["H1_dimension"],c["Z1_dimension"]) for c in modular_v5]==[(1,2,6),(1,1,5)]
        original_images=[tuple(catalogue["elements"][i]) for i in mapping["arc_images"]]
        modular_v4=cocycle_dimensions(original_images,row["presentation_words"],5)
        assert (modular_v4["invariants"],modular_v4["H1_dimension"],modular_v4["Z1_dimension"])==(1,1,4)
        results.append({"name":row["name"],"orbit_index":row["orbit_index"],"H1_H_torsion_factors":torsion,
                        "maximum_gain_all_primes_except2_3":row["maximum_gain_all_primes_except2_3"],
                        "deck_element_PSL_index":deck[0],"field_checks":checks,
                        "V5_nonsplit_characteristics":modular_v5,"V4_characteristic5":modular_v4})
    controls=[]
    for row,expected_free,expected_gain in zip(report["controls"],[7,3,7],[5,3,7]):
        images=[tuple(p) for p in row["generator_images"]]
        matrix,columns=fox_cover_matrix(images,row["presentation_words"],row["tree_edges"])
        diagonal=replay(matrix,row["smith"])
        assert columns-len(diagonal)==expected_free and not [d for d in diagonal if d>1]
        assert cohomology(images,row["presentation_words"],5)[0]==expected_gain
        assert row["maximum_gain_all_primes_except2_3"]==expected_gain
        controls.append({"group":row["group"],"cover_free_rank":expected_free,"gain":expected_gain})
    assert Counter(r["maximum_gain_all_primes_except2_3"] for r in results)=={0:8,1:8}
    inputs=list(report["input_sha256"])+["six_point_cover_results.json","verify_six_point_covers.py",
            "verify_core_subgroup_homology.py","verify_integer_lifts.py"]
    receipt={"status":"PASS","covers":16,"coefficient_group":"PSL2(F5) acting on P1(F5)",
             "A5_to_PSL2_identification":identification,"Sylow_action_to_projective_sheet_bijection":sheet_map,
             "results":results,"controls":controls,"maximum_gain":1,"new_sharp_cases":0,
             "V4_V5_cover_cocycle_rank_bound_at_most2_all_fields":True,
             "orbifold_nonprincipal_classes":[{"name":r["name"],"orbit_index":r["orbit_index"],"prime":c["prime"],"dimension":c["orbifold_H1_V5"]}
                                               for r in results for c in r["field_checks"] if c["orbifold_H1_V5"]],
             "scope":"All A5-induced degree6 covers only; A6 and other monodromy not excluded. The Q cocycle calculation tests descent to meridian-involution knot representations, not integrability or varying peripheral curves. Shared PD parser and Smith replay disclosed.",
             "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in dict.fromkeys(inputs)}}
    (HERE/"six_point_cover_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({k:receipt[k] for k in ["status","covers","maximum_gain","new_sharp_cases","orbifold_nonprincipal_classes"]},indent=2))


if __name__=="__main__":
    main()
