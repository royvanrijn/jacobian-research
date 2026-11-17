#!/usr/bin/env python3
"""Replay cover Smith forms and independently build augmentation cocycles.

The Fox-prefix covering matrix is separate from the producer's edge walk.
The four-dimensional cocycle matrix is a third path to the field dimensions.
The PD parser and elementary Smith replay remain shared and disclosed.
"""

from collections import Counter
import hashlib
from itertools import permutations
import json
from math import gcd
from pathlib import Path

from parabolic_lifts import oriented_diagram
from verify_core_subgroup_homology import fox_cover_matrix,generated_group,inverse,multiply
from verify_integer_lifts import nullspace,replay

HERE=Path(__file__).resolve().parent


def prime_factors(value):
    factors=set()
    for divisor in range(2,value+1):
        if divisor*divisor>value:
            break
        if value%divisor==0:
            factors.add(divisor)
            while value%divisor==0:
                value//=divisor
    if value>1:
        factors.add(value)
    return factors


def augmentation_action(permutation):
    # Basis e_i-e_4, i=0,...,3. This does not divide by5.
    return [[int(j==permutation[i])-int(j==permutation[4]) for i in range(4)] for j in range(4)]


def cocycle_matrix(images,words):
    rows=[]
    for word in words:
        prefix=tuple(range(5))
        block=[[0]*(4*len(images)) for _ in range(4)]
        for letter in word:
            i=abs(letter)-1
            if letter<0:
                prefix=multiply(prefix,inverse(images[i]))
            action=augmentation_action(prefix)
            sign=1 if letter>0 else -1
            for a in range(4):
                for b in range(4):
                    block[a][4*i+b]+=sign*action[a][b]
            if letter>0:
                prefix=multiply(prefix,images[i])
        assert prefix==tuple(range(5))
        rows.extend(block)
    principal=[]
    for p in images:
        action=augmentation_action(p)
        principal.extend([[action[a][b]-int(a==b) for b in range(4)] for a in range(4)])
    assert all(sum(row[j]*principal[j][b] for j in range(len(principal)))==0
               for row in rows for b in range(4))
    return rows,principal


def verify_homology(row):
    images=[tuple(p) for p in row["generator_images"]]
    assert len(generated_group(images))==60
    words=row["presentation_words"]
    matrix,columns=fox_cover_matrix(images,words,row["tree_edges"])
    diagonal=replay(matrix,row["smith"])
    free=columns-len(diagonal)
    torsion=[d for d in diagonal if d>1]
    assert free==row["H1_H_free_rank"] and torsion==row["H1_H_torsion_factors"]
    assert columns==row["columns"] and len(matrix)==row["rows"]
    bf,bt=row["H1_G_free_rank"],row["H1_G_torsion_factors"]
    exceptional=sorted(set().union(*(prime_factors(d) for d in torsion+bt))-{5})
    assert [r["prime"] for r in row["exceptional_primes"]]==exceptional
    generic=row["all_other_primes_except5"]
    assert generic["sample_prime"] not in exceptional and generic["sample_prime"]!=5
    assert generic["b1_G"]==bf and generic["b1_H"]==free and generic["gain"]==free-bf
    cocycles,principal=cocycle_matrix(images,words)
    checks=[]
    for expected in row["exceptional_primes"]+[dict(generic,prime=generic["sample_prime"])]:
        prime=expected["prime"]
        assert prime>1 and all(prime%d for d in range(2,prime))
        bg=bf+sum(d%prime==0 for d in bt)
        bh=free+sum(d%prime==0 for d in torsion)
        assert expected["b1_G"]==bg and expected["b1_H"]==bh and expected["gain"]==bh-bg>=0
        z=len(nullspace(cocycles,4*len(images),prime))
        invariants=len(nullspace(principal,4,prime))
        assert invariants==0
        h1=z-(4-invariants)
        assert h1==bh-bg
        assert len(nullspace(matrix,columns,prime))==bh
        checks.append({"prime":prime,"Z1_dimension":z,"principal_dimension":4,"H1_dimension":h1,
                       "b1_G":bg,"b1_H":bh,"gain":bh-bg})
    maximum=max(c["gain"] for c in checks)
    assert maximum==row["maximum_gain_all_primes_except5"]
    assert row["sharp_four_via_augmentation"]==(maximum>=5)
    return {"H1_H_free_rank":free,"H1_H_torsion_factors":torsion,
            "maximum_gain_all_primes_except5":maximum,"field_checks":checks,
            "sharp_four":maximum>=5}


def all_index_five_gate(catalogue,base):
    # Every index5 action has trivial sign since H1(G;Z) has odd order.
    # A proper nontrivial A5 image would supply a forbidden abelian quotient
    # at a prime2,3,5, as independently checked for every subgroup below.
    elements=[tuple(p) for p in catalogue["elements"]]
    index={p:i for i,p in enumerate(elements)}
    for record in catalogue["subgroups"]:
        subgroup=set(record["elements"])
        if len(subgroup) in [1,60]:
            continue
        commutators=[]
        for i in subgroup:
            for j in subgroup:
                comm=multiply(multiply(multiply(elements[i],elements[j]),inverse(elements[i])),inverse(elements[j]))
                commutators.append(index[comm])
        derived={0}
        pending=[0]
        while pending:
            a=pending.pop()
            for b in commutators:
                c=index[multiply(elements[a],elements[b])]
                if c not in derived:
                    derived.add(c)
                    pending.append(c)
        assert derived<subgroup and len(subgroup)%len(derived)==0
        assert prime_factors(len(subgroup)//len(derived))<={2,3,5}
    assert all(r["degrees"][0]["H1_free_rank"]==0 and
               len(r["degrees"][0]["H1_torsion_factors"])==1 and
               gcd(r["degrees"][0]["H1_torsion_factors"][0],30)==1
               for r in base["all_degree_field_panel"])
    return True


def main():
    report=json.loads((HERE/"augmentation_cover_results.json").read_text())
    plan=json.loads((HERE/"augmentation_cover_plan.json").read_text())
    catalogue=json.loads((HERE/"joint_a5_catalogue.json").read_text())
    panel=json.loads((HERE/"joint_a5_panel.json").read_text())
    kernel=json.loads((HERE/"joint_a5_results.json").read_text())
    base=json.loads((HERE/"branched_cover_verification.json").read_text())
    for receipt in [report,plan]:
        for path,digest in receipt["input_sha256"].items():
            assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest
    assert report["status"]=="COMPLETE" and len(report["results"])==16
    assert plan["targets"]==[r["name"] for r in panel["cases"][:16]]
    expected={(r["name"],m["orbit_index"]):(case,r,m) for case,r in zip(panel["cases"][:16],kernel["results"][:16])
              for m in r["epimorphism_representatives"]}
    assert set(expected)=={(r["name"],r["orbit_index"]) for r in report["results"]}
    a5=[p for p in permutations(range(5)) if sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))%2==0]
    assert a5==[tuple(p) for p in catalogue["elements"]]
    base_by_name={r["name"]:r["degrees"][0] for r in base["all_degree_field_panel"]}
    results=[]
    for row in report["results"]:
        case,old,mapping=expected[(row["name"],row["orbit_index"])]
        assert mapping["source_group"]==row["source_group"]==old["source_group"]
        images=[a5[i] for i in mapping["arc_images"]]
        assert [list(p) for p in images]==row["generator_images"]
        root=old["normalized_arc"]
        assert images[root]==tuple(range(5))
        cr=oriented_diagram(case["pd"])
        assert row["presentation_words"]==[[o+1,-(a+1),o+1,-(b+1)] for o,a,b,_ in cr]+[[root+1]]
        assert row["H1_G_free_rank"]==base_by_name[row["name"]]["H1_free_rank"]==0
        assert row["H1_G_torsion_factors"]==base_by_name[row["name"]]["H1_torsion_factors"]
        checked=verify_homology(row)
        assert checked["H1_H_free_rank"]==0
        assert checked["H1_H_torsion_factors"]==[3*row["H1_G_torsion_factors"][0]]
        assert checked["maximum_gain_all_primes_except5"]==1
        assert all(c["gain"]==int(c["prime"]==3) for c in checked["field_checks"])
        results.append({"name":row["name"],"orbit_index":row["orbit_index"],"source_group":row["source_group"],**checked})
    assert len({r["name"] for r in results})==12 and len(report["skipped_no_A5_quotient"])==4
    control_results=[]
    for row in report["controls"]:
        checked=verify_homology(row)
        control_results.append({"group":row["group"],**checked})
    assert [r["maximum_gain_all_primes_except5"] for r in control_results]==[4,2,5]
    assert control_results[0]["H1_H_free_rank"]==6
    assert control_results[2]["H1_H_free_rank"]==2 and control_results[2]["H1_H_torsion_factors"]==[3]*6
    # Failure of the splitting at characteristic5 is an actual control.
    f2=report["controls"][0]
    matrix,principal=cocycle_matrix([tuple(p) for p in f2["generator_images"]],f2["presentation_words"])
    invariants5=len(nullspace(principal,4,5))
    h1v5=len(nullspace(matrix,8,5))-(4-invariants5)
    assert invariants5==1 and h1v5==5 and 6-2==4
    complete_index5=all_index_five_gate(catalogue,base)
    inputs=list(report["input_sha256"])+["augmentation_cover_results.json","verify_augmentation_covers.py",
            "verify_core_subgroup_homology.py","verify_integer_lifts.py"]
    receipt={"status":"PASS","knots":16,"A5_kernel_classes":16,"knots_with_quotients":12,
             "all_index_five_subgroups_covered":complete_index5,
             "skipped_no_index_five_subgroup":report["skipped_no_A5_quotient"],
             "maximum_gain_all_primes_except5":1,"new_sharp_cases":0,"results":results,"controls":control_results,
             "characteristic5_countercontrol":{"V4_invariant_dimension":invariants5,"H1_G_V4":h1v5,"b1_H_minus_b1_G":4},
             "scope":"All connected index5 covers of these16 cover groups are represented: odd H1 kills sign and coprime-to30 H1 excludes every proper A5 image. Other degrees and targets untested. Shared PD parser and Smith replay disclosed; independent Fox-prefix and module cocycle paths.",
             "input_sha256":{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in dict.fromkeys(inputs)}}
    (HERE/"augmentation_cover_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps({k:receipt[k] for k in ["status","knots","A5_kernel_classes","all_index_five_subgroups_covered","maximum_gain_all_primes_except5","new_sharp_cases"]},indent=2))


if __name__=="__main__":
    main()
