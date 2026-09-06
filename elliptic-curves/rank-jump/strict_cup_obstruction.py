#!/usr/bin/env python3
"""Finite cup-action verification and retained strict CT detection of ideal obstructions."""
import argparse
from itertools import permutations
from pathlib import Path
import retrospective as r
import local_collision as lc
import torsion_difference as td
import strict_deformation_solubility as strict
import strict_Sha_Artin as art
import affine_selmer as af
from cubic_bridge import Cubic

PROTOCOL=Path(__file__).with_name("STRICT_CUP_OBSTRUCTION_PROTOCOL.json")
OUTPUT=r.OUT/"rank_jump_strict_cup_obstruction_v1.json"


def embed(v):
    return r.pack(td.weil(t,v) for t in (1,2,3))


def projection(v):
    return v^(7 if v.bit_count()%2 else 0)


def finite_check():
    checks=[]
    for perm in permutations(range(3)):
        g=(perm[0]+1,perm[1]+1)
        for signs in range(8):
            for v in range(4):
                gv=td.act(g,v)
                actual=embed(td.act(td.sign_matrix(signs),gv))
                expected=projection(signs&embed(gv))
                assert actual==expected
                checks.append({"permutation":list(perm),"gamma_sign_mask":signs,
                               "input_V":v,"cup_product_projection":expected,
                               "Jacobian_connecting_action":actual})
    assert len(checks)==192
    for v in range(4):
        assert embed(v).bit_count()%2==0
        for w in range(4):
            assert (embed(v)&embed(w)).bit_count()%2==td.weil(v,w)
    return checks


def calculate():
    local=r.read(lc.INPUT);s=r.read(strict.OUTPUT)
    A,B=map(r.F,local["anchor"]["short_model_ainvariants"][3:])
    K=Cubic(A,B);cases=[]
    for u in r.read(PROTOCOL)["parameters"]:
        old=next(x for x in s["single_deformations"] if x["u"]==u)
        raw=next(x for x in r.read(af.INPUT)["cases"] if x["u"]==u)
        gamma=K.sub(K.one,K.scale(K.theta,u));D=1+A*u*u+B*u**3
        assert K.norm(gamma)==D
        remaining=abs(D.numerator);support=[]
        for p in [x["place"] for x in raw["local"] if x["place"]!="infinity"]:
            exponent=0
            while remaining%p==0:remaining//=p;exponent+=1
            if exponent:support.append([p,exponent])
        assert remaining==1 and D.denominator==1
        CT=next(x["matrix"] for x in local["ct"] if x["u"]==u)
        words=old["CT_cross_report"]["strict_basis_in_inherited_coordinates"]
        form=[[lc.pairing(v,w,CT) for w in words] for v in words]
        packed=list(map(r.pack,form));rank=r.rank(packed)
        assert all(form[i][i]==0 for i in range(len(form)))
        assert form==list(map(list,zip(*form))) and rank%2==0
        radical=lc.orthogonal(packed,len(words))
        cases.append({"u":u,"strict_anchor_basis":old["strict_anchor_basis"],
                      "gamma_coordinates":list(map(str,gamma)),"norm_gamma":str(D),"norm_support":support,
                      "strict_CT_matrix":form,"detected_cup_image_dimension_lower_bound":rank,
                      "cup_kernel_in_retained_space_dimension_upper_bound":len(words)-rank,
                      "retained_character_annihilator_masks":radical,
                      "retained_character_annihilator_anchor_masks":[lc.lift(v,old["strict_anchor_basis"]) for v in radical],
                      "full_CT_cross_rank":old["CT_cross_report"]["cross_pairing_rank"]})
    u1=next(x for x in cases if x["u"]==-1)
    ar=r.read(art.OUTPUT)["result"]
    selected=ar["selected_character_indices"]
    dual=ar["dual_half_ideal_words"]
    projected_words=[]
    for ctrow in u1["strict_CT_matrix"]:
        coordinates=r.pack(ctrow[i] for i in selected)
        projected_words.append(lc.lift(coordinates,dual))
    projection_rank=r.rank(projected_words)
    matrix=ar["Artin_matrix_rows"]
    for i,word in enumerate(projected_words):
        assert [sum(matrix[j][k] for k in range(5) if word>>k&1)%2 for j in selected]==[
            u1["strict_CT_matrix"][i][j] for j in selected]
    return {"cases":cases,"u_minus_one_elementary_factor_projection":{
                "selected_character_indices":selected,"factor_dual_half_ideal_words":dual,
                "projected_cup_half_ideal_words":projected_words,"projected_image_dimension":projection_rank,
                "cup_image_in_complement_dimension_lower_bound":
                    u1["detected_cup_image_dimension_lower_bound"]-projection_rank,
                "coordinate_source":"Inferred through the proved cup/CT identity and retained CT values. No independent norm-witness evaluation.",
                "left_unresolved_strict_class_anchor_mask":u1["retained_character_annihilator_anchor_masks"][0]},
            "boundary":"Cup-image ranks are lower bounds using retained strict characters. Vanishing against those characters does not prove the cup ideal class zero. Explicit norm witnesses and independent recovery of the nine target CT bits remain uncomputed."}


def build(check=False):
    paths=(Path(__file__),PROTOCOL,td.OUTPUT,strict.OUTPUT,art.OUTPUT,lc.INPUT,af.INPUT)
    report={"schema":"rank-jump.strict-cup-obstruction.v1",
            "bindings":{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths},
            "finite_action_checks":finite_check(),**calculate()}
    if check:
        assert r.read(OUTPUT)==report;print("PASS 192 cup-action identities and retained obstruction ranks")
    else:r.write_new(OUTPUT,report)
    for row in report["cases"]:
        print(row["u"],"strict",len(row["strict_anchor_basis"]),
              "cup image >=",row["detected_cup_image_dimension_lower_bound"],
              "kernel <=",row["cup_kernel_in_retained_space_dimension_upper_bound"])
    print(report["u_minus_one_elementary_factor_projection"])


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mode",choices=("build","check"))
    build(p.parse_args().mode=="check")
