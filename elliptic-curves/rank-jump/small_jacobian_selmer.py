#!/usr/bin/env python3
"""Labelled small Jacobian Selmer calculation, with finite local certificates."""
import argparse
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import retrospective as r
import verify_nonscalar_cup_control as arithmetic
import nonscalar_cup_orbit as orbit
import small_quotient_block as elliptic
import verify_small_quotient_block as verified

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/"SMALL_JACOBIAN_SELMER_PROTOCOL.json"
OUTPUT=r.OUT/"rank_jump_small_jacobian_selmer_v1.json"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,orbit.OUTPUT,elliptic.OUTPUT,verified.OUTPUT,
             HERE/"JACOBIAN_LOCAL_CONDITIONS_AND_CT.md",HERE/"RELATIVE_FULL_SELMER_THEOREM.md")}


def canonical(word):return min(word,word^63)


def act(word,pair_flip=0,rotate=False):
    out=0
    for i in range(6):
        if word>>i&1:
            pair,side=divmod(i,2)
            side^=(pair_flip>>pair)&1
            if rotate:pair=(pair+1)%3
            out^=1<<(2*pair+side)
    return canonical(out)


def intersection(A,B):
    span=lambda base:{r.pack([])} if not base else {
        __import__('functools').reduce(int.__xor__,(b for i,b in enumerate(base) if mask>>i&1),0)
        for mask in range(1<<len(base))}
    return r.basis(sorted(span(A)&span(B)))


def build(check=False):
    old=r.read(elliptic.OUTPUT);proof=r.read(verified.OUTPUT)
    assert proof["status"]=="PASS" and proof["analysis_sha256"]==r.digest(elliptic.OUTPUT.read_bytes())
    tau=list(map(F,r.read(orbit.OUTPUT)["automorphism_theta"]))
    u0=list(map(F,[-1,-1,0]));u1=list(map(F,[1,12,-1]));theta=list(map(F,[0,1,0]))
    assert arithmetic.mul(tau,list(map(F,[1,1,0])))==[-1,0,0]
    assert orbit.substitute(u0,tau)==u1
    assert arithmetic.mul(u0,tau)==[1,0,0]
    assert arithmetic.mul(u0,u1)==theta
    squares={tuple(int(v)%8 for v in arithmetic.mul(z,z)) for z in product(range(8),repeat=3)}
    square_rows=[{"unit":list(map(str,u)),"residue_mod8":list(map(lambda x:int(x)%8,u)),
                  "has_square_root_mod8":tuple(int(x)%8 for x in u) in squares} for u in (u0,u1,theta)]
    assert not any(x["has_square_root_mod8"] for x in square_rows)
    words=sorted({canonical(w) for w in range(64) if w.bit_count()%2==0});assert len(words)==16
    # Global coordinates beta0,beta1,u0,u1; strict classes have zero local images.
    S0=[1,2,4];Splus=[1,2,4];Sminus=[1,2,8];rows=[]
    for sign,S1 in [(1,Splus),(-1,Sminus)]:
        SC=intersection(S0,S1);SD=r.basis(S0+S1)
        expected_common=3 if sign==1 else 2
        assert len(SC)==expected_common and len(SD)==(3 if sign==1 else 4)
        global_flips=(1,2,4) if sign==1 else (3,5)
        global_fixed=[w for w in words if act(w,rotate=True)==w and all(act(w,pair_flip=f)==w for f in global_flips)]
        assert global_fixed==[0]
        real_flip=1 if sign==1 else 6
        real_fixed=[w for w in words if act(w,pair_flip=real_flip)==w]
        real_torsion_dim=(len(real_fixed)).bit_length()-1
        real_middle_dim=real_torsion_dim-2
        local=[]
        for place in ("2","163","infinity"):
            ell_dim=0 if place=="163" else 1
            common=ell_dim if sign==1 else 0
            summed=2*ell_dim-common
            middle=2 if place=="2" else (real_middle_dim if place=="infinity" else 0)
            connecting=summed+common-middle
            assert connecting>=0
            if place!="infinity":assert connecting==0
            local.append({"place":place,"elliptic_local_dimension":ell_dim,"common_dimension":common,
                          "sum_dimension":summed,"connecting_rank":connecting,"Jacobian_local_dimension":middle})
        difference=[[0,1,0],[1,0,0],[0,0,0]] if sign==1 else [[0,0],[0,0]]
        ct_rank=r.rank(list(map(r.pack,difference)));radical=expected_common-ct_rank
        selmer_dim=len(SD)+radical
        mw_rank=sum(old["exact_MW_ranks"]) if sign==1 else 2*old["exact_MW_ranks"][0]
        rational_image=len(intersection([4],Splus if sign==1 else [8]))
        assert rational_image==(1 if sign==1 else 0)
        rows.append({"sign":sign,"second_quotient":"minus-one twist" if sign==1 else "original curve",
                     "pullback_of_standard_second_Selmer_basis":S1,"common_Selmer_basis":SC,"sum_Selmer_basis":SD,
                     "global_Jacobian_2_torsion_classes":global_fixed,"real_2_torsion_dimension":real_torsion_dim,
                     "local_conditions":local,"difference_CT_matrix":difference,"difference_radical_dimension":radical,
                     "Jacobian_2_Selmer_dimension":selmer_dim,"Jacobian_exact_rank":mw_rank,
                     "Jacobian_Sha_2_dimension":selmer_dim-mw_rank,
                     "rational_image_dimension_in_common_Selmer":rational_image,
                     "strict_image_dimension_from_Jacobian_Selmer":0 if sign==1 else 2,
                     "strict_image_dimension_from_rational_Jacobian_points":0,
                     "index_of_product_rational_image":2**rational_image})
    out={"schema":"rank-jump.small-jacobian-selmer.v1","bindings":bindings(),"status":"PASS",
         "dyadic_residue_count":512,"dyadic_unit_checks":square_rows,"torsion_pullback":"tau(theta)=-1/(theta+1)",
         "rows":rows,"positive_sign_extra_divisor":{
             "divisor":"infinity_plus - (0,1)","projection_to_E0":["-1","-1"],
             "projection_to_Eplus":["0","1"],"common_Kummer_class":"[-1-theta]",
             "identity":"2j=Phi(-P,R), P=(-1,1), R=(0,1)"},
         "boundary":"Exact for the two fixed genus-two controls. Jacobian Selmer lifts need not be rational; higher Sha and production ranks not inferred."}
    if check:assert out==r.read(OUTPUT)
    else:r.write_new(OUTPUT,out)
    print("PASS labelled Jacobian Selmer/rank/Sha2 dimensions:",[(v["sign"],v["Jacobian_2_Selmer_dimension"],v["Jacobian_exact_rank"],v["Jacobian_Sha_2_dimension"]) for v in rows])


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--check",action="store_true")
    build(parser.parse_args().check)
