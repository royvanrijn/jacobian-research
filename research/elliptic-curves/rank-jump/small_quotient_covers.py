#!/usr/bin/env python3
"""Explicit paired two-cover quadrics for the two strict classes."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import retrospective as r
import small_quotient_block as ex
import verify_small_quotient_block as proof
import verify_nonscalar_cup_control as a
from nonscalar_unramified_lifts import inverse

OUTPUT=r.OUT/"rank_jump_small_quotient_covers_v1.json"
PAIRS=[(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]


def evaluate(form,z):
    return sum(c*z[i]*z[j] for c,(i,j) in zip(form,PAIRS))


def build(check=False):
    data=r.read(ex.OUTPUT);verified=r.read(proof.OUTPUT);source=r.read(ex.scalar_cup.CONTROL)
    assert verified["status"]=="PASS" and verified["analysis_sha256"]==r.digest(ex.OUTPUT.read_bytes())
    records=[]
    for witness in data["records"][1]["strict_class_rational_witnesses"]:
        index=witness["beta_index"];beta=list(map(F,source["norm_witnesses"][index]["beta"]))
        columns=[]
        for i,j in PAIRS:
            ei=[F(k==i) for k in range(3)];ej=[F(k==j) for k in range(3)]
            column=a.mul(beta,a.mul(ei,ej))
            columns.append([v*(1 if i==j else 2) for v in column])
        forms=[[column[i] for column in columns] for i in range(3)]
        z=inverse(list(map(F,witness["square_root_beta_over_Kummer"])))
        X,Y=map(F,witness["point"])
        if 25*a.norm(z)!=Y:z=[-v for v in z]
        assert [evaluate(form,z) for form in forms]==[X,1,0]
        assert 25*a.norm(z)==Y
        records.append({"beta_index":index,"beta":list(map(str,beta)),"monomials":["u^2","u*v","u*w","v^2","v*w","w^2"],
                        "Q0":list(map(str,forms[0])),"Q1":list(map(str,forms[1])),"Q2":list(map(str,forms[2])),
                        "original_cover":["Q2=0","Q1+h^2=0"],"twist_cover":["Q2=0","Q1-h^2=0"],
                        "map":"x=Q0/h^2, y=25*N(u+v*theta+w*theta^2)/h^3",
                        "twist_rational_point":list(map(str,z))+["1"],"elliptic_image":witness["point"],
                        "original_rational_solubility":"NO: nonzero Sha[2] class, exact rank/CT certificate",
                        "twist_rational_solubility":"YES: explicit projective point",
                        "local_solubility":"YES on both curves: retained strict Selmer classes"})
    out={"schema":"rank-jump.small-quotient-covers.v1","status":"PASS","records":records,
         "block_sha256":r.digest(ex.OUTPUT.read_bytes()),"rank_verification_sha256":r.digest(proof.OUTPUT.read_bytes()),
         "producer_sha256":r.digest(Path(__file__).read_bytes()),
         "boundary":"Solubility of these four fixed two-covers; not a general rule that a zero cup product implies rational solubility."}
    if check:assert out==r.read(OUTPUT)
    else:r.write_new(OUTPUT,out)
    print("PASS paired explicit two-covers: original locally soluble but pointless; twist explicitly soluble")
    for rec in records:print(rec["beta_index"],rec["Q1"],rec["Q2"],rec["twist_rational_point"])


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--check",action="store_true")
    build(parser.parse_args().check)
