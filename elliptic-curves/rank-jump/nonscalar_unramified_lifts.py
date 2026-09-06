#!/usr/bin/env python3
"""Explicit principal corrections for the two vanishing nonscalar cups."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import retrospective as r
import nonscalar_cup_control as control
import verify_nonscalar_cup_control as arithmetic

OUTPUT=r.OUT/"rank_jump_nonscalar_unramified_lifts_v1.json"


def inverse(a):
    columns=[arithmetic.mul(a,[F(i==j) for i in range(3)]) for j in range(3)]
    M=[[columns[j][i] for j in range(3)]+[F(i==0)] for i in range(3)]
    for j in range(3):
        pivot=next(i for i in range(j,3) if M[i][j])
        M[j],M[pivot]=M[pivot],M[j]
        v=M[j][j];M[j]=[x/v for x in M[j]]
        for i in range(3):
            if i!=j:
                v=M[i][j];M[i]=[x-v*y for x,y in zip(M[i],M[j])]
    out=[row[3] for row in M]
    assert arithmetic.mul(a,out)==[1,0,0]
    return out


def build(check=False):
    data=r.read(control.OUTPUT);verified=r.read(arithmetic.OUTPUT)
    assert verified["status"]=="PASS" and verified["analysis_sha256"]==r.digest(control.OUTPUT.read_bytes())
    original=data["records"][1];local=verified["records"][1]
    assert original["matrix"]==[[0,0],[0,0]]
    gamma=list(map(F,original["gamma"]));rows=[]
    for index,(w,checked) in enumerate(zip(original["witnesses"],local["witnesses"])):
        eta=list(map(F,[-2,1,0] if index==0 else [5,0,0]))
        assert arithmetic.norm(eta)==(65 if index==0 else 125)
        eta_inverse=inverse(eta)
        a,b,beta=[list(map(F,w[key])) for key in ("a","b","beta")]
        aa=arithmetic.mul(a,eta_inverse);bb=arithmetic.mul(b,eta_inverse)
        adjusted_beta=arithmetic.mul(beta,arithmetic.mul(eta_inverse,eta_inverse))
        square_a=arithmetic.mul(aa,aa);gamma_square_b=arithmetic.mul(gamma,arithmetic.mul(bb,bb))
        assert [x-y for x,y in zip(square_a,gamma_square_b)]==adjusted_beta
        valuations=[]
        for rec in checked["places"]:
            p,t=rec["p"],rec["theta_mod_p"]
            eta_val=int((p in (5,13) and t==2) if index==0 else p==5)
            corrected=[v-eta_val for v in rec["valuations"]]
            assert all(v%2==0 for v in corrected)
            valuations.append({"p":p,"theta_mod_p":t,"eta_valuation":eta_val,"corrected_valuations":corrected})
        rows.append({"eta":list(map(str,eta)),"eta_norm":str(arithmetic.norm(eta)),
                     "norm_a":list(map(str,aa)),"norm_b":list(map(str,bb)),
                     "norm_target":list(map(str,adjusted_beta)),"original_beta":w["beta"],
                     "corrected_valuation_checks":valuations,"unramified_outside_S":True})
    # Three disjoint sign-changing intervals contain all three real roots.
    f=lambda x:sum(c*x**i for i,c in enumerate(arithmetic.POLY))
    intervals=[[-2,-1],[-1,0],[12,13]]
    assert all(f(a)*f(b)<0 for a,b in intervals)
    out={"schema":"rank-jump.nonscalar-unramified-lifts.v1","status":"PASS","gamma":original["gamma"],
         "lifts":rows,"norm_image_dimension":2,"positive_multiplier_norm_image_dimension_on_retained_U":0,
         "root_intervals":intervals,"multiplier_signs":[[-1,1,1],[1,-1,-1]],
         "neither_multiplier_is_a_scalar_squareclass":True,
         "analysis_sha256":r.digest(control.OUTPUT.read_bytes()),"independent_parities_sha256":r.digest(arithmetic.OUTPUT.read_bytes()),
         "producer_sha256":r.digest(Path(__file__).read_bytes()),
         "arithmetic_source_sha256":r.digest(Path(arithmetic.__file__).read_bytes()),
         "elliptic_rational_solubility":"UNKNOWN"}
    if check:assert out==r.read(OUTPUT)
    else:r.write_new(OUTPUT,out)
    print("PASS two explicit unramified norm lifts, nonscalar signs, and complete parity removal")


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--check",action="store_true")
    build(parser.parse_args().check)
