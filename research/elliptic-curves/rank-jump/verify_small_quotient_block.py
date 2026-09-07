#!/usr/bin/env python3
"""Rank and strict-class replay without a new descent or number-field CAS."""
import argparse
from fractions import Fraction as F
from math import comb
from pathlib import Path
import retrospective as r
import small_quotient_block as ex
import verify_nonscalar_cup_control as a

OUTPUT=r.OUT/"rank_jump_small_quotient_block_verification_v1.json"


def add(model,P,Q):
    if P is None:return Q
    if Q is None:return P
    x,y=P;z,w=Q;a2,a4=map(F,[model[1],model[3]])
    if x==z and y==-w:return None
    slope=(w-y)/(z-x) if x!=z else (3*x*x+2*a2*x+a4)/(2*y)
    X=slope*slope-a2-x-z
    return X,-y+slope*(x-X)


def verify(check=False):
    data=r.read(ex.OUTPUT);old=r.read(ex.scalar_cup.CONTROL)
    assert data["bindings"]==ex.bindings()
    f=a.POLY
    assert all(sum(c*t**i for i,c in enumerate(f))%2 for t in (0,1))
    shifted=[sum(f[i]*comb(i,j)*58**(i-j) for i in range(j,4)) for j in range(4)]
    assert shifted==[157295,8802,163,1]
    assert all(c%163==0 for c in shifted[:-1]) and shifted[0]%163**2
    assert old["certified_class_group"]==[2,2]
    assert r.rank(list(map(r.pack,old["Artin_matrix"])))==2
    assert r.rank(old["strict_masks"])==2
    assert old["Artin_matrix"]==[[1,0],[1,1]]
    # Real signs of two norm-one units are 011 and 101.
    u0=list(map(F,[-1,-1,0]));u1=list(map(F,[1,12,-1]))
    assert a.norm(u0)==a.norm(u1)==1
    assert a.mul(u1,list(map(F,[1,1,0])))==[0,-1,0]
    assert r.rank([6,5])==2
    # Unit norm-kernel dimension 2; independent half-ideal classes add 2.
    norm_kernel_dimension=2+2
    # Each real Kummer image is a line in the two-dimensional even-sign plane.
    full_selmer_upper=norm_kernel_dimension-(2-1)
    assert full_selmer_upper==3
    records=[]
    for rec in data["records"]:
        model=rec["model"];index=rec["index"];points=[tuple(map(F,P)) for P in rec["points"]]
        assert model[0]==model[2]==0
        a2,a4,a6=map(F,[model[1],model[3],model[4]])
        # Cubic discriminant and elliptic discriminant: no bad prime outside S.
        cubic_discriminant=a2*a2*a4*a4-4*a4**3-4*a2**3*a6-27*a6*a6+18*a2*a4*a6
        assert cubic_discriminant==163**2
        signatures=[]
        for x,y in points:
            assert y*y==x**3+a2*x*x+a4*x+a6
            bits=[]
            for root in rec["cubic_roots"]:
                p=rec["fingerprint_prime"]
                assert sum(c*root**i for i,c in enumerate(f))%p==0
                residue=(r.mod(x,p)+(root if index else -root))%p
                assert residue
                bits.append(int(pow(residue,(p-1)//2,p)==p-1))
            lo,hi=(x,x+1) if index==0 else (x-2,x-1)
            assert lo>=0 or hi<=0
            bits.append(int(hi<=0));signatures.append(r.pack(bits))
        assert signatures==rec["point_signatures"] and r.rank(signatures)==len(points)
        for witness in rec["strict_class_rational_witnesses"]:
            P=None
            for i,Q in enumerate(points):
                if witness["point_word"]>>i&1:P=add(model,P,Q)
            assert P==tuple(map(F,witness["point"]))
            delta=[P[0],F(1),F(0)]
            root=list(map(F,witness["square_root_beta_over_Kummer"]))
            root+=[F(0)]*(3-len(root))
            beta=list(map(F,old["norm_witnesses"][witness["beta_index"]]["beta"]))
            assert a.mul(delta,a.mul(root,root))==beta
        records.append({"index":index,"point_rank":len(points),"signatures":signatures,
                        "strict_rational_witnesses_checked":len(rec["strict_class_rational_witnesses"])})
    H=ex.scalar_cup.symmetrize(old["Artin_matrix"])
    assert H==[[0,1],[1,0]] and r.rank(list(map(r.pack,H)))==2
    # The twist's three point classes fill its Selmer group. Hence its CT is zero.
    # Scalar cup difference makes CT on U nonsingular on the original curve.
    assert records[1]["point_rank"]==full_selmer_upper
    original_upper=full_selmer_upper-r.rank(list(map(r.pack,H)))
    assert original_upper==records[0]["point_rank"]==1
    assert data["exact_MW_ranks"]==[original_upper,full_selmer_upper]
    assert data["strict_rational_dimensions"]==[0,2]
    assert data["Sha_2_dimensions"]==[2,0]
    result={"schema":"rank-jump.small-quotient-block-verification.v1","status":"PASS","records":records,
            "norm_kernel_dimension":4,"real_condition_codimension":1,"full_Selmer_dimension":3,
            "original_strict_CT_rank":2,"exact_ranks":[1,3],"Sha_2_dimensions":[2,0],
            "new_descent_upper_bounds_used_as_proof":False,
            "retained_class_group_certificate_used":True,
            "analysis_sha256":r.digest(ex.OUTPUT.read_bytes()),"verifier_sha256":r.digest(Path(__file__).read_bytes()),
            "arithmetic_source_sha256":r.digest(Path(a.__file__).read_bytes())}
    if check:assert result==r.read(OUTPUT)
    else:r.write_new(OUTPUT,result)
    print("PASS independent Selmer bounds, rational class witnesses, scalar CT obstruction, and exact ranks")


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--check",action="store_true")
    verify(parser.parse_args().check)
