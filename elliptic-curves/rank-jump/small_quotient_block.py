#!/usr/bin/env python3
"""Exact point/Kummer witnesses for the fixed small strict block."""
import argparse
from pathlib import Path
import retrospective as r
import small_quotient_solubility as desc
import nonscalar_cup_control as cup
import nonscalar_cup_orbit as orbit
import scalar_cup

OUTPUT=r.OUT/"rank_jump_small_quotient_block_v1.json"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),desc.OUTPUT,cup.OUTPUT,orbit.OUTPUT,scalar_cup.CONTROL)}


def build(check=False):
    from sage.all import QQ,PolynomialRing,NumberField,EllipticCurve
    old=r.read(scalar_cup.CONTROL);raw=r.read(desc.OUTPUT)
    R=PolynomialRing(QQ,"t");f=R(old["polynomial_ascending"]);K=NumberField(f,"theta");theta=K.gen()
    assert f.is_irreducible() and f.discriminant()==163**2
    p=37;roots=[i for i in range(p) if int(f(i))%p==0];assert len(roots)==3
    beta=[K(list(map(QQ,w["beta"]))) for w in old["norm_witnesses"]]
    records=[]
    for index,rawrow in enumerate(raw["records"]):
        assert rawrow["status"]=="COMPUTED"
        E=EllipticCurve(rawrow["model"]);points=[E(list(map(QQ,P))) for P in rawrow["points"]]
        assert len(points)<=r.read(desc.PROTOCOL)["limits"]["maximum_point_word_dimension"]
        signatures=[]
        for P in points:
            x=P[0];chars=[]
            for root in roots:
                residue=(r.mod(str(x),p)+(root if index else -root))%p
                assert residue
                chars.append(int(pow(residue,(p-1)//2,p)==p-1))
            # Use the middle embedding for E0, the first for its twist.
            # Their theta intervals are (-1,0) and (-2,-1), respectively.
            lo,hi=((x,x+1) if index==0 else (x-2,x-1))
            assert lo>=0 or hi<=0
            chars.append(int(hi<=0))
            signatures.append(r.pack(chars))
        assert r.rank(signatures)==len(points)
        matches=[]
        if index==1:
            for j,b in enumerate(beta):
                found=None
                for word in range(1,1<<len(points)):
                    P=E(0)
                    for i,Q in enumerate(points):
                        if word>>i&1:P+=Q
                    delta=K(P[0])+theta;ratio=b/delta
                    if ratio.is_square():
                        root=ratio.sqrt();assert delta*root*root==b
                        found={"beta_index":j,"point_word":word,"point":list(map(str,P[:2])),
                               "square_root_beta_over_Kummer":list(map(str,root.list()))}
                        break
                assert found is not None
                matches.append(found)
        records.append({"index":index,"model":rawrow["model"],"points":rawrow["points"],
                        "fingerprint_prime":p,"cubic_roots":roots,"point_signatures":signatures,
                        "certified_point_rank":len(points),"strict_class_rational_witnesses":matches})
    result={"schema":"rank-jump.small-quotient-block.v1","bindings":bindings(),"records":records,
            "full_2_Selmer_dimensions":[3,3],"exact_MW_ranks":[1,3],"Sha_2_dimensions":[2,0],
            "strict_rational_dimensions":[0,2],"strict_CT_matrices":[[[0,1],[1,0]],[[0,0],[0,0]]],
            "proof_dependencies":"Real-sign upper bound on the four-dimensional norm-kernel class space; independent point fingerprints; scalar cup difference; see canonical proof.",
            "scope":"Two fixed small control curves. No production specialization rank claim."}
    if check:assert result==r.read(OUTPUT)
    else:r.write_new(OUTPUT,result)
    print("PASS small block: exact ranks 1/3, strict rational dimensions 0/2, Sha[2] dimensions 2/0")
    print(records[1]["strict_class_rational_witnesses"])


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--check",action="store_true")
    build(parser.parse_args().check)
