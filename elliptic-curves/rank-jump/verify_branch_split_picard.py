#!/usr/bin/env python3
"""Independent fibre point counts for the branch-character Picard certificate."""
import argparse
from pathlib import Path
import hashlib
import struct
import retrospective as r
import branch_split_picard as experiment


def verify(prime):
    from sage.all import GF,PolynomialRing,EllipticCurve,ZZ,QQ,pari
    from sage.version import version
    data=r.read(experiment.INPUT);assert data["bindings"]==experiment.bindings()
    raw=next(row for row in data["primes"] if row["p"]==prime)
    assert raw["status"]=="COUNTED"
    anchor=r.read(experiment.lc.INPUT)["anchor"]
    A,B=map(int,anchor["short_model_ainvariants"][3:]);p=prime
    retained={row["degree"]:row for row in raw["fields"]}
    decoded={degree:experiment.decode(row) for degree,row in retained.items()}
    Fp=GF(p);R=PolynomialRing(Fp,"x");n=retained[2]["quadratic_nonresidue"]
    Fq=GF(p*p,name="w",modulus=R([-n,0,1]));w=Fq.gen()
    roots=retained[1]["roots"];assert roots==retained[2]["roots"]
    assert sorted(map(int,R([B,A,0,1]).roots(multiplicities=False)))==roots
    independent={1:[],2:[None]*(p*p)};curves=0;singular=0

    def trace(F,u):
        nonlocal curves,singular
        rs=[F(t)+u*F(t*t) for t in roots]
        if len(set(rs))<3:
            double=next(x for x in rs if rs.count(x)==2)
            single=next(x for x in rs if rs.count(x)==1)
            singular+=1
            return -1 if (double-single).is_square() else 1
        E=EllipticCurve(F,[0,2*F(A)*u,0,F(A)+3*F(B)*u+F(A*A)*u*u,
                              F(B)+F(A*B)*u*u-F(B*B)*u**3])
        curves+=1
        return int(E.cardinality(algorithm="pari")-F.cardinality()-1)

    for a in range(p):
        v=trace(Fp,Fp(a));independent[1].append(v)
        # Independent count in degree two, including the rational base parameters.
        independent[2][a]=trace(Fq,Fq(a))
    for b in range(1,(p-1)//2+1):
        for a in range(p):
            value=trace(Fq,Fq(a)+Fq(b)*w)
            independent[2][a+p*b]=value
            independent[2][a+p*(p-b)]=value
    assert independent==decoded
    records=[]
    for degree,F in ((1,Fp),(2,Fq)):
        traces=independent[degree];values=[0,0,0]
        for index,v in enumerate(traces):
            u=F(index) if degree==1 else F(index%p)+F(index//p)*w
            gammas=[1-u*F(t) for t in roots]
            for k in range(3):
                other=[i for i in range(3) if i!=k]
                d=gammas[other[0]]*gammas[other[1]]
                chi=0 if d==0 else (1 if d.is_square() else -1)
                values[k]+=chi*v
        assert values==retained[degree]["residual_H2_traces"]
        records.append({"degree":degree,"base_parameters":len(traces),
                        "trace_sha256":hashlib.sha256(struct.pack("<"+"h"*len(traces),*traces)).hexdigest(),
                        "residual_H2_traces":values})
    report=r.read(experiment.OUTPUT)
    for row in [x for x in report["reductions"] if x["p"]==p and x["status"]=="RHO_20_REDUCTION"]:
        t=-row["quadratic_factor_ascending"][1]
        sf=int(ZZ(t*t-4*p*p).squarefree_part())
        assert sf==row["NS_discriminant_squareclass"]
        # Over F_(p^2), the residual quadratic is X^2-(t^2-2p^2)X+p^4.
        normalized=QQ(4)-QQ(t*t)/(p*p)
        assert normalized>0 and not t%p==0
        assert int(ZZ((-normalized/(p*p)).numerator()*(-normalized/(p*p)).denominator()).squarefree_part())==sf
    out={"prime":p,"status":"PASS","software":{"sage":version,"pari":str(pari.version())},
         "smooth_fibre_cardinality_calls":curves,"singular_fibre_checks":singular,"fields":records,
         "input_sha256":r.digest(experiment.INPUT.read_bytes()),
         "analysis_sha256":r.digest(experiment.OUTPUT.read_bytes()),
         "verifier_sha256":r.digest(Path(__file__).read_bytes())}
    print(out)
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--prime",type=int,required=True)
    p.add_argument("--destination",type=Path)
    args=p.parse_args();result=verify(args.prime)
    if args.destination:r.write_new(args.destination,result)
