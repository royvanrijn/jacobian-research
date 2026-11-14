#!/usr/bin/env sage -python
"""Construct explicit 2-isogeny/Jacobian endpoints for the frozen15 paired bases.

Exact arithmetic only, no point search, rank computation or Selmer assertion.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import resource
from sage.all import EllipticCurve, PolynomialRing, QQ

ROOT=Path(__file__).resolve().parents[2]
PILOT=ROOT/"artifacts/generated-results/hill1152-transfer-20261113-v2"
DATA=ROOT/"elkies-k3/data/external/hill2026/data"


def digest(path):return sha256(path.read_bytes()).hexdigest()


def invariants(f):
    e,d,c,b,a=[f[i] for i in range(5)]
    return 12*a*e-3*b*d+c*c,72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c*c*c


def rows(f):return [str(f[i]) for i in range(max(0,int(f.degree()))+1)]


def run():
    resource.setrlimit(resource.RLIMIT_CPU,(30,35))
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    result=json.loads((PILOT/"result.json").read_text())
    cert=json.loads((PILOT/"four-generic-covers.json").read_text())
    c=QQ(result["common_chart"]["c"])
    R=PolynomialRing(QQ,"U");u=R.gen();K=R.fraction_field()
    covers=result["quadratic_covers"]
    maps={}
    for f in cert["families"]:
        maps[f["index"]]=tuple(R(f[k]) for k in ["N","D","W"])
    for index,published in [(13,4339),(14,36773)]:
        d=json.loads((DATA/f"cover-{published}.json").read_text())["quadratic_basechange"]
        oldN,oldD,oldW=[R(d[k]) for k in ["parameter_numerator","parameter_denominator","branch_root_numerator"]]
        maps[index]=(oldD,oldN-c*oldD,oldW)
    records=[]
    pairs=result["locally_surviving_disjoint_pairs"]
    assert len(pairs)==15
    for i,j in pairs:
        q,r=R(covers[i]["q"]),R(covers[j]["q"])
        assert q.degree()==r.degree()==2 and q.is_squarefree() and r.is_squarefree() and q.gcd(r).degree()==0
        a=q[2]*r[0]+q[0]*r[2]-q[1]*r[1]/2
        discq=q[1]**2-4*q[2]*q[0];discr=r[1]**2-4*r[2]*r[0]
        b=discq*discr/16
        resultant=q.resultant(r)
        assert a*a-4*b==resultant
        assert b and resultant
        E=EllipticCurve(QQ,[0,a,0,b,0])
        N,D,W=maps[i]
        assert N.gcd(D).degree()==0
        assert q[2]*N*N+q[1]*N*D+q[0]*D*D==W*W
        F=r[2]*N*N+r[1]*N*D+r[0]*D*D
        assert F.degree() in (3,4) and F.is_squarefree()
        I,J=invariants(F)
        A,B=b-a*a/3,2*a**3/27-a*b/3
        assert A and B
        k2=((-27*J)/B)/((-27*I)/A)
        assert k2.is_square()
        k=k2.sqrt()
        assert k and k**4*A==-27*I and k**6*B==-27*J
        # Hence (X,Y)=(k^2*(x+a/3),k^3*y) is an exact Q-isomorphism
        # E -> Jac(Y^2=F), with no rational point on the torsor required.
        product=q*r
        Ip,Jp=invariants(product)
        assert Ip==a*a+12*b and Jp==2*a*(36*b-a*a)
        quotient=EllipticCurve(QQ,[0,-2*a,0,resultant,0])
        x=u;X=x+a+b/x
        assert (x**3+a*x*x+b*x)*(1-b/x**2)**2==X**3-2*a*X*X+resultant*X
        records.append({
            "pair":[i,j],"jacobian_a_invariants":["0",str(a),"0",str(b),"0"],
            "quadratic_discriminants":[str(discq),str(discr)],
            "radicand_resultant":str(resultant),
            "conic_map":{"N":rows(N),"D":rows(D),"W":rows(W)},
            "paired_base_quartic":rows(F),"quartic_I":str(I),"quartic_J":str(J),
            "jacobian_to_quartic_scale":str(k),
            "quotient_a_invariants":["0",str(-2*a),"0",str(resultant),"0"],
            "product_quartic_I":str(Ip),"product_quartic_J":str(Jp),
            "phi_kernel":"O,(0,0)",
            "phi_formula":"X=x+a+b/x; Y=y*(1-b/x^2)",
            "quotient_to_product_jacobian":"X=9*x-6*a; Y=27*y",
            "rational_points_on_paired_base":"UNKNOWN",
            "complete_selmer_bound":"UNKNOWN"})
    # A twist-sensitive comparison against the published paired-base Jacobian.
    published=json.loads((DATA/"pair.json").read_text())
    E0=EllipticCurve(QQ,list(map(QQ,published["auxiliary_a_invariants"])))
    last=next(r for r in records if r["pair"]==[13,14])
    E=EllipticCurve(QQ,list(map(QQ,last["jacobian_a_invariants"])))
    assert E.is_isomorphic(E0)
    return {"schema":"hill1152-paired-jacobians-2isogenies-v1","status":"PASS",
        "pilot_result_sha256":digest(PILOT/"result.json"),
        "cover_certificate_sha256":digest(PILOT/"four-generic-covers.json"),
        "published_control_source_sha256":digest(DATA/"pair.json"),
        "source_hashes":{f"cover-{i}.json":digest(DATA/f"cover-{i}.json") for i in [4339,36773]},
        "pairs":records,"published_pair_Q_isomorphism_checked":True,
        "boundary":"Classical Jacobian/isogeny construction instantiated on15 frozen bases. No torsor point, full rank, complete Selmer space or global-solubility conclusion."}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",required=True,type=Path)
    a=p.parse_args();result=run();a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open("x") as f:json.dump(result,f,indent=2,sort_keys=True);f.write("\n")
    print(json.dumps({"status":"PASS","paired_jacobians":len(result["pairs"]),"degree2_isogenies":len(result["pairs"])}))
