#!/usr/bin/env sage -python
"""Exact marked-fibration distinction using intrinsic j-map critical values.

A fixed finite witness; no Mobius search, basis enumeration or fibre search.
The shared F197 point count is only a comparison, never a surface identity.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import resource
from sage.all import GF, PolynomialRing, QQ, matrix

ROOT=Path(__file__).resolve().parents[2]
SOURCES={
 "hill":"elkies-k3/data/external/hill2026/data/parent-rank17.json",
 "record":"elkies-k3/data/fibrations/icarm4239_x1152_generic_parent.json",
}


def run():
    resource.setrlimit(resource.RLIMIT_CPU,(10,12))
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    prime=197
    R=PolynomialRing(GF(prime),"t");t=R.gen()
    rows=[]
    chi=[0]+[1 if pow(x,(prime-1)//2,prime)==1 else -1 for x in range(1,prime)]
    for label,source in SOURCES.items():
        path=ROOT/source;data=json.loads(path.read_text())
        if label=="hill":
            aa,bb=[data["jacobian_model"][k] for k in ["a4_low_to_high","a6_low_to_high"]]
        else:
            aa,bb=data["short_A_low_to_high"],data["short_B_low_to_high"]
        A,B=R(list(map(QQ,aa))),R(list(map(QQ,bb)))
        C=3*A.derivative()*B-2*A*B.derivative()
        Delta=4*A**3+27*B**2
        assert [A.degree(),B.degree(),C.degree(),Delta.degree()]==[8,12,18,24]
        assert all(f.is_squarefree() for f in [A,B,C,Delta])
        assert A.gcd(B).degree()==0 and C.gcd(A*B*Delta).degree()==0
        c=C.monic()
        value=(6912*A**3*Delta.inverse_mod(c))%c
        multiplication=matrix(GF(prime),18,18,lambda i,j:((value*t**j)%c)[i])
        trace=int(multiplication.trace())
        characteristic=multiplication.charpoly()
        assert (-characteristic[17])%prime==trace
        count=0
        for parameter in list(range(prime))+[None]:
            a,b=(int(A(parameter)),int(B(parameter))) if parameter is not None else (int(A[8]),int(B[12]))
            count+=prime+1+sum(chi[(x**3+a*x+b)%prime] for x in range(prime))
        rows.append({"label":label,"source":source,"source_sha256":sha256(path.read_bytes()).hexdigest(),
            "critical_polynomial_monic":[int(x) for x in c.list()],
            "j_in_critical_algebra":[int(x) for x in value.list()],
            "critical_value_trace":trace,
            "critical_value_characteristic":[int(x) for x in characteristic.list()],
            "surface_Fp_points":count,"surface_H2_trace":count-1-prime**2})
    assert rows[0]["critical_value_trace"]!=rows[1]["critical_value_trace"]
    return {"schema":"x1152-critical-value-distinction-v1","prime":prime,"status":"PASS",
        "models":rows,
        "conclusion":"No PGL2 base change over Qbar identifies the two j-maps; consequently the marked elliptic fibrations are not isomorphic.",
        "surface_identity":"UNKNOWN",
        "boundary":"A fibration distinction is not a K3 surface distinction. Equal F197 surface counts supply no surface identity, Picard bound or basis transport."}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",required=True,type=Path)
    a=p.parse_args();result=run();a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open("x") as f:json.dump(result,f,indent=2,sort_keys=True);f.write("\n")
    print(json.dumps({"status":"PASS","critical_value_traces":[r["critical_value_trace"] for r in result["models"]],
                      "surface_points":[r["surface_Fp_points"] for r in result["models"]]}))
