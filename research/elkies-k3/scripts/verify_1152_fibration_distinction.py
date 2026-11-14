#!/usr/bin/env python3
"""Independent finite-polynomial replay of the intrinsic j-critical-value test.

No Sage, Mobius fitting, point search or cached arithmetic. The written
characteristic-zero invariant argument is in the canonical note.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from math import isqrt
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
SOURCES={
 "hill":"elkies-k3/data/external/hill2026/data/parent-rank17.json",
 "record":"elkies-k3/data/fibrations/icarm4239_x1152_generic_parent.json",
}
DEFAULT=ROOT/"artifacts/generated-results/x1152-fibration-distinction-20261114-v1.json"


def require(ok,message):
    if not ok:raise ValueError(message)


def trim(f):
    f=list(f)
    while f and not f[-1]:f.pop()
    return f


def add(a,b,p):
    return trim([((a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0))%p
                 for i in range(max(len(a),len(b)))])


def scale(a,s,p):return trim([s*x%p for x in a])


def mul(a,b,p):
    if not a or not b:return []
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%p
    return trim(c)


def divmod_poly(a,b,p):
    require(b,"zero divisor")
    r=trim(a);q=[0]*max(0,len(r)-len(b)+1);inv=pow(b[-1],-1,p)
    while len(r)>=len(b):
        i=len(r)-len(b);s=r[-1]*inv%p;q[i]=s
        for j,x in enumerate(b):r[i+j]=(r[i+j]-s*x)%p
        r=trim(r)
    return trim(q),r


def rem(a,b,p):return divmod_poly(a,b,p)[1]


def gcd_poly(a,b,p):
    while b:a,b=b,rem(a,b,p)
    return scale(a,pow(a[-1],-1,p),p) if a else []


def inverse(a,modulus,p):
    r0,r1=modulus,rem(a,modulus,p);s0,s1=[],[1]
    while r1:
        q,r=divmod_poly(r0,r1,p)
        r0,r1=r1,r;s0,s1=s1,add(s0,scale(mul(q,s1,p),-1,p),p)
    require(len(r0)==1,"noninvertible critical-algebra denominator")
    answer=rem(scale(s0,pow(r0[0],-1,p),p),modulus,p)
    require(rem(mul(a,answer,p),modulus,p)==[1],"inverse witness")
    return answer


def derivative(a,p):return trim([i*x%p for i,x in enumerate(a)][1:])


def reduction(cs,p):
    values=[Q(c) for c in cs]
    require(all(x.denominator%p for x in values),"nonintegral source at witness prime")
    out=trim([x.numerator*pow(x.denominator,-1,p)%p for x in values])
    require(len(out)==len(values),"source degree drop")
    return out


def evaluate(a,t,p):
    s=0
    for c in reversed(a):s=(s*t+c)%p
    return s


def critical_trace(A,B,p):
    A2=mul(A,A,p);A3=mul(A2,A,p);B2=mul(B,B,p)
    Delta=add(scale(A3,4,p),scale(B2,27,p),p)
    C=add(scale(mul(derivative(A,p),B,p),3,p),
          scale(mul(A,derivative(B,p),p),-2,p),p)
    require([len(A)-1,len(B)-1,len(C)-1,len(Delta)-1]==[8,12,18,24],"degree gates")
    for f in [A,B,C,Delta]:
        require(gcd_poly(f,derivative(f,p),p)==[1],"multiple roots")
    require(gcd_poly(A,B,p)==[1] and gcd_poly(C,mul(mul(A,B,p),Delta,p),p)==[1],"critical strata overlap")
    C=scale(C,pow(C[-1],-1,p),p)
    value=rem(mul(scale(A3,6912,p),inverse(Delta,C,p),p),C,p)
    trace=0
    for i in range(18):
        col=rem([0]*i+value,C,p)
        trace=(trace+(col[i] if i<len(col) else 0))%p
    return C,value,trace


def surface_count(A,B,p):
    chi=[0]*p
    for x in range(1,p):chi[x]=1 if pow(x,(p-1)//2,p)==1 else -1
    total=0
    for t in range(p+1):
        a,b=(evaluate(A,t,p),evaluate(B,t,p)) if t<p else (A[8],B[12])
        total+=p+1+sum(chi[(x*x*x+a*x+b)%p] for x in range(p))
    return total


def verify(result,check_sources=True):
    require(result["schema"]=="x1152-critical-value-distinction-v1","schema")
    p=result["prime"]
    require(p>24 and all(p%d for d in range(2,isqrt(p)+1)),"prime")
    require([r["label"] for r in result["models"]]==["hill","record"],"source coverage")
    traces=[];counts=[]
    for row in result["models"]:
        require(row["source"]==SOURCES[row["label"]],"source identity")
        path=ROOT/row["source"]
        require(sha256(path.read_bytes()).hexdigest()==row["source_sha256"],"source hash")
        data=json.loads(path.read_text())
        if row["label"]=="hill":
            aa,bb=[data["jacobian_model"][k] for k in ["a4_low_to_high","a6_low_to_high"]]
        else:aa,bb=data["short_A_low_to_high"],data["short_B_low_to_high"]
        A,B=reduction(aa,p),reduction(bb,p)
        C,f,trace=critical_trace(A,B,p)
        require(C==row["critical_polynomial_monic"] and f==row["j_in_critical_algebra"],"critical algebra")
        require(trace==row["critical_value_trace"],"critical-value trace")
        # Newton identities recover all coefficients from exact quotient traces.
        # Since p>24, divisions by1,...,18 are valid.
        powers=[];g=[1]
        for k in range(1,19):
            g=rem(mul(g,f,p),C,p)
            tr=0
            for i in range(18):
                col=rem([0]*i+g,C,p)
                tr=(tr+(col[i] if i<len(col) else 0))%p
            powers.append(tr)
        coefficients=[1]
        for k in range(1,19):
            coefficients.append(-sum(coefficients[k-i]*powers[i-1]
                                     for i in range(1,k+1))*pow(k,-1,p)%p)
        require(list(reversed(coefficients))==row["critical_value_characteristic"],
                "complete critical-value characteristic polynomial")
        count=surface_count(A,B,p)
        require(count==row["surface_Fp_points"] and count-1-p*p==row["surface_H2_trace"],"surface count")
        traces.append(trace);counts.append(count)
    require(traces[0]!=traces[1],"no distinction supplied")
    require(result["surface_identity"]=="UNKNOWN","surface identity overclaim")
    return {"schema":"x1152-critical-value-independent-replay-v2","status":"PASS",
        "prime":p,"critical_value_traces":traces,"surface_Fp_points":counts,
        "geometric_marked_fibrations_distinct":True,
        "underlying_surface_identity":"UNKNOWN",
        "boundary":"Intrinsic j-critical divisor invariant only. No full rank17 record basis, surface isomorphism/nonisomorphism or Picard proof."}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--certificate",type=Path,default=DEFAULT)
    p.add_argument("--output",type=Path);a=p.parse_args()
    result=verify(json.loads(a.certificate.read_text()))
    if a.output:
        with a.output.open("x") as f:json.dump(result,f,indent=2,sort_keys=True);f.write("\n")
    print(json.dumps(result,sort_keys=True))
