#!/usr/bin/env python3
"""Exact finite admission controls for the all-orders Painleve-I Long construction."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from math import comb, factorial
import importlib.util
import hashlib
import json
import sympy as s

HERE=Path(__file__).resolve().parent
# Local PBW exponent order: (R, X, B, h), with [B,X]=3h X^2 and R central.


def mon(a=0,b=0,c=0,d=0,k=1):
    return {(a,b,c,d):F(k)} if k else {}


def add(*polys):
    out=defaultdict(F)
    for p in polys:
        for e,c in p.items():
            out[e]+=c
    return {e:c for e,c in out.items() if c}


def scale(p,k):
    return {e:c*k for e,c in p.items() if c*k}


def mul(p,q,max_h=None):
    out=defaultdict(F)
    for (a,b,c,d),v in p.items():
        for (A,B,C,D),w in q.items():
            rising=1
            for j in range(c+1):
                if j:
                    rising*=B+j-1
                if max_h is not None and d+D+j>max_h:
                    continue
                value=v*w*comb(c,j)*3**j*rising
                if value:
                    out[(a+A,b+B+j,c+C-j,d+D+j)]+=value
    return {e:c for e,c in out.items() if c}


def power(p,n):
    out=mon()
    for _ in range(n):
        out=mul(out,p)
    return out


def require(value,label):
    if not value:
        raise AssertionError(label)


def equal(p,q,label):
    require(not add(p,scale(q,-1)),label)


def trunc(p,n):
    return {e:c for e,c in p.items() if e[3]<=n}


def comm(p,q,n=None):
    return add(mul(p,q,n),scale(mul(q,p,n),-1))


def globalize(p):
    require(all(e[0]>=0 for e in p),"no negative R powers")
    out=defaultdict(F)
    for (a,b,c,d),coef in p.items():
        for j in range(a+1):
            out[(a+b+j,j,c,d)]+=coef*comb(a,j)*2**(a-j)*(-3)**j
    out={e:c for e,c in out.items() if c}
    require(all(e[0]>=0 for e in out),"no negative original X powers")
    return out


def rows(p):
    return [[list(e),str(c)] for e,c in sorted(p.items())]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    R=mon(a=1)
    X=mon(b=1)
    BB=mon(c=1)
    h=mon(d=1)
    a0=mon(a=1,b=-1)
    w=add(power(a0,3),scale(power(a0,2),-5),scale(a0,9),mon(k=-2),
          mon(a=1,b=1,c=1,k=-1),mon(a=1,b=2,d=1,k=-F(3,2)))
    xx=scale(a0,8)
    U0=add(scale(power(w,2),4),scale(xx,-3))
    V0=add(scale(power(w,3),8),scale(add(mul(xx,w),mul(w,xx)),-F(9,2)))
    Fn=[mon(a=2,b=4,d=2,k=6),mon(a=3,b=9,d=4,k=F(147,4))]
    U=add(U0,*Fn)
    V=add(V0,*[scale(add(mul(f,w),mul(w,f)),F(3,2)) for f in Fn])
    SS=mul(add(mon(k=16),scale(U,-1)),mon(a=-1,k=F(1,12)))
    TT=mul(add(V,scale(U,6),mon(k=-32)),mon(a=-2,k=F(1,432)))
    CK=add(mon(b=-2,c=1,k=-F(1,2)),mon(a=1,b=-3,c=1,k=F(1,6)),
        mon(a=-1,b=-1,c=1,k=-F(1,15)),mon(b=-4,k=-F(3,4)),
        mon(a=1,b=-5,k=F(1,5)),mon(a=-1,b=-3,k=F(3,5)),
        mon(a=-2,b=-2,k=F(1,5)),mon(a=-1,d=1,k=F(1,10)))
    AA=mul(add(scale(mul(R,add(mul(SS,TT),mul(TT,SS))),2),
               scale(TT,-32),scale(mul(SS,SS),-1)),mon(a=-2,k=F(1,20)))
    DD=add(CK,AA)
    glob={name:globalize(p) for name,p in [("S",SS),("T",TT),("D_minus_D0",DD)]}

    # The original polynomial D0 gives this exact derivation of the localized C algebra.
    dB=add(mon(b=-1,c=1,k=3),mon(a=1,b=-2,c=1,k=-F(3,2)))
    cache={0:{}}
    def delta_Bpower(n):
        if n not in cache:
            cache[n]=add(*[mul(mul(power(BB,j),dB),power(BB,n-1-j))
                           for j in range(n)])
        return cache[n]
    def delta(p):
        out={}
        for (a,b,c,d),coef in p.items():
            posdelta=add(mon(a=a-1,b=b,d=d,k=a),
                         mon(a=a,b=b-1,d=d,k=F(3*b,2)),
                         mon(a=a+1,b=b-2,d=d,k=-F(b,2)))
            term=add(mul(posdelta,power(BB,c)),
                     mul(mon(a=a,b=b,d=d),delta_Bpower(c)))
            out=add(out,scale(term,coef))
        return out

    equal(add(mul(h,delta(w)),comm(CK,w)),
          mul(mon(a=-1,d=1,k=F(3,5)),w),"exact Euler action on w")
    equal(add(mul(h,delta(xx)),comm(CK,xx)),
          mul(mon(a=-1,d=1,k=F(6,5)),xx),"exact Euler action on input x")
    # [R,S]=[R,T]=0 is exact because R is central in C.
    equal(comm(SS,TT,6),h,"pair relation through h6")
    equal(trunc(add(mul(h,delta(SS)),comm(DD,SS,6)),6),{},
          "full D-S relation through h6")
    equal(trunc(add(mul(h,delta(TT)),comm(DD,TT,6)),6),{},
          "full D-T relation through h6")
    equal(delta(R),mon(),"D-R relation")
    require(max(e[2] for e in glob["S"])==2,"S momentum order")
    require(max(e[2] for e in glob["T"])==3,"T momentum order")
    require(max(e[2] for e in glob["D_minus_D0"])==5,"D momentum order")

    # A single generic string identity for arbitrary W(t), not a high-order sweep.
    t,q,ee=s.symbols("t q ee")
    W=s.Function("W")(t)
    LP=4*q**2+W
    LQ=8*q**3+3*W*q
    bracket=s.diff(LQ,q)*s.diff(LP,t)-s.diff(LQ,t)*s.diff(LP,q)
    cubic=s.diff(LQ,q,3)*s.diff(LP,t,3)
    require(s.simplify(bracket-3*W*s.diff(W,t))==0,"generic string first term")
    require(s.simplify(cubic/24-2*s.diff(W,t,3))==0,"generic string dispersive term")
    def generic_star(f,g):
        degree=int(s.degree(f,q)+s.degree(g,q))
        out=0
        for n in range(degree+1):
            term=sum((-1)**j*comb(n,j)*s.diff(f,q,n-j,t,j)
                     *s.diff(g,t,n-j,q,j) for j in range(n+1))
            out+=ee**n*term/(2**n*factorial(n))
        return s.expand(out)
    cubic_energy=s.expand(generic_star(LQ,LQ)-generic_star(generic_star(LP,LP),LP))
    require(s.degree(cubic_energy,q)==2,"generic energy remains quadratic")
    require(s.simplify(cubic_energy.coeff(q,2)+3*W**2+4*ee**2*s.diff(W,t,2))==0,
            "exact quadratic energy coefficient")
    coeff=[F(-3)]
    for n in (1,2):
        conv=sum((coeff[j]*coeff[n-j] for j in range(1,n)),F(0))
        coeff.append(conv/6+F(2,9)*(25*(n-1)**2-1)*coeff[n-1])
    require(coeff==[F(-3),F(2,3),F(98,27)],"forced first PI coefficients")
    for n in (1,2):
        const=8*coeff[n]*F(9,8)**n
        equal(Fn[n-1],mon(a=n+1,b=5*n-1,d=2*n,k=const),
              "global scaling clears PI denominators")

    # Compare the whole classical tuple with the actual retained Long formulas.
    source=HERE.parent/"dixmier-audit"/"verify_long_pbw.py"
    require(sha(source)=="266cfed8261bdec371bd26e2bb3a46cab9f588030adcc92c4e75dbd6b150bca2",
            "retained Long source hash")
    spec=importlib.util.spec_from_file_location("retained_long",source)
    audit=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
    xo,qo,po,zo=(audit.var(i) for i in range(4))
    Bclass=audit.add(audit.scale(audit.mul(audit.power(xo,2),po),3),
        audit.scale(audit.mul(audit.add(audit.const(1),
            audit.scale(audit.mul(xo,qo),-3)),zo),2))
    def classical(p):
        out={}
        for (a,b,c,d),coef in p.items():
            if d==0:
                term=audit.mul(audit.mul(audit.power(xo,a),audit.power(qo,b)),
                               audit.power(Bclass,c))
                out=audit.add(out,audit.scale(term,coef))
        return out
    D0=audit.add(audit.scale(audit.mul(
        audit.add(audit.const(1),audit.scale(audit.mul(xo,qo),3)),po),F(1,2)),
        audit.scale(audit.mul(audit.power(qo,2),zo),-3))
    actual=audit.long_images()
    require(classical(glob["S"])==actual[3],"complete classical S")
    require(classical(glob["T"])==actual[1],"complete classical T")
    require(audit.add(D0,classical(glob["D_minus_D0"]))==actual[2],
            "complete corrected classical D")

    # Check the advertised exact S coefficients; all n>=2 follow the written formula.
    s2={e:c for e,c in glob["S"].items() if e[3]==2}
    expected2=globalize(mon(a=1,b=4,d=2,k=-F(17,4)))
    equal(s2,expected2,"S2 low-order exception")
    s4={e:c for e,c in glob["S"].items() if e[3]==4}
    equal(s4,globalize(mon(a=2,b=9,d=4,k=-F(49,16))),"S4 general formula")

    packet={"claim":"QPL1","scope":"finite admission seed, containing the base and forced PI n=1,2 terms",
            "local_exponent_order":["R","X","B","h"],
            "global_exponent_order":["X","q","B","h"],
            "local":{name:rows(p) for name,p in [("S",SS),("T",TT),("D_minus_D0",DD)]},
            "global":{name:rows(p) for name,p in glob.items()}}
    (HERE/"painleve_long_seed.json").write_text(json.dumps(packet,indent=2)+"\n")
    locked=["PAINLEVE_GLOBAL_LONG_LIFT.md","check_painleve_long.py",
            "painleve_long_seed.json","ELLIPTIC_ISOGENY_AND_QUANTUM_ENERGY.md",
            "ELLIPTIC_QUANTUM_ENERGY_MOTION.md"]
    out={"claim":"QPL1","status":"PASS",
        "admission_seed_global_term_counts":{name:len(p) for name,p in glob.items()},
        "full_classical_retained_Long_tuple":"PASS",
        "all_six_full_relations_through_h6":"PASS",
        "exact_Euler_connection_identities":"PASS",
        "generic_PI_string_identity":"PASS",
        "forced_PI_coefficients":[str(a) for a in coeff],
        "uniform_momentum_orders":{"R":0,"S":2,"T":3,"D":5},
        "formal_cubic_energy_scalar_order":2,
        "formal_cubic_energy_leading_coefficient":"-54*t",
        "S_coefficient_degree_formula":"8*n-1 for every n>=1; nonzero coefficients",
        "S2_exception":"-(17/4)*R*X^4",
        "S4":"-(49/16)*R^2*X^9",
        "hashes":{name:sha(HERE/name) for name in locked},
        "retained_Long_source_sha256":sha(source),
        "coverage_boundary":"Exact finite admission seed and symbolic identities only. The all-orders PI construction, positivity, global admission beyond the seed, full relations in all orders, exact degree growth and divergence use written proofs. This constructs an unbounded full formal Long lift with bounded momentum orders, not a finite Weyl endomorphism or all-representative nonlifting theorem.",
        "independent_replay":False,"formal_verification":False,"external_review":False}
    (HERE/"painleve_long_controls.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:out[k] for k in ["claim","status","admission_seed_global_term_counts",
        "full_classical_retained_Long_tuple","uniform_momentum_orders"]},indent=2))


if __name__=="__main__":
    main()
