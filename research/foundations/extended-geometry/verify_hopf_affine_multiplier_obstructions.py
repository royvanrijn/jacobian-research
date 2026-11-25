#!/usr/bin/env python3
"""Exact finite identity checks for two affine Hopf-multiplier exclusions.

The unrestricted statements follow from the written necessary-equation proofs
and separate endpoint converses. This program checks the displayed finite
polynomials and certificates; it does not infer an infinite moment tail.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from math import comb, factorial
from pathlib import Path
import json
import platform

def sparse_mul(p, q):
    out=defaultdict(F)
    for k,c in p.items():
        for ell,d in q.items():
            out[tuple(a+b for a,b in zip(k,ell))]+=c*d
    return {k:v for k,v in out.items() if v}

def trim(p):
    p=list(map(F,p))
    while len(p)>1 and not p[-1]:p.pop()
    return p or [F(0)]

def add(p,q):
    return trim([(p[i] if i<len(p) else 0)+(q[i] if i<len(q) else 0)
                 for i in range(max(len(p),len(q)))])

def scale(p,c):
    return trim([v*c for v in p])

def mul(p,q):
    out=[F(0)]*(len(p)+len(q)-1)
    for i,c in enumerate(p):
        for j,d in enumerate(q):out[i+j]+=c*d
    return trim(out)

def power(p,n):
    out=[F(1)]
    for _ in range(n):out=mul(out,p)
    return out

def check(condition, label):
    if not condition:raise ArithmeticError(label)

# Keys: exponents (x,t,a,b,c,d). These are literal Laurent inputs.
P4={(-1,0,0,0,0,0):F(1),(0,0,0,0,0,0):F(1),
    (-1,2,0,0,0,0):F(-1),(0,2,0,0,0,0):F(-3),
    (1,2,0,0,0,0):F(-3),(2,2,0,0,0,0):F(-1)}
UT={(0,1,0,0,0,0):F(1),(1,1,0,0,0,0):F(1)}
P6=sparse_mul(P4,UT)
ELL={(0,0,1,0,0,0):F(1),(1,0,0,1,0,0):F(1),
     (-1,0,0,0,1,0):F(1,2),(-1,2,0,0,1,0):F(-1,2),
     (0,1,0,0,0,1):F(1)}

def raw_moments(seed, wanted):
    f=sparse_mul(seed,ELL)
    p={(0,0,0,0,0,0):F(1)}
    out={}
    for m in range(1,max(wanted)+1):
        p=sparse_mul(p,f)
        if m not in wanted:continue
        moment=defaultdict(F)
        for (x,t,a,b,c,d),v in p.items():
            if x==0 and t%2==0:moment[(a,b,c,d)]+=v/F(t+1)
        out[m]={k:v for k,v in moment.items() if v}
    return out

def angular_closed(m,i,j,h,sextic):
    parity=(m+h) if sextic else h
    if parity%2:return F(0)
    lower=m-i+j
    total=F(0)
    for v in range(m+1):
        upper=(2*m if sextic else m)+2*v
        if not 0<=lower<=upper:continue
        height=sum((F((-1)**w*comb(j,w),
                      (m if sextic else 0)+2*v+h+2*w+1)
                    for w in range(j+1)),F(0))
        total+=(-1)**v*comb(m,v)*comb(upper,lower)*height
    return total/F(2**j)

def closed_moment(m,sextic):
    out={}
    for b in range(m+1):
        for c in range(m-b+1):
            for d in range(m-b-c+1):
                a=m-b-c-d
                coefficient=F(factorial(m),
                    factorial(a)*factorial(b)*factorial(c)*factorial(d))
                coefficient*=angular_closed(m,b,c,d,sextic)
                if coefficient:out[(a,b,c,d)]=coefficient
    return out

Pq=list(map(F,[14254,56680,32175]))
Qq=list(map(F,[-643,1802,79220,60775]))
Uq=list(map(F,[290136844604,-1247587548025,-1531111624900]))
Vq=list(map(F,[1031834055265,810588507300]))
QUARTIC_BEZOUT=3472141285450021
Ps=list(map(F,[6731307,167187968,-485630508,152847024,7148449616,
               -13319044992,-5200227648,2971558656]))
Qs=list(map(F,[5677522608208,-446879603799373,2571951562119874,
               -19554621670633104,-13510803549968736,423512030238845856,
               -786969297304179264,-1083906565862605056,
               4002560438537396736,804612135850617600,
               -2119204298356231680,596608834727854080]))
Us=list(map(F,[24,81,76,19,51,70,47,30,93,1,69]))
Vs=list(map(F,[68,96,31,6,17,98,9]))

def check_integer_bezout(rhs):
    check(add(mul(Uq,Pq),mul(Vq,Qq))==[F(rhs)],
          "quartic integer Bezout identity")

def check_mod_bezout(p,q,u,v):
    check(all(x.denominator==1 for x in p+q+u+v),"integer modular inputs")
    check(p[-1].numerator%101!=0 and q[-1].numerator%101!=0,
          "leading-coefficient modulo-101 units")
    lhs=add(mul(u,p),mul(v,q))
    check([x.numerator%101 for x in lhs]==[1]+[0]*(len(lhs)-1),
          "sextic modulo-101 Bezout identity")

def expect_rejection(name, action):
    try:action()
    except ArithmeticError:return name
    raise ArithmeticError("corrupted control was accepted: "+name)

def q_normalize(poly):
    out=[F(0)]
    v=list(map(F,[F(-27,25),F(-130,25)]))
    for (a,b,c,d),co in poly.items():
        check(d%2==0,"quartic even d parity")
        term=mul([F(0)]*a+[F(1)],power(v,d//2))
        out=add(out,scale(term,co*F(3,10)**b))
    return out

def s_clear_denominator(poly,m):
    # c=1,b=z,a=-(132z^2-44z+53)/(44(6z+1)),d=0.
    numerator=list(map(F,[-53,44,-132]))
    denominator=list(map(F,[44,264]))
    out=[F(0)]
    for (a,b,c,d),co in poly.items():
        if d:continue
        term=mul(power(numerator,a),power(denominator,m-a))
        term=[F(0)]*b+term
        out=add(out,scale(term,co))
    return out

def run():
    q=raw_moments(P4,{1,2,3,4})
    s=raw_moments(P6,{1,2,4,6})
    for m,p in q.items():check(p==closed_moment(m,False),
                              "quartic raw/closed moment "+str(m))
    for m,p in s.items():check(p==closed_moment(m,True),
                              "sextic raw/closed moment "+str(m))
    check(q[1]=={(0,1,0,0):F(2,3),(0,0,1,0):F(-1,5)},"quartic M1")
    check(q[2]=={(1,1,0,0):F(16,15),(1,0,1,0):F(8,105),
                 (0,2,0,0):F(8,15),(0,1,1,0):F(-8,105),
                 (0,0,2,0):F(2,35),(0,0,0,2):F(8,105)},"quartic M2")
    check(q_normalize(q[3])==scale(Pq,F(8,625625)),"quartic reduced M3")
    check(q_normalize(q[4])==scale(Qq,F(256,31906875)),"quartic reduced M4")
    check(s[1]=={(0,0,0,1):F(-2,15)},"sextic M1")
    s2={k:v for k,v in s[2].items() if not k[3]}
    check(s2=={(1,1,0,0):F(16,105),(1,0,1,0):F(8,315),
               (0,2,0,0):F(8,105),(0,1,1,0):F(-8,315),
               (0,0,2,0):F(106,3465)},"sextic M2 at d=0")
    # Denominator guard at z=-1/6 and the c=0, b=-2a branch.
    check(F(132,36)+F(44,6)+53==64,"sextic denominator guard")
    branch=sum((co*(-2)**b for (a,b,c,d),co in s[4].items()
                if c==d==0),F(0))
    check(branch==F(-2048,15015),"sextic c=0 fourth-moment exclusion")
    w=[F(1),F(6)]
    check(s_clear_denominator(s[4],4)==scale(mul(w,Ps),
          F(-4*44**4,6455143695)),"sextic reduced M4")
    check(s_clear_denominator(s[6],6)==scale(mul(w,Qs),
          F(2*44**6,48450701243889945)),"sextic reduced M6")
    check(Ps[-1].numerator%101==84 and Qs[-1].numerator%101==63,
          "displayed leading residues")
    check_integer_bezout(QUARTIC_BEZOUT)
    check_mod_bezout(Ps,Qs,Us,Vs)
    bad_u=Us.copy();bad_u[0]+=1
    bad_p=Ps.copy();bad_p[-1]=F(101)
    rejected=[
        expect_rejection("wrong quartic Bezout constant",
                         lambda:check_integer_bezout(QUARTIC_BEZOUT+1)),
        expect_rejection("corrupted sextic Bezout multiplier",
                         lambda:check_mod_bezout(Ps,Qs,bad_u,Vs)),
        expect_rejection("nonunit leading coefficient",
                         lambda:check_mod_bezout(bad_p,Qs,Us,Vs))]
    def encode(p):
        return [{"exponents":list(k),"coefficient":str(v)}
                for k,v in sorted(p.items())]
    return {"status":"PASS","quartic_orders":sorted(q),"sextic_orders":sorted(s),
            "quartic_moments":{str(m):encode(p) for m,p in q.items()},
            "sextic_moments":{str(m):encode(p) for m,p in s.items()},
            "corrupted_controls_rejected":rejected,
            "scope":"Exact finite rational moment identities, necessary-equation reductions and Bezout certificates. Unrestricted purity converses and all-field/degree statements are written proofs.",
            "formal_verification":False,"independent_replay":False,
            "external_review":False}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output",type=Path)
    args=ap.parse_args()
    receipt=run()
    receipt["checker_sha256"]=sha256(Path(__file__).read_bytes()).hexdigest()
    receipt["python_version"]=platform.python_version()
    if args.output:
        args.output.write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n")
    print("PASS: quartic orders 1,2,3,4; sextic orders 1,2,4,6; both necessary-equation certificates; three corrupted controls rejected.")

if __name__=="__main__":main()
