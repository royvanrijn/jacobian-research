#!/usr/bin/env python3
"""Exact finite certificate controls for EOP1; Python standard library."""
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt, prod
from pathlib import Path
import json
import platform

HERE=Path(__file__).resolve().parent
CHECKS=[]


def check(label,condition):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def valuation(q,p):
    q=F(q)
    if not q:
        return float("inf")
    a,b=abs(q.numerator),q.denominator
    v=0
    while a%p==0:
        a//=p
        v+=1
    while b%p==0:
        b//=p
        v-=1
    return v


def modq(q,m):
    q=F(q)
    return q.numerator*pow(q.denominator,-1,m)%m


def prime(n):
    if n<2:
        return False
    if n%2==0:
        return n==2
    return all(n%d for d in range(3,isqrt(n)+1,2))


def square(q):
    q=F(q)
    return q>=0 and isqrt(q.numerator)**2==q.numerator and isqrt(q.denominator)**2==q.denominator


def on_curve(P,a,b):
    return P is None or P[1]**2==P[0]**3+a*P[0]**2+b*P[0]


def add(P,Q,a):
    if P is None:
        return Q
    if Q is None:
        return P
    x,y=P
    u,v=Q
    if x==u and y==-v:
        return None
    if P==Q:
        slope=(3*x*x+2*a*x+27091730)/(2*y)
    else:
        slope=(v-y)/(u-x)
    xx=slope*slope-a-x-u
    return xx,-y+slope*(x-xx)


def mul(P,n,a):
    out=None
    while n:
        if n%2:
            out=add(out,P,a)
        P=add(P,P,a)
        n//=2
    return out


def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:
        p.pop()
    return p


def remainder(f,g):
    f=trim(f)
    g=trim(g)
    while len(f)>=len(g) and any(f):
        c=f[-1]/g[-1]
        offset=len(f)-len(g)
        for i,v in enumerate(g):
            f[i+offset]-=c*v
        f=trim(f)
    return f


def gcd_degree(f,g):
    f,g=trim(f),trim(g)
    while any(g):
        f,g=g,remainder(f,g)
    return len(trim(f))-1


def discriminant(a,b,c):
    return a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c


def factors_at(P):
    x,y=P
    v=F(3,6400)*(x+432)
    s0=F(1816347,819200)
    q=[v,F(9,10),F(1)]
    h=[F(10,9)*v*v-F(10,9)*s0*v+F(125,1024),s0-v,-F(12,5),F(1)]
    D=q[1]**2-4*q[0]
    Delta=discriminant(h[2],h[1],h[0])
    return q,h,D,Delta


def depth(P,p):
    return valuation(-(P[0]-3552)/P[1],p)


def replay():
    cert=json.loads((HERE/"elliptic_exception_one_prime_certificate.json").read_text())
    a,b=cert["curve"]["a2"],cert["curve"]["a4"]
    check("fixed integral model",(a,b,cert["curve"]["a6"])==(-10656,27091730,0))
    points={k:tuple(F(q) for q in P) for k,P in cert["points"].items()}
    for name,P in points.items():
        check("exact curve point "+name,on_curve(P,a,b))
    points["R"]=add(points["P"],points["T"],a)
    points["R_minus_S"]=add(points["R"],(points["S"][0],-points["S"][1]),a)
    for name in ["R","R_minus_S"]:
        check("exact derived curve point "+name,on_curve(points[name],a,b))
    check("rational two-torsion",add(points["T"],points["T"],a) is None)
    for name,alpha in cert["classes"].items():
        if name=="T":
            check("two-torsion squareclass",alpha==b)
        else:
            check("exact global squareclass "+name,square(points[name][0]/alpha))

    fourP=mul(points["P"],4,a)
    check("four-P exact curve equation",on_curve(fourP,a,b))
    fac=cert["four_P_x_denominator_factors"]
    check("distinct complete denominator factors",len({p for p,e in fac})==len(fac))
    for p,e in fac:
        check("denominator factor is prime %d"%p,prime(p) and e>0)
    product=prod(p**e for p,e in fac)
    check("exact complete four-P denominator factorization",fourP[0].denominator==product)
    exceptions=sorted(p for p,e in fac if p>=11 and p%3==2)
    check("complete inert denominator exceptions",exceptions==cert["inert_good_exceptions"])
    q0=cert["q0"]
    check("q0 exact prime",q0==2709173 and prime(q0))
    discE=16*b*b*(a*a-4*b)
    check("complete elliptic discriminant factorization",discE==2**9*5**2*7**3*1889*q0*q0)
    for p in exceptions:
        check("exception prime is good %d"%p,discE%p!=0)
    check("additional inert bad primes",cert["additional_inert_bad_primes"]==[1889,q0])
    check("1889 exact prime",prime(1889) and 1889%3==2)
    for p in exceptions+[1889,q0]:
        for name in ["P","R"]:
            check("exception witness unit %s at %d"%(name,p),valuation(128*points[name][0]-25865,p)==0)

    for name in ["P","Q","R","R_minus_S"]:
        q,h,D,Delta=factors_at(points[name])
        check("witness quadratic discriminant identity "+name,D==-F(3,1600)*points[name][0])
        check("witness quadratic separable "+name,D!=0)
        check("witness cubic separable "+name,Delta!=0)
        check("witness factors coprime "+name,gcd_degree(q,h)==0)
        check("witness matching discriminant classes "+name,square(Delta/D))

    check("class-two two-adic translated depth",depth(points["R"],2)==4)
    check("class-two two-adic positive coordinate",valuation(points["P"][0],2)==9)
    check("class-two three-adic residue and depth",modq(points["Q"][0]-3552,3)==2 and valuation(points["Q"][1],3)==2)
    check("class-two five-adic integral coordinate",valuation(points["P"][0],5)>=0)
    check("class-five-q two-adic x valuation",valuation(points["R_minus_S"][0],2)==-12)
    check("class-five-q two-adic formal depth",depth(points["R_minus_S"],2)==6)
    tr=-(points["R"][0]-3552)/points["R"][1]
    check("class-five-q three-adic formal depth",valuation(tr,3)==3)
    check("class-five-q three-adic cube unit",modq(tr/27,9)==8)
    check("class-five-q five-adic integral coordinate",valuation(points["R"][0],5)>=0)
    corrupted=fac.copy()
    corrupted[0]=[2,13]
    check("corrupted factor exponent rejected",prod(p**e for p,e in corrupted)!=fourP[0].denominator)
    wrong=(points["P"][0]+1,points["P"][1])
    check("corrupted rational point rejected",not on_curve(wrong,a,b))
    return {"passed":True,"exact_controls":len(CHECKS),"checks":CHECKS,
            "four_P_x_denominator":str(product),
            "derived_points":{k:[str(q) for q in points[k]] for k in ["R","R_minus_S"]},
            "complete_exceptional_inert_primes":exceptions+[1889,q0],
            "unrestricted_proof":"EOP1 three-residue argument, EOP2 profinite coset argument, and EOP3 canonical-height reduction-image argument",
            "bounded_checks_establish_all_primes_by_extrapolation":False,
            "all_finite_sets_in_other_two_classes_claimed":False,
            "new_all_places_Hasse_seed_claimed":False,
            "formal_verification":False,"independent_replay":False,"external_review":False,
            "python":platform.python_version(),"arithmetic":"standard-library exact integers and fractions"}


if __name__=="__main__":
    report=replay()
    names=["ELLIPTIC_EXCEPTION_ONE_PRIME_POSITIVITY.md",Path(__file__).name,
           "elliptic_exception_one_prime_certificate.json"]
    report["input_sha256"]={name:sha256((HERE/name).read_bytes()).hexdigest() for name in names}
    (HERE/"elliptic_exception_one_prime_positivity_checks.json").write_text(
        json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"passed":True,"exact_controls":report["exact_controls"],
                      "input_sha256":report["input_sha256"]},indent=2))
