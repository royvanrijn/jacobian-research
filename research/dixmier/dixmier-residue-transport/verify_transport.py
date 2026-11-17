#!/usr/bin/env python3
"""Exact symbolic-lambda transport of the retained strong order-five residue.

Sparse arithmetic below uses Fraction and four exponents (X,Q,Z,lambda).
The source constructor and 42-dimensional lower-lift solve are retained code.
No missing artifact is rebuilt and no parameter search is performed.
"""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OLD=ROOT/"research/foundations/scripts"
sys.path.insert(0,str(OLD))
from sympy.polys.domains import QQ
from explore_degree_five_quantum_residue import degree_five_family, third_order_family, laurent_monomials, solve_affine, column_rank
from verify_degree_five_laurent_quantum_obstruction import SUPPORT, COEFFICIENTS

def clean(p): return {m:c for m,c in p.items() if c}
def add(*ps):
    out=defaultdict(F)
    for p in ps:
        for m,c in p.items(): out[m]+=c
    return clean(out)
def scale(p,c): return clean({m:c*v for m,v in p.items()})
def mul(p,q):
    out=defaultdict(F)
    for m,c in p.items():
        for n,d in q.items(): out[tuple(x+y for x,y in zip(m,n))]+=c*d
    return clean(out)
def dz(p,n=1):
    for _ in range(n): p={(i,j,k-1,l):c*k for (i,j,k,l),c in p.items() if k}
    return p
def delta(p,n=1):
    for _ in range(n):
        out=defaultdict(F)
        for (i,j,k,l),c in p.items():
            out[i+1,j,k,l]+=c*(3*i-6*j)
            if j: out[i,j-1,k,l]+=c*2*j
        p=clean(out)
    return p
def pi(p,q,n=1):
    out={}
    for j in range(n+1):
        out=add(out,scale(mul(delta(dz(p,n-j),j),dz(delta(q,n-j),j)),(-1)**j*comb(n,j)))
    return out
def convert(p): return {m+(0,):F(str(c)) for m,c in p.items()}
def shear(p,sign=1):
    # tau: Z -> Z + 3*lambda*X^4.
    out=defaultdict(F)
    for (i,j,k,l),c in p.items():
        for r in range(k+1):
            out[i+4*r,j,k-r,l+r]+=c*comb(k,r)*(sign*3)**r
    return clean(out)
A2={(6,0,0,1):F(45,2)}
A4={(8,0,0,1):F(1701,16)}
def alpha2(p): return shear(mul(A2,dz(p,3)))
def alpha4(p):
    return shear(add(mul(A4,dz(p,5)),scale(mul(mul(A2,A2),dz(p,6)),F(1,2))))
def q2(seed,lower): return add(shear(lower),alpha2(seed))
def q4(seed,lower): return add(alpha2(lower),alpha4(seed))
def defect(S,T,s,t):
    return add(pi(s,t),scale(pi(s,T,3),F(1,24)),scale(pi(S,t,3),F(1,24)),scale(pi(S,T,5),F(1,1920)))
def lin(S,T,s,t): return add(pi(s,T),pi(S,t))
def transport_dual(coeffs,zcap=4):
    out=defaultdict(F)
    # Pull back each coefficient by tau^-1. This is complete for ALL X
    # exponents (including negative ones) and 0<=Z-degree<=zcap.
    for (i,j,k),c in coeffs.items():
        for K in range(k,zcap+1):
            r=K-k
            out[i-4*r,j,K,r]+=c*comb(K,k)*(-3)**r
    return clean(out)
def pair(dual,p):
    out=defaultdict(F)
    index=defaultdict(list)
    for (i,j,k,l),c in dual.items(): index[i,j,k].append((l,c))
    for (i,j,k,l),c in p.items():
        for r,d in index[i,j,k]: out[l+r]+=c*d
    return clean(out)
def serial(p): return [[list(m) if isinstance(m,tuple) else m,str(c)] for m,c in sorted(p.items())]
def digest(p): return hashlib.sha256(json.dumps(serial(p),separators=(",",":")).encode()).hexdigest()
def expect(c,msg):
    if not c: raise AssertionError(msg)

def main():
    start=time.monotonic()
    S0,T0=degree_five_family(QQ,-QQ(1)/2,QQ(1))
    family=third_order_family(S0,T0,QQ)
    expect((family.column_count,family.operator_rank,len(family.kernel))==(1569,1527,42),"full lower-lift solve")
    S,T=convert(S0),convert(T0)
    s,t=map(convert,family.base)
    ks=[tuple(map(convert,row)) for row in family.kernel]
    coeffs={m:F(*c) for m,c in zip(SUPPORT,COEFFICIENTS)}
    dual0={m+(0,):c for m,c in coeffs.items()}
    dual=transport_dual(coeffs)
    expected={0:-F(47547660815739,190658)}
    expect(pi(S,T)=={(0,0,0,0):F(1)},"canonical seed")
    expect((max(m[2] for m in S),max(m[2] for m in T))==(5,4),"source Z bounds")
    expect(lin(S,T,s,t)==scale(pi(S,T,3),-F(1,24)),"lower equation")
    expect(shear(shear(S),-1)==S,"symbolic inverse")
    # Exact coefficient checks of the entire bounded affine lower-lift scheme.
    O=defect(S,T,s,t)
    expect(pair(dual0,O)==expected,"retained scalar value")
    expect(pair(dual,shear(O))==expected,"transported scalar value")
    nonlinear_checks=0
    for u,v in ks:
        diag=pi(u,v)
        linear=add(pi(s,v),pi(u,t),scale(pi(u,T,3),F(1,24)),scale(pi(S,v,3),F(1,24)))
        for col in (linear,diag):
            expect(not pair(dual0,col),"old lower-lift variation")
            expect(not pair(dual,shear(col)),"transported lower-lift variation")
            nonlinear_checks+=1
    for i,j in combinations_with_replacement(range(42),2):
        if i==j: continue
        col=add(pi(ks[i][0],ks[j][1]),pi(ks[j][0],ks[i][1]))
        expect(not pair(dual0,col),"old cross variation")
        expect(not pair(dual,shear(col)),"transported cross variation")
        nonlinear_checks+=1
    # Current corrections: transport every column, rather than imposing the
    # same numerical degree cutoffs again in the new chart.
    counts=[]
    leakage=[]
    for bound,zcap,which in ((21,1,0),(17,0,1)):
        mons=laurent_monomials(bound,zcap,0,3)
        counts.append(len(mons))
        for m in mons:
            p={m+(0,):F(1)}
            col=pi(p,T) if which==0 else pi(S,p)
            expect(not pair(dual0,col),"all old current columns")
            expect(not pair(dual,shear(col)),"all transported current columns")
            if any(i+j+3*k>bound for i,j,k,l in shear(p)):
                leakage.append((which,m))
    # Independently expand the actual transformed hbar equations, including
    # the mandatory affine higher-order shift from the quantum automorphism.
    Sh,Th=shear(S),shear(T)
    # Separate test: reimposing the old numerical box is NOT the transport.
    # Record whether the transported cocycle happens to kill that box too.
    native_failures=[]
    for bound,zcap,which in ((21,1,0),(17,0,1)):
        for m in laurent_monomials(bound,zcap,0,3):
            p={m+(0,):F(1)}
            col=pi(p,Th) if which==0 else pi(Sh,p)
            value=pair(dual,col)
            if value: native_failures.append({"component":which,"monomial":list(m),"value":serial(value)})
    # Complete relevance rectangle for the original current-order operator
    # with only its Z-order restriction retained (no Bernstein-degree bound).
    # Negative X inputs cannot reach X>=14 (shift<=12); Q-input<=8 and
    # X-input<=19 cover every nonnegative input meeting dual support.
    wider_failures=[]
    wider_count=0
    for which,zcap in ((0,1),(1,0)):
        for i in range(20):
            for j in range(9):
                for k in range(zcap+1):
                    p={(i,j,k,0):F(1)}
                    value=pair(dual0,pi(p,T) if which==0 else pi(S,p))
                    wider_count+=1
                    if value: wider_failures.append({"component":which,"monomial":[i,j,k],"value":serial(value)})
    sh,th=q2(S,s),q2(T,t)
    expect(lin(Sh,Th,sh,th)==scale(pi(Sh,Th,3),-F(1,24)),"quantum lower equation")
    raw=defect(Sh,Th,sh,th)
    shift=lin(Sh,Th,q4(S,s),q4(T,t))
    expect(add(raw,shift)==shear(O),"actual transformed order-five identity")
    # Genuine lower-lift direction; repeat coefficientwise at a second point.
    u,v=ks[0]
    ss,tt=add(s,u),add(t,v)
    raw2=defect(Sh,Th,q2(S,ss),q2(T,tt))
    shift2=lin(Sh,Th,q4(S,ss),q4(T,tt))
    expect(add(raw2,shift2)==shear(defect(S,T,ss,tt)),"second affine-lift quantum identity")
    # Gauge map d0 and d1 commute with tau, checked on a complete small
    # polynomial basis; the proof supplies arbitrary Hamiltonian coverage.
    gauges=laurent_monomials(4,4,0,1)
    for m in gauges:
        H={m+(0,):F(1)}
        u,v=pi(H,S),pi(H,T)
        expect(not lin(S,T,u,v),"gauge complex")
        expect(pi(shear(H),Sh)==shear(u) and pi(shear(H),Th)==shear(v),"gauge transport")
    # This specific native box misses explicit transported directions.
    witness={(18,0,1,0):F(1)}
    expect(shear(witness)=={(18,0,1,0):F(1),(22,0,0,1):F(3)},"support leakage witness")
    expect(native_failures==[{"component":0,"monomial":[12,6,1],"value":[[1,"-1296"]]}],"exact native counterexample")
    expect(pair(dual0,pi({(16,6,0,0):F(1)},T))=={0:F(432)},"independent inverse-support explanation")
    # Does the single new direction repair ALL order-five equations at the
    # selected lower lift? Vanishing of the scalar period alone is insufficient.
    original_current=[]
    for bound,zcap,which in ((21,1,0),(17,0,1)):
        for m in laurent_monomials(bound,zcap,0,3):
            pp={m+(0,):F(1)}
            original_current.append(pi(pp,T) if which==0 else pi(S,pp))
    original_current.append(pi({(16,6,0,0):F(1)},T))
    def toQQ(p): return {m:QQ(c.numerator,c.denominator) for m,c in p.items()}
    try:
        extra_base,extra_kernel,extra_rank=solve_affine(list(map(toQQ,original_current)),toQQ(scale(O,-1)),QQ)
        residual=add(O,*[scale(original_current[i],F(str(c))) for i,c in extra_base.items()])
        expect(not residual,"extended solution direct verification")
        one_extra={"status":"SOLVABLE","rank":extra_rank,"kernel_dimension":len(extra_kernel),
                   "extra_coefficient":str(extra_base.get(614,0)),"solution":[[i,str(c)] for i,c in sorted(extra_base.items())]}
    except ValueError:
        ranks=[column_rank(list(map(toQQ,original_current))),column_rank(list(map(toQQ,original_current))+[toQQ(O)])]
        expect(ranks[1]==ranks[0]+1,"exact enlarged rank jump")
        one_extra={"status":"INCONSISTENT","ranks":ranks,"scope":"all order-five equations at the selected bounded lower lift, with just X^16 Q^6 adjoined"}
    # Keep the finite-Laurent domain explicit even if this dual table happens
    # to be supported at nonnegative X exponents.
    source_names=["explore_degree_five_quantum_residue.py","explore_degree_five_a2_subprincipal.py","verify_degree_five_laurent_quantum_obstruction.py"]
    receipt={"status":"PASS","date":"2026-11-14","arithmetic":"exact rational polynomials in symbolic lambda",
             "seed":{"kappa":0,"tau":1,"a":"-1/2","S_terms":len(S),"T_terms":len(T)},
             "lower_family":{"columns":1569,"rank":1527,"dimension":42,"variation_coefficients_checked":nonlinear_checks},
             "current_correction_columns":counts,"native_support_leakage":len(leakage),
             "native_box_failures":native_failures,"unbounded_degree_relevance_checks":wider_count,
             "one_extra_current_direction":one_extra,
             "unbounded_degree_failures":wider_failures,
             "gauge_controls":len(gauges),"period":str(expected[0]),
             "transported_dual_terms":len(dual),"transported_dual":serial(dual),
             "actual_transformed_order5_identity":True,"two_lower_lift_points_checked":True,
             "raw_period_without_affine_shift":serial(pair(dual,raw)),
             "affine_shift_period":serial(pair(dual,shift)),
             "seed_h2_shift_S":serial(alpha2(S)),"seed_h2_shift_T":serial(alpha2(T)),
             "support_leakage_witness":serial(shear(witness)),
             "hashes":{n:hashlib.sha256((OLD/n).read_bytes()).hexdigest() for n in source_names},
             "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             "source_pair_sha256":[digest(S),digest(T)],
             "elapsed_seconds":round(time.monotonic()-start,4),
             "scope":"All-degree shear transport is a written automorphism proof. Replay covers symbolic cubic potential, all 614 current columns and 945 lower-family coefficients. Full finite-Laurent strength is inherited from the retained written source, not independently replayed here."}
    (HERE/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    # A small standalone input packet for independent arithmetic replay.
    packet={"S":serial(S),"T":serial(T),"s2":serial(s),"t2":serial(t),
            "lambda":serial(dual0),"transported_lambda":serial(dual)}
    (HERE/"inputs.json").write_text(json.dumps(packet,indent=2)+"\n")
    receipt["inputs_sha256"]=hashlib.sha256((HERE/"inputs.json").read_bytes()).hexdigest()
    (HERE/"receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print("PASS: symbolic shear, quantum affine shifts and actual order-five equations")
    print("PASS:",counts,"current columns;",nonlinear_checks,"lower-family coefficients;",len(gauges),"gauge controls")
    print("PERIOD",expected[0],"dual terms",len(dual),"native leakage",len(leakage))
    print("RAW_PERIOD",pair(dual,raw),"AFFINE_PERIOD",pair(dual,shift))
    print("NATIVE_FAILURES",len(native_failures),"WIDER_FAILURES",len(wider_failures),wider_failures[:2])
    print("ONE_EXTRA_DIRECTION",one_extra["status"])
    print("SECONDS",round(time.monotonic()-start,3))

if __name__=="__main__": main()
