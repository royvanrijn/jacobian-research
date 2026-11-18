#!/usr/bin/env python3
"""Exact controls for the known weighted optimum and uniform equality proof."""
from __future__ import annotations
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb,factorial
from pathlib import Path
import platform

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent

def require(ok,message):
    if not ok:
        raise AssertionError(message)

def index(a,r):
    return factorial(a*r)//(factorial(a)**r*factorial(r))

def partitions(n,upper=None):
    if n==0:
        yield ()
        return
    upper=min(n,upper or n)
    for p in range(upper,0,-1):
        for tail in partitions(n-p,p):
            yield (p,)+tail

def cover(n,exhaust=False):
    require(n>=6 and n%4==2,"Covering theorem domain")
    b=n//2
    costs=[comb(n,k) for k in range(1,b,2)]+[comb(n,b)//2]
    require(sum(costs)==2**(n-2),"Odd-subset/half-partition rank")
    types=0
    if exhaust:
        for cyc in partitions(n):
            if (n-len(cyc))%2:
                continue
            types+=1
            subs={0}
            for d in cyc:
                subs|={v+d for v in list(subs)}
            require(any(k in subs for k in range(1,b,2)) or b in subs
                    or all(d%2==0 for d in cyc),"Uncovered even cycle type")
    return {"degree":n,"rank":sum(costs),"factor_degrees":costs,
            "even_cycle_types_checked":types if exhaust else None}

def digest(p):
    return sha256(p.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-receipt",type=Path,
        default=HERE/"alternating_weighted_cover_priority_verification.json")
    args=parser.parse_args()
    require(index(3,3)==280>2**7,"a3r3 induction base")
    require(index(2,4)==105>2**6,"a2r4 induction base")
    require(index(2,3)==15<2**4,"Necessary n6 exception")
    imprimitive=[]
    for a in range(2,51):
        for r in range(3,51):
            if a*r<=200 and a*r>=8:
                I=index(a,r)
                require(I>2**(a*r-2),"Uniform imprimitive finite control")
                imprimitive.append((a,r,I))
            if r<=20 and a<=20:
                ratio=Fraction(index(a,r+1),index(a,r))
                require(ratio==Fraction(comb((r+1)*a,a),r+1),
                        "r-step exact ratio")
                require(comb((r+1)*a,a)>=(r+1)**a,
                        "Labeled-group binomial injection bound")
                require(ratio>=2**a,"r-step growth")
        require(Fraction(index(a+1,3),index(a,3))==
                Fraction(3*(3*a+1)*(3*a+2),(a+1)**2),"a-step ratio")
        require(19*a*a+11*a-2>0,"a-step strict growth")
    primitive=[]
    for n in range(14,202,2):
        require(factorial(n)>2**(n-1)*3**n,"Primitive cutoff finite control")
        require((n+1)*(n+2)>36,"Primitive cutoff induction growth")
        primitive.append(n)
    controls=[cover(n,True) for n in (6,10,14,18,22,26,30)]
    admissible=[cover(m*m+1,False) for m in (3,5,7,9,11)]
    paths=[
        HERE/"ALTERNATING_WEIGHTED_COVERS_PRIORITY_AND_OPTIMAL_HASSE.md",
        HERE/"sources/maroti-covering-symmetric-2004.pdf",
        HERE/"sources/maroti-primitive-orders-2002.pdf",
        HERE/"sources/bubboloni-praeger-normal-1011.4368v1.pdf",
        HERE/"sources/bubboloni-sonn-intersective-1507.08593v1.pdf",
        HERE/"sources/sonn-intersective-two-remarks-2009.pdf",
        HERE/"sources/lenstra-chebotarev-density-notes.pdf",
        HERE/"ARALA_LAURENT_PENCILS_AND_KELLER_FIBRES.md"]
    result={"schema_version":1,"result":"PASS","python":platform.python_version(),
        "checker_sha256":digest(Path(__file__)),
        "source_sha256":{str(p.relative_to(ROOT)):digest(p) for p in paths},
        "imported_covering_theorem":{
          "source":"Maroti, Covering the symmetric groups with proper subgroups",
          "locator":"Theorem4.1, printedp9; explicitcover/proofpp10--11",
          "domain":"n>=6,n=2mod4","ordinary_covering_minimum":"2^(n−2)",
          "checker_reproves_import":False},
        "imported_primitive_order_bound":{
          "source":"Maroti, On the Orders of Primitive Groups",
          "locator":"Corollary1.2, printedp3",
          "statement":"Primitive H<=Sn not containingAn has |H|<3^n.",
          "classification_dependent":True,"checker_reproves_import":False},
        "imprimitive_control_count":len(imprimitive),
        "imprimitive_control_domain":"2<=a<=50,3<=r<=50,8<=ar<=200",
        "imprimitive_bases":{"I3_3":280,"I2_4":105,"excluded_I2_3":15},
        "primitive_cutoff_even_degrees":primitive,
        "small_degree_upper_cover_controls":controls,
        "admissible_degree_rank_controls":admissible,
        "uniform_equality_domain":"n>=14,n=2mod4",
        "numerical_optimum_priority":"Known-theorem consequence of Maroti2005.",
        "equality_priority":"Not resolved; no novelty claim.",
        "scope":"Finite controls supplement written uniform index/equality proof. "
        "All-degree numerical minimum imports MarotiTheorem4.1. The Arala Hasse "
        "family consequently has globally minimum rank in every admissible degree.",
        "formal_verification":False,"independent_external_replay":False,
        "external_human_review":False}
    args.write_receipt.write_text(json.dumps(result,indent=2)+"\n")
    print("PASS: known weighted-optimum and uniform equality arithmetic controls.")
if __name__=="__main__":
    main()

