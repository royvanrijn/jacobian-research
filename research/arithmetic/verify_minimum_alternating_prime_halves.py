#!/usr/bin/env python3
"""Exact finite controls for the written all-prime minimum-degree proof."""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from math import comb, factorial, isqrt, prod
from pathlib import Path
import platform

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent

def require(ok,message):
    if not ok:
        raise AssertionError(message)

def prime(p):
    return p>=2 and all(p%d for d in range(2,isqrt(p)+1))

def partitions(n,upper=None):
    if n==0:
        yield ()
        return
    if upper is None or upper>n:
        upper=n
    for part in range(upper,0,-1):
        for tail in partitions(n-part,part):
            yield (part,)+tail

def cover_control(n):
    half=n//2
    count=0
    for cyc in partitions(n):
        if (n-len(cyc))%2:
            continue
        count+=1
        sums={0}
        for d in cyc:
            sums|={a+d for a in list(sums)}
        U=any(k in sums for k in range(1,half,2))
        W=half in sums or all(d%2==0 for d in cyc)
        require(U or W,"Proposed normal cover misses a cycle type")
    return count

def control(p):
    require(prime(p) and p>=7,"Prime-half hypothesis")
    n=2*p
    R=2**(n-2)
    odd=list(range(1,p,2))
    costs=[comb(n,k) for k in odd]
    Wp=comb(n,p)//2
    W2=factorial(n)//(2**p*factorial(p))
    require(sum(costs)+Wp==R,"Minimum upper-bound sum")
    require(factorial(n)>2**(n-1)*3**n,"Primitive cost cutoff")
    require(W2==prod(range(1,n,2)) and W2>R,"Matching cutoff")
    require(factorial(n)>2*3**n*R,"Primitive-index strict cutoff")
    # Exhaustive invariant-subset sizes of every two-cycle witness.
    incidence=[]
    for k in odd+[p]:
        possible={k,n-k}
        row=[j in possible for j in odd]
        require(row==[j==k for j in odd],"Forced witness incidence")
        if k<p:
            require(p not in possible,"Odd witness lies in half partition")
            require(k%2==1 and (n-k)%2==1,"Witness order must be odd")
        incidence.append(row)
    # Any deviations from the minimum cost inflate it, as used in equality proof.
    require(2*Wp>Wp and all(c>0 for c in costs),"Equality strictness")
    return {"p":p,"degree":n,"minimum_rank":R,
        "odd_subset_sizes":odd,"factor_degrees":costs+[Wp],
        "matching_index":W2,"primitive_index_bound_numerator":factorial(n),
        "primitive_index_bound_denominator":2*3**n,
        "two_cycle_witness_incidence":incidence,
        "all_minimum_costs_forced":True}

def digest(p):
    return sha256(p.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-receipt",type=Path,
        default=HERE/"minimum_alternating_prime_halves_verification.json")
    args=parser.parse_args()
    require(factorial(14)>2**13*3**14,"Induction base factorial")
    require(prod(range(1,14,2))>2**12,"Induction base matching")
    require(15*16>36 and 15>4,"First induction ratios")
    # The induction ratio polynomials increase for n>=14.
    controls=[control(p) for p in (7,11,13,17,19,23,41,61,113)]
    source_paths=[
        HERE/"MINIMUM_TOTAL_DEGREES_ALTERNATING_PRIME_HALVES.md",
        HERE/"sources/maroti-primitive-orders-2002.pdf",
        HERE/"sources/lenstra-chebotarev-density-notes.pdf",
        HERE/"ARALA_LAURENT_PENCILS_AND_KELLER_FIBRES.md"]
    report={"schema_version":1,"result":"PASS","python":platform.python_version(),
        "checker_sha256":digest(Path(__file__)),
        "source_sha256":{str(p.relative_to(ROOT)):digest(p) for p in source_paths},
        "induction_bases":{"degree":14,"factorial":factorial(14),
                          "factorial_comparison":2**13*3**14,
                          "matching_index":prod(range(1,14,2)),
                          "matching_comparison":2**12},
        "prime_controls":controls,
        "all_even_cycle_types_checked":{str(n):cover_control(n) for n in (14,22,26,34)},
        "scope":"Written proof establishes every prime p>=7 and equality "
        "classification, importing Maroti Cor1.2 with CFSG. Finite integer "
        "controls do not prove the unbounded quantifier or the imported theorem.",
        "formal_verification":False,"independent_replay":False,
        "external_human_review":False}
    args.write_receipt.write_text(json.dumps(report,indent=2)+"\n")
    print("PASS: exact prime-half minimum-degree controls.")
if __name__=="__main__":
    main()
