#!/usr/bin/env python3
"""Finite exact ordering diagnostic; no correction-support exclusion."""
from functools import lru_cache
from itertools import combinations
from math import factorial
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import time
import verify_long_pbw as a

def compositions(n):
    for i in range(n+1):
        for j in range(n-i+1):
            for k in range(n-i-j+1):
                yield i,j,k,n-i-j-k

def bidifferential(left, right, n):
    # Pi = d_p tensor d_x - d_x tensor d_p
    #    + d_z tensor d_q - d_q tensor d_z.
    @lru_cache(None)
    def derivative(side, orders):
        p = left if side == 0 else right
        for i, count in enumerate(orders):
            for _ in range(count):
                p = a.diff(p,i)
        return p
    out = {}
    for i,j,k,l in compositions(n):
        dl = derivative(0,(j,l,i,k))
        dr = derivative(1,(i,k,j,l))
        coeff = F(factorial(n), factorial(i)*factorial(j)*factorial(k)*factorial(l))
        out = a.add(out,a.scale(a.mul(dl,dr),coeff*(-1)**(j+l)))
    return out

def main():
    start=time.monotonic()
    a.require(a.scale(bidifferential(a.power(a.var(2),3),a.power(a.var(0),3),3),F(1,24))
              == a.const(F(3,2)),"known p^3,x^3 Moyal normalization")
    ims=a.long_images()
    names=["R","T","D","S"]
    rows=[]
    for i,j in combinations(range(4),2):
        layers=[]
        # A contraction consumes one momentum derivative on one of the
        # two factors, so total momentum degree bounds ALL higher layers.
        cap=max(e[2]+e[3] for e in ims[i])+max(e[2]+e[3] for e in ims[j])
        for n in range(3,cap+1,2):
            p=a.scale(bidifferential(ims[i],ims[j],n),F(1,2**(n-1)*factorial(n)))
            if (i,j,n)==(0,1,3):
                x,q=a.var(0),a.var(1)
                compact=a.scale(a.mul(a.power(x,8),a.mul(
                    a.add(a.scale(a.mul(q,x),3),a.const(-2)),
                    a.add(a.scale(a.mul(q,x),3),a.const(-1)))), -F(3,2))
                a.require(p==compact,"compact symmetric RT defect")
            if p:
                first=next(iter(sorted(p.items())))
                layers.append({"h_order":n,"terms":len(p),"sha256":a.digest(p),
                               "witness":[list(first[0]),str(first[1])]})
        rows.append({"pair":[names[i],names[j]],"momentum_bound":cap,"layers":layers})
    a.require(not rows[2]["layers"],"Weyl ordered R,S commute exactly")
    a.require(any(row["layers"] for row in rows),"some symmetric-order defects survive")
    out={"date":"2026-11-14","status":"PASS","method":"complete Moyal commutator; Pi^n/(2^(n-1)*n!) for odd n",
         "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         "source_sha256":hashlib.sha256(Path(a.__file__).read_bytes()).hexdigest(),
         "compact_RT_h3":"-3*x^8*(3*q*x-2)*(3*q*x-1)/2",
         "relations":rows,"elapsed_seconds":round(time.monotonic()-start,4),
         "boundary":"Only the literal symmetric Weyl quantizations are tested; corrections and other images remain open."}
    Path(__file__).with_name("weyl_ordering_receipt.json").write_text(json.dumps(out,indent=2)+"\n")
    for r in rows:
        print(r["pair"],[(x["h_order"],x["terms"]) for x in r["layers"]])

if __name__=="__main__":
    main()
