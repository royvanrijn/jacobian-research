#!/usr/bin/env python3
"""Exact direct expansion controls for delayed nonpurity."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import comb,factorial,prod
from pathlib import Path

if not __debug__:
    raise SystemExit("Refusing -O: exact assertions must remain enabled.")


def raw_moment(a,N,epsilon,m):
    if a*m%2 or m%N:return Q(0)
    ell=m//N
    return sum((Q(comb(m,j)*comb(j,ell)*(-1)**(j-ell),a*m+2*j+1)
                *epsilon**ell for j in range(ell,m+1)),Q(0))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path(__file__).with_name(
        "two_height_delayed_nonpurity_verification.json"))
    args=parser.parse_args()
    controls=raw=0;rows=[]
    for a in range(9):
        for N in (2,4,6,8,10):
            for eps in (Q(1),Q(-1),Q(2,3)):
                values=[raw_moment(a,N,eps,m) for m in range(1,N+1)]
                assert all(v==0 for v in values[:-1])
                exact=eps*Q(2**(N-1)*factorial(N),
                            prod(a*N+3+2*j for j in range(N)))
                assert values[-1]==exact and exact!=0
                assert raw_moment(a,N,2*eps,N)==2*exact
                raw+=N+1;controls+=1
                rows.append({"a":a,"N":N,"epsilon":str(eps),
                             "phase_degree":N+a+1,"first_bad_moment":N,
                             "moment_value":str(exact)})
    proof=Path(__file__).with_name("ALL_START_TWO_HEIGHT_EFFECTIVE_PURITY_PREFIX.md")
    out={"status":"PASS","scope":"Finite exact direct-expansion controls; the unrestricted first-bad-order statement is the written angular-weight proof.",
         "proof_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
         "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         "fixtures":controls,"direct_moment_controls":raw,"results":rows,
         "pure_counterexample_claimed":False,"external_or_formal_review_claimed":False}
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in out.items() if k!="results"},sort_keys=True))


if __name__=="__main__":main()
