#!/usr/bin/env python3
"""Exact controls for the uniform BGV critical-genus and rational-correction proof.

Finite controls do not prove the all-prime argument or surface criterion.
Sparse arithmetic is reused from the existing BGV lifting checker.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path
import platform
import re
import shutil
import subprocess
import verify_bgv_odd_prime_lifting as w

HERE=Path(__file__).resolve().parent
X,Y,ONE=w.X,w.Y,w.ONE

def shift(a,i,j):
    return w.clean({(m+i,n+j):c for (m,n),c in a.items()})

def substitute(a,xx,yy):
    out={}
    for (i,j),c in a.items():
        assert i>=0 and j>=0
        out=w.add(out,w.scale(w.mul(w.power(xx,i),w.power(yy,j)),c))
    return out

def univariate_gcd(a,b):
    def remainder(f,g):
        f=dict(f)
        dg=max(j for (i,j) in g)
        cg=g[(0,dg)]
        while f:
            df=max(j for (i,j) in f)
            if df<dg:break
            cf=f[(0,df)]
            f=w.sub(f,w.scale(shift(g,0,df-dg),Fraction(cf,cg)))
        return f
    while b:a,b=b,remainder(a,b)
    return a

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write-receipt",type=Path)
    args=ap.parse_args()
    rows=[]
    for p in (3,5,7):
        z,P,Q=w.family(p)
        J=w.jac(P,Q)
        # Compute J(x,(z-x)/x^p) independently from the expanded source.
        Jxz={}
        for (i,j),c in J.items():
            Jxz=w.add(Jxz,w.scale(shift(w.power(w.sub(Y,X),j),i-p*j,0),c))
        lhs=shift(Jxz,1,0)
        compact=w.add(w.scale(X,1-p),w.scale(Y,p))
        for i,j,c in [(p,p,p*(p+1)),(p+1,p-1,-p*(p-1)),
                      (p*p,p*p-p+1,-p*(p-1)),(p*p+p,p*p,p)]:
            compact=w.add(compact,{(i,j):c})
        assert lhs==compact
        # x=u/z. After multiplying xJ by z^p, get u*z^(p-1)*J.
        K=w.clean({(i,j-i+p):c for (i,j),c in compact.items()})
        expected={}
        for i,j,c in [(1,p-1,1-p),(0,p+1,p),(p,p,p*(p+1)),
                      (p+1,p-2,-p*(p-1)),(p*p,1,-p*(p-1)),
                      (p*p+p,0,p)]:
            expected=w.add(expected,{(i,j):c})
        assert K==expected
        H=w.add(w.add(w.power(Y,p+1),w.scale(w.power(Y,p),p+1)),
                w.add(w.scale(Y,-(p-1)),ONE))
        T=w.add(w.add(w.scale(w.power(Y,2),p-1),
                     w.scale(Y,(p+1)*(p-2))),w.scale(ONE,-(p+1)))
        N=2*p-1
        # z=u^p*v gives the exact Kummer equation.
        Knew=w.clean({(i+p*j,j):c for (i,j),c in K.items()})
        numerator=w.scale(w.mul(w.power(Y,p-2),w.add(Y,w.scale(ONE,p))),p-1)
        kummer=w.sub(w.scale(shift(H,N,0),p),numerator)
        assert Knew==shift(kummer,p*p-p+1,0)
        bad=w.add(kummer,ONE)
        assert Knew!=shift(bad,p*p-p+1,0)
        # Exact repeated-root reduction.
        identity=w.sub(w.scale(w.mul(w.add(Y,w.scale(ONE,p)),H),p+1),
                       w.mul(w.mul(Y,w.add(Y,w.scale(ONE,p+1))),w.diff(H,1)))
        assert identity==w.scale(T,-p)
        shiftedT=substitute(T,X,w.sub(Y,ONE))
        assert shiftedT==w.add(w.add(w.scale(w.power(Y,2),p-1),
                                    w.scale(Y,p*(p-3))),w.scale(ONE,-p*(p-1)))
        shiftedH=substitute(H,X,w.sub(Y,ONE))
        assert w.mod(shiftedH,p)==w.power(Y,p+1)
        assert shiftedH.get((0,0),0)==0
        assert shiftedH[(0,1)]==p*(p-1)
        assert max(j for (i,j) in univariate_gcd(H,w.diff(H,1)))==0
        evaluate=lambda a,n:sum(c*n**j for (i,j),c in a.items())
        assert evaluate(H,0)==1
        assert evaluate(H,-p)==p*p-p+1-p**p!=0
        d0=gcd(N,p-2)
        different=(N-d0)+(p+3)*(N-1)
        genus=(different-2*N+2)//2
        assert 2*genus-2==different-2*N
        assert genus==p*p+p-(3+d0)//2>0
        rows.append(dict(prime=p,source_jacobian_terms=len(J),
                         birational_critical_identities=True,
                         corrupted_identity_rejected=True,
                         squarefree_controls=True,kummer_degree=N,
                         zero_branch_gcd=d0,critical_genus=genus,
                         noncontraction_order=p,ramification_index=2))
    singular=shutil.which("Singular")
    assert singular, "Singular is required for the independent p=3 genus control"
    singular_input='''LIB "normal.lib";
ring R=0,(u,v),dp;
poly H=v^4+4*v^3-2*v+1;
ideal I=3*u^5*H-2*v*(v+3);
print("BGV_GENUS_BEGIN");
genus(I);
print("BGV_GENUS_END");
ideal I0=v^2-u^3;
print("RATIONAL_GENUS_BEGIN");
genus(I0);
print("RATIONAL_GENUS_END");
quit;
'''
    proc=subprocess.run([singular,"-q"],input=singular_input,text=True,
                        capture_output=True,timeout=120,check=True)
    assert re.search(r"BGV_GENUS_BEGIN\s+10\s+BGV_GENUS_END",proc.stdout)
    assert re.search(r"RATIONAL_GENUS_BEGIN\s+0\s+RATIONAL_GENUS_END",proc.stdout)
    assert "error" not in proc.stdout.lower() and "error" not in proc.stderr.lower()
    version=subprocess.check_output([singular,"--version"],text=True).splitlines()[0]
    lib=Path("/usr/share/singular/LIB/normal.lib")
    assert lib.is_file()
    # Rational genus-zero control: (x^2,y/(2x)) has Jacobian one.
    rationalP=w.power(X,2)
    rationalQ={(-1,1):Fraction(1,2)}
    assert w.jac(rationalP,rationalQ)==ONE
    paths=["BGV_CRITICAL_GENUS_AND_RATIONAL_CORRECTIONS.md",
           "verify_bgv_odd_prime_lifting.py",
           "bgv_critical_genus_source_ledger.json",
           "sources/brosnan-surfaces-elimination-indeterminacy.html",
           "../../hessian/hc4-triangular-graphs/sources/stacks-0C1B-riemann-hurwitz.html"]
    report=dict(schema_version=1,result="PASS",kind="internal_exact_replay",
        python=platform.python_version(),
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256={q:hashlib.sha256((HERE/q).read_bytes()).hexdigest() for q in paths},
        controls=rows,
        independent_genus_control=dict(command="Singular -q",version=version,
            normal_lib_sha256=hashlib.sha256(lib.read_bytes()).hexdigest(),
            input=singular_input,stdout=proc.stdout,stderr=proc.stderr,
            result="genus 10 for BGV p=3 and genus 0 for a rational cusp"),
        scope="Finite exact identities at p=3,5,7 and a separate Singular genus calculation. The written proof establishes the uniform critical genus and the general rational pre/post-correction obstruction. No independent whole-proof replay, formal verification, external review or JC2 solution.",
        formal_verification=False,external_review=False,independent_implementation=False)
    out=json.dumps(report,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:args.write_receipt.write_text(out)
    print(out,end="")

if __name__=="__main__":main()
