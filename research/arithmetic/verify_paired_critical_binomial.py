#!/usr/bin/env python3
"""Exact controls for the normalized binomial-critical obstruction.

The all-b proof is written in PAIRED_CRITICAL_BINOMIAL_OBSTRUCTION.md.
These exact degree10/26 controls verify the full discriminant identities
and squarefree nonzero squared-critical-value polynomials.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sympy as sp

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt",type=Path)
    args=parser.parse_args()
    assert sp.__version__=="1.14.0"
    u,z,s,y=sp.symbols("u z s y")
    records=[]
    for b in (5,13):
        n=2*b;d=b-1
        k=-sp.Rational(1,(n-1)*(n-3))
        rho=sp.Rational(n-1,n-3)
        N=u**b-rho*u**d+(u+1)/(n-3)
        S=(n-1)*(u**d-k)
        assert sp.expand(2*u*sp.diff(N,u)-N-(u-1)*S)==0
        assert N.subs(u,1)==0 and sp.diff(N,u).subs(u,1)==0
        assert S.subs(u,0)!=0 and S.subs(u,1)!=0
        assert sp.gcd(S,sp.diff(S,u))==1
        reduced=sp.rem(N,u**d-k,u)
        alpha=sp.Poly(reduced,u).nth(1)
        beta=reduced.subs(u,0)
        assert beta/alpha==rho
        collision=sp.factor((alpha**2*u+2*alpha*beta+beta**2/u)
                          -(alpha**2*y+2*alpha*beta+beta**2/y))
        assert sp.cancel(collision-(u-y)*(alpha**2-beta**2/(u*y)))==0
        assert k*k<1<rho**(2*d)
        assert abs((-rho)**d)>1>abs(k)
        R=u**d-k
        chi=sp.cancel(sp.resultant(R,u*y-reduced**2,u)/sp.resultant(R,u,u))
        assert sp.Poly(chi,y).LC()==1 and sp.Poly(chi,y).degree()==d
        assert chi.subs(y,0)!=0 and sp.gcd(chi,sp.diff(chi,y))==1
        H=sp.expand(N.subs(u,z*z)-N.subs(u,0))
        a=N.subs(u,0)
        discriminant=sp.discriminant(H-s*z+a,z)
        expected=(-1)**(b+1)*(n-1)**(n-1)*s**2*chi.subs(y,s*s)
        assert sp.expand(discriminant-expected)==0
        # Removing the zero branch leaves a squarefree polynomial of degree n-2.
        residual=sp.cancel(discriminant/s**2)
        assert sp.Poly(residual,s).degree()==n-2
        assert sp.gcd(residual,sp.diff(residual,s))==1
        assert sp.gcd(discriminant,sp.diff(discriminant,s))==s
        assert sp.expand(discriminant-expected+s)!=0
        records.append({
            "b":b,"degree":n,"H":str(H),"a":str(a),"k":str(k),"rho":str(rho),
            "squared_critical_value_polynomial_degree":d,
            "squared_critical_value_polynomial":str(chi),
            "nonzero_squared_values":True,"all_squared_values_distinct":True,
            "full_discriminant_identity":True,
            "finite_branch_passport":"one double transposition at zero; n-2 separate transpositions",
            "monodromy_with_written_primitive_transposition_argument":f"S{n}",
        })
    report={
        "schema_version":1,"kind":"internal_exact_symbolic_controls",
        "python":platform.python_version(),"sympy":sp.__version__,
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "cases":records,
        "scope":"Finite b=5,13 controls. The written proof excludes the normalized binomial S ansatz for every b>=3 over C; no general degree26 impossibility follows.",
        "independent_replay":False,"formal_verification":False,
        "external_review":False,"result":"PASS",
    }
    output=json.dumps(report,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:args.write_receipt.write_text(output)
    print(output,end="")
if __name__=="__main__":main()

