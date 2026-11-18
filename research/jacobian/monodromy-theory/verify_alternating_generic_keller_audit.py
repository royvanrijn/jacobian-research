#!/usr/bin/env python3
"""Exact small-prime controls for the uniform primitive-monomial counterexample.

The all-exponent/all-prime theorem is proved in the accompanying note. These
checks verify the substitutions, determinant, integral presentation, and
bounded monomial remainders, and reject a constant-monomial control.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from pathlib import Path
import sympy as sp


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt",type=Path)
    args=parser.parse_args()
    assert sp.__version__=="1.14.0"
    x,y,u,v=sp.symbols("x y u v")
    rows=[]
    for ell in (3,5,7,11):
        U=x+y**(ell-1)
        V=x*y+2*y**ell
        assert sp.expand(V-(y**ell+U*y))==0
        jac=sp.det(sp.Matrix([[sp.diff(f,z) for z in (x,y)] for f in (U,V)]))
        assert sp.expand(jac-x-(ell+1)*y**(ell-1))==0
        assert not jac.is_number
        relation=sp.Poly(y**ell+u*y-v,y,domain=sp.QQ.frac_field(u,v))
        assert relation.is_irreducible
        cases=0
        for i in range(13):
            for j in range(13):
                remainder=sp.rem(sp.Poly((u-y**(ell-1))**i*y**j,y,
                                         domain=sp.QQ.frac_field(u,v)),relation)
                if i==j==0:
                    assert remainder.degree()==0
                else:
                    assert remainder.degree()>0,(ell,i,j)
                    cases+=1
        plus=sp.rem(sp.Poly(u-y**(ell-1)+y,y,domain=sp.QQ.frac_field(u,v)),relation)
        assert plus.degree()==ell-1
        # Modular weights force u exponents to be multiples of ell.
        assert all(a%ell==0 for a in range(4*ell+1) for b in range(4*ell+1)
                   if ((ell-1)*a+ell*b)%ell==0)
        rows.append({"prime":ell,"monomial_controls":cases,
                     "xy_sum_remainder_degree":plus.degree(),
                     "jacobian":str(jac),"extension_degree":relation.degree()})
    report={
        "schema_version":1,"kind":"internal_exact_replay","result":"PASS",
        "command":".venv/bin/python research/jacobian/monodromy-theory/verify_alternating_generic_keller_audit.py",
        "python":platform.python_version(),"sympy":sp.__version__,
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "prime_controls":rows,"monomial_exponents_tested":"0<=i,j<=12, excluding (0,0)",
        "scope":"Finite controls supplement the written uniform weighted proof. The counterexample refutes the stated field-theoretic auxiliary criterion, not a Keller conclusion. These maps have nonconstant determinant.",
        "independent_replay":False,"formal_verification":False,"external_review":False,
    }
    output=json.dumps(report,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:args.write_receipt.write_text(output)
    print(output,end="")


if __name__=="__main__":main()
