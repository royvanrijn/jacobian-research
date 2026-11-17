#!/usr/bin/env python3
"""Analyze every GL2 transverse completion at the 49-term collision pair.

Diagonal column scalings and a column swap preserve support, so the two
unipotent parameters cover all invertible transverse bases.
"""
import hashlib,json
from pathlib import Path
import sympy as sp
import two_parameter as core
import specialize as fibre
HERE=Path(__file__).resolve().parent
OUT=HERE/"full-transverse"
x,y,z,w,t=sp.symbols("x y z w t")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    OUT.mkdir(exist_ok=True)
    source=HERE.parent/"sparse-nine/normalized-49-witness.json"
    data=json.loads(source.read_text())
    H=sp.Matrix([sum(sp.Rational(v["coefficient"])*sp.prod(
        a**b for a,b in zip((x,y,z),v["exponent"])) for v in row) for row in data["coefficients"]])
    K=H.subs(z,z+w*y);K[2]-=w*K[1]
    K=K.subs(y,y+t*z);K[1]-=t*K[2]
    coefficients=[];factors={}
    for i,entry in enumerate(K):
        for exponent,c in sp.Poly(sp.expand(entry),x,y,z).terms():
            if sum(exponent)<2:continue
            polynomial=sp.Poly(c,w,t,domain=sp.QQ)
            coefficients.append(dict(row=i,exponent=list(exponent),coefficient=str(polynomial.as_expr())))
            for f,m in sp.factor_list(polynomial)[1]:
                f=f.monic();factors[str(f.as_expr())]=f
    family=dict(
        schema="full-transverse-frame-family.v1",source_sha256=sha(Path(__file__)),
        predecessor_sha256=sha(source),generic_support=len(coefficients),coefficients=coefficients,
        matrix=[["1","0","0"],["0","1","t"],["0","w","1+w*t"]],
        factors=[dict(polynomial=key,degree_w=f.degree(w),degree_t=f.degree(t),
                      vanishing_coefficients=[i for i,row in enumerate(coefficients)
                        if sp.Poly(row["coefficient"],w,t,domain=sp.QQ).rem(f).is_zero])
                 for key,f in sorted(factors.items())])
    core.save(OUT/"two-parameter-family.json",family)
    core.HERE=OUT
    core.project(family)
    fibre.HERE=OUT
    fibre.main()
if __name__=="__main__":main()

