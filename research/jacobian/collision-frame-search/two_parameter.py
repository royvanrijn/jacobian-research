#!/usr/bin/env python3
"""Exact lower-triangular transverse completions of every boundary collision frame."""
from pathlib import Path
from fractions import Fraction
import argparse,hashlib,json,time
import sympy as sp
HERE=Path(__file__).resolve().parent
PREV=HERE.parent/"sparse-nine"
w,t=sp.symbols("w t")
x,y,z=sp.symbols("x y z")
def save(path,payload):
    temporary=path.with_suffix(path.suffix+".research/.scratch")
    temporary.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    temporary.replace(path)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def compile_family():
    old=json.loads((PREV/"frame-family-certificate.json").read_text())
    H=sp.Matrix([x,y,z])
    for row in old["coefficient_rows"]:
        H[row["row"]]+=sp.sympify(row["coefficient"])*sp.prod(
            variable**power for variable,power in zip((x,y,z),row["exponent"]))
    K=H.subs(z,z+t*y);K[2]-=t*K[1]
    data=[];factors={}
    for i,entry in enumerate(K):
        for exponent,c in sp.Poly(sp.expand(entry),x,y,z).terms():
            if sum(exponent)<2:continue
            polynomial=sp.Poly(c,w,t,domain=sp.QQ)
            data.append(dict(row=i,exponent=list(exponent),coefficient=str(polynomial.as_expr())))
            for f,multiplicity in sp.factor_list(polynomial)[1]:
                f=f.monic();factors[str(f.as_expr())]=f
    result=dict(schema="two-parameter-frame-family.v1",source_sha256=sha(Path(__file__)),
                predecessor_sha256=sha(PREV/"frame-family-certificate.json"),
                coefficients=data,generic_support=len(data),
                factors=[dict(polynomial=key,degree_w=f.degree(w),degree_t=f.degree(t),
                              vanishing_coefficients=[i for i,row in enumerate(data)
                                  if sp.Poly(row["coefficient"],w,t).rem(f).is_zero])
                         for key,f in sorted(factors.items())])
    save(HERE/"two-parameter-family.json",result)
    return result

def project(family):
    path=HERE/"projection-checkpoint.json"
    polynomials=[sp.Poly(row["polynomial"],t,w,domain=sp.QQ) for row in family["factors"]]
    vertical=[f for f in polynomials if f.degree(t)==0]
    moving=[f for f in polynomials if f.degree(t)>0]
    sources=[]
    for f in vertical:sources.append(("vertical",str(f.as_expr())))
    for f in moving:
        leading=sp.Poly(f.as_expr(),t).LC()
        sources.append(("leading",str(leading)))
    sources.extend(("pair",str(f.as_expr()),str(g.as_expr()))
                   for i,f in enumerate(moving) for g in moving[:i])
    if path.exists():
        state=json.loads(path.read_text())
        assert state["source_sha256"]==sha(Path(__file__))
        assert state["family_sha256"]==sha(HERE/"two-parameter-family.json")
    else:state=dict(source_sha256=sha(Path(__file__)),family_sha256=sha(HERE/"two-parameter-family.json"),
                    done=0,projection_factors={},resultants=[])
    for index in range(state["done"],len(sources)):
        source=sources[index]
        if source[0]=="pair":
            expression=sp.resultant(sp.sympify(source[1]),sp.sympify(source[2]),t)
            assert expression!=0
        else:expression=sp.sympify(source[1])
        polynomial=sp.Poly(expression,w,domain=sp.QQ)
        if polynomial.degree()>0:
            for f,multiplicity in sp.factor_list(polynomial)[1]:
                f=f.monic()
                state["projection_factors"][str(f.as_expr())]=f.degree()
        state["resultants"].append(dict(source=list(source),polynomial=str(polynomial.as_expr())))
        state["done"]=index+1
        if (index+1)%25==0:
            save(path,state)
            print(json.dumps(dict(stage="projection",done=index+1,total=len(sources),
                                  factors=len(state["projection_factors"]))),flush=True)
    state["status"]="complete";save(path,state)
    save(HERE/"projection-factors.json",dict(
        source_sha256=sha(Path(__file__)),family_sha256=sha(HERE/"two-parameter-family.json"),
        status="complete",factors=[dict(polynomial=f,degree=d)
            for f,d in sorted(state["projection_factors"].items(),key=lambda item:(item[1],item[0]))],
        counts_by_degree={str(d):list(state["projection_factors"].values()).count(d)
                          for d in sorted(set(state["projection_factors"].values()))}))
    print(json.dumps(dict(stage="projection",status="complete",sources=len(sources),
                          factors=len(state["projection_factors"]))),flush=True)

def main():
    family=compile_family()
    print(json.dumps(dict(stage="family",coefficients=len(family["coefficients"]),
                          factors=len(family["factors"]))),flush=True)
    project(family)
if __name__=="__main__":main()

