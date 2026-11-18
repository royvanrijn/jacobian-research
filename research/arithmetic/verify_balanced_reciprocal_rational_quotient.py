#!/usr/bin/env python3
"""Uniform balanced sign-quotient and degree26 finite rational search system."""
from __future__ import annotations
import argparse,hashlib,importlib.util,itertools,json,platform
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
SOURCE=HERE/"verify_reciprocal_parameter_finiteness.py"
spec=importlib.util.spec_from_file_location("reciprocal_raw_source",SOURCE)
finite=importlib.util.module_from_spec(spec);spec.loader.exec_module(finite)
w=finite.w;r=finite.r
eta=sp.symbols("eta")

def quotient_case(k):
    data=finite.raw_system(2*k,k)
    pc=data["variables"][:k];qc=data["variables"][k:2*k]
    x=sp.symbols("x0:"+str(k));z=sp.symbols("z0:"+str(k))
    variables=list(x)+list(z)+[eta]
    signs=[(-1)**(k+i)for i in range(k)]
    substitution={}
    tau={r:-r}
    for i in range(k):
        substitution[pc[i]]=(x[i]+r*z[i])/2
        substitution[qc[i]]=signs[i]*(x[i]-r*z[i])/2
        tau[pc[i]]=signs[i]*qc[i]
        tau[qc[i]]=signs[i]*pc[i]
    def parts(expr):
        polynomial=sp.Poly(sp.expand(expr.subs(substitution,simultaneous=True)),r)
        E=O=sp.Integer(0)
        for(j,),co in polynomial.terms():
            if j%2:O+=co*eta**((j-1)//2)
            else:E+=co*eta**(j//2)
        return sp.expand(E),sp.expand(O)
    # Exact involution on the original rational-polynomial reconstruction.
    for j,co in enumerate(data["F"]):
        assert sp.expand(co.subs(tau,simultaneous=True)-(-1)**(j+1)*co)==0
    assert sp.expand(data["X"].subs(tau,simultaneous=True)-data["Y"].subs(w,-w))==0
    assert sp.expand(data["Y"].subs(tau,simultaneous=True)-data["X"].subs(w,-w))==0
    assert sp.expand(data["G"].subs(tau,simultaneous=True)-data["G"])==0

    raw=[];strong=[]
    for target,power in [(raw,1),(strong,2)]:
        plus=sp.Poly(sp.rem(data["X"],data["P"]**power,w),w)
        minus=sp.Poly(sp.rem(data["Y"],data["Q"]**power,w),w)
        for j in range(power*k):
            E,O=parts(plus.nth(j))
            ME,MO=parts(minus.nth(j))
            assert sp.expand(ME-(-1)**j*E)==0
            assert sp.expand(MO+(-1)**j*O)==0
            target.extend([E,O])
    GE,GO=parts(data["G"])
    assert GO==0 and GE.subs(eta,0)==8*k+1
    raw.append(GE);strong.append(GE)
    numerator=[]
    b=4*k+1
    for j,co in enumerate(data["F"]):
        E,O=parts(co)
        if j%2:
            assert O==0
            numerator.append(sp.cancel(E*eta**((j-b)//2)))
        else:
            assert E==0
            numerator.append(sp.cancel(O*eta**((j-b+1)//2)))
    assert numerator[-1]==1
    assert numerator[0]==(8*k+1)/eta**(2*k)
    assert sp.cancel(sum(numerator)-2*GE/eta**(2*k))==0
    assert sp.cancel(sum(j*co for j,co in enumerate(numerator))-GE/eta**(2*k))==0

    def integers(equations):
        out=[]
        for expression in equations:
            _,poly=sp.Poly(expression,*variables,domain=sp.QQ).clear_denoms()
            _,primitive=poly.primitive()
            out.append(str(primitive.as_expr()))
        return out
    raw_length=int(4*k*sp.binomial(4*k,k)**2)
    return dict(schema_version=1,k=k,h=2*k,degree=8*k+2,
       variables=list(map(str,variables)),
       raw_rational_equations=list(map(str,raw)),
       raw_primitive_integer_equations=integers(raw),
       strong_rational_equations=list(map(str,strong)),
       strong_primitive_integer_equations=integers(strong),
       raw_equation_count=len(raw),strong_equation_count=len(strong),
       raw_labelled_length=raw_length,raw_quotient_length=raw_length//2,
       coefficient_reconstruction_constant_to_leading=list(map(str,numerator)),
       eta_nonzero_forced_by_equations=True,
       eta_squareclass="same as lambda; lambda=1/eta",
       raw_quotient_rational_point_list_computed=False)

def low_degree_exact_control(case):
    variables=list(map(sp.Symbol,case["variables"]))
    polys=list(map(sp.sympify,case["strong_primitive_integer_equations"]))
    GB=sp.groebner(polys,*variables,domain=sp.QQ,order="lex")
    expected=sp.Poly.from_list(finite.mixed.E18_COEFFICIENTS,eta,domain=sp.QQ)
    assert GB.is_zero_dimensional and len(GB.polys)==3
    assert sp.expand(GB.polys[-1].as_expr()-expected.as_expr()/20480)==0
    assert sp.degree(GB.polys[0].as_expr(),variables[0])==1
    assert sp.degree(GB.polys[1].as_expr(),variables[1])==1
    LM=[poly.monoms()[0]for poly in GB.polys]
    bounds=[min(ex[i]for ex in LM if ex[i]and
            all(ex[j]==0 for j in range(len(variables))if j!=i))
            for i in range(len(variables))]
    count=sum(1 for ex in itertools.product(*(range(B)for B in bounds))
        if not any(all(x>=y for x,y in zip(ex,lm))for lm in LM))
    assert count==18
    return dict(k=1,degree=10,coefficient_field="QQ",
       strong_quotient_length=18,known_eta_eliminant_exactly_recovered=True,
       complete_lex_basis=[str(poly.as_expr())for poly in GB.polys])

def run(write_system=None):
    cases=[quotient_case(k)for k in [1,2,3]]
    case26=cases[-1]
    assert case26["raw_labelled_length"]==580800
    assert case26["raw_quotient_length"]==290400
    assert case26["raw_equation_count"]==7 and case26["strong_equation_count"]==13
    control=low_degree_exact_control(cases[0])
    system=HERE/"degree26_balanced_rational_quotient_system.json"
    if write_system is not None:
        write_system.write_text(json.dumps(case26,indent=2,sort_keys=True)+"\n")
    if system.exists():
        assert json.loads(system.read_text())==case26
    paths=[
       "research/arithmetic/MIXED_RECIPROCAL_COUPLED_ODE.md",
       SOURCE.relative_to(ROOT).as_posix(),
       "research/arithmetic/RECIPROCAL_PARAMETER_FINITENESS_AND_INTEGRALITY.md",
       "research/arithmetic/BALANCED_RECIPROCAL_RATIONAL_QUOTIENT.md",
       Path(__file__).relative_to(ROOT).as_posix()]
    if system.exists():paths.append(system.relative_to(ROOT).as_posix())
    return dict(schema_version=1,status="PASS",
       python=platform.python_version(),sympy=sp.__version__,
       uniform_symbolic_controls=[{
           key:case[key]for key in [
               "k","h","degree","raw_equation_count","strong_equation_count",
               "raw_labelled_length","raw_quotient_length"]}
           for case in cases],
       low_degree_exact_QQ_control=control,
       degree26_system_path=system.relative_to(ROOT).as_posix(),
       uniform_quotient_proof_required=True,
       degree26_rational_points_computed=False,
       degree26_rational_existence_or_nonexistence_proved=False,
       sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                    for path in paths})

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-receipt",type=Path)
    parser.add_argument("--write-degree26-system",type=Path)
    args=parser.parse_args()
    result=run(args.write_degree26_system)
    text=json.dumps(result,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:args.write_receipt.write_text(text)
    print(text,end="")
if __name__=="__main__":main()
