#!/usr/bin/env python3
"""Global Q_61 obstruction for the degree-26 mixed (2,10) equations.

Two independent complete finite-field implementations verify the affine
and weighted-leading censuses. The companion written valuation lemma
upgrades them to exclusion of every Q_61 point, with arbitrary denominators.
Default replay needs SymPy and g++; the C++ binary is built in a temp folder.
"""
from __future__ import annotations
import argparse,hashlib,json,platform,shutil,subprocess,tempfile
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
SYSTEM=HERE/"degree26_mixed_2_10_system.json"
CPP=HERE/"verify_degree_twenty_six_mixed_2_10_finite_field.cpp"
PRIME=61
WEIGHTS=(1,2,1,1)
DEGREES=(7,6,11,12)
a,b,c,r=sp.symbols("a b c r")

def load():
    data=json.loads(SYSTEM.read_text())
    assert data["h"]==6 and data["degree"]==26
    assert data["variables"]==["a","b","c","r"]
    assert data["sign_degrees"]==[2,10]
    assert data["top_down_divisors"]==[21,17,13,9,5]
    polys=[sp.Poly(sp.sympify(text),a,b,c,r,domain=sp.ZZ)
           for text in data["equations"]]
    assert len(polys)==4
    assert all(abs(int(co))<2**63 and max(ex)<16
               for poly in polys for ex,co in poly.terms())
    for poly in polys[:2]:
        assert poly.degree(c)<=1 and poly.degree(r)<=1
    assert [poly.total_degree()for poly in polys]==[6,6,11,12]
    leading=[]
    for poly,expected in zip(polys,DEGREES):
        degree=max(sum(x*y for x,y in zip(exponent,WEIGHTS))
                   for exponent,coefficient in poly.terms())
        assert degree==expected
        leading.append(sp.Poly.from_dict(
             {exponent:coefficient for exponent,coefficient in poly.terms()
              if sum(x*y for x,y in zip(exponent,WEIGHTS))==degree},
             (a,b,c,r),domain=sp.ZZ))
    return polys,leading

def python_census(polys):
    """Solve both first linear c-equations, keeping every singular case."""
    p=PRIME
    terms=[[(tuple(exponent),int(coefficient)%p)
             for exponent,coefficient in poly.terms()]for poly in polys]
    points=[];candidates=0
    for av in range(p):
        A=[pow(av,j,p)for j in range(16)]
        for bv in range(p):
            B=[pow(bv,j,p)for j in range(16)]
            rows=[]
            for j in [0,1]:
                rows.append({
                    (ce,re):sum(co*A[ae]*B[be]
                       for(ae,be,cc,rr),co in terms[j]
                       if cc==ce and rr==re)%p
                    for ce in [0,1]for re in [0,1]})
            C0,C1=rows
            for rv in range(p):
                ca=(C0[1,1]*rv+C0[1,0])%p
                cb=(C0[0,1]*rv+C0[0,0])%p
                da=(C1[1,1]*rv+C1[1,0])%p
                db=(C1[0,1]*rv+C1[0,0])%p
                if ca==0:
                    if cb:continue
                    if da==0:
                        if db:continue
                        choices=range(p)
                    else:choices=[-db*pow(da,-1,p)%p]
                else:
                    cv=-cb*pow(ca,-1,p)%p
                    if(da*cv+db)%p:continue
                    choices=[cv]
                R=[pow(rv,j,p)for j in range(16)]
                for cv in choices:
                    candidates+=1
                    C=[pow(cv,j,p)for j in range(16)]
                    if all(sum(co*A[ae]*B[be]*C[ce]*R[re]
                         for(ae,be,ce,re),co in poly)%p==0
                         for poly in terms[2:]):
                        points.append([av,bv,cv,rv])
    return dict(first_two_candidates=candidates,points=points)

def cpp_censuses(polys):
    compiler=shutil.which("g++")
    if compiler is None:raise RuntimeError("g++ is required for the independent census")
    lines=[str(PRIME)]
    for poly in polys:
        lines.append(str(len(poly.terms())))
        for exponent,coefficient in poly.terms():
            lines.append(" ".join([str(coefficient)]+list(map(str,exponent))))
    data="\n".join(lines)+"\n"
    with tempfile.TemporaryDirectory(prefix="mixed26-mod61-")as folder:
        binary=Path(folder)/"census"
        subprocess.run([compiler,"-O2","-std=c++17",str(CPP),"-o",str(binary)],
                       check=True,capture_output=True,text=True)
        process=subprocess.run([str(binary)],input=data,check=True,
                               capture_output=True,text=True)
    reports=[json.loads(line)for line in process.stdout.splitlines()]
    assert len(reports)==2 and [item.pop("mode")for item in reports]==[0,1]
    return reports

def run():
    polys,leading=load()
    affine=python_census(polys)
    infinity=python_census(leading)
    assert affine==dict(first_two_candidates=3758,points=[])
    assert infinity==dict(first_two_candidates=10921,points=[[0,0,0,0]])
    cpp=cpp_censuses(polys)
    assert cpp==[affine,infinity]
    # The half-integral valuation case has only b nonzero on reduction.
    only_b=sp.Poly(leading[1].as_expr().subs({a:0,c:0,r:0}),b)
    assert only_b.as_expr()==124215*b**3
    assert 124215%PRIME==19
    paths=[
       SYSTEM.relative_to(ROOT).as_posix(),
       "research/arithmetic/MIXED_RECIPROCAL_COUPLED_ODE.md",
       "research/arithmetic/verify_mixed_reciprocal_coupled_ode.py",
       "research/arithmetic/DEGREE_TWENTY_SIX_MIXED_2_10_GLOBAL_OBSTRUCTION.md",
       Path(__file__).relative_to(ROOT).as_posix(),
       CPP.relative_to(ROOT).as_posix()]
    return dict(schema_version=1,status="PASS",
       python=platform.python_version(),sympy=sp.__version__,
       prime=PRIME,variables=["a","b","c","r"],weights=list(WEIGHTS),
       weighted_degrees=list(DEGREES),
       affine_python=affine,weighted_leading_python=infinity,
       affine_cpp=cpp[0],weighted_leading_cpp=cpp[1],
       independent_enumerations_agree=True,
       weighted_leading_forms=[str(poly.as_expr())for poly in leading],
       half_integral_chart=dict(only_b_polynomial=str(only_b.as_expr()),
                                coefficient_mod61=19),
       valuation_argument="The written weighted residue lemma excludes nonintegral Q61 coordinates after a totally ramified quadratic rescaling; integral coordinates are excluded by the affine census.",
       covers_arbitrary_denominators=True,
       saturation_or_good_reduction_hypothesis=False,
       conclusion="The entire four-equation affine closure has no Q61 point, hence no rational point. No nondegenerate rational degree26 reciprocal model with sign degrees2+10 in either orientation exists.",
       remaining_cases=["squarelambda4+8","squarelambda6+6","nonsquarelambda6+6"],
       sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                  for path in paths})

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-receipt",type=Path)
    args=parser.parse_args()
    result=run();text=json.dumps(result,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:args.write_receipt.write_text(text)
    print(text,end="")
if __name__=="__main__":main()
