#!/usr/bin/env python3
"""Balanced degree26 local-square exclusion by a complete F37 census.

Default replay derives the 13 exact strong equations, binds C++ coefficient
reconstruction to symbolic F,X,Y,G, then runs two complete candidate methods.
The uniform finiteness note provides denominator control for all Q37 points.
"""
from __future__ import annotations
import argparse,hashlib,importlib.util,json,platform,shutil,subprocess,tempfile
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
CPP=HERE/"verify_degree_twenty_six_balanced_mod37.cpp"
SOURCE=HERE/"verify_mixed_reciprocal_coupled_ode.py"
spec=importlib.util.spec_from_file_location("balanced_mixed_source",SOURCE)
mixed=importlib.util.module_from_spec(spec);spec.loader.exec_module(mixed)
w=mixed.w;r=mixed.r
p2,p1,p0,q2,q1,q0=sp.symbols("p2 p1 p0 q2 q1 q0")
VARIABLES=(p2,p1,p0,q2,q1,q0,r)

def derive():
    P=w**3+p2*w*w+p1*w+p0
    Q=w**3+q2*w*w+q1*w+q0
    critical,F,X,Y=mixed.from_critical(6,sp.expand(P*Q))
    XP=sp.Poly(sp.rem(X,P**2,w),w)
    YQ=sp.Poly(sp.rem(Y,Q**2,w),w)
    plus=[sp.expand(XP.nth(j))for j in range(6)]
    minus=[sp.expand(YQ.nth(j))for j in range(6)]
    G=sp.cancel(sum(F[j]*r**j for j in range(len(F)))/(2*r))
    assert G.subs(r,0)==25
    for expression in plus:
        assert sp.Poly(expression,q2,q1,q0).total_degree()<=1
        assert sp.Poly(expression,r).degree()<=1
    for expression in minus:
        assert sp.Poly(expression,p2,p1,p0).total_degree()<=1
        assert sp.Poly(expression,r).degree()<=1
    polys=[sp.Poly(expr,*VARIABLES,domain=sp.QQ)
           for expr in plus+minus+[G]]
    denominators=set()
    for poly in polys:
        for coefficient in poly.coeffs():
            denominators.update(map(int,sp.factorint(coefficient.q)))
    assert all(prime<=25 for prime in denominators)
    return P,Q,F,X,Y,G,polys,denominators

def mod37(value):
    value=sp.Rational(value)
    return int(value.p)*pow(int(value.q),-1,37)%37

def bind_all_coefficients(binary,X,Y,G):
    """Unisolvent tensor controls, with exact symbolic degree bounds.

    Both the C++ coefficient construction and these symbolic polynomials
    are affine in each cubic coefficient group; XY are linear in r and
    G has r degree at most12. The C++ multiply/operator code has these
    bounds directly: PQ occurs once and the G loop uses powers r^0..r^12.
    Zero and three unit vectors determine an affine polynomial, and
    thirteen distinct r values determine a degree<=12 polynomial in F37.
    """
    vectors=[[sp.Poly(expr.coeff(w,j),*VARIABLES,domain=sp.QQ)
              for j in range(13)] for expr in (X,Y)]
    gp=sp.Poly(G,*VARIABLES,domain=sp.QQ)
    for poly in vectors[0]+vectors[1]+[gp]:
        for group in (VARIABLES[:3],VARIABLES[3:6]):
            assert sp.Poly(poly.as_expr(),*group).total_degree()<=1
    for poly in vectors[0]+vectors[1]:assert poly.degree(r)<=1
    assert gp.degree(r)<=12
    def terms(poly):return [(ex,mod37(co))for ex,co in poly.terms()]
    vectors=[[terms(poly)for poly in row]for row in vectors];gt=terms(gp)
    def evaluate(poly,point):
        total=0
        for ex,co in poly:
            value=co
            for x,e in zip(point,ex):value=value*pow(x,int(e),37)%37
            total=(total+value)%37
        return total
    process=subprocess.run([str(binary),"--basis-controls"],check=True,
                           capture_output=True,text=True)
    controls=[json.loads(line)for line in process.stdout.splitlines()]
    expected=[]
    for pi in range(-1,3):
        for qi in range(-1,3):
            for rv in range(13):
                point=[0]*7
                if pi>=0:point[pi]=1
                if qi>=0:point[3+qi]=1
                point[6]=rv;expected.append(point)
    assert [row['tuple']for row in controls]==expected
    for control in controls:
        point=control['tuple']
        assert control['X']==[evaluate(poly,point)for poly in vectors[0]]
        assert control['Y']==[evaluate(poly,point)for poly in vectors[1]]
        assert control['G']==evaluate(gt,point)
    return dict(result="PASS",prime=37,controls=len(controls),
        tensor_basis="zero or a coefficient unit vector in each cubic group; r=0,...,12",
        degree_bounds=dict(P_group=1,Q_group=1,r_in_XY=1,r_in_G=12),
        all_X_Y_G_coefficient_polynomials_bound=True)

def run():
    P,Q,F,X,Y,G,polys,denominators=derive()
    compiler=shutil.which("g++")
    if compiler is None:raise RuntimeError("g++ is required for the complete census")
    with tempfile.TemporaryDirectory(prefix="balanced26-mod37-")as folder:
        binary=Path(folder)/"census"
        subprocess.run([compiler,"-O2","-std=c++17",str(CPP),"-o",str(binary)],
                        check=True,capture_output=True,text=True)
        controls_process=subprocess.run([str(binary),"--controls"],check=True,
                         capture_output=True,text=True)
        controls=[json.loads(line)for line in controls_process.stdout.splitlines()]
        assert len(controls)==4
        for control in controls:
            substitution=dict(zip(VARIABLES,control["tuple"]))
            expectedX=[mod37(X.coeff(w,j).subs(substitution))for j in range(13)]
            expectedY=[mod37(Y.coeff(w,j).subs(substitution))for j in range(13)]
            expectedG=mod37(G.subs(substitution))
            assert control["X"]==expectedX
            assert control["Y"]==expectedY
            assert control["G"]==expectedG
        coefficient_binding=bind_all_coefficients(binary,X,Y,G)
        census_process=subprocess.run([str(binary)],check=True,
                         capture_output=True,text=True)
        censuses=[json.loads(line)for line in census_process.stdout.splitlines()]
    assert len(censuses)==2
    primary,secondary=censuses
    assert primary==dict(mode="strong_gauss",outer=1823508,
        candidates=1925562,points=[],consistent_ranks=[0,5,4021,1769940])
    assert secondary==dict(mode="weak_cramer",outer=1823508,
        candidates=12275266,points=[],consistent_ranks=None)
    assert primary["outer"]==37**3*36
    assert primary["candidates"]==sum(
        count*37**(3-rank)for rank,count in enumerate(primary["consistent_ranks"]))
    integer_equations=[]
    for poly in polys:
        denominator,integer=poly.clear_denoms()
        assert denominator%37
        _,primitive=integer.primitive()
        integer_equations.append(str(primitive.as_expr()))
    paths=[
       "research/arithmetic/MIXED_RECIPROCAL_COUPLED_ODE.md",
       SOURCE.relative_to(ROOT).as_posix(),
       "research/arithmetic/RECIPROCAL_PARAMETER_FINITENESS_AND_INTEGRALITY.md",
       "research/arithmetic/DEGREE_TWENTY_SIX_BALANCED_LOCAL_SQUARE_OBSTRUCTION.md",
       Path(__file__).relative_to(ROOT).as_posix(),
       CPP.relative_to(ROOT).as_posix()]
    return dict(schema_version=1,status="PASS",
       python=platform.python_version(),sympy=sp.__version__,prime=37,
       variables=list(map(str,VARIABLES)),strong_equations=integer_equations,
       coefficient_denominator_primes=sorted(denominators),
       polynomial_control_bindings=controls,censuses=censuses,
       unisolvent_coefficient_binding=coefficient_binding,
       whole_labelled_strong_closure_has_Q37_point=False,
       arbitrary_denominators_covered_by_uniform_integrality=True,
       rational_square_lambda_model_exists=False,
       nonsquare_lambda_necessary_condition="lambda is a unit at37 and Legendre(lambda,37)=-1",
       all_nonsquare_lambda_excluded=False,
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
