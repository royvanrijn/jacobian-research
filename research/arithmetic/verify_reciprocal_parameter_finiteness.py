#!/usr/bin/env python3
"""Exact finite controls for reciprocal parameter finiteness/integrality.

The written uniform proof supplies all degrees. This replay verifies the
weighted degeneration on actual raw equations and computes full small raw
quotient lengths independently by finite-field Groebner bases.
"""
from __future__ import annotations
import argparse,hashlib,importlib.util,itertools,json,platform
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[2]
SOURCE=Path(__file__).with_name("verify_mixed_reciprocal_coupled_ode.py")
spec=importlib.util.spec_from_file_location("retained_mixed_ode",SOURCE)
mixed=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mixed)
w=mixed.w;r=mixed.r

def raw_system(h,p):
    q=h-p
    pc=sp.symbols("p0:"+str(p));qc=sp.symbols("q0:"+str(q))
    P=w**p+sum(pc[j]*w**j for j in range(p))
    Q=w**q+sum(qc[j]*w**j for j in range(q))
    variables=list(pc)+list(qc)+[r]
    weights=[p-j for j in range(p)]+[q-j for j in range(q)]+[1]
    s=sp.expand(P*Q)
    _,F,X,Y=mixed.from_critical(h,s)
    Xr=sp.Poly(sp.rem(X,P,w),w)
    Yr=sp.Poly(sp.rem(Y,Q,w),w)
    G=sp.cancel(sum(F[j]*r**j for j in range(len(F)))/(2*r))
    equations=[Xr.nth(j)for j in range(p)]+[
               Yr.nth(j)for j in range(q)]+[G]
    degrees=[2*h-j for j in range(p)]+[
             2*h-j for j in range(q)]+[2*h]
    return dict(h=h,p=p,q=q,P=P,Q=Q,s=s,variables=variables,
                weights=weights,F=F,X=X,Y=Y,G=G,equations=equations,
                degrees=degrees)

def weighted_top(expr,variables,weights,degree):
    poly=sp.Poly(sp.expand(expr),*variables,domain=sp.QQ)
    assert all(sum(i*j for i,j in zip(ex,weights))<=degree
               for ex,co in poly.terms())
    return sp.expand(sum(co*sp.prod(var**e for var,e in zip(variables,ex))
              for ex,co in poly.terms()
              if sum(i*j for i,j in zip(ex,weights))==degree))

def symbolic_control(h,p):
    d=raw_system(h,p);L=4*h+1
    P,Q,s=d["P"],d["Q"],d["s"]
    variables,weights=d["variables"],d["weights"]
    Z=sp.expand(L*sum(s.coeff(w,k)*(
           w**(h+k)/sp.Integer(2*h+2*k+1)-
           r*w**(h+k-1)/sp.Integer(2*h+2*k-1))
           for k in range(h+1)))
    for T in [d["X"],d["Y"]]:
        for j in range(2*h+1):
            assert weighted_top(T.coeff(w,j),variables,weights,2*h-j)==Z.coeff(w,j)
    gt=weighted_top(d["G"],variables,weights,2*h)
    assert sp.expand(2*gt-Z.subs(w,r))==0
    assert d["G"].subs(r,0)==L
    assert sp.expand(Z+2*w*sp.diff(Z,w)-
                     L*w**(h-1)*(w-r)*s)==0
    t=sp.symbols("t")
    g=t*Z.subs(w,t*t)
    assert sp.expand(sp.diff(g,t)-
                     L*t**(2*h-2)*(t*t-r)*s.subs(w,t*t))==0
    assert sp.Poly(g,t).LC()==1 and sp.Poly(g,t).degree()==L
    assert g.subs(t,0)==0
    expected=[sp.Poly(sp.rem(Z,P,w),w).nth(j)for j in range(p)]+[
              sp.Poly(sp.rem(Z,Q,w),w).nth(j)for j in range(h-p)]+[gt]
    actual=[weighted_top(expr,variables,weights,degree)
             for expr,degree in zip(d["equations"],d["degrees"])]
    assert all(sp.expand(x-y)==0 for x,y in zip(actual,expected))
    assert all(expr!=0 for expr in actual)
    denominator_primes=set()
    for expr in d["equations"]:
        for co in sp.Poly(expr,*variables,domain=sp.QQ).coeffs():
            denominator_primes.update(map(int,sp.factorint(co.q)))
    assert all(prime<=L for prime in denominator_primes)
    length=2*h*sp.binomial(2*h,p)*sp.binomial(2*h,h-p)
    assert sp.Rational(sp.prod(d["degrees"]),sp.prod(weights))==length
    return dict(h=h,p=p,q=h-p,variables=list(map(str,variables)),
       weights=weights,equation_weighted_degrees=d["degrees"],
       numerator_denominator_primes=sorted(denominator_primes),
       exact_raw_length=int(length),top_equations=[str(x)for x in actual],
       Z=str(Z),weighted_limit_and_single_critical_value_identities=True,
       r_zero_constraint=L)

def groebner_length_control(h,p):
    d=raw_system(h,p);variables=d["variables"]
    integers=[]
    for expression in d["equations"]:
        denominator,poly=sp.Poly(expression,*variables,domain=sp.QQ).clear_denoms()
        assert denominator%29
        integers.append(poly.as_expr())
    GB=sp.groebner(integers,*variables,modulus=29,order="lex")
    assert GB.is_zero_dimensional
    monomials=[poly.monoms()[0]for poly in GB.polys]
    bounds=[]
    for i in range(len(variables)):
        bounds.append(min(ex[i]for ex in monomials if ex[i]and
                      all(ex[j]==0 for j in range(len(variables))if j!=i)))
    count=sum(1 for exponent in itertools.product(*(range(x)for x in bounds))
            if not any(all(x>=y for x,y in zip(exponent,lm))
                       for lm in monomials))
    expected=int(2*h*sp.binomial(2*h,p)*sp.binomial(2*h,h-p))
    assert count==expected
    return dict(h=h,p=p,q=h-p,prime=29,quotient_length=count,
         standard_monomial_box=bounds,
         full_lex_basis=[str(poly.as_expr())for poly in GB.polys])

def run():
    controls=[symbolic_control(h,p)for h,p in [
        (1,0),(1,1),(2,0),(2,1),(3,1),(6,1),(6,3)]]
    lengths=[groebner_length_control(h,p)for h,p in [
        (1,0),(1,1),(2,0),(2,1),(2,2)]]
    degree26={str((p,6-p)):int(12*sp.binomial(12,p)*sp.binomial(12,6-p))
              for p in [1,2,3]}
    assert degree26=={"(1, 5)":114048,"(2, 4)":392040,"(3, 3)":580800}
    paths=[
       "research/arithmetic/MIXED_RECIPROCAL_COUPLED_ODE.md",
       SOURCE.relative_to(ROOT).as_posix(),
       "research/arithmetic/DEGREE_TWENTY_SIX_MIXED_2_10_GLOBAL_OBSTRUCTION.md",
       "research/arithmetic/RECIPROCAL_PARAMETER_FINITENESS_AND_INTEGRALITY.md",
       Path(__file__).relative_to(ROOT).as_posix(),
       "research/arithmetic/sources/stacks_02JN_regular_sequence_2026_11_18.html",
       "research/arithmetic/sources/stacks_00N7_cohen_macaulay_2026_11_18.html"]
    return dict(schema_version=1,status="PASS",python=platform.python_version(),
        sympy=sp.__version__,symbolic_controls=controls,
        independent_small_quotient_lengths=lengths,
        degree26_raw_lengths=degree26,
        written_uniform_proof_required=True,
        raw_lengths_include_degenerate_points_and_multiplicities=True,
        good_reduction_claim=False,rational_existence_claim=False,
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
