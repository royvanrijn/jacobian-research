"""Bounded exact controls for safe recomposition and original filtered images.

The general classification and sharp filtered/jet results are written
proofs. Independent controls below contract the original three-term
cubic directly, without using the inverse shear or period transform.
"""
from __future__ import annotations
import argparse
from functools import lru_cache
import hashlib
import json
from math import factorial
from pathlib import Path
import platform
import sympy as sp

from gaussian_safe_recomposition import (
    U,Z,W,Y,T,recomposition_data,recompose_polynomial_period,
    original_period,filtered_basis,marker_lift,moment_jets,jet_bound,
)

if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")


def zero(value):
    assert sp.expand(value)==0,value


@lru_cache(maxsize=None)
def raw_cubic_monomial(i,j,k,n,functional):
    """L(Z^i W^j Y^k Pstar^n), by original multinomial Wick contraction."""
    answer=0
    for c in range(n+1):
        numerator_a=n+i-j+2*c
        numerator_b=n-i+j-4*c
        if numerator_a%2 or numerator_b%2:
            continue
        a,b=numerator_a//2,numerator_b//2
        if min(a,b)<0 or a+b+c!=n or (k+b)%2:
            continue
        value=(sp.Rational(factorial(n),factorial(a)*factorial(b)*factorial(c))
               *(-sp.Rational(1,2))**c*sp.factorial(j+a)*sp.factorial2(k+b-1))
        if functional=="sphere_3":
            value/=sp.factorial2(2*(j+a)+k+b+1)
        answer+=value
    return sp.cancel(answer)


@lru_cache(maxsize=None)
def raw_recomposed_monomial(i,j,k,m,polynomial,functional):
    answer=0
    for (n,),coefficient in sp.Poly(sp.expand(polynomial**m),U).terms():
        answer+=coefficient*raw_cubic_monomial(i,j,k,n,functional)
    return sp.cancel(answer)


def raw_recomposed_moment(marker,m,polynomial,functional):
    return sp.expand(sum(coefficient*raw_recomposed_monomial(i,j,k,m,polynomial,functional)
                         for (i,j,k),coefficient in sp.Poly(sp.expand(marker),Z,W,Y).terms()))


def raw_period(marker,polynomial,functional,order):
    return sp.expand(sum(raw_recomposed_moment(marker,m,polynomial,functional)*T**m/sp.factorial(m)
                         for m in range(order+1)))


def parameter_truncate(value,parameters,order):
    return sp.expand(sum(coefficient*sp.prod(p**n for p,n in zip(parameters,powers))
                         for powers,coefficient in sp.Poly(sp.expand(value),*parameters).terms()
                         if sum(powers)<=order))


def raw_perturbed_jet(H,parameters,order,polynomial,functional,m):
    answer=0
    for j in range(min(order,m)+1):
        marker=parameter_truncate(H**j,parameters,order)
        answer+=sp.binomial(m,j)*raw_recomposed_moment(marker,m-j,polynomial,functional)
    return parameter_truncate(answer,parameters,order)


def must_reject(callback):
    try:
        callback()
    except ValueError:
        return
    raise AssertionError("Invalid backend input was accepted")


def checks():
    f=U**2+7*U**5
    monomials=[Z**i*W**j*Y**k for i in range(4) for j in range(4-i) for k in range(4-i-j)]
    fixtures={
        "gaussian":[(1,1),(Y,2*T),(Z*W,1-6*T*T),(W**2,-4*T+30*T**3),
                    (Y*W,-420*T),(Y*W**2,-36*T*T+420*T**4)],
        "sphere_3":[(1,1),(Y,2*T/15),(Z*W,sp.Rational(1,3)-2*T*T/315),
                    (W**2,-4*T/105+2*T**3/9009),(Y*W,-4*T/99),
                    (Y*W**2,-4*T*T/1155+4*T**4/328185)],
    }
    a,b=sp.symbols("a b")
    H0=W**3+Y*W**2
    H=a*H0+b*(Y-Z**2)+b*b*Z
    for functional in fixtures:
        for marker in monomials:
            degree=sp.Poly(marker,Z,W,Y).total_degree()
            actual=original_period(marker,f,functional)
            zero(actual-raw_period(marker,f,functional,8))
            assert sp.Poly(actual,T).degree()<=3*degree//2
        for marker,expected in fixtures[functional]:
            zero(original_period(marker,f,functional)-expected)
            zero(raw_period(marker,f,functional,8)-expected)
        for polynomial in [f,3*U**3-2*U**4,U,2*U]:
            value,e=recomposition_data(polynomial)
            for degree in range(4):
                basis=filtered_basis(polynomial,degree,functional)
                assert len(basis)==3*degree//e+1
                for m,entry in enumerate(basis):
                    assert entry["time_degree"]==m
                    assert entry["marker_degree"]==(e*m+2)//3
                    expected=((-sp.Rational(1,2))**(e*m//3)*sp.factorial(e*m)
                              *value.nth(e)**m/sp.factorial(m))
                    if functional=="sphere_3":
                        expected/=sp.factorial2(2*e*m+1)
                    zero(entry["leading_coefficient"]-expected)
            for m in range(5):
                marker=marker_lift(T**m,polynomial,functional)
                assert sp.Poly(marker,Z,W,Y).total_degree()==(e*m+2)//3
                zero(original_period(marker,polynomial,functional)-T**m)
        # Independent original contractions of finite parameter jets, only order2.
        jets=moment_jets(H,(a,b),2,f,functional)
        bound=jet_bound(f,3,2)
        assert bound==11
        for m in range(1,bound+2):
            zero(jets.get(m,0)-raw_perturbed_jet(H,(a,b),2,f,functional,m))
        for r in [1,2]:
            B=jet_bound(f,3,r)
            n=B-r
            j=2*n
            c=(-sp.Rational(1,2))**(3*r) if r%2==0 else r*(-sp.Rational(1,2))**(3*r-1)
            expected=sp.binomial(B,r)*sp.factorial(j)*c
            if functional=="sphere_3":
                expected/=sp.factorial2(2*j+1)
            assert expected!=0
            actual=moment_jets(a*H0,(a,),r,f,functional)
            zero(sp.expand(actual[B]).coeff(a,r)-expected)
            assert max(actual)==B

    # Composition order and factorial normalization are tested exactly.
    g=1+2*T-3*T**4+5*T**10
    h=3*U**2+U**5
    q=U**2-2*U**3
    zero(recompose_polynomial_period(recompose_polynomial_period(g,h),q)
         -recompose_polynomial_period(g,q.subs(U,h)))
    assert recomposition_data(q.subs(U,h))[1]==4
    zero(recompose_polynomial_period(g,U)-g)
    zero(recompose_polynomial_period(g,2*U)-g.subs(T,2*T))
    short_source=1+T**2+T**4
    zero(recompose_polynomial_period(short_source,U**2+7*U**1000)
         -recompose_polynomial_period(short_source,U**2))
    for polynomial in [0,1,U+1]:
        must_reject(lambda polynomial=polynomial:recomposition_data(polynomial))
    must_reject(lambda:moment_jets(1+a*Z,(a,),2,f))
    for parameters in [(T,),(U,),(a,a),(sp.Integer(1),)]:
        must_reject(lambda parameters=parameters:moment_jets(a*Z,parameters,2,f))
    # A coefficient-7 corruption changes an actual original marked period.
    assert original_period(Y*W,U**2,"gaussian")!=original_period(Y*W,f,"gaussian")
    try:
        zero(original_period(Y*W,U**2,"gaussian")+420*T)
    except AssertionError:
        pass
    else:
        raise AssertionError("The deliberately corrupted recomposition fixture was accepted")
    return {
        "status":"passed","functional_variants":["gaussian","sphere_3"],
        "benchmark":"f(U)=U^2+7U^5, spatial phase degree15, exact filtered cutoff floor(3D/2)",
        "original_wick_controls":{"monomials":len(monomials),"spatial_degree_max":3,"moment_orders":[0,8],
                                  "method":"Direct multinomial contraction of original Pstar; no shear or transform in the reference calculation."},
        "hand_fixtures":{name:[{"marker":str(marker),"period":str(period)} for marker,period in values]
                         for name,values in fixtures.items()},
        "filtered_basis_controls":{"degree_range":[0,3],"recompositions":["U^2+7U^5","3U^3-2U^4","U","2U"]},
        "minimum_marker_degree":"ceil(e*m/3), monomial-period lifts m=0..4",
        "independent_perturbation_jets":{"parameters":2,"order":2,"moment_orders":[1,12],"budget":11},
        "sharp_jet_benchmark":{"direction":"W^3+YW^2","jet_orders":[1,2],"budgets":[5,11]},
        "composition_identity":"R_f composed with R_h = R_(f composed with h)",
        "high_degree_truncation":"f=U^2+7U^1000 agrees with U^2 on source period degree4",
        "negative_controls":["zero recomposition rejected by scoped backend","constant recomposition rejected",
                             "nonzero f(0) rejected","nonvanishing parameter-origin perturbation rejected",
                             "reserved-time/formal-variable parameter aliases rejected",
                             "duplicate and nonsymbol coefficient parameters rejected",
                             "coefficient-7 fixture corruption rejected"],
        "scope":"Bounded exact original controls. General safe-recomposition classification, all-field filtered surjectivity, minimal marker degrees and sharp all-r moment-jet budgets use the written proof.",
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=checks()
    here=Path(__file__).resolve()
    files=[here,here.with_name("gaussian_safe_recomposition.py"),
           here.with_name("moment_jet_tail_certificates.py"),
           here.with_name("GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md")]
    receipt={"status":"passed","checked_on":"2026-11-17",
             "environment":{"python":platform.python_version(),"sympy":sp.__version__},
             "sha256":{path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in files},
             "checks":result,"independent_replay":False,"formal_verification":False,"external_review":False}
    if args.output:
        args.output.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2))


if __name__=="__main__":
    main()
