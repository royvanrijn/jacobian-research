#!/usr/bin/env python3
"""Exact reciprocal degree-twelve obstruction and weighted-Rolle controls.

The written note supplies the uniform proof. This script derives the two
degree-twelve equations, their exact resultant, a projective mod-31
obstruction, and independent Sturm/interval controls over the real numbers.
"""
from __future__ import annotations
import argparse,hashlib,json,platform
from pathlib import Path
import sympy as sp

ROOT=Path(__file__).resolve().parents[2]
a,c,w,v=sp.symbols("a c w v")
E_COEFFICIENTS=[
1125,-225,-16335,-43362,207009,-100551,2204367,-15624026,
48834247,-83352841,67437755,-3843328,2647888]
MOD31_VALUES=[23,2,1,3,18,13,15,1,25,8,3,13,7,14,26,18,
11,15,10,11,12,14,11,19,14,5,15,5,13,19,8]

def inverse_operator(rhs,sign):
    t=sp.symbols("t")
    shifted=sp.Poly(sp.expand(rhs.subs(w,t-2*sign)),t)
    return sp.expand(sum(co*(w+2*sign)**j/sp.Integer(2*j+1)
                         for (j,),co in shifted.terms()))

def primitive(expr):
    num,den=sp.fraction(sp.cancel(expr))
    assert not den.free_symbols
    return sp.Poly(num,a,c,domain=sp.QQ).primitive()[1]

def interval_horner(poly,interval):
    lo=hi=sp.Integer(0)
    for coefficient in poly.all_coeffs():
        products=[lo*interval[0],lo*interval[1],
                  hi*interval[0],hi*interval[1]]
        lo=min(products)+coefficient
        hi=max(products)+coefficient
    return lo,hi

def general_controls():
    """Symbolic operator controls for three even critical degrees."""
    records=[]
    for k in [2,3,4]:
        degree=2*k
        coefficients=sp.symbols("f0:"+str(degree))
        F=v**degree+sum(coefficients[j]*v**j for j in range(degree))
        reciprocal=sp.expand(v**degree*F.subs(v,1/v))
        # Direct Laurent definitions, represented without Laurent division.
        Xv=sp.cancel((F+v*F.subs(v,1/v))/(v+1))
        Yv=sp.cancel((F-v*F.subs(v,1/v))/(v-1))
        C=[sp.Integer(1),w-1];D=[sp.Integer(1),w+1]
        for j in range(1,degree-1):
            C.append(sp.expand(w*C[-1]-C[-2]))
            D.append(sp.expand(w*D[-1]-D[-2]))
        coeff=list(coefficients)+[sp.Integer(1)]
        X=sp.expand(coeff[0]+sum(coeff[j]*C[j-1]for j in range(1,degree+1)))
        Y=sp.expand(-coeff[0]+sum(coeff[j]*D[j-1]for j in range(1,degree+1)))
        assert sp.cancel(X.subs(w,v+1/v)-Xv)==0
        assert sp.cancel(Y.subs(w,v+1/v)-Yv)==0
        assert sp.Poly(X,w).degree()==degree-1
        assert sp.Poly(Y,w).degree()==degree-1
        R=sp.expand(2*v*sp.diff(F,v)-F)
        LX=X+2*(w+2)*sp.diff(X,w)
        LY=Y+2*(w-2)*sp.diff(Y,w)
        assert sp.cancel((v-1)*LX.subs(w,v+1/v)+
                         (v+1)*LY.subs(w,v+1/v)-2*R)==0
        records.append(dict(k=k,critical_degree=degree,
                            operator_identity=True))
    return records

def run():
    q1,q0=sp.symbols("q1 q0")
    C2=w*w-w-1;D2=w*w+w-1
    Q=w*w+q1*w+q0
    error=sp.Poly(sp.expand(4*(w-2)*(w+a)*sp.diff(Q,w)+
             (3*w+a-4)*Q-11*C2*(w-c)),w)
    solution={}
    divisors=[]
    for degree,variable in [(2,q1),(1,q0)]:
        coefficient=sp.expand(error.nth(degree).subs(solution))
        divisor=sp.diff(coefficient,variable)
        divisors.append(int(divisor))
        solution[variable]=sp.cancel(-coefficient.subs(variable,0)/divisor)
    assert divisors==[7,3]
    Q=sp.expand(Q.subs(solution));s=sp.expand((w-c)*Q)
    remaining=sp.Poly(error.as_expr().subs(solution),w)
    assert remaining.degree()<=0
    X=inverse_operator(11*D2*s,1)
    Y=inverse_operator(11*C2*s,-1)
    e0=primitive(remaining.nth(0))
    e1=primitive(sp.cancel(X.subs(w,c)/(c+2)))
    assert sp.Poly(X.subs(w,c),c).rem(sp.Poly(c+2,c)).is_zero
    # Bind both operators to the actual monic degree-six polynomial.
    Sigma=sp.Poly(sp.expand(v**3*s.subs(w,v+1/v)),v)
    assert Sigma.degree()==6 and Sigma.LC()==Sigma.nth(0)==1
    F=sp.expand(sum(11*Sigma.nth(j)*v**j/sp.Integer(2*j-1)
                    for j in range(7)))
    assert sp.Poly(F,v).LC()==1 and F.subs(v,0)==-11
    assert sp.expand(2*v*sp.diff(F,v)-F-11*Sigma.as_expr())==0
    assert sp.cancel((F+v*F.subs(v,1/v))/(v+1)-
                     X.subs(w,v+1/v))==0
    assert sp.cancel((F-v*F.subs(v,1/v))/(v-1)-
                     Y.subs(w,v+1/v))==0
    # The residual Y-equation is precisely the sole constant compatibility.
    assert sp.rem(Y-Q**2*(w+a),e0.as_expr(),a)==0
    eliminant=sp.Poly.from_list(E_COEFFICIENTS,c,domain=sp.ZZ)
    resultant=sp.resultant(e0.as_expr(),e1.as_expr(),a)
    assert sp.expand(resultant+57153600*eliminant.as_expr())==0
    values=[int(eliminant.eval(j))%31 for j in range(31)]
    assert values==MOD31_VALUES and all(values)
    assert int(eliminant.LC())%31==9
    assert sp.gcd(eliminant,eliminant.diff()).degree()==0

    # Independent complete lex basis: every real candidate has one a-coordinate.
    groebner=sp.groebner([e0.as_expr(),e1.as_expr()],a,c,
                        order="lex",domain=sp.QQ)
    assert len(groebner.polys)==2 and groebner.is_zero_dimensional
    assert sp.expand(groebner.polys[1].as_expr()-
                     eliminant.as_expr()/1125)==0
    first=groebner.polys[0].as_expr()
    assert sp.diff(first,a)==1
    A=sp.expand(a-first)
    assert a not in A.free_symbols
    Q2=sp.Poly(sp.expand(Q.subs({a:A,w:2})),c,
                domain=sp.QQ).rem(eliminant)
    real_intervals=sp.polys.polytools.intervals(
        eliminant,eps=sp.Rational(1,10)**30)
    assert len(real_intervals)==2
    assert eliminant.count_roots(-sp.oo,sp.oo)==2
    real_records=[]
    for interval,multiplicity in real_intervals:
        assert multiplicity==1 and interval[0]>2
        low,high=interval_horner(Q2,interval)
        assert high<0
        real_records.append(dict(
            c_interval=[str(x)for x in interval],
            Q_at_two_integer_bounds=[int(sp.floor(low)),int(sp.ceiling(high))],
            Q_at_two_strictly_negative=True))
    assert [record["Q_at_two_integer_bounds"]for record in real_records]==[[-10,-9],[-6,-5]]

    # A shorter, purely polynomial proof that no candidate has c<=0.
    x=sp.symbols("x")
    quartic=1125*x**4-16335*x**2+207009
    rest=sp.Poly(sp.expand(eliminant.as_expr().subs(c,-x)-x**8*quartic),x)
    assert all(coefficient>=0 for coefficient in rest.all_coeffs())
    assert rest.nth(0)>0
    deficit=4*1125*207009-16335**2
    assert deficit==664708275 and deficit>0

    paths=[
      "research/arithmetic/PAIRED_CRITICAL_RECIPROCAL_OBSTRUCTION.md",
      "research/arithmetic/MIXED_RECIPROCAL_COUPLED_ODE.md",
      "research/arithmetic/RECIPROCAL_DEGREE_TWELVE_OBSTRUCTION.md",
      "research/arithmetic/verify_reciprocal_degree_twelve.py"]
    pins={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
          for path in paths}
    return dict(schema_version=1,status="PASS",
       python=platform.python_version(),sympy=sp.__version__,
       uniform_controls=general_controls(),
       critical_degree=6,cover_degree=12,top_down_divisors=divisors,
       Q=str(Q),equations=[str(e0.as_expr()),str(e1.as_expr())],
       F=str(F),X=str(X),Y=str(Y),
       eliminant_coefficients=E_COEFFICIENTS,
       resultant_multiplier=-57153600,
       projective_mod31=dict(affine_values=values,value_at_infinity=9),
       complete_lex_basis=[str(p.as_expr())for p in groebner.polys],
       real_candidate_count=2,real_candidates=real_records,
       negative_axis_certificate=dict(
         quartic=str(quartic),quadratic_discriminant_deficit=deficit,
         remainder=str(rest.as_expr())),
       written_uniform_proof_required=True,
       rational_good_reciprocal_models_exist=False,
       all_real_good_reciprocal_models_exist=False,
       general_nonreciprocal_existence_claim=False,
       sha256=pins)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-receipt",type=Path)
    args=parser.parse_args()
    result=run()
    text=json.dumps(result,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:args.write_receipt.write_text(text)
    print(text,end="")
if __name__=="__main__":main()
