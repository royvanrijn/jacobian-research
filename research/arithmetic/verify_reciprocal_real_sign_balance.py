#!/usr/bin/env python3
"""Exact finite controls for the uniform reciprocal real-sign balance proof."""
from __future__ import annotations
import argparse,hashlib,json,platform
from pathlib import Path
import sympy as sp
ROOT=Path(__file__).resolve().parents[2]
w,v,u,r=sp.symbols("w v u r")

def chebyshev(maximum):
    C=[sp.Integer(1),w-1];D=[sp.Integer(1),w+1]
    for j in range(1,maximum):
        C.append(sp.expand(w*C[-1]-C[-2]))
        D.append(sp.expand(w*D[-1]-D[-2]))
    return C[:maximum+1],D[:maximum+1]

def open_positive_root_count(expression):
    poly=sp.Poly(expression,w,domain=sp.QQ)
    endpoint_multiplicity=0
    while poly.eval(2)==0:
        poly=poly.exquo(sp.Poly(w-2,w))
        endpoint_multiplicity+=1
    return int(poly.count_roots(2,sp.oo)),endpoint_multiplicity

def operator_controls():
    C,D=chebyshev(12)
    records=[]
    for h in [1,2,3,6]:
        coefficients=sp.symbols("s0:"+str(h))
        s=w**h+sum(coefficients[j]*w**j for j in range(h))
        Sigma=sp.Poly(sp.expand(v**h*s.subs(w,v+1/v)),v)
        L=4*h+1
        F=sp.expand(sum(L*((Sigma.nth(j-1)if j else 0)-
                   r*(Sigma.nth(j)if j<=2*h else 0))*v**j/
                   sp.Integer(2*j-1)for j in range(2*h+2)))
        fc=sp.Poly(F,v)
        X=sp.expand(fc.nth(0)+sum(fc.nth(j)*C[j-1]
                                  for j in range(1,2*h+2)))
        Y=sp.expand(-fc.nth(0)+sum(fc.nth(j)*D[j-1]
                                   for j in range(1,2*h+2)))
        assert sp.expand(2*v*sp.diff(F,v)-F-
                         L*(v-r)*Sigma.as_expr())==0
        assert sp.expand(X+2*(w+2)*sp.diff(X,w)-
                         L*(D[h]-r*D[h-1])*s)==0
        assert sp.expand(Y+2*(w-2)*sp.diff(Y,w)-
                         L*(C[h]-r*C[h-1])*s)==0
        records.append(dict(h=h,actual_coefficient_operator_identity=True))
    return records

def run():
    C,D=chebyshev(12)
    records=[]
    t=sp.symbols("t")
    for h in range(1,13):
        threshold=sp.Rational(2*h+1,2*h-1)
        assert C[h].subs(w,2)==C[h-1].subs(w,2)==1
        assert D[h].subs(w,2)==2*h+1
        assert D[h-1].subs(w,2)==2*h-1
        ratio_records=[]
        for label,Z in [("C",C),("D",D)]:
            # Positive numerator of (Z_h/Z_(h-1))' for w>=2.
            W=sp.Poly(sp.expand((sp.diff(Z[h],w)*Z[h-1]-
                    Z[h]*sp.diff(Z[h-1],w)).subs(w,t+2)),t)
            assert W.nth(0)>0
            assert all(x>=0 for x in W.all_coeffs())
            ratio_records.append(dict(sequence=label,
                    derivative_shifted_coefficients=list(map(str,W.all_coeffs()))))
        tests=[]
        sample_parameters=[sp.Rational(1,2),sp.Integer(1),
                           (1+threshold)/2,threshold,threshold+1]
        for parameter in sample_parameters:
            beta,endpoint_beta=open_positive_root_count(C[h]-parameter*C[h-1])
            alpha,endpoint_alpha=open_positive_root_count(D[h]-parameter*D[h-1])
            assert beta==int(bool(parameter>1))
            assert alpha==int(bool(parameter>threshold))
            assert endpoint_beta==int(bool(parameter==1))
            assert endpoint_alpha==int(bool(parameter==threshold))
            tests.append(dict(r=str(parameter),beta_roots_above_two=beta,
                  alpha_roots_above_two=alpha,beta_endpoint=endpoint_beta,
                  alpha_endpoint=endpoint_alpha))
        records.append(dict(h=h,threshold=str(threshold),tests=tests,
                            positive_derivative_controls=ratio_records))
    integer_records=[]
    for h in range(1,41):
        for eps_beta,eps_alpha in [(0,0),(1,0),(1,1)]:
            permitted=[p for p in range(h+1)
                         if h-p<=p+eps_beta and p<=h-p+1+eps_alpha]
            if h%2==0:
                expected=[h//2]
                if eps_alpha:expected.append(h//2+1)
            else:
                expected=[(h+1)//2]
                if eps_beta:expected.insert(0,(h-1)//2)
            assert permitted==expected
            integer_records.append(dict(h=h,epsilon_beta=eps_beta,
                 epsilon_alpha=eps_alpha,permitted_p=permitted))
    degree26={}
    for label,eb,ea in [("r<=1",0,0),("1<r<=13/11",1,0),("r>13/11",1,1)]:
        degree26[label]=[[p,6-p]for p in range(7)
                        if 6-p<=p+eb and p<=7-p+ea]
    assert degree26=={"r<=1":[[3,3]],"1<r<=13/11":[[3,3]],
                     "r>13/11":[[3,3],[4,2]]}

    # Sharpness: the already-certified rational degree-ten model.
    N=sum(sp.Rational(numerator,denominator)*u**j
       for j,(numerator,denominator)in enumerate([
           (16,5625),(1804,5625),(-14006,5625),(3129,625),(-96,25),(1,1)]))
    R=sp.Poly(2*u*sp.diff(N,u)-N,u)
    S,remainder=R.div(sp.Poly(u-1,u))
    assert remainder.is_zero and N.subs(u,1)==sp.diff(N,u).subs(u,1)==0
    assert sp.gcd(R,R.diff()).degree()==0
    assert R.count_roots(0,sp.oo)==5
    rho=sp.Rational(2,15)
    Sigma=sp.Poly(sp.expand(S.as_expr().subs(u,rho*v)/(9*rho**4)),v)
    assert Sigma.LC()==Sigma.nth(0)==1
    P=w*w+Sigma.nth(3)*w+Sigma.nth(2)-2
    assert sp.expand(v*v*P.subs(w,v+1/v)-Sigma.as_expr())==0
    assert sp.Poly(P,w).count_roots(2,sp.oo)==2
    F=sp.Poly(sp.expand(N.subs(u,rho*v)/rho**5),v)
    X=sp.expand(F.nth(0)+sum(F.nth(j)*C[j-1]for j in range(1,6)))
    assert sp.expand(X-P**2)==0
    sharpness=dict(h=2,r="15/2",p=2,q=0,
          all_five_critical_u_roots_positive=True,
          normalized_residual_P=str(P),plus_square_identity=True)

    paths=[
      "research/arithmetic/MIXED_RECIPROCAL_COUPLED_ODE.md",
      "research/arithmetic/RECIPROCAL_DEGREE_TWELVE_OBSTRUCTION.md",
      "research/arithmetic/RATIONAL_DEGREE_TEN_ALTERNATING_PENCIL.md",
      "research/arithmetic/RECIPROCAL_REAL_SIGN_BALANCE.md",
      "research/arithmetic/verify_reciprocal_real_sign_balance.py"]
    return dict(schema_version=1,status="PASS",python=platform.python_version(),
        sympy=sp.__version__,threshold_controls=records,
        actual_operator_controls=operator_controls(),
        integer_balance_controls=integer_records,degree26_counts=degree26,
        sharpness=sharpness,written_uniform_proof_required=True,
        sufficiency_claim=False,general_rational_degree26_claim=False,
        sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                for path in paths})

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-receipt",type=Path)
    args=parser.parse_args()
    result=run();text=json.dumps(result,sort_keys=True,indent=2)+"\n"
    if args.write_receipt:args.write_receipt.write_text(text)
    print(text,end="")
if __name__=="__main__":main()
