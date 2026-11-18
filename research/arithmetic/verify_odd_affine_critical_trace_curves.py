#!/usr/bin/env python3
"""Exact controls for odd-numerator affine critical-trace curves.

The note supplies the uniform operator proof. The default checker rebuilds
the complete h2/h6 rational curve data, checks all coefficients, replays the
rational-slope exception certificate, and verifies a good F29 control.
"""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
from math import isqrt
import sympy as sp
from sympy.polys.rings import ring
from verify_paired_critical_trace_chart import Algebra,charpoly

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
u,r,s,t=sp.symbols("u r s t")

def require(ok,message):
    if not ok:
        raise AssertionError(message)

def gammas(n):
    g=[sp.Integer(0),sp.Integer(1)]
    for j in range(2,n+1):
        g.append(sp.expand(t*g[-1]-g[-2]))
    return g

def operator_columns(h):
    L=4*h+1
    g=gammas(2*h+1)
    k=g[2*h+1]-g[2*h]
    K=L+k*r
    e=s+(t+2)*r
    A,B=sp.Integer(1),sp.Integer(0)
    columns=[]
    for j in range(h+1):
        J=A*((t+2)*r*(r-1)-2*s)+r*r*B*(t*s+(t+2)*(t-1)*r-t-2)
        op=sp.Poly(sp.expand(L*J-((t+2)*r*K-2*e*L)*r**j
                      -4*e*K*r*(j*r**(j-1) if j else 0)),r)
        require(op.nth(0)==op.nth(1)==0,"Operator not divisible by R²")
        kj=g[2*j+1]-g[2*j]
        require(sp.expand(op.nth(j+2)-(t+2)*(L*kj-(4*j+1)*k))==0,
                "Incorrect triangular diagonal")
        columns.append(sp.Poly(op.as_expr()/r**2,r))
        A,B=sp.expand(-r*A-t*r*r*B),sp.expand(t*A+(t*s+(t*t-1)*r)*B)
    return columns

def norm_and_numerator(h,coeff,T=t):
    """Return S,N for q=sum coeff[j]R^j over the chosen coefficient ring."""
    rems=[(sp.Integer(1),sp.Integer(0)),(sp.Integer(0),sp.Integer(1))]
    for j in range(2,h+1):
        a,b=rems[-1]
        rems.append((sp.expand(-(u*u-s*u)*b),sp.expand(a+T*u*b)))
    A=sum(coeff[j]*rems[j][0] for j in range(h+1))
    B=sum(coeff[j]*rems[j][1] for j in range(h+1))
    S=sp.Poly(sp.expand(A*A+T*u*A*B+(u*u-s*u)*B*B),u)
    require(S.degree()==2*h and S.LC()==1,"Norm polynomial degree")
    RN=sp.Poly((4*h+1)*(u-1)*S.as_expr(),u)
    N=sp.Poly(sum(RN.nth(j)*u**j/(2*j-1) for j in range(2*h+2)),u)
    return S,N

def universal_identity_control(h):
    qs=sp.symbols("q0:"+str(h))
    coeff=list(qs)+[sp.Integer(1)]
    q=sum(c*r**j for j,c in enumerate(coeff))
    S,N=norm_and_numerator(h,coeff)
    m=2*h+1;L=4*h+1
    C=[sp.Integer(0),sp.Integer(1)]
    T=[sp.Integer(2),s+t*r]
    for j in range(2,m+1):
        C.append(sp.expand((s+t*r)*C[-1]-r*r*C[-2]))
        T.append(sp.expand((s+t*r)*T[-1]-r*r*T[-2]))
    D=sp.expand(sum(N.nth(j)*(T[j]-(s+(t-2)*r)*C[j])
                    for j in range(m+1)))
    g=gammas(m);k=g[m]-g[m-1];K=L+k*r;e=s+(t+2)*r
    E=sp.expand(D-2*K*q*q)
    cols=operator_columns(h)
    op=sp.expand(r*r*sum(coeff[j]*cols[j].as_expr() for j in range(h+1)))
    require(sp.expand(e*r*sp.diff(E,r)+((t+2)*r/2-e)*E-q*op)==0,
            "Universal Euler/operator identity")
    require(sp.expand(sp.Poly(E,r).nth(2*h).subs(t,-2)
                      +sp.Rational(4,3)*L*h*(h-1)*s)==0,
            "Slope minus-two coefficient")
    require(sp.expand(N.nth(0)-L*qs[0]**2)==0,"Constant norm")
    return {"h":h,"universal_identity":True,
            "slope_minus_two_top_coefficient":str(-sp.Rational(4,3)*L*h*(h-1)*s)}

def encode(poly):
    return {str(j[0]):{"numerator":str(c.numer.as_expr()),
                      "denominator":str(c.denom.as_expr())}
            for j,c in poly.items()}

def build_curve(h):
    """Construct the monic kernel and accessory in Q(t)[s], exactly."""
    cols=operator_columns(h)
    K=sp.QQ.frac_field(t)
    P,x=ring("s",K)
    c=[P.zero]*h+[P.one]
    for j in range(h-1,-1,-1):
        num=sum((P.from_expr(cols[k].nth(j))*c[k] for k in range(j+1,h+1)),P.zero)
        c[j]=-num.quo_ground(K.from_sympy(cols[j].nth(j)))
        require(c[j].degree()<=h-j,"q coefficient degree in intercept")
    for i in range(h+1):
        require(sum((P.from_expr(cols[j].nth(i))*c[j] for j in range(h+1)),P.zero)==0,
                "Full rational-function matrix kernel")
    R,ss,uu=ring("s,u",K)
    cc=[R.from_expr(ci.as_expr()) for ci in c]
    theta=uu*K.from_sympy(t)
    norm=uu*uu-ss*uu
    rems=[(R.one,R.zero),(R.zero,R.one)]
    for j in range(2,h+1):
        a,b=rems[-1];rems.append((-norm*b,a+theta*b))
    A=sum((cc[j]*rems[j][0] for j in range(h+1)),R.zero)
    B=sum((cc[j]*rems[j][1] for j in range(h+1)),R.zero)
    S=A*A+theta*A*B+norm*B*B
    RN=(4*h+1)*(uu-1)*S
    accessory=P.zero
    for (i,j),v in RN.items():
        accessory+=x**i*P.ground_new(v/K.from_sympy(sp.Integer(2*j-1)))
    require(accessory.degree()<=2*h,"Accessory degree")
    return {"h":h,"q_coefficients_ascending":[encode(ci) for ci in c],
            "accessory_s_coefficients":encode(accessory),
            "kernel_checked":True}

def fixed_slope(h,T):
    cols=operator_columns(h)
    c=[sp.Integer(0)]*h+[sp.Integer(1)]
    for j in range(h-1,-1,-1):
        c[j]=sp.cancel(-sum(cols[k].nth(j).subs(t,T)*c[k] for k in range(j+1,h+1))
                       /cols[j].nth(j).subs(t,T))
    require(all(sp.expand(sum(cols[j].nth(i).subs(t,T)*c[j]
                            for j in range(h+1)))==0 for i in range(h+1)),
            "Fixed-slope kernel")
    S,N=norm_and_numerator(h,c,T)
    return c,S,N

def prime_divisors(n):
    out=[];p=2
    while p*p<=n:
        if n%p==0:
            out.append(p)
            while n%p==0:n//=p
        p+=1
    if n>1:out.append(n)
    return out

def valuation(n,p):
    v=0
    while n%p==0:v+=1;n//=p
    return v

def gamma_integer(n,T):
    a,b=0,1
    for _ in range(n):a,b=b,T*b-a
    return a

def diagonal_certificate(h):
    L=4*h+1;rows=[]
    def kj(j,T):return gamma_integer(2*j+1,T)-gamma_integer(2*j,T)
    for j in range(h):
        Lj=4*j+1;step=2*(h-j)
        den=[]
        for p in prime_divisors(Lj):
            gap=valuation(Lj,p)-valuation(L,p)
            require(gap<=0 or gap%step!=0,"Possible rational denominator")
            den.append({"prime":p,"valuation_gap":gap,"required_positive_multiple":step})
        const=L*(-1)**j-Lj*(-1)**h
        require(const!=0,"Integer-root constant zero")
        candidates=sorted({sgn*d for d in range(1,abs(const)+1) if const%d==0
                           for sgn in (-1,1)})
        values=[[T,L*kj(j,T)-Lj*kj(h,T)] for T in candidates]
        require([T for T,v in values if v==0]==[-2],"Rational diagonal exception")
        rows.append({"j":j,"leading_coefficient":-Lj,"constant":const,
                     "denominator_prime_tests":den,"integer_candidate_values":values,
                     "rational_roots":[-2]})
    return rows

def good29_control():
    h=6;p=29;sv=22
    c,S0,N0=fixed_slope(h,2)
    def red(x):
        x=sp.Rational(x)
        return int(x.p)*pow(int(x.q),-1,p)%p
    q=[red(ci.subs(s,sv)) for ci in c]
    N=[red(N0.nth(j).subs(s,sv)) for j in range(14)]
    S=[red(S0.nth(j).subs(s,sv)) for j in range(13)]
    B=Algebra(q,p);rr=B.reduce([0,1])
    theta=B.add(B.scale(rr,2),B.scale(B.one,sv))
    d=B.add(theta,B.scale(rr,-2))
    delta=B.add(B.mul(theta,theta),B.scale(B.mul(rr,rr),-4))
    C=[B.zero,B.one]
    for j in range(2,14):
        C.append(B.add(B.mul(theta,C[-1]),B.mul(B.mul(rr,rr),C[-2]),-1))
    a=B.zero
    for j,nj in enumerate(N):a=B.add(a,B.scale(C[j],nj))
    M=B.mul(B.mul(a,a),d)
    cm=charpoly(B.columns(M),p)
    def disc(cs):
        x=sp.Symbol("x")
        return int(sp.discriminant(sum(v*x**j for j,v in enumerate(cs)),x))%p
    def norm(v):
        return int(charpoly(B.columns(v),p)[0])%p
    units={"q0":q[0],"Discq":disc(q),"DiscS":disc(S),"S1":sum(S)%p,
           "norm_secant":norm(a),"norm_d":norm(d),"norm_Delta":norm(delta),
           "norm_M":norm(M),"disc_value_polynomial":disc(cm),
           "accessory_s_derivative":red(sp.diff(N0.eval(1),s).subs(s,sv))}
    require(all(units.values()),"Good-open or Hensel unit failed")
    require(any(a[1:]),"Secant accidentally scalar")
    require(B.mul(rr,rr)!=B.scale(B.one,B.mul(rr,rr)[0]),"Pair product scalar")
    require(sum(N)%p==0,"Distinguished zero value")
    # Full discriminant in F29(w), with field subtraction normalizing units.
    z,w=sp.symbols("z w")
    field=sp.GF(p).frac_field(w)
    F=sp.Poly(sum(v*z**(2*j) for j,v in enumerate(N))-w*z,z,domain=field)
    disc_actual=F.rep.discriminant()
    value_poly=sum(v*w**(2*j) for j,v in enumerate(cm))
    expected=field.from_sympy(25**25*(w*value_poly)**2)
    require(disc_actual-expected==field.zero,"Full finite-field square identity")
    return {"h":h,"prime":p,"s":sv,"t":2,
            "q_coefficients_ascending":q,"N_coefficients_ascending":N,
            "S_coefficients_ascending":S,"secant_coefficients_ascending":a,
            "squared_value_polynomial_ascending":cm,"good_units":units,
            "full_discriminant_identity":True,"secant_nonconstant":True,
            "paired_product_nonconstant":True}

def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-receipt",type=Path,
        default=HERE/"odd_affine_critical_trace_curves_verification.json")
    args=parser.parse_args()
    retained=json.loads((HERE/"odd_affine_critical_trace_curve_inputs.json").read_text())
    curves={str(h):build_curve(h) for h in (2,6)}
    require(curves==retained["curves"],"Explicit curve input mismatch")
    positive=good29_control()
    require(positive==retained["good29_control"],"Good control input mismatch")
    # Uniform nonzero-accessory specialization, controlled in retained degrees.
    for h in (2,6):
        c,S,N=fixed_slope(h,2);L=4*h+1
        require(all(sp.expand(ci.subs(s,0)-sp.binomial(h,j)*L**(h-j))==0
                    for j,ci in enumerate(c)),"Degenerate closure kernel")
        require(sp.expand(N.as_expr().subs(s,0)-(u+L)**(2*h+1))==0,
                "Nonzero accessory specialization")
    identities=[universal_identity_control(h) for h in (1,2,3)]
    certificates={str(h):diagonal_certificate(h) for h in (2,6)}
    paths=[
        HERE/"ODD_AFFINE_CRITICAL_TRACE_CURVES.md",
        HERE/"odd_affine_critical_trace_curve_inputs.json",
        HERE/"verify_paired_critical_trace_chart.py",
        HERE/"PAIRED_CRITICAL_TRACE_CHART.md",
        HERE/"UNIFORM_REAL_AFFINE_TRACE_OBSTRUCTION.md",
        HERE/"MIXED_RECIPROCAL_COUPLED_ODE.md",
        HERE/"sources/jones-beauville-1304.5444v1.pdf",
        HERE/"ALTERNATING_HURWITZ_DISCRIMINANT_BRIDGE.md"]
    result={"schema_version":1,"result":"PASS","python":platform.python_version(),
        "sympy":sp.__version__,"checker_sha256":digest(Path(__file__)),
        "source_sha256":{str(p.relative_to(ROOT)):digest(p) for p in paths},
        "symbolic_universal_controls":identities,
        "curve_accessory_s_degrees":{"2":4,"6":12},
        "rational_slope_exception_certificates":certificates,
        "good29_control":positive,
        "scope":"Written uniform operator equivalence plus complete good-open "
        "conditions. Exact h2/h6 curve data rebuilt, rational slope -2 excluded, "
        "good h6 F29 point and simple accessory root establish a Q29/numberfield "
        "point. No rational Q point, global obstruction or exhaustive Q-point "
        "determination is asserted.",
        "formal_verification":False,"independent_external_replay":False,
        "external_human_review":False}
    args.write_receipt.write_text(json.dumps(result,indent=2)+"\n")
    print("PASS: exact odd affine-trace curve, exception and good-open controls.")
if __name__=="__main__":
    main()

