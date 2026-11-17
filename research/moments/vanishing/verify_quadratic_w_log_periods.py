#!/usr/bin/env python3
"""Algebra for nondegenerate quadratic W coefficients; analytic lemma is written."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import time,json,hashlib
import sympy as S
from verify_moment_transfers import add,mul,scale,gaussian


def cut(p,n):return {m:c for m,c in p.items() if m[0]<=n and c}
def powp(p,k,n):
    v={(0,0):Q(1)}
    for _ in range(k):v=cut(mul(v,p),n)
    return v

def series(a,b,c,d,e,f,G,n):
    one={(0,0):Q(1)};T={(0,1):Q(1)};g={}
    def shift(p):return {(m[0]+1,m[1]):v for m,v in p.items() if m[0]+1<=n}
    for _ in range(n+1):
        h=add({(0,0):a,(0,1):c,(0,2):f},scale(g,b),scale(mul(g,g),d),scale(mul(g,T),e));g=cut(shift(h),n)
    ht=add({(0,0):b,(0,1):e},scale(g,2*d));q=shift(ht);Dinv=one.copy();qk=one
    for j in range(1,n+1):qk=cut(mul(qk,q),n);Dinv=add(Dinv,qk)
    phase={}
    for (i,j),co in G.items():phase=add(phase,scale(shift(mul(powp(g,i,n),powp(T,j,n))),co))
    exp=one.copy();pk=one
    for j in range(1,n+1):pk=cut(mul(pk,phase),n);exp=add(exp,scale(pk,Q(1,factorial(j))))
    result=cut(mul(exp,Dinv),n);out=[Q(0)]*(n+1)
    for (m,k),co in result.items():
        if k%2==0:out[m]+=co*Q(factorial(k),2**(k//2)*factorial(k//2))
    return out

def wick(a,b,c,d,e,f,G,n):
    P={(0,1,0):a,(1,1,0):b,(0,1,1):c,(2,1,0):d,(1,1,1):e,(0,1,2):f};P=add(P,{(i,0,j):v for (i,j),v in G.items()});P={m:v for m,v in P.items() if v};v={(0,0,0):Q(1)};out=[]
    for _ in range(n):v=mul(v,P);out.append(gaussian(v))
    return out


def algebra():
    t,x,z,T,a,b,c,f,rho,sigma=S.symbols('t x z T a b c f rho sigma',nonzero=True);eta=f*(rho-sigma);e=f*(rho+sigma);d=f*rho*sigma
    u=(rho*c-b+1/t)/eta;v=(b-1/t-sigma*c)/eta;k=u*v-a/f
    U=T+rho*z+u;V=T+sigma*z+v;h=a+b*z+c*T+d*z*z+e*z*T+f*T*T;D=1-t*S.diff(h,z)
    assert S.cancel(h-z/t-f*(U*V-k))==0
    assert S.cancel(S.expand(t*t*k).subs(t,0)+1/eta**2)==0
    Z=(x-k/x-u+v)/(rho-sigma);Ts=x-u-rho*Z
    assert S.cancel((h-z/t).subs({z:Z,T:Ts},simultaneous=True))==0
    hz=S.diff(h,z).subs({z:Z,T:Ts},simultaneous=True);ht=S.diff(h,T).subs({z:Z,T:Ts},simultaneous=True)
    assert S.cancel(S.diff(Z,x)-ht/(eta*x))==0
    assert S.cancel(S.diff(Ts,x)-(1-t*hz)/(t*eta*x))==0
    assert S.cancel(D+rho*t*S.diff(h,T)-t*eta*U)==0
    # f=0 with the actual pole x=z+c/e retained.
    ee,dd=S.symbols('ee dd',nonzero=True);zz=x-c/ee;A=a+b*zz+dd*zz*zz;F=(zz/t-A)/(ee*x);hh=a+b*z+c*T+dd*z*z+ee*z*T
    assert S.cancel((hh-z/t).subs({z:zz,T:F},simultaneous=True))==0
    Ds=(1-t*S.diff(hh,z)).subs({z:zz,T:F},simultaneous=True)
    assert S.cancel(S.diff(F,x)-Ds/(t*ee*x))==0
    T0=(1/t-b)/ee;assert S.expand(1-t*(b+ee*(x+T0)))==-ee*t*x
    # The formal derivative flux identity holds for an arbitrary coordinate X.
    X=S.Function('X')(t,T);kap=S.symbols('kap',nonzero=True);w=S.diff(X,T)/(t*kap*X);R=S.diff(X,t)/S.diff(X,T)
    assert S.simplify(S.diff(w,t)+w/t-S.diff(w*R,T))==0
    # Null rotation removes the rank-one quadratic cross coefficient.
    E,Fc,Db=S.symbols('E Fc Db');lam=-E/(2*Fc)
    assert S.cancel(E+2*Fc*lam)==0
    assert S.cancel((Db+E*lam+Fc*lam**2)-(4*Db*Fc-E*E)/(4*Fc))==0
    return {'smooth_conic':'PASS: constraint, nonzero k, logarithmic Jacobian and formal X_T identity',
            'WZT_pole_face':'PASS: c!=0 or a!=0; x a Laurent unit',
            'trapped_component':'PASS: g=0 retained separately, D=-e*t*x, x=-1/(e*t)+O(1)',
            'formal_flux_identity':'PASS for arbitrary X(t,T)',
            'rank_one_null_rotation':'PASS'}


def connection_check(phi,t,x,m,label):
    K=S.QQ.frac_field(t);phi=S.expand(phi)
    exprs={int(term.as_powers_dict().get(x,0)):None for term in S.Add.make_args(phi)}
    low=min(exprs);high=max(exprs);q=max(0,-low);p=max(0,high)
    if p and q:basis=list(range(-q,p))
    elif p:basis=list(range(p))
    elif q:basis=list(range(0,-q,-1))
    else:basis=[0]
    def parse(expr):
        expr=S.expand(expr);powers=[int(term.as_powers_dict().get(x,0)) for term in S.Add.make_args(expr)];offset=max(0,-min(powers))
        poly=S.Poly(S.expand(expr*x**offset),x,domain=K)
        return {mon[0]-offset:co for mon,co in poly.rep.to_dict().items() if co}
    coeff=parse(phi);L={j:K.convert(j)*v for j,v in coeff.items() if j and v};valid=set(basis)
    def combine(out,j,value):
        out[j]=out.get(j,K.zero)+value
        if not out[j]:del out[j]
    def exact(qpoly):
        out={}
        for j,value in qpoly.items():
            combine(out,j,K.convert(j)*value)
            for k,co in L.items():combine(out,j+k,co*value)
        return out
    def reduce(poly):
        out=poly.copy();witness={}
        while any(j not in valid for j in out):
            j=max(out) if max(out)>max(valid) else min(out)
            if j>max(valid) and p:shift=j-p;lc=L[p]
            elif j<min(valid) and q:shift=j+q;lc=L[-q]
            else:shift=j;lc=K.convert(j);assert j!=0
            factor=out[j]/lc;combine(witness,shift,factor)
            for k,co in exact({shift:factor}).items():combine(out,k,-co)
        difference=poly.copy()
        for j,value in out.items():combine(difference,j,-value)
        for j,value in exact(witness).items():combine(difference,j,-value)
        assert not difference
        return out
    H=t**(m+1)
    if p:H*=S.cancel(t**m*K.to_sympy(L[p]))**p
    if q:H*=S.cancel(t**m*K.to_sympy(L[-q]))**q
    H=S.expand(H);rank=len(basis);bound=(m+1)*(rank+1);maxdeg=0
    for j in basis:
        remainder=reduce(parse(x**j*(S.diff(phi,t)-1/t)))
        for value in remainder.values():
            num=S.cancel(H*K.to_sympy(value));assert S.denom(num)==1
            deg=S.degree(num,t)
            if deg!=-S.oo:maxdeg=max(maxdeg,int(deg));assert deg<=bound
    assert S.degree(H,t)<=bound
    return {'case':label,'Laurent_orders':[p,q],'rank':rank,'H_degree':int(S.degree(H,t)),'max_numerator_degree':maxdeg,'bound':bound,'all_exact_reduction_witnesses':'PASS'}


def degree_support_checks():
    def product(A,B):return {(i+k,j+l) for i,j in A for k,l in B}
    records=[]
    for degree in [3,4,5]:
      for pole_order,m in [(2,2*degree-1),(1,degree-1)]:
          support={(i,j) for i in range(-pole_order,1) for j in [-1,0,1]}
          phases=set()
          for iz in range(degree+1):
              for jt in range(degree+1-iz):
                  terms={(0,0)}
                  for _ in range(iz+jt):terms=product(terms,support)
                  phases.update((i+1,j) for i,j in terms)
          phases.update(product(support,support))
          assert all(-m<=i<=1 and -degree<=j<=degree for i,j in phases)
          records.append({'coordinate_coefficient_t_pole_order':pole_order,'phase_t_powers':[-m,1],
                          'phase_x_orders':[-degree,degree],'G_degree':degree,'all_G_monomials':'PASS','cleared_coefficient_degree':m+1})
    return records


def main():
    start=time.monotonic();root=Path(__file__).parent;proof=algebra();t,x=S.symbols('t x')
    models=[('two_poles_m5',t*(x**3+x**-3)+x*x/t**2+x**-2/t**5,5),('two_poles_m2',t*(x**3+x**-3)+x*x/t+x**-2/t**2,2),
            ('regular_zero',t*x**3+x/t**2,2),('regular_infinity',t/x**3+1/(x*t**2),2),('constant_phase',t+1/t,2)]
    matrices=[connection_check(phi,t,x,m,label) for label,phi,m in models]
    G={(0,0):Q(-1),(0,1):Q(2),(0,2):Q(-1),(0,3):Q(1),(1,0):Q(2),(1,1):Q(-2),(1,2):Q(1),(2,0):Q(1),(2,1):Q(-1),(3,0):Q(1)}
    cases=[('f_nonzero',(1,2,3,-1,0,1),G),('f_nonzero_a_c_zero',(0,1,0,2,3,1),G),
           ('WZT_c_nonzero',(1,-1,2,1,1,0),G),('WZT_c_zero_a_nonzero',(2,3,0,0,1,0),G),
           ('trapped_component',(0,2,0,3,1,0),G),('trapped_empty_G',(0,0,0,0,1,0),{}),
           ('excluded_boundary_control',(0,0,0,0,0,1),{})]
    samples=[]
    for label,params,g in cases:
        av,bv,cv,dv,ev,fv=map(Q,params);mom=wick(av,bv,cv,dv,ev,fv,g,8);srs=series(av,bv,cv,dv,ev,fv,g,8)
        assert all(mom[m-1]==factorial(m)*srs[m] for m in range(1,9))
        disc=ev*ev-4*dv*fv
        if label=='excluded_boundary_control':assert disc==0 and not any(mom)
        else:assert disc!=0 and any(mom)
        samples.append({'case':label,'discriminant':str(disc),'moments_through_8':list(map(str,mom))})
    gates=[{'r':r,'WZT_prefix':3*(r+1)*r*(r+1)//2+r,'general_conic_prefix':6*(r+1)*r*(r+1)//2+r} for r in range(1,7)]
    assert gates[-1]['WZT_prefix']==447 and gates[-1]['general_conic_prefix']==888
    inputs=['verify_moment_transfers.py','verify_mixed_period_infinity.py']
    rec={'scope':'Universal conic and formal-flux algebra, exact Laurent-module witnesses, representative connection bounds and Wick/GL1 regressions; logarithmic growth and comparison are written analytic proofs.',
         'universal_phase_support':degree_support_checks(),'algebra':proof,'connection_models':matrices,'Wick_GL1_cases':samples,'all_degree_prefixes':[{'D':D,'WZT':D*D*(2*D+1)**2+2*D,'general':2*D*D*(2*D+1)**2+2*D} for D in [3,4,5]],'finite_degree_table':gates,
         'proof_boundary':'No computed contour estimate, external review, formalization or global ternary frame theorem.',
         'sympy_version':S.__version__,'elapsed_seconds':time.monotonic()-start,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'input_hashes':{n:hashlib.sha256((root/n).read_bytes()).hexdigest() for n in inputs}}
    (root/'quadratic_w_log_verification.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))


if __name__=='__main__':main()
