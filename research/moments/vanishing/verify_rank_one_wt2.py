#!/usr/bin/env python3
"""Exact algebra for the rank-one WT2 classification; analytic proof is written."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import json,time,hashlib
import sympy as S
from verify_quadratic_w_log_periods import wick,series
from verify_moment_transfers import mul,gaussian


def algebra_checks():
    t,z,T,a,b,c,f=S.symbols('t z T a b c f',nonzero=True);q=1-b*t;A=a+c*T+f*T*T;g=t*A/q
    h=A+b*z;assert S.cancel(g-t*h.subs(z,g))==0
    coeff=S.symbols('r0 u0 v0 s r1 u1 v1 r2 u2 r3');r0,u0,v0,s,r1,u1,v1,r2,u2,r3=coeff
    G=r0+u0*T+v0*T*T+s*T**3+r1*z+u1*z*T+v1*z*T*T+r2*z*z+u2*z*z*T+r3*z**3
    L=s*T/f+v1*z/f+v0/f-s*c/f**2;R=S.expand(G-(b*z+A)*L)
    assert S.degree(R,T)<=1 and S.degree(R.subs(T,0),z)<=3 and S.degree(S.diff(R,T),z)<=2
    assert S.expand(G-R-(b*z+A)*L)==0
    degrees=[]
    for mon in [1,T,z,z*T,z*z,z*z*T,z**3]:degrees.append(int(S.degree(S.sympify(mon).subs(z,-A/b),T)))
    assert degrees==list(range(7))
    U=0
    for i in range(4):
        for j in range(4-i):U+=S.symbols('c'+str(i)+str(j))*t**(i+1)*A**i*T**j*q**(3-i)
    U=S.expand(U-q**3*T*T/2);assert S.degree(U,t)<=4 and S.degree(U,T)<=6
    source=S.expand(S.diff(U,t)*q+3*b*U+b*q**3);assert S.degree(source,t)<=4
    assert S.cancel(S.diff(U/q**3,t)+b/q-source/q**4)==0
    # Flattened Gaussian: exact phase coefficient and logarithmic-derivative obstruction.
    rr=S.symbols('rr');gamma=-b/(2*f);Phi=S.cancel(t*rr+gamma*g-T*T/2)
    K=S.diff(Phi,T,2)/2;C=S.diff(Phi,T).subs(T,0);R0=Phi.subs(T,0)
    assert S.cancel(K+1/(2*q))==0
    H=t*rr-b*a*t/(2*f*q)+b*b*c*c*t*t/(8*f*f*q)
    assert S.cancel(R0+q*C*C/2-H)==0
    logderivative=S.diff(H,t)+b/(2*q)
    assert S.residue(logderivative,t,1/b)==-S.Rational(1,2)
    # Highest-t-degree chart table is a universal monomial-degree calculation.
    degree_cases=[]
    for j,dgj in [(3,0),(2,0),(2,1),(1,0),(1,1),(1,2)]:
        leading_degree=2*j+dgj;lower_max=max([2]+[2*k+3-k for k in range(j)])
        escape=leading_degree<lower_max
        assert escape==(j==1 and dgj==0)
        degree_cases.append({'highest_Z_power':j,'coefficient_T_degree':dgj,'leading_T_degree':leading_degree,'lower_T_degree_bound':lower_max,'possible_escape':escape})
    return {'division':'PASS: quotient linear, remainder R0(Z)+T*R1(Z)',
            'restriction_leading_degrees':degrees,'q3_phase_numerator_degree':4,'q4_source_numerator_degree':4,
            'flattened_quadratic_coefficient':'-1/(2q)','flattened_log_derivative_residue':'-1/2',
            'b_zero_degree_cases':degree_cases}


def frame_criterion_check():
    z,w,T,vz,vw,vT=S.symbols('z w T vz vw vT');xs=[z,w,T];vs=[vz,vw,vT]
    pieces={d:0 for d in range(4)}
    for i in range(4):
      for j in range(4-i):
       for k in range(4-i-j):
        d=i+j+k;pieces[d]+=S.symbols('p'+str(i)+str(j)+str(k))*z**i*w**j*T**k
    P=sum(pieces.values());direction=lambda expr:sum(v*S.diff(expr,x) for x,v in zip(xs,vs))
    substitution=dict(zip(xs,vs))
    expected=2*sum(S.diff(pieces[3],x).subs(substitution,simultaneous=True)*x for x in xs)+2*pieces[2].subs(substitution,simultaneous=True)
    assert S.expand(direction(direction(P))-expected)==0
    # In an actual circular W-linear frame the W direction meets every equation.
    generic=S.expand(P).subs({S.symbols('p030'):0,S.symbols('p020'):0,S.symbols('p120'):0,S.symbols('p021'):0})
    assert S.diff(generic,w,2)==0
    return 'PASS: universal D_v^2 identity and the four omitted W-power coefficients'


def escape_checks():
    x,v,k=S.symbols('x v k',nonzero=True);cases=[(2,1,2,3),(1,0,1,2),(1,-1,1,2)]
    cases.extend((1,0,n,m) for m in [2,3,4] for n in range(1,m))
    rows=[]
    for p,r,n,m in sorted(set(cases)):
        F=v*x**n+k*x**m;rel=x**(m-n)+S.Rational(n,m)*v/k
        action=(1-S.Rational(n,m))*v*x**n;hessian=n*(n-m)*v*x**(n-2)
        for expr in [S.diff(F,x),F-action,S.diff(F,x,2)-hessian]:
            num,den=S.fraction(S.cancel(expr));assert S.rem(num,rel,x)==0 and S.gcd(den,rel,x)==1
        assert S.gcd(S.fraction(action)[0],rel,x)==1
        assert S.gcd(rel,S.diff(rel,x),x)==1
        beta=S.Rational(p-r,m-n);power=p+n*beta
        rows.append({'p_r_n_m':[p,r,n,m],'escape_scale':str(beta),'highest_action_power':str(power),'Hessian_nonzero':True,'bounded_plus_escape_dimension':(n-1)+(m-n)})
    return rows


def matrix_check(Phi,b,t,T,label):
    K=S.QQ.frac_field(t);q=1-b*t;Pz=S.Poly(S.diff(Phi,T),T,domain=K);d=Pz.degree()+1
    A=S.cancel(q**3*Pz.LC());assert S.denom(A)==1 and S.degree(A,t)<=4
    H=S.expand(q**4*A**d);bound=4*d+4;maxdeg=0
    for j in range(d-1):
        source=S.Poly(T**j*(S.diff(Phi,t)+b/q),T,domain=K);out=source;witness=S.Poly(0,T,domain=K)
        while out.degree()>=d-1:
            m=out.degree()-(d-1);factor=out.LC()/Pz.LC();quot=S.Poly(T**m,T,domain=K).mul_ground(factor)
            exact=quot.diff()+Pz*quot;out-=exact;witness+=quot
        assert S.cancel((source-out-witness.diff()-Pz*witness).as_expr())==0
        for co in out.all_coeffs():
            numerator=S.cancel(H*co);assert S.denom(numerator)==1
            deg=S.degree(numerator,t)
            if deg!=-S.oo:maxdeg=max(maxdeg,int(deg));assert deg<=bound
    assert S.degree(H,t)<=bound
    return {'case':label,'phase_degree':d,'H_degree':int(S.degree(H,t)),'max_N_degree':maxdeg,'bound':bound,'exact_reduction_witnesses':'PASS'}


def main():
    start=time.monotonic();p=Path(__file__).parent;algebra=algebra_checks();escapes=escape_checks();t,T=S.symbols('t T');q=1-t;gg=t*T*T/q
    models=[('b0_Z3',t**4*T**6-T*T/2,0),('b0_Z2',t**3*T**4-T*T/2,0),
            ('b0_exceptional_escape',t*t*T*T+t*T**3-T*T/2,0),('b_nonzero_alpha',gg*gg-T*T/2,1),
            ('b_nonzero_beta',gg*T-T*T/2,1),('flattened_constant',-gg/2-T*T/2,1),
            ('flattened_linear',t*T-gg/2-T*T/2,1),('safe_null_boundary',-T*T/2,0)]
    matrices=[matrix_check(S.cancel(Phi),b,t,T,label) for label,Phi,b in models]
    dense={(0,0):Q(1),(0,1):Q(-2),(0,2):Q(1),(0,3):Q(1),(1,0):Q(-1),(1,1):Q(2),(1,2):Q(1),(2,0):Q(-1),(2,1):Q(1),(3,0):Q(1)}
    cases=[('b0_dense',(1,0,2,1),dense),('b1_dense',(1,1,-1,2),dense),('b0_escape',(0,0,0,1),{(1,0):Q(1),(0,3):Q(1)}),
           ('flattened_constant',(1,2,-2,3),{(1,0):Q(-2,3),(0,0):Q(1,2)-Q(1,3),(0,1):Q(2,3),(0,2):Q(-1)}),
           ('flattened_linear',(0,1,0,1),{(0,1):Q(1),(1,0):Q(-1,2),(0,2):Q(-1,2)}),
           ('safe_WT2',(0,0,0,1),{}),('safe_general_product',(2,0,-3,2),{})]
    samples=[]
    for label,params,G in cases:
        a,b,c,f=map(Q,params);mom=wick(a,b,c,Q(0),Q(0),f,G,8);srs=series(a,b,c,Q(0),Q(0),f,G,8)
        assert all(mom[m-1]==factorial(m)*srs[m] for m in range(1,9))
        if b==0 and not G:assert not any(mom)
        else:assert any(mom)
        samples.append({'case':label,'moments_through_8':list(map(str,mom))})
    # Independent closed Gaussian series check in the cancellation stratum.
    av,bv,cv,fv,rv=map(S.Rational,[1,2,-2,3,S.Rational(1,2)]);qv=1-bv*t
    Hv=t*rv-bv*av*t/(2*fv*qv)+bv*bv*cv*cv*t*t/(8*fv*fv*qv)
    exact=S.series(S.exp(Hv)/S.sqrt(qv),t,0,9).removeO().expand()
    G={(0,0):Q(1,2)-Q(1,3),(1,0):Q(-2,3),(0,1):Q(2,3),(0,2):Q(-1)}
    mom=wick(Q(1),Q(2),Q(-2),Q(0),Q(0),Q(3),G,8)
    assert all(S.Rational(mom[m-1].numerator,mom[m-1].denominator)==factorial(m)*exact.coeff(t,m) for m in range(1,9))
    # Mixed circular-weight boundary and sharpness in P=W*T^2.
    P={(0,1,2):Q(1)};mixed=0
    for i,j,k in [(0,0,0),(1,0,1),(2,1,0),(3,0,2),(1,2,1)]:
        v={(i,j,k):Q(1)}
        for m in range(1,i+j+k+3):
            v=mul(v,P)
            if i!=m+j:assert gaussian(v)==0
            if m>i+j+k:assert gaussian(v)==0
            mixed+=1
    for n in range(1,6):
        value=gaussian({(n,n,2*n):Q(1)});assert value==Q(factorial(n)*factorial(2*n),2**n*factorial(n))
    gates=[{'d':d,'B':4*d+4,'prefix':(4*d+4)*(d-1)*d//2+d-1} for d in range(2,7)];assert gates[-1]['prefix']==425
    inputs=['verify_quadratic_w_log_periods.py','verify_moment_transfers.py','verify_mixed_period_infinity.py']
    rec={'scope':'Rank-one WT2 algebra and finite regressions supporting the written all-W-linear cubic theorem; analytic proof is not computed.',
         'frame_criterion':frame_criterion_check(),'algebra':algebra,'escape_models':escapes,'connection_models':matrices,'Wick_GL1_cases':samples,
         'flattened_closed_Gaussian_series':'PASS through 8','mixed_weight_checks':mixed,'sharp_product_cases':5,'rank_one_gate_table':gates,
         'full_W_linear_pure_gate':888,'full_W_linear_mixed_cutoff':'3*deg Q, sharp by existing shear theorem',
         'proof_boundary':'Written shrinking-period observation and imported Hien--Roucairol comparison; no external review, formal verification or global frame normalization.',
         'sympy_version':S.__version__,'elapsed_seconds':time.monotonic()-start,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'input_hashes':{n:hashlib.sha256((p/n).read_bytes()).hexdigest() for n in inputs}}
    (p/'rank_one_wt2_verification.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))


if __name__=='__main__':main()
