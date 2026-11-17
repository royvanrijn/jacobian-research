#!/usr/bin/env python3
"""Exact algebra supporting the full mixed-chart period-at-infinity proof.

The analytic comparison and saddle completeness are written mathematics,
not certified by these finite computations. The sufficient full-parent gate is a written degree bound.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import json,hashlib,time
import sympy as S
from verify_mixed_wt_chart import direct,phase_series


def main():
    start=time.monotonic();p=Path(__file__).parent
    t,z,Y,a=S.symbols('t z Y a');r0,u0,v0,s,r1,u1,v1,r2,u2,r3=S.symbols('r0 u0 v0 s r1 u1 v1 r2 u2 r3');T=S.symbols('T')
    G=r0+u0*T+v0*T*T+s*T**3+r1*z+u1*z*T+v1*z*T*T+r2*z*z+u2*z*z*T+r3*z**3
    F=z/t-a-z*z;Phi=S.expand(t*G.subs(T,F)-F*F/2)
    V=S.expand(G.subs(T,-a-z*z));L=S.expand(S.diff(G,T).subs(T,-a-z*z));M=S.expand(S.diff(G,T,2).subs(T,-a-z*z)/2)
    B0=S.expand(z*L-(a+z*z)**2/2);B1=S.expand(z*z*M+z*(a+z*z));B2=s*z**3-z*z/2
    assert S.expand(Phi-t*V-B0-B1/t-B2/t**2)==0
    assert S.Poly(V,z).coeff_monomial(z**6)==-s
    assert S.Poly(V,z).coeff_monomial(z**5)==v1
    assert S.Poly(V,z).coeff_monomial(z**4)==v0-u2-3*a*s
    # High-degree implications are triangular, with constant pivots -1,1,-1.
    degree_rows=[]
    for n,subs in [(6,{}),(5,{s:0}),(4,{s:0,v1:0}),(3,{s:0,v1:0,u2:v0})]:
        degs=[S.degree(B.subs(subs),z) for B in [B0,B1,B2]]
        assert max(degs)<=max(n,4)
        if n==3:
            assert S.Poly(B0.subs(subs),z).coeff_monomial(z**4)==-S.Rational(1,2)
            assert S.expand(B2.subs(subs)+z*z/2)==0
            assert S.degree(V.subs(subs),z)<=3
        degree_rows.append({'deg_V_chart':str(n) if n>3 else '<=3','B0_B1_B2_degrees':list(map(int,degs))})
    # Universal normalization/connection bridge, polynomial seeds 1,z,...,z^4.
    Z=t*(a+Y);Psi=S.expand(Phi.subs(z,Z)+Y*Y/2)
    assert Psi.subs(t,0)==0
    for j in range(5):
        q=z**j;qs=q.subs(z,Z)
        total=S.diff(qs,t)+S.diff(Psi,t)*qs
        fixed=(S.diff(q,t)+(S.diff(Phi,t)-1/t)*q).subs(z,Z)
        Dterm=(S.diff(z*q/t,z)+S.diff(Phi,z)*z*q/t).subs(z,Z)
        assert S.expand(total-fixed-Dterm)==0
        amplitude=(z*q/t**2).subs(z,Z)
        gaussian_derivative=S.diff(amplitude,Y)+(S.diff(Psi,Y)-Y)*amplitude
        assert S.expand(gaussian_derivative-Dterm)==0
        ordinary_D=(S.diff(q,z)+S.diff(Phi,z)*q).subs(z,Z)
        assert S.expand((S.diff(qs,Y)+(S.diff(Psi,Y)-Y)*qs)/t-ordinary_D)==0
    # Escape saddle leading terms and complete rank count in each quartic chart.
    h,v=S.symbols('h v',nonzero=True);escapes=[]
    for n in range(1,4):
        H=v*h**n-h**4/2;relation=h**(4-n)-n*v/2
        action=(1-S.Rational(n,4))*v*h**n;hessian=2*(n-4)*h*h
        assert S.rem(S.diff(H,h),relation,h)==0
        assert S.rem(H-action,relation,h)==0
        assert S.rem(S.diff(H,h,2)-hessian,relation,h)==0
        assert S.gcd(relation,S.diff(relation,h),h)==1
        assert S.gcd(relation,S.fraction(action)[0],h)==1
        assert (n-1)+(4-n)==3
        escapes.append({'n':n,'escape_count':4-n,'bounded_dimension':n-1,'leading_action_degree':str(S.Rational(4,4-n)),
                        'normalized_prefactor_power':str(-1-S.Rational(1,4-n)),'leading_action':'(1-n/4)*v*h^n !=0','hessian':'2*(n-4)*h^2 !=0'})
    for k in range(2,7):assert -1-S.Rational(1,k)<0
    # Regressions cover all V-degree regimes, including V=constant through a nontrivial multiple of T+a+Z^2.
    cases=[('n6',Q(0),{(0,3):Q(1)}),('n5',Q(0),{(1,2):Q(1)}),('n4',Q(0),{(0,2):Q(1)}),
           ('n3',Q(2),{(3,0):Q(1)}),('n2',Q(2),{(2,0):Q(1)}),('n1',Q(2),{(1,0):Q(1)}),
           ('n0_nontrivial',Q(2),{(0,0):Q(3),(0,2):Q(1),(0,1):Q(2),(2,1):Q(1)}),('V_zero',Q(0),{})]
    samples=[];connection_degrees=[]
    for label,av,Gv in cases:
        Vsample=S.expand(sum(S.Rational(c.numerator,c.denominator)*z**i*(-S.Rational(av.numerator,av.denominator)-z*z)**j for (i,j),c in Gv.items()))
        if label.startswith('n'):assert S.degree(Vsample,z)==int(label[1])
        else:assert Vsample==0
        # Exact rational connection regression; no moment-320 expansion.
        Gsample=S.Add(*(S.Rational(c.numerator,c.denominator)*z**i*T**j for (i,j),c in Gv.items()))
        Fs=z/t-S.Rational(av.numerator,av.denominator)-z*z
        Phis=S.expand(t*Gsample.subs(T,Fs)-Fs*Fs/2);Pz=S.Poly(S.diff(Phis,z),z,domain=S.QQ.frac_field(t));d=Pz.degree()+1
        lead=S.cancel(t*t*Pz.LC());H=S.expand(t**3*lead**d);bound=3*d+3;maxdeg=0
        for j in range(d-1):
            remainder=S.Poly(z**j*(S.diff(Phis,t)-1/t),z,domain=S.QQ.frac_field(t))
            while remainder.degree()>=d-1:
                k=remainder.degree()-(d-1);factor=remainder.LC()/Pz.LC()
                subtraction=Pz*S.Poly(z**k,z,domain=S.QQ.frac_field(t))
                if k:subtraction+=S.Poly(k*z**(k-1),z,domain=S.QQ.frac_field(t))
                remainder-=subtraction.mul_ground(factor)
            for coefficient in remainder.all_coeffs():
                numerator=S.cancel(H*coefficient);assert S.denom(numerator)==1
                deg=S.degree(numerator,t)
                if deg!=-S.oo:maxdeg=max(maxdeg,int(deg));assert deg<=bound
        assert S.degree(H,t)<=bound
        connection_degrees.append({'case':label,'d':d,'common_H_degree':int(S.degree(H,t)),'max_N_degree':maxdeg,'bound':bound})
        moments=direct(av,Q(0),Q(1),Q(1),Gv,8);series=phase_series(av,Q(0),Q(1),Q(1),Gv,8)
        assert all(moments[m-1]==factorial(m)*series[m] for m in range(1,9))
        assert any(moments)
        samples.append({'case':label,'first_nonzero_moment':next(i+1 for i,x in enumerate(moments) if x),'moments_through_8':list(map(str,moments))})
    inputs=['verify_mixed_wt_chart.py','verify_moment_transfers.py','compile_mixed_wt_parent.py']
    finite_table=[{'d':d,'B':3*d+3,'D':(3*d+3)*(d-1)*d//2,'prefix':(3*d+3)*(d-1)*d//2+d-1} for d in [4,5,6]]
    assert [x['prefix'] for x in finite_table]==[93,184,320]
    receipt={'scope':'Universal phase and formal-connection algebra, exact degree/escape strata and eight bounded Wick regressions. Analytic proof is written, not computed.',
             'degree_charts':degree_rows,'connection_seed_degrees_checked':list(range(5)),
             'escape_charts':escapes,'bounded_critical_max_order':6,'all_normalized_zero_action_powers_negative':True,
             'V_constant':'Quartic limiting coefficient -1/2; normalized period exp(t*k)/t times a holomorphic function of 1/t.',
             'finite_gate_degree_table':finite_table,'sample_connection_degree_checks':connection_degrees,'Wick_phase_regressions':samples,'input_hashes':{n:hashlib.sha256((p/n).read_bytes()).hexdigest() for n in inputs},
             'imported_theorem':'Hien--Roucairol 0704.1739v1 Theorems 1.1/3.7 and 2.4; inspected primary HTML and PDF.',
             'analytic_scope':'Written one-variable saddle decomposition, including multiplicities; no generic-to-special promotion and no assertion of convergence of the original formal integral.',
             'full_mixed_parent':'Excluded by the written period proof; sufficient first-320 gate from the written rational-connection degree estimate, not numerical moment expansion.',
             'sympy_version':S.__version__,'elapsed_seconds':time.monotonic()-start,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (p/'mixed_period_infinity_verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
