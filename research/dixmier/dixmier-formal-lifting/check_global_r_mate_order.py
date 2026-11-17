#!/usr/bin/env python3
"""Global first-order R-mate obstruction and admission facts for the nonsymmetric perturbation."""
from __future__ import annotations
import argparse
from functools import lru_cache
from math import comb
from pathlib import Path
import hashlib,json,platform,tempfile
import sympy as s
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    X,Q,h,r,v,H,lam,nu=s.symbols('X Q h r v H lambda nu',nonzero=True)
    def clean(op):return {i:s.factor(c) for i,c in op.items() if s.factor(c)!=0}
    def add(*ops):
        out={}
        for op in ops:
            for i,c in op.items():out[i]=out.get(i,0)+c
        return clean(out)
    def scale(op,c):return clean({i:c*a for i,a in op.items()})
    def same(a,b):return add(a,scale(b,-1))=={}
    def gm(a,b):
        out={}
        for (i,j),ca in a.items():
            for (k,l),cb in b.items():
                for u in range(i+1):
                    for w in range(j+1):
                        key=(i+k-u,j+l-w)
                        out[key]=out.get(key,0)+ca*comb(i,u)*comb(j,w)*h**(u+w)*s.diff(cb,X,u,Q,w)
        return clean(out)
    def gc(a,b):return add(gm(a,b),scale(gm(b,a),-1))
    pos=lambda c:{(0,0):c}
    R=X*(2-3*X*Q)
    B={(1,0):3*X*X,(0,1):2*(1-3*X*Q)}
    D0={(1,0):(1+3*X*Q)/2,(0,1):-3*Q*Q}
    assert gc(D0,pos(R))==pos(h)
    A0=add(pos(Q+3*X*Q*Q),scale(gm(pos(X),B),-s.Rational(1,3)))
    C=add(scale(A0,R),pos(-s.Rational(2,3)))
    perturbation=scale(gm(pos(X*X),C),nu*h)
    assert gc(C,pos(R))=={} and gc(perturbation,pos(R))=={}
    assert max(i+j for i,j in perturbation)==1
    assert all(s.Poly(a,X,Q,h).is_multivariate for a in perturbation.values())
    @lru_cache(None)
    def delta(f,k):
        if k==0:return f
        return s.factor(s.diff(delta(f,k-1),v)/v)
    def mul(a,b):
        out={}
        for i,ca in a.items():
            for j,cb in b.items():
                for k in range(i+1):
                    n=i+j-k
                    out[n]=out.get(n,0)+ca*comb(i,k)*h**k*delta(cb,k)
        return clean(out)
    def power(op,n):
        ans={0:s.Integer(1)}
        for _ in range(n):ans=mul(ans,op)
        return ans
    A={1:s.Integer(1)}
    centC={1:r,0:-r*v}
    W=scale(power(centC,2),v**-3)
    qfull=add(A,scale(W,lam*h))
    Slambda=add(scale(mul(qfull,qfull),-3*r),{0:2*v+4/(3*r)})
    Snew=add(Slambda,scale(centC,nu*h*v**-2))
    p4=-3*lam**2*h*h*r**5*v**-6
    assert Slambda[4]==p4 and Snew[4]==p4
    assert same(clean({i:a.subs(h,0) for i,a in Snew.items()}),{2:-3*r,0:2*v+4/(3*r)})
    f=5*v*v/(14*r)
    coefficient=s.diff(p4,r)+f*delta(p4,1)-4*p4*delta(f,1)
    assert s.factor(coefficient)==0
    assert s.factor(s.diff(p4,r)/r**4)==-15*lam**2*h*h*v**-6
    assert s.denom(s.cancel(f)).has(r)
    fhom=v**(-s.Rational(3,2))
    assert s.factor(-6*fhom/v**2-4*s.diff(fhom,v)/v)==0
    # Generic first-order commutator formula, checked with a free rational test coefficient.
    ftest=r*r*v**3+r*v**-2
    qtest=r*v*v+h/v
    first={1:ftest,0:qtest}
    comm=add(mul(first,Snew),scale(mul(Snew,first),-1))
    assert s.factor(comm[4]/h-(ftest*delta(p4,1)-4*p4*delta(ftest,1)))==0
    # Critical face extracted from the actual admitted deformation's weighted support.
    critical={}
    for j,a in Snew.items():
        if j==0:a=a-4/(3*r)
        for term in s.Add.make_args(s.expand(a)):
            powers=term.as_powers_dict()
            hn=int(powers.get(h,0));vm=int(powers.get(v,0))
            weight=5*hn+j+2*vm-2
            assert weight>=0
            if weight==0:
                critical[j]=critical.get(j,0)+term/h**hn
    critical=clean(critical)
    beta,sigma,g,mu=s.symbols('beta sigma g mu',nonzero=True)
    q={2:beta*v**-3,1:s.Integer(1)}
    base=add(scale(mul(q,q),g),{0:mu*v})
    crit=add(base,{1:sigma*v**-2})
    coupled={beta:lam*r*r,sigma:nu*r,g:-3*r,mu:s.Integer(2),h:s.Integer(1)}
    assert same(critical,clean({i:a.subs(coupled) for i,a in crit.items()}))
    xi,t=s.symbols('xi t')
    weighted=[];canonical=[]
    Cj=lambda j:s.prod(H-2*l*h for l in range(j))
    for j,a in crit.items():
        for term in s.Add.make_args(s.expand(a)):
            vm=int(term.as_powers_dict().get(v,0))
            weighted.append((vm+3*j,term*xi**j))
            canonical.append((vm-2*j,s.factor(term/v**vm)*Cj(j)))
    assert max(a for a,b in weighted)==6
    top=s.factor(sum(b for a,b in weighted if a==6))
    assert s.factor(top-g*(xi+beta*v**-3*xi**2)**2)==0
    assert s.expand(top).coeff(xi,3)==2*beta*g*v**-3
    assert min(a for a,b in canonical)==-14
    bottom=s.factor(sum(b for a,b in canonical if a==-14))
    expected=g*beta*beta*H*(H-2*h)*(H-7*h)*(H-9*h)
    assert s.factor(bottom-expected)==0
    assert s.rem((1+t*t)*(1+t**7),s.cyclotomic_poly(14,t),t)==0
    # Only the q-fixing adjoint is checked; this is not a mate obstruction for sigma !=0.
    adjoint_shift=-3*h*v**-2-v**3/beta
    Astar={1:-s.Integer(1),0:adjoint_shift}
    def dagger(op):
        ans={}
        for j,a in op.items():ans=add(ans,mul(power(Astar,j),{0:a}))
        return ans
    assert same(dagger(q),q) and same(dagger(base),base)
    diff=add(dagger(crit),scale(crit,-1))
    assert diff.get(1)==-2*sigma*v**-2 and diff!={}
    # Negative control: forgetting the r denominator fails the leading equation.
    wrong=5*v*v/14
    bad=s.factor(s.diff(p4,r)+wrong*delta(p4,1)-4*p4*delta(wrong,1))
    assert bad!=0
    output=args.output or Path(tempfile.mkdtemp(prefix='global-rmate-'))/'global_r_mate_order_controls.json'
    output.parent.mkdir(parents=True,exist_ok=True)
    encode=lambda op:{str(i):str(a) for i,a in op.items()}
    receipt={
        'scope':'No globally admitted first-order R-mate for mixed S_lambda and S_lambda+nu*h*X^2*C; new critical admission and symbol facts only',
        'status':'PASS','environment':{'python':platform.python_version(),'sympy':s.__version__},
        'domain':'D in A2(K((h))) of original momentum order <=1; f,q in K((h))[r,v,v^-1], no r denominators',
        'leading_coefficient':str(p4),'necessary_rational_f':str(f),
        'r_divisor_derivative_coefficient':str(s.factor(s.diff(p4,r)/r**4)),
        'new_global_perturbation':'nu*h*X^2*C; globally admitted first-order, classical reduction zero',
        'critical_coupling':{'g':'-3r','beta':'lambda*r^2','mu':'2','sigma':'nu*r','eta':'1'},
        'critical_operator':encode(crit),'top_weight':6,'top_symbol':str(top),'bottom_weight':-14,
        'bottom_polynomial':str(bottom),'primitive14_character_zero':True,
        'q_adjoint_shift':str(adjoint_shift),'newterm_adjoint_difference':encode(diff),
        'negative_control':'Removing required r denominator from f makes top commutator nonzero',
        'hashes':{name:sha(HERE/name) for name in ['GLOBAL_R_MATE_ORDER_OBSTRUCTION.md','check_global_r_mate_order.py',
                'EXACT_POSITION_CENTRALIZER_AND_COMPLETION.md','GLOBAL_ADMISSION_AND_QUARTIC_CRITICAL_FACE.md',
                'ORDERED_QUARTIC_INFINITY_OBSTRUCTION.md']},
        'coverage_boundary':'Exact finite identities; universal polynomial-r no-go is written divisor proof. No higher-order R-mate claim, finite scalar mate calculation, or bounded Long lift.',
        'independent_replay':False,'formal_verification':False,'external_review':False}
    output.write_text(json.dumps(receipt,indent=2)+'\n')
    print('Global first-order R-mate controls PASS')
    print('Verification output:',output)
if __name__=='__main__':main()
