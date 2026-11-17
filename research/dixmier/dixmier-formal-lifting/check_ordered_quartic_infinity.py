#!/usr/bin/env python3
"""Exact finite controls for the all-order ordered-quartic infinity proof."""
from __future__ import annotations

import argparse
from pathlib import Path
import hashlib
import json
import tempfile

import sympy as sp

HERE=Path(__file__).resolve().parent
u,w,g,mu,s,eta=sp.symbols('u w g mu s eta')


def clean(op):
    return {key:sp.expand(c) for key,c in op.items() if sp.expand(c)!=0}


def add(left,right,factor=1):
    return clean({key:left.get(key,0)+factor*right.get(key,0)
                  for key in set(left)|set(right)})


def scale(op,c):
    return clean({key:c*value for key,value in op.items()})


def mul(left,right,minimum=0,quantum=eta):
    out={}
    for (i,e),a in left.items():
        for (j,f),b in right.items():
            falling=sp.Integer(1)
            for k in range(max(-1,i+j-minimum)+1):
                if i>=0 and k>i:
                    break
                key=(i+j-k,e+f-2*k)
                out[key]=out.get(key,0)+sp.binomial(i,k)*falling*quantum**k*a*b
                falling*=f-2*k
                if falling==0:
                    break
    return clean(out)


def power(op,n,quantum=eta):
    out={(0,0):sp.Integer(1)}
    for _ in range(n):
        out=mul(out,op,quantum=quantum)
    return out


def comm(left,right,quantum=eta):
    return add(mul(left,right,quantum=quantum),
               mul(right,left,quantum=quantum),-1)


def build(ordering_u,ordering_w,G=g,Mu=mu,S=s,Eta=eta):
    L={(2,-3):S,(1,-5):(2-ordering_u)*S*Eta,
       (0,-7):ordering_w*S*Eta**2}
    L=clean(L)
    P=add(add(mul(L,L,quantum=Eta),{(2,0):G}),{(0,1):Mu})
    return L,P


def half_positive(L,P,n,S=s,G=g,Eta=eta):
    assert n%2==1
    R=add(L,{(0,3):G/(2*S)})
    K=2*n-2
    for k in range(1,K+1):
        defect=add(mul(R,R,minimum=2-k,quantum=Eta),
                   {key:c for key,c in P.items() if key[0]>=2-k},-1)
        correction={(-k,e+3):-c/(2*S)
                    for (i,e),c in defect.items() if i==2-k}
        R=add(R,correction)
    defect=add(mul(R,R,minimum=2-K,quantum=Eta),
               {key:c for key,c in P.items() if key[0]>=2-K},-1)
    assert not defect,defect
    B=mul(power(P,(n-1)//2,quantum=Eta),R,quantum=Eta)
    return B


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    controls=[]
    cases=[('all_ordering_M2',1,u,w,g,mu,s,eta),
           ('all_ordering_M6',3,u,w,g,mu,s,eta),
           ('u5_w3_M10',5,5,3,g,mu,sp.Integer(1),sp.Integer(1))]
    for label,n,U,W,G,Mu,S,Eta in cases:
        L,P=build(U,W,G,Mu,S,Eta)
        B=half_positive(L,P,n,S,G,Eta)
        alpha=sp.binomial(sp.Rational(n,2),n)
        assert sp.expand(B.get((0,3*n),0)-alpha*G**n/S**n)==0
        for (j,e),c in B.items():
            assert e+3*j<=3*n
            assert (3*n-e-3*j)%5==0
            if j%2:
                assert e+3*j<=3*n-5
        C=comm(P,B,quantum=Eta)
        assert set(i for i,_ in C)<={0,1,2}
        target=sp.factor(C.get((1,3*n-2),0))
        expected=6*n*Eta*alpha*G**(n+1)/S**n
        assert sp.expand(target-expected)==0 and target!=0
        assert max(e for (j,e),c in C.items() if j==1)==3*n-2
        if n==3:
            assert sp.expand(target.subs({u:5,w:3})-expected)==0
            assert sp.expand(target.subs(w,(u*u-16)/4)-expected)==0
            bad=dict(B);bad[(0,9)]+=G**3/(16*S**3)
            assert sp.expand(comm(P,bad,quantum=Eta).get((1,7),0)-expected)!=0
        controls.append({'case':label,'partner_order':2*n,
                         'ordering_u':str(U),'ordering_w':str(W),
                         'binomial':str(alpha),'target':'V^'+str(3*n-2)+'*A',
                         'target_coefficient':str(target),'positive_part_terms':len(B),
                         'commutator_terms':len(C),'degree_bounds_pass':True})

    # Centralizer enlargement when the required g is zero.
    L0,P0=build(u,w,sp.Integer(0),sp.Integer(0))
    assert comm(P0,L0)=={}
    assert max(j for j,_ in L0)==2 and max(j for j,_ in P0)==4

    # A strict rational infinity-cone example beyond the Laurent family.
    V=sp.symbols('V',nonzero=True)
    delta=lambda f:sp.diff(f,V)/V
    def rational_mul(left,right):
        out={}
        for i,a in left.items():
            for j,b in right.items():
                for k in range(i+1):
                    out[i+j-k]=out.get(i+j-k,0)+sp.binomial(i,k)*a*b
                    b=delta(b)
        return {i:sp.cancel(c) for i,c in out.items() if sp.cancel(c)!=0}
    def rational_add(left,right,factor=1):
        return {i:sp.cancel(left.get(i,0)+factor*right.get(i,0))
                for i in set(left)|set(right)
                if sp.cancel(left.get(i,0)+factor*right.get(i,0))!=0}
    rho=(1+1/V)*V**-3
    Pr={4:rho**2,3:V**-4,2:1+V**-1,1:V**2,0:V**5}
    Rr={2:rho}
    for j in [1,0]:
        partial=rational_mul(Rr,Rr)
        Rr[j]=sp.cancel((Pr[j+2]-partial.get(j+2,0))/(2*rho))
    Cr=rational_add(rational_mul(Pr,Rr),rational_mul(Rr,Pr),-1)
    assert sp.limit(Cr[1]/V,V,sp.oo)==3
    nonsquare_p4=V**-6*(1+1/V)
    assert sp.limit(nonsquare_p4*V**6,V,sp.oo)==1
    assert sp.limit(nonsquare_p4/(V+1),V,-1)!=0

    # Exact global ordered-square coefficients and their infinity degrees.
    r,h,epsilon,tau=sp.symbols('r h epsilon tau')
    C={(1,0):r,(0,1):-r,(0,0):-sp.Rational(2,3)}
    X3={(0,-3):sp.Integer(1)}
    Wglobal=add(add(mul(X3,mul(C,C,quantum=h),quantum=h),
                   mul({(0,-5):(2-u)*h*r},C,quantum=h)),
                {(0,-7):w*h**2*r**2})
    Sb={(2,0):-3*r,(1,0):sp.Integer(4),(0,1):sp.Integer(2)}
    Sglobal=add(Sb,scale(mul(Wglobal,Wglobal,quantum=h),epsilon*h**2))
    leading={}
    for j in range(5):
        e=max(E for (J,E) in Sglobal if J==j)
        leading[j]=(e,sp.factor(Sglobal[(j,e)]))
    assert leading[4]==(-6,epsilon*h**2*r**4)
    assert leading[3]==(-5,-4*epsilon*h**2*r**4)
    assert leading[2]==(0,-3*r)
    assert leading[1]==(0,sp.Integer(4))
    assert leading[0]==(1,sp.Integer(2))
    remainders={2:Sb[(2,0)],1:Sb[(1,0)],0:Sb[(0,1)]}
    for j,removed in remainders.items():
        subtract_key=(j,0) if j else (0,1)
        rest=dict(Sglobal);rest[subtract_key]-=removed;rest=clean(rest)
        assert max(e for (J,e) in rest if J==j)=={2:-4,1:-3,0:-2}[j]
    centered={(1,0):sp.Integer(1),(0,1):-tau}
    Wscaled=scale(add(add(mul(X3,mul(centered,centered,quantum=1),quantum=1),
                          mul({(0,-5):2-u},centered,quantum=1)),{(0,-7):w}),r**2)
    critical=clean({key:c.subs(tau,0) for key,c in Wscaled.items()})
    expected=scale(build(u,w,sp.Integer(0),sp.Integer(0),sp.Integer(1),sp.Integer(1))[0],r**2)
    assert add(critical,expected,-1)=={}

    names=['ORDERED_QUARTIC_INFINITY_OBSTRUCTION.md','check_ordered_quartic_infinity.py',
           'GLOBAL_ADMISSION_AND_QUARTIC_CRITICAL_FACE.md']
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    receipt={'claim':'QOI1','status':'PASS','fractional_power_controls':controls,
             'exceptional_orderings_u5_and_discriminant16':True,
             'corrupted_M6_scalar_term_rejected':True,
             'g_mu_zero_centralizer_boundary':True,
             'strict_rational_infinity_cone_M2_coefficient':3,
             'nonsquare_leading_divisor_boundary':True,
             'global_ordered_square_leading_coefficients':
                 {str(j):{'V_degree':e,'coefficient':str(c)} for j,(e,c) in leading.items()},
             'global_ordered_square_remainder_degrees':{'p2':-4,'p1':-3,'p0':-2},
             'global_critical_rescaling_matches_L_F':True,
             'hashes':{name:sha(HERE/name) for name in names},
             'coverage_boundary':(
                 'Exact positive parts and commutators for orders2,6 with all '
                 'ordering parameters and order10 at u5,w3; weighted degree '
                 'bounds; a strict rational infinity example; a nonsquare '
                 'divisor control; and full global ordered-square coefficients. '
                 'All orders, polynomial centralizer exhaustion, trace descent '
                 'and the full infinity cone use the written weighted-symbol '
                 'proof, not extrapolation from these finite inputs.'),
             'independent_replay':False,'formal_verification':False,'external_review':False}
    output=args.output or Path(tempfile.mkdtemp(prefix='qoi1-check-'))/'ordered_quartic_infinity_controls.json'
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'claim':'QOI1','status':'PASS','receipt':str(output)},indent=2))


if __name__=='__main__':
    main()
