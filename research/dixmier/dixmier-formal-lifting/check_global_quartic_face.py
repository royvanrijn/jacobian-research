#!/usr/bin/env python3
"""Global exact-R admission and local negative-pole face controls; no mate construction."""
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
    X,Q,h,r,v,H,t,eps=s.symbols('X Q h r v H t epsilon')
    def clean(op):
        return {i:s.factor(c) for i,c in op.items() if s.factor(c)!=0}
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
    pos=lambda f:{(0,0):f}
    RR=X*(2-3*X*Q)
    B={(1,0):3*X*X,(0,1):2*(1-3*X*Q)}
    A0=add(pos(Q+3*X*Q*Q),scale(gm(pos(X),B),-s.Rational(1,3)))
    a=add(A0,pos(1/X))
    assert gc(B,pos(RR))=={} and gc(A0,pos(RR))=={}
    assert gc(a,pos(1/X))==pos(h*X)
    F=3/X-s.Rational(5,3)*RR/X**2+s.Rational(1,3)*RR**2/X**3
    assert s.factor(F-(1/X+Q+3*X*Q*Q))==0
    seed=add(scale(gm(A0,A0),-3*RR),scale(A0,18*X*Q-8),pos(9*Q-3*h*RR*X))
    seedchart=add(scale(a,4),scale(gm(a,a),-3*RR),pos(2/X))
    assert same(seed,seedchart) and gc(seed,pos(RR))=={}
    assert all(s.denom(s.cancel(c))==1 or s.denom(s.cancel(c)).is_Integer for c in seed.values())
    C=add(scale(A0,RR),pos(-s.Rational(2,3)))
    assert gc(C,pos(RR))=={}

    @lru_cache(None)
    def delta(f,k):
        if k==0:return f
        return s.factor(s.diff(delta(f,k-1),v)/v)
    def cm(a,b):
        out={}
        for i,ca in a.items():
            for j,cb in b.items():
                for k in range(i+1):
                    n=i+j-k
                    out[n]=out.get(n,0)+ca*comb(i,k)*h**k*delta(cb,k)
        return clean(out)
    def cp(a,n):
        ans={0:s.Integer(1)}
        for _ in range(n):ans=cm(ans,a)
        return ans
    def cc(a,b):return add(cm(a,b),scale(cm(b,a),-1))
    z={1:s.Integer(1)}
    centC={1:r,0:-r*v}
    original_a=add(z,{0:2/(3*r)})
    assert same(add(scale(original_a,r),{0:-r*v-s.Rational(2,3)}),centC)
    assert cc(z,{0:v})=={0:h/v}
    def support(op):
        out=[]
        for j,c in op.items():
            for term in s.Add.make_args(s.expand(c)):
                powers=term.as_powers_dict()
                n=int(powers.get(h,0));m=int(powers.get(v,0))
                coeff=s.factor(term/(h**n*v**m))
                assert not coeff.has(v,h)
                out.append((n,m,j,coeff,5*n+j+2*m-2))
        return out
    quartic=scale(cp(centC,4),eps*h*h*v**-6)
    cubic=scale(cp(centC,3),h*v**-3)
    samples=[]
    for n,m,j in [(2,-6,4),(1,-3,3),(2,-5,4),(2,-7,4),(3,-8,5)]:
        op=scale(cp(centC,j),h**n*v**m)
        data=support(op); low=min(x[4] for x in data)
        expected=5*n+j+2*m-2
        assert low==expected
        bottom=[x for x in data if x[4]==low]
        assert len(bottom)==1 and bottom[0][:3]==(n,m,j)
        samples.append({'n':n,'m':m,'j':j,'least_weight':low,'unique_least_monomial':True})
    qsupport=support(quartic)
    assert min(x[4] for x in qsupport)==0
    assert [(n,m,j) for n,m,j,c,w in qsupport if w==0]==[(2,-6,4)]
    assert max(quartic)==4
    assert s.factor(quartic[4]-eps*h*h*r**4*v**-6)==0
    assert all(n>=2 for n,m,j,c,w in qsupport if j>=3)
    # Euler polynomial representations in the canonical v valuation.
    Cj=lambda j:s.prod(H-2*i*h for i in range(j))
    for j in range(1,7):
        hj={0:s.Integer(1)}
        euler={1:v*v}
        # Verify v^-2j C_j(H) = a^j by direct finite differential multiplication.
        pol={0:s.Integer(1)}
        for i in range(j):pol=cm(pol,add(euler,{0:-2*i*h}))
        assert same(scale(pol,v**(-2*j)),cp(z,j))
    C4=s.expand(Cj(4))
    phi14=s.cyclotomic_poly(14,t)
    good=1+t*t+t**4+t**6
    bad=(1+t*t)*(1+t**7)
    assert s.rem(good,phi14,t)!=0
    assert s.rem(bad,phi14,t)==0
    # Corruption control: changing the second pair of roots by seven kills the criterion.
    assert good!=bad and s.rem(good,phi14,t)!=s.rem(bad,phi14,t)
    monomial_controls=[]
    for j in range(2,8):
        for m in range(-5,0):
            N=2*j-m
            char=sum(t**(2*i) for i in range(j))
            assert s.rem(char,s.cyclotomic_poly(N,t),t)!=0
            monomial_controls.append((m,j,N))
    # Written degree comparison checked at representative polynomial degrees and shifts.
    degree_controls=[]
    for j in range(2,6):
        cpoly=s.prod(H-(i+1)*h for i in range(j))
        for e in range(0,5):
            dpoly=H**e+3 if e else s.Integer(1)
            for N in (3,7,14):
                boundary=s.expand(cpoly.subs(H,H+N*h)*dpoly-dpoly.subs(H,H-N*h)*cpoly)
                coeff=s.expand(boundary).coeff(H,j+e-1)
                assert s.factor(coeff-N*h*(j+e))==0
                degree_controls.append((j,e,N))
    # Leading coefficient of a k-th power solves the required bottom difference identity.
    power_controls=[]
    for k in range(1,5):
        D=s.prod(C4.subs(H,H-14*i*h) for i in range(k))
        assert s.factor(C4.subs(H,H-14*k*h)*D-D.subs(H,H-14*h)*C4)==0
        assert s.degree(D,H)==4*k
        power_controls.append(k)
    # Actual square ordering gives the distinct, vanishing-character polynomial.
    W=scale(cp(centC,2),v**-3)
    sq=scale(cm(W,W),eps*h*h)
    assert max(sq)==4
    normalized_sq=scale(cm(scale(cp(z,2),v**-3),scale(cp(z,2),v**-3)),1)
    expected_bottom=s.expand(Cj(2).subs(H,H-7*h)*Cj(2))
    euler_poly={0:s.Integer(1)}
    for root in (0,2,7,9):euler_poly=cm(euler_poly,add({1:v*v},{0:-root*h}))
    assert same(normalized_sq,scale(euler_poly,v**-14))
    assert s.factor(expected_bottom-H*(H-2*h)*(H-7*h)*(H-9*h))==0
    # The mixed cubic-quartic candidate: exact admission and symbol facts only.
    lam,xi=s.symbols('lambda xi',nonzero=True)
    assert cc(centC,{0:v})=={0:r*h/v}
    vW=cm({0:v},W);Wv=cm(W,{0:v})
    expected_vW=scale(cp(centC,2),v**-2)
    expected_Wv=add(expected_vW,scale(centC,2*r*h*v**-4),{0:-r*r*h*h*v**-6})
    assert same(vW,expected_vW) and same(Wv,expected_Wv)
    E=add(scale(cp(centC,2),2*v**-2),scale(centC,2*r*h*v**-4),{0:-r*r*h*h*v**-6})
    assert same(E,add(vW,Wv))
    seed_normalized={2:-3*r,0:2*v}
    mixed=add(seed_normalized,scale(add(scale(E,r),cm(centC,W),cm(W,centC)),-3*lam*h),
                        scale(cm(W,W),-3*r*lam**2*h*h))
    shifted=add(z,scale(W,lam*h))
    completed=add(scale(cm(shifted,shifted),-3*r),{0:2*v})
    assert same(mixed,completed)
    assert max(mixed)==4
    assert s.factor(mixed[4]+3*lam**2*h*h*r**5*v**-6)==0
    data=support(mixed)
    assert min(item[4] for item in data)==0
    critical={}
    for n,m,j,coefficient,weight in data:
        if weight==0:critical[j]=critical.get(j,0)+coefficient*v**m
    critical=clean(critical)
    K=scale(cp(z,2),r*r*v**-3)
    critical_root=add(z,scale(K,lam))
    critical_expected=add(scale(cm(critical_root,critical_root),-3*r),{0:2*v})
    critical_expected=clean({i:a.subs(h,1) for i,a in critical_expected.items()})
    assert same(critical,critical_expected)
    weighted_terms=[]
    for j,coefficient in critical.items():
        for term in s.Add.make_args(s.expand(coefficient)):
            exponent=int(term.as_powers_dict().get(v,0))
            weighted_terms.append((exponent+3*j,term*xi**j))
    assert max(weight for weight,term in weighted_terms)==6
    top=s.factor(sum(term for weight,term in weighted_terms if weight==6))
    expected_top=-3*r*(xi+lam*r*r*v**-3*xi**2)**2
    assert s.factor(top-expected_top)==0
    cubic_top=s.expand(top).coeff(xi,3)
    assert s.factor(cubic_top+6*lam*r**3*v**-3)==0
    # Positive control in ordinary b coordinates: P=b^j+v has mate -b.
    def ordinary_mul(a,b):
        out={}
        for i,ca in a.items():
            for j,cb in b.items():
                for k in range(i+1):
                    n=i+j-k
                    out[n]=out.get(n,0)+ca*comb(i,k)*h**k*s.diff(cb,v,k)
        return clean(out)
    for j in range(2,6):
        PP={j:s.Integer(1),0:v};T={1:s.Integer(-1)}
        assert add(ordinary_mul(PP,T),scale(ordinary_mul(T,PP),-1))=={0:h}
        assert s.rem(sum(t**i for i in range(j)),s.cyclotomic_poly(j,t),t)==0
    output=args.output or Path(tempfile.mkdtemp(prefix='global-quartic-'))/'global_quartic_face_controls.json'
    output.parent.mkdir(parents=True,exist_ok=True)
    receipt={
        'scope':'Exact global admission, negative-pole exclusions and closed square/mixed quartic families; mixed exclusion uses the separate written adjoint/residue proofs',
        'status':'PASS','environment':{'python':platform.python_version(),'sympy':s.__version__},
        'global_coordinates':{'A0':'Q+3XQ^2-XB/3=a-v','C':'R*A0-2/3=r*(centered_a-v)',
                              'seed':'-3R*A0^2+(18XQ-8)A0+9Q-3hRX=4a-3ra^2+2v'},
        'support_controls':samples,'quartic_unique_critical_term':'epsilon*h^2*r^4*v^-6*centered_a^4',
        'quartic_bottom':{'N':14,'C':str(C4),'character':str(good),'primitive14_character_nonzero':True},
        'general_monomial_controls':monomial_controls,'degree_controls':degree_controls,
        'power_difference_controls':power_controls,
        'closed_square_ordering':{'candidate':'S_b+epsilon*h^2*(X^3*C^2)^2','bottom_roots_in_h_units':[0,2,7,9],
                                   'character':str(bad),'primitive14_character_zero':True,'mate_constructed':False,
                                   'closure':'All-order square-face residue proof and QOI1 exclude every finite mate'},
        'admitted_mixed_family':{'candidate':'S_b-3lambda*h*(R*E+C*W+W*C)-3R*lambda^2*h^2*W^2',
              'W':'X^3*C^2','E':'2X^2*C^2+2R*h*X^4*C-R^2*h^2*X^6',
              'Wv_identity':True,'completed_square_identity':True,'critical_coefficients':{str(i):str(a) for i,a in critical.items()},
              'top_weight':6,'top_symbol':str(top),'top_cubic_coefficient':str(cubic_top),
              'outside_QOI1_strict_quartic_cone':True,'mate_constructed':False,
              'mate_excluded_by_separate_written_proof':'QAS1 twisted adjoint; QAS2 universal residue for the specified nonsymmetric perturbation'},
        'negative_control':'Character criterion deliberately fails for actual pair b^j+v,-b; corrupted quartic divisor rejected',
        'hashes':{name:sha(HERE/name) for name in ['GLOBAL_ADMISSION_AND_QUARTIC_CRITICAL_FACE.md','check_global_quartic_face.py',
                    'EXACT_POSITION_CENTRALIZER_AND_COMPLETION.md','CRITICAL_LAURENT_FACE_OBSTRUCTION.md',
                    'SQUARE_ORDERED_QUARTIC_FACE.md','ORDERED_QUARTIC_INFINITY_OBSTRUCTION.md']},
        'coverage_boundary':'Finite exact Ore, mixed-admission and critical-symbol identities plus representative divisor controls; universal exclusions rely on written proofs. Mixed admission calculations do not establish mate exclusion; that separate written adjoint/residue proof is explicitly linked. No actual Long mate or global completion.',
        'independent_replay':False,'formal_verification':False,'external_review':False}
    output.write_text(json.dumps(receipt,indent=2)+'\n')
    print('Global admission/quartic face controls PASS')
    print('Verification output:',output)
if __name__=='__main__':main()
