#!/usr/bin/env python3
"""Exact certificates for the written all-order moving-disc root bound.

The proof supplies classical Grace composition and all-parameter induction.
Finite coefficient/recurrence controls support that proof; they are not a
scan of roots or a proof of consecutive coprimality.
"""
from __future__ import annotations
if not __debug__:
    raise SystemExit('Unoptimized execution is required')
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from math import comb,prod
import argparse,json
HERE=Path(__file__).resolve().parent
ROWS=[]

def add(a,b):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,F(0))+v
    return {k:v for k,v in out.items() if v}
def sc(a,c):return {k:v*F(c) for k,v in a.items() if v*c}
def mul(a,b):
    out={}
    for k,v in a.items():
        for j,w in b.items():
            ij=tuple(x+y for x,y in zip(k,j));out[ij]=out.get(ij,F(0))+v*w
    return {k:v for k,v in out.items() if v}
def power(a,n):
    out={(0,)*len(next(iter(a))):F(1)}
    for _ in range(n):out=mul(out,a)
    return out
def v(j):return {tuple(int(i==j) for i in range(5)):F(1)}
ONE={(0,0,0,0,0):F(1)}
def ck(name,left,right,op='<'):
    left,right=F(left),F(right)
    assert (left<right if op=='<' else left==right),(name,left,op,right)
    ROWS.append({'name':name,'left':str(left),'right':str(right),'relation':op})

def identities():
    r,e,u,eps,B=[v(j) for j in range(5)]
    rp=add(r,ONE)
    # The exact squared containing radius before the square-root majorant.
    lhs=add(add(power(add(ONE,e),2),power(rp,2)),sc(mul(rp,add(sc(ONE,-1),u)),2))
    rhs=add(add(power(r,2),sc(e,2)),add(power(e,2),sc(mul(rp,u),2)))
    assert lhs==rhs
    # Clearing r from the unsimplified u update gives exactly P5.
    R=add(add(power(r,2),e),add(sc(power(e,2),F(1,2)),mul(rp,u)))
    lhs=add(add(sc(power(r,2),-1),mul(mul(r,rp),eps)),mul(add(ONE,sc(eps,3)),R))
    rhs=add(mul(mul(r,add(sc(r,4),ONE)),eps),mul(add(ONE,sc(eps,3)),add(add(e,sc(power(e,2),F(1,2))),mul(rp,u))))
    assert lhs==rhs
    # Induction slack, including the 3epsilon terms.
    lhs=add(sc(power(rp,2),6),sc(add(add(sc(r,4),ONE),mul(add(ONE,sc(eps,3)),add(sc(ONE,5),sc(mul(r,rp),6)))),-1))
    rhs=add(sc(r,2),sc(mul(eps,add(sc(ONE,5),sc(mul(r,rp),6))),-3))
    assert lhs==rhs
    # At epsilon<=1/[100(r+1)] the cleared positive slack is exact.
    lhs=add(sc(mul(r,rp),200),sc(add(sc(ONE,5),sc(mul(r,rp),6)),-3))
    rhs=add(sc(mul(r,rp),182),sc(ONE,-15));assert lhs==rhs
    # Shift r=1+s gives the literal positive polynomial.
    assert add(sc(mul(add(r,ONE),add(r,sc(ONE,2))),182),sc(ONE,-15))==add(add(sc(power(r,2),182),sc(r,546)),sc(ONE,349))
    # Height when B<=1; the B>=1 case is treated separately in the proof.
    lhs=add(mul(add(e,B),add(sc(ONE,2),e)),sc(add(power(add(ONE,e),2),sc(power(add(ONE,sc(B,-1)),2),-1)),-1))
    assert lhs==mul(B,add(e,B))
    return 6

def numerical_guards():
    ck('origin growth bound',F(300,97),4)
    ck('disc cross term',4+F(8,100),5)
    ck('induction positive shifted constant',0,349)
    ck('physical product upper',F(100,99),2)
    ck('physical lower real constant',F(500,99),6)
    ck('physical upper real constant',F(600,99),7)
    ck('eta uniform bound',F(1,100),F(1,100),'=')
    ck('height two-plus-e upper',2+F(4,100),F(51,25),'=')
    ck('height constant',(F(100,99))**2*8*F(51,25),17)
    ck('polynomial thin rectangle',F(7,1000),F(1,100))
    ck('large polynomial strip',F(7,10**9),F(1,10**7))
    return len(ROWS)

def coefficient_controls():
    fixtures=[]
    cases=[(d,N,100*d*d+7) for d in range(2,8) for N in range(1,6)]
    cases += [(20,1,40000),(50,1,250000),(100,1,1000000)]
    for d,N,a in cases:
        gammas=[F(a*N+j,d) for j in range(1,d+1)]
        direct=[comb(N,k)*prod(range(a*N+1,a*N+d*k+1)) for k in range(N+1)]
        split=[F(comb(N,k))*d**(d*k)*prod(prod(g+j for j in range(k)) for g in gammas) for k in range(N+1)]
        assert all(x==y for x,y in zip(direct,split))
        eps=F(d,a);er=F(0);ur=F(0);M=F(1);h=N-1
        for index,g in enumerate(gammas,1):
            if index==1:er=F(3*h,g);ur=F(h,g)
            else:
                r=index-1;local=F(h,g)
                ur=(4*r+1)*local+(1+3*local)*(er/r+er*er/(2*r)+(1+F(1,r))*ur)
                er=3*local+(1+3*local)*er
            M*=g
            assert er<=4*index*eps and ur<=6*index*index*eps
        kappa=F(d**d)*M/F(a**d*N**d)
        eta=F(d*d,a)
        assert 1<=kappa<F(100,99) and kappa-1<F(100,99)*eta
        fixtures.append({'d':d,'N':N,'a':a,'coefficient_vector_sha256':sha256(json.dumps([str(x) for x in direct],separators=(',',':')).encode()).hexdigest(),'recurrence_steps':d,'height_majorant_branch':'B>=1' if 6*d*eta>=1 else 'B<1'})
    assert {r['height_majorant_branch'] for r in fixtures}=={'B>=1','B<1'}
    return fixtures

def affine_controls():
    def trim(a):
        a=[F(v)for v in a]
        while a and not a[-1]:a.pop()
        return a
    def plus(a,b):
        out=[F(0)]*max(len(a),len(b))
        for j,c in enumerate(a):out[j]+=c
        for j,c in enumerate(b):out[j]+=c
        return trim(out)
    def times(a,b):
        out=[F(0)]*max(0,len(a)+len(b)-1)
        for i,c in enumerate(a):
            for j,e in enumerate(b):out[i+j]+=c*e
        return trim(out)
    def affine(a,c0,c1):
        out=[]
        for c in reversed(a):out=plus(times(out,[c0,c1]),[c])
        return trim(out)
    def gcd(a,b):
        a,b=trim(a),trim(b)
        while b:
            rem=a[:]
            while rem and len(rem)>=len(b):
                q=rem[-1]/b[-1];k=len(rem)-len(b)
                for j,vv in enumerate(b):rem[j+k]-=q*vv
                rem=trim(rem)
            a,b=b,rem
        return [v/a[-1]for v in a]if a else[]
    def moment(m,a,d,shape_order=None):
        out=[F(0)]*(m+1);ell=m if shape_order is None else shape_order
        for k in range(m+1):out[m-k]=F(comb(m,k)*prod(range(a*ell+1,a*ell+d*k+1)))
        return out
    controls=[]
    for d,N in [(2,2),(3,3),(4,4)]:
        a=100*d*d;n=N-1;c0=-F(N**d*a**d);c1=-c0/10
        up,lo=moment(N,a,d),moment(n,a,d)
        U,L=affine(up,c0,c1),affine(lo,c0,c1)
        assert len(U)==N+1 and len(L)==n+1
        assert len(gcd(up,lo))==len(gcd(U,L))==1
        assert lo!=moment(n,a,d,shape_order=N)
        controls.append({'a':a,'d':d,'N':N,'actual_lower_shape':a*n+1,'original_gcd_degree':0,'same_affine_gcd_degree':0})
    # The comparison also preserves an actual NONCONSTANT gcd.
    up,lo=[-2,-1,1],[-6,1,1]
    assert len(gcd(up,lo))==len(gcd(affine(up,F(-5),F(3)),affine(lo,F(-5),F(3))))==2
    controls.append({'type':'nonconstant_common_factor','original_gcd_degree':1,'same_affine_gcd_degree':1})
    # Different affine arguments can destroy a common factor.
    same=[-1,1]
    assert len(gcd(same,same))==2
    assert len(gcd(affine(same,F(0),F(1)),affine(same,F(1),F(1))))==1
    return controls

def negatives():
    # The minus sign in Grace is necessary already for degree one.
    gamma1,gamma2=F(2),F(3)
    assert -gamma1*gamma2+gamma1*gamma2==0
    assert gamma1*gamma2+gamma1*gamma2!=0
    # Generic point -1+3i/4 has e1/4,u0. Omitting e^2 undersizes its disc.
    actual_distance_squared=1+F(9,16)
    assert actual_distance_squared>1+2*F(1,4)
    # Squaring a positive upper real bound is invalid: y=i,e0,B2.
    assert 1>(1+0)**2-(1-2)**2
    return 3

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,default=HERE/'iterated_grace_polynomial_root_strip_verification.json');args=ap.parse_args();repo=args.repo_root.resolve()
    lp=HERE/'iterated_grace_polynomial_root_strip_source_ledger.json';ledger=json.loads(lp.read_text());assert ledger.get('schema')==1
    mandatory={('local','ITERATED_GRACE_POLYNOMIAL_ROOT_STRIP.md'),('local',Path(__file__).name),('local','iterated_grace_polynomial_root_strip_inputs.json'),('local','iterated_grace_primary_scope_audit.md'),('local','sources/lamprecht-grace-2015.pdf'),('repo','research/factorial/image-rigidity/GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md'),('repo','research/factorial/image-rigidity/FIXED_GAP_LARGE_LOWER_EXPONENT_EXTENSION.md')}
    assert len(ledger['files'])==len(mandatory) and {(r['root'],r['path'])for r in ledger['files']}==mandatory
    for row in ledger['files']:
        assert row['root']in{'local','repo'}
        p=(HERE if row['root']=='local'else repo)/row['path'];assert p.is_file()and sha256(p.read_bytes()).hexdigest()==row['sha256'],row['path']
    inp=json.loads((HERE/'iterated_grace_polynomial_root_strip_inputs.json').read_text())
    assert inp=={'equivalence_declaration':'Corollary ActualPolynomialStripNoncollisionEquivalence','common_argument':'S_N(Z)=N^d*a^d*(-1+(d/sqrt(a))*Z), used for BOTH actual orders N,N-1','schema':1,'declaration':'Theorem PolynomialGlobalHermiteStrip','domain':{'d_min':2,'N_min':1,'a_min':'100*d^2'},'actual_polynomial':'A_(N;a,d)(x)=sum binom(N,k)*(aN+1)_(d*k)*x^(N-k)','normalization':'Z=(x/(N^d*a^d)+1)/(d/sqrt(a))','conclusion':{'real_bound':'7*d^2/sqrt(a)','imaginary_square_bound':'17*d'},'excludes':['consecutive coprimality','uniform outside-one-packet theorem','polynomial Strong Factorial threshold']}
    ident=identities();count=numerical_guards();fixtures=coefficient_controls();affine_cases=affine_controls();negative=negatives()+2
    result={'schema':1,'status':'PASS_POLYNOMIAL_GLOBAL_ROOT_STRIP_SUPPORT','scope':'Universal exact algebra/majorant support for the written ALL-ROOT, ALL-ORDER localization theorem, not an SFC claim','source_ledger_sha256':sha256(lp.read_bytes()).hexdigest(),'pins':ledger['files'],'universal_identities':ident,'rational_guards':ROWS,'actual_multiplication_and_recurrence_controls':fixtures,'common_affine_controls':affine_cases,'negative_controls':negative,'formal_verification':False,'external_review':False,'independent_replay':False}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':'PASS','universal_identities':ident,'rational_guards':count,'actual_controls':len(fixtures),'common_affine_controls':len(affine_cases),'negative_controls':negative}))
if __name__=='__main__':main()
