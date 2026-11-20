#!/usr/bin/env python3
"""Exact algebra/contour guards for the existential large-lower gap-two theorem."""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import argparse, hashlib, json
if not __debug__:
    raise SystemExit('Refusing -O/-OO: exact checks are required.')
HERE=Path(__file__).resolve().parent
ROWS=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def st(x):
    x=F(x)
    return str(x.numerator) if x.denominator==1 else f'{x.numerator}/{x.denominator}'
def ck(name,left,right,relation='<'):
    left,right=F(left),F(right)
    assert (left<right if relation=='<' else left==right),(name,st(left),relation,st(right))
    ROWS.append({'name':name,'left':st(left),'relation':relation,'right':st(right)})
def eb(x):
    x=F(x)
    if x<0:
        if -x>100:return F(0),1/eb(F(100))[0]
        lo,hi=eb(-x);return 1/hi,1/lo
    term=total=F(1)
    for k in range(1,257):term*=x/k;total+=term
    assert x/258<1
    return total,total+term*x/257/(1-x/258)
def eu(x):return eb(x)[1]
def el(x):return eb(x)[0]
def sb(x):
    x=F(x);q=10**12;k=isqrt(x.numerator*q*q//x.denominator)
    assert F(k*k,q*q)<=x<F((k+1)**2,q*q)
    return F(k,q),F(k+1,q)
def pl(a,b):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,F(0))+v
    return {k:v for k,v in out.items() if v}
def ps(a,c):return{k:v*c for k,v in a.items()if v*c}
def pm(a,b):
    out={}
    for i,u in a.items():
        for j,v in b.items():
            k=tuple(x+y for x,y in zip(i,j));out[k]=out.get(k,F(0))+u*v
    return{k:v for k,v in out.items()if v}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,default=HERE/'gap_two_endpoint_contour_identities_verification.json');args=ap.parse_args()
    repo=args.repo_root.resolve();proof=HERE/'GAP_TWO_LARGE_LOWER_EXPONENT_ENDPOINT_THEOREM.md'
    ledger_path=HERE/'gap_two_large_lower_endpoint_source_ledger.json';ledger=json.loads(ledger_path.read_text())
    mandatory={('local',proof.name),('local',Path(__file__).name),('local','gap_two_large_lower_endpoint_primary_audit.md'),('repo','research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex'),('repo','research/factorial/image-rigidity/GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md'),('repo','research/factorial/image-rigidity/COALESCING_BINOMIAL_HERMITE_BULK.md')}
    mandatory|={('local','sources/dlmf-'+tag+'-endpoint.html')for tag in ['2-4','9-2','9-7','9-9']}
    mandatory|={('repo','research/factorial/image-rigidity/'+name)for name in ['COALESCING_BINOMIAL_SADDLE_BRIDGE.md','COALESCING_HERMITE_BULK_CONTOUR_EXCLUSION.md','GAP_TWO_EXTENDED_BULK_ENDPOINT_REDUCTION.md','GAP_TWO_OUTSIDE_ENDPOINT_ONE_PACKET.md']}
    assert ledger.get('schema')==1
    assert len(mandatory)==14 and len(ledger['files'])==14
    assert mandatory=={(q['root'],q['path'])for q in ledger['files']}
    pins=[]
    for q in ledger['files']:
        path=(HERE if q['root']=='local'else repo)/q['path'];assert path.is_file()and sha(path)==q['sha256'],q['path'];pins.append(q)
    # Universal polynomial identity for the bent vertical cubic contour.
    z={(0,0,0,0):F(1)}
    def mono(e):return{e:F(1)}
    P,Q,X,Y=[mono(tuple(int(i==j)for i in range(4)))for j in range(4)]
    def cm(a,b):return(pl(pm(a[0],b[0]),ps(pm(a[1],b[1]),-1)),pl(pm(a[0],b[1]),pm(a[1],b[0])))
    pair=(X,Y);square=cm(pair,pair);cube0=cm(square,pair);sqpart=cm((P,Q),square);re=pl(sqpart[0],ps(cube0[0],F(1,3)))
    target={(1,0,2,0):1,(1,0,0,2):-1,(0,1,1,1):-2,(0,0,3,0):F(1,3),(0,0,1,2):-1}
    assert re==target
    ck('positive contour cross coefficient',F(2,6),F(1,3),'=')
    ck('positive contour cubic loss',F(8,9)-F(1,3),F(5,9),'=')
    ck('positive contour envelope constant squared',8*4,6**2)
    ck('cos5pi12 exceeds one quarter squared reduction',48,49)
    # Exact rational identities for the negative-axis curve, after clearing y.
    def upoly(a):return{(k,):F(v)for k,v in enumerate(a)if v}
    one=upoly([1]);y=upoly([0,1]);ym=upoly([-1,1]);yp=upoly([2,1]);two=upoly([1,2])
    numer=pm(pm(ym,ym),yp)
    imag=pl(numer,pl(ps(pm(y,pm(y,y)),-1),ps(y,3)))
    assert imag==upoly([2])
    real=pl(numer,pl(ps(pm(y,pm(y,y)),-9),ps(y,9)))
    assert real==ps(pm(ym,pm(two,two)),-2)
    # One-variable Gaussian majorants giving the literal constant50.
    ck('negative middle square-root lower',F(16,25),F(7,9))
    ck('negative middle D factor',F(4,9),F(2*4,1)/(9*F(3,2))*F(4,5))
    ck('negative middle arc derivative',F(4,3)+F(1,2)/(3*F(1,4)*F(4,5))+1,4)
    ck('negative middle Gaussian constant squared',36*4*F(3,2),15**2)
    ck('negative upper D quadratic',F(4,81)*F(3,2),F(2,27),'=')
    ck('negative upper arc derivative',1+F(2,3)/F(3,2)+1,3)
    ck('negative upper Gaussian constant squared',F(243,8)*4*F(3,2),F(27,2)**2,'=')
    ck('negative lower arc derivative',F(1,2)+F(2,3)+F(1,2),2)
    ck('negative lower Gaussian constant squared',72*4*F(3,2),21**2)
    ck('negative full envelope',15+F(27,2)+21,50)
    # The explicit critical phase is quartically lower on the horizontal.
    # d/dt[t-log(1+t)-t^2/(2(1+t))]=t^2/(2(1+t)^2).
    t=upoly([0,1]);den=pm(pl(one,t),pl(one,t))
    derivative_numer=pl(ps(pm(t,pl(one,t)),2),ps(pl(ps(t,2),pm(t,t)),-1))
    assert derivative_numer==pm(t,t)
    # Polynomial divided amplitude: even and odd interpolation through all monomials.
    T={(1,0):F(1)};Z={(0,1):F(1)};one2={(0,0):F(1)};divisor=pl(pm(T,T),ps(Z,-1));division_cases=0
    for degree in range(18):
        h={(degree,0):F(1)};j=degree//2
        A={(0,j):F(1)}if degree%2==0 else{}
        B={(0,j):F(1)}if degree%2 else{}
        G={}
        for k in range(j):G[(2*(j-1-k)+(degree%2),k)]=F(1)
        assert h==pl(pl(A,pm(B,T)),pm(divisor,G));division_cases+=1
    # Actual gap-two stationary polynomial from the exact rational phase derivative.
    V0={(1,0,0):F(1)};Z0={(0,1,0):F(1)};Delta={(0,0,1):F(1)};one3={(0,0,0):F(1)}
    h0=pl(V0,ps(pm(Delta,pm(V0,V0)),F(1,4)));hp0=pl(one3,ps(pm(Delta,V0),F(1,2)))
    stationary=pl(pm(V0,pl(Z0,h0)),ps(pm(hp0,hp0),-1))
    expected={(2,0,0):1,(1,1,0):1,(0,0,0):-1,(3,0,1):F(1,4),(1,0,1):-1,(2,0,2):-F(1,4)}
    assert stationary==expected
    # Exact Robin continuation identity using Ai''=z Ai.
    # Variables t,lambda,z,U=Ai(z),V=Ai'(z); denominator=1-t^2lambda^2z.
    def var(j):return{tuple(int(i==j)for i in range(5)):F(1)}
    t,l,z,U,V=[var(j)for j in range(5)];one5={(0,0,0,0,0):F(1)}
    D=pl(one5,ps(pm(pm(t,t),pm(pm(l,l),z)),-1));H=pl(U,pm(pm(t,l),V))
    lhs=pl(pl(ps(pm(l,V),-1),pm(pm(t,l),ps(pm(pm(l,z),U),-1))),pm(pm(l,V),D))
    rhs=ps(pm(pm(pm(t,pm(l,l)),z),H),-1)
    assert lhs==rhs
    H0=pl(U,pm(l,V));Hz=pl(V,pm(pm(l,z),U));Hz_identity=pl(pm(V,pl(one5,ps(pm(pm(l,l),z),-1))),pm(pm(l,z),H0));assert Hz==Hz_identity
    ck('Robin denominator bootstrap',1-F(1,4)-2*F(1,16),F(5,8),'=')
    ck('Robin speed bootstrap',F(8,5),2)
    ck('endpoint drift coefficient',F(3,4)*15-F(15,4),F(1,2)*15,'=')
    # Exact actual radial-shift lower bound rho>2/N.
    N=upoly([0,1]);Nm=upoly([-1,1]);rho_numer=pl(pm(N,upoly([-1,2])),ps(pm(Nm,Nm),-2))
    assert rho_numer==upoly([-2,3])
    ck('minimum terminal positive rho gap',0,3*2-2)
    # Exact Gaussian-rational endpoint Taylor coefficients.
    def ga(a,b):return(a[0]+b[0],a[1]+b[1])
    def gm(a,b):return(a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
    def gs(a,c):return(a[0]*c,a[1]*c)
    I=(F(0),F(1));mi=(F(0),F(-1));power=(F(1),F(0));phase=[(F(0),F(0))for _ in range(9)]
    phase[0]=(F(1,2),F(0));phase[1]=I;phase[2]=(F(-1,2),F(0))
    for k in range(1,9):power=gm(power,mi);phase[k]=ga(phase[k],gs(power,F((-1)**(k+1),k)))
    assert phase[1]==phase[2]==(0,0) and phase[3]==(0,F(1,3))and phase[4]==(F(-1,4),0)
    # t=-iu+u^2/4 gives t^3/3=i*u^3/3-u^4/4 through degree4.
    tc={1:mi,2:(F(1,4),F(0))};cube={}
    for i,a in tc.items():
        for j,b in tc.items():
            for k,c in tc.items():cube[i+j+k]=ga(cube.get(i+j+k,(F(0),F(0))),gs(gm(gm(a,b),c),F(1,3)))
    assert cube[3]==phase[3]and cube[4]==phase[4]
    # Corruptions and scope guards.
    bad=dict(target);bad[(0,1,1,1)]=2;assert bad!=re
    bad_div=pl(pl(A,pm(B,T)),pm(ps(divisor,-1),G));assert bad_div!=h
    assert not 1-F(1,4)-2*F(1,2)>0 # Dropping the small-lambda cubic guard breaks the bootstrap.
    result={'schema':1,'scope':'universal algebra and rational contour guards; the written argument supplies analytic uniformity and the nonnumerical threshold','proof_sha256':sha(proof),'checker_sha256':sha(Path(__file__)),'source_ledger_sha256':sha(ledger_path),'pins':pins,'comparisons':ROWS,'universal_identity_count':10,'polynomial_division_cases':division_cases,'negative_controls':3}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','comparisons':len(ROWS),'universal_identities':10,'division_cases':division_cases,'negative_controls':3,'output':str(args.output)}))
if __name__=='__main__':main()
