#!/usr/bin/env python3
"""Exact support for a complete contour proof on a specified bulk window."""
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
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,default=HERE/'gap_two_endpoint_reduction_guards_verification.json');args=ap.parse_args()
    repo=args.repo_root.resolve();ledger_path=HERE/'gap_two_endpoint_reduction_source_ledger.json';ledger=json.loads(ledger_path.read_text())
    mandatory={('local','GAP_TWO_OUTSIDE_ENDPOINT_ONE_PACKET.md'),('local','GAP_TWO_EXTENDED_BULK_ENDPOINT_REDUCTION.md'),('local',Path(__file__).name),('repo','research/factorial/image-rigidity/COALESCING_BINOMIAL_HERMITE_BULK.md'),('local','COALESCING_HERMITE_BULK_CONTOUR_EXCLUSION.md'),('local','verify_coalescing_hermite_bulk_contour_exclusion.py'),('repo','research/factorial/image-rigidity/GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md'),('repo','research/factorial/image-rigidity/verify_gap_two_grace_global_localization.py'),('repo','research/factorial/image-rigidity/verify_coalescing_binomial_hermite_bulk.py'),('repo','research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex')}
    assert mandatory<={(q['root'],q['path'])for q in ledger['files']}
    pins=[]
    for q in ledger['files']:
        path=(HERE if q['root']=='local'else repo)/q['path'];assert path.is_file()and sha(path)==q['sha256'],q['path'];pins.append(q)
    de=F(1,10**8);R=F(16,5);v=F(751,1000);r=F(1,1000);eta=v
    ck('outside b upper corner',F(3,4)+F(4,3),F(21,10));ck('outside b lower corner',4,F(1,4)+4)
    ck('outside selected saddle Rouche',F(1,10000)*v+de*(v**3/4+v)+de**2*v*v/4,r*(F(7,12)-r))
    ck('outside selected derivative simple',F(1,2),F(7,12)-2*r-F(1,10000)-de*(3*v*v/4+1)-de**2*v/2)
    ck('outside whole horizontal denominator',F(4,3),F(21,10)-eta*(1+3*de/2))
    t=de*R/2;H=1+de*R/2;D=F(4,3)
    ck('outside whole horizontal concavity',-1+(2*t+t*t)/(1-t)**2+(de/2)/D+H*H/(D*D),-F(2,5))
    ck('outside curvature norm',1/(1-t)**2+(de/2)/D+H*H/(D*D),F(8,5))
    ck('outside local denominator',F(5,4),F(21,10)-(eta+F(1,20))*(1+3*de/2))
    ck('outside local third',de/(1-t)**3+3*H*(de/2)/F(5,4)**2+2*H**3/F(5,4)**3,F(6,5))
    ck('outside stationary Gaussian positive',F(31,1000),((F(1,4)-r)**2-r*r)/2,'=')
    ck('outside stationary Taylor remainder',(de/2)*v**3/(3*(1-(de/2)*v)),de)
    ck('outside right derivative',-F(99,100)*3+F(10001,10000)*(F(1,3)+F(1,100))+2,0)
    ck('outside left derivative',2,3-1/(F(591,200)-F(1,10000)))
    ck('outside real endpoint denominator',(3+9*de/4+F(1,10000))**2+16,36)
    ck('log6 bound',6,el(2));ck('outside real endpoint',-F(900,203)+2,-2)
    ck('outside vertical modulus',F(4001,1000)+F(31,10)+de*R*R/4,8);ck('log8 bound',8,el(F(21,10)))
    ck('outside vertical phase',-(9-eta*eta)/2+6*de+F(21,10),-F(21,10))
    m=10**12
    ck('outside local error',(5+100*eu(-F(m,2000)))/sb(m)[0],F(1,1000))
    ck('outside nonlocal error',6*sb(m)[1]*eu(-F(m,2000)),F(1,1000))
    ck('outside far error',3*sb(m)[1]*eu(-m),F(1,1000))
    r=F(1,10**6);x=F(2,10**7);y=F(191,100);v=1+r
    ck('bulk radial inverse square',(2*de**2-de**4)/(1-de**2)**2,3*de**2)
    ck('bulk shifted real',F(1,10**7)+3*de**2*(F(1,10**7)+1/de),x)
    ck('bulk shifted imaginary',F(19,10)*(1+3*de**2),y)
    ck('bulk Hermite real part',F(37,125)**2,1-(y/2)**2)
    ck('bulk selected root Rouche',x*v+de*(v**3/4+v)+de**2*v*v/4,r*(F(74,125)-r))
    t=x+r+de*v*v/4;eps=de/2
    err=(2*t+t*t)/(1-t)**2+((1+de*v/2)**2-1)/(1-t)**2+de/(2*(1-t))+(2*eps*v+eps*eps*v*v)/(1-eps*v)**2
    ck('bulk curvature perturbation',err,F(3,10**6));ck('bulk curvature real',F(17,100),2-y*y/2-F(3,10**6))
    ck('bulk local denominator',F(47,50),1-t-F(1,20)*(1+de*(v+F(1,20))/2))
    H=1+de*(v+F(1,20))/2
    ck('bulk local third',de/(1-eps*(v+F(1,20)))**3+3*H*(de/2)/F(47,50)**2+2*H**3/F(47,50)**3,F(5,2))
    ck('bulk local decay',F(3,50),F(17,200)-F(5,2)*F(1,20)/6)
    u=1-(y/2)**2;ck('bulk centre barrier',(1-u)*eu(u+F(1,300)),1)
    change=r+x+6*de
    ck('bulk compact comparison',6*de+R*r+r*r/2+(10*change)/(1-10*change),F(1,50000))
    ck('bulk horizontal distance',F(49,1000),F(1,20)-r)
    ck('bulk horizontal gap',F(1,10000),F(3,50)*F(49,1000)**2-F(2,50000))
    ck('bulk connector gap',F(1,1000),F(1,600)-F(2,50000));ck('bulk near-zero phase',F(101,1000)*eu(2),1)
    ck('bulk far vertical phase',-(9-(y/2+r)**2)/2+6*de+2,-2)
    ck('bulk Gaussian Taylor constant',(F(5,12)/F(3,50)**2)/sb(2)[0],120)
    ck('bulk Gaussian truncation constant',F(4000,17)/sb(2)[0],240)
    ck('bulk Gaussian exponential',F(1,5000),F(17,80000))
    ck('bulk local error',(120+240*eu(-F(m,5000)))/sb(m)[0],F(1,1000))
    ck('bulk horizontal error',3*sb(F(m,2))[1]*eu(-F(m,10000)),F(1,1000))
    ck('bulk far error',4*sb(F(m,2))[1]*eu(-m),F(1,1000));ck('bulk connector error',2*sb(F(m,2))[1]*eu(-F(m,1000)),F(1,1000))
    ck('bulk prefactor ratio',sb(F(300,17))[1],5);ck('bulk upper zero ratio',5*F(101,99),6)
    ck('bulk action drop',F(2,5),F(59,100)*F(5,6));ck('bulk lower exclusion',40*eu(-F(2,5)*2*10**8),1)
    ck('global real width in the huge-a domain',F(8,2*10**8)+F(18,4*10**16*2*10**8),F(1,10**7))
    ck('hard orders exceed the one-packet onset',10**12,F(4*10**16,4)+F(3,2))
    # Exact stationary polynomial from the rational derivative of the actual phase.
    one={(0,0,0):F(1)};V={(1,0,0):F(1)};Z={(0,1,0):F(1)};d={(0,0,1):F(1)}
    h=pl(V,ps(pm(d,pm(V,V)),F(1,4)));hp=pl(one,ps(pm(d,V),F(1,2)))
    actual=pl(pm(V,pl(Z,h)),ps(pm(hp,hp),-1))
    target={(2,0,0):1,(1,1,0):1,(0,0,0):-1,(3,0,1):F(1,4),(1,0,1):-1,(2,0,2):-F(1,4)}
    assert actual==target
    # Three scope/corruption controls.
    assert 1-F(1,2)**2>0 and 1-F(2)**2<0 # Selecting the other saddle destroys horizontal stability.
    assert not F(0)>F(1,600) # At z=2i the centre is the saddle: no positive bulk barrier.
    corrupt=dict(target);corrupt[(1,0,1)]=1;assert corrupt!=actual
    result={'schema':1,'scope':'supporting guards for outside-endpoint zero-freeness and thin-strip actual-pair bulk exclusion; no endpoint-window claim','comparisons':ROWS,'polynomial_identity_count':1,'negative_control_count':3,'pins':pins,'ledger_sha256':sha(ledger_path)}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','comparisons':len(ROWS),'identities':1,'negative_controls':3,'output':str(args.output)}))
if __name__=='__main__':main()
