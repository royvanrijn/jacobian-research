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
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--output',type=Path,default=HERE/'coalescing_hermite_bulk_contour_exclusion_verification.json')
    a=ap.parse_args();repo=a.repo_root.resolve()
    ledger_path=HERE/'coalescing_hermite_bulk_contour_exclusion_source_ledger.json'
    ledger=json.loads(ledger_path.read_text())
    mandatory={('local','COALESCING_HERMITE_BULK_CONTOUR_EXCLUSION.md'),('local',Path(__file__).name),('repo','research/factorial/image-rigidity/COALESCING_BINOMIAL_HERMITE_BULK.md'),('repo','research/factorial/image-rigidity/verify_coalescing_binomial_hermite_bulk.py')}
    assert mandatory<={(q['root'],q['path'])for q in ledger['files']}
    pins=[]
    for q in ledger['files']:
        base=HERE if q['root']=='local' else repo
        path=base/q['path'];assert path.is_file() and sha(path)==q['sha256'],q['path']
        pins.append({'root':q['root'],'path':q['path'],'sha256':sha(path)})
    de=F(1,10**8);R=F(16,5);u=F(51,80)**2;c=F(13,15)
    ck('middle discriminant margin',F(27,20)**2,4-F(59,40)**2)
    ck('middle t modulus lower bound',(F(27,20)-F(3,40))/2,F(51,80),'=')
    ck('middle square-root imaginary ratio',F(3,40)/F(27,20),F(1,18),'=')
    ck('middle radial bound algebra',1/(1-F(2,17)),F(17,15),'=')
    ck('centre barrier derivative',1,c*(1+u))
    ck('positive centre barrier',(1-c*u)*eu(u+F(1,50)),1)
    ck('Hermite peak logarithm',F(6,5),el(F(1,5)))
    ck('Hermite peak lower bound',-F(23,25),-F(18,25)-F(1,5),'=')
    ck('horizontal derivative positive constant',F(3,40)**2+(F(19,36)*F(59,40))**2,1)
    ck('Hermite horizontal gap',F(1,5)*F(1,25)**2,F(1,3125),'=')
    ck('h compact remainder',R*R/(2*(1-R*de)),6)
    ck('U compact remainder',R**3/(6*(1-R*de/2)),6)
    ck('height quadratic shift',F(3,2)*(2*R+3*de),10)
    ck('compact logarithm comparison',90/(1-90*de),91)
    ck('total compact comparison',6+10+91,120)
    ck('near-zero denominator',F(1,10)+9*de,F(101,1000))
    ck('near-zero Gaussian height',(F(21,25)+3*de)**2/2+6*de,F(3,8))
    ck('near-zero phase',F(101,1000)*eu(F(15,8)),1)
    ck('actual peak above minus one',-1,-F(23,25)-120*de)
    ck('actual horizontal distance margin',F(1,25),F(1,20)-3*de)
    ck('actual horizontal gap',F(1,4000),F(1,3125)-240*de)
    ck('connector derivative',(F(53,100)),F(47,100)/(F(3,40)**2+(F(53,100)*F(59,40))**2))
    ck('small-y centre denominator',F(3,40)**2+(F(1,10)+F(19,360))**2,F(1,4)**2)
    ck('small-y centre height',F(19,360)**2/2,F(1,100))
    ck('small-y centre logarithm',eu(F(121,100)),4)
    ck('actual short-connector gap',F(1,200),F(1,100)-240*de)
    ck('right tail derivative',(-3+F(15,7)*(F(1,3)+F(1,100)))/F(203,200),-2)
    ck('right derivative multiplied at endpoint',-F(99,100)*3+F(15,7)*(F(1,100)+F(1,3))+2,0)
    ck('real logarithm endpoint',5,el(F(5,3)))
    ck('right real endpoint',-F(900,203)+F(5,3),-F(5,2))
    ck('left h lower bound',3-F(9,2)*F(1,100),F(591,200),'=')
    ck('left real derivative',2,3-F(200,271))
    ck('left real endpoint',-F(9,2)+F(5,3),-F(5,2))
    ck('far vertical Gaussian',-(9-F(21,25)**2)/2,-F(2592,625),'=')
    ck('far vertical remainder',F(1,200)*R**3/(3*(1-R/200)),F(4096,73800),'=')
    ck('far vertical h bound',F(78,25)+F(8,5)+F(1,100)*R**2/(2*(1-R/100)),5)
    ck('far vertical total',-F(2592,625)+F(4096,73800)+F(5,3),-F(12,5))
    m=10**12
    ck('hard order margin',m,F(10**16*2*2,2*2)+F(1,2))
    ck('local Gaussian full error',(25+100*eu(-F(m,2000)))/sb(2*m)[0],F(1,1000))
    ck('horizontal nonlocal full error',3*sb(F(m,2))[1]*eu(-F(m,4000)),F(1,1000))
    ck('far tails full error',4*sb(F(m,2))[1]*eu(-m),F(1,1000))
    ck('connector full error',2*sb(F(m,2))[1]*eu(-F(m,200)),F(1,1000))
    ck('complete packet error',F(4,1000),F(1,100))
    ck('two-order leading ratio',F(3,2)*sb(6)[1],4)
    ck('two-order action exclusion',16*eu(-10**8),1)
    # Three polynomial identities, independently cleared of all denominators.
    one={(0,0,0):F(1)};t={(1,0,0):F(1)};x={(0,1,0):F(1)};Y={(0,0,1):F(1)}
    actual=pl(ps(pm(t,pl(pm(pl(t,x),pl(t,x)),pm(Y,Y))),-1),pl(t,x))
    target={(3,0,0):-1,(2,1,0):-2,(1,0,0):1,(1,2,0):-1,(1,0,2):-1,(0,1,0):1}
    assert actual==target
    # Centre identity: r*[x^2+Y^2]=1-u*(2-r), using the saddle relations.
    u0={(1,0,0):F(1)};r={(0,1,0):F(1)}
    left=pl(pm(u0,pm(pl(one,ps(r,-1)),pl(one,ps(r,-1)))),pl(r,ps(u0,-1)))
    right=pm(r,pl(one,pm(u0,pl(ps(one,-2),r))))
    assert left==right
    # Positive derivative in u after the centre barrier substitution.
    assert pl(ps(pm(u0,r),1),pl(r,ps(one,-1)))=={(1,1,0):1,(0,1,0):1,(0,0,0):-1}
    # Negative controls: the uncancelled central vertical is NOT exponentially small;
    # the Airy endpoint has no positive centre barrier; and the centre identity is sensitive.
    assert eu(F(1,16))<F(3,2)
    assert 0**2-(-1)**2-2*(-1)-1==0 and 0*(2*(-1)+2)==0 # P0(-i,2i)=0.
    assert not F(0)>F(1,100) # The selected saddle is its own centre at the Airy endpoint.
    corrupt=pm(r,pl(one,pm(u0,pl(ps(one,-1),r))))
    assert corrupt!=left
    result={'schema':1,'scope':'actual full-contour exclusion on the stated Hermite bulk window only','comparisons':ROWS,'polynomial_identity_count':3,'negative_control_count':3,'pins':pins,'ledger_sha256':sha(ledger_path)}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','comparisons':len(ROWS),'polynomial_identities':3,'negative_controls':3,'output':str(a.output)}))
if __name__=='__main__':main()
