#!/usr/bin/env python3
"""Matched endpoint data, exact-critical-point full inverse and conditional Robin budget; no inverse uniform over a nonzero parameter window, actual C_R bound or numerical SFC cutoff."""
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
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,default=HERE/'quantitative_gap_two_endpoint_data_verification.json');args=ap.parse_args();repo=args.repo_root.resolve()
    ledger_path=HERE/'quantitative_gap_two_endpoint_data_source_ledger.json';ledger=json.loads(ledger_path.read_text())
    mandatory={('local','QUANTITATIVE_GAP_TWO_ENDPOINT_DATA.md'),('local',Path(__file__).name),('repo','research/factorial/image-rigidity/GAP_TWO_LARGE_LOWER_EXPONENT_ENDPOINT_THEOREM.md')}
    mandatory|={('local','HERMITE_CRITICAL_CUBIC_INVERSE_DISK.md'),('local','sources/olde-daalhuis-temme-1994-uniform-airy.pdf'),('repo','research/factorial/image-rigidity/gap_two_large_lower_endpoint_primary_audit.md')}
    assert ledger.get('schema')==1
    assert len(mandatory)==6 and len(ledger['files'])==6
    assert mandatory=={(v['root'],v['path'])for v in ledger['files']};pins=[]
    for row in ledger['files']:
        p=(HERE if row['root']=='local'else repo)/row['path'];assert p.is_file()and sha(p)==row['sha256'],row['path'];pins.append(row)
    r=F(1,3);w=r/2
    def cosh_upper(x):return(eu(x)+eu(-x))/2
    def sinh_upper(x):return(eu(x)-el(-x))/2
    def sinh_lower(x):return(el(x)-eu(-x))/2
    ck('Hermite u Rouche',F(101,2000)+cosh_upper(r)-1-r*r/2,r*r/2)
    S=(F(3,4)/r**3)*(sinh_upper(2*r)-2*r)
    ck('Hermite critical-action cube root input',S-1,F(1,40))
    ck('principal cube root transfer',F(1,3)*F(40,39)*F(1,40),F(1,117),'=')
    T=(cosh_upper(w)-sinh_lower(w)/w)/(2-cosh_upper(w))
    ck('Hermite tanh quotient',T,F(1,100))
    ck('Hermite sinh quotient',sinh_upper(r)/r-1,F(1,50))
    ck('Hermite A square',(F(1,117)+F(1,100))/(1-F(1,100)),F(1,50))
    ck('Hermite A modulus',F(989,1000)**2,F(49,50))
    ck('Hermite zeta modulus',r*r*F(118,117)**2,F(1,8))
    ck('Hermite zeta derivative',(F(1,50)+F(1,117))/(1-F(1,117)),F(3,100))
    de=F(1,1000);v=F(3,2)
    ck('outer selected-pair Rouche',F(11,100)*v+de*(v**3/4+v)+de*de*v*v/4,F(1,4))
    RW=F(11,100)+F(1,2)+de*v*v/4;ck('outer denominator radius',RW,F(9769,16000),'=')
    Lerr=de*de/4+de*(v+2*v)/4;ck('outer quotient numerator',Lerr,F(4501,4000000),'=')
    Qerr=(Lerr+de*v/2)/(1-de*v/2);ck('outer multiplicative quotient error',Qerr,F(1,500))
    ck('outer K right half-plane',F(3,5),1/(1+RW)-F(1,500)/(1-RW))
    ck('outer K modulus',F(501,500)/(1-RW),F(13,5))
    logK=(de/4)/(1-Lerr)+(de/2)/(1-de*v/2)+(1+de*v/2)/(1-RW)
    ck('outer logarithmic derivative',logK,3)
    ck('J cube-root lower',F(4,5)**3,F(3,5));ck('J cube-root upper',F(13,5),F(7,5)**3)
    ck('stationary auxiliary H modulus squared',F(7,3),4)
    ck('stationary auxiliary B modulus',3*F(5,4),4)
    ck('outer matched zeta',F(1,4)*F(7,5)**2,F(49,100),'=')
    small=F(1,10**8);ratio=small/(de-small);ck('Cauchy transfer ratio',ratio,F(1,99999),'=')
    ck('Cauchy A change',2*ratio,F(21,10**6));ck('Cauchy B change',4*ratio,F(41,10**6))
    ck('Cauchy zeta change',F(49,100)*ratio,F(1,200000))
    ck('actual A modulus',F(49,50),F(989,1000)-F(21,10**6))
    ck('actual B modulus',F(500,989)+F(41,10**6),F(13,25))
    ck('actual zeta modulus',r*r*F(118,117)**2+F(1,200000),F(1,8))
    ck('outer zeta z-derivative',F(49,100)/F(9,1000),F(490,9),'=')
    ck('Cauchy zeta derivative change',F(490,9)*ratio,F(3,5000))
    ck('actual zeta derivative',F(3,100)+F(3,5000),F(1,4))
    ck('critical full inverse boundary',F(51,100)**3*F(8,5),F(3,5)**3)
    ck('critical true Jacobian Cauchy bound',F(3,5)/(F(1,2)-F(2,5)),6,'=')
    ck('critical true divided-amplitude derivative',F(2,5)*6/((F(2,5)-F(3,10))**2*F(2,5)**2),1500,'=')
    # Conditional huge-order Robin budget, explicitly not a remainder certificate.
    ck('huge order two-thirds bound',(4*10**10)**3,(10**16)**2)
    K=(F(7,10)+F(1,4))/(F(9,10)-F(1,10**6))
    ck('conditional effective Robin K',K,F(11,10))
    ck('conditional Robin product',F(11,10)**2/F(8),F(1,4))
    ck('conditional Robin cubic',F(11,10)**3/F(10**16),F(1,16))
    ck('conditional actual radial separation',12*F(11,10),2*10**8)
    # Universal rational identities underlying the Hermite amplitude formulas.
    one={(0,0):F(1)};E={(1,0):F(1)};b={(0,1):F(1)}
    E1=pl(E,one);Em=pl(E,ps(one,-1));E2m=pl(pm(E,E),ps(one,-1))
    assert pm(E1,E2m)==pm(pm(E1,E1),Em)
    assert pm(Em,E2m)==pm(pm(Em,Em),E1)
    assert ps(pm(b,pm(E1,Em)),-2)==ps(pm(b,E2m),-2) # AB numerator reduces to -2b(E²-1).
    # Actual stationary cubic quotient: verify its leading terms without assuming root positions.
    def var(j):return{tuple(int(i==j)for i in range(4)):F(1)}
    V,S,P,D=[var(j)for j in range(4)];o={(0,0,0,0):F(1)}
    quad=pl(pl(pm(V,V),ps(pm(S,V),-1)),P)
    L=pl(pl(o,ps(pm(D,D),-F(1,4))),ps(pm(D,pl(V,S)),F(1,4)))
    product=pm(quad,L)
    leading={k:x for k,x in product.items()if k[0]>=2}
    assert leading=={(3,0,0,1):F(1,4),(2,0,0,0):1,(2,0,0,2):-F(1,4)}
    # Exact weighted-average moments: odd powers cancel; the total mass is1.
    weight_cases=0
    for k in range(18):
        raw=F(0)if k%2 else F(3,4)*(F(2,k+1)-F(2,k+3))
        want=F(0)if k%2 else F(3,(k+1)*(k+3))
        assert raw==want;weight_cases+=1
    assert F(3,4)*(2-F(2,3))==1
    # Gaussian-rational Taylor control: true inverse differs from auxiliary H at orderu².
    def ga(a,b):return(a[0]+b[0],a[1]+b[1])
    def gm(a,b):return(a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
    def gs(a,c):return(a[0]*c,a[1]*c)
    mi=(F(0),F(-1));power=(F(1),F(0));phase=[(F(0),F(0))for _ in range(7)]
    phase[1]=(F(0),F(1));phase[2]=(F(-1,2),F(0))
    for k in range(1,7):power=gm(power,mi);phase[k]=ga(phase[k],gs(power,F((-1)**(k+1),k)))
    tc={1:mi,2:(F(1,4),0),3:(0,F(11,80))};cube={}
    for i,a in tc.items():
        for j,bv in tc.items():
            for k,c in tc.items():cube[i+j+k]=ga(cube.get(i+j+k,(0,0)),gs(gm(gm(a,bv),c),F(1,3)))
    assert all(cube[k]==phase[k]for k in[3,4,5])
    # Reciprocal t_v coefficient: (−i+u/2+33iu²/80)*(i+u/2+13iu²/80)=1+O(u³).
    tv={0:(0,F(-1)),1:(F(1,2),0),2:(0,F(33,80))};jac={0:(0,F(1)),1:(F(1,2),0),2:(0,F(13,80))};prod={}
    for i,a in tv.items():
        for j,bv in jac.items():prod[i+j]=ga(prod.get(i+j,(0,0)),gm(a,bv))
    assert prod[0]==(1,0)and prod[1]==prod[2]==(0,0)
    assert F(13,80)-F(1,8)==F(3,80)
    # Negative controls: auxiliary is not the true Jacobian; removing z-margin destroys derivative control;
    # and the conditional C_R budget must not be mistaken for an established remainder bound.
    assert F(13,80)!=F(1,8)
    assert F(49,100)/F(1,10**12)*ratio>F(1,4)
    badK=(F(7,10)+F(10**12,4*10**10))/(F(9,10)-F(10**12,10**16));assert badK**2/F(8)>F(1,4)
    result={'schema':1,'scope':'Matched endpoint data, exact-critical-point full inverse and conditional Robin budget; no inverse uniform over a nonzero parameter window, actual C_R bound or numerical SFC cutoff','comparisons':ROWS,'universal_algebra_identity_count':5,'weighted_average_moment_cases':weight_cases,'exact_taylor_identity_count':3,'negative_controls':3,'pins':pins,'source_ledger_sha256':sha(ledger_path)}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','rational_guards':len(ROWS),'universal_identities':5,'weighted_moments':weight_cases,'taylor_identities':3,'negative_controls':3,'output':str(args.output)}))
if __name__=='__main__':main()
