#!/usr/bin/env python3
"""Reproducible standard-library finite exact uniform diagnostics.

Run without Python -O from any directory. No output file is written unless
--output is supplied. These controls do not verify the analytic proof.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
from math import gcd,isqrt
from functools import lru_cache
import argparse,hashlib,json,sys
if not __debug__:
    raise RuntimeError("Exact controls require assertions; run without Python -O.")
parser=argparse.ArgumentParser(description="Finite exact asymmetric arithmetic controls; not an analytic proof verifier.")
parser.add_argument("--output",type=Path,help="Optional JSON receipt; otherwise print only.")
args=parser.parse_args()
ROOT=Path(__file__).resolve().parents[3]
@lru_cache(None)
def factors(n):
    out={};p=2
    while p*p<=n:
        while n%p==0:out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=out.get(n,0)+1
    return out

def weight(p,e,shift,side):
    if e==0:return F(1)
    if side=='A':return F(shift+e+1,shift+1)**2 if p%12==1 else F(0)
    return F(shift+e+1,shift+1) if p%3==1 else F(0)
def series(p,shift,side):
    x=F(1,p);q=1-x
    if side=='A':
        return 1/q+2*x/((shift+1)*q*q)+x*(1+x)/((shift+1)**2*q**3) if p%12==1 else F(1)
    return 1/q+x/((shift+1)*q*q) if p%3==1 else F(1)
def I(p,n,s,t):
    x=F(1,p);q=1-x;A=lambda e:weight(p,e,s,'A');B=lambda e:weight(p,e,t,'B')
    ta=series(p,s,'A')-sum(A(i)*x**i for i in range(n+1))
    tb=series(p,t,'B')-sum(B(i)*x**i for i in range(n+1))
    return q*sum(A(i)*B(i)*x**i for i in range(n))+x**n*(1-2*x)*A(n)*B(n)+q*(B(n)*ta+A(n)*tb)
def local_formula(p,a,b,n,s,t):
    common=min(a,b)
    if n<common:return F(0)
    scale=p**common;a-=common;b-=common;n-=common;x=F(1,p);q=1-x
    if not a and not b:ans=I(p,n,s,t)
    elif a:
        if n<a:ans=weight(p,n,t,'B')*q*series(p,s,'A')
        else:
            multiplier=weight(p,a,t,'B')
            ans=multiplier*I(p,n-a,s,t+a) if multiplier else F(0)
    else:
        if n<b:ans=weight(p,n,s,'A')*q*series(p,t,'B')
        else:
            multiplier=weight(p,b,s,'A')
            ans=multiplier*I(p,n-b,s+b,t) if multiplier else F(0)
    return scale*ans

def valuation(n,p):
    e=0
    while n and n%p==0:e+=1;n//=p
    return e
@lru_cache(None)
def residue_weights(p,h,shift,side):
    mod=p**h;x=F(1,p);q=1-x
    tail=series(p,shift,side)-sum(weight(p,i,shift,side)*x**i for i in range(h))
    return [q*x**(-h)*tail]+[weight(p,valuation(n,p),shift,side) for n in range(1,mod)]
def exact_residue_integral(p,a,b,n,s,t):
    h=n+1;mod=p**h
    wa=residue_weights(p,h,s,'A');wb=residue_weights(p,h,t,'B')
    if a<b:a,b=b,a;wa,wb=wb,wa
    common=p**min(b,h);reduced=mod//common
    inv=pow((p**b)//common,-1,reduced) if reduced>1 else 0
    total=F(0)
    for xx in range(mod):
        if not wa[xx]:continue
        rhs=(p**n-p**a*xx)%mod
        if rhs%common:continue
        yy=(rhs//common*inv)%reduced if reduced>1 else 0
        subtotal=sum(wb[yy+j*reduced] for j in range(common))
        total+=wa[xx]*subtotal
    return total/mod

lc=zero=commoncases=0
for p,s,t in [(13,2,3),(7,0,3),(5,0,0)]:
    cases=list(product(range(3),range(3),range(3)))
    if p in [7,13]:cases +=[(2,0,3),(0,2,3),(2,1,3),(1,2,3)]
    for a,b,n in cases:
        got=exact_residue_integral(p,a,b,n,s,t);want=local_formula(p,a,b,n,s,t)
        assert got==want,(p,a,b,n,got,want)
        lc+=1;zero+=int(want==0);commoncases+=int(min(a,b)>0)

@lru_cache(None)
def shifted(n,shift,side):
    fs=factors(shift);answer=F(1)
    for p,e in factors(n).items():answer*=weight(p,e,fs.get(p,0),side)
    return answer
@lru_cache(None)
def count(B,C,u,v,N):
    return sum((shifted(m,B,'A')*shifted((N-u*m)//v,C,'B') for m in range(1,(N-1)//u+1) if (N-u*m)%v==0),F(0))
reductions=obstructions=0
for N,B,C,u,v in product(list(range(2,501,6))+[1352], [1,169],[1,49,13],[1,13,169],[1,7,13,91]):
    original=count(B,C,u,v,N);r=gcd(u,v)
    if N%r:
        reduced=F(0);obstructions+=1
    else:
        uu=u//r;vv=v//r;NN=N//r;d=gcd(uu,NN);e=gcd(vv,NN)
        if any(p%12==7 for p in factors(e)):
            reduced=F(0);obstructions+=1
        else:
            reduced=shifted(e,B,'A')*shifted(d,C,'B')*count(B*e,C*d,uu//d,vv//e,NN//(d*e))
            assert gcd((uu//d)*(vv//e),NN//(d*e))==1
    assert original==reduced,(B,C,u,v,N,original,reduced)
    reductions+=1

def chi3(n):return 0 if n%3==0 else (1 if n%3==1 else -1)
def base(n):
    if gcd(n,6)>1:return 0
    return sum(chi3(d) for d in range(1,n+1) if n%d==0)
def kernel(n,C):
    value=F(1);fc=factors(C)
    for p,e in factors(n).items():
        if p in fc:
            value*= -F(fc[p],fc[p]+1) if e==1 else 0
        elif p>=5 and p%3==2:
            value*= -1 if e==2 else 0
        else:value=F(0)
    return value
kernels=0
for C in [7,49,91]:
    for n in range(1,501):
        reconstructed=sum((kernel(k,C)*base(n//k) for k in range(1,n+1) if n%k==0),F(0))
        assert reconstructed==shifted(n,C,'B')
        if kernel(n,C):assert n%6==1
        kernels+=1
assert kernel(7,7)==-F(1,2) and 7%12!=1
smallmod=[]
for v,N in product([1,7],[2,8]):
    density=12*sum(F(1,2) for b in [1,7] if (1+v*b-N)%12==0)
    assert density==6;smallmod.append([v,N,str(density)])
paths=['research/arithmetic/ASYMMETRIC_UNIFORM_COEFFICIENTS.md','research/arithmetic/ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md','research/arithmetic/FIXED_PERIODIC_DIVISOR_BOXES.md','research/arithmetic/sources/blomer-grimmelt-li-rydin-myerson-2111.01601v5.pdf','research/arithmetic/sources/heath-brown-d3-progressions-1986.pdf']
result={'scope':'Exact local, coefficient reduction, p7 shifted-kernel and fixed-modulus diagnostics; not a proof of dispersion','exact_local_cases':lc,'exact_zero_local_cases':zero,'common_coefficient_local_cases':commoncases,'global_coefficient_reduction_cases':reductions,'global_obstruction_cases':obstructions,'p7_shifted_kernel_cases':kernels,'negative_control_mod12_support':'e_G7(7)=-1/2, so support 1 mod12 is false','small_modulus_density':smallmod,'source_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}}
result['checker']={'path':Path(__file__).relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version.split()[0]}
result['result']='PASS'
serialized=json.dumps(result,indent=2)+'\n'
if args.output is not None:
    args.output.write_text(serialized)
print(serialized,end='')
