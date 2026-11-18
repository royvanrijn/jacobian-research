#!/usr/bin/env python3
"""Portable exact finite controls for the asymmetric U*V arithmetic identities.

Run from any directory without Python -O. The global theorem is analytic:
this checks selected local rational identities, bounded residue integrals,
shift polynomials, and small gcd/kernel/hyperbola identities only. It does
not verify the pointwise asymptotic, its uniform error, or cofactor families.
No file is written unless --output is explicitly supplied.
"""
from fractions import Fraction as F
from math import comb, gcd, isqrt
from pathlib import Path
import argparse
import hashlib
import json
import sys

if not __debug__:
    raise RuntimeError("These finite controls require assertions; run without Python -O.")

parser=argparse.ArgumentParser(
    description="Exact finite arithmetic controls for asymmetric C12+A3 sums; not a global asymptotic verifier."
)
parser.add_argument("--output",type=Path,help="Optional JSON receipt path; otherwise print only.")
args=parser.parse_args()
ROOT=Path(__file__).resolve().parents[3]
SOURCE_PATHS=(
    "research/arithmetic/ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md",
    "research/arithmetic/FIXED_PERIODIC_DIVISOR_BOXES.md",
)
source_hashes={}
source_files_missing=[]
for relative in SOURCE_PATHS:
    source=ROOT/relative
    if source.is_file():
        source_hashes[relative]=hashlib.sha256(source.read_bytes()).hexdigest()
    else:
        source_files_missing.append(relative)
if SOURCE_PATHS[0] not in source_hashes:
    raise FileNotFoundError(ROOT/SOURCE_PATHS[0])
checker=Path(__file__).resolve()
checker_hash=hashlib.sha256(checker.read_bytes()).hexdigest()

def aa(p,e):
    return (e+1)**2 if p%12==1 else int(e==0)
def bb(p,e):
    return e+1 if p%3==1 else int(e==0)
def ha(p,x,v=0):
    if p%12!=1:
        assert v==0
        return F(1)
    q=1-x
    return 1/q + F(2,v+1)*x/q**2 + F(1,(v+1)**2)*x*(1+x)/q**3
def hb(p,x,v=0):
    if p%3!=1:
        assert v==0
        return F(1)
    q=1-x
    return 1/q+F(1,v+1)*x/q**2
def local_integral(p,r):
    x=F(1,p); q=1-x
    ta=ha(p,x)-sum(aa(p,i)*x**i for i in range(r+1))
    tb=hb(p,x)-sum(bb(p,i)*x**i for i in range(r+1))
    return (q*sum(aa(p,i)*bb(p,i)*x**i for i in range(r))
      +x**r*(1-2*x)*aa(p,r)*bb(p,r)
      +q*(bb(p,r)*ta+aa(p,r)*tb))
def strata_integral(p,r):
    x=F(1,p); q=1-x; ans=F(0)
    for i in range(r+1):
        a=aa(p,i); b=bb(p,i)
        if not a*b:
            continue
        part=q if i<r else q*ha(p,x,i)+q*hb(p,x,i)-1
        ans+=x**i*a*b*part
    return ans
def vp(n,p):
    e=0
    while n and n%p==0:
        e+=1; n//=p
    return e
def residue_integral(p,r,L=3):
    assert r<L
    x=F(1,p); q=1-x; modulus=p**L; target=p**r
    ta=ha(p,x)-sum(aa(p,i)*x**i for i in range(L))
    tb=hb(p,x)-sum(bb(p,i)*x**i for i in range(L))
    total=F(0)
    for a in range(modulus):
        b=(target-a)%modulus
        if a==0:
            total+=q*ta*bb(p,r)
        elif b==0:
            total+=q*tb*aa(p,r)
        else:
            total+=F(aa(p,vp(a,p))*bb(p,vp(b,p)),modulus)
    return total

local_cases=0
for p in [5,7,11,13,17,19,23,37,61,73,97]:
    x=F(1,p); q=1-x; ma=q*ha(p,x); mb=q*hb(p,x)
    assert local_integral(p,0)==ma+mb-1
    chi=1 if p%3==1 else -1
    eb=q*(1-chi*x)*hb(p,x)
    assert 1+((eb/(1-chi*x))-1)/ma==(ma+mb-1)/ma
    for r in range(13):
        val=local_integral(p,r)
        assert val==strata_integral(p,r), (p,r)
        assert val>0
        if p%12!=1 and r:
            assert val==q
        local_cases+=1
residue_cases=0
for p in [5,7,11,13]:
    for r in range(3):
        assert residue_integral(p,r)==local_integral(p,r), (p,r)
        residue_cases+=1

shift_cases=0
for v in [0,1,2,10,1000]:
    h=[F((v+i+1)**2,(v+1)**2) for i in range(13)]
    c=[sum((-1)**j*comb(4,j)*h[i-j] for j in range(min(4,i)+1)) for i in range(13)]
    assert all(z==0 for z in c[4:])
    assert all(abs(z)<=16 for z in c)
    h2=[F(v+i+1,v+1) for i in range(13)]
    c2=[sum((-1)**j*comb(2,j)*h2[i-j] for j in range(min(2,i)+1)) for i in range(13)]
    assert all(z==0 for z in c2[2:])
    assert all(abs(z)<=1 for z in c2)
    shift_cases+=1

def factors(n):
    fs={}; p=2
    while p*p<=n:
        while n%p==0:
            fs[p]=fs.get(p,0)+1; n//=p
        p+=1
    if n>1:
        fs[n]=fs.get(n,0)+1
    return fs
def U(n):
    fs=factors(n)
    if any(p%12!=1 for p in fs):
        return 0
    out=1
    for e in fs.values():
        out*=(e+1)**2
    return out
def V(n):
    fs=factors(n)
    if any(p%3!=1 for p in fs):
        return 0
    out=1
    for e in fs.values():
        out*=e+1
    return out
def divs(n):
    return [d for d in range(1,n+1) if n%d==0]
def chi(n):
    return 1 if n%6==1 else -1 if n%6==5 else 0
def ge(n):
    return sum(chi(d) for d in divs(n)) if gcd(n,6)==1 else 0
def ge_hyperbola(n):
    assert n%6==1
    s=isqrt(n)
    return 2*sum(chi(q) for q in divs(n) if q*q<n)+(chi(s) if s*s==n else 0)
def kernel(g,M):
    e=[F(0) for _ in range(M+1)]
    for n in range(1,M+1):
        e[n]=F(V(g*n),V(g))-sum(e[d]*ge(n//d) for d in divs(n) if d<n)
        if e[n]:
            assert n%12==1
    return e
finite_records=[]
for N in list(range(2,201,6))+[338]:
    direct=sum(U(a)*V(N-a) for a in range(1,N))
    strata=F(0); expanded=F(0)
    for g in divs(N):
        if not U(g):
            continue
        M=N//g; e=kernel(g,M)
        for m in range(1,M):
            if gcd(m,M)!=1:
                continue
            f=F(U(g*m),U(g)); G=F(V(g*(M-m)),V(g))
            strata+=U(g)*V(g)*f*G
            if not f:
                continue
            n=M-m
            assert n%6==1
            ex=F(0)
            for k in divs(n):
                if not e[k]:
                    continue
                assert gcd(k,M)==1
                t=n//k
                assert ge_hyperbola(t)==ge(t)
                ex+=e[k]*ge_hyperbola(t)
            assert ex==G, (N,g,m)
            expanded+=U(g)*V(g)*f*ex
    assert direct==strata==expanded, N
    finite_records.append({"N":N,"weighted_count":direct})
# Negative controls reject three concrete omissions in the exact reduction.
# Without the inert-square kernel correction, gE6 would be mistaken for V.
inert_wrong=ge(25)
inert_correct=V(25)
inert_kernel=kernel(1,25)
assert inert_wrong==1 and inert_correct==0
assert inert_kernel[25]==-1
assert inert_wrong+inert_kernel[25]*ge(1)==inert_correct
assert inert_wrong!=inert_correct

# Dropping the factor 2 below the square boundary loses one at 49.
n=49
hyperbola_wrong=sum(chi(q) for q in divs(n) if q*q<n)+chi(isqrt(n))
hyperbola_correct=ge(n)
assert hyperbola_correct==3 and hyperbola_wrong==2
assert hyperbola_wrong!=hyperbola_correct

# At N=26, the nontrivial g=13 stratum supplies the entire count.
target=26
gcd_correct=sum(U(a)*V(target-a) for a in range(1,target))
gcd_wrong=sum(U(a)*V(target-a) for a in range(1,target) if gcd(a,target)==1)
assert gcd_correct==8 and gcd_wrong==0
assert gcd_wrong!=gcd_correct
negative_controls=[
    {"omission":"inert-square correction for V at 25","incorrect":inert_wrong,
     "correct":inert_correct,"rejected":True},
    {"omission":"hyperbola factor 2 at 49","incorrect":hyperbola_wrong,
     "correct":hyperbola_correct,"rejected":True},
    {"omission":"nontrivial gcd strata at N=26","incorrect":gcd_wrong,
     "correct":gcd_correct,"rejected":True},
]
out={
 "scope":"Exact rational local identities, bounded residue integrals, normalized shift polynomials, small gcd/kernel/hyperbola identities, and concrete negative controls only; no global asymptotic, uniform analytic estimate, or cofactor theorem certified.",
 "provenance":{
     "checker":checker.relative_to(ROOT).as_posix(),
     "checker_sha256":checker_hash,
     "source_sha256":source_hashes,
     "optional_source_files_missing":source_files_missing,
     "python_version":sys.version.split()[0],
 },
 "local_integral_and_strata_cases":local_cases,
 "residue_integral_cases":residue_cases,
 "shift_polynomial_cases":shift_cases,
 "finite_identity_cases":finite_records,
 "negative_controls":negative_controls,
 "result":"PASS"
}
serialized=json.dumps(out,indent=2)+"\n"
if args.output is not None:
    args.output.write_text(serialized)
print(serialized,end="")
