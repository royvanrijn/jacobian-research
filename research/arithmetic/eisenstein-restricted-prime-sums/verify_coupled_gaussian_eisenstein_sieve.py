#!/usr/bin/env python3
"""Finite exact diagnostics only; no eventual coverage or UBV claim."""
from fractions import Fraction
from pathlib import Path
from math import gcd,isqrt
import argparse,hashlib,json

if not __debug__:
 raise RuntimeError("Run these exact finite diagnostics without Python -O.")

ROOT=Path(__file__).resolve().parents[3]
ARITH=ROOT/"research/arithmetic"
SOURCE_PATHS=['research/arithmetic/COUPLED_GAUSSIAN_EISENSTEIN_SIEVE.md','research/arithmetic/ASYMMETRIC_SUBGROUP_SPLITTING.md','research/arithmetic/sources/henriot-erratum-2014.pdf','research/arithmetic/sources/granville-koukoulopoulos-1710.01389v3.pdf','research/arithmetic/sources/loughran-matthiesen-1904.12845v2.pdf']
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output",type=Path,help="Optional receipt path; otherwise print only.")
args=parser.parse_args()
if args.output is not None:
 protected={(ROOT/name).resolve() for name in SOURCE_PATHS}
 protected.add(Path(__file__).resolve())
 if args.output.resolve() in protected:
  raise ValueError("Output path would overwrite a proof source or this checker.")

LIMIT=2400
spf=list(range(LIMIT+1))
for p in range(2,isqrt(LIMIT)+1):
 if spf[p]==p:
  for m in range(p*p,LIMIT+1,p):
   if spf[m]==m: spf[m]=p

def fac(n):
 out={}
 while n>1:
  p=spf[n];out[p]=out.get(p,0)+1;n//=p
 return out
F=[{}]+[fac(n) for n in range(1,LIMIT+1)]
def allowed(n,q): return n>0 and all(p%q==1 for p in F[n])
C4=[False]+[allowed(n,4) for n in range(1,LIMIT+1)]
C3=[False]+[allowed(n,3) for n in range(1,LIMIT+1)]
C=[False]+[allowed(n,12) for n in range(1,LIMIT+1)]
counters={}
def add(k,n=1): counters[k]=counters.get(k,0)+n

def check(cond,label):
 if not cond: raise AssertionError(label)
 add(label)

Ds=[1,5,7,17,19,35,85,95,119,133,323]
for N in range(2,1201,12):
 pairs=[(a,N-a) for a in range(1,N) if C4[a] and C3[N-a]]
 for a,b in pairs:
  check(a%12==b%12==1,'separate_residues')
  check(sum(e for p,e in F[a].items() if p%12==5)%2==0 and sum(e for p,e in F[b].items() if p%12==7)%2==0,'separate_even_parities')
 for d in Ds:
  d5=d7=1
  for p,e in fac(d).items():
   check(e==1 and p%12 in (5,7),'supported_squarefree_input')
   if p%12==5:d5*=p
   else:d7*=p
  lhs=sum(a*b%d==0 for a,b in pairs)
  rhs=0
  for x in range(1,(N-1)//d5+1):
   rem=N-d5*x
   if rem>0 and rem%d7==0 and C4[x] and C3[rem//d7]:rhs+=1
  check(lhs==rhs,'exact_coefficient_identity')
  if gcd(d,N)>1:check(lhs==0,'target_prime_zero')
 for z in sorted(set([isqrt(N),max(3,round(N**(1/3))+1)])):
  if z**3<=N or z<3:continue
  for a,b in pairs:
   if any(p<z and p%12 in (5,7) for p in set(F[a])|set(F[b])):continue
   for v,cl in [(a,5),(b,7)]:
    bad=[p for p,e in F[v].items() if p%12==cl for _ in range(e)]
    if bad:
     check(len(bad)==2 and min(bad)>=z,'two_large_prime_structure')
     c=v//(bad[0]*bad[1])
     check(C[c] and c*z*z<=N,'good_cofactor_bound')

# Content drop is exact for every good-prime multiplier.
for h in range(1,301):
 if not C[h]:continue
 for n in range(1,min(101,LIMIT//h+1)):
  for z in (4,12,30):
   for arr,cl in [(C3,7),(C4,5)]:
    def sift(v):return arr[v] and not any(p<z and p%12==cl for p in F[v])
    check(sift(h*n)==sift(n),'exact_good_content_drop')

# Haar ratios at the missing primes, counted in Z/p^2Z.
for p in (5,7,17,19):
 for r in (0,1,2,3):
  N=1 if r==0 else p**r
  total=event=0
  for a in range(p*p):
   b=(N-a)%(p*p)
   if p%12==5:
    ok=b%p!=0;hit=a%p==0
   else:
    ok=a%p!=0;hit=b%p==0
   total+=ok;event+=ok and hit
  check(Fraction(event,total)==(Fraction(1,p-1) if r==0 else 0),'local_g0_all_target_depths')

# Primitive-form small-prime integrals. At good ell, all Q2 valuations
# are accepted; target/leading valuations can be arbitrary after gcd drop.
for ell in (13,37):
 for r in range(4):
  for j in range(4):
   h=min(r,j);M=ell**(r-h);L=ell**(j-h)
   for offset in (5,7):
    mod=ell*ell
    count=sum((12*t+offset)%ell!=0 for t in range(mod))
    check(Fraction(count,mod)==Fraction(ell-1,ell),'good_local_all_valuation_depths')
for ell in (5,7,11,17,19,23):
 for r in range(4):
  M=1 if r==0 else ell**r
  for L in (1,2,3):
   for offset in (5,7):
    mod=ell*ell
    count=sum((12*t+offset)%ell!=0 and (M-L*(12*t+offset))%ell!=0 for t in range(mod))
    expect=Fraction(ell-(2 if r==0 else 1),ell)
    check(Fraction(count,mod)==expect,'bad_local_exact_target_correction')

# Exact rectangular/domain and target-correction identities supporting the
# averaged drift proof. These do NOT test an asymptotic or its rate.
for X in (80,160,320,640):
 targets=list(range(X+1,2*X+1))
 targets=[N for N in targets if N%12==2]
 RN={N:sum(C4[a] and C3[N-a] for a in range(1,N)) for N in targets}
 M=sum(RN.values())
 rect=sum(1 for a in range(1,2*X) if C4[a] and a%12==1 for b in range(1,2*X-a+1) if C3[b] and b%12==1 and X<a+b<=2*X)
 check(M==rect,'exact_averaged_rectangle')
 for p in (5,17,29,41):
  AP={N:sum(C4[a] and C3[N-a] and a%p==0 for a in range(1,N)) for N in targets}
  scaled=sum(1 for m in range(1,(2*X-1)//p+1) if C4[m] and m%12==5 for b in range(1,2*X-p*m+1) if C3[b] and b%12==1 and X<p*m+b<=2*X)
  check(sum(AP.values())==scaled,'exact_power_cofactor_domain')
  signed=sum(Fraction(AP[N])-Fraction(0 if N%p==0 else RN[N],p-1) for N in targets)
  correction=sum(Fraction(RN[N],p-1) for N in targets if N%p==0)
  check(signed==sum(AP.values())-Fraction(M,p-1)+correction,'exact_averaged_target_correction')
  check(correction>=0,'target_correction_nonnegative')

# Negative controls are rejected deliberately, not counted as proof.
check(not C[25] and not C[49] and C4[25] and C3[49] and 25+49==74,'residue_does_not_imply_intersection')
p=5
check(Fraction(1,p-1)!=Fraction(2*p-1,p*p-p+1),'weighted_density_is_different')
check(not C3[5] and C3[1],'bad_content_drop_rejected')
# An unsifted component may contain four forbidden factors, so the
# two-prime conclusion genuinely needs z^3>N and the lower-prime sieve.
check(sum(F[625].values())==4 and C4[625],'unsifted_four_prime_control')

pins={}
for name in SOURCE_PATHS:
 path=ROOT/name;pins[name]=hashlib.sha256(path.read_bytes()).hexdigest()
result={'result':'PASS','scope':'finite exact diagnostics only; no analytic asymptotic, UBV, eventual coverage, cutoff or external review certified','limit':LIMIT,'checks':counters,'total_checks':sum(counters.values()),'source_sha256':pins,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
if args.output is not None:
 args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2))
