"""Bounded exact controls of the all-prime proof; no noncollision claim."""
from pathlib import Path
from math import comb
import json,hashlib
HERE = Path(__file__).resolve().parent
COUNTS={}
def check(k,ok):
 if not ok:raise AssertionError(k)
 COUNTS[k]=COUNTS.get(k,0)+1
def vp(n,p):
 if not n:return 10**9
 v=0
 while n%p==0:n//=p;v+=1
 return v
def vpf(n,p):
 a=0
 while n:n//=p;a+=n
 return a
def sp(n,p):
 s=0
 while n:s+=n%p;n//=p
 return s
def raise_poly(f,c):
 g=[0]*(len(f)+1)
 for k,x in enumerate(f):g[k]+=(c-k)*x;g[k+1]+=x
 return g
def binom(n,k):return comb(n,k) if n>=0 else (-1)**k*comb(k-n-1,k)
profiles=[]
for p in (2,3,5,7,11):
 q=p-1
 for N in sorted(set((1,2,p-1,p,p+1,2*p-1,p*p-1))):
  C=N+q-1;e=[(sp(N,p)-a+q-2)//q for a in range(q)]
  f=[1]
  for j in range(N):f=raise_poly(f,C+q*j)
  for r in range(13):
   den=vpf(N,p)+vpf(r,p)
   for a in range(q):
    for j in range(a,len(f),q):
     check('actual_integer_valued_binomial_column_content',
           vp(f[j],p)+(j-a)//q-den>=e[a])
   f=raise_poly(f,C-r)
  profiles.append({'p':p,'N':N,'binomial_degrees':13})
 for j in range(1,1001):
  ell,r=divmod(j,q)
  for eps in (0,1):
   val=ell+eps-vpf(j,p)
   check('formal_gauge_coefficient_integral',val>=0)
   check('formal_gauge_digit_formula',val==eps+(sp(j,p)-r)//q)
 for depth in range(1,6):
  j=p**depth
  if p==2:check('dyadic_nonrestricted_tail',j-vpf(j,p)==1)
  else:check('odd_prime_unit_tail',(j-1)//q-vpf(j,p)==0)
 for digit_sum in range(1,151):
  e=[(digit_sum-a+q-2)//q for a in range(q)]
  weights=[1+e[-1]-e[0]]+[e[a-1]-e[a]for a in range(1,q)]
  check('exact_one_weighted_link',sum(weights)==1 and set(weights)<={0,1})
 # Explicit class-indicator approximants, including negative integer inputs.
 for r in (1,2):
  modulus=p**r
  for a in sorted(set((0,1,p-1,p%modulus,modulus-1))):
   for x in sorted(set((-modulus-1,-1,0,1,a,a+modulus,a-modulus,modulus+3))):
    rho=1
    for i in range(r):
     digit=(a//p**i)%p
     rho*=1-(binom(x,p**i)-digit)**(p-1)
    check('Lucas_class_indicator',rho%p==int(x%modulus==a))
    for m in range(3):
     if x%modulus==a:
      check('uniform_on_class_projector_bound',pow(rho,p**m,p**(m+1))==1)
     else:
      check('uniform_off_class_projector_bound',pow(rho,p**m,p**(p**m))==0)
# Distinguish the dyadic tail valuation1 from an incorrect unit-tail extension.
check('reject_wrong_dyadic_unit_tail',2**5-vpf(2**5,2)!=0)
files=[Path(__file__).resolve(),HERE/'PHYSICAL_PRIME_RAISING_CALCULUS.tex',
 HERE/'ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex']
out={'status':'PASS','controls':COUNTS,'total_controls':sum(COUNTS.values()),
 'profiles':profiles,
 'scope':'Bounded normalization/formal-coefficient/projector controls; unbounded all-prime calculus proof in TeX. No noncollision or finite spectral closure claim.',
 'source_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in files},
 'independent_whole_proof_replay':False,'formal_verification':False,'external_review':False}
(HERE/'physical_prime_raising_calculus_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'PASS','total_controls':sum(COUNTS.values()),'controls':COUNTS},indent=2))
