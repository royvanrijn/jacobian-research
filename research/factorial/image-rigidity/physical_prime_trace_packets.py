#!/usr/bin/env python3
"""Exact supplementary controls for the written generic physical trace-packet proof."""
from pathlib import Path
import hashlib,json,math,platform
HERE=Path(__file__).resolve().parent
COUNTS={}
def check(key,ok=True,n=1):
 assert ok,key
 COUNTS[key]=COUNTS.get(key,0)+n
def trim(a,m):
 a=[int(x)%m for x in a]
 while len(a)>1 and not a[-1]:a.pop()
 return a or [0]
def add(a,b,m,scale=1):
 c=[0]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=scale*x
 return trim(c,m)
def mul(a,b,m):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return trim(c,m)
def rem(a,q,m):
 a=trim(a,m)
 assert q[-1]%m==1
 while len(a)>=len(q):
  d=len(a)-len(q);z=a[-1]
  for j,c in enumerate(q):a[d+j]=(a[d+j]-z*c)%m
  a=trim(a,m)
 return a
def power(a,k,q,m):
 out=[1]
 while k:
  if k&1:out=rem(mul(out,a,m),q,m)
  a=rem(mul(a,a,m),q,m);k//=2
 return out
def chi(p,z,m):
 out=[math.comb(p,j)*(-1)**(p-j)%m for j in range(p+1)]
 out[p-1]=(out[p-1]-z)%m
 return out
def order(a,p):
 t=1
 for d in range(1,p):
  t=t*a%p
  if t==1:return d
 raise AssertionError((a,p))
def prime(p):return p>=2 and all(p%d for d in range(2,math.isqrt(p)+1))
PRIMES=[p for p in range(3,102,2) if prime(p)]
check("all_odd_primes_through_101",PRIMES==[3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101])
# Every residue class and one principal-unit perturbation at the first five odd primes.
# These replay the finite-free algebra step, including the a=1 Artin-Schreier case.
period_cases=[]
for p in (3,5,7,11,13):
 for residue in range(1,p):
  a=(-residue)%p;d=order(a,p) if a!=1 else p
  for lift in (0,1):
   z=residue+p*lift;m=p*p;q=chi(p,z,m)
   x=[0,1];lhs=power(x,p**(d+1),q,m);rhs=power(x,p,q,m)
   check("full_mod_p_squared_unit_period_identity",lhs==rhs)
   residue_q=chi(p,z,p)
   check("residue_affine_Frobenius_identity",power(x,p**d,residue_q,p)==x)
   # Adding p must be visible modulo p^2 even though it is invisible modulo p.
   check("corrupted_precision_rejected",rem(add(add(lhs,rhs,m,-1),[p],m),q,m)!=[0])
  period_cases.append({"p":p,"Z_residue":residue,"depth_period":d})
def traces_and_derivatives(p,z,last,m):
 # Newton sums and their derivatives with respect to z. Only c_1=-p-z varies.
 coefficients=list(reversed(chi(p,z,m)[:-1]))
 vals=[p%m];derivs=[0]
 for n in range(1,last+1):
  upto=min(n-1,p)
  value=-sum(coefficients[k-1]*vals[n-k] for k in range(1,upto+1))
  deriv=-sum(coefficients[k-1]*derivs[n-k] for k in range(1,upto+1))
  if n>=2:deriv+=vals[n-1]
  if n<=p:
   value-=n*coefficients[n-1]
   if n==1:deriv+=1
  vals.append(value%m);derivs.append(deriv%m)
 return vals,derivs
constant_cases=[]
for p in PRIMES:
 m=p*p;last=p+(p-1)*(p-1)
 vals,derivs=traces_and_derivatives(p,1,last,m)
 # A whole family of constants is checked by independently generated Newton sums.
 for b in range(1,p):
  n=p+b*(p-1);binomial=(p+1)*math.comb(p,b)%m
  unit=(p*((-1)**(b-1)*pow(b,-1,p)%p))%m
  check("exact_uniform_binomial_trace_constants",vals[n]==binomial==unit)
  check("special_trace_derivative_vanishes",derivs[n]%p==0)
  check("corrupted_prime_digit_constant_rejected",vals[n]!=(unit+p)%m)
 constant_cases.append({"p":p,"digits_checked":p-1})
 # Direct Taylor-invariance evaluations, independent of the derivative calculation.
 for e in (1,2):
  lifted,_=traces_and_derivatives(p,1+p*e,last,m)
  for b in range(1,p):
   n=p+b*(p-1)
   check("direct_principal_unit_trace_invariance",lifted[n]==vals[n])
def actual_H_point(p,n,y):
 # Direct integer coefficient definition, not the trace congruence.
 total=p;product=1;ppower=1
 for k in range(1,n+1):
  product*=math.prod(range(n+(p-1)*(k-1),n+(p-1)*k))
  numerator=math.comb(n,k)*product
  check("actual_physical_coefficient_integrality",numerator%ppower==0)
  coefficient=numerator//ppower
  total+=coefficient*y**k
  if k==n:check("actual_physical_leading_unit",coefficient%p!=0)
  ppower*=p
 return total
actual_cases=[]
for p in (3,5,7):
 for b in range(1,p):
  for r in (0,1,2):
   n=b*p**r+p-b;m=p*p
   value=actual_H_point(p,n,-1)
   vals,_=traces_and_derivatives(p,-math.factorial(p-1),n,m)
   check("direct_H_and_trace_packet_agree",value%m==vals[n])
   if r%2:
    expected=p*((-1)**(b-1)*pow(b,-1,p)%p)
    check("odd_packet_has_exact_valuation_one",value%m==expected and expected%p==0 and expected%m!=0)
   else:check("even_packet_is_unit",value%p==1)
   actual_cases.append({"p":p,"b":b,"r":r,"N":n,"H_at_minus_one_mod_p_squared":value%m})
# The valuation-one condition cannot be weakened merely to divisibility by p.
# (Y+1+p)^2 has a double root in the residue disk and value p^2 at -1.
for p in (3,5,7,11):
 F=[(1+p)**2,2*(1+p),1]
 value=F[0]-F[1]+F[2];root=-1-p
 check("valuation_two_repeated_root_guard",value==p*p and
       F[0]+F[1]*root+F[2]*root*root==0 and F[1]+2*F[2]*root==0)
imports={
 "PHYSICAL_PRIME_TOTAL_DEGREE_REDUCTION.tex":"42c5fc615400a49dd3ae4f19a9a832381607e69ce4f6ff6d255b57ca4cfd156f",
 "GAP_SIX_SMALL_OFFSET_FAMILIES.tex":"1b83ed6d1ea139e40de35f2b236d179e743bbe2202b555deea5dad547db3620b",
}
for name,digest in imports.items():
 check("frozen_mathematical_imports",hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest)
names=["PHYSICAL_PRIME_TRACE_PACKETS.tex",Path(__file__).name]+list(imports)
out={"passed":True,"kind":"internal exact supplementary replay; all quantifiers are proved in the written source",
 "scope":{"odd_primes":"every odd prime in written theorems; finite checker primes through101",
  "unit_point_depth_period":True,"uniform_prime_digit_constant":True,"packet":"N=b*p^r+p-b,1<=b<p,r>=0",
  "entire_residue_disk_including_ramified_points":True,"local_exclusion_only":True,
  "global_SFC":False,"Gamma_shape_extension":False,"indispensable_finite_certificate":False},
 "controls":COUNTS,"total_controls":sum(COUNTS.values()),"finite_period_profiles":period_cases,
 "finite_prime_digit_profiles":constant_cases,"direct_physical_packet_profiles":actual_cases,
 "source_sha256":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
 "python":platform.python_version(),"arithmetic":"pure Python integer lists modulo p and p^2; no polynomial factorization or Groebner input",
 "independent_replay":False,"formal_verification":False,"external_review":False,
 "novelty_claim":False,"priority_claim":False,
 "command":".venv/bin/python research/factorial/image-rigidity/physical_prime_trace_packets.py"}
(HERE/"physical_prime_trace_packets_verification.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps({"passed":True,"total_controls":sum(COUNTS.values()),"controls":COUNTS},indent=2))
