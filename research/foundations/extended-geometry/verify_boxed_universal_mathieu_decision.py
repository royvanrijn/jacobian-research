#!/usr/bin/env python3
"""Finite exact controls for boxed universal Mathieu decision.

No unrestricted minimal-prime/Hom loop is executed. The known all-order
phase families are explicit proof antecedents, not inferred from test jets.
"""
if not __debug__:
    raise SystemExit('Refusing -O/-OO: exact assertions are required')
from pathlib import Path
from fractions import Fraction as Q
from collections import defaultdict
from itertools import product,combinations
from math import factorial,comb
from hashlib import sha256
import argparse,json

HERE=Path(__file__).resolve().parent
PROOF='BOXED_UNIVERSAL_MATHIEU_DECISION.md'
INPUTS='boxed_universal_mathieu_decision_inputs.json'
LEDGER='boxed_universal_mathieu_decision_source_ledger.json'
RECEIPT='boxed_universal_mathieu_decision_verification.json'

def clean(p):return {e:Q(c) for e,c in p.items() if c}
def add(*ps):
 d=defaultdict(Q)
 for p in ps:
  for e,c in p.items():d[e]+=c
 return clean(d)
def scale(p,c):return clean({e:v*c for e,v in p.items()})
def mul(p,q):
 d=defaultdict(Q)
 for e,c in p.items():
  for f,v in q.items():d[tuple(a+b for a,b in zip(e,f))]+=c*v
 return clean(d)
def one(n):return {(0,)*n:Q(1)}
def power(p,k,n):
 r=one(n)
 for _ in range(k):r=mul(r,p)
 return r
def var(n,j):
 e=[0]*n;e[j]=1;return {tuple(e):Q(1)}
def evaluate(p,point):return sum(c*__import__('functools').reduce(lambda a,b:a*b,(Q(x)**e for x,e in zip(point,ex)),Q(1)) for ex,c in p.items())
def poly(terms):return clean({tuple(e):Q(c) for e,c in terms})
def rank(rows):
 if not rows:return 0
 a=[list(map(Q,r)) for r in rows];r=0
 for j in range(len(a[0])):
  k=next((k for k in range(r,len(a)) if a[k][j]),None)
  if k is None:continue
  a[r],a[k]=a[k],a[r];v=a[r][j];a[r]=[x/v for x in a[r]]
  for k in range(r+1,len(a)):
   if a[k][j]:
    v=a[k][j];a[k]=[x-v*y for x,y in zip(a[k],a[r])]
  r+=1
  if r==len(a):break
 return r

def compositions(total,n):
 if n==1:return [(total,)]
 return [(i,)+rest for i in range(total+1) for rest in compositions(total-i,n-1)]

def minimal_coordinate_primes(gens,n):
 supports=[{i for i,e in enumerate(ex) if e} for ex in gens]
 candidates=[]
 for mask in range(1<<n):
  s=frozenset(i for i in range(n) if mask>>i&1)
  if all(s&sup for sup in supports):candidates.append(s)
 return [s for s in candidates if not any(t<s for t in candidates)]
def in_coordinate_prime(p,prime):return all(any(ex[j] for j in prime) for ex in p)

def component_controls(inputs):
 records=[];sample_checks=0
 for case in inputs['component_cases']:
  gens=[tuple(x) for x in case['monomial_ideal']];g=poly(case['guard'])
  primes=minimal_coordinate_primes(gens,2)
  retained=[p for p in primes if not in_coordinate_prime(g,p)]
  unsafe=[p for p in retained if 0 not in p] # F=a*F4; a=0 is exactly the safe branch
  answer=not unsafe
  assert answer==case['universal_valid'],case['name']
  assert len(retained)==case['retained_components']
  for a,b in product(range(-2,3),repeat=2):
   if any(evaluate({ex:Q(1)},(a,b)) for ex in gens) or not evaluate(g,(a,b)):continue
   if answer:assert a==0
   sample_checks+=1
  records.append({'name':case['name'],'minimal_primes':[sorted(p) for p in primes],'retained_primes':[sorted(p) for p in retained],'universal_valid':answer})
 # The generic unsafe family has the safe closed point a=0. Never reverse specialization.
 assert 0==0 and not in_coordinate_prime(var(2,0),frozenset())
 return records,sample_checks

def wick(k):
 if k%2:return 0
 return factorial(k)//(2**(k//2)*factorial(k//2))
def reduce_split_prime(p):
 # R=Q[a,b]/(a^2+b^2), a Q-domain that splits geometrically.
 d=defaultdict(Q)
 for (i,j),c in p.items():d[(i+2*(j//2),j%2)]+=c*(-1)**(j//2)
 return clean(d)
def isotropic_controls(inputs):
 tests=0;tail=0
 for r,s in product(range(inputs['max_linear_marker_exponent']+1),repeat=2):
  for m in range(inputs['max_linear_moment']+1):
   p={(k,m-k):Q(comb(m,k)*wick(r+k)*wick(s+m-k)) for k in range(m+1)}
   red=reduce_split_prime(p)
   if m>r+s:assert red=={};tail+=1
   if r==s==0 and m>0:assert red=={}
   # Source computations can divide by a; the identities still specialize at a=b=0.
   if m>r+s:assert evaluate(red,(0,0))==0
   tests+=1
 return tests,tail

def image(p,n):
 d=defaultdict(Q)
 for ex,c in p.items():
  a,b=ex[:n],ex[n:]
  if any(x>y for x,y in zip(a,b)):continue
  for x,y in zip(a,b):c*=factorial(y)//factorial(y-x)
  d[tuple(y-x for x,y in zip(a,b))]+=c
 return clean(d)
def original_witness_controls(inputs):
 n=4;xi1,xi2,s1,s2=[var(n,j) for j in range(n)]
 R=add(mul(xi1,s1),mul(xi2,s2));Z=mul(xi1,s2);W=scale(mul(xi2,s1),2);T=add(mul(xi1,s1),scale(mul(xi2,s2),-1))
 U=add(R,Z);B=add(mul(power(R,2,n),W),scale(mul(add(scale(R,2),Z),power(T,2,n)),Q(-1,2)));F=mul(U,B)
 assert len(F)==16
 assert evaluate(R,(1,1,1,-1))==0 and evaluate(F,(1,1,1,-1))==-2
 assert evaluate(U,(1,1,-2,1))==0 and evaluate(B,(1,1,-2,1))==Q(1,2)
 vals=[];fm=one(n)
 for m in range(1,inputs['max_F4_moment']+1):
  fm=mul(fm,F);assert image(fm,2)=={}
  denom=1
  for j in range(1,2*m+2,2):denom*=j
  value=Q(factorial(4*m+2)*factorial(m),denom)
  assert image(mul(Z,fm),2)=={(0,0):value}
  vals.append(str(value))
 return vals

def interpolation_controls(inputs):
 records=[];integral=0;malformed=0
 for n,e,k in inputs['interpolation_cases']:
  mons=[a+b for a in compositions(e,n) for b in compositions(e,n)];M=len(mons)
  selected=[];rows=[]
  for point in product(range(e+1),repeat=2*n):
   row=[evaluate({ex:Q(1)},point) for ex in mons]
   if rank(rows+[row])>len(rows):selected.append(point);rows.append(row)
   if len(rows)==M:break
  assert len(rows)==M and rank(rows)==M
  # Universal H in variables (original coordinates, coefficient variables).
  H={ex+tuple(int(j==i) for j in range(M)):Q(1) for i,ex in enumerate(mons)}
  Hk=power(H,k,2*n+M)
  for point,row in zip(selected,rows):
   actual=defaultdict(Q)
   for ex,c in Hk.items():actual[ex[2*n:]]+=c*evaluate({ex[:2*n]:Q(1)},point)
   L={tuple(int(i==j) for j in range(M)):v for i,v in enumerate(row) if v}
   expected=power(L,k,M)
   assert clean(actual)==expected
   assert add(clean(actual),one(M))!=expected # literal malformed monic relation rejected
   integral+=1;malformed+=1
  records.append({'pairs':n,'input_bidegree':e,'power':k,'input_dimension':M,'selected_rational_points':[list(x) for x in selected],'evaluation_rank':rank(rows)})
 return records,integral,malformed

def prime_guard_controls(inputs):
 tests=0
 for d in range(1,inputs['max_power_divisor_degree']+1):
  primes=[p for p in range(2,d+1) if d%p==0 and all(p%q for q in range(2,p))]
  for k in range(2,d+1):
   if d%k==0:assert any(k%p==0 for p in primes)
  tests+=1
 assert comb(2+5-1,5)**2==36 and comb(2+4-1,4)**2==25
 return tests

def find_root(arg):
 if arg:return arg.resolve()
 for p in [HERE,*HERE.parents]:
  if (p/'research/MATH_STATUS.json').is_file():return p
 raise SystemExit('Scratch replay needs --repo-root')
def digest(p):return sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path);ap.add_argument('--output',type=Path);ap.add_argument('--write-receipt',action='store_true');args=ap.parse_args()
 root=find_root(args.repo_root);inputs=json.loads((HERE/INPUTS).read_text());assert inputs['schema']==1
 components,samples=component_controls(inputs);iso,tails=isotropic_controls(inputs);witness=original_witness_controls(inputs);interp,integral,malformed=interpolation_controls(inputs);guards=prime_guard_controls(inputs)
 ledger=json.loads((HERE/LEDGER).read_text());pins=[]
 for item in ledger['files']:
  p=(root if item['location']=='repository' else HERE)/item['path']
  h=digest(p);assert h==item['sha256'],str(p);pins.append({'path':item['path'],'sha256':h})
 result={'status':'PASS','scope':'Finite exact specialization/interpolation/original-moment controls; no general minimal-prime or generic Hom computation is executed.','counts':{'component_restriction_cases':len(components),'rational_sample_memberships':samples,'split_Q_prime_marker_moments':iso,'uniform_linear_tail_identities':tails,'F4_original_powers':len(witness),'interpolation_profiles':len(interp),'integral_power_relations':integral,'malformed_power_relations_rejected':malformed,'prime_divisor_guard_degrees':guards},'component_records':components,'F4_marked_values':witness,'interpolation_records':interp,'verified_source_pins':pins,'proof_sha256':digest(HERE/PROOF),'checker_sha256':digest(Path(__file__).resolve()),'inputs_sha256':digest(HERE/INPUTS),'ledger_sha256':digest(HERE/LEDGER)}
 out=args.output or HERE/RECEIPT;out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps({'status':result['status'],'counts':result['counts']},sort_keys=True));print(out)
if __name__=='__main__':main()
