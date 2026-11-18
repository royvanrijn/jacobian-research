#!/usr/bin/env python3
"""Reproducible standard-library finite exact cofactor diagnostics.

Run without Python -O from any directory. No output file is written unless
--output is supplied. These controls do not verify the analytic proof.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
from math import comb,gcd,isqrt,prod
import argparse,hashlib,json,sys
if not __debug__:
    raise RuntimeError("Exact controls require assertions; run without Python -O.")
parser=argparse.ArgumentParser(description="Finite exact asymmetric arithmetic controls; not an analytic proof verifier.")
parser.add_argument("--output",type=Path,help="Optional JSON receipt; otherwise print only.")
args=parser.parse_args()
ROOT=Path(__file__).resolve().parents[3]
def primes(n):
    flags=bytearray([1])*(n+1);flags[:2]=b'\0\0'
    for p in range(2,isqrt(n)+1):
        if flags[p]:flags[p*p:n+1:p]=b'\0'*((n-p*p)//p+1)
    return [p for p in range(2,n+1) if flags[p]]
ps=primes(1500)
def factor(n):
    out={}
    for p in ps:
        if p*p>n:break
        while n%p==0:out[p]=out.get(p,0)+1;n//=p
    if n>1:out[n]=out.get(n,0)+1
    return out
facts=[{}]+[factor(n) for n in range(1,1501)]
def U(n):return prod((e+1)**2 for e in facts[n].values()) if all(p%12==1 for p in facts[n]) else 0
def V(n):return prod(e+1 for e in facts[n].values()) if all(p%3==1 for p in facts[n]) else 0
local=0
for p in primes(300):
    if p<5:continue
    x=F(1,p);q=1-x;E=F(p-1,p-2);a=F(p-1,p*(p-2))
    assert (q*x)/(1-2*x)==a
    if p%12==1:
        for e in range(1,31):
            previous=1 if e==1 else E
            total=E+x*previous
            assert total>=1+2*a
            assert total-(1+2*a)==(0 if e==1 else F(1,p*(p-2)))
            assert x*previous/total<=x
            local+=1
    if p%12==7:
        assert E-(1+a)==F(1,p*(p-2));local+=1
    assert F(2,3)/(F(1,3)*F(1,3))==6

def upper(n,h):return sum((-1)**k*comb(n,k) for k in range(min(h,n)+1))
signcases=0
for n in product(range(11),repeat=3):
    u=[upper(n[j],2*(1+j)) for j in range(3)]
    l=[upper(n[j],2*(1+j)+1) for j in range(3)]
    value=prod(u)-sum((u[j]-l[j])*prod(u[i] for i in range(3) if i!=j) for j in range(3))
    event=int(all(k==0 for k in n))
    assert all(k>=0 for k in u) and value<=event<=prod(u)
    signcases+=1
assert upper(4,3)**2==1 # Deliberately invalid product of lower truncations.

def elementary(gs):
    es=[F(1)]
    for g in gs:
        es.append(F(0))
        for k in range(len(es)-1,0,-1):es[k]+=g*es[k-1]
    return es
blocks=[[F(1,10+j)]*20 for j in range(3)]
models=[]
for j,gs in enumerate(blocks):
    es=elementary(gs);h=2*(4+j)
    eu=sum((-1)**k*es[k] for k in range(min(h,len(gs))+1))
    el=eu-(es[h+1] if h+1<len(es) else 0)
    true=prod(1-g for g in gs);models.append((eu,el,true))
pred=prod(t[0] for t in models)-sum((models[j][0]-models[j][1])*prod(models[i][0] for i in range(3) if i!=j) for j in range(3))
true=prod(t[2] for t in models)
assert pred>=true/2

# Generate all admitted tuples and residual pairs directly, then check independent recovery.
z=17;w=13;small=[p for p in ps if 5<=p<=z]
seen={};fibres=0;high_target=0
for N in range(2,1501,6):
  for d in [1,13]:
    if N%d:continue
    M=N//d
    for u in [1,13]:
      for v in [1,7,13,91]:
        if gcd(d,u*v)>1 or gcd(u,v)>1 or gcd(u*v,N)>1:continue
        directives=set(facts[u])|set(facts[v])
        for x in range(1,M):
          y=M-x
          if not U(x) or not V(y):continue
          if any(facts[x].get(p,0)!=1 for p in facts[u]):continue
          if any(facts[y].get(p,0)!=1 for p in facts[v]):continue
          if any(x%p==0 or y%p==0 for p in small if p not in directives):continue
          a=d*x;b=d*y
          key=(a,b)
          assert key not in seen,(key,seen.get(key),(d,u,v))
          seen[key]=(d,u,v)
          rd=ru=rv=1
          for p in [7,13]:
            ea=facts[a].get(p,0);eb=facts[b].get(p,0)
            if ea and eb:
              assert ea==eb==1 and p%12==1
              rd*=p
            elif ea:
              assert ea==1 and p%12==1
              ru*=p
            elif eb:
              assert eb==1 and p%3==1
              rv*=p
          assert (rd,ru,rv)==(d,u,v)
          assert F(U(x)*V(y),4**len(facts[u])*2**len(facts[v]))==U(x//u)*V(y//v)
          assert all(p>z for p in facts[x//u]) and all(p>z for p in facts[y//v])
          if N%169==0 and d==13:high_target+=1
          fibres+=1
assert high_target>0
paths=['research/arithmetic/ASYMMETRIC_SEMIGROUP_SCALE.md','research/arithmetic/ASYMMETRIC_UNIFORM_COEFFICIENTS.md','research/arithmetic/ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md']
result={'scope':'Finite exact local compensation, block-sign and fibre diagnostics only; no uniform analytic estimate or unrestricted count certified','local_compensation_cases':local,'block_sign_cases':signcases,'negative_control_product_of_lower_polynomials':1,'bernoulli_lower_model_at_least_half':True,'fibre_cases':fibres,'high_target_13_squared_common_d_cases':high_target,'source_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}}
result['checker']={'path':Path(__file__).relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version.split()[0]}
result['result']='PASS'
serialized=json.dumps(result,indent=2)+'\n'
if args.output is not None:
    args.output.write_text(serialized)
print(serialized,end='')
