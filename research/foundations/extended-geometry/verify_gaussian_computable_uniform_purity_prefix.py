#!/usr/bin/env python3
"""Finite exact controls for the Gaussian uniform-prefix argument."""
if not __debug__:raise RuntimeError('Refusing -O: exact assertions are required')
import argparse,json,hashlib
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
from collections import defaultdict
HERE=Path(__file__).resolve().parent

def clean(p):return {e:c for e,c in p.items() if c}
def add(*ps):
 out=defaultdict(Q)
 for p in ps:
  for e,c in p.items():out[e]+=c
 return clean(out)
def mul(p,q):
 out=defaultdict(Q)
 for a,c in p.items():
  for b,d in q.items():out[tuple(x+y for x,y in zip(a,b))]+=c*d
 return clean(out)
def scale(p,c):return clean({e:v*c for e,v in p.items()})
def power(p,n,dim):
 out={(0,)*dim:Q(1)}
 for _ in range(n):out=mul(out,p)
 return out
def diff(p,i):
 out={}
 for e,c in p.items():
  if e[i]:
   v=list(e);v[i]-=1;out[tuple(v)]=e[i]*c
 return out
def wick(e):
 a=1
 for v in e:
  if v%2:return 0
  a*=factorial(v)//(2**(v//2)*factorial(v//2))
 return a

def controls():
 symbolic=0
 for n in range(1,5):
  lin={};quad={}
  for j in range(n):
   e=[0]*(2*n);e[j]=e[n+j]=1;lin[tuple(e)]=Q(1)
   v=[0]*n;v[j]=2;quad[tuple(v)]=Q(1)
  for m in range(13):
   actual=defaultdict(Q)
   for e,c in power(lin,m,2*n).items():actual[e[:n]]+=c*wick(e[n:])
   expected={} if m%2 else scale(power(quad,m//2,n),wick((m,)))
   assert clean(actual)==expected
   symbolic+=1
 # Literal source identities A(D_t+P)1 = -nabla(B).
 t={(1,0):Q(1)};x={(0,1):Q(1)};one={(0,0):Q(1)}
 source=0
 for deg,B in [(1,one),(2,x)]:
  P=power(x,deg,2)
  nabla=add(diff(B,1),scale(mul(x,B),-1),mul(mul(t,diff(P,1)),B))
  if deg==1:left=add(P,scale(t,-1))
  else:left=add(mul(add(one,scale(t,-2)),P),scale(one,-1))
  assert left==scale(nabla,-1);source+=1
 # Actual formal Gaussian EGF jets and exact scalar ODEs.
 jets=0
 for degree in (1,2):
  seq=[Q(wick((degree*m,)),factorial(m)) for m in range(33)]
  for k in range(30):
   if degree==1: lhs=(k+1)*seq[k+1]-(seq[k-1] if k else 0)
   else:lhs=(k+1)*seq[k+1]-(2*k+1)*seq[k]
   assert lhs==0;jets+=1
  assert seq[0]==1
  if degree==1:assert seq[1]==0 and seq[2]==Q(1,2)
  else:assert seq[1]==1
 # Isotropic complex Gaussian, performed without numerical complex arithmetic.
 for m in range(2,32,2):
  value=sum(Q(factorial(m),factorial(k)*factorial(m-k))*(-1)**(k//2)*wick((k,m-k)) for k in range(0,m+1,2))
  assert value==0
 return {'symbolic_centered_linear_moments':symbolic,'literal_unit_source_identities':source,'formal_egf_ode_coefficients':jets,'isotropic_complex_moments':15,'computed_numerical_unrestricted_prefix':False}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,default=None);ap.add_argument('--output',type=Path,default=HERE/'gaussian_computable_uniform_purity_prefix_verification.json');args=ap.parse_args()
 
 if args.repo_root is None:
  args.repo_root=next((p for p in [HERE,*HERE.parents]if(p/'research/MATH_STATUS.json').is_file()),None)
  if args.repo_root is None:raise SystemExit('Use --repo-root when copying this checker outside the repository.')
 counts=controls();ledgerpath=HERE/'shifted_gaussian_primary_source_ledger.json';ledger=json.loads(ledgerpath.read_text())
 checked=[]
 for r in ledger['files']:
  path=(args.repo_root if r['location']=='repository' else HERE)/r['path']
  actual=hashlib.sha256(path.read_bytes()).hexdigest();assert actual==r['sha256'],str(path)
  checked.append({'path':r['path'],'sha256':actual})
 out={'status':'PASS','scope':'Finite formula and source-identity controls supplement the unrestricted termination proof; no universal prime-stratification run or numerical B_G is claimed.','counts':counts,'verified_source_pins':checked,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'ledger_sha256':hashlib.sha256(ledgerpath.read_bytes()).hexdigest()}
 args.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
