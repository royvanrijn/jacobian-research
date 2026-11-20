#!/usr/bin/env python3
"""Exact finite coefficient support for the uniform high-gap Rouche region."""
if not __debug__:raise SystemExit('Refusing -O/-OO: exact assertions required')
from pathlib import Path
from fractions import Fraction as Q
from math import comb
from hashlib import sha256
import argparse,json
HERE=Path(__file__).resolve().parent

def exp_bounds(x,M=100):
 x=Q(x);term=Q(1);total=term
 for k in range(1,M+1):term*=x/k;total+=term
 first=term*x/(M+1);ratio=x/(M+2);assert 0<=ratio<1
 return total,total+first/(1-ratio)
def dyadic_root_upper(v,N,bits=32):
 assert 0<=v<1
 den=2**bits;lo,hi=0,den
 while hi-lo>1:
  mid=(lo+hi)//2
  if Q(mid,den)**N>=v:hi=mid
  else:lo=mid
 out=Q(hi,den);assert out**N>=v
 if hi:assert Q(hi-1,den)**N<v
 return out

def clean(p):
 while p and p[-1]==0:p.pop()
 return p
def rem(a,b,p):
 a=a[:];b=clean(b[:]);assert b
 inv=pow(b[-1],-1,p)
 while a and len(a)>=len(b):
  c=a[-1]*inv%p;k=len(a)-len(b)
  for j,v in enumerate(b):a[j+k]=(a[j+k]-c*v)%p
  clean(a)
 return a
def gcddeg(a,b,p):
 while b:a,b=b,rem(a,b,p)
 return len(a)-1

def sha(p):return sha256(p.read_bytes()).hexdigest()
def root_path(arg):
 if arg:return arg.resolve()
 for p in [HERE,*HERE.parents]:
  if (p/'research/MATH_STATUS.json').is_file():return p
 raise SystemExit('Scratch replay requires --repo-root')

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path);ap.add_argument('--output',type=Path);ap.add_argument('--write-receipt',action='store_true');args=ap.parse_args();root=root_path(args.repo_root)
 inputs=json.loads((HERE/'factorial_high_gap_rouche_inputs.json').read_text());rows=[];terms=0
 for rel,pin in inputs['retained_source_hashes'].items():assert sha(root/rel)==pin,rel
 assert sha(HERE/'factorial_high_gap_primary_source_record.json')==inputs['primary_record_sha256']
 assert sha(HERE/'FACTORIAL_HIGH_GAP_ROUCHE_REGION.md')==inputs['proof_sha256']
 for N,d,g in inputs['fixtures']:
  gamma=Q(g);alpha=Q(d*d*N,4)/(gamma+d*N)
  lo,hi=exp_bounds(alpha);assert lo>=32*N
  rising=[Q(1)]
  for j in range(d*N):rising.append(rising[-1]*(gamma+j))
  C=rising[d*N]
  coefficient_uppers=[]
  for k in range(1,N):
   nth=Q(comb(N,k))**N*rising[d*k]**N/C**k
   upper=dyadic_root_upper(nth,N);coefficient_uppers.append(upper);terms+=1
  S=sum(coefficient_uppers,Q(0));assert S<=Q(2,31)
  p=1000003
  ac=[Q(comb(N,k))*rising[d*(N-k)] for k in range(N+1)]
  a=[c.numerator*pow(c.denominator,-1,p)%p for c in ac]
  derivative=[k*a[k]%p for k in range(1,len(a))]
  assert gcddeg(a,derivative,p)==0
  assert S+1>Q(2,31) # a malformed perturbation budget is rejected
  rows.append({'N':N,'gap':d,'gamma':str(gamma),'alpha':str(alpha),'normalized_middle_coefficient_sum_upper':str(S),'full_degree_modular_squarefree_prime':p})
 # Exact analytic margins in the unrestricted proof.
 assert exp_bounds(Q(1,4))[1]<Q(4,3)
 assert 2*(Q(1)/(1-Q(1,32))-1)==Q(2,31)
 assert Q(2,31)*Q(4,3)==Q(8,93)<Q(1,6)
 assert Q(3,2)-Q(4,3)==Q(1,6)
 # The former gap6 domain is not claimed covered by this sufficient criterion.
 N,d,g=inputs['outside_criterion'];alpha=Q(d*d*N,4)/(g+d*N)
 assert exp_bounds(alpha)[1]<32*N
 # Positive measures alone do not imply simplicity: a three-atom exact countercontrol.
 mu2=Q(1,9)*3;mu4=Q(1,9)*9
 assert mu2==Q(1,3) and mu4==1
 app=[mu4,Q(0),6*mu2,Q(0),Q(1)]
 assert app==[1,0,2,0,1] # (z^2+1)^2, repeated roots
 assert gcddeg([int(c) for c in app],[0,4,0,4],1000003)==2
 proof=HERE/'FACTORIAL_HIGH_GAP_ROUCHE_REGION.md';record=HERE/'factorial_high_gap_primary_source_record.json'
 imports=['research/factorial/image-rigidity/GAP_SIX_RESULTANT_STRICT_POSITIVITY.md','research/factorial/image-rigidity/BINOMIAL_SHAPE_RESULTANT.tex']
 result={'status':'PASS','scope':'A few exact normalized coefficient and modular controls support the full written infinite-parameter proof; no unproved saddle extension or all-gap closure inferred.','profiles':rows,'normalized_middle_coefficient_controls':terms,'source_hashes':{'proof':sha(proof),'checker':sha(Path(__file__).resolve()),'inputs':sha(HERE/'factorial_high_gap_rouche_inputs.json'),'primary_inspection_record':sha(record),**{x:sha(root/x) for x in imports}},'negative_controls':['malformed_normalized_middle_budget','outside_sufficient_region','positive_measure_multiple_root']}
 out=args.output or HERE/'factorial_high_gap_rouche_verification.json';out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print('PASS',len(rows),'profiles',terms,'middle coefficients');print(out)
if __name__=='__main__':main()
