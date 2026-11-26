#!/usr/bin/env python3
"""Supplementary exact floor/supporting-line controls, not a finite-prefix universal proof."""
if not __debug__:raise SystemExit('Ordinary Python execution required.')
import argparse
from fractions import Fraction as F
from pathlib import Path
import json,hashlib

def vf(n,p):
 out=0
 while n:n//=p;out+=n
 return out

def v(n,p):
 out=0
 while n%p==0:n//=p;out+=1
 return out

def prefixes(q,p):
 ds=[];x=q;s=0
 while x:
  d=x%p
  if d:ds.append((s,d))
  x//=p;s+=1
 H=0;out=[]
 for s,d in reversed(ds):
  out.append((s,H,H+d*p**s));H+=d*p**s
 return out

def line(q,p,s,H,j):return -F(vf(H,p))-F(p**s-1,(p-1)*p**s)*(j-H)
def term(q,p,j,k):return k*v(q,p)+vf(2*q-j-k,p)-vf(q,p)-vf(q-j,p)-vf(k,p)-vf(j-k,p)
def require(c,msg):
 if not c:raise ValueError(msg)

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--output',type=Path,help='Write a fresh supplementary receipt; never overwrite an existing file.')
args=ap.parse_args()
floors=0;negT=0;rejects=[]
for P in [3,5,7,9,11,13,25]:
 for u in range(P):
  for a in range(P):
   for b in range(P):
    T=(2*u-a-b)//P+(a>u)+(b>a)
    require(T>=-1,'T lower bound')
    require(T!=-1 or a+b>P,'negative T boundary')
    require(a>u or T>=0,'equal high quotient boundary')
    require(F(T)+F(a+b,P)>=0,'low fractional defect')
    floors+=1;negT+=T==-1
for p in [3,5,7,11,13]:
 for k in range(2,101):
  require(F(k-k//p-1)-F(k,p*(p-1))>0,'low surplus pays entire tail')
  require(k-k//p-1>=0,'discarded low levels are nonnegative')
for P in [3,5,9,25,125]:
 for a in range(P):
  T=(-a-1)//P+(a>0)+(1>a)
  require(T==0,'k1 exact low T')
# Direct factorial-floor valuation controls, including a case outside the old strict height gate.
full_cases=[];terms=0
for q,p in [(3,3),(12,3),(39,3),(55,5),(135,5),(635,5)]:
 t=v(q,p)
 for j in range(q+1):
  for k in range(1,j+1):
   z=term(q,p,j,k)
   for s,H,H1 in prefixes(q,p):require(F(z)>line(q,p,s,H,j),'strict ordinary facet inequality')
   terms+=1
 full_cases.append({'q':q,'p':p,'terms':q*(q+1)//2,'facet_positions':[s for s,H,H1 in prefixes(q,p)]})
# Huge counterfamily values are handled by exact floor divisions only, without giant factorials.
q=5**19+10;p=5
huge=[]
for j,k in [(q-11,21),(q-11,22),(q-11,23),(5**19,22),(q,1),(q,2),(q,q)]:
 z=term(q,p,j,k);margins=[F(z)-line(q,p,s,H,j) for s,H,H1 in prefixes(q,p)]
 require(min(margins)>0,'huge family ordinary facet')
 huge.append({'j':j,'k':k,'ordinary_term_valuation':z,'minimum_facet_margin':str(min(margins))})
# False variants are actually rejected by finite witnesses, not merely listed.
for remove,label in [(0,'deleted_v_gt_u_indicator'),(1,'deleted_w_gt_v_indicator')]:
 witness=None
 for P in [3,5,7]:
  for u in range(P):
   for a in range(P):
    for b in range(P):
     T=(2*u-a-b)//P+(0 if remove==0 else a>u)+(0 if remove==1 else b>a)
     if T<-1 or (T==-1 and a+b<=P):witness=[P,u,a,b,T];break
    if witness:break
   if witness:break
  if witness:break
 require(witness is not None,'missing genuine corruption witness');rejects.append({'corruption':label,'witness':witness})
require(term(3,3,3,3)-3*v(3,3)<=line(3,3,1,0,3),'missing-kt corruption must fail')
rejects.append({'corruption':'deleted_kt_valuation','q':3,'p':3,'j':3,'k':3})
out={'status':'PASS','scope':'Supplementary finite exact floor and supporting-line controls; universal proof is the separately written all-range argument. No huge factorial expanded.','floor_remainder_cases':floors,'negative_T_cases':negT,'direct_factorial_floor_terms':terms,'full_cases':full_cases,'huge_family_floor_only_cases':huge,'rejected_corruptions':rejects,'actual_engine_cases':0,'canonical_mutations':False,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
encoded=json.dumps(out,indent=2)+'\n'
if args.output:
 with args.output.open('x',encoding='utf-8') as f:f.write(encoded)
print(encoded,end='')