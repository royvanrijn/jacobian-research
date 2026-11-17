#!/usr/bin/env python3
"""Seek rank-three tangent cubics p^2*q in S, deleting four proportional old columns."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,combinations_with_replacement
from collections import defaultdict
from math import gcd,lcm
from functools import reduce
import hashlib,json,time
from flint import fmpq_mat
from power_line_compress import cubic
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def plus(a,b,c):
 out=dict(a)
 for e,v in b.items():out[e]=out.get(e,Q(0))+c*v
 return {e:v for e,v in out.items() if v}
def normalize(v):
 den=lcm(*(c.denominator for c in v));w=[int(c*den) for c in v];g=reduce(gcd,w,0)
 if not g:return None,Q(0)
 if next(c for c in w if c)<0:g=-g
 return tuple(c//g for c in w),Q(g,den)
def solve(M,rhs,d):
 rr,rank=fmpq_mat([[str(c) for c in row]+[str(b)] for row,b in zip(M,rhs)]).rref();x=[Q(0)]*d
 for i in range(rr.nrows()):
  pivot=next((j for j in range(d) if rr[i,j]),None)
  if pivot is None:
   if rr[i,d]:return None
  else:x[pivot]=Q(str(rr[i,d]))
 return x
def main():
 out=HERE/'power-tangent-search';out.mkdir(exist_ok=True);source=HERE/'power-basis-anneal/candidate.json';d=json.loads(source.read_text());poolpath=HERE/'power-two-exchange/pool.json';pool=json.loads(poolpath.read_text())['forms'];B=[[Q(c) for c in row] for row in d['B']];n=89
 cubes=[cubic(f) for f in d['forms']];mons=sorted(set().union(*(set(f) for f in cubes)));mi={e:i for i,e in enumerate(mons)}
 R,_=fmpq_mat([[str(f.get(e,0)) for f in cubes]+[str(int(i==j)) for j in range(len(mons))] for i,e in enumerate(mons)]).rref();RR=fmpq_mat([[str(R[i,n+j]) for j in range(len(mons))] for i in range(len(mons))])
 groups=defaultdict(list)
 for j in range(n):
  v=tuple(row[j] for row in B);c=next(a for a in v if a);groups[tuple(a/c for a in v)].append((j,c))
 groups={v:ids for v,ids in groups.items() if len(ids)>=4};assert sum(len(list(combinations(ids,4))) for ids in groups.values())==253
 stats=dict(directions=0,nonzero_tangent_spaces=0,rank1_only=0,linear_systems=0,best=89);dimensions=defaultdict(int);started=time.monotonic();answer=None;complete=True
 for index,p in enumerate(pool):
  if time.monotonic()-started>60:complete=False;break
  stats['directions']+=1;square={}
  for i,j in combinations_with_replacement([a for a,c in enumerate(p) if c],2):
   e=[0]*18;e[i]+=1;e[j]+=1;square[tuple(e)]=(1 if i==j else 2)*p[i]*p[j]
  columns=[]
  for j in range(18):
   f={}
   for e,c in square.items():ee=list(e);ee[j]+=1;f[tuple(ee)]=c
   columns.append(f)
  V=RR*fmpq_mat([[str(f.get(e,0)) for f in columns] for e in mons])
  residual=[[Q(str(V[i,j])) for j in range(18)] for i in range(n,len(mons))]
  outside=sorted(set().union(*(set(f) for f in columns))-set(mi));residual.extend([[Q(f.get(e,0)) for f in columns] for e in outside])
  residual=[row for row in residual if any(row)]
  if residual:
   rr,rank=fmpq_mat([[str(c) for c in row] for row in residual]).rref();pivots=[next(j for j in range(18) if rr[i,j]) for i in range(rank)]
  else:rank=0;pivots=[];rr=None
  free=[j for j in range(18) if j not in pivots];dim=len(free);dimensions[dim]+=1
  if not dim:continue
  stats['nonzero_tangent_spaces']+=1;N=[[Q(0)]*dim for _ in range(18)]
  for a,j in enumerate(free):
   N[j][a]=1
   for k,pivot in enumerate(pivots):N[pivot][a]=-Q(str(rr[k,j]))
  if dim==1:
   q0=[row[0] for row in N];nonzero=next(j for j,c in enumerate(p) if c);ratio=q0[nonzero]/Q(p[nonzero])
   if all(a==ratio*b for a,b in zip(q0,p)):stats['rank1_only']+=1;continue
  A=[[sum((Q(str(V[j,k]))*N[k][a] for k in range(18)),Q(0)) for a in range(dim)] for j in range(n)]
  for direction,ids in groups.items():
   possible=[(j,c) for j,c in ids if any(A[j])]
   for selected in combinations(possible,4):
    stats['linear_systems']+=1;theta=solve([A[j] for j,c in selected],[c for j,c in selected],dim)
    if theta is None:continue
    q=[sum((c*t for c,t in zip(row,theta)),Q(0)) for row in N];a=[sum((c*t for c,t in zip(row,theta)),Q(0)) for row in A]
    removed=[j for j,c in ids if a[j]==c];assert len(removed)>=4
    byform={tuple(f):[row[j]-direction[o]*a[j] for o,row in enumerate(B)] for j,f in enumerate(d['forms'])}
    rationalforms=[[x+Q(y) for x,y in zip(q,p)],[x-Q(y) for x,y in zip(q,p)],q]
    for f,weight in zip(rationalforms,[Q(1,6),Q(1,6),Q(-1,3)]):
     form,factor=normalize(f)
     if form is None:continue
     coeff=weight*factor**3
     if form not in byform:byform[form]=[Q(0)]*18
     byform[form]=[x+coeff*y for x,y in zip(byform[form],direction)]
    newforms=sorted(f for f,col in byform.items() if any(col));newB=[[byform[f][o] for f in newforms] for o in range(18)];size=len(newforms);assert size<89
    for o in range(18):
     expected={};actual={}
     for c,f in zip(B[o],cubes):expected=plus(expected,f,c)
     for c,f in zip(newB[o],map(cubic,newforms)):actual=plus(actual,f,c)
     assert actual==expected
    tau=(0,)*17+(1,);answer=dict(status='exact-primary-candidate',source_sha256=sha(Path(__file__)),input_sha256=sha(source),pool_sha256=sha(poolpath),size=size,forms=newforms,B=[[str(c) for c in row] for row in newB],padding_index=newforms.index(tau),p=p,q=list(map(str,q)),removed=removed,scope='Exact smaller cube identity from a tangent cubic; no smaller Keller witness until full matrix and collision construction and independent replay.')
    (out/'candidate.json').write_text(json.dumps(answer,indent=2,sort_keys=True)+'\n');stats['best']=size;print('BEST',size,removed,flush=True);break
   if answer:break
  if answer:break
 receipt=dict(status='candidate' if answer else ('completed-directions-no-improvement' if complete else 'bounded-no-improvement'),source_sha256=sha(Path(__file__)),input_sha256=sha(source),pool_sha256=sha(poolpath),**stats,tangent_dimensions=dict(dimensions),seconds=round(time.monotonic()-started,3),scope='For930 selected rational p directions, all q with p^2*q in S are solved linearly and every four-column cancellation subset of each proportional old output direction is checked. New q and the three cube forms are not restricted to the old pool. This is not a search over arbitrary p or all rank-three cubics.')
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
