#!/usr/bin/env python3
"""Exact three-new-form search using rank-two residual circuits in the930-form pool."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,product
from collections import defaultdict
import hashlib,json,time
from flint import fmpq_mat
from power_line_compress import cubic
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def plus(a,b,c):
 out=dict(a)
 for i,v in b.items():out[i]=out.get(i,Q(0))+c*v
 return {i:v for i,v in out.items() if v}
def main():
 out=HERE/"power-three-rank2";out.mkdir(exist_ok=True);source=HERE/"power-basis-anneal/candidate.json";d=json.loads(source.read_text());poolpath=HERE/"power-two-exchange/pool.json";pool=list(map(tuple,json.loads(poolpath.read_text())['forms']));B=[[Q(c) for c in row] for row in d['B']];n=89
 cubes=[cubic(f) for f in d['forms']];mons=sorted(set().union(*(set(f) for f in cubes)));mi={e:i for i,e in enumerate(mons)}
 R,_=fmpq_mat([[str(f.get(e,0)) for f in cubes]+[str(int(i==j)) for j in range(len(mons))] for i,e in enumerate(mons)]).rref();RR=fmpq_mat([[str(R[i,n+j]) for j in range(len(mons))] for i in range(len(mons))])
 grouping={};directions=[];groups=[]
 for index,f in enumerate(pool):
  p=cubic(f);v=RR*fmpq_mat([[str(p.get(e,0))] for e in mons]);res={('i',i):Q(str(v[i,0])) for i in range(n,len(mons)) if v[i,0]};res.update({('o',e):Q(c) for e,c in p.items() if e not in mi})
  if not res:continue
  scalar=res[min(res)];normal={e:c/scalar for e,c in res.items()};key=tuple(sorted(normal.items()));coords={i:Q(str(v[i,0])) for i in range(n) if v[i,0]}
  if key not in grouping:grouping[key]=len(directions);directions.append(normal);groups.append([])
  groups[grouping[key]].append((index,scalar,coords))
 assert len(directions)==445
 circuits=[]
 for i,r in enumerate(directions):
  pivot=min(r);quotients=defaultdict(list)
  for j in range(i+1,len(directions)):
   s=directions[j];a=s.get(pivot,Q(0));q=plus(s,r,-a);assert q;t=q[min(q)];key=tuple((e,c/t) for e,c in sorted(q.items()));quotients[key].append((j,a,t))
  for group in quotients.values():
   for (j,aj,tj),(k,ak,tk) in combinations(group,2):
    ci=tj*ak-tk*aj;assert ci;circuits.append((i,j,k,ci,tk,-tj))
 print('DIRECTIONS',len(directions),'CIRCUITS',len(circuits),flush=True)
 started=time.monotonic();tested=0;best=89;bestremoved=0;answer=None;complete=True
 for i,j,k,ci,cj,ck in circuits:
  for (p,sp,ap),(q,sq,aq),(r,sr,ar) in product(groups[i],groups[j],groups[k]):
   if time.monotonic()-started>60:complete=False;break
   tested+=1;weights=[ci/sp,cj/sq,ck/sr];a=plus(plus({h:weights[0]*c for h,c in ap.items()},aq,weights[1]),ar,weights[2])
   cancel=defaultdict(list)
   for h,c in a.items():
    if not B[17][h]:cancel[tuple(row[h]/c for row in B)].append(h)
   removed=max(cancel.values(),key=len,default=[]);bestremoved=max(bestremoved,len(removed))
   if len(removed)>=4 and n-len(removed)+3<best:
    best=n-len(removed)+3;ref=removed[0];b=[row[ref]/a[ref] for row in B];keep=[h for h in range(n) if h not in removed]
    newforms=[d['forms'][h] for h in keep]+[pool[p],pool[q],pool[r]]
    newB=[[row[h]-b[o]*a.get(h,Q(0)) for h in keep]+[b[o]*w for w in weights] for o,row in enumerate(B)]
    for o in range(18):
     expected={};actual={}
     for c,f in zip(B[o],cubes):expected=plus(expected,f,c)
     for c,f in zip(newB[o],map(cubic,newforms)):actual=plus(actual,f,c)
     assert expected==actual
    answer=dict(status='exact-primary-candidate',source_sha256=sha(Path(__file__)),input_sha256=sha(source),pool_sha256=sha(poolpath),size=best,forms=newforms,B=[[str(c) for c in row] for row in newB],padding_index=keep.index(d['padding_index']),removed=removed,new_pool_indices=[p,q,r],weights=list(map(str,weights)),scope='Exact smaller cube identity. A smaller Keller matrix and collision require separate full construction and verification.')
    (out/'candidate.json').write_text(json.dumps(answer,indent=2,sort_keys=True)+'\n');print('BEST',best,removed,flush=True)
  if not complete:break
 receipt=dict(status='candidate' if answer else ('finite-pool-no-improvement' if complete else 'bounded-no-improvement'),source_sha256=sha(Path(__file__)),input_sha256=sha(source),pool_sha256=sha(poolpath),residual_directions=len(directions),direction_circuits=len(circuits),triples_tested=tested,max_canceled=bestremoved,best=best,seconds=round(time.monotonic()-started,3),scope='Three outside candidate forms whose scalar residues have rank2, forcing proportional output columns. Does not cover residual rank1 or unrestricted forms. No optimality claim.')
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
