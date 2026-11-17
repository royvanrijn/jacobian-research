#!/usr/bin/env python3
"""Three-new-form search with rank-one residuals and freely re-solved output columns."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
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
 out=HERE/'power-three-rank1';out.mkdir(exist_ok=True);source=HERE/'power-basis-anneal/candidate.json';d=json.loads(source.read_text());poolpath=HERE/'power-two-exchange/pool.json';pool=list(map(tuple,json.loads(poolpath.read_text())['forms']));B=[[Q(c) for c in row] for row in d['B']];n=89;old=set(map(tuple,d['forms']))
 cubes=[cubic(f) for f in d['forms']];mons=sorted(set().union(*(set(f) for f in cubes)));mi={e:i for i,e in enumerate(mons)}
 R,_=fmpq_mat([[str(f.get(e,0)) for f in cubes]+[str(int(i==j)) for j in range(len(mons))] for i,e in enumerate(mons)]).rref();RR=fmpq_mat([[str(R[i,n+j]) for j in range(len(mons))] for i in range(len(mons))])
 groups=defaultdict(list);inside=[]
 for index,f in enumerate(pool):
  p=cubic(f);v=RR*fmpq_mat([[str(p.get(e,0))] for e in mons]);a={i:Q(str(v[i,0])) for i in range(n) if v[i,0]};res={('i',i):Q(str(v[i,0])) for i in range(n,len(mons)) if v[i,0]};res.update({('o',e):Q(c) for e,c in p.items() if e not in mi})
  if not res:
   if f not in old:inside.append((index,a))
   continue
  scalar=res[min(res)];key=tuple((e,c/scalar) for e,c in sorted(res.items()));groups[key].append((index,scalar,a))
 weights=[1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093,1097,1103,1109]
 h=[sum((Q(w)*row[j] for w,row in zip(weights,B)),Q(0)) for j in range(n)]
 stats=dict(tasks=0,too_small_support=0,full_fits=0,max_canceled=0,best=89);started=time.monotonic();answer=None;complete=True
 def fit(u,v):
  T=sorted((set(u)|set(v))-{d['padding_index']})
  if len(T)<4:stats['too_small_support']+=1;return None
  rows={i:(u.get(i,Q(0)),v.get(i,Q(0))) for i in T};buckets={};parallel=defaultdict(list)
  for i,(a,b) in rows.items():
   c=a if a else b;parallel[(a/c,b/c,h[i]/c)].append(i)
  choices=[]
  for ids in parallel.values():
   if len(ids)>=4:
    j=ids[0];a,b=rows[j]
    choices.append(([row[j]/a for row in B],[Q(0)]*18) if a else ([Q(0)]*18,[row[j]/b for row in B]))
  for p,q in combinations(T,2):
   ap,bp=rows[p];aq,bq=rows[q];det=ap*bq-aq*bp
   if not det:continue
   x=(h[p]*bq-h[q]*bp)/det;y=(h[q]*ap-h[p]*aq)/det;key=(x,y)
   if key not in buckets:buckets[key]=[set(),[]]
   bucket=buckets[key];bucket[0].update([p,q]);bucket[1].append((p,q,det))
  for ids,pairs in buckets.values():
   if len(ids)<4:continue
   for p,q,det in pairs:
    ap,bp=rows[p];aq,bq=rows[q]
    choices.append(([(row[p]*bq-row[q]*bp)/det for row in B],[(row[q]*ap-row[p]*aq)/det for row in B]))
  for X,Y in choices:
   stats['full_fits']+=1;removed=[j for j in T if all(row[j]==rows[j][0]*x+rows[j][1]*y for row,x,y in zip(B,X,Y))]
   stats['max_canceled']=max(stats['max_canceled'],len(removed))
   if len(removed)>=4:return X,Y,removed
  return None
 def attempt(u,v,indices,W):
  nonlocal answer
  stats['tasks']+=1
  if stats['tasks']%100==0 and time.monotonic()-started>60:return False
  result=fit(u,v)
  if result:
   X,Y,removed=result;keep=[j for j in range(n) if j not in removed]
   newcols=[[x*W[0][k]+y*W[1][k] for x,y in zip(X,Y)] for k in range(3)]
   active=[k for k,col in enumerate(newcols) if any(col)];newforms=[d['forms'][j] for j in keep]+[pool[indices[k]] for k in active]
   newB=[[row[j]-X[o]*u.get(j,Q(0))-Y[o]*v.get(j,Q(0)) for j in keep]+[newcols[k][o] for k in active] for o,row in enumerate(B)]
   for o in range(18):
    expected={};actual={}
    for c,p in zip(B[o],cubes):expected=plus(expected,p,c)
    for c,p in zip(newB[o],map(cubic,newforms)):actual=plus(actual,p,c)
    assert expected==actual
   size=len(newforms);assert size<89
   answer=dict(status='exact-primary-candidate',source_sha256=sha(Path(__file__)),input_sha256=sha(source),pool_sha256=sha(poolpath),size=size,forms=newforms,B=[[str(c) for c in row] for row in newB],padding_index=keep.index(d['padding_index']),removed=removed,new_pool_indices=indices,scope='Exact smaller cube identity; requires full integer map and collision construction and independent replay.')
   (out/'candidate.json').write_text(json.dumps(answer,indent=2,sort_keys=True)+'\n');stats['best']=size;print('BEST',size,removed,flush=True);return False
  return True
 for group in groups.values():
  for (p,sp,ap),(q,sq,aq),(r,sr,ar) in combinations(group,3):
   if not attempt(plus(ap,ar,-sp/sr),plus(aq,ar,-sq/sr),[p,q,r],[[Q(1),Q(0),-sp/sr],[Q(0),Q(1),-sq/sr]]):complete=False;break
  if not complete:break
  for (p,sp,ap),(q,sq,aq) in combinations(group,2):
   u=plus(ap,aq,-sp/sq)
   for r,ar in inside:
    if not attempt(u,ar,[p,q,r],[[Q(1),-sp/sq,Q(0)],[Q(0),Q(0),Q(1)]]):complete=False;break
   if not complete:break
  if not complete:break
 receipt=dict(status='candidate' if answer else ('finite-pool-no-improvement' if complete else 'bounded-no-improvement'),source_sha256=sha(Path(__file__)),input_sha256=sha(source),pool_sha256=sha(poolpath),**stats,
   outside_residual_groups=len(groups),new_inside_forms=len(inside),seconds=round(time.monotonic()-started,3),
   scope='Padding-retaining three-new-form exchanges in the930-form pool with residual rank1. Includes three outside forms or two outside plus a genuinely new inside form. Scalar projection only filters necessary fits; full18-output equalities validate every surviving fit. Time bound60 seconds; no unrestricted three-form exclusion.')
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
