#!/usr/bin/env python3
"""Exact two-new-form exchange in the union of previous candidate pools."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product,combinations
from collections import defaultdict
from math import gcd
from functools import reduce
import hashlib,json,time
from flint import fmpq_mat
from power_line_compress import cubic
from power_search import terms,decomposition,DB
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def primitive(v):
 g=reduce(gcd,v,0)
 if not g:return None
 if next(x for x in v if x)<0:g=-g
 return tuple(x//g for x in v)
def main():
 started=time.monotonic();out=HERE/"power-two-exchange";out.mkdir(exist_ok=True);source=HERE/"power-basis-anneal/candidate.json";d=json.loads(source.read_text());n=d["size"];assert n==89
 B=[[Q(c) for c in row] for row in d["B"]];forms=list(map(tuple,d["forms"]));pool=set(forms)
 for _,_,factors in terms(True):
  for mode in range(5):pool.update(decomposition(factors,mode))
  for weights in product(range(-2,3),repeat=3):
   f=primitive(tuple(sum(w*v[j] for w,v in zip(weights,factors)) for j in range(18)))
   if f:pool.add(f)
 pencil=HERE/"power-basis-pencils/pool.json";pool.update(map(tuple,json.loads(pencil.read_text())["forms"]));pool=sorted(pool)
 (out/"pool.json").write_text(json.dumps(dict(forms=pool),indent=2)+'\n')
 cubes=[cubic(f) for f in forms];mons=sorted(set().union(*(set(p) for p in cubes)));mi={e:i for i,e in enumerate(mons)}
 joint=fmpq_mat([[str(p.get(e,0)) for p in cubes]+[str(int(i==j)) for j in range(len(mons))] for i,e in enumerate(mons)])
 rr,rank=joint.rref();assert rank==len(mons) and all(rr[i,i]==1 for i in range(n))
 R=fmpq_mat([[str(rr[i,n+j]) for j in range(len(mons))] for i in range(len(mons))])
 groups=defaultdict(list);inside=0
 for index,f in enumerate(pool):
  p=cubic(f);v=R*fmpq_mat([[str(p.get(e,0))] for e in mons])
  residual={('i',i):Q(str(v[i,0])) for i in range(n,len(mons)) if v[i,0]}
  residual.update({('o',e):Q(c) for e,c in p.items() if e not in mi})
  if not residual:inside+=1;continue
  scalar=residual[min(residual)];key=tuple((e,c/scalar) for e,c in sorted(residual.items()));coordinate=[Q(str(v[i,0])) for i in range(n)]
  groups[key].append((index,scalar,coordinate))
 tested=0;best=89;answer=None;histogram=defaultdict(int)
 for group in groups.values():
  for (p,cp,ap),(q,cq,aq) in combinations(group,2):
   tested+=1;ratio=-cp/cq;a=[x+ratio*y for x,y in zip(ap,aq)];cancel=defaultdict(list)
   for j,c in enumerate(a):
    if c and not B[17][j]:cancel[tuple(row[j]/c for row in B)].append(j)
   removed=max(cancel.values(),key=len,default=[]);histogram[len(removed)]+=1
   if len(removed)<3:continue
   size=n-len(removed)+2
   if size>=best:continue
   best=size;j=removed[0];b=[row[j]/a[j] for row in B];keep=[h for h in range(n) if h not in removed]
   newB=[[row[h]-b[o]*a[h] for h in keep]+[b[o],b[o]*ratio] for o,row in enumerate(B)]
   newforms=[forms[h] for h in keep]+[pool[p],pool[q]];newcubes=[cubic(f) for f in newforms]
   for o in range(18):
    expected={};actual={}
    for co,cube in zip(B[o],cubes):
     for e,c in cube.items():expected[e]=expected.get(e,Q(0))+co*c
    for co,cube in zip(newB[o],newcubes):
     for e,c in cube.items():actual[e]=actual.get(e,Q(0))+co*c
    assert {e:c for e,c in expected.items() if c}=={e:c for e,c in actual.items() if c}
   answer=dict(status='exact-identity-candidate',source_sha256=sha(Path(__file__)),input_sha256=sha(source),
     size=size,forms=newforms,B=[[str(c) for c in row] for row in newB],padding_index=keep.index(d['padding_index']),
     removed=removed,new_pool_indices=[p,q],ratio=str(ratio),
     scope='Exact primary smaller cube identity; no smaller Keller map until full integer factorization, determinant and collision replay.')
   (out/'candidate.json').write_text(json.dumps(answer,indent=2,sort_keys=True)+'\n');print('BEST',best,removed,flush=True)
 receipt=dict(status='candidate' if answer else 'finite-pool-no-improvement',source_sha256=sha(Path(__file__)),input_sha256=sha(source),
   pool_sha256=sha(out/'pool.json'),pencil_pool_sha256=sha(pencil),candidate_forms=len(pool),inside_saved_span=inside,
   nonzero_residual_directions=len(groups),pairs_tested=tested,max_canceled=max(histogram,default=0),
   cancellation_histogram=dict(histogram),best=best,seconds=round(time.monotonic()-started,3),
   source_inputs={name:sha(DB/name) for name in ['verify.py','polarize.py','polarization.json']},
   scope='Exhaustive two-new-candidate exchange for this finite pool, allowing all retained coefficients to change. Both-inside-span pairs are excluded by the separately checked fixed-span theorem. Not a global two-form exclusion.')
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
