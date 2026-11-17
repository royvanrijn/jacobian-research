#!/usr/bin/env python3
"""Independent exact quotient-space replay of the930-form two-new-cube test."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from collections import defaultdict
import hashlib,json,time
from verify import add,mul
from verify_power_single_compress import lin,scale
from verify_power_span_optimality import rows_rref
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 started=time.monotonic();out=HERE/"power-two-exchange";source=HERE/"power-basis-anneal/candidate.json";d=json.loads(source.read_text());B=[[Q(c) for c in row] for row in d["B"]];n=d["size"];assert n==89
 span=json.loads((HERE/"power-span-split/validation.json").read_text());assert span["status"]=="pass" and span["source_sha256"]==sha(HERE/"verify_power_span_optimality.py")
 assert span["input_sha256"]==sha(source)
 poolpath=out/"pool.json";pool=json.loads(poolpath.read_text())["forms"];assert len(pool)==930 and len({tuple(f) for f in pool})==930
 basis={};cube=lambda f:mul(mul(f,f),f)
 for j,f in enumerate(d["forms"]):
  p=cube(lin(f));coord={j:Q(1)}
  for e,(a,ac) in sorted(basis.items()):
   c=p.get(e,Q(0))
   if c:p=add(p,a,-c);coord=add(coord,ac,-c)
  assert p;e=min(p);c=p[e];basis[e]=(scale(p,1/c),scale(coord,1/c))
 def project(p):
  p=dict(p);coord={}
  for e,(a,ac) in sorted(basis.items()):
   c=p.get(e,Q(0))
   if c:p=add(p,a,-c);coord=add(coord,ac,c)
  return p,coord
 groups=defaultdict(list);inside=0;allcubes=[]
 for f in pool:
  assert len(f)==18 and all(type(c) is int for c in f)
  p=cube(lin(f));allcubes.append(p);residual,coords=project(p)
  if not residual:inside+=1;continue
  value=residual[min(residual)];signature=tuple((e,c/value) for e,c in sorted(residual.items()))
  groups[signature].append((value,coords))
 tested=0;hist=defaultdict(int)
 for group in groups.values():
  for (x,a),(y,b) in combinations(group,2):
   tested+=1;g=add(a,b,-x/y);canceled=defaultdict(list)
   for j,c in g.items():
    if c and not B[17][j]:canceled[tuple(row[j]/c for row in B)].append(j)
   maximum=max(map(len,canceled.values()),default=0);hist[maximum]+=1;assert maximum<=2
 # Independent rank of the whole finite cube pool, using sparse elimination.
 from verify_power_single_compress import echelon
 poolrank=len(echelon(allcubes));assert poolrank==154
 primary=json.loads((out/'checkpoint.json').read_text());assert primary['status']=='finite-pool-no-improvement' and primary['source_sha256']==sha(HERE/'power_two_exchange.py')
 assert primary['pool_sha256']==sha(poolpath) and primary['input_sha256']==sha(source)
 assert (primary['inside_saved_span'],primary['nonzero_residual_directions'],primary['pairs_tested'])==(inside,len(groups),tested)
 assert {int(k):v for k,v in primary['cancellation_histogram'].items()}==dict(hist)
 # A further cheap raw triple test: every proportional-column triple has source rank3,
 # except the five exact tau^2 Wi identities. None can collapse to two cubes without
 # changing retained coefficient columns.
 outputgroups=defaultdict(list)
 for j in range(n):
  v=tuple(row[j] for row in B);c=next(a for a in v if a);outputgroups[tuple(a/c for a in v)].append((j,c))
 rank3=0;binary=0
 for group in outputgroups.values():
  for triple in combinations(group,3):
   ids=[j for j,c in triple];rank=len(rows_rref([d['forms'][j] for j in ids],18))
   if rank==3:rank3+=1;continue
   assert rank==2;binary+=1;g={}
   for j,c in triple:g=add(g,cube(lin(d['forms'][j])),c)
   assert len(g)==1;e,c=next(iter(g.items()));assert e[17]==2 and sum(e)==3 and any(e[i]==1 for i in range(12,17))
 assert (rank3,binary)==(209,5)
 receipt=dict(status='pass',source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/'verify.py'),helper_sha256=sha(HERE/'verify_power_single_compress.py'),
  primary_sha256=sha(out/'checkpoint.json'),input_sha256=sha(source),pool_sha256=sha(poolpath),fixed_span_validation_sha256=sha(HERE/'power-span-split/validation.json'),
  candidate_forms=len(pool),cube_pool_rank=poolrank,inside_saved_span=inside,nonzero_residual_directions=len(groups),pairs_tested=tested,
  max_canceled=max(hist),cancellation_histogram=dict(hist),raw_rank3_triples=rank3,raw_tangent_binary_triples=binary,seconds=round(time.monotonic()-started,3),
  scope='No shortening by adding at most two forms from this930-form pool while re-solving all retained coefficients. Pair quotients checked exactly; both-inside-span cases use the separate fixed-span theorem. This is not an unrestricted two-new-form obstruction or global rank minimum.')
 (out/'validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
