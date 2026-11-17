#!/usr/bin/env python3
"""Independent Fraction replay of the arbitrary single-cube compression obstruction."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,combinations_with_replacement
from collections import defaultdict
import hashlib,json,time
from verify import add,mul,variable
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def scale(p,c):return {e:v*c for e,v in p.items() if v*c}
def decode(terms,n):
 out={}
 for e,c in terms:
  assert len(e)==n and all(type(a) is int and a>=0 for a in e)
  e=tuple(e);assert e not in out;v=Q(c)
  assert v;out[e]=v
 return out
def lin(f):return {tuple(int(i==j) for i in range(len(f))):Q(c) for j,c in enumerate(f) if c}
def echelon(polys):
 basis={}
 for p0 in polys:
  p=reduce(p0,basis)
  if p:
   e=min(p);basis[e]=scale(p,1/p[e])
 return basis
def reduce(p,basis):
 p=dict(p)
 for e,b in sorted(basis.items()):
  c=p.get(e,Q(0))
  if c:p=add(p,b,-c)
 return p
def main():
 started=time.monotonic();out=HERE/"power-single-compress";source=HERE/"power-basis-anneal/candidate.json";d=json.loads(source.read_text());n=d["size"];assert n==89
 construction=json.loads((HERE/"power89/validation.json").read_text())
 assert construction["status"]=="pass" and construction["candidate_sha256"]==sha(source)
 assert construction["source_sha256"]==sha(HERE/"verify_power89.py")
 B=[[Q(c) for c in row] for row in d["B"]];forms=d["forms"]
 assert len(B)==18 and all(len(row)==n for row in B)
 lins=[lin(f) for f in forms];cubes=[mul(mul(f,f),f) for f in lins]
 basis={}
 for j,cube in enumerate(cubes):
  p=dict(cube);coord={j:Q(1)}
  for e,(b,bc) in sorted(basis.items()):
   c=p.get(e,Q(0))
   if c:p=add(p,b,-c);coord=add(coord,bc,-c)
  assert p
  e=min(p);c=p[e];basis[e]=(scale(p,1/c),scale(coord,1/c))
 def project(p):
  p=dict(p);coord={}
  for e,(b,bc) in sorted(basis.items()):
   c=p.get(e,Q(0))
   if c:p=add(p,b,-c);coord=add(coord,bc,c)
  return p,coord
 for j,cube in enumerate(cubes):assert project(cube)==({}, {j:Q(1)})
 groups=defaultdict(list)
 for j in range(n):
  v=tuple(row[j] for row in B);c=next(c for c in v if c);groups[tuple(a/c for a in v)].append(j)
 pairs=sorted((j,k) for ids in groups.values() for j,k in combinations(ids,2))
 primary=json.loads((out/"checkpoint.json").read_text());certification=json.loads((out/"certification.json").read_text())
 assert primary["status"]=="primary-all-closed" and primary["source_sha256"]==sha(HERE/"power_single_compress.py")
 assert certification["status"]=="exact-primary-identities" and certification["source_sha256"]==sha(HERE/"certify_power_single_compress.py")
 assert certification["primary_checkpoint_sha256"]==sha(out/"checkpoint.json")
 cases=[json.loads((out/("case-%02d.json"%i)).read_text()) for i in range(15)]
 masks=[];adj=[0]*18
 for f in forms:
  for i,j in combinations([i for i,c in enumerate(f) if c],2):adj[i]|=1<<j;adj[j]|=1<<i
 for case in cases:
  S=case["support"];assert S==sorted(set(S));assert all(0<=i<18 for i in S)
  assert all(adj[i]&(1<<j) for i,j in combinations(S,2));masks.append(sum(1<<i for i in S))
 assert len(set(masks))==15
 # Exhaustive subset check is independent of the producer's maximal-clique algorithm.
 clique=bytearray(1<<18);clique[0]=1;covered=0
 for mask in range(1,1<<18):
  bit=mask&-mask;i=bit.bit_length()-1;rest=mask^bit
  if clique[rest] and rest&~adj[i]==0:
   clique[mask]=1;covered+=1;assert any(mask&~c==0 for c in masks),mask
 identities=0;max_degree=0;case_receipts=[]
 for case in cases:
  i=case["case"];S=case["support"];m=len(S);constraints=[decode(f,m) for f in case["constraints"]];coords=[decode(f,m) for f in case["coordinates"]];assert len(coords)==n
  assert all(sum(e)==3 for f in constraints+coords for e in f)
  # Reconstruct membership equations with an independently built sparse elimination basis.
  actual_coords=[{} for _ in range(n)];residual_polys={}
  for ids in combinations_with_replacement(range(m),3):
   e=[0]*18;u=[0]*m
   for j in ids:e[S[j]]+=1;u[j]+=1
   e=tuple(e);u=tuple(u);mult=1 if ids[0]==ids[2] else (3 if ids[0]==ids[1] or ids[1]==ids[2] else 6)
   residual,coordinate=project({e:Q(mult)})
   for a,c in coordinate.items():actual_coords[a][u]=c
   for exponent,c in residual.items():residual_polys.setdefault(exponent,{})[u]=c
  I=echelon(constraints);J=echelon(residual_polys.values())
  assert len(I)==len(J)
  assert all(not reduce(f,I) for f in residual_polys.values())
  assert all(not reduce(f,J) for f in constraints)
  assert all(not reduce(add(a,b,Q(-1)),I) for a,b in zip(coords,actual_coords))
  # Coordinates zero as polynomials can never cancel the corresponding nonzero old coefficient column.
  selected=[]
  for j,k in pairs:
   if not coords[j] or not coords[k]:continue
   v=next(o for o in range(18) if B[o][j]);selected.append((j,k,B[v][k]/B[v][j]))
  assert set((j,k,str(r)) for j,k,r in selected)==set(tuple(p) for p in case["pairs"])
  assert len(selected)==len(case["pairs"])
  selected=[(j,k,Q(r)) for j,k,r in case["pairs"]]
  cp=out/("cert-%02d.json"%i);certificate=json.loads(cp.read_text())
  assert certificate["case_sha256"]==sha(out/("case-%02d.json"%i))
  assert certificate["source_sha256"]==sha(HERE/"certify_power_single_compress.py")
  assert certificate["singular_source_sha256"]==sha(out/("cert-%02d.sing"%i))
  assert certification["cases"][i]["certificate_sha256"]==sha(cp)
  assert len(certificate["certificates"])==len(selected)
  one={(0,)*m:Q(1)}
  for (j,k,ratio),cert in zip(selected,certificate["certificates"]):
   assert cert["pair"]==[j,k] and Q(cert["ratio"])==ratio
   generators=constraints+[add(coords[j],one,Q(-1)),add(coords[k],one,-ratio)]
   multipliers=[decode(f,m) for f in cert["multipliers"]];assert len(multipliers)==len(generators)
   result={}
   for a,b in zip(generators,multipliers):result=add(result,mul(a,b))
   assert result==one,(i,j,k)
   max_degree=max(max_degree,max((sum(e) for p in multipliers for e in p),default=0));identities+=1
  case_receipts.append(dict(case=i,support=S,membership_rank=len(I),identities=len(selected),certificate_sha256=sha(cp)))
 assert identities==certification["identities"]==86
 # Positive control: four independent binary cubes really can collapse to one new cube.
 a=variable(0,2);b=variable(1,2);cube=lambda f:mul(mul(f,f),f)
 control=add(add(scale(cube(a),Q(-3)),scale(cube(b),Q(6))),add(scale(cube(add(a,b)),Q(3)),cube(add(a,b,Q(-1)))))
 assert control==cube(add(a,b,Q(2)))
 receipt=dict(status="pass",source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/"verify.py"),
  input_sha256=sha(source),construction_validation_sha256=sha(HERE/"power89/validation.json"),
  primary_checkpoint_sha256=sha(out/"checkpoint.json"),certification_sha256=sha(out/"certification.json"),
  independent_cube_rank=89,coordinate_support_cases=15,covered_nonempty_cliques=covered,
  proportional_coefficient_pairs=len(pairs),unit_identities=identities,max_multiplier_degree=max_degree,
  all_membership_constraints_reconstructed=True,all_coordinate_maps_verified_mod_membership=True,
  all_coordinate_supports_covered=True,all_cancellation_pairs_covered=True,positive_compression_control=True,
  cases=case_receipts,seconds=round(time.monotonic()-started,3),
  scope="The saved89-form padded companion identity cannot be shortened by removing at least two selected forms, adding one arbitrary new linear cube, and re-solving all retained coefficients. All characteristic-zero extensions are covered. Other decompositions using two or more new forms and arbitrary smaller Keller constructions remain open.")
 (out/"validation.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
 print(json.dumps({k:v for k,v in receipt.items() if k!="cases"}))
if __name__=="__main__":main()
