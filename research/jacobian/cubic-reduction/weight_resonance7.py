#!/usr/bin/env python3
"""Complete finite weight-ray frontend; no coefficient solve or bounded weight scan."""
from pathlib import Path
from itertools import combinations,permutations,product
from math import gcd
import json,hashlib,time
HERE=Path(__file__).resolve().parent;OUT=HERE/'weight-resonance7';OUT.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def primitive(v):
 g=gcd(*v);v=tuple(z//g for z in v)
 return tuple(-z for z in v) if next(z for z in v if z)<0 else v

def canonical(v):return min(primitive(tuple(v[j] for j in p)) for p in permutations(range(3)))
def cross(v,w):return (v[1]*w[2]-v[2]*w[1],v[2]*w[0]-v[0]*w[2],v[0]*w[1]-v[1]*w[0])
def cyclic(edges):
 return any((i,j) in edges and (j,i) in edges for i in range(3) for j in range(i+1,3)) or any(all((p[j],p[(j+1)%3]) in edges for j in range(3)) for p in permutations(range(3)))
def matching(support):
 rows=[{e for i,e in support if i==j and sum(e)>3} for j in range(3)]
 return max((len(assign) for k in range(4) for rr in combinations(range(3),k) for assign in product(*(rows[i] for i in rr)) if len(set(assign))==len(assign)),default=0)
def mixed(w):
 if 0 in w or not min(w)<0<max(w):return None
 if sum(z>0 for z in w)==2:w=tuple(-z for z in w)
 return (next(z for z in w if z>0),)+tuple(sorted([-z for z in w if z<0]))
def prior(m):
 if m is None:return None
 p,q,r=m
 if q==r and (p,q)!=(1,1):return 'KCUB32'
 if p>=2 and p in [q,r]:return 'KCUB33'
 if p==q==1 and r==2:return 'KCUB15:raw13-minimum-only'
 if p==q==1 and r==3:return 'KCUB18'
 if p==q==1 and r>=4:return 'KCUB34'
 if (p,q,r)==(1,2,3):return 'KCUB16'
 if (p,q,r)==(2,1,3):return 'KCUB17'
 if p==2 and q==1 and r>=4:return 'KCUB35'
 return None

def main():
 start=time.monotonic();atoms=[]
 for i in range(3):
  for e in product(range(8),repeat=3):
   if 2<=sum(e)<=7:atoms.append((i,e,tuple(e[j]-int(j==i) for j in range(3))))
 relations=sorted({v for i,e,v in atoms});planes=sorted({primitive(v) for v in relations});rays={}
 for v,w in combinations(relations,2):
  z=cross(v,w)
  if any(z):rays.setdefault(canonical(z),(v,w))
 records=[];cp={'status':'running','source_sha256':sha(Path(__file__)),'degree_cap':7,'labelled_atoms':len(atoms),'distinct_relation_vectors':len(relations),'primitive_relation_planes':len(planes),'candidate_rays':len(rays),'processed':0};(OUT/'checkpoint.json').write_text(json.dumps(cp,indent=2)+'\n')
 for w,(v0,v1) in sorted(rays.items()):
  support=[(i,e) for i,e,v in atoms if sum(a*b for a,b in zip(v,w))==0];edges=sorted({(j,i) for i,e in support for j in range(3) if j!=i and e[j]});maxima=[max([1]+[sum(e) for i,e in support if i==j]) for j in range(3)];m=mixed(w);gate='acyclic' if not cyclic(edges) else 'GLDC10:high-rank<3' if matching(support)<3 else prior(m) or 'unresolved'
  rec={'weight':list(w),'mixed_positive_negative_magnitudes':list(m) if m else None,'nonlinear_slots':len(support),'per_output_slots':[sum(i==j for i,e in support) for j in range(3)],'max_output_degrees':maxima,'high_degree_structural_rank':matching(support),'cross_edges':[list(e) for e in edges],'gate':gate,'support':[{'output':i,'exponent':list(e)} for i,e in support]};records.append(rec)
  if len(records)%100==0:cp['processed']=len(records);(OUT/'checkpoint.json').write_text(json.dumps(cp,indent=2)+'\n')
 rank1=[]
 for v in planes:
  supp=[(i,e) for i,e,r in atoms if not any(cross(v,r))];rank1.append({'relation_line':list(v),'support':[{'output':i,'exponent':list(e)} for i,e in supp],'boundary':'elementary-monomial-shear' if min(v)<0 else 'identity-by-positive-leading-scalar'})
 gates={}
 for r in records:gates[r['gate']]=gates.get(r['gate'],0)+1
 out={'status':'complete-primary-census','source_sha256':sha(Path(__file__)),'degree_cap':7,'labelled_atoms':len(atoms),'distinct_relation_vectors':len(relations),'primitive_relation_planes':len(planes),'rank2_weight_rays_up_to_permutation_and_sign':len(records),'gate_counts':gates,'rank1_planes':rank1,'records':records,'coverage_boundary':'Rank-zero affine; universal written rank-one automorphism lemma; rank-two kernel rays from complete native monomial relations. Gates are theorem imports/structural necessary bounds, not a new global exclusion.'};op=OUT/'census.json';op.write_text(json.dumps(out,indent=2)+'\n');cp.update(status='complete',processed=len(records),census_sha256=sha(op),elapsed_seconds=time.monotonic()-start);(OUT/'checkpoint.json').write_text(json.dumps(cp,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['rank1_planes','records']},indent=2));remaining=[r for r in records if r['gate']=='unresolved'];print('REMAINING',len(remaining));print([(r['mixed_positive_negative_magnitudes'] or r['weight'],r['nonlinear_slots'],r['max_output_degrees']) for r in sorted(remaining,key=lambda x:(x['nonlinear_slots'],x['weight']))])
if __name__=='__main__':main()
