#!/usr/bin/env python3
"""Independent sequential-substitution replay and full determinant checks."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib
from verify import decode,variable,compose,add,scale,evaluate,derivative,identity,mmul,madd,determinant3,mul
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 path=HERE/'cubic12-integer-witness.json';raw=json.loads(path.read_text());F=[decode(p,12) for p in raw['coefficients']];X=[variable(j,12) for j in range(12)]
 folder=HERE/'carrier-subset-frontier12';primary=json.loads((folder/'checkpoint.json').read_text())
 assert primary['status']=='complete' and primary['source_sha256']==sha(HERE/'carrier_subset_frontier12.py') and primary['input_sha256']==sha(path)
 assert len(primary['cases'])==512
 rows=[];best={};witness_hashes={};detrecords=[]
 for mask in range(512):
  deleted=[j for j in range(3,12) if mask&(1<<(j-3))];keep=[j for j in range(12) if j not in deleted];current=list(F)
  # Unlike the producer's full embedding composition, eliminate the live map one row at a time.
  for j in deleted:
   graph=add(X[j],current[j],Q(-1));assert all(e[j]==0 for e in graph)
   substitution=list(X);substitution[j]=graph
   current=[compose(p,substitution,12) for p in current];assert not current[j]
  assert all(not current[j] for j in deleted)
  assert all(not any(e[j] for j in deleted) for i in keep for e in current[i])
  G=[{tuple(e[j] for j in keep):c for e,c in current[i].items()} for i in keep];n=len(keep)
  degree=max(sum(e) for p in G for e in p);count=sum(sum(e)>=2 for p in G for e in p)
  record={'mask':mask,'dimension':n,'degree':degree,'nonlinear_total':count,'deleted':[j+1 for j in deleted]};assert record==primary['cases'][mask];rows.append(record)
  points=[[Q(v) for j,v in enumerate(point) if j in keep] for point in raw['points']];images=[[evaluate(p,point) for p in G] for point in points]
  assert points[0]!=points[1] and images[0]==images[1]
  if n not in best or (degree,count)<(best[n]['degree'],best[n]['nonlinear_total']):best[n]=record
  if mask not in [primary['best_by_dimension'][str(a)]['mask'] for a in (3,4,5,7,12)]:continue
  wp=folder/('best-n%02d.json'%n);d=json.loads(wp.read_text())
  assert all(d[k]==v for k,v in record.items())
  assert d['source_sha256']==primary['source_sha256'] and d['input_sha256']==sha(path)
  assert [decode(p,n) for p in d['coefficients']]==G
  if n==4:
   a,b,c,u=[variable(j,4) for j in range(4)];aa=mul(a,a);bb=mul(b,b);dd=add(scale(b,3),mul(a,c))
   displayed=[add(a,mul(aa,dd),Q(-1)),add(add(add(b,scale(mul(a,c),3)),scale(mul(mul(a,b),add(scale(b,2),mul(a,c))),12)),scale(mul(mul(aa,bb),dd),12)),add(add(add(add(add(c,scale(bb,8)),scale(mul(a,u),8)),scale(mul(mul(a,b),c),6)),scale(mul(mul(a,bb),dd),12)),scale(mul(mul(aa,u),dd),-8)),add(u,mul(bb,b))]
   assert displayed==G and [sum(sum(e)>=2 for e in p) for p in G]==[2,5,7,1]
  assert [[Q(v) for v in point] for point in d['points']]==points
  assert list(map(Q,d['common_image']))==images[0]
  assert d['source_variables']==[raw['variables'][j] for j in keep]
  # Reconstruct the source embedding from the independently substituted live carrier equations.
  E=[decode(p,12) for p in d['source_embedding']];assert len(E)==12
  assert all(E[j]==X[j] for j in keep)
  assert all(not compose(F[j],E,12) for j in deleted)
  assert all(compose(F[j],E,12)==current[j] for j in keep)
  J=[[derivative(p,j) for j in range(n)] for p in G];r=n-3;one={(0,)*n:Q(1)}
  assert all(J[i][j].get((0,)*n,0)==Q(i==j) for i in range(n) for j in range(n))
  if r:
   D=[row[3:] for row in J[3:]];I=identity(r,n);N=madd(D,I,Q(-1))
   assert all(not N[i][j] for i in range(r) for j in range(i,r))
   inv=I;power=I
   for k in range(1,r):power=mmul(power,N);inv=madd(inv,power,Q((-1)**k))
   assert mmul(power,N)==[[{} for j in range(r)] for i in range(r)]
   assert mmul(D,inv)==I and mmul(inv,D)==I
   schur=madd([row[:3] for row in J[:3]],mmul(mmul([row[3:] for row in J[:3]],inv),[row[:3] for row in J[3:]]),Q(-1))
  else:schur=J
  assert determinant3(schur,n)==one
  assert all(c.denominator==1 for p in G for c in p.values())
  witness_hashes[str(n)]=sha(wp);detrecords.append({'dimension':n,'degree':degree,'nonlinear_total':count,'determinant_one':True,'full_symbolic_schur_complement':True,'distinct_integer_collision':True,'integral_coefficients':True,'deleted':record['deleted']})
  print('WITNESS PASS',n,degree,count,flush=True)
 assert {str(n):v for n,v in best.items()}==primary['best_by_dimension']
 triples=sorted(set((a['dimension'],a['degree'],a['nonlinear_total']) for a in rows));pareto=[p for p in triples if not any(q!=p and all(a<=b for a,b in zip(q,p)) for q in triples)]
 assert pareto==list(map(tuple,primary['pareto_triples']))==[(3,7,13),(4,6,15),(5,5,17),(7,4,23),(12,3,36)]
 assert min(a['nonlinear_total'] for a in rows)==13 and min(a['dimension'] for a in rows if a['degree']==3)==12
 receipt={'status':'pass','source_sha256':sha(Path(__file__)),'arithmetic_sha256':sha(HERE/'verify.py'),'proof_sha256':sha(HERE/'CARRIER_SUBSET_FRONTIER.md'),'input_sha256':sha(path),'producer_sha256':sha(HERE/'carrier_subset_frontier12.py'),'checkpoint_sha256':sha(folder/'checkpoint.json'),'witness_hashes':witness_hashes,'all_subsets_replayed':512,'all_collisions_checked':512,'pareto_triples':pareto,'determinants':sorted(detrecords,key=lambda x:x['dimension']),'scope':'Exhaustion only of512 literal zero-carrier subsets of KCUB3, no target repair. Five full polynomial determinant-one collision witnesses; no global dimension, degree or support minimum.'}
 (folder/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS all512 restrictions and five exact determinants; scoped frontier',pareto,flush=True)
if __name__=='__main__':main()
