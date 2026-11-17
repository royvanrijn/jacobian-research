#!/usr/bin/env python3
"""Compare complete literal zero-carrier subset frontiers of three saved maps."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,sys
from verify import decode,variable,compose,add,scale
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(F):return [[{'exponent':list(e),'coefficient':str(c)} for e,c in sorted(p.items())] for p in F]
def main():
 sys.path.insert(0,str(ROOT/'archive/non-elliptic/scripts'))
 from audit_k12_coordinate_pair_frontier import build_k12
 import sympy as sp
 z,expressions=build_k12();old=[{tuple(e):Q(str(c)) for e,c in sp.Poly(p,*z).terms()} for p in expressions]
 maps={'BCR2':old}
 for label,name in [('KCUB3','cubic12-integer-witness.json'),('KCUB2','cubic13-sparse-witness.json')]:
  d=json.loads((HERE/name).read_text());maps[label]=[decode(p,d['dimension']) for p in d['coefficients']]
 folder=HERE/'carrier-subset-comparison';folder.mkdir(exist_ok=True);results={}
 for label,F in maps.items():
  n=len(F);X=[variable(i,n) for i in range(n)];nonlin=[add(F[j],X[j],Q(-1)) for j in range(n)]
  order=[];pending=set(range(3,n))
  while pending:
   ready=[j for j in sorted(pending) if not any(e[k] for e in nonlin[j] for k in pending)]
   assert ready,'nontriangular carrier graph';order.extend(ready);pending.difference_update(ready)
  rows=[];best={}
  for mask in range(1<<(n-3)):
   deleted=[j for j in order if mask&(1<<(j-3))];keep=[j for j in range(n) if j not in deleted];E=list(X)
   for j in deleted:E[j]=scale(compose(nonlin[j],E,n),Q(-1))
   full=[compose(p,E,n) for p in F];assert all(not full[j] for j in deleted)
   G=[{tuple(e[j] for j in keep):c for e,c in full[i].items()} for i in keep]
   assert all(not any(e[j] for j in deleted) for i in keep for e in full[i])
   degree=max(sum(e) for p in G for e in p);count=sum(sum(e)>=2 for p in G for e in p);dim=len(keep)
   rec={'mask':mask,'dimension':dim,'degree':degree,'nonlinear_total':count,'deleted':[j+1 for j in deleted]};rows.append(rec)
   if dim not in best or (degree,count)<(best[dim]['degree'],best[dim]['nonlinear_total']):best[dim]=rec
  triples=sorted(set((r['dimension'],r['degree'],r['nonlinear_total']) for r in rows));pareto=[p for p in triples if not any(q!=p and all(a<=b for a,b in zip(q,p)) for q in triples)]
  results[label]={'carrier_order':[j+1 for j in order],'cases':rows,'best_by_dimension':{str(k):v for k,v in sorted(best.items())},'pareto_triples':pareto}
  (folder/(label+'-input.json')).write_text(json.dumps({'coefficients':enc(F),'source_sha256':sha(Path(__file__))},indent=2)+'\n')
  print(label,'PARETO',pareto,flush=True)
 result={'status':'complete-primary-comparison','source_sha256':sha(Path(__file__)),'families':results,'scope':'Only literal carrier-zero restrictions of these three maps; no global or literature optimum. Independent replay required.'}
 (folder/'checkpoint.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
