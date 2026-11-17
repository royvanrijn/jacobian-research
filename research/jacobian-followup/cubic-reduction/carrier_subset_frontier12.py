#!/usr/bin/env python3
"""Exact complete512 literal carrier-subset restrictions of the saved cubic."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
import json,hashlib,time
from verify import decode,variable,compose,add,scale,evaluate
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(F):return [[{'exponent':list(e),'coefficient':str(c)} for e,c in sorted(p.items())] for p in F]
def main():
 folder=HERE/'carrier-subset-frontier12';folder.mkdir(exist_ok=True);path=HERE/'cubic12-integer-witness.json';raw=json.loads(path.read_text());F=[decode(p,12) for p in raw['coefficients']];X=[variable(i,12) for i in range(12)];points=[[Q(v) for v in p] for p in raw['points']];rows=[];best={};start=time.monotonic()
 for mask in range(512):
  deleted=[j for j in range(3,12) if mask&(1<<(j-3))];keep=[i for i in range(12) if i not in deleted];E=list(X)
  for j in deleted:E[j]=scale(compose(add(F[j],X[j],Q(-1)),E,12),Q(-1))
  full=[compose(p,E,12) for p in F];assert all(not full[j] for j in deleted)
  assert all(not any(e[j] for j in deleted) for i in keep for e in full[i])
  G=[{tuple(e[j] for j in keep):c for e,c in full[i].items()} for i in keep]
  degree=max(sum(e) for p in G for e in p);count=sum(sum(e)>=2 for p in G for e in p)
  pp=[[pt[j] for j in keep] for pt in points];im=[[evaluate(p,pt) for p in G] for pt in pp];assert pp[0]!=pp[1] and im[0]==im[1]
  n=len(keep);rec={'mask':mask,'dimension':n,'degree':degree,'nonlinear_total':count,'deleted':[j+1 for j in deleted]};rows.append(rec)
  if n not in best or (degree,count)<(best[n]['degree'],best[n]['nonlinear_total']):
   best[n]=rec;payload={'status':'graph-inherited-witness-candidate','source_sha256':sha(Path(__file__)),'input_sha256':sha(path),**rec,'source_variables':[raw['variables'][j] for j in keep],'source_embedding':enc(E),'coefficients':enc(G),'points':[[str(v) for v in p] for p in pp],'common_image':[str(v) for v in im[0]],'scope':'Exact primary carrier-zero restriction. Keller inheritance by triangular graph/source block proof; selected determinants and all coefficients require independent replay before a new bound is accepted.'};(folder/('best-n%02d.json'%n)).write_text(json.dumps(payload,indent=2)+'\n')
  if (mask+1)%64==0:print('SUBSETS',mask+1,'BEST',[(n,best[n]['degree'],best[n]['nonlinear_total']) for n in sorted(best)],flush=True)
 unique=sorted(set((a['dimension'],a['degree'],a['nonlinear_total']) for a in rows));pareto=[p for p in unique if not any(q!=p and all(a<=b for a,b in zip(q,p)) for q in unique)]
 result={'status':'complete','source_sha256':sha(Path(__file__)),'input_sha256':sha(path),'subsets':512,'cases':rows,'best_by_dimension':{str(n):v for n,v in sorted(best.items())},'pareto_triples':pareto,'minimum_nonlinear_total':min(a['nonlinear_total'] for a in rows),'seconds':round(time.monotonic()-start,3),'scope':'Complete512 literal zero-carrier subsets only; no target repair or arbitrary frame/graph/counterexample optimum.'};(folder/'checkpoint.json').write_text(json.dumps(result,indent=2)+'\n');print('PARETO',pareto,flush=True)
if __name__=='__main__':main()
