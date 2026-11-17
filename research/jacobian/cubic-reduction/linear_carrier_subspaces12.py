#!/usr/bin/env python3
"""Complete echelon-chart screen of arbitrary linear carrier subspaces."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from collections import defaultdict
import hashlib,json,argparse,time
from verify import decode,variable,compose,add,scale,mul
HERE=Path(__file__).resolve().parent
WEIGHTS=[1,-1,-2,0,-3,2,-1,-2,-2,-1,-1,0]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p):return [{'exponent':list(e),'coefficient':str(c)} for e,c in sorted(p.items())]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--retained',type=int,default=3);ap.add_argument('--degree',type=int,default=4);ap.add_argument('--seconds',type=float,default=60);args=ap.parse_args()
 assert 0<=args.retained<=9 and args.degree>=3
 path=HERE/'cubic12-integer-witness.json';raw=json.loads(path.read_text());F=[decode(p,12) for p in raw['coefficients']]
 folder=HERE/('linear-subspaces-n%02d-d%d'%(args.retained+3,args.degree));folder.mkdir(exist_ok=True);rows=[];start=time.monotonic();sets=list(combinations(range(3,12),args.retained))
 def checkpoint(status):
  d={'status':status,'source_sha256':sha(Path(__file__)),'input_sha256':sha(path),'retained_carrier_dimension':args.retained,'target_degree':args.degree,'charts_total':len(sets),'charts_completed':len(rows),'cases':rows,'seconds':round(time.monotonic()-start,3),'scope':'Complete highest-pivot echelon charts for simultaneous linear zero relations without weight restrictions; no target repair, base-output relations, nonlinear graphs or global lower bound.'};(folder/'checkpoint.json').write_text(json.dumps(d,indent=2)+'\n')
 for case,kept in enumerate(sets):
  if time.monotonic()-start>=args.seconds:checkpoint('time-cap');print('TIME CAP',len(rows),flush=True);return
  deleted=[j for j in range(3,12) if j not in kept];pairs=[(j,k) for j in deleted for k in kept if k<j];n=12+len(pairs);X=[variable(j,n) for j in range(n)];E=X[:12];embed=lambda p:{e+(0,)*len(pairs):c for e,c in p.items()}
  for j in deleted:
   correction=compose(embed(add(F[j],variable(j,12),Q(-1))),E,n)
   for z,(a,k) in enumerate(pairs):
    if a==j:correction=add(correction,mul(X[12+z],compose(embed(F[k]),E,n)))
   assert all(e[j]==0 for e in correction);E[j]=scale(correction,-1)
  full=[compose(embed(p),E,n) for p in F];retained=[0,1,2]+list(kept);G={i:full[i] for i in retained}
  # Every deleted row equals minus its indicated linear combination of retained carriers.
  for j in deleted:
   relation=full[j]
   for z,(a,k) in enumerate(pairs):
    if a==j:relation=add(relation,mul(X[12+z],full[k]))
   assert not relation
  high={}
  for i,p in G.items():
   groups=defaultdict(dict)
   for e,c in p.items():
    if sum(e[:12])>args.degree:groups[e[:12]][e[12:]]=c
   for e,q in groups.items():high[(i,e)]=q
  constant=next(((i,e,q) for (i,e),q in high.items() if set(q)=={(0,)*len(pairs)}),None)
  status='constant-obstruction' if constant else ('universal-degree-family' if not high else 'parameter-equations-unresolved')
  row={'case':case,'kept':[j+1 for j in kept],'deleted':[j+1 for j in deleted],'parameter_pairs':[[j+1,k+1] for j,k in pairs],'parameter_count':len(pairs),'high_coefficient_equations':len(high),'generic_degree':max(sum(e[:12]) for p in G.values() for e in p),'status':status}
  if constant:
   i,e,q=constant;row['constant_obstruction']={'output':i+1,'source_exponent':list(e),'coefficient':str(q[(0,)*len(pairs)])}
  model={'input_sha256':sha(path),'source_sha256':sha(Path(__file__)),**row,'source_embedding':[enc(p) for p in E],'outputs':{str(i+1):enc(p) for i,p in G.items()},'high_coefficients':[{'output':i+1,'source_exponent':list(e),'parameter_polynomial':enc(q)} for (i,e),q in sorted(high.items())]}
  (folder/('case-%03d.json'%case)).write_text(json.dumps(model,indent=2)+'\n');rows.append(row);checkpoint('running')
  if not constant:print('NONTRIVIAL',row,flush=True)
 checkpoint('complete');counts={s:sum(r['status']==s for r in rows) for s in sorted(set(r['status'] for r in rows))};print('COMPLETE',counts,'seconds',round(time.monotonic()-start,3),flush=True)
if __name__=='__main__':main()
