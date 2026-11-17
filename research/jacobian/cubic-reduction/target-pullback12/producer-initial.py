#!/usr/bin/env python3
"""Nonlinear degree<=3 target polynomials whose pullback by K remains cubic."""
from pathlib import Path
from itertools import combinations_with_replacement
from fractions import Fraction as Q
import json,hashlib,time
from verify import decode,mul
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 out=HERE/'target-pullback12';out.mkdir(exist_ok=True);path=HERE/'cubic12-integer-witness.json';raw=json.loads(path.read_text());F=[decode(row,12) for row in raw['coefficients']];start=time.monotonic()
 words=[ix for d in (2,3) for ix in combinations_with_replacement(range(12),d)];cols=[];products={}
 for ix in words:
  if len(ix)==2:p=mul(F[ix[0]],F[ix[1]]);products[ix]=p
  else:p=mul(products[ix[:2]],F[ix[2]])
  cols.append({e:int(c) for e,c in p.items() if sum(e)>3})
 mons=sorted(set().union(*(set(p) for p in cols)));prime=1000000007;pivots={};dependent=[]
 for j,col in enumerate(cols):
  v={e:c%prime for e,c in col.items() if c%prime}
  while v:
   e=max(v)
   if e not in pivots:
    inv=pow(v[e],-1,prime);pivots[e]={m:c*inv%prime for m,c in v.items()};break
   co=v[e]
   for m,c in pivots[e].items():
    value=(v.get(m,0)-co*c)%prime
    if value:v[m]=value
    elif m in v:del v[m]
  if not v:dependent.append(j)
  if (j+1)%100==0:print('COLS',j+1,'RANK',len(pivots),flush=True)
 data={'status':'modular-screen-complete','source_sha256':sha(Path(__file__)),'input_sha256':sha(path),'target_degree_bound':3,'columns':len(cols),'source_high_monomials':len(mons),'prime':prime,'modular_rank':len(pivots),'dependent_columns':dependent,'words':[[j+1 for j in ix] for ix in words],'seconds':round(time.monotonic()-start,3),'scope':'Complete nonlinear target degree2/3 vector space in the fixed K frame. Modular rank is a rational lower bound; a dependent modular column is not yet an exact polynomial relation. No arbitrary target degree or sequence coverage.'}
 (out/'checkpoint.json').write_text(json.dumps(data,indent=2)+'\n');print(json.dumps({k:v for k,v in data.items() if k not in ['words','dependent_columns']}),flush=True)
 if dependent:
  from flint import fmpq_mat
  M=fmpq_mat([[col.get(e,0) for col in cols] for e in mons]);rr,rank=M.rref();pp=[next(j for j in range(len(cols)) if rr[i,j]) for i in range(rank)];free=[j for j in range(len(cols)) if j not in pp];basis=[]
  for j in free:
   vec={j:Q(1)}
   for row,k in enumerate(pp):
    co=-Q(str(rr[row,j]))
    if co:vec[k]=co
   basis.append([{'word':[a+1 for a in words[k]],'coefficient':str(c)} for k,c in sorted(vec.items())])
  data.update(status='exact-kernel-complete',rational_rank=rank,kernel=basis,kernel_dimension=len(basis));(out/'checkpoint.json').write_text(json.dumps(data,indent=2)+'\n');print('KERNEL',len(basis),json.dumps(basis)[:10000],flush=True)
if __name__=='__main__':main()
