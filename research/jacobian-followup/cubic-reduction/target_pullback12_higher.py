#!/usr/bin/env python3
"""Complete higher target degree4/5 screens, checkpointed and capped at60s."""
from pathlib import Path
from itertools import combinations_with_replacement
from fractions import Fraction as Q
import hashlib,json,time
from verify import decode,mul
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 folder=HERE/'target-pullback12-higher';folder.mkdir(exist_ok=True);path=HERE/'cubic12-integer-witness.json';raw=json.loads(path.read_text());F=[decode(r,12) for r in raw['coefficients']];prime=1000000007;start=time.monotonic()
 products={(i,):f for i,f in enumerate(F)};piv={};dependent=[];words=[];degree_records=[];terminal=False
 for degree in range(2,6):
  for ix in combinations_with_replacement(range(12),degree):
   if time.monotonic()-start>60:terminal=True;break
   value=mul(products[ix[:-1]],F[ix[-1]]);products[ix]=value;high={e:int(c)%prime for e,c in value.items() if sum(e)>3 and c%prime}
   words.append(ix);v=high
   while v:
    e=max(v)
    if e not in piv:
     inv=pow(v[e],-1,prime);piv[e]={m:c*inv%prime for m,c in v.items()};break
    co=v[e]
    for mon,c in piv[e].items():
     w=(v.get(mon,0)-co*c)%prime
     if w:v[mon]=w
     elif mon in v:del v[mon]
   if not v:dependent.append(len(words)-1)
   if len(words)%500==0:print('TARGET DEG',degree,'COLUMNS',len(words),'RANK',len(piv),'SEC',round(time.monotonic()-start,2),flush=True)
  rec={'degree':degree,'columns':len(words),'modular_rank':len(piv),'dependent_columns':len(dependent),'seconds':round(time.monotonic()-start,3)};degree_records.append(rec)
  data={'status':'time-budget-exhausted' if terminal else ('complete-through5' if degree==5 else 'running'),'source_sha256':sha(Path(__file__)),'input_sha256':sha(path),'prime':prime,'time_budget_seconds':60,'records':degree_records,'words':[[i+1 for i in ix] for ix in words],'dependent_columns':dependent,'pivot_monomials':[list(e) for e in piv],'scope':'Complete nonlinear target vector spaces through each completed total degree; a modular dependency is not an exact characteristic-zero relation. No arbitrary target degree or automorphism extension claim.'}
  (folder/'checkpoint.json').write_text(json.dumps(data,indent=2)+'\n');print('DEGREE',rec,flush=True)
  if terminal:break
 print('TERMINAL',data['status'],flush=True)
if __name__=='__main__':main()
