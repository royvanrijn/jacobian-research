#!/usr/bin/env python3
"""Extract polynomial dual pairings and exact exceptional parameter roots."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib
import sympy as s
from flint import fmpq_mat
from verify import decode,variable,mul,add,derivative
from graph_pencil_derivations12 import cap,at,rank_rows
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 folder=HERE/'graph-pencils12/pivot-10-D7';model=json.loads((folder/'model.json').read_text());D=7
 C={int(a)-1:decode(row,13) for a,row in model['derivation_coefficients'].items()};b=decode(model['defect'],13);cols=[]
 for ex in model['source_monomials']:
  q={tuple(ex)+(0,):Q(1)};out={}
  for a,p in C.items():out=add(out,cap(mul(p,derivative(q,a)),D))
  cols.append(out)
 S=set()
 for point in map(Q,['0','1','-1','2','-2','3']):
  exact=[at(c,point) for c in cols]+[at(b,point)];rows,piv=rank_rows(exact);assert len(piv)==353
  A=fmpq_mat([[str(p.get(e,0)) for p in exact] for e in rows]);v=A.transpose().solve(fmpq_mat([[int(i==352)] for i in range(353)]))
  S.update(e for i,e in enumerate(rows) if v[i,0])
 S=sorted(S);index={e:i for i,e in enumerate(S)};t=s.Symbol('t');records=[]
 for degree in [1,2,3]:
  variables=len(S)*(degree+1);eqs={}
  for j,c in enumerate(cols):
   for e,co in c.items():
    if e[:12] not in index:continue
    row=index[e[:12]]
    for k in range(degree+1):
     rec=eqs.setdefault((j,e[12]+k),{});idx=row*(degree+1)+k;rec[idx]=rec.get(idx,Q(0))+co
  M=fmpq_mat([[str(eq.get(j,0)) for j in range(variables)] for eq in eqs.values()]);rr,rank=M.rref();piv=[next(j for j in range(variables) if rr[i,j]) for i in range(rank)];free=[j for j in range(variables) if j not in piv];pairs=[]
  for f in free:
   v=[Q(0)]*variables;v[f]=1
   for i,j in enumerate(piv):v[j]=-Q(str(rr[i,f]))
   pairing={}
   for e,c in b.items():
    if e[:12] not in index:continue
    row=index[e[:12]]
    for k in range(degree+1):
     if v[row*(degree+1)+k]:pairing[e[12]+k]=pairing.get(e[12]+k,Q(0))+c*v[row*(degree+1)+k]
   pairing={k:c for k,c in pairing.items() if c}
   if not pairing:continue
   weights=[]
   for i,e in enumerate(S):
    terms={str(k):str(v[i*(degree+1)+k]) for k in range(degree+1) if v[i*(degree+1)+k]}
    if terms:weights.append({'source_exponent':list(e),'parameter_coefficients':terms})
   pairs.append({'weights':weights,'pairing':{str(k):str(c) for k,c in pairing.items()}})
  gcd=s.Poly(0,t)
  for pair in pairs:gcd=s.gcd(gcd,s.Poly(sum(s.Rational(c)*t**int(k) for k,c in pair['pairing'].items()),t))
  result={'parameter_degree':degree,'unknowns':variables,'rank':rank,'nonzero_pairings':len(pairs),'gcd':str(gcd.monic().as_expr()) if gcd else '0','pairs':pairs};records.append(result);print('DEG',degree,'GCD',result['gcd'],'PAIRINGS',len(pairs),flush=True)
  if gcd and gcd.degree()<=1:break
 payload={'status':'homogeneous-dual-candidates','source_sha256':sha(Path(__file__)),'model_sha256':sha(folder/'model.json'),'records':records,'scope':'All identities require independent coefficient replay. A pairing polynomial excludes only its nonzero locus until its roots are handled.'};(folder/'dual-gcd.json').write_text(json.dumps(payload,indent=2)+'\n')
 if gcd and gcd.degree()==1:
  root=-gcd.nth(0)/gcd.nth(1);point=Q(str(root));rr,pc=rank_rows([at(c,point) for c in cols]+[at(b,point)]);print('EXCEPTION',root,'RANK',sum(j<352 for j in pc),'AUG',len(pc),flush=True)
if __name__=='__main__':main()
