#!/usr/bin/env python3
"""Close the uncovered t=3/8 locus by a polynomial dual, no denominator drop."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
from flint import fmpq_mat
from verify import decode,add,mul,derivative
from graph_pencil_derivations12 import cap,at,rank_rows
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();folder=HERE/'graph-pencils12/pivot-10-D7';model=json.loads((folder/'model.json').read_text());D=7
 C={int(a)-1:decode(row,13) for a,row in model['derivation_coefficients'].items()};b=decode(model['defect'],13);cols=[]
 for ex in model['source_monomials']:
  q={tuple(ex)+(0,):Q(1)};p={}
  for a,c in C.items():p=add(p,cap(mul(c,derivative(q,a)),D))
  cols.append(p)
 S=set();probes=[]
 for point in map(Q,['0','1','-1','2','-2','3','3/8']):
  exact=[at(c,point) for c in cols]+[at(b,point)];rows,pp=rank_rows(exact);assert len(pp)==353
  M=fmpq_mat([[str(p.get(e,0)) for p in exact] for e in rows]);v=M.transpose().solve(fmpq_mat([[int(i==352)] for i in range(353)]));support=[e for i,e in enumerate(rows) if v[i,0]];S.update(support);probes.append({'parameter':str(point),'rank':353,'dual_support':len(support)})
 S=sorted(S);index={e:i for i,e in enumerate(S)};grouped=[]
 for p in cols+[b]:
  grouped.append({(index[e[:12]],e[12]):c for e,c in p.items() if e[:12] in index})
 log=[]
 for degree in range(7):
  n=len(S)*(degree+1);eqs={}
  for j,p in enumerate(grouped):
   for (row,power),c in p.items():
    for k in range(degree+1):
     rec=eqs.setdefault((j,power+k),{});idx=row*(degree+1)+k;rec[idx]=rec.get(idx,Q(0))+c
  eqs.setdefault((352,0),{});keys=sorted(eqs)
  M=fmpq_mat([[str(eqs[key].get(j,0)) for j in range(n)]+[str(int(key==(352,0)))] for key in keys]);rr,rank=M.rref();piv=[next(j for j in range(n+1) if rr[i,j]) for i in range(rank)]
  row={'parameter_degree':degree,'equation_rows':len(keys),'unknowns':n,'augmented_rank':rank,'consistent':n not in piv};log.append(row);print('ADAPTIVE',row,flush=True)
  if n in piv:continue
  sol=[Q(0)]*n
  for i,j in enumerate(piv):sol[j]=Q(str(rr[i,n]))
  weights=[]
  for i,e in enumerate(S):
   terms={str(k):str(sol[i*(degree+1)+k]) for k in range(degree+1) if sol[i*(degree+1)+k]}
   if terms:weights.append({'source_exponent':list(e),'parameter_coefficients':terms})
  out={'status':'exact-polynomial-dual-candidate','source_sha256':sha(Path(__file__)),'model_sha256':sha(folder/'model.json'),'parameter_degree':degree,'weights':weights,'probe_cover':probes,'attempts':log,'seconds':round(time.monotonic()-start,3),'scope':'The t=3/8 row support supplements an earlier generic pairing. This global polynomial dual requires independent reconstruction and multiplication before acceptance.'};(folder/'dual-adaptive.json').write_text(json.dumps(out,indent=2)+'\n');print('FOUND',len(weights),'DEG',degree,flush=True);break
 else:
  (folder/'dual-adaptive.json').write_text(json.dumps({'status':'no-dual-within-bound','source_sha256':sha(Path(__file__)),'model_sha256':sha(folder/'model.json'),'attempts':log,'probes':probes},indent=2)+'\n')
if __name__=='__main__':main()
