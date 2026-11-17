#!/usr/bin/env python3
"""Univariate nonliteral carrier graphs and cubic-source slice equations."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement
import hashlib,json,time,sys
from verify import decode,variable,add,mul,scale,compose,derivative
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p):return [{'exponent':list(e),'coefficient':str(c)} for e,c in sorted(p.items())]
def cap(p,D):return {e:c for e,c in p.items() if sum(e[:12])<=D}
def build(pivot,other,bad,D):
 path=HERE/'cubic12-integer-witness.json';raw=json.loads(path.read_text());F=[decode(row,12) for row in raw['coefficients']];invpath=HERE/'target-derivation12/inverse.json';V=json.loads(invpath.read_text());V=[[decode(p,12) for p in row] for row in V['inverse']]
 X=[variable(i,13) for i in range(12)];t=variable(12,13)
 embed=lambda p:{e+(0,):c for e,c in p.items()}
 graph=scale(add(add(embed(F[pivot]),X[pivot],Q(-1)),mul(t,embed(F[other]))),Q(-1));assert all(e[pivot]==0 for e in graph)
 sub=list(X);sub[pivot]=graph;G={j:compose(f,sub,13) for j,f in enumerate(F) if j!=pivot}
 C={a:cap(compose({e:c for e,c in V[a][bad].items() if sum(e)<=D},sub,13),D) for a in range(12) if a!=pivot}
 def apply(f):
  out={}
  for a in C:out=add(out,cap(mul(C[a],derivative(f,a)),D))
  return out
 one={(0,)*13:Q(1)}
 for j,f in G.items():
  assert cap(apply(f),D)==(one if j==bad else {}),(pivot,j)
  if j!=bad:assert max(sum(e[:12]) for e in f)<=3
 b=add(one,C[bad],Q(-1));assert not b.get((0,)*13,0)
 ix=[j for j in range(12) if j!=pivot];mons=[];cols=[]
 for d in (2,3):
  for word in combinations_with_replacement(ix,d):
   p=one
   for j in word:p=mul(p,X[j])
   mons.append(next(iter(p))[:12]);cols.append(apply(p))
 assert len(cols)==352
 payload={'source_sha256':sha(Path(__file__)),'input_sha256':sha(path),'inverse_sha256':sha(invpath),'pivot':pivot+1,'other':other+1,'bad':bad+1,'source_order_cap':D,'parameter':'t','graph':enc(graph),'raw_outputs':{str(j+1):enc(f) for j,f in G.items()},'derivation_coefficients':{str(a+1):enc(p) for a,p in C.items()},'source_monomials':[list(e) for e in mons],'defect':enc(b),'scope':'Cubic f can be normalized to Xi_bad+homogeneous quadratic/cubic terms by subtracting the other ten already-cubic raw outputs and a constant. Therefore a slice requires the352-column source system D(q)=1-D(Xi_bad). Truncation supplies necessary equations only.'}
 return payload,cols,b

def at(p,value):
 out={}
 for e,c in p.items():
  ex=e[:12];co=c*value**e[12]
  if co:out[ex]=out.get(ex,Q(0))+co
 return {e:c for e,c in out.items() if c}
def rank_rows(columns,prime=1000000007):
 piv={};pc=[]
 for j,col in enumerate(columns):
  v={e:(c.numerator*pow(c.denominator,-1,prime))%prime for e,c in col.items()};v={e:c for e,c in v.items() if c}
  while v:
   e=max(v)
   if e not in piv:
    inv=pow(v[e],-1,prime);piv[e]={m:c*inv%prime for m,c in v.items()};pc.append(j);break
   sc=v[e]
   for m,c in piv[e].items():
    co=(v.get(m,0)-sc*c)%prime
    if co:v[m]=co
    elif m in v:del v[m]
 return list(piv),pc

def main():
 from flint import fmpq_mat
 D=int(sys.argv[1]) if len(sys.argv)>1 else 5
 for pivot,other,bad in [(9,6,1),(8,7,2)]:
  start=time.monotonic();folder=HERE/'graph-pencils12'/('pivot-%02d-D%d'%(pivot+1,D));folder.mkdir(parents=True,exist_ok=True)
  payload,cols,b=build(pivot,other,bad,D);(folder/'model.json').write_text(json.dumps(payload,indent=2)+'\n');support=set();probes=[]
  for value in map(Q,['0','1','-1','2','-2','3']):
   exact=[at(p,value) for p in cols]+[at(b,value)];rr,pc=rank_rows(exact)
   rcols=sum(j<352 for j in pc);rec={'t':str(value),'column_rank':rcols,'augmented_rank':len(pc)};probes.append(rec);print('PROBE',pivot+1,D,rec,flush=True)
   if len(pc)==353:
    A=fmpq_mat([[str(p.get(e,0)) for p in exact] for e in rr]);rhs=fmpq_mat([[int(i==352)] for i in range(353)])
    sol=A.transpose().solve(rhs);functional={e:Q(str(sol[j,0])) for j,e in enumerate(rr) if sol[j,0]};support.update(functional)
    rec['dual_support']=len(functional)
  checkpoint={'status':'probes-complete','source_sha256':sha(Path(__file__)),'model_sha256':sha(folder/'model.json'),'probes':probes,'dual_support_union':len(support),'seconds':round(time.monotonic()-start,3)}
  (folder/'checkpoint.json').write_text(json.dumps(checkpoint,indent=2)+'\n')
  if any(p['column_rank']==p['augmented_rank'] for p in probes):print('NECESSARY CUBIC SLICE CANDIDATE',pivot+1,flush=True);continue
  # Seek one polynomial-in-t dual separator on the small discovered support.
  S=sorted(support);allcols=cols+[b];grouped=[]
  for p in allcols:
   by={}
   for e,c in p.items():
    if e[:12] in support:by[(S.index(e[:12]),e[12])]=c
   grouped.append(by)
  for degree in range(4):
   variables=len(S)*(degree+1);equations={}
   for col,by in enumerate(grouped):
    for (row,power),co in by.items():
     for k in range(degree+1):
      key=(col,power+k);eq=equations.setdefault(key,{})
      unknown=row*(degree+1)+k;eq[unknown]=eq.get(unknown,Q(0))+co
   equations.setdefault((352,0),{})
   keys=sorted(equations);M=fmpq_mat([[str(equations[key].get(j,0)) for j in range(variables)]+[str(Q(key==(352,0)))] for key in keys]);rr,rank=M.rref();piv=[next(j for j in range(variables+1) if rr[i,j]) for i in range(rank)]
   print('DUAL',pivot+1,'DEG',degree,'ROWS',len(keys),'VARS',variables,'RANK',rank,'SOLVABLE',variables not in piv,flush=True)
   if variables in piv:continue
   sol=[Q(0)]*variables
   for i,j in enumerate(piv):sol[j]=Q(str(rr[i,variables]))
   weights=[]
   for i,e in enumerate(S):
    terms={k:sol[i*(degree+1)+k] for k in range(degree+1) if sol[i*(degree+1)+k]}
    if terms:weights.append({'source_exponent':list(e),'parameter_coefficients':{str(k):str(c) for k,c in terms.items()}})
   result={'status':'exact-polynomial-dual-candidate','model_sha256':sha(folder/'model.json'),'source_sha256':sha(Path(__file__)),'parameter_degree':degree,'weights':weights,'scope':'Requires independent reconstruction/multiplication of every column identity before any universal graph exclusion is accepted.'}
   (folder/'dual.json').write_text(json.dumps(result,indent=2)+'\n');checkpoint.update(status='dual-candidate',parameter_degree=degree,functional_support=len(weights));print('DUAL FOUND',pivot+1,'TERMS',len(weights),flush=True);break
  checkpoint['seconds']=round(time.monotonic()-start,3);(folder/'checkpoint.json').write_text(json.dumps(checkpoint,indent=2)+'\n')
if __name__=='__main__':main()
