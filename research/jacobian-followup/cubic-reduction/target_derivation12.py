#!/usr/bin/env python3
"""Polynomial inverse-Jacobian derivations on the complete source cubic space."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement
import json,hashlib,time
from verify import decode,derivative,add,mul,scale,identity,madd,mmul,determinant3,variable
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def inverse_jacobian(F):
 n=12;J=[[derivative(f,j) for j in range(n)] for f in F]
 A=[r[:3] for r in J[:3]];B=[r[3:] for r in J[:3]];C=[r[:3] for r in J[3:]];D=[r[3:] for r in J[3:]]
 eye=identity(9,n);N=madd(D,eye,Q(-1));N2=mmul(N,N);assert not any(p for r in mmul(N2,N) for p in r)
 Di=madd(madd(eye,N,Q(-1)),N2);E=mmul(Di,C);T=mmul(B,Di);S=madd(A,mmul(B,E),Q(-1));assert determinant3(S,n)=={(0,)*n:Q(1)}
 P=[]
 for i in range(3):
  row=[]
  for j in range(3):
   rr=[a for a in range(3) if a!=j];cc=[a for a in range(3) if a!=i]
   co=add(mul(S[rr[0]][cc[0]],S[rr[1]][cc[1]]),mul(S[rr[0]][cc[1]],S[rr[1]][cc[0]]),Q(-1))
   row.append(scale(co,Q((-1)**(i+j))))
  P.append(row)
 EP=mmul(E,P);PT=mmul(P,T);bottom=madd(Di,mmul(EP,T))
 inv=[P[i]+[scale(p,-1) for p in PT[i]] for i in range(3)]+[[scale(p,-1) for p in EP[i]]+bottom[i] for i in range(9)]
 assert mmul(J,inv)==identity(n,n) and mmul(inv,J)==identity(n,n)
 return J,inv

def main():
 folder=HERE/'target-derivation12';folder.mkdir(exist_ok=True);path=HERE/'cubic12-integer-witness.json';raw=json.loads(path.read_text());F=[decode(r,12) for r in raw['coefficients']];start=time.monotonic();J,inv=inverse_jacobian(F)
 print('INVERSE',sum(len(p) for row in inv for p in row),'TERMS','MAXDEG',max(sum(e) for row in inv for p in row for e in p),flush=True)
 witness={'input_sha256':sha(path),'source_sha256':sha(Path(__file__)),'inverse':[[[{'exponent':list(e),'coefficient':str(c)} for e,c in sorted(p.items())] for p in row] for row in inv]};(folder/'inverse.json').write_text(json.dumps(witness,indent=2)+'\n')
 words=[ix for d in range(4) for ix in combinations_with_replacement(range(12),d)];exps=[]
 for ix in words:
  ex=[0]*12
  for j in ix:ex[j]+=1
  exps.append(tuple(ex))
 data={'status':'running','input_sha256':sha(path),'source_sha256':sha(Path(__file__)),'source_degree_bound':3,'source_columns':len(words),'inverse_sha256':sha(folder/'inverse.json'),'prime':1000000007,'axes':[],'scope':'Kernel of each inverse-Jacobian column derivation on the complete degree<=3 source space. Kernel is only a necessary condition for membership in the other-output polynomial algebra. Full-rank upper bounds plus known independent raw output invariants can prove equality; no arbitrary-frame claim.'}
 (folder/'checkpoint.json').write_text(json.dumps(data,indent=2)+'\n');prime=data['prime']
 for i in range(12):
  piv={};dep=[];cols=[];pivot_columns=[]
  for k,ex in enumerate(exps):
   col={}
   for a,ea in enumerate(ex):
    if not ea:continue
    ee=list(ex);ee[a]-=1
    for power,c in inv[a][i].items():
     exponent=tuple(u+v for u,v in zip(ee,power));col[exponent]=col.get(exponent,Q(0))+ea*c
   col={e:c for e,c in col.items() if c};cols.append(col)
   v={e:(c.numerator*pow(c.denominator,-1,prime))%prime for e,c in col.items() if c%prime}
   while v:
    e=max(v)
    if e not in piv:
     sc=pow(v[e],-1,prime);piv[e]={m:c*sc%prime for m,c in v.items()};pivot_columns.append(k);break
    sc=v[e]
    for m,c in piv[e].items():
     co=(v.get(m,0)-sc*c)%prime
     if co:v[m]=co
     elif m in v:del v[m]
   if not v:dep.append(k)
  rec={'axis':i+1,'modular_rank':len(piv),'kernel_upper_bound':len(words)-len(piv),'known_kernel_dimension':12,'dependent_column_count':len(dep),'pivot_monomials':[list(e) for e in piv],'pivot_columns':pivot_columns,'seconds':round(time.monotonic()-start,3)}
  data['axes'].append(rec);(folder/'checkpoint.json').write_text(json.dumps(data,indent=2)+'\n');print('AXIS',i+1,'RANK',len(piv),'KERNEL<=',len(words)-len(piv),flush=True)
  if len(piv)<443:
   from flint import fmpq_mat
   mons=sorted(set().union(*(set(p) for p in cols)));M=fmpq_mat([[str(p.get(e,0)) for p in cols] for e in mons]);rr,rank=M.rref();pp=[next(j for j in range(len(cols)) if rr[a,j]) for a in range(rank)];free=[j for j in range(len(cols)) if j not in pp];basis=[]
   for j in free:
    vec={j:Q(1)}
    for a,k in enumerate(pp):
     co=-Q(str(rr[a,j]))
     if co:vec[k]=co
    basis.append([{'exponent':list(exps[k]),'coefficient':str(c)} for k,c in vec.items()])
   rec.update(rational_rank=rank,rational_kernel_dimension=len(basis));(folder/('axis-%02d-kernel.json'%(i+1))).write_text(json.dumps({'axis':i+1,'basis':basis},indent=2)+'\n');data['axes'][-1]=rec
 data.update(status='complete',seconds=round(time.monotonic()-start,3));(folder/'checkpoint.json').write_text(json.dumps(data,indent=2)+'\n')
 print('COMPLETE',[(r['axis'],r['kernel_upper_bound']) for r in data['axes']],flush=True)
if __name__=='__main__':main()
