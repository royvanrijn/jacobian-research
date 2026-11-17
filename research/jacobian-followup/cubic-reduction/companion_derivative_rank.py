#!/usr/bin/env python3
"""Exact Laurent-unit Gaussian rank certificates for one-parameter tensors."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,sys
from flint import fmpq_mat
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def plus(a,b,c=Q(1)):
 out=dict(a)
 for k,v in b.items():out[k]=out.get(k,Q(0))+c*v
 return {k:v for k,v in out.items() if v}
def times(a,b):
 out={}
 for k,v in a.items():
  for l,w in b.items():out[k+l]=out.get(k+l,Q(0))+v*w
 return {k:v for k,v in out.items() if v}
def derivative(p,j):
 out={}
 for e,c in p.items():
  if e[j]:a=list(e);a[j]-=1;out[tuple(a)]=c*e[j]
 return out

def run(tag):
 folder=HERE/tag;ip=folder/'input.json';d=json.loads(ip.read_text());H=[[{tuple(e):Q(c) for e,c in row} for row in d[key]] for key in ['H0','H1']];derivatives=[(i,j,derivative(H[0][i],j),derivative(H[1][i],j)) for i in range(19) for j in range(19)];mons=sorted(set().union(*(set(p) for i,j,a,b in derivatives for p in [a,b])));initial=[];labels=[]
 for i,j,a,b in derivatives:
  row={}
  for k,e in enumerate(mons):
   p={0:a[e]} if e in a else {}
   if e in b:p[1]=b[e]
   if p:row[k]=p
  if row:initial.append(row);labels.append([i,j])
 rows=[{k:dict(p) for k,p in row.items()} for row in initial];used=set();pivots=[]
 for step in range(min(len(rows),len(mons))):
  choices=[(abs(next(iter(p))),len(rows[i]),i,k) for i in range(step,len(rows)) for k,p in rows[i].items() if k not in used and len(p)==1]
  if not choices:break
  _,_,i,k=min(choices);rows[step],rows[i]=rows[i],rows[step];power,c=next(iter(rows[step][k].items()));pivots.append(dict(row=i,column=k,power=power,coefficient=str(c)));rows[step]={j:{e-power:v/c for e,v in p.items()} for j,p in rows[step].items()};assert rows[step][k]=={0:Q(1)}
  for j in range(step+1,len(rows)):
   factor=rows[j].get(k)
   if factor:
    for col,p in rows[step].items():
     value=plus(rows[j].get(col,{}),times(factor,p),Q(-1))
     if value:rows[j][col]=value
     else:rows[j].pop(col,None)
  used.add(k)
 rank=len(pivots);zero=all(not row for row in rows[rank:]);numeric=fmpq_mat([[str(row.get(k,{}).get(0,Q(0))) for k in range(len(mons))] for row in initial]);payload=dict(status='exact-unit-rank-certificate' if zero else 'nonunit-residual-unresolved',source_sha256=sha(Path(__file__)),input_sha256=sha(ip),parameter='sigma',excluded_parameter_values=['0'],monomials=[list(e) for e in mons],row_labels=labels,pivots=pivots,rank=rank,zero_residual=zero,rank_at_zero_primary=int(numeric.rank()),scope='Univariate Laurent-unit row operations: every pivot is a rational monomial in sigma, so the complete zero residual certifies the exact derivative-space rank for every nonzero sigma. Independent derivative reconstruction and arithmetic replay required.');(folder/'derivative-rank-certificate.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print(tag,payload['status'],'NONZERO RANK',rank,'ZERO RANK',payload['rank_at_zero_primary'],flush=True)
def main():
 for tag in sys.argv[1:] or ['companion-pencil19','companion-w19']:run(tag)
if __name__=='__main__':main()
