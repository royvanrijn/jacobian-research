#!/usr/bin/env python3
"""Coupled carrier-zero translation centres; exact primary exploration."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
from math import comb
import hashlib,json,time
from verify import decode,evaluate,derivative,add,scale
from optimize12 import serialize
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def support(F):return sum(sum(e)>1 for f in F for e in f)
def inverse(A):
 n=len(A);B=[list(row)+[Q(i==j) for j in range(n)] for i,row in enumerate(A)]
 for j in range(n):
  k=next(k for k in range(j,n) if B[k][j]);B[j],B[k]=B[k],B[j]
  c=B[j][j];B[j]=[a/c for a in B[j]]
  for k in range(n):
   if k!=j and B[k][j]:
    c=B[k][j];B[k]=[a-c*b for a,b in zip(B[k],B[j])]
 return [row[n:] for row in B]
def shift(f,centre):
 out={};n=len(centre)
 for ex,co in f.items():
  terms={(0,)*n:co}
  for j,k in enumerate(ex):
   if not k:continue
   new={}
   for ee,v in terms.items():
    for power in range(k+1):
     vv=v*comb(k,power)*centre[j]**(k-power)
     if vv:
      zz=list(ee);zz[j]=power;new[tuple(zz)]=vv
   terms=new
  out=add(out,terms)
 out.pop((0,)*n,None)
 return out
def make_map(F,J,centre):
 A=[[evaluate(f,centre) for f in row] for row in J];inv=inverse(A)
 assert all(sum(A[i][k]*inv[k][j] for k in range(12))==Q(i==j) for i in range(12) for j in range(12))
 shifted=[shift(f,centre) for f in F];G=[]
 for row in inv:
  p={}
  for co,f in zip(row,shifted):p=add(p,f,co)
  G.append(p)
 assert all(f.get(tuple(int(j==i) for j in range(12)),0)==1 for i,f in enumerate(G))
 return G,inv

def main():
 folder=HERE/'joint-translation12';folder.mkdir(exist_ok=True);path=HERE/'cubic12-integer-witness.json';raw=json.loads(path.read_text());F=[decode(row,12) for row in raw['coefficients']];J=[[derivative(f,j) for j in range(12)] for f in F];start=time.monotonic()
 vals=list(map(Q,['-1','-1/2','0','1/2','1']));bases=set(product(vals,vals,map(Q,['-8','-4','0','4'])))
 bases.update(tuple(map(Q,p[:3])) for p in raw['points']);bases.add((Q(-1,2),Q(1),Q(-2)))
 rows=[];candidates=[]
 for abc in sorted(bases):
  centre=list(abc)+[Q(0)]*9
  for i in range(3,12):centre[i]=-evaluate(F[i],centre)
  assert all(evaluate(F[i],centre)==0 for i in range(3,12))
  if sum(bool(v) for v in centre)<=1:continue # Previously completed coordinate lines.
  G,inv=make_map(F,J,centre);score=support(G)
  row={'base':list(map(str,abc)),'centre':list(map(str,centre)),'nonlinear_total':score};rows.append(row)
  candidates.append((score,centre,G,inv))
  if score<36:
   witness={'status':'unverified-smaller-candidate','input_sha256':sha(path),'source_sha256':sha(Path(__file__)),'centre':list(map(str,centre)),'matrix':[[str(v) for v in row] for row in inv],'coefficients':serialize(G),'nonlinear_total':score,'scope':'Primary exact affine transform; independent determinant/collision/conjugacy replay required.'}
   (folder/'candidate.json').write_text(json.dumps(witness,indent=2)+'\n')
  if len(rows)%10==0:
   print('CENTRES',len(rows),'BEST',min(x['nonlinear_total'] for x in rows),flush=True)
   checkpoint={'status':'running','source_sha256':sha(Path(__file__)),'input_sha256':sha(path),'centres':rows,'scope':'Finite rational carrier-zero graph centres, excluding old coordinate axes; not all translations or affine frames.'};(folder/'checkpoint.json').write_text(json.dumps(checkpoint,indent=2)+'\n')
 candidates.sort(key=lambda x:(x[0],sum(v!=0 for v in x[1]),tuple(x[1])))
 selected=[];scores={}
 for score,centre,G,inv in candidates:
  scores[score]=scores.get(score,0)+1
  if len(selected)<4:
   item={'base':list(map(str,centre[:3])),'centre':list(map(str,centre)),'matrix':[[str(v) for v in row] for row in inv],'nonlinear_total':score,'coefficients':serialize(G)};selected.append(item)
 checkpoint={'status':'finite-centres-complete','source_sha256':sha(Path(__file__)),'input_sha256':sha(path),'centres':rows,'score_histogram':scores,'best':min(scores),'selected':selected,'seconds':round(time.monotonic()-start,3),'scope':'Finite rational carrier-zero graph centres, excluding old coordinate axes; no arbitrary-translation or global affine optimum.'};(folder/'checkpoint.json').write_text(json.dumps(checkpoint,indent=2)+'\n')
 print('COMPLETE',len(rows),'HISTOGRAM',scores,'SELECTED',[(a['base'],a['nonlinear_total']) for a in selected],flush=True)
if __name__=='__main__':main()
