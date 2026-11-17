"""Measure bounded exact determinant replay before choosing full-grid cap."""
from pathlib import Path
import json,time,sympy as S
from flint import fmpz_poly,fmpz_mat
HERE=Path(__file__).resolve().parent
src=json.loads((HERE/'generic-triple-input.json').read_text());delta=json.loads((HERE/'residual-discriminant-candidate.json').read_text());d,a,s=S.symbols('d a s');ct=S.Poly(S.sympify(src['equations'][1]['polynomial']),d,a,s).terms()
C=[{} for _ in range(9)]
for (i,j,k),c in ct:C[i].setdefault(j,{})[k]=int(c)
def unpack(rows):
 out=[]
 for row in rows:
  out.append([fmpz_poly([row.get(j,{}).get(k,0) for k in range(max(row.get(j,{0:0}))+1)]) for j in range(max(row)+1)])
 return out
C=unpack(C);D={}
for (i,j,k),c in delta['terms']:D.setdefault(j,{})[k]=int(c)
D=unpack([D])[0]
def bareiss(M):
 W=[r[:] for r in M];prev=1;sign=1
 for k in range(len(W)-1):
  v=next((i for i in range(k,len(W)) if W[i][k]),None)
  if v is None:return 0
  if v!=k:W[k],W[v]=W[v],W[k];sign=-sign
  pivot=W[k][k]
  for i in range(k+1,len(W)):
   left=W[i][k]
   for j in range(k+1,len(W)):
    num=W[i][j]*pivot-left*W[k][j];q,r=divmod(num,prev);assert r==0;W[i][j]=q
   W[i][k]=0
  prev=pivot
 return sign*W[-1][-1]
def sylv(f):
 g=[(j+1)*f[j+1] for j in range(8)];M=[[0]*15 for _ in range(15)]
 for i in range(7):
  for j,c in enumerate(reversed(f)):M[i][i+j]=c
 for i in range(8):
  for j,c in enumerate(reversed(g)):M[7+i][i+j]=c
 return M
start=time.monotonic();timings={'native':0,'plain':0}
for sv in [-202,-1,1,203]:
 ca=[fmpz_poly([int(P(sv)) for P in row]) for row in C];da=fmpz_poly([int(P(sv)) for P in D])
 for av in range(-10,15):
  f=[int(P(av)) for P in ca];M=sylv(f);expect=f[-1]*int(da(av));t=time.monotonic();n=int(fmpz_mat(M).det());timings['native']+=time.monotonic()-t;assert n==expect
  t=time.monotonic();b=bareiss(M);timings['plain']+=time.monotonic()-t;assert b==expect
out={'samples':100,'timings':timings,'total_seconds':time.monotonic()-start,'extrapolated_43036_seconds':{k:v*430.36 for k,v in timings.items()},'theorem':False}
(HERE/'discriminant-replay-profile.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
