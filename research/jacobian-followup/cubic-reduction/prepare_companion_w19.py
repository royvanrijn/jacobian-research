#!/usr/bin/env python3
"""A genuinely changed rank6 companion pencil, with exact finite cube coordinates."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement,product
from math import gcd,lcm
from functools import reduce
from collections import defaultdict
import hashlib,json,time
import sympy as S
from flint import fmpq_mat
HERE=Path(__file__).resolve().parent
N=19
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def primitive(v):
 v=list(map(Q,v));den=lcm(*(q.denominator for q in v));w=[int(q*den) for q in v];g=reduce(gcd,w,0)
 if not g:return (0,)*len(v),Q(0)
 if next(a for a in w if a)<0:g=-g
 return tuple(a//g for a in w),Q(g,den)**3
def cubes(v):
 d={}
 for i,j,k in combinations_with_replacement([a for a,c in enumerate(v) if c],3):
  e=[0]*N
  for a in (i,j,k):e[a]+=1
  mult=1 if i==k else (3 if i==j or j==k else 6);d[tuple(e)]=Q(mult*v[i]*v[j]*v[k])
 return d
def plus(*vs):return tuple(map(sum,zip(*vs)))
def scaled(v,k):return tuple(k*a for a in v)
def decompose(fs,mode):
 out=defaultdict(Q)
 def put(v,c):
  v,scale=primitive(v)
  if scale:out[v]+=Q(c)*scale
 if fs[0]==fs[2]:put(fs[0],1)
 elif len(set(fs))==2:
  a=next(v for v in fs if fs.count(v)==2);b=next(v for v in fs if fs.count(v)==1);put(plus(a,b),Q(1,6));put(plus(a,scaled(b,-1)),Q(-1,6));put(b,Q(-1,3))
 elif mode==1:
  a,b,c=fs
  for sb,sc in product([1,-1],repeat=2):put(plus(a,scaled(b,sb),scaled(c,sc)),Q(sb*sc,24))
 else:
  a,b,c=fs;sb,sc=[(1,1),(1,-1),(-1,1),(-1,-1)][0 if mode==0 else mode-1];a,b,c=a,scaled(b,sb),scaled(c,sc);sg=sb*sc;put(plus(a,b,c),Q(sg,6))
  for u,v in [(a,b),(a,c),(b,c)]:put(plus(u,v),Q(-sg,6))
  for v in [a,b,c]:put(v,Q(sg,6))
 return {v:c for v,c in out.items() if c}
def encode(p):return [[list(e),str(c)] for e,c in sorted(p.items()) if c]
def main():
 start=time.monotonic();out=HERE/'companion-w19';out.mkdir(exist_ok=True);base=HERE/'cubic12-integer-witness.json';d=json.loads(base.read_text());X=S.symbols('x1:20');sigma=S.Symbol('sigma');K=[sum(S.Rational(a['coefficient'])*S.prod(x**k for x,k in zip(X,a['exponent'])) for a in row) for row in d['coefficients']];E=list(X[:12]);E[7]+=sigma*X[2]*X[3];E[11]+=sigma*X[3]**2;G=[S.expand(f.subs(dict(zip(X[:12],E)),simultaneous=True)) for f in K];assert all(S.Poly(f,*X).total_degree()<=3 for f in G);Q2=[];C3=[]
 for f in G:
  pp=S.Poly(f,*X);Q2.append(sum(c*S.prod(x**k for x,k in zip(X,e)) for e,c in pp.terms() if sum(e)==2));C3.append(sum(c*S.prod(x**k for x,k in zip(X,e)) for e,c in pp.terms() if sum(e)==3))
 active=[0,1,2,3,4,9];assert all(not f for i,f in enumerate(C3) if i not in active);H=[X[18]*f for f in Q2]
 for j,i in enumerate(active):H[i]-=X[18]**2*X[12+j]
 H+= [C3[i] for i in active]+[X[18]**3];H=list(map(S.expand,H));H0=[S.expand(f.subs(sigma,0)) for f in H];H1=[S.expand(S.diff(f,sigma)) for f in H];assert all(S.expand(f-a-sigma*b)==0 for f,a,b in zip(H,H0,H1));P0=[{e:Q(c) for e,c in S.Poly(f,*X).terms() if c} for f in H0];P1=[{e:Q(c) for e,c in S.Poly(f,*X).terms() if c} for f in H1]
 oldpath=HERE/'power89/witness.json';old=json.loads(oldpath.read_text());V=S.eye(18);V[8,8]=S.Rational(1,2);V[8,7]=S.Rational(3,2);V[10,9]=3;raw=S.Matrix(old['raw_B']);transport=V.inv()*raw;forms=[];coef=[]
 def lift(v):return tuple(v[:17])+(0,v[17])
 L=S.Matrix(old['forms'])*V
 for j in range(L.rows):
  form,scale=primitive(lift([Q(c) for c in L.row(j)]));forms.append(form);column=[Q(c)*scale for c in transport.col(j)];coef.append(column[:17]+[Q(0),column[17]])
 assert len(set(forms))==len(forms)==89
 # Directly check every transported cube against the old tensor embedded in19.
 reconstructed=[{} for i in range(N)]
 for f,col in zip(forms,coef):
  for e,c in cubes(f).items():
   for i,b in enumerate(col):
    if b:reconstructed[i][e]=reconstructed[i].get(e,Q(0))+b*c
 reconstructed=[{e:c for e,c in row.items() if c} for row in reconstructed];oldtarget=[dict(p) for p in P0];extra_e=tuple([0]*17+[1,2]);assert oldtarget[9].pop(extra_e)==-1;assert reconstructed==oldtarget
 oldpool=HERE/'power-basis-pencils/pool.json';pool=set(forms)
 for f in json.loads(oldpool.read_text())['forms']:
  vec=(S.Matrix(1,18,f)*V).tolist()[0];pool.add(primitive(lift(vec))[0])
 seed_extra=[];terms=[];unit=[tuple(int(i==j) for i in range(N)) for j in range(N)]
 delta=P1+[{extra_e:Q(-1)}]
 for row,p in enumerate(delta):
  for e,c in p.items():
   fs=[unit[j] for j,k in enumerate(e) for _ in range(k)];assert len(fs)==3;terms.append(dict(row=row,coefficient=str(c),exponent=list(e)))
   for mode in range(5):pool.update(decompose(fs,mode))
   seed_extra+=list(decompose(fs,1))
 pool.discard((0,)*N);forms=forms+sorted(pool-set(forms));m=len(forms);positions={f:i for i,f in enumerate(forms)};priority=list(range(89))+[positions[f] for f in seed_extra if positions[f]>=89];priority=list(dict.fromkeys(priority));priority+=[i for i in range(m) if i not in priority];polys=[cubes(f) for f in forms];monomials=sorted(set().union(*(set(p) for p in polys+P0+P1)));C=fmpq_mat([[str(p.get(e,0)) for p in polys] for e in monomials]);T=fmpq_mat([[str(p.get(e,0)) for p in P0+P1] for e in monomials]);ordered=fmpq_mat([[str(C[i,j]) for j in priority] for i in range(len(monomials))]);rr,rank=ordered.rref();basis=[priority[next(j for j in range(m) if rr[i,j])] for i in range(rank)];assert basis[:89]==list(range(89));joint=fmpq_mat([[str(C[i,j]) for j in basis]+[str(C[i,j]) for j in range(m)]+[str(T[i,j]) for j in range(38)] for i in range(len(monomials))]);rr,r=joint.rref();assert r==rank
 A=[{str(j):str(rr[i,rank+j]) for j in range(m) if rr[i,rank+j]} for i in range(rank)];B=[{str(j):str(rr[i,rank+m+j]) for j in range(38) if rr[i,rank+m+j]} for i in range(rank)];padding=positions[(0,)*18+(1,)];assert all(Q(row.get('18','0'))==int(basis[i]==padding) and not Q(row.get('37','0')) for i,row in enumerate(B))
 payload=dict(source_sha256=sha(Path(__file__)),base_sha256=sha(base),old89_sha256=sha(oldpath),old_pool_sha256=sha(oldpool),variables=list(map(str,X)),dimension=19,axis='w',parameter='sigma',excluded_parameter_values=['0'],active_cubic_outputs=[i+1 for i in active],source_shear=list(map(str,E)),base_family=list(map(str,G)),H0=[encode(p) for p in P0],H1=[encode(p) for p in P1],forms=forms,padding_index=padding,rank=rank,monomials=len(monomials),difference_terms=terms,scope='Changed rank6 companion tensor H_sigma, sigma!=0, from two commuting quadratic source shears of KCUB3. The old89 identity is transported and independently checked at preparation, then extended by polarization cubes. Finite-pool and primary exact matrices only; no improved Keller witness or minimum claim.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');(out/'basis.json').write_text(json.dumps(dict(input_sha256=sha(out/'input.json'),basis=basis,A=A,B=B),indent=2,sort_keys=True)+'\n');print('PREPARED forms',m,'rank',rank,'monomials',len(monomials),'seconds',round(time.monotonic()-start,3),flush=True)
if __name__=='__main__':main()
