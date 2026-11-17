#!/usr/bin/env python3
from pathlib import Path
from functools import lru_cache
from math import comb
import json,time
import sympy as S
p=Path(__file__).resolve().parent;start=time.monotonic()
a,b,c,d,e,f=S.symbols('a b c d e f');vars=[a,b,c,d,e,f];cov=[[a,d,e],[d,b,f],[e,f,c]]
@lru_cache(None)
def wick(v):
 if not any(v):return S.Integer(1)
 if sum(v)%2:return S.Integer(0)
 i=next(j for j,n in enumerate(v) if n);q=list(v);q[i]-=1;out=0
 for j,n in enumerate(q):
  if n:
   w=q[:];w[j]-=1;out+=n*cov[i][j]*wick(tuple(w))
 return S.Poly(out,*vars).as_expr()
moments=[]
for n in [2,4,6,8,10]:
 v=S.Poly(sum(comb(n,j)*wick((2*j,3*(n-j),j)) for j in range(n+1)),*vars)
 den,z=v.clear_denoms();cont,prim=z.primitive();expr=prim.as_expr();moments.append(str(expr));print('MOMENT',n,'TERMS',len(v.terms()),'DEGREE',v.total_degree(),flush=True)
(p/'cuspidal_moments_ten_input.json').write_text(json.dumps({'normal_form':'X^2 Z + Y^3','variables':list(map(str,vars)),'even_moment_indices':[2,4,6,8,10],'primitive_moments':moments,'nonzero_integer_scalars':[str(S.Poly(sum(comb(n,j)*wick((2*j,3*(n-j),j)) for j in range(n+1)),*vars).content()) for n in [2,4,6,8,10]],'elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
