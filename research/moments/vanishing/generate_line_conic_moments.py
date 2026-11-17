#!/usr/bin/env python3
from pathlib import Path
from math import comb,factorial
import json,time
import sympy as S
p=Path(__file__).resolve().parent;start=time.monotonic()
vs=S.symbols('a b c d e f')
def moment(n,kind):
 terms={}
 for j in range(n+1):
  exps=(n+j,2*(n-j),j) if kind=='tangent' else (3*n-2*j,j,j)
  u,v,w=exps
  for xy in range(min(u,v)+1):
   for xz in range(min(u-xy,w)+1):
    for yz in range(min(v-xy,w-xz)+1):
     rem=[u-xy-xz,v-xy-yz,w-xz-yz]
     if any(k<0 or k%2 for k in rem):continue
     loops=[k//2 for k in rem];den=2**sum(loops)*factorial(xy)*factorial(xz)*factorial(yz)
     for k in loops:den*=factorial(k)
     num=comb(n,j)*factorial(u)*factorial(v)*factorial(w);assert num%den==0
     key=tuple(loops+[xy,xz,yz]);terms[key]=terms.get(key,0)+num//den
 return S.Poly.from_dict(terms,vs)
record={}
for kind in ['tangent','transverse']:
 polys=[];scalars=[]
 for n in [2,4,6,8,10]:
  full=moment(n,kind);content,prim=full.primitive();polys.append(str(prim.as_expr()));scalars.append(str(content));print(kind,n,len(full.terms()),flush=True)
 record[kind]={'form':'X^2*Z+X*Y^2' if kind=='tangent' else 'X*Y*Z+X^3','primitive_moments':polys,'scalars':scalars}
(p/'line_conic_moments_input.json').write_text(json.dumps({'variables':list(map(str,vs)),'indices':[2,4,6,8,10],'families':record,'elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
