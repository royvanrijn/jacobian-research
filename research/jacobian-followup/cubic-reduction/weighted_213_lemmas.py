#!/usr/bin/env python3
"""Finite parity/degree data for the mixed-weight boundary proof."""
from pathlib import Path
import hashlib,json
import sympy as sp
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def main():
 source=HERE/'weighted-sparse7/2_-1_-3/input.json';d=json.loads(source.read_text());W=sp.symbols('c5:17');by=dict(zip(map(str,W),W));u,v,t,l,sigma=sp.symbols('u v t l sigma');g=u;h=v
 for a in d['atoms']:
  if a['row']>0:
   b,c=a['exponent'][1:]
   if a['row']==1:g+=by[a['coefficient']]*u**b*v**c
   else:h+=by[a['coefficient']]*u**b*v**c
 for mode in ['linearv','parallel']:
  coords={u:l,v:t} if mode=='linearv' else {u:t,v:l-sigma*t};gh=[sp.expand(z.subs(coords,simultaneous=True)) for z in [g,h]];P=[] if mode=='linearv' else [sigma];V=tuple(P)+W;linear=[]
  for z,limit in zip(gh,[1,3]):
   for (k,),c in sp.Poly(z,t).terms():
    if k>limit:
     assert sp.Poly(c,l).degree()<2
     linear += [co for e,co in sp.Poly(c,l).terms() if co]
  M,b=sp.linear_eq_to_matrix(linear,W);red,pivots=M.row_join(b).rref();assert all(i<len(W) for i in pivots);Z=tuple(W[i] for i in range(len(W)) if i not in pivots);solve={W[p]:sp.factor(red[k,len(W)]-sum(red[k,j]*W[j] for j in range(len(W)) if j not in pivots)) for k,p in enumerate(pivots)};payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),parameters=list(map(str,P)),variables=list(map(str,W)),linear_equations=[enc(f,V) for f in linear],solution={str(k):str(v) for k,v in solve.items()},free_tail_variables=list(map(str,Z)),linear_rank=len(pivots),scope='Necessary global transverse degree bounds: boundary g*Hprime-3gprime*H=nonzero constant and degree<=5 force degG<=1,degH<=3. Above-threshold coefficients have normal degree<2 and vanish at at least two roots. Parity and univariate degree/divisibility proof remain written mathematics.');out=HERE/'weighted-213/lemmas'/mode;out.mkdir(parents=True,exist_ok=True);(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print(mode,len(pivots),len(Z))
if __name__=='__main__':main()
