#!/usr/bin/env python3
"""Finite source-degree data for the two written univariate-boundary lemmas."""
from pathlib import Path
import hashlib,json
import sympy as sp
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def main():
 source=HERE/'weighted-zero7/both-zero/input.json';d=json.loads(source.read_text());W=sp.symbols(' '.join(n for n in d['variables'] if n[0]!='f'));by=dict(zip(map(str,W),W));u,w,t,l,sigma=sp.symbols('u w t l sigma');g,h=[sp.sympify(f,locals=by|{'u':u,'w':w}) for f in d['quotients'][1:]]
 for mode in ['parallel','univariate-u']:
  gh=[sp.expand(z.subs({u:t,w:l-sigma*t},simultaneous=True)) for z in [g,h]] if mode=='parallel' else [sp.expand(z.subs({u:l,w:t},simultaneous=True)) for z in [g,h]];parameters=[sigma] if mode=='parallel' else [];V=tuple(parameters)+W;linear=[]
  for z,limit in zip(gh,[1,2]):
   for (k,),c in sp.Poly(z,t).terms():
    if k>limit:
     assert sp.Poly(c,l).degree()<2
     linear += [co for e,co in sp.Poly(c,l).terms() if co]
  M,b=sp.linear_eq_to_matrix(linear,W);red,pivots=M.row_join(b).rref();assert all(i<len(W) for i in pivots);Z=tuple(W[i] for i in range(len(W)) if i not in pivots);solve={W[p]:sp.factor(red[k,len(W)]-sum(red[k,j]*W[j] for j in range(len(W)) if j not in pivots)) for k,p in enumerate(pivots)};payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),mode=mode,parameters=list(map(str,parameters)),variables=list(map(str,W)),linear_equations=[enc(z,V) for z in linear],linear_solution={str(k):str(v) for k,v in solve.items()},free_tail_variables=list(map(str,Z)),linear_rank=len(pivots),scope='Necessary global degree constraints when a squarefree univariate boundary has at least two roots. Above-threshold coefficients have degree<2 in the normal coordinate, so vanishing on all boundary lines makes them identically zero. Written divisibility argument then excludes a Keller completion.');out=HERE/'weighted-zero7/lemmas'/mode;out.mkdir(parents=True,exist_ok=True);(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print(mode,len(pivots),len(Z))
if __name__=='__main__':main()
