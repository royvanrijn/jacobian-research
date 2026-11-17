#!/usr/bin/env python3
"""Nonzero quadratic-w tail branch for f=1-3u, after highest-w consequences."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def dec(terms,V):return sp.Add(*(sp.Rational(c)*sp.prod(v**a for v,a in zip(V,e)) for e,c in terms))
def main():
 mode=sys.argv[1];assert mode in ['generic','normalized'];source=HERE/'weighted-zero7/nu-zero/linear/input.json';d=json.loads(source.read_text());W=sp.symbols(' '.join(d['variables']));by=dict(zip(d['variables'],W));full=json.loads((HERE/'weighted-zero7/nu-zero/input.json').read_text());u,w=sp.symbols('u w');h=sp.sympify(full['quotients'][2],locals=by|{'u':u,'w':w});g=sp.sympify(full['quotients'][1],locals=by|{'u':u,'w':w});H=sp.Poly(h,w);G=sp.Poly(g,w);K=sp.expand(G.coeff_monomial(w**2).subs(u,0));C=sp.expand(H.coeff_monomial(w**2).subs(u,0));linear=[dec(f,W) for f in d['boundary_equations']];linear += [c for (k,),c in sp.Poly(H.coeff_monomial(w**2),u).terms() if k>=1]+[c for (k,),c in sp.Poly(H.coeff_monomial(w),u).terms() if k>=2]
 if mode=='normalized':linear.append(C-1)
 M,b=sp.linear_eq_to_matrix(linear,W);red,pivots=M.row_join(b).rref();assert all(i<len(W) for i in pivots);Z=tuple(W[i] for i in range(len(W)) if i not in pivots);sub={W[p]:sp.expand(red[k,len(W)]-sum(red[k,j]*W[j] for j in range(len(W)) if j not in pivots)) for k,p in enumerate(pivots)};inv=sp.Symbol('branch_inv');V=Z+(inv,);eqs=[sp.expand(dec(f,W).subs(sub,simultaneous=True)) for f in d['equations']];eqs=[f for f in eqs if f];nonzero=sp.expand((K*C).subs(sub,simultaneous=True));eqs.append(1-inv*nonzero);out=HERE/'weighted-zero7/nu-zero/high'/mode;out.mkdir(parents=True,exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),variables=list(map(str,V)),linear_equations=[enc(f,W) for f in linear],substitution={str(v):str(sub.get(v,v)) for v in W},nonzero=str(nonzero),K=str(K),C=str(C),equations=[enc(f,V) for f in eqs],scope='K!=0 highest-w branch. Boundary gives G2=K*v,H3=0; top equations force H2 constant nonzero and deg_v H1<=1. This includes all compatible tails after those necessary consequences. Normalized case uses diagonal c scaling to set H2=1, sending K toK/H2^2; no coefficient values sampled.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r=0,('+','.join(map(str,V))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-high-w-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),mode=mode,variables=len(V),seconds=round(time.monotonic()-started,3));(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True);print(text[:2000],flush=True)
if __name__=='__main__':main()
