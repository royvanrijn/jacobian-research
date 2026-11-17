#!/usr/bin/env python3
"""Determinant-only residual tail strata, without fixed collision conditions."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def dec(terms,V):return sp.Add(*(sp.Rational(c)*sp.prod(v**a for v,a in zip(V,e)) for e,c in terms))
def main():
 mode=sys.argv[1];assert mode in ['Aone','Bzero','raw13'];source=HERE/'weighted-unmarked7/linear/input.json';d=json.loads(source.read_text());W=sp.symbols(' '.join(d['variables']),seq=True);by=dict(zip(d['variables'],W));full=json.loads((HERE/'weighted-unmarked7/input.json').read_text());u,w=sp.symbols('u w');allby=by|{n:sp.Integer(0) for n in ['f1','f2','f3','f4']}|{'u':u,'w':w};g=sp.sympify(full['quotients'][1],locals=allby);G=sp.Poly(sp.expand(g.subs(w,1-3*u)),u);A=G.coeff_monomial(u);B=G.coeff_monomial(1);linear=[dec(f,W) for f in d['boundary_equations']]
 if mode=='Aone':linear.append(A-1)
 elif mode=='Bzero':linear.append(B)
 else:
  raw=json.loads((HERE/'raw13-specialization/input.json').read_text());raw_atoms={(a['row'],tuple(a['exponent'])) for a in raw['atoms']}
  linear += [by[a['coefficient']] for a in full['atoms'] if a['row']>0 and (a['row'],tuple(a['exponent'])) not in raw_atoms]
 M,b=sp.linear_eq_to_matrix(linear,W);R,pivots=M.row_join(b).rref();assert all(i<len(W) for i in pivots);Z=tuple(W[i] for i in range(len(W)) if i not in pivots);sub={W[p]:sp.expand(R[k,len(W)]-sum(R[k,j]*W[j] for j in range(len(W)) if j not in pivots)) for k,p in enumerate(pivots)};eqs=[sp.expand(dec(f,W).subs(sub,simultaneous=True)) for f in d['equations']];eqs=[f for f in eqs if f];out=HERE/'weighted-unmarked7/tails'/mode;out.mkdir(parents=True,exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),variables=list(map(str,Z)),linear_equations=[enc(f,W) for f in linear],substitution={str(v):str(sub.get(v,v)) for v in W},equations=[enc(f,Z) for f in eqs],boundary_slope=str(A),boundary_intercept=str(B),scope='Fixed f=1-3u-w and determinant1 only. Boundary slope/intercept exceptional strata or raw13 support. A proper basis requires explicit coefficient point extraction, full map determinant/collision replay and support recount.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r=0,('+','.join(map(str,Z))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-tail-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),mode=mode,variables=len(Z),seconds=round(time.monotonic()-started,3));(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True);print(text[:7000],flush=True)
if __name__=='__main__':main()
