#!/usr/bin/env python3
"""Eliminate necessary boundary-linear equations before the weighted tail search."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def expr(terms,V):return s.Add(*(s.Rational(c)*s.prod(x**a for x,a in zip(V,e)) for e,c in terms))
def encode(f,V):return [[list(e),str(c)] for e,c in s.Poly(f,*V).terms() if c]
def main():
 mode=sys.argv[1];assert mode in ['generic','marked'];out=HERE/'weighted-tail7'/mode;out.mkdir(exist_ok=True);source=HERE/'weighted-tail7/input.json';d=json.loads(source.read_text());V=s.symbols(' '.join(d['variables']));byname=dict(zip(d['variables'],V));boundary=[expr(f,V) for f in d['boundary_equations']]
 if mode=='marked':
  u,w=s.symbols('u w');g=u+sum(byname[a['coefficient']]*u**a['quotient_exponent'][0]*w**a['quotient_exponent'][1] for a in d['tail_atoms'] if a['row']==1);G=s.Poly(s.expand(g.subs(w,1-3*u)),u)
  boundary += [G.coeff_monomial(u)-4,G.coeff_monomial(1)-3,byname['beta']-1,byname['sat']-1]
 M,b=s.linear_eq_to_matrix(boundary,V);R,pivots=M.row_join(b).rref();assert all(i<len(V) for i in pivots)
 free=[i for i in range(len(V)) if i not in pivots];W=tuple(V[i] for i in free);sub={V[pivot]:s.expand(R[k,len(V)]-sum(R[k,j]*V[j] for j in free)) for k,pivot in enumerate(pivots)}
 eqs=[];seen=set()
 for terms in d['equations']:
  f=s.expand(expr(terms,V).subs(sub,simultaneous=True))
  if not f:continue
  f=s.Poly(f,*W).monic().as_expr()
  if str(f) not in seen:seen.add(str(f));eqs.append(f)
 payload=dict(mode=mode,source_sha256=sha(Path(__file__)),input_sha256=sha(source),variables=list(map(str,W)),substitution={str(v):str(sub.get(v,v)) for v in V},equations=[encode(f,W) for f in eqs],linear_rank=len(pivots),scope='Necessary boundary equations explicitly eliminated; marked branch additionally fixes boundary slope4,constant3 and collision beta1.')
 (out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r=0,('+','.join(map(str,W))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=str(e.stdout or '')+'\n'+str(e.stderr or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,mode=mode,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),variables=len(W),equations=len(eqs),seconds=round(time.monotonic()-started,3),scope='Weighted tail coefficient calculation; a basis requires parameter extraction and exact full-map verification before a smaller construction is claimed.')
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
