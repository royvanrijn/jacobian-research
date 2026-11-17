#!/usr/bin/env python3
"""Explore the remaining nonzero g(w^2) tail branch after structural reductions."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def expr(terms,V):return s.Add(*(s.Rational(c)*s.prod(x**a for x,a in zip(V,e)) for e,c in terms))
def encode(f,V):return [[list(e),str(c)] for e,c in s.Poly(f,*V).terms() if c]
def main():
 mode='generic';out=HERE/'weighted-tail7/mid-w';out.mkdir(parents=True,exist_ok=True);source=HERE/'weighted-tail7/input.json';d=json.loads(source.read_text());V=s.symbols(' '.join(d['variables']));names=dict(zip(d['variables'],V));lookup={(a['row'],tuple(a['quotient_exponent'])):names[a['coefficient']] for a in d['tail_atoms']}
 linear=[expr(f,V) for f in d['boundary_equations']]+[lookup[key] for key in [(1,(1,2)),(2,(0,3)),(2,(1,2)),(2,(4,0)),(2,(3,1))]]
 K=lookup[1,(0,2)];T=lookup[2,(0,3)]
 if mode=='unit-slice':linear += [K-1,T-1]
 M,b=s.linear_eq_to_matrix(linear,V);R,pivots=M.row_join(b).rref();free=[i for i in range(len(V)) if i not in pivots];W=tuple(V[i] for i in free);sub={V[p]:s.expand(R[k,len(V)]-sum(R[k,j]*V[j] for j in free)) for k,p in enumerate(pivots)}
 eqs=[];seen=set()
 for f in [expr(terms,V) for terms in d['equations']]:
  f=s.expand(f.subs(sub,simultaneous=True))
  if not f:continue
  f=s.Poly(f,*W).monic().as_expr()
  if str(f) not in seen:seen.add(str(f));eqs.append(f)
 if mode=='generic':
  t=s.Symbol('extra_inv');W=W+(t,);eqs.append(s.expand(1-t*K.subs(sub)))
 payload=dict(mode=mode,source_sha256=sha(Path(__file__)),input_sha256=sha(source),variables=list(map(str,W)),substitution={str(v):str(sub.get(v,v)) for v in V},equations=[encode(f,W) for f in eqs],scope='Nonzero g(w^2), with g(u*w^2),h(w^3),h(u*w^2),h(u^4),h(u^3*w) zero as forced by the high-v equations.')
 (out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r=0,('+','.join(map(str,W))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=str(e.stdout or '')+'\n'+str(e.stderr or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),mode=mode,variables=len(W),equations=len(eqs),seconds=round(time.monotonic()-started,3),scope='Exploratory high-w branch. Proper basis requires actual coefficient extraction and full-map determinant/collision replay.')
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
