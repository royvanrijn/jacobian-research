#!/usr/bin/env python3
"""Test the at-most12-term small-third-output branch of the enlarged weighted model."""
from pathlib import Path
import hashlib,json,subprocess,time
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def expr(terms,V):return s.Add(*(s.Rational(c)*s.prod(v**a for v,a in zip(V,e)) for e,c in terms))
def encode(f,V):return [[list(e),str(c)] for e,c in s.Poly(f,*V).terms() if c]
def main():
 out=HERE/'weighted-tail7/small-h';out.mkdir(exist_ok=True);source=HERE/'weighted-tail7/input.json';d=json.loads(source.read_text());V=s.symbols(' '.join(d['variables']));names=dict(zip(d['variables'],V));linear=[expr(f,V) for f in d['boundary_equations']]
 eliminated=[a['coefficient'] for a in d['tail_atoms'] if a['row']==2 and a['quotient_exponent'] not in [[2,0],[1,1]]];linear += [names[a] for a in eliminated]
 M,b=s.linear_eq_to_matrix(linear,V);R,pivots=M.row_join(b).rref();free=[i for i in range(len(V)) if i not in pivots];W=tuple(V[i] for i in free);sub={V[p]:s.expand(R[k,len(V)]-sum(R[k,j]*V[j] for j in free)) for k,p in enumerate(pivots)}
 eqs=[];seen=set()
 for terms in d['equations']:
  f=s.expand(expr(terms,V).subs(sub,simultaneous=True))
  if not f:continue
  f=s.Poly(f,*W).monic().as_expr()
  if str(f) not in seen:seen.add(str(f));eqs.append(f)
 payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),variables=list(map(str,W)),substitution={str(v):str(sub.get(v,v)) for v in V},equations=[encode(f,W) for f in eqs],eliminated_tail_atoms=eliminated,nonlinear_bound=12,scope='Fixed first output a-3a^2b-a^3c, every weight-compatible tail atom through degree7 in output2, and output3=c+I*b^2+J*a*b*c; beta!=0 collision family retained.')
 (out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r=0,('+','.join(map(str,W))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=str(e.stdout or '')+'\n'+str(e.stderr or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),variables=len(W),equations=len(eqs),seconds=round(time.monotonic()-started,3),scope='Proper ideals require actual point extraction and exact full-map replay. No new12-term map claimed from primary basis status alone.')
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
