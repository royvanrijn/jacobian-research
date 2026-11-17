#!/usr/bin/env python3
"""Necessary polynomial-boundary Wronskian constraints for two residual lines."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in s.Poly(f,*V).terms() if c]
def main():
 mode=sys.argv[1];assert mode in ['cubicgraph','parabolic'];out=HERE/'weighted-newfirst7/parameter-field'/(mode+'-boundary');out.mkdir(exist_ok=True);source=HERE/'weighted-newfirst7/input.json';d=json.loads(source.read_text());V=s.symbols(' '.join(d['variables']));by=dict(zip(d['variables'],V));P,Q,T=[by[n] for n in ['f2','f3','f4']];t=s.Symbol('t');u,w=s.symbols('u w');W=tuple(v for v in V if v not in [P,Q,T]);sub={Q:0,T:0} if mode=='cubicgraph' else {P:16*T,Q:8*T};parameter=P if mode=='cubicgraph' else T
 U=t if mode=='cubicgraph' else -1+t-T*t**2;B=1-3*t+P*t**2+P*t**3 if mode=='cubicgraph' else t-4*U
 gh=[]
 for row in [1,2]:
  expr=(u if row==1 else w)+sum(s.sympify(d['substitution'].get(a['coefficient'],a['coefficient']),locals=by)*u**a['quotient_exponent'][0]*w**a['quotient_exponent'][1] for a in d['atoms'] if a['row']==row)
  gh.append(s.expand(expr.subs(sub,simultaneous=True).subs({u:U,w:B},simultaneous=True)))
 boundary=[]
 for f,limit in zip(gh,[1,2]):
  boundary += [c for (k,),c in s.Poly(f,t).terms() if k>limit and c]
 eqs=[s.expand(sum(s.Rational(c)*s.prod(v**a for v,a in zip(V,e)) for e,c in terms).subs(sub,simultaneous=True)) for terms in d['equations']];eqs=[f for f in eqs if f]
 M,b=s.linear_eq_to_matrix(boundary,W);aug=M.row_join(b);R,pivots=aug.rref();assert all(i<len(W) for i in pivots);free=[i for i in range(len(W)) if i not in pivots];Z=tuple(W[i] for i in free);solve={W[p]:s.factor(R[k,len(W)]-sum(R[k,j]*W[j] for j in free)) for k,p in enumerate(pivots)};reduced=[s.factor(f.subs(solve,simultaneous=True)) for f in eqs];reduced=[f for f in reduced if f]
 payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),mode=mode,parameters=[str(parameter)],variables=list(map(str,W)),equations=[enc(f,(parameter,)+W) for f in eqs+boundary],boundary_equations=[enc(f,(parameter,)+W) for f in boundary],boundary_g=str(gh[0]),boundary_h=str(gh[1]),boundary_coordinates=[str(U),str(B)],linear_rank=len(pivots),free_variables=list(map(str,Z)),solution={str(k):str(v) for k,v in solve.items()},scope='Boundary Wronskian lemma: g h_prime-2g_prime h=1 forces deg(g)<=1,deg(h)<=2. Additional linear constraints are necessary; they are not asserted to be ideal memberships without that lemma.')
 (out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r=(0,'+str(parameter)+'),('+','.join(map(str,Z))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in reduced)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-boundary-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),linear_rank=len(pivots),free_variables=list(map(str,Z)),seconds=round(time.monotonic()-started,3));(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True);print(text[:5000],flush=True)
if __name__=='__main__':main()
