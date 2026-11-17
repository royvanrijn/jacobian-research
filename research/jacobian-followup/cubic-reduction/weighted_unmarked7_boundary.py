#!/usr/bin/env python3
"""Necessary boundary degrees for unmarked weighted Keller completions."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def dec(terms,V):return sp.Add(*(sp.Rational(c)*sp.prod(v**a for v,a in zip(V,e)) for e,c in terms))
def main():
 mode=sys.argv[1];assert mode in ['cubicgraph','quadraticgraph','parabola','parallel'];source=HERE/'weighted-unmarked7/input.json';d=json.loads(source.read_text());V=sp.symbols(' '.join(d['variables']));by=dict(zip(d['variables'],V));R,P,Q,T=[by['f'+str(i)] for i in range(1,5)];W=tuple(v for v in V if str(v)[0]!='f');t,l,kappa=sp.symbols('t l kappa');u,w=sp.symbols('u w');f,g,h=[sp.sympify(f,locals=by|{'u':u,'w':w}) for f in d['quotients']]
 if mode=='cubicgraph':sub={Q:0,T:0};parameters=[R,P];U=t;B=1-3*t+P*t**2+R*t**3
 elif mode=='quadraticgraph':sub={R:0,Q:0,T:0};parameters=[P];U=t;B=1-3*t+P*t**2
 elif mode=='parabola':sub={R:0,P:(kappa+3)**2*T,Q:2*(kappa+3)*T};parameters=[kappa,T];U=(-1+t-T*t**2)/kappa;B=t-(kappa+3)*U
 else:sub={R:0,P:9*T,Q:6*T};parameters=[T];U=t;B=l-3*t
 gh=[sp.cancel(z.subs(sub,simultaneous=True).subs({u:U,w:B},simultaneous=True)) for z in [g,h]];boundary=[]
 for z,limit in zip(gh,[1,2]):
  num=sp.fraction(z)[0]
  for (k,),c in sp.Poly(num,t).terms():
   if k<=limit or not c:continue
   if mode=='parallel':
    remainder=sp.rem(c,1-l+T*l**2,l);boundary += [sp.fraction(sp.cancel(c))[0] for e,c in sp.Poly(remainder,l).terms() if c]
   else:boundary.append(c)
 eqs=[sp.expand(dec(z,V).subs(sub,simultaneous=True)) for z in d['equations']];eqs=[z for z in eqs if z];out=HERE/'weighted-unmarked7'/mode;out.mkdir(exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),parameters=list(map(str,parameters)),variables=list(map(str,W)),substitution={str(a):str(b) for a,b in sub.items()},boundary_coordinates=[str(U),str(B)],boundary_equations=[enc(z,tuple(parameters)+W) for z in boundary],equations=[enc(z,tuple(parameters)+W) for z in eqs+boundary],scope='Necessary boundary degree constraints; cubic/quadraticgraph have Wr=1, parabola has Wr=1/kappa for kappa!=0. Parallel lines require distinct roots of1-l+T*l^2, hence T!=1/4; degree conditions are imposed modulo that quadratic. Parameter zero/denominator exceptional strata are not discarded.')
 (out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r=(0,'+','.join(map(str,parameters))+'),('+','.join(map(str,W))+'),dp;','short=0;','ideal I='+','.join(sing(z) for z in eqs+boundary)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-boundary-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),mode=mode,boundary_equations=len(boundary),seconds=round(time.monotonic()-started,3));(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True);print(text[:1000],flush=True)
if __name__=='__main__':main()
