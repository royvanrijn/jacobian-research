#!/usr/bin/env python3
"""Polynomial-boundary constraints on the remaining zero-coefficient branches."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def dec(terms,V):return sp.Add(*(sp.Rational(c)*sp.prod(v**a for v,a in zip(V,e)) for e,c in terms))
def main():
 case,mode=sys.argv[1:];assert case in ['mu-zero','nu-zero'];assert mode in ['cubicgraph','quadraticgraph','linear','parabola'];assert case=='mu-zero' or mode in ['linear','parabola'];source=HERE/'weighted-zero7'/case/'input.json';d=json.loads(source.read_text());V=sp.symbols(' '.join(d['variables']));by=dict(zip(d['variables'],V));R,P,Q,T=[by['f'+str(i)] for i in range(1,5)];W=tuple(v for v in V if str(v)[0]!='f');u,w,t,sigma=sp.symbols('u w t sigma');mu,nu=d['fixed_first_coefficients'];g,h=[sp.sympify(f,locals=by|{'u':u,'w':w}) for f in d['quotients'][1:]]
 if mode=='cubicgraph':sub={Q:0,T:0};parameters=[R,P];U=t;B=1+P*t**2+R*t**3;required=sp.Integer(1)
 elif mode=='quadraticgraph':sub={R:0,Q:0,T:0};parameters=[P];U=t;B=1+P*t**2;required=sp.Integer(1)
 elif mode=='linear':
  sub={R:0,P:0,Q:0,T:0};parameters=[];U=t if case=='mu-zero' else sp.Rational(1,3);B=sp.Integer(1) if case=='mu-zero' else t;required=sp.Integer(1)
 else:
  sub={R:0,P:sigma**2*T,Q:2*sigma*T};parameters=[sigma,T];kappa=mu-nu*sigma;U=(-1-nu*t-T*t**2)/kappa;B=t-sigma*U;required=kappa*T
 gh=[sp.cancel(f.subs(sub,simultaneous=True).subs({u:U,w:B},simultaneous=True)) for f in [g,h]];boundary=[]
 for z,limit in zip(gh,[1,2]):boundary += [c for (k,),c in sp.Poly(sp.fraction(z)[0],t).terms() if k>limit and c]
 eqs=[sp.expand(dec(z,V).subs(sub,simultaneous=True)) for z in d['equations']];eqs=[z for z in eqs if z];out=HERE/'weighted-zero7'/case/mode;out.mkdir(exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),parameters=list(map(str,parameters)),variables=list(map(str,W)),substitution={str(a):str(b) for a,b in sub.items()},boundary_coordinates=[str(U),str(B)],boundary_equations=[enc(z,tuple(parameters)+W) for z in boundary],equations=[enc(z,tuple(parameters)+W) for z in eqs+boundary],required_nonzero=str(required),scope='Necessary polynomial-boundary Wronskian degree constraints. Parameter and coordinate denominators require exact nonzero conditions and exceptional analysis; no collision assumption.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r='+('(0,'+','.join(map(str,parameters))+')' if parameters else '0')+',('+','.join(map(str,W))+'),dp;','short=0;','ideal I='+','.join(sing(z) for z in eqs+boundary)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-boundary-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),case=case,mode=mode,boundary_equations=len(boundary),seconds=round(time.monotonic()-started,3));(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True);print(text[:1000],flush=True)
if __name__=='__main__':main()
