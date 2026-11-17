#!/usr/bin/env python3
"""Degree-seven negative-weight completions without any collision assumptions."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in s.Poly(f,*V).terms() if c]
def dec(terms,V):return s.Add(*(s.Rational(c)*s.prod(v**a for v,a in zip(V,e)) for e,c in terms))
def main():
 mode=sys.argv[1];out=HERE/'weighted-unmarked7';out.mkdir(exist_ok=True)
 if mode=='prepare':
  source=HERE/'weighted-newfirst7/input.json';d=json.loads(source.read_text());V=s.symbols(' '.join(d['substitution']));by=dict(zip(d['substitution'],V));u,w=s.symbols('u w');fgh=[1-3*u-w,u,w]
  for atom in d['atoms']:
   if 'quotient_exponent' not in atom:continue
   b,c=atom['quotient_exponent'];fgh[atom['row']]+=by[atom['coefficient']]*u**b*w**c
  f,g,h=fgh;J=s.expand(g*(s.diff(f,u)*s.diff(h,w)-s.diff(f,w)*s.diff(h,u))+2*h*(s.diff(g,u)*s.diff(f,w)-s.diff(g,w)*s.diff(f,u))+f*(s.diff(g,u)*s.diff(h,w)-s.diff(g,w)*s.diff(h,u))-1);eqs=[c for e,c in s.Poly(J,u,w).terms() if c]
  payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),variables=list(map(str,V)),atoms=d['atoms'],quotients=list(map(str,fgh)),equations=[enc(f,V) for f in eqs],scope='All23 weight-compatible degree<=7 nonlinear atoms, identity linear part, fixed first a^2b/a^3c coefficients-3,-1. Only detJF=1; no collision or sparsity equations. A completion must be extracted and fully checked; no global lower bound follows from a timeout.')
  (out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print('PREPARED',len(V),len(eqs),flush=True);return
 assert mode in ['graph','conic'];source=out/'input.json';d=json.loads(source.read_text());V=s.symbols(' '.join(d['variables']));by=dict(zip(d['variables'],V));R,P,Q,T=[by['f'+str(i)] for i in range(1,5)];sub={T:0} if mode=='graph' else {R:0};parameters=[R,P,Q] if mode=='graph' else [P,Q,T];W=tuple(v for v in V if str(v)[0]!='f');eqs=[s.expand(dec(terms,V).subs(sub,simultaneous=True)) for terms in d['equations']];eqs=[f for f in eqs if f];out=out/mode;out.mkdir(exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),parameters=list(map(str,parameters)),variables=list(map(str,W)),substitution={str(a):str(b) for a,b in sub.items()},equations=[enc(f,tuple(parameters)+W) for f in eqs],scope='Generic coefficient-field test without collision. Denominator exceptional strata require separate exact analysis.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r=(0,'+','.join(map(str,parameters))+'),('+','.join(map(str,W))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-parameter-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),mode=mode,parameters=list(map(str,parameters)),variables=len(W),equations=len(eqs),seconds=round(time.monotonic()-started,3));(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
