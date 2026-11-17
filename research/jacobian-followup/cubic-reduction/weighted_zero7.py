#!/usr/bin/env python3
"""Complete degree-seven weighted determinant models on zero first coefficients."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
CASES={'mu-zero':(0,-1),'nu-zero':(-3,0),'both-zero':(0,0)}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def dec(terms,V):return sp.Add(*(sp.Rational(c)*sp.prod(v**a for v,a in zip(V,e)) for e,c in terms))
def main():
 case,mode=sys.argv[1:];assert case in CASES;out=HERE/'weighted-zero7'/case;out.mkdir(parents=True,exist_ok=True)
 if mode=='prepare':
  source=HERE/'weighted-unmarked7/input.json';d=json.loads(source.read_text());V=sp.symbols(' '.join(d['variables']));by=dict(zip(d['variables'],V));u,w=sp.symbols('u w');mu,nu=CASES[case];fgh=[1+mu*u+nu*w,u,w];atoms=[]
  for atom in d['atoms']:
   atom=dict(atom)
   if 'quotient_exponent' not in atom:atom['coefficient']=str(mu if atom['exponent']==[2,1,0] else nu)
   else:
    b,c=atom['quotient_exponent'];fgh[atom['row']]+=by[atom['coefficient']]*u**b*w**c
   atoms.append(atom)
  f,g,h=fgh;J=sp.expand(g*(sp.diff(f,u)*sp.diff(h,w)-sp.diff(f,w)*sp.diff(h,u))+2*h*(sp.diff(g,u)*sp.diff(f,w)-sp.diff(g,w)*sp.diff(f,u))+f*(sp.diff(g,u)*sp.diff(h,w)-sp.diff(g,w)*sp.diff(h,u))-1);eqs=[c for e,c in sp.Poly(J,u,w).terms() if c];payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),case=case,fixed_first_coefficients=[mu,nu],variables=list(map(str,V)),atoms=atoms,quotients=list(map(str,fgh)),equations=[enc(f,V) for f in eqs],scope='All weight-compatible degree<=7 atoms with identity linear part and declared zero/nonzero first coefficients; only determinant1, no collision or sparsity condition. Proper ideals require full coefficient extraction and determinant/collision replay before a counterexample claim.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print(case,'PREPARED',len(V),len(eqs),flush=True);return
 assert mode in ['graph','conic'];source=out/'input.json';d=json.loads(source.read_text());V=sp.symbols(' '.join(d['variables']));by=dict(zip(d['variables'],V));R,P,Q,T=[by['f'+str(i)] for i in range(1,5)];sub={T:0} if mode=='graph' else {R:0};parameters=[R,P,Q] if mode=='graph' else [P,Q,T];W=tuple(v for v in V if str(v)[0]!='f');eqs=[sp.expand(dec(terms,V).subs(sub,simultaneous=True)) for terms in d['equations']];eqs=[f for f in eqs if f];out=out/mode;out.mkdir(exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),case=case,parameters=list(map(str,parameters)),variables=list(map(str,W)),substitution={str(a):str(b) for a,b in sub.items()},equations=[enc(f,tuple(parameters)+W) for f in eqs],scope='Generic parameter-field test only. Every denominator exceptional branch remains to be checked. No collision assumptions.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r=(0,'+','.join(map(str,parameters))+'),('+','.join(map(str,W))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-parameter-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),case=case,mode=mode,seconds=round(time.monotonic()-started,3));(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
