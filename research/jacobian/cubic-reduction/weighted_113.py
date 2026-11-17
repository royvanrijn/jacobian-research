#!/usr/bin/env python3
"""Unreduced first-parameter models for the remaining20-slot weight class."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
SUBS={'full':{},'qzero':{'c5':0},'univariate':{'c5':0,'c3':0},'both':{'c1':-1,'c3':-1,'c5':0},'mu-zero':{'c1':0,'c3':-1,'c5':0},'both-linear':{'c1':-1,'c2':0,'c3':-1,'c4':0,'c5':0},'mu-zero-linear':{'c1':0,'c2':0,'c3':-1,'c4':0,'c5':0}}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def dec(terms,V):return sp.Add(*(sp.Rational(c)*sp.prod(v**a for v,a in zip(V,e)) for e,c in terms))
def main():
 mode=sys.argv[1];out=HERE/'weighted-113';out.mkdir(exist_ok=True)
 if mode=='prepare':
  source=HERE/'weighted-sparse7/1_-1_-3/input.json';d=json.loads(source.read_text());names=[a['coefficient'] for a in d['atoms']];V=sp.symbols(' '.join(names),seq=True);by=dict(zip(names,V));u,w=sp.symbols('u w');fgh=[sp.Integer(1),u,w]
  for a in d['atoms']:
   b,c=a['exponent'][1:];fgh[a['row']]+=by[a['coefficient']]*u**b*w**c
  f,g,h=fgh;J=sp.expand(g*(sp.diff(f,u)*sp.diff(h,w)-sp.diff(f,w)*sp.diff(h,u))+3*h*(sp.diff(g,u)*sp.diff(f,w)-sp.diff(g,w)*sp.diff(f,u))+f*(sp.diff(g,u)*sp.diff(h,w)-sp.diff(g,w)*sp.diff(h,u))-1);eqs=[c for e,c in sp.Poly(J,u,w).terms() if c];payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),variables=names,atoms=d['atoms'],quotients=list(map(str,fgh)),equations=[enc(f,V) for f in eqs],scope='All20 degree<=7 weight(1,-1,-3) nonlinear slots, identity linear part, det1 only. No collision/sparsity restriction; proper solutions require exact extraction and replay.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print('PREPARED',len(V),len(eqs));return
 assert mode in SUBS;source=out/'input.json';d=json.loads(source.read_text());V=sp.symbols(' '.join(d['variables']),seq=True);by=dict(zip(d['variables'],V));sub={by[k]:sp.Integer(v) for k,v in SUBS[mode].items()};P=tuple(v for v in V[:5] if v not in sub);W=V[5:];eqs=[sp.expand(dec(f,V).subs(sub,simultaneous=True)) for f in d['equations']];eqs=[f for f in eqs if f];out=out/mode;out.mkdir(exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),mode=mode,parameters=list(map(str,P)),variables=list(map(str,W)),substitution={str(k):str(v) for k,v in sub.items()},equations=[enc(f,P+W) for f in eqs],scope='Changed first-parameter coefficient models. Generic units require every denominator stratum separately; proper ideals require actual coefficient extraction/determinant/collision replay. No global lower bound from a bounded solve.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r='+('(0,'+','.join(map(str,P))+')' if P else '0')+',('+','.join(map(str,W))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-first-parameter-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),mode=mode,parameters=list(map(str,P)),seconds=round(time.monotonic()-started,3));(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True);print(text[:3000],flush=True)
if __name__=='__main__':main()
