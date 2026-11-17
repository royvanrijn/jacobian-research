#!/usr/bin/env python3
"""Normalized first-output slices for the13-slot weight(1,-2,-3) model."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
CASES={'both':(-1,-1),'mu-zero':(0,-1),'nu-zero':(-1,0),'both-zero':(0,0)}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def main():
 case,mode=sys.argv[1:];assert case in CASES and mode in ['zero','field'];source=HERE/'weighted-sparse7/1_-2_-3/input.json';d=json.loads(source.read_text());V=sp.symbols(' '.join(d['variables']));by=dict(zip(d['variables'],V));mu,nu=CASES[case];sub={by['c1']:mu,by['c2']:nu};parameters=[by['c3']]
 if mode=='zero':sub[by['c3']]=0;parameters=[]
 W=tuple(v for v in V if v not in [by['c1'],by['c2'],by['c3']]);eqs=[]
 for terms in d['equations']:
  f=sp.expand(sum(sp.Rational(c)*sp.prod(v**a for v,a in zip(V,e)) for e,c in terms).subs(sub,simultaneous=True))
  if f:eqs.append(f)
 out=HERE/'weighted-sparse7/1_-2_-3/slices'/case/mode;out.mkdir(parents=True,exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),case=case,mode=mode,parameters=list(map(str,parameters)),variables=list(map(str,W)),substitution={str(k):str(v) for k,v in sub.items()},equations=[enc(f,tuple(parameters)+W) for f in eqs],scope='All ten tail coefficients retained. Diagonal normalization sets the specified nonzero first coefficients to-1. Zero mode eliminates the third first atom, so any completion has<=12 terms and a nonlinear first component except the both-zero case. Field mode is generic only; all denominator exceptions remain to be resolved. Proper ideals require extraction and full exact replay.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r='+('(0,'+','.join(map(str,parameters))+')' if parameters else '0')+',('+','.join(map(str,W))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-slice-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),case=case,mode=mode,seconds=round(time.monotonic()-started,3));(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True);print(text[:5000],flush=True)
if __name__=='__main__':main()
