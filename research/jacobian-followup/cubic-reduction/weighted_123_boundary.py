#!/usr/bin/env python3
"""Necessary low-degree (2,3)-Wronskian boundary constraints for normalized f."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def main():
 case=sys.argv[1];assert case in ['both','mu-zero'];source=HERE/'weighted-sparse7/1_-2_-3/slices'/case/'field/input.json';d=json.loads(source.read_text());P=sp.Symbol('c3');W=sp.symbols(' '.join(d['variables']),seq=True);by=dict(zip(map(str,W),W));u,w,t=sp.symbols('u w t');full=json.loads((HERE/'weighted-sparse7/1_-2_-3/input.json').read_text());gh=[u,w]
 for a in full['atoms']:
  if a['row']>0:
   row=a['row'];b,c=a['exponent'][1:];gh[row-1]+=by[a['coefficient']]*u**b*w**c
 U=t;B=1-t+P*t**2 if case=='both' else 1+P*t**2;boundary=[]
 for f,limit in zip(gh,[2,3]):boundary += [c for (k,),c in sp.Poly(sp.expand(f.subs({u:U,w:B},simultaneous=True)),t).terms() if k>limit and c]
 eqs=[sp.Add(*(sp.Rational(c)*sp.prod(v**a for v,a in zip((P,)+W,e)) for e,c in terms)) for terms in d['equations']];out=source.parent.parent/'boundary';out.mkdir(exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),parameters=['c3'],variables=list(map(str,W)),boundary_equations=[enc(f,(P,)+W) for f in boundary],equations=[enc(f,(P,)+W) for f in eqs+boundary],scope='On f=0,2G*Hprime-3Gprime*H=1 with degreesG,H<=4. Leading degrees force degG<=2,degH<=3. Thus these three boundary equations are necessary for every coefficient solution; no generic exception is removed.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r=(0,c3),('+','.join(map(str,W))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs+boundary)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-boundary-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),case=case,seconds=round(time.monotonic()-started,3));(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True);print(text[:2000],flush=True)
if __name__=='__main__':main()
