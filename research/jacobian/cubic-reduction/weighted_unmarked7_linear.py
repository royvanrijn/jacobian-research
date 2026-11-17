#!/usr/bin/env python3
"""Fixed first output, det1 only, necessary polynomial-boundary degrees."""
from pathlib import Path
import hashlib,json,subprocess,time
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def main():
 source=HERE/'weighted-unmarked7/input.json';d=json.loads(source.read_text());V=sp.symbols(' '.join(d['variables']));by=dict(zip(d['variables'],V));W=tuple(v for v in V if str(v)[0]!='f');sub={by['f'+str(i)]:0 for i in range(1,5)};u,w,t=sp.symbols('u w t');g,h=[sp.sympify(f,locals=by|{'u':u,'w':w}).subs(sub).subs({u:t,w:1-3*t},simultaneous=True) for f in d['quotients'][1:]];boundary=[]
 for f,limit in [(g,1),(h,2)]:boundary += [c for (k,),c in sp.Poly(sp.expand(f),t).terms() if k>limit and c]
 eqs=[sp.expand(sum(sp.Rational(c)*sp.prod(v**a for v,a in zip(V,e)) for e,c in terms).subs(sub)) for terms in d['equations']];eqs=[f for f in eqs if f];out=HERE/'weighted-unmarked7/linear';out.mkdir(exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),parameters=[],variables=list(map(str,W)),equations=[enc(f,W) for f in eqs+boundary],boundary_equations=[enc(f,W) for f in boundary],scope='Fixed f=1-3u-w, all17 degree7 tails, detJF1 and necessary boundary degrees only. No collision assumed; a proper ideal is not a checked map.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r=0,('+','.join(map(str,W))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs+boundary)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=__import__('subprocess').run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-boundary-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),seconds=round(time.monotonic()-started,3));(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True);print(text[:5000],flush=True)
if __name__=='__main__':main()
