#!/usr/bin/env python3
"""Changed first-coefficient models for the open16-slot mixed-weight search."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as sp
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
FIXED={'u2':(1,0,0,0),'uv':(0,1,0,0),'v2':(0,0,0,1),'diagonal':(1,0,0,1),'rank1':(1,2,0,1),'mixed':(1,1,0,1),'u2uv':(1,1,0,0),'uvv2':(0,1,0,1)}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(f,V):return [[list(e),str(c)] for e,c in sp.Poly(f,*V).terms() if c]
def main():
 mode=sys.argv[1];assert mode in ['full','conic','linearv']+list(FIXED);source=HERE/'weighted-sparse7/2_-1_-3/input.json';d=json.loads(source.read_text());V=sp.symbols(' '.join(d['variables']),seq=True);by=dict(zip(d['variables'],V));first=[by['c'+str(i)] for i in range(1,5)];W=tuple(v for v in V if v not in first);sub={}
 if mode=='conic':sub[first[2]]=0
 elif mode=='linearv':sub[first[3]]=0
 elif mode in FIXED:sub=dict(zip(first,map(sp.Integer,FIXED[mode])))
 P=tuple(v for v in first if v not in sub);eqs=[]
 for terms in d['equations']:
  f=sp.expand(sum(sp.Rational(c)*sp.prod(v**a for v,a in zip(V,e)) for e,c in terms).subs(sub,simultaneous=True))
  if f:eqs.append(f)
 out=HERE/'weighted-213'/mode;out.mkdir(parents=True,exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),mode=mode,parameters=list(map(str,P)),variables=list(map(str,W)),substitution={str(a):str(b) for a,b in sub.items()},equations=[enc(f,P+W) for f in eqs],scope='Changed coefficient models. A generic field unit excludes no denominator roots without further proof. Fixed cases are construction attempts, not exhaustive parameter coverage. Proper ideals need actual coefficient extraction and independent determinant/support/collision replay. No12-term or global claim from basis status.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r='+('(0,'+','.join(map(str,P))+')' if P else '0')+',('+','.join(map(str,W))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-changed-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),mode=mode,parameters=list(map(str,P)),seconds=round(time.monotonic()-started,3));(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True);print(text[:4000],flush=True)
if __name__=='__main__':main()
