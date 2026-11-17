#!/usr/bin/env python3
"""Parameter-field completion tests for new degree-seven first outputs."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 mode=sys.argv[1];assert mode in ['full','conic','graph','cubicgraph'];out=HERE/'weighted-newfirst7/parameter-field'/mode;out.mkdir(parents=True,exist_ok=True);source=HERE/'weighted-newfirst7/input.json';d=json.loads(source.read_text());V=s.symbols(' '.join(d['variables']));byname=dict(zip(d['variables'],V));P,Q,T=[byname[n] for n in ['f2','f3','f4']]
 sub={};parameters=[P,Q,T]
 if mode=='conic':sub[P]=4*Q-16*T;parameters=[Q,T]
 if mode=='graph':sub[T]=0;parameters=[P,Q]
 if mode=='cubicgraph':sub[Q]=sub[T]=0;parameters=[P]
 W=tuple(v for v in V if str(v) not in ['f2','f3','f4']);eqs=[]
 for terms in d['equations']:
  f=s.expand(sum(s.Rational(c)*s.prod(v**a for v,a in zip(V,e)) for e,c in terms).subs(sub,simultaneous=True))
  if f:eqs.append(f)
 inputdata=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),mode=mode,parameters=list(map(str,parameters)),variables=list(map(str,W)),substitution={str(a):str(b) for a,b in sub.items()},equations=[[[list(e),str(c)] for e,c in s.Poly(f,*(tuple(parameters)+W)).terms() if c] for f in eqs],scope='Generic parameter-field test only; any exceptional denominator hypersurfaces require separate analysis.')
 (out/'input.json').write_text(json.dumps(inputdata,indent=2,sort_keys=True)+'\n');script=out/'model.sing';script.write_text('\n'.join(['ring r=(0,'+','.join(map(str,parameters))+'),('+','.join(map(str,W))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-parameter-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=str(e.stdout or '')+'\n'+str(e.stderr or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),mode=mode,parameters=list(map(str,parameters)),variables=len(W),equations=len(eqs),seconds=round(time.monotonic()-started,3),scope='A generic unit is not a global exclusion. A proper ideal requires actual coefficient extraction and full-map replay. No smaller map is claimed from basis status.')
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
