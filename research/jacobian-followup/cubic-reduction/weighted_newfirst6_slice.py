#!/usr/bin/env python3
"""Exact rational slices of the new first-output coefficient in the degree-six search."""
from pathlib import Path
import hashlib,json,subprocess,time,sys
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 theta=s.Rational(sys.argv[1]);tag=str(theta).replace('/','p');out=HERE/'weighted-newfirst6'/('theta-'+tag);out.mkdir(exist_ok=True);source=HERE/'weighted-newfirst6/input.json';d=json.loads(source.read_text());V=s.symbols(' '.join(d['variables']));byname=dict(zip(d['variables'],V));assert 'f2' in byname
 W=tuple(v for v in V if str(v)!='f2');eqs=[]
 for p in d['equations']:
  f=sum(s.Rational(c)*s.prod(v**a for v,a in zip(V,e)) for e,c in p);f=s.expand(f.subs(byname['f2'],theta))
  if f:eqs.append(f)
 script=out/'model.sing';script.write_text('\n'.join(['ring r=0,('+','.join(map(str,W))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=5);text=r.stdout+r.stderr;status='primary-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=str(e.stdout or '')+'\n'+str(e.stderr or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(source),theta=str(theta),variables=len(W),seconds=round(time.monotonic()-started,3),scope='One rational new-first-output slice in the degree-six model. A proper ideal needs coefficient extraction and independent collision replay; no global degree-six exclusion.')
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
