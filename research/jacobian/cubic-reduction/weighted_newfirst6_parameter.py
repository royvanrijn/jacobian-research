#!/usr/bin/env python3
"""Treat the new first-output coefficient as a parameter field to expose exceptional slices."""
from pathlib import Path
import hashlib,json,subprocess,time
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 out=HERE/'weighted-newfirst6/parameter-field';out.mkdir(exist_ok=True);source=HERE/'weighted-newfirst6/input.json';d=json.loads(source.read_text());V=s.symbols(' '.join(d['variables']));theta=s.Symbol('theta');old={str(v):v for v in V};assert 'f2' in old;W=tuple(v for v in V if str(v)!='f2')
 eqs=[s.expand(sum(s.Rational(c)*s.prod(v**a for v,a in zip(V,e)) for e,c in terms).subs(old['f2'],theta)) for terms in d['equations']]
 script=out/'model.sing';script.write_text('\n'.join(['ring r=(0,theta),('+','.join(map(str,W))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=r.stdout+r.stderr;status='primary-parameter-basis' if r.returncode==0 and text.startswith('DIM:') and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=str(e.stdout or '')+'\n'+str(e.stderr or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(source),seconds=round(time.monotonic()-started,3),scope='Generic coefficient field Q(theta). A unit ideal only excludes generic theta; a lifted certificate and its exceptional denominator roots are required for global parameter coverage.')
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
