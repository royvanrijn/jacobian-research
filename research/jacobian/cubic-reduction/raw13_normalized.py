#!/usr/bin/env python3
"""Torus-normalized raw13 coefficient model, with collision parameter beta retained."""
from pathlib import Path
import hashlib,json,subprocess,time
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 out=HERE/'raw13-normalized';out.mkdir(exist_ok=True);source=HERE/'raw13-specialization/input.json';data=json.loads(source.read_text());X=s.symbols('x y z');A=s.symbols('a1:14');beta,t=s.symbols('beta sat');values=list(A);values[0]=-1;values[1]=-3
 assert data['atoms'][0]['row']==0 and data['atoms'][0]['exponent']==[3,0,1]
 assert data['atoms'][1]['row']==0 and data['atoms'][1]['exponent']==[2,1,0]
 F=list(X)
 for a,atom in zip(values,data['atoms']):F[atom['row']]+=a*s.prod(x**e for x,e in zip(X,atom['exponent']))
 coefficients=[c for e,c in s.Poly(s.expand(s.det(s.Matrix(F).jacobian(X))-1),*X).terms() if c]
 point=dict(zip(X,[-1,beta,1+3*beta]));assert s.expand(F[0].subs(point))==0
 eqs=coefficients+[s.expand(F[1].subs(point)-beta),1-t*beta];V=tuple(A[2:])+(beta,t)
 payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),variables=list(map(str,V)),fixed_coefficients={'a1':'-1','a2':'-3'},equations=[[[list(e),str(c)] for e,c in s.Poly(f,*V).terms() if c] for f in eqs],scope='Normalize the two nonzero first-output coefficients by diagonal source/target conjugation, retaining the transformed common b-coordinate beta!=0. Only this13-atom pattern and collision family are tested.')
 (out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
 script=out/'normalized.sing';script.write_text('\n'.join(['ring r=0,('+','.join(map(str,V))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in eqs)+';','ideal G=std(I);','print("DIM:"+string(dim(G)));','print(G);','quit;'])+'\n');started=time.monotonic()
 try:
  r=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=10);text=r.stdout+r.stderr;status='primary-basis' if r.returncode==0 and '?' not in text and text.startswith('DIM:') else 'error'
 except subprocess.TimeoutExpired as exc:text=str(exc.stdout or '')+'\n'+str(exc.stderr or '');status='timeout'
 (out/'normalized.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),seconds=round(time.monotonic()-started,3),scope='Primary normalized coefficient ideal only; no raw12-term exclusion or smaller witness without extraction and independent replay.')
 (out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
