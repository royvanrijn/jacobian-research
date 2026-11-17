#!/usr/bin/env python3
"""Complete base-monomial carrier extensions, composed through old carrier inverse."""
from pathlib import Path
import hashlib,json,subprocess,time
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 source=HERE/'cubic12-integer-witness.json';d=json.loads(source.read_text());X=s.symbols('x1:13');F=[sum(s.Rational(t['coefficient'])*s.prod(x**k for x,k in zip(X,t['exponent'])) for t in row) for row in d['coefficients']];atoms=[(i,e,c) for i in range(3,12) for e,c in s.Poly(F[i]-X[i],*X).terms() if c]
 weights=[1,-1,-2,0,-3,2,-1,-2,-2,-1,-1,0];seen={(i,e) for i,e,c in atoms}
 for i in range(3,12):
  for a in range(4):
   for b in range(4):
    for c in range(4):
     e=(a,b,c)+(0,)*9
     if 2<=a+b+c<=3 and a-b-2*c==weights[i] and (i,e) not in seen:atoms.append((i,e,s.Integer(1)));seen.add((i,e))
 assert len(atoms)==22;params=s.symbols('z1:23');new=list(X);labels=[]
 for z,(i,e,c) in zip(params,atoms):new[i]+=z*c*s.prod(x**k for x,k in zip(X,e));labels.append(dict(parameter=str(z),row=i,exponent=list(e),base_coefficient=str(c)))
 E=list(X)
 for i in range(3,12):E[i]=s.expand(new[i]-(F[i]-X[i]).subs(dict(zip(X,E)),simultaneous=True));assert not E[i].has(X[i]) or s.diff(E[i],X[i])==1
 G=[s.expand(f.subs(dict(zip(X,E)),simultaneous=True)) for f in F];assert G[3:]==[s.expand(f) for f in new[3:]];high=[];low=[]
 for row,f in enumerate(G):
  for e,c in s.Poly(f,*X).terms():
   if sum(e)>3:high.append(s.factor(c))
   elif sum(e)>=2:low.append(dict(row=row,exponent=list(e),coefficient=str(s.factor(c))))
 unique=list(dict.fromkeys(high));out=HERE/'carrier-reparam22';out.mkdir(exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),variables=list(map(str,params)),atoms=labels,source_substitution=list(map(str,E)),carrier_outputs=list(map(str,new[3:])),high_equations=[[[list(e),str(c)] for e,c in s.Poly(f,*params).terms() if c] for f in unique],low_coefficients=low,scope='All original carrier atoms plus every base(a,b,c)-only weight-compatible monomial of degrees2,3 varied independently; source automorphism E=C_old^-1 C_new fixes base variables and is triangular with identity diagonal. G=F E is Keller with transported collision; degree3 requires every high equation. Only actual extracted/replayed candidate improves bounds.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print('PARAMETERS',len(params),'HIGH',len(unique),'LOW',len(low),flush=True);script=out/'model.sing';script.write_text('\n'.join(['ring r=0,('+','.join(map(str,params))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in unique)+';','ideal B=std(I);','print("DIM:"+string(dim(B)));','print(B);','quit;'])+'\n');started=time.monotonic()
 try:
  p=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=p.stdout+p.stderr;status='primary-basis' if p.returncode==0 and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (out/'model.log').write_text(text);(out/'checkpoint.json').write_text(json.dumps(dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),seconds=round(time.monotonic()-started,3)),indent=2,sort_keys=True)+'\n');print(status,text,flush=True)
if __name__=='__main__':main()
