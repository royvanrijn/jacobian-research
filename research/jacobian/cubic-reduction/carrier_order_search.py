#!/usr/bin/env python3
"""Complete weighted degree2/3 carrier models in changed triangular orders.
Each is Keller by E=C_old^-1 C_new. Seek an identity carrier, which would
permit a cubic dimension deletion with the inherited collision.
"""
from pathlib import Path
from itertools import combinations_with_replacement
import hashlib,json,subprocess,time,sys
import sympy as s
from power_single_compress import sing
HERE=Path(__file__).resolve().parent
W=[1,-1,-2,0,-3,2,-1,-2,-2,-1,-1,0]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(order):
 source=HERE/'cubic12-integer-witness.json';d=json.loads(source.read_text());X=s.symbols('x1:13');F=[sum(s.Rational(t['coefficient'])*s.prod(x**k for x,k in zip(X,t['exponent'])) for t in row) for row in d['coefficients']];position={j:k for k,j in enumerate(order)};allowed={i:sorted([0,1,2]+order[:position[i]]) for i in order};original={(i,e):c for i in range(3,12) for e,c in s.Poly(F[i]-X[i],*X).terms() if c};atoms=[];seen=set()
 for (i,e),c in original.items():
  if all(not k or j in allowed[i] for j,k in enumerate(e)):atoms.append((i,e,c));seen.add((i,e))
 for i in range(3,12):
  for degree in [2,3]:
   for js in combinations_with_replacement(allowed[i],degree):
    e=tuple(js.count(a) for a in range(12))
    if sum(W[j] for j in js)==W[i] and (i,e) not in seen:atoms.append((i,e,s.Integer(1)));seen.add((i,e))
 z=s.symbols('z1:'+str(len(atoms)+1));new=list(X);labels=[]
 for zz,(i,e,c) in zip(z,atoms):new[i]+=zz*c*s.prod(x**k for x,k in zip(X,e));labels.append(dict(parameter=str(zz),row=i,exponent=list(e),base_coefficient=str(c)))
 E=list(X)
 for i in range(3,12):E[i]=s.expand(new[i]-(F[i]-X[i]).subs(dict(zip(X,E)),simultaneous=True))
 G=[s.expand(f.subs(dict(zip(X,E)),simultaneous=True)) for f in F];assert G[3:]==[s.expand(f) for f in new[3:]];high=[];low=[]
 for i,f in enumerate(G):
  for e,c in s.Poly(f,*X).terms():
   if sum(e)>3:high.append(s.factor(c))
   elif sum(e)>=2:low.append(dict(row=i,exponent=list(e),coefficient=str(s.factor(c))))
 high=list(dict.fromkeys(high));slug='_'.join(str(j+1) for j in order);out=HERE/'carrier-order-search'/slug;out.mkdir(parents=True,exist_ok=True);payload=dict(source_sha256=sha(Path(__file__)),input_sha256=sha(source),order=[i+1 for i in order],weights=W,variables=list(map(str,z)),atoms=labels,source_substitution=list(map(str,E)),carrier_outputs=list(map(str,new[3:])),high_equations=[[[list(e),str(c)] for e,c in s.Poly(f,*z).terms() if c] for f in high],low_coefficients=low,scope='Complete weighted degree2/3 carrier monomials triangular in the stated order, linear partI. E=C_old^-1 C_new gives G=K E. An identity carrier row would give an eleven-variable cubic counterexample after restriction, requiring independent extraction/determinant/collision replay.');(out/'input.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');lines=['ring r=0,('+','.join(map(str,z))+'),dp;','short=0;','ideal I='+','.join(sing(f) for f in high)+';','ideal B=std(I);','print("MODEL_DIM:"+string(dim(B)));'];branches=[]
 for i in range(3,12):
  coeff=[zz for zz,a in zip(z,labels) if a['row']==i];branches.append(dict(row=i+1,zero_parameters=list(map(str,coeff))));lines+=['ideal J'+str(i)+'=B,'+','.join(map(str,coeff))+';','ideal T'+str(i)+'=std(J'+str(i)+');','print("CARRIER:'+str(i+1)+':DIM:"+string(dim(T'+str(i)+')));','if(T'+str(i)+'[1]!=1){print("PROPER:'+str(i+1)+'");print(T'+str(i)+');}']
 lines+=['quit;'];script=out/'model.sing';script.write_text('\n'.join(lines)+'\n');started=time.monotonic();print('PREPARED',slug,'ATOMS',len(z),'HIGH',len(high),flush=True)
 try:
  result=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20);text=result.stdout+result.stderr;status='primary-bases-complete' if result.returncode==0 and '?' not in text else 'error'
 except subprocess.TimeoutExpired as e:text=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (out/'model.log').write_text(text);receipt=dict(status=status,source_sha256=sha(Path(__file__)),input_sha256=sha(out/'input.json'),order=[j+1 for j in order],identity_branches=branches,seconds=round(time.monotonic()-started,3),scope='Primary coefficient bases only; dimension/proper ideals are not extracted Keller witnesses or independent nonexistence certificates. Timeouts remain unresolved.');(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(slug,status,text[:5000],flush=True)
def main():
 indices=list(map(int,sys.argv[1:])) or [3,4,5,6,7]
 for j in indices:
  order=list(range(3,12));order[j],order[j+1]=order[j+1],order[j];run(order)
if __name__=='__main__':main()
