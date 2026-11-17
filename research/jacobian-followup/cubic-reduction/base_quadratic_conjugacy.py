#!/usr/bin/env python3
"""Full weighted homogeneous quadratic elementary conjugacies in each base axis.

New multi-atom family; does not rerun the saved monomial-line screen.
Primary Singular bases are not independent certificates or global optima.
"""
from pathlib import Path
from itertools import combinations_with_replacement
import hashlib,json,subprocess,time,sys
import sympy as s
HERE=Path(__file__).resolve().parent
W=[1,-1,-2,0,-3,2,-1,-2,-2,-1,-1,0]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
def run(i):
 folder=HERE/'base-quadratic-conjugacy'/('axis-%02d'%(i+1));folder.mkdir(parents=True,exist_ok=True)
 source=HERE/'cubic12-integer-witness.json';raw=json.loads(source.read_text());X=s.symbols('x1:13')
 F=[sum(s.Rational(t['coefficient'])*s.prod(x**e for x,e in zip(X,t['exponent'])) for t in row) for row in raw['coefficients']]
 pairs=[(j,k) for j,k in combinations_with_replacement([a for a in range(12) if a!=i],2) if W[j]+W[k]==W[i]]
 Z=s.symbols('z1:'+str(len(pairs)+1));q=sum(z*X[j]*X[k] for z,(j,k) in zip(Z,pairs));start=time.monotonic()
 checkpoint={'status':'preparing','axis':i+1,'source_sha256':sha(Path(__file__)),'input_sha256':sha(source),'atoms':[[j+1,k+1] for j,k in pairs],'parameter_count':len(Z),'solver_timeout_seconds':20,'scope':'Complete weight-compatible homogeneous quadratic q, independent of this one base axis, in H^-1 K H. Primary experiment only; no full affine-orbit/dimension result.'}
 dump(folder/'checkpoint.json',checkpoint);print('PREPARE',i+1,len(Z),flush=True)
 # Exact elementary inverse is H^-1_i=xi-q, since q is independent of xi.
 FH=[s.Poly(s.expand(f.subs(X[i],X[i]+q)),*X) for f in F]
 target=s.Poly(0,*X)
 for z,(j,k) in zip(Z,pairs):target+=z*FH[j]*FH[k]
 G=FH[:];G[i]=G[i]-target
 high=[];seen=set();low=[]
 for ri,row in enumerate(G):
  for ex,co in row.terms():
   if not co:continue
   co=s.Poly(s.expand(co),*Z,domain=s.QQ)
   if sum(ex)>3:
    co=co.monic();key=str(co.as_expr())
    if key not in seen:seen.add(key);high.append(co.as_expr())
   elif sum(ex)>=2:low.append({'row':ri+1,'exponent':list(ex),'coefficient':str(co.as_expr())})
 inp={'axis':i+1,'weights':W,'variables':list(map(str,Z)),'atoms':[[j+1,k+1] for j,k in pairs],'q':str(q),'high_equations':[str(h) for h in high],'low_coefficients':low,'input_sha256':sha(source),'producer_sha256':sha(Path(__file__))}
 dump(folder/'input.json',inp)
 def sing(e):
  terms=[]
  for ex,co in s.Poly(e,*Z).terms():
   mon='*'.join(str(z)+'^'+str(a) for z,a in zip(Z,ex) if a) or '1'
   terms.append('('+str(co)+')*'+mon)
  return '+'.join(terms) or '0'
 lines=['ring r=0,('+','.join(map(str,Z))+'),dp;','short=0;','ideal I='+','.join(sing(e) for e in high)+';','ideal B=std(I);','print("BASIS_BEGIN");','print(B);','print("BASIS_END");','print("DIM:"+string(dim(B)));']
 for z in Z:lines+=['print("REDUCE_'+str(z)+':"+string(reduce('+str(z)+',B)));']
 lines+=['quit;'];script=folder/'model.sing';script.write_text('\n'.join(lines)+'\n')
 checkpoint.update(status='solver-running',high_equations=len(high),low_slots=len(low),preparation_seconds=round(time.monotonic()-start,3),model_sha256=sha(script));dump(folder/'checkpoint.json',checkpoint);print('SOLVE',i+1,len(high),len(low),flush=True)
 try:
  result=subprocess.run(['Singular','-q',str(script)],capture_output=True,text=True,timeout=20)
  log=result.stdout+result.stderr;status='primary-basis-complete' if result.returncode==0 and '?' not in log else 'error'
 except subprocess.TimeoutExpired as e:
  log=(e.stdout or b'').decode() if isinstance(e.stdout,bytes) else str(e.stdout or '');status='timeout'
 (folder/'model.log').write_text(log);checkpoint.update(status=status,seconds=round(time.monotonic()-start,3));dump(folder/'checkpoint.json',checkpoint)
 print('END',i+1,status,log[:8000],flush=True)
def main():
 for i in map(int,sys.argv[1:] or ['0','1','2']):run(i)
if __name__=='__main__':main()
